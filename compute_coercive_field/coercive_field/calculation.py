"""Coercive-field calculations."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .branch import build_branch_images
from .constants import EV_PER_EANG_MV_CM, UC_CM2_PER_EANG_A3
from .directions import resolve_direction
from .models import PathResult
from .vasp_io import parse_structure, read_raw_images


def resolve_final_polarization(text: str, initial_polarization: np.ndarray) -> np.ndarray:
    if text == "opposite-initial":
        return -initial_polarization
    return np.array([float(part.strip()) for part in text.split(",")], dtype=float)


def target_polarization(result: PathResult, image_index: int) -> np.ndarray:
    image_count = len(result.raw_images)
    progress = image_index / (image_count - 1) if image_count > 1 else 0.0
    return result.target_initial_polarization_uC_cm2 + progress * (
        result.target_final_polarization_uC_cm2 - result.target_initial_polarization_uC_cm2
    )


def switching_axis(result: PathResult) -> np.ndarray:
    axis = result.target_initial_polarization_uC_cm2
    norm = float(np.linalg.norm(axis))
    if norm == 0.0:
        return result.branch_direction
    return axis / norm


def target_delta_dipoles(result: PathResult) -> np.ndarray:
    initial_target = target_polarization(result, 0)
    rows = []
    for image_index, raw in enumerate(result.raw_images):
        delta_p = target_polarization(result, image_index) - initial_target
        rows.append(delta_p * raw.volume_A3 / UC_CM2_PER_EANG_A3)
    return np.array(rows, dtype=float)


def corrected_relative_energies(result: PathResult, field_MV_cm: float) -> np.ndarray:
    rel = np.array([img.rel_energy_eV for img in result.raw_images])
    projected_delta = target_delta_dipoles(result) @ result.field_direction
    return rel - EV_PER_EANG_MV_CM * field_MV_cm * projected_delta


def barrier_at_field(result: PathResult, field_MV_cm: float) -> float:
    corrected = corrected_relative_energies(result, field_MV_cm)
    if result.barrier_mode == "to-saddle":
        window = corrected[1 : result.saddle_image + 1]
    elif result.barrier_mode == "intermediate":
        window = corrected[1:-1]
    else:
        window = corrected[1:]
    if len(window) == 0:
        return 0.0
    return float(np.max(window) - corrected[0])


def analyze_path(
    label: str,
    path: Path,
    energy_source: str,
    field_direction_arg: str,
    field_coordinates: str,
    branch_direction_arg: str,
    branch_coordinates: str,
    branch_mode: str,
    barrier_mode: str,
    initial_polarization: np.ndarray,
    final_polarization_arg: str,
    continuity_threshold_uC_cm2: float,
    branch_search_radius: int,
) -> PathResult:
    raw_images, notes = read_raw_images(label, path, energy_source)
    reference_lattice, _ = parse_structure(raw_images[0].structure_path)

    if field_direction_arg == "reverse-initial":
        raw_field_direction = -initial_polarization
        field_direction = raw_field_direction / np.linalg.norm(raw_field_direction)
    else:
        field_direction = resolve_direction(field_direction_arg, field_coordinates, reference_lattice)

    if branch_direction_arg == "field":
        branch_direction = field_direction
    elif branch_direction_arg == "reverse-initial":
        raw_branch_direction = -initial_polarization
        branch_direction = raw_branch_direction / np.linalg.norm(raw_branch_direction)
    else:
        branch_direction = resolve_direction(branch_direction_arg, branch_coordinates, reference_lattice)

    rel_energies = np.array([img.rel_energy_eV for img in raw_images])
    saddle_image = int(np.argmax(rel_energies))
    final_polarization = resolve_final_polarization(final_polarization_arg, initial_polarization)
    branch_images = build_branch_images(
        raw_images=raw_images,
        branch_mode=branch_mode,
        branch_direction=branch_direction,
        center_image=saddle_image,
        continuity_threshold_uC_cm2=continuity_threshold_uC_cm2,
        search_radius=branch_search_radius,
        initial_polarization=initial_polarization,
        final_polarization=final_polarization,
    )

    target_result = SimpleNamespace(
        raw_images=raw_images,
        target_initial_polarization_uC_cm2=initial_polarization,
        target_final_polarization_uC_cm2=final_polarization,
        branch_direction=branch_direction,
    )
    projected_delta = target_delta_dipoles(target_result) @ field_direction
    if barrier_mode == "to-saddle":
        candidate_indices = range(1, saddle_image + 1)
    elif barrier_mode == "intermediate":
        candidate_indices = range(1, len(raw_images) - 1)
    elif barrier_mode == "full-path":
        candidate_indices = range(1, len(raw_images))
    else:
        raise ValueError(f"unknown barrier mode: {barrier_mode}")

    requirements: list[tuple[float, int]] = []
    impossible: list[int] = []
    for image in candidate_indices:
        de = float(rel_energies[image])
        dm = float(projected_delta[image])
        if de <= 0:
            continue
        if dm <= 0:
            impossible.append(image)
        else:
            requirements.append((de / (EV_PER_EANG_MV_CM * dm), image))

    if impossible:
        coercive_field = math.inf
        limiting_image = None
        notes.append(
            "no finite field for images "
            + ",".join(str(image) for image in impossible)
            + f" in {barrier_mode} mode"
        )
    elif requirements:
        coercive_field, limiting_image = max(requirements, key=lambda item: item[0])
        coercive_field = float(coercive_field)
        limiting_image = int(limiting_image)
    else:
        coercive_field = 0.0
        limiting_image = None

    discontinuous = [img.image for img in branch_images if not img.component_continuity_ok]
    if discontinuous:
        notes.append(
            "branch component jump exceeds threshold at images "
            + ",".join(str(image) for image in discontinuous)
        )

    return PathResult(
        label=label,
        path=path,
        raw_images=raw_images,
        branch_images=branch_images,
        field_direction=field_direction,
        branch_direction=branch_direction,
        target_initial_polarization_uC_cm2=initial_polarization,
        target_final_polarization_uC_cm2=final_polarization,
        branch_mode=branch_mode,
        barrier_mode=barrier_mode,
        continuity_threshold_uC_cm2=continuity_threshold_uC_cm2,
        saddle_image=saddle_image,
        saddle_energy_eV=float(rel_energies[saddle_image]),
        coercive_field_MV_cm=coercive_field,
        limiting_image=limiting_image,
        endpoint_rel_energy_eV=float(rel_energies[-1]),
        endpoint_projected_dipole_eA=float(projected_delta[-1]),
        notes=notes,
    )
