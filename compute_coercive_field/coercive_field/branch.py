"""Berry-branch correction logic."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np

from .constants import UC_CM2_PER_EANG_A3
from .models import BranchImage, RawImage


def fold_to_target(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    coeff = np.linalg.solve(quantum_vectors.T, target - raw_polarization)
    center = np.rint(coeff).astype(int)
    best_score = math.inf
    best_pol: np.ndarray | None = None
    best_index: np.ndarray | None = None
    for shift in np.ndindex(3, 3, 3):
        index = center + np.array(shift, dtype=int) - 1
        candidate = raw_polarization + index @ quantum_vectors
        score = float(np.linalg.norm(candidate - target))
        if score < best_score:
            best_score = score
            best_pol = candidate
            best_index = index
    assert best_pol is not None and best_index is not None
    return best_pol, best_index


def candidate_branches(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    center_index: np.ndarray,
    radius: int,
) -> Iterable[tuple[np.ndarray, np.ndarray]]:
    for shift in np.ndindex(*([2 * radius + 1] * 3)):
        index = center_index + np.array(shift, dtype=int) - radius
        yield raw_polarization + index @ quantum_vectors, index


def fold_near_previous(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    previous_polarization: np.ndarray,
    branch_direction: np.ndarray,
    side: str,
    search_radius: int,
) -> tuple[np.ndarray, np.ndarray]:
    coeff = np.linalg.solve(quantum_vectors.T, previous_polarization - raw_polarization)
    center_index = np.rint(coeff).astype(int)
    previous_projection = float(previous_polarization @ branch_direction)
    best: tuple[float, np.ndarray, np.ndarray] | None = None
    fallback: tuple[float, np.ndarray, np.ndarray] | None = None
    for candidate, index in candidate_branches(
        raw_polarization,
        quantum_vectors,
        center_index,
        search_radius,
    ):
        delta = candidate - previous_polarization
        distance = float(np.linalg.norm(delta))
        max_component_jump = float(np.max(np.abs(delta)))
        # Prefer component continuity; Euclidean distance breaks ties.
        score = 10.0 * max_component_jump + distance
        projection = float(candidate @ branch_direction)
        if fallback is None or score < fallback[0]:
            fallback = (score, candidate, index)
        monotonic = projection >= previous_projection - 1e-8 if side == "forward" else projection <= previous_projection + 1e-8
        if monotonic and (best is None or score < best[0]):
            best = (score, candidate, index)
    chosen = best if best is not None else fallback
    assert chosen is not None
    return chosen[1], chosen[2]


def fold_continuous_near_target(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    previous_polarization: np.ndarray,
    target: np.ndarray,
    search_radius: int,
) -> tuple[np.ndarray, np.ndarray]:
    coeff = np.linalg.solve(quantum_vectors.T, previous_polarization - raw_polarization)
    center_index = np.rint(coeff).astype(int)
    best: tuple[float, np.ndarray, np.ndarray] | None = None
    for candidate, index in candidate_branches(
        raw_polarization,
        quantum_vectors,
        center_index,
        search_radius,
    ):
        delta = candidate - previous_polarization
        max_component_jump = float(np.max(np.abs(delta)))
        continuity_distance = float(np.linalg.norm(delta))
        target_distance = float(np.linalg.norm(candidate - target))
        score = 1000.0 * max_component_jump + continuity_distance + 0.01 * target_distance
        if best is None or score < best[0]:
            best = (score, candidate, index)
    assert best is not None
    return best[1], best[2]


def build_branch_images(
    raw_images: list[RawImage],
    branch_mode: str,
    branch_direction: np.ndarray,
    center_image: int,
    continuity_threshold_uC_cm2: float,
    search_radius: int,
    initial_polarization: np.ndarray,
    final_polarization: np.ndarray,
) -> list[BranchImage]:
    raw_pols = [img.raw_polarization_uC_cm2 for img in raw_images]
    quanta = [img.quantum_vectors_uC_cm2 for img in raw_images]
    image_count = len(raw_images)

    branch_pols: list[np.ndarray | None] = [None] * image_count
    branch_indices: list[np.ndarray | None] = [None] * image_count

    if branch_mode == "center-zero":
        center_pol, center_index = fold_to_target(raw_pols[center_image], quanta[center_image], np.zeros(3))
        branch_pols[center_image] = center_pol
        branch_indices[center_image] = center_index
        for image in range(center_image - 1, -1, -1):
            assert branch_pols[image + 1] is not None
            branch_pols[image], branch_indices[image] = fold_near_previous(
                raw_pols[image],
                quanta[image],
                branch_pols[image + 1],
                branch_direction,
                "backward",
                search_radius,
            )
        for image in range(center_image + 1, image_count):
            assert branch_pols[image - 1] is not None
            branch_pols[image], branch_indices[image] = fold_near_previous(
                raw_pols[image],
                quanta[image],
                branch_pols[image - 1],
                branch_direction,
                "forward",
                search_radius,
            )
    elif branch_mode == "endpoint-target":
        initial_target = initial_polarization
        branch_pols[0], branch_indices[0] = fold_to_target(raw_pols[0], quanta[0], initial_target)
        for image in range(1, image_count):
            assert branch_pols[image - 1] is not None
            progress = image / (image_count - 1) if image_count > 1 else 0.0
            target = initial_polarization + progress * (final_polarization - initial_polarization)
            branch_pols[image], branch_indices[image] = fold_continuous_near_target(
                raw_pols[image],
                quanta[image],
                branch_pols[image - 1],
                target,
                search_radius,
            )
    else:
        raise ValueError(f"unknown branch mode: {branch_mode}")

    branch_images: list[BranchImage] = []
    previous_pol: np.ndarray | None = None
    previous_projection: float | None = None
    for raw, pol, index in zip(raw_images, branch_pols, branch_indices):
        assert pol is not None and index is not None
        branch_dipole = pol * raw.volume_A3 / UC_CM2_PER_EANG_A3
        projection = float(pol @ branch_direction)
        if previous_pol is None:
            delta = np.zeros(3)
            delta_projection = 0.0
            max_jump = 0.0
            monotonic = True
            continuous = True
        else:
            delta = pol - previous_pol
            delta_projection = projection - float(previous_projection)
            max_jump = float(np.max(np.abs(delta)))
            monotonic = delta_projection >= -1e-8
            continuous = max_jump <= continuity_threshold_uC_cm2
        branch_images.append(
            BranchImage(
                image=raw.image,
                branch_index=index,
                branch_polarization_uC_cm2=pol,
                branch_dipole_eA=branch_dipole,
                branch_projection_uC_cm2=projection,
                delta_branch_polarization_from_prev_uC_cm2=delta,
                delta_branch_projection_from_prev_uC_cm2=delta_projection,
                max_abs_component_jump_from_prev_uC_cm2=max_jump,
                monotonic_along_branch=monotonic,
                component_continuity_ok=continuous,
            )
        )
        previous_pol = pol
        previous_projection = projection
    return branch_images
