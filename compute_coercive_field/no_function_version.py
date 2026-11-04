from __future__ import annotations

import argparse
import itertools
import math
import re
from pathlib import Path

import numpy as np

from coercive_field.constants import EV_PER_EANG_MV_CM, UC_CM2_PER_EANG_A3
from coercive_field.models import BranchImage, PathResult
from coercive_field.outputs import write_outputs
from coercive_field.vasp_io import read_raw_images


parser = argparse.ArgumentParser(description="No-custom-function coercive-field workflow.")
parser.add_argument("--root", type=Path, default=Path("Pmn2_1") / "xyz-101-000")
parser.add_argument("--sample", action="append", default=[])
parser.add_argument("--energy-source", choices=["oszicar", "neb", "auto"], default="oszicar")
parser.add_argument("--initial-polarization", default="auto")
parser.add_argument("--final-polarization", default="auto")
parser.add_argument("--continuity-threshold", default="auto")
parser.add_argument("--field-direction", default="0,0,-1")
parser.add_argument("--field-coordinates", choices=["lattice", "cartesian"], default="lattice")
parser.add_argument("--barrier-mode", choices=["to-saddle", "intermediate", "full-path"], default="full-path")
parser.add_argument("--coercive-criterion", choices=["barrier", "downhill"], default="downhill")
parser.add_argument("--field-step", type=float, default=0.1)
parser.add_argument("--field-max", type=float, default=None)
parser.add_argument("--output-dir", type=Path, default=Path("results_no_function"))
args = parser.parse_args()

if args.sample:
    samples = []
    for item in args.sample:
        if "=" not in item:
            raise argparse.ArgumentTypeError("--sample must be LABEL=PATH/to/sample")
        label_text, path_text = item.split("=", 1)
        samples.append((label_text, Path(path_text)))
else:
    root_name = args.root.name
    samples = [
        ("HfO2", args.root),
        ("a-HZO", args.root / f"a-{root_name}"),
        ("b-HZO", args.root / f"b-{root_name}"),
        ("c-HZO", args.root / f"c-{root_name}"),
    ]

results = []

for label, path in samples:
    raw_images, notes = read_raw_images(label, path, args.energy_source)
    image_count = len(raw_images)
    center_image = image_count // 2
    reference_lattice = raw_images[0].lattice_A

    raw_polarizations = []
    for image in raw_images:
        raw_polarizations.append(image.raw_polarization_uC_cm2)

    quantum_vectors = []
    for image in raw_images:
        quantum_vectors.append(image.quantum_vectors_uC_cm2)

    phase_names = {"pca21", "pca2_1", "pmn21", "pmn2_1", "fe1", "fe2"}
    phase = None
    state = None
    path_parts = path.parts
    for part_index, part in enumerate(path_parts):
        if part.lower() not in phase_names:
            continue
        phase = part
        for candidate in path_parts[part_index + 1 :]:
            if re.search(r"[xyz]{3}-[01]{3}-[01]{3}", candidate):
                state = re.search(r"([xyz]{3}-[01]{3}-[01]{3})", candidate).group(1)
                break
        break
    if phase is None or state is None:
        raise ValueError(f"could not identify phase and state from {path}")

    max_quantum = 0.0
    for raw in raw_images:
        for coeff_tuple in itertools.product((-1, 0, 1), repeat=3):
            if coeff_tuple == (0, 0, 0):
                continue
            quantum = np.array(coeff_tuple, dtype=float) @ raw.quantum_vectors_uC_cm2
            max_quantum = max(max_quantum, float(np.linalg.norm(quantum)))

    if str(args.continuity_threshold).strip().lower() == "auto":
        continuity_threshold = max_quantum / (2.0 * image_count)
        notes.append(
            f"auto branch continuity threshold {continuity_threshold:.6g} uC/cm^2 "
            "from max(|integer-combined polarization quantum|)/(2*image_count)"
        )
    else:
        continuity_threshold = float(args.continuity_threshold)

    direction_vector = np.array([float(part.strip()) for part in args.field_direction.split(",")], dtype=float)
    if len(direction_vector) != 3 or float(np.linalg.norm(direction_vector)) == 0.0:
        raise ValueError("--field-direction must be a nonzero x,y,z vector")
    if args.field_coordinates == "cartesian":
        field_raw = direction_vector
    else:
        field_raw = direction_vector @ reference_lattice
    field_direction = field_raw / float(np.linalg.norm(field_raw))

    branch_pols = [None] * image_count
    branch_indices = [None] * image_count
    coeff = np.linalg.solve(quantum_vectors[center_image].T, np.zeros(3) - raw_polarizations[center_image])
    branch_indices[center_image] = np.floor(coeff + 0.5).astype(int)
    branch_pols[center_image] = raw_polarizations[center_image] + branch_indices[center_image] @ quantum_vectors[center_image]

    for image in range(center_image - 1, -1, -1):
        coeff = np.linalg.solve(quantum_vectors[image].T, branch_pols[image + 1] - raw_polarizations[image])
        branch_indices[image] = np.floor(coeff + 0.5).astype(int)
        branch_pols[image] = raw_polarizations[image] + branch_indices[image] @ quantum_vectors[image]

    for image in range(center_image + 1, image_count):
        coeff = np.linalg.solve(quantum_vectors[image].T, branch_pols[image - 1] - raw_polarizations[image])
        branch_indices[image] = np.floor(coeff + 0.5).astype(int)
        branch_pols[image] = raw_polarizations[image] + branch_indices[image] @ quantum_vectors[image]

    initial_text = str(args.initial_polarization).strip().lower()
    if initial_text in {"auto", "vasp", "from-vasp"}:
        if phase.lower() in {"pca21", "pca2_1", "fe1"}:
            initial_lattice_direction = np.array([0.0, 0.0, 1.0])
        elif phase.lower() in {"pmn21", "pmn2_1", "fe2"}:
            initial_lattice_direction = np.array([1.0, 0.0, 1.0])
        else:
            raise ValueError(f"unknown phase: {phase}")
        initial_direction_raw = initial_lattice_direction @ reference_lattice
        initial_direction = initial_direction_raw / float(np.linalg.norm(initial_direction_raw))
        initial_magnitude = float(np.linalg.norm(branch_pols[0]))
        initial_target = initial_magnitude * initial_direction
        coeff = np.linalg.solve(quantum_vectors[0].T, initial_target - raw_polarizations[0])
        initial_index = np.floor(coeff + 0.5).astype(int)
        initial_p = raw_polarizations[0] + initial_index @ quantum_vectors[0]
    else:
        initial_p = np.array([float(part.strip()) for part in args.initial_polarization.split(",")], dtype=float)
        initial_norm = float(np.linalg.norm(initial_p))
        initial_direction = initial_p / initial_norm if initial_norm != 0.0 else np.zeros(3)

    notes.append(
        "initial polarization "
        + ",".join(f"{value:.6g}" for value in initial_p)
        + " uC/cm^2 from direction "
        + ",".join(f"{value:.6g}" for value in initial_direction)
    )

    final_text = str(args.final_polarization).strip().lower()
    if final_text == "auto":
        prefix, suffix, translation = state.split("-")
        del translation
        mirror_matrix = np.diag([-1 if value == "1" else 1 for value in suffix])
        permutation_matrix = np.zeros((3, 3), dtype=int)
        for row_index, axis_name in enumerate(prefix):
            permutation_matrix[row_index, ord(axis_name) - ord("x")] = 1
        rotation_matrix = mirror_matrix @ permutation_matrix
        final_p = rotation_matrix @ initial_p
    elif final_text == "opposite-initial":
        final_p = -initial_p
    else:
        final_p = np.array([float(part.strip()) for part in args.final_polarization.split(",")], dtype=float)

    initial_norm = float(np.linalg.norm(initial_p))
    final_norm = float(np.linalg.norm(final_p))
    non_180_switch = False
    if initial_norm != 0.0 and final_norm != 0.0:
        cosine = float((initial_p @ final_p) / (initial_norm * final_norm))
        non_180_switch = not np.isclose(cosine, -1.0, atol=1e-3)

    if non_180_switch:
        coeff = np.linalg.solve(quantum_vectors[-1].T, final_p - raw_polarizations[-1])
        desired_final_index = np.floor(coeff + 0.5).astype(int)
        desired_final_pol = raw_polarizations[-1] + desired_final_index @ quantum_vectors[-1]
        desired_score = float(np.linalg.norm(desired_final_pol - final_p))
        current_score = float(np.linalg.norm(branch_pols[-1] - final_p))
        offset = desired_final_index - branch_indices[-1]
        if desired_score < current_score and np.any(offset != 0):
            proposed_pols = list(branch_pols)
            proposed_indices = list(branch_indices)
            for image in range(center_image + 1, image_count):
                proposed_pols[image] = proposed_pols[image] + offset @ quantum_vectors[image]
                proposed_indices[image] = proposed_indices[image] + offset
            offset_is_continuous = True
            for image in range(center_image + 1, image_count):
                jump = float(np.max(np.abs(proposed_pols[image] - proposed_pols[image - 1])))
                if jump > continuity_threshold:
                    offset_is_continuous = False
                    break
            if offset_is_continuous:
                branch_pols = proposed_pols
                branch_indices = proposed_indices
        notes.append(
            "non-180 switch detected; final polarization "
            + ",".join(f"{value:.6g}" for value in final_p)
            + " uC/cm^2 used as endpoint branch reference"
        )

    branch_images = []
    previous_pol = None
    previous_projection = None
    for raw, pol, index in zip(raw_images, branch_pols, branch_indices):
        branch_dipole = pol * raw.volume_A3 / UC_CM2_PER_EANG_A3
        projection = float(pol @ field_direction)
        if previous_pol is None:
            delta_pol = np.zeros(3)
            delta_projection = 0.0
            max_jump = 0.0
            monotonic = True
            continuous = True
        else:
            delta_pol = pol - previous_pol
            delta_projection = projection - float(previous_projection)
            max_jump = float(np.max(np.abs(delta_pol)))
            monotonic = delta_projection >= -1e-8
            continuous = max_jump <= continuity_threshold
        branch_images.append(
            BranchImage(
                image=raw.image,
                branch_index=index,
                branch_polarization_uC_cm2=pol,
                branch_dipole_eA=branch_dipole,
                branch_projection_uC_cm2=projection,
                delta_branch_polarization_from_prev_uC_cm2=delta_pol,
                delta_branch_projection_from_prev_uC_cm2=delta_projection,
                max_abs_component_jump_from_prev_uC_cm2=max_jump,
                monotonic_along_branch=monotonic,
                component_continuity_ok=continuous,
            )
        )
        previous_pol = pol
        previous_projection = projection

    rel_energies = np.array([image.rel_energy_eV for image in raw_images], dtype=float)
    saddle_image = int(np.argmax(rel_energies))
    first_branch_dipole = branch_images[0].branch_dipole_eA
    projected_delta = np.array(
        [(image.branch_dipole_eA - first_branch_dipole) @ field_direction for image in branch_images],
        dtype=float,
    )

    discontinuous_images = [image.image for image in branch_images if not image.component_continuity_ok]
    if discontinuous_images:
        notes.append(
            "branch component jump exceeds threshold at images "
            + ",".join(str(image) for image in discontinuous_images)
            + "; coercive field is not reliable until branch correction is fixed"
        )
        coercive_field = math.inf
        limiting_image = None
    else:
        if args.barrier_mode == "to-saddle":
            candidate_indices = range(1, saddle_image + 1)
        elif args.barrier_mode == "intermediate":
            candidate_indices = range(1, image_count - 1)
        else:
            candidate_indices = range(1, image_count)

        requirements = []
        impossible = []
        if args.coercive_criterion == "barrier":
            for image in candidate_indices:
                de = float(rel_energies[image] - rel_energies[0])
                dm = float(projected_delta[image] - projected_delta[0])
                if de <= 0:
                    continue
                if dm <= 0:
                    impossible.append(image)
                else:
                    requirements.append((de / (EV_PER_EANG_MV_CM * dm), image))
        else:
            for image in candidate_indices:
                de = float(rel_energies[image] - rel_energies[image - 1])
                dm = float(projected_delta[image] - projected_delta[image - 1])
                if de <= 0:
                    continue
                if dm <= 0:
                    impossible.append(image)
                else:
                    requirements.append((de / (EV_PER_EANG_MV_CM * dm), image))

        if impossible:
            if args.coercive_criterion == "downhill":
                impossible_text = ",".join(f"{image - 1}->{image}" for image in impossible)
            else:
                impossible_text = ",".join(str(image) for image in impossible)
            notes.append(
                "no finite field for "
                + args.coercive_criterion
                + " constraints "
                + impossible_text
                + f" in {args.barrier_mode} mode"
            )
            coercive_field = math.inf
            limiting_image = None
        elif requirements:
            coercive_field, limiting_image = max(requirements, key=lambda item: item[0])
            coercive_field = float(coercive_field)
        else:
            coercive_field = 0.0
            limiting_image = None

    results.append(
        PathResult(
            label=label,
            path=path,
            raw_images=raw_images,
            branch_images=branch_images,
            field_direction=field_direction,
            branch_direction=field_direction,
            target_initial_polarization_uC_cm2=initial_p,
            target_final_polarization_uC_cm2=final_p,
            branch_mode="center-zero",
            barrier_mode=args.barrier_mode,
            coercive_criterion=args.coercive_criterion,
            continuity_threshold_uC_cm2=continuity_threshold,
            saddle_image=saddle_image,
            saddle_energy_eV=float(rel_energies[saddle_image]),
            coercive_field_MV_cm=coercive_field,
            limiting_image=limiting_image,
            endpoint_rel_energy_eV=float(rel_energies[-1]),
            endpoint_projected_dipole_eA=float(projected_delta[-1]),
            notes=notes,
        )
    )

write_outputs(results, args.output_dir, args.field_step, args.field_max)

print("label, zero-field barrier (eV), coercive field (MV/cm), limiting constraint")
for result in results:
    ec_text = "inf" if not math.isfinite(result.coercive_field_MV_cm) else f"{result.coercive_field_MV_cm:.4f}"
    limiting = result.limiting_image
    if result.coercive_criterion == "downhill" and result.limiting_image is not None:
        limiting = f"{result.limiting_image - 1}->{result.limiting_image}"
    print(f"{result.label}, {result.saddle_energy_eV:.6f}, {ec_text}, {limiting}")
    for note in result.notes:
        print(f"  note: {note}")
print(f"wrote outputs to {args.output_dir.resolve()}")
