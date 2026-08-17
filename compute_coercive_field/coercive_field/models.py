from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class RawImage:
    label: str
    sample_path: Path
    image: int
    energy_path: Path
    structure_path: Path
    dipole_outcar_path: Path
    energy_eV: float
    rel_energy_eV: float
    lattice_A: np.ndarray
    volume_A3: float
    ion_dipole_eA: np.ndarray
    electronic_dipole_eA: np.ndarray
    raw_dipole_eA: np.ndarray
    raw_polarization_uC_cm2: np.ndarray
    quantum_vectors_uC_cm2: np.ndarray


@dataclass
class BranchImage:
    image: int
    branch_index: np.ndarray
    branch_polarization_uC_cm2: np.ndarray
    branch_dipole_eA: np.ndarray
    branch_projection_uC_cm2: float
    delta_branch_polarization_from_prev_uC_cm2: np.ndarray
    delta_branch_projection_from_prev_uC_cm2: float
    max_abs_component_jump_from_prev_uC_cm2: float
    monotonic_along_branch: bool
    component_continuity_ok: bool


@dataclass
class PathResult:
    label: str
    path: Path
    raw_images: list[RawImage]
    branch_images: list[BranchImage]
    field_direction: np.ndarray
    branch_direction: np.ndarray
    target_initial_polarization_uC_cm2: np.ndarray
    target_final_polarization_uC_cm2: np.ndarray
    branch_mode: str
    barrier_mode: str
    coercive_criterion: str
    continuity_threshold_uC_cm2: float
    saddle_image: int
    saddle_energy_eV: float
    coercive_field_MV_cm: float
    limiting_image: int | None
    endpoint_rel_energy_eV: float
    endpoint_projected_dipole_eA: float
    notes: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class PathSegment:
    segment_index: int

    start_index: int
    maximum_index: int
    end_index: int

    start_energy_eV: float
    maximum_energy_eV: float
    end_energy_eV: float

    forward_obstacle_eV: float
    reverse_obstacle_eV: float
    energy_change_eV: float

    segment_type: str
    maximum_is_internal: bool