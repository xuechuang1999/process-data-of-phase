#!/usr/bin/env python3
import re
from pathlib import Path
from dataclasses import dataclass
from typing import Tuple, Union

import numpy as np


# 1 e / Å^2 = 1602.176634 μC / cm^2
E_PER_A2_TO_UC_CM2 = 1602.176634

PathLike = Union[str, Path]


@dataclass
class PolarizationResult:
    directory: Path

    lattice: np.ndarray          # shape (3, 3), rows are a, b, c, unit Å
    volume: float                # Å^3

    p_ion: np.ndarray            # e Å
    p_elc: np.ndarray            # e Å
    p_tot: np.ndarray            # e Å
    p_shifted: np.ndarray        # e Å

    P_raw: np.ndarray            # μC/cm^2
    P_shifted: np.ndarray        # μC/cm^2

    Q: np.ndarray                # shape (3, 3), rows are Qa, Qb, Qc, unit μC/cm^2
    branch_n: np.ndarray         # [na, nb, nc]


def read_lattice_from_poscar(poscar: PathLike) -> np.ndarray:
    """
    Read lattice vectors from POSCAR/CONTCAR.

    Returns
    -------
    lattice : np.ndarray
        Shape (3, 3), rows are lattice vectors a, b, c in Å.
    """
    poscar = Path(poscar)
    lines = poscar.read_text().splitlines()

    scale_values = [float(x) for x in lines[1].split()]
    lattice_raw = np.array(
        [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)],
        dtype=float,
    )

    if len(scale_values) == 1:
        scale = scale_values[0]

        # VASP convention:
        # If scale < 0, abs(scale) is the desired cell volume.
        if scale < 0:
            target_volume = abs(scale)
            raw_volume = abs(np.linalg.det(lattice_raw))
            scale = (target_volume / raw_volume) ** (1.0 / 3.0)

        lattice = lattice_raw * scale

    elif len(scale_values) == 3:
        # Three scale factors. Usually rare.
        # Here they are applied to a, b, c respectively.
        scales = np.array(scale_values, dtype=float)
        lattice = lattice_raw * scales[:, None]

    else:
        raise ValueError(f"Cannot parse scale line in {poscar}: {lines[1]}")

    return lattice


def cell_volume(lattice: np.ndarray) -> float:
    """
    Calculate cell volume from lattice vectors.

    lattice rows are a, b, c.
    """
    a, b, c = lattice
    return abs(np.dot(a, np.cross(b, c)))


def read_vasp_dipole_from_outcar(outcar: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read ionic and electronic dipole moments from OUTCAR.

    Requires LCALCPOL = .TRUE.

    Returns
    -------
    p_ion : np.ndarray
        Ionic dipole moment in e Å.
    p_elc : np.ndarray
        Total electronic dipole moment in e Å.
    """
    outcar = Path(outcar)
    text = outcar.read_text(errors="ignore")

    ion_pattern = re.compile(
        r"Ionic dipole moment:\s*p\[ion\]=\(\s*"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)"
    )

    elc_pattern = re.compile(
        r"Total electronic dipole moment:\s*p\[elc\]=\(\s*"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)"
    )

    ion_matches = ion_pattern.findall(text)
    elc_matches = elc_pattern.findall(text)

    if not ion_matches:
        raise RuntimeError(
            f"Cannot find 'Ionic dipole moment' in {outcar}. "
            "Check whether LCALCPOL = .TRUE. was used."
        )

    if not elc_matches:
        raise RuntimeError(
            f"Cannot find 'Total electronic dipole moment' in {outcar}. "
            "Check whether LCALCPOL = .TRUE. was used."
        )

    # Use the last occurrence, safer for OUTCAR with multiple ionic steps.
    p_ion = np.array([float(x) for x in ion_matches[-1]], dtype=float)
    p_elc = np.array([float(x) for x in elc_matches[-1]], dtype=float)

    return p_ion, p_elc


def polarization_quantum(lattice: np.ndarray, volume: float) -> np.ndarray:
    """
    Calculate polarization quantum basis vectors.

    Q_a = e a / Ω
    Q_b = e b / Ω
    Q_c = e c / Ω

    Returns
    -------
    Q : np.ndarray
        Shape (3, 3).
        Rows are Qa, Qb, Qc.
        Columns are Cartesian x, y, z.
        Unit: μC/cm^2.
    """
    return lattice / volume * E_PER_A2_TO_UC_CM2


def calculate_vasp_polarization(
    directory: PathLike,
    outcar_name: str = "OUTCAR",
    structure_name: str = "CONTCAR",
    branch_n: Tuple[int, int, int] = (0, 0, 0),
) -> PolarizationResult:
    """
    Calculate polarization from VASP OUTCAR and CONTCAR/POSCAR.

    Parameters
    ----------
    directory : str or Path
        VASP calculation directory.
    outcar_name : str
        OUTCAR filename.
    structure_name : str
        CONTCAR or POSCAR filename.
    branch_n : tuple[int, int, int]
        Polarization branch shift integers: (na, nb, nc).

    Returns
    -------
    PolarizationResult
    """
    directory = Path(directory)

    outcar = directory / outcar_name
    structure = directory / structure_name

    if not outcar.exists():
        raise FileNotFoundError(f"OUTCAR not found: {outcar}")

    if not structure.exists():
        raise FileNotFoundError(f"Structure file not found: {structure}")

    lattice = read_lattice_from_poscar(structure)
    volume = cell_volume(lattice)

    p_ion, p_elc = read_vasp_dipole_from_outcar(outcar)
    p_tot = p_ion + p_elc

    Q = polarization_quantum(lattice, volume)

    branch_n_array = np.array(branch_n, dtype=int)

    # In dipole moment unit e Å:
    # p' = p + na*a + nb*b + nc*c
    p_shifted = p_tot + branch_n_array @ lattice

    # Convert e Å / Å^3 to μC/cm^2.
    P_raw = p_tot / volume * E_PER_A2_TO_UC_CM2
    P_shifted = p_shifted / volume * E_PER_A2_TO_UC_CM2

    return PolarizationResult(
        directory=directory,
        lattice=lattice,
        volume=volume,
        p_ion=p_ion,
        p_elc=p_elc,
        p_tot=p_tot,
        p_shifted=p_shifted,
        P_raw=P_raw,
        P_shifted=P_shifted,
        Q=Q,
        branch_n=branch_n_array,
    )