"""Readers for VASP NEB and LCALCPOL outputs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .constants import UC_CM2_PER_EANG_A3
from .models import RawImage


def parse_structure(path: Path) -> tuple[np.ndarray, float]:
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 8:
        raise ValueError(f"{path} is too short to be a VASP structure file")
    scale = float(lines[1].split()[0])
    lattice = np.array(
        [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)],
        dtype=float,
    )
    lattice *= scale
    volume = abs(float(np.linalg.det(lattice)))
    return lattice, volume


def image_structure_path(path: Path, image: int) -> Path:
    image_dir = path / f"{image:02d}"
    contcar = image_dir / "CONTCAR"
    if contcar.exists():
        return contcar
    poscar = image_dir / "POSCAR"
    if poscar.exists():
        return poscar
    raise FileNotFoundError(f"could not find CONTCAR or POSCAR in {image_dir}")


def parse_oszicar_energy(path: Path) -> float:
    text = path.read_text(errors="ignore")
    matches = re.findall(r"E0=\s*([-+0-9.Ee]+)", text)
    if not matches:
        raise ValueError(f"could not find E0 in {path}")
    return float(matches[-1])


def parse_neb_dat(path: Path) -> np.ndarray:
    values: list[float] = []
    for line in path.read_text(errors="ignore").splitlines():
        fields = line.split()
        if len(fields) >= 3:
            values.append(float(fields[2]))
    if not values:
        raise ValueError(f"could not parse relative energies from {path}")
    return np.array(values, dtype=float)


def parse_lcalcpol_dipoles(outcar: Path) -> tuple[np.ndarray, np.ndarray]:
    text = outcar.read_text(errors="ignore")
    ion_matches = re.findall(
        r"Ionic dipole moment: p\[ion\]=\(\s*([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)",
        text,
    )
    elc_matches = re.findall(
        r"Total electronic dipole moment: p\[elc\]=\(\s*([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)",
        text,
    )
    if not ion_matches or not elc_matches:
        raise ValueError(f"could not find LCALCPOL dipoles in {outcar}")
    ion = np.array([float(x) for x in ion_matches[-1]], dtype=float)
    electronic = np.array([float(x) for x in elc_matches[-1]], dtype=float)
    return ion, electronic


def polarization_quanta(lattice: np.ndarray, volume_A3: float) -> np.ndarray:
    return lattice * UC_CM2_PER_EANG_A3 / volume_A3


def discover_image_count(path: Path) -> int:
    image_dirs = sorted(p for p in path.iterdir() if p.is_dir() and re.fullmatch(r"\d\d", p.name))
    if not image_dirs:
        raise ValueError(f"could not find image directories in {path}")
    return len(image_dirs)


def read_relative_energies(path: Path, source: str, image_count: int) -> tuple[np.ndarray, list[str]]:
    notes: list[str] = []
    oszicar_energies = np.array(
        [parse_oszicar_energy(path / f"{i:02d}" / "OSZICAR") for i in range(image_count)],
        dtype=float,
    )
    oszicar_rel = oszicar_energies - oszicar_energies[0]

    neb_path = path / "neb.dat"
    neb_rel: np.ndarray | None = None
    if neb_path.exists():
        neb_rel = parse_neb_dat(neb_path)
        if len(neb_rel) == image_count:
            max_diff = float(np.max(np.abs(neb_rel - oszicar_rel)))
            if max_diff > 1e-3:
                notes.append(f"neb.dat differs from OSZICAR by up to {max_diff:.4f} eV")
        else:
            notes.append(f"neb.dat has {len(neb_rel)} points, expected {image_count}")
            neb_rel = None

    if source == "neb":
        if neb_rel is None:
            raise ValueError(f"{path} has no usable neb.dat")
        return neb_rel, notes
    if source == "auto" and neb_rel is not None:
        return neb_rel, notes
    return oszicar_rel, notes


def read_raw_images(label: str, path: Path, energy_source: str) -> tuple[list[RawImage], list[str]]:
    image_count = discover_image_count(path)
    rel_energies, notes = read_relative_energies(path, energy_source, image_count)
    raw_images: list[RawImage] = []
    for image in range(image_count):
        structure_path = image_structure_path(path, image)
        lattice, volume = parse_structure(structure_path)
        ion_dipole, electronic_dipole = parse_lcalcpol_dipoles(
            path / "1-dipol" / f"{image:02d}" / "OUTCAR"
        )
        raw_dipole = ion_dipole + electronic_dipole
        raw_polarization = raw_dipole * UC_CM2_PER_EANG_A3 / volume
        raw_images.append(
            RawImage(
                label=label,
                sample_path=path,
                image=image,
                energy_path=path / f"{image:02d}" / "OSZICAR",
                structure_path=structure_path,
                dipole_outcar_path=path / "1-dipol" / f"{image:02d}" / "OUTCAR",
                energy_eV=float(rel_energies[image] + rel_energies[0]),
                rel_energy_eV=float(rel_energies[image]),
                lattice_A=lattice,
                volume_A3=volume,
                ion_dipole_eA=ion_dipole,
                electronic_dipole_eA=electronic_dipole,
                raw_dipole_eA=raw_dipole,
                raw_polarization_uC_cm2=raw_polarization,
                quantum_vectors_uC_cm2=polarization_quanta(lattice, volume),
            )
        )
    return raw_images, notes

