"""Readers for VASP NEB and LCALCPOL outputs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .constants import UC_CM2_PER_EANG_A3
from .models import RawImage


def parse_structure(path: Path) -> tuple[np.ndarray, float]:
    """Read lattice vectors and cell volume from POSCAR/CONTCAR.

    Formula:
        Omega = |det(a, b, c)|

    Meaning of symbols:
        a,b,c: lattice rows read from POSCAR/CONTCAR after applying scale.
        Omega: returned volume in Angstrom^3.
    """

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
    """Return the relaxed structure file for one NEB image.

    CONTCAR is preferred because it is the relaxed image structure used for
    volume and polarization conversion. POSCAR is only a fallback.
    返回一个NEB image的结构文件路径，优先使用CONTCAR文件，如果不存在则使用POSCAR文件
    """

    image_dir = path / f"{image:02d}"
    contcar = image_dir / "CONTCAR"

    if contcar.exists():
        return contcar
    dipol_dir = path / "1-dipol" / f"{image:02d}"
    if dipol_dir.exists():
        dipol_poscar = dipol_dir / "POSCAR"
        if dipol_poscar.exists():
            return dipol_poscar
    poscar = image_dir / "POSCAR"
    if poscar.exists():
        return poscar
    raise FileNotFoundError(f"could not find CONTCAR or POSCAR in {image_dir} or {dipol_dir}")


def parse_oszicar_energy(path: Path) -> float:
    """Read the last E0 value from OSZICAR.

    Meaning:
        E0 is the zero-field electronic/ionic energy printed by VASP for the
        final ionic step in this image directory. Relative energies are built
        later as U_i = E_i - E_0.
    """

    text = path.read_text(errors="ignore")
    matches = re.findall(r"E0=\s*([-+0-9.Ee]+)", text)
    if not matches:
        raise ValueError(f"could not find E0 in {path}")
    return float(matches[-1])


def parse_neb_dat(path: Path) -> np.ndarray:
    """
    由于当运行了nebresults.pl脚本后，neb.dat文件中会包含每个image的相对能量信息，所以这里直接从neb.dat文件的第三列中读取相对能量值。
    """

    values: list[float] = []
    for line in path.read_text(errors="ignore").splitlines():
        fields = line.split()
        if len(fields) >= 3:
            values.append(float(fields[2]))
    if not values:
        raise ValueError(f"could not parse relative energies from {path}")
    return np.array(values, dtype=float)


def parse_lcalcpol_dipoles(outcar: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read ionic and electronic dipoles printed by VASP LCALCPOL.

    VASP prints both terms in e*Angstrom. The total raw dipole used later is
    p_ion + p_elc, before any Berry-branch translation is applied.

    Formula:
        mu_raw_i = p_ion_i + p_elc_i

    Meaning of symbols:
        p_ion_i: returned ion vector, ionic dipole of image i.
        p_elc_i: returned electronic vector, electronic dipole of image i.
        mu_raw_i: raw_dipole_eA computed in read_raw_images. 总的原始偶极矩是离子偶极矩和电子偶极矩的和。
    
    根据传入的outcar路径，读取OUTCAR文件中的离子和电子偶极矩，并返回两个numpy数组，分别表示离子偶极矩和电子偶极矩，单位为e*Angstrom。
    """

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
    """Return the three polarization quantum vectors in uC/cm^2.

    Translating the dipole by one lattice vector changes polarization by
    UC_CM2_PER_EANG_A3*a/Omega, UC_CM2_PER_EANG_A3*b/Omega, or
    UC_CM2_PER_EANG_A3*c/Omega.

    Formula:
        Q_a = UC_CM2_PER_EANG_A3*a/Omega
        Q_b = UC_CM2_PER_EANG_A3*b/Omega
        Q_c = UC_CM2_PER_EANG_A3*c/Omega

    Meaning of symbols:
        a,b,c: lattice vectors in Angstrom.
        Omega: volume_A3, cell volume in Angstrom^3.
        Q_a,Q_b,Q_c: returned polarization quantum vectors in uC/cm^2.
    """

    return lattice * UC_CM2_PER_EANG_A3 / volume_A3


def discover_image_count(path: Path) -> int:
    """Count NEB image directories named 00, 01, ..."""

    image_dirs = sorted(p for p in path.iterdir() if p.is_dir() and re.fullmatch(r"\d\d", p.name))
    # image_dirs = [f"{i:02d}" for i in range(2, 7)]
    if not image_dirs:
        raise ValueError(f"could not find image directories in {path}")
    return len(image_dirs)


def read_relative_energies(path: Path, source: str, image_count: int) -> tuple[np.ndarray, list[str]]:
    """Return zero-field relative energies for all images.

    OSZICAR is always parsed so that neb.dat can be checked against it. If
    source is "neb", neb.dat must exist and have the expected length. If
    source is "auto", neb.dat is used when valid; otherwise OSZICAR is used.

    Formula:
        U_i = E_i - E_0

    Meaning of symbols:
        E_i: absolute OSZICAR E0 value of image i.
        E_0: OSZICAR E0 value of image 00.
        U_i: returned relative zero-field energy, in eV.
    """

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
        return neb_rel, notes, oszicar_energies[0]
    if source == "auto" and neb_rel is not None:
        return neb_rel, notes, oszicar_energies[0]
    return oszicar_rel, notes, oszicar_energies[0]


def read_raw_images(label: str, path: Path, energy_source: str) -> tuple[list[RawImage], list[str]]:
    """Read all raw per-image data needed by the coercive-field workflow.

    Main formulas per image i:
        mu_raw_i = p_ion_i + p_elc_i
        P_raw_i = UC_CM2_PER_EANG_A3 * mu_raw_i / Omega_i
        Q_i = UC_CM2_PER_EANG_A3 * lattice_i / Omega_i

    Meaning of symbols:
        mu_raw_i: raw total dipole in e*Angstrom.
        P_raw_i: raw polarization before Berry-branch correction, in uC/cm^2.
        Omega_i: image cell volume in Angstrom^3.
        Q_i: 3x3 matrix of polarization quantum vectors for image i.
    """

    image_count = discover_image_count(path)
    rel_energies, notes, reference_energy = read_relative_energies(path, energy_source, image_count) 

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
                image=image, # image为image的序号，从0开始计数
                energy_path=path / f"{image:02d}" / "OSZICAR",
                structure_path=structure_path,
                dipole_outcar_path=path / "1-dipol" / f"{image:02d}" / "OUTCAR",
                energy_eV=float(rel_energies[image] + reference_energy),
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
