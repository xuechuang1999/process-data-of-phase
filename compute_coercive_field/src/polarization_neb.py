#!/usr/bin/env python3
from __future__ import annotations

import re
import itertools
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


E_PER_A2_TO_UC_CM2 = 1602.176634  # e/Å^2 -> μC/cm^2


@dataclass
class ImagePolarization:
    image: int
    image_name: str
    image_dir: Path

    volume: float
    lattice: np.ndarray          # shape (3, 3), rows are a, b, c

    p_ion: np.ndarray            # e Å
    p_elc: np.ndarray            # e Å
    p_tot: np.ndarray            # e Å

    P_raw: np.ndarray            # μC/cm^2
    Q: np.ndarray                # shape (3, 3), rows are Qa, Qb, Qc, μC/cm^2

    branch: np.ndarray | None = None       # [na, nb, nc]
    P_corr: np.ndarray | None = None       # μC/cm^2
    jump_norm: float | None = None         # |P_i - P_{i-1}| after correction


def read_lattice_from_poscar(poscar: str | Path) -> np.ndarray:
    """
    Read lattice vectors from POSCAR/CONTCAR.

    Returns
    -------
    lattice : np.ndarray
        shape (3, 3), rows are lattice vectors a, b, c, unit Å.
    """
    poscar = Path(poscar)
    lines = poscar.read_text(errors="ignore").splitlines()

    scale_values = [float(x) for x in lines[1].split()]
    lattice_raw = np.array(
        [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)],
        dtype=float,
    )

    if len(scale_values) == 1:
        scale = scale_values[0]

        # VASP convention:
        # negative scale means target cell volume.
        if scale < 0:
            target_volume = abs(scale)
            raw_volume = abs(np.linalg.det(lattice_raw))
            scale = (target_volume / raw_volume) ** (1.0 / 3.0)

        lattice = lattice_raw * scale

    elif len(scale_values) == 3:
        # Rare POSCAR form: separate scale factors for a, b, c
        scales = np.array(scale_values, dtype=float)
        lattice = lattice_raw * scales[:, None]

    else:
        raise ValueError(f"Cannot parse scale line in {poscar}: {lines[1]}")

    return lattice


def cell_volume(lattice: np.ndarray) -> float:
    """
    Cell volume in Å^3.
    """
    a, b, c = lattice
    return abs(np.dot(a, np.cross(b, c)))


def find_structure_file(image_dir: str | Path) -> Path:
    """
    Prefer CONTCAR, then POSCAR.
    """
    image_dir = Path(image_dir)

    for name in ["CONTCAR", "POSCAR"]:
        path = image_dir / name
        if path.exists() and path.stat().st_size > 0:
            return path

    raise FileNotFoundError(f"No CONTCAR or POSCAR found in {image_dir}")


def read_vasp_dipole_from_outcar(outcar: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read ionic and total electronic dipole moments from VASP OUTCAR.

    Returns
    -------
    p_ion : np.ndarray
        Ionic dipole moment, unit e Å.
    p_elc : np.ndarray
        Total electronic dipole moment, unit e Å.
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
            f"Cannot find Ionic dipole moment in {outcar}. "
            "Please check LCALCPOL = .TRUE."
        )

    if not elc_matches:
        raise RuntimeError(
            f"Cannot find Total electronic dipole moment in {outcar}. "
            "Please check LCALCPOL = .TRUE."
        )

    # Use the last occurrence in case OUTCAR contains several ionic/electronic steps.
    p_ion = np.array([float(x) for x in ion_matches[-1]], dtype=float)
    p_elc = np.array([float(x) for x in elc_matches[-1]], dtype=float)

    return p_ion, p_elc


def polarization_quantum(lattice: np.ndarray, volume: float) -> np.ndarray:
    """
    Calculate polarization quantum basis.

    Q_a = e a / Ω
    Q_b = e b / Ω
    Q_c = e c / Ω

    Returns
    -------
    Q : np.ndarray
        shape (3, 3).
        Rows are Qa, Qb, Qc.
        Columns are Cartesian x, y, z.
        Unit: μC/cm^2.
    """
    return lattice / volume * E_PER_A2_TO_UC_CM2


def read_single_image_polarization(image_dir: str | Path) -> ImagePolarization:
    """
    Read one image directory, such as:

        1-dipol/00
        1-dipol/01
        ...
        1-dipol/08

    The image directory should contain OUTCAR and CONTCAR/POSCAR.
    """
    image_dir = Path(image_dir)
    outcar = image_dir / "OUTCAR"

    if not outcar.exists():
        raise FileNotFoundError(f"OUTCAR not found: {outcar}")

    structure = find_structure_file(image_dir)

    lattice = read_lattice_from_poscar(structure)
    volume = cell_volume(lattice)

    p_ion, p_elc = read_vasp_dipole_from_outcar(outcar)
    p_tot = p_ion + p_elc

    P_raw = p_tot / volume * E_PER_A2_TO_UC_CM2
    Q = polarization_quantum(lattice, volume)

    try:
        image = int(image_dir.name)
    except ValueError:
        image = -1

    return ImagePolarization(
        image=image,
        image_name=image_dir.name,
        image_dir=image_dir,
        volume=volume,
        lattice=lattice,
        p_ion=p_ion,
        p_elc=p_elc,
        p_tot=p_tot,
        P_raw=P_raw,
        Q=Q,
    )


def unwrap_polarization_branch(
    images: Sequence[ImagePolarization],
    start_branch: Sequence[int] = (0, 0, 0),
    search_range: int = 2,
    metric: str = "norm",
    direction: Sequence[float] | None = None,
) -> list[ImagePolarization]:
    """
    Automatically unwrap polarization branches along a NEB path.

    For each image i, choose integer branch:

        branch_i = [na, nb, nc]

    so that:

        P_corr_i = P_raw_i + branch_i @ Q_i

    is as continuous as possible with respect to image i-1.

    Parameters
    ----------
    images
        List of ImagePolarization sorted by image number.
    start_branch
        Branch of the first image. Usually use (0, 0, 0).
    search_range
        Search around previous branch.
        If previous branch is [0,0,0] and search_range=2,
        then next branch candidates are from -2 to +2 in each lattice direction.
    metric
        "norm" or "direction".
        "norm": minimize full 3D vector jump.
        "direction": minimize jump projected along given direction.
    direction
        Required if metric="direction".
        For example [1,0,0], [0,1,0], [0,0,1], or your switching direction.

    Returns
    -------
    images
        Same list, with branch, P_corr, jump_norm filled.
    """
    if len(images) == 0:
        return []

    images = sorted(images, key=lambda x: x.image)

    start_branch = np.array(start_branch, dtype=int)

    images[0].branch = start_branch
    images[0].P_corr = images[0].P_raw + start_branch @ images[0].Q
    images[0].jump_norm = 0.0

    if metric not in {"norm", "direction"}:
        raise ValueError("metric should be 'norm' or 'direction'.")

    if metric == "direction":
        if direction is None:
            raise ValueError("direction must be provided when metric='direction'.")
        direction_vec = np.array(direction, dtype=float)
        norm = np.linalg.norm(direction_vec)
        if norm < 1e-12:
            raise ValueError("direction vector is zero.")
        direction_vec = direction_vec / norm
    else:
        direction_vec = None

    local_shifts = list(
        itertools.product(
            range(-search_range, search_range + 1),
            repeat=3,
        )
    )

    for i in range(1, len(images)):
        prev = images[i - 1]
        cur = images[i]

        assert prev.P_corr is not None
        assert prev.branch is not None

        best_score = None
        best_branch = None
        best_P_corr = None
        best_jump = None

        # Search around previous branch, not only around zero.
        # This allows cumulative branch changes along a long path.
        for dn in local_shifts:
            dn = np.array(dn, dtype=int)
            branch = prev.branch + dn

            P_corr = cur.P_raw + branch @ cur.Q
            jump_vec = P_corr - prev.P_corr
            jump_norm = float(np.linalg.norm(jump_vec))

            if metric == "norm":
                score = jump_norm
            else:
                # minimize projected jump along selected switching direction
                score = abs(float(np.dot(jump_vec, direction_vec)))

            if best_score is None or score < best_score:
                best_score = score
                best_branch = branch
                best_P_corr = P_corr
                best_jump = jump_norm

        cur.branch = best_branch
        cur.P_corr = best_P_corr
        cur.jump_norm = best_jump

    return images


def read_dipol_path(
    dipol_dir: str | Path,
    images: Iterable[int] = range(9),
    start_branch: Sequence[int] = (0, 0, 0),
    search_range: int = 2,
    metric: str = "norm",
    direction: Sequence[float] | None = None,
) -> list[ImagePolarization]:
    """
    Read and unwrap one 1-dipol path.

    Example
    -------
    dipol_dir = "xyz-000-011/1-dipol"
    """
    dipol_dir = Path(dipol_dir)

    image_results: list[ImagePolarization] = []

    for i in images:
        image_dir = dipol_dir / f"{i:02d}"
        if not image_dir.exists():
            raise FileNotFoundError(f"Image directory not found: {image_dir}")

        image_results.append(read_single_image_polarization(image_dir))

    image_results = unwrap_polarization_branch(
        image_results,
        start_branch=start_branch,
        search_range=search_range,
        metric=metric,
        direction=direction,
    )

    return image_results


def images_to_dataframe(
    image_results: Sequence[ImagePolarization],
    system: str,
    path_name: str,
    case_dir: str | Path,
) -> pd.DataFrame:
    """
    Convert corrected polarization results to a pandas DataFrame.
    """
    rows = []

    for r in image_results:
        if r.branch is None or r.P_corr is None:
            raise RuntimeError("Polarization branch is not unwrapped yet.")

        row = {
            "system": system,
            "path": path_name,
            "case_dir": str(case_dir),
            "image": r.image,
            "image_name": r.image_name,
            "volume_A3": r.volume,

            "branch_a": int(r.branch[0]),
            "branch_b": int(r.branch[1]),
            "branch_c": int(r.branch[2]),
            "jump_norm_uC_cm2": r.jump_norm,

            "p_ion_x_eA": r.p_ion[0],
            "p_ion_y_eA": r.p_ion[1],
            "p_ion_z_eA": r.p_ion[2],

            "p_elc_x_eA": r.p_elc[0],
            "p_elc_y_eA": r.p_elc[1],
            "p_elc_z_eA": r.p_elc[2],

            "p_tot_x_eA": r.p_tot[0],
            "p_tot_y_eA": r.p_tot[1],
            "p_tot_z_eA": r.p_tot[2],

            "P_raw_x_uC_cm2": r.P_raw[0],
            "P_raw_y_uC_cm2": r.P_raw[1],
            "P_raw_z_uC_cm2": r.P_raw[2],

            "P_corr_x_uC_cm2": r.P_corr[0],
            "P_corr_y_uC_cm2": r.P_corr[1],
            "P_corr_z_uC_cm2": r.P_corr[2],
        }

        # Save polarization quantum matrix Q.
        # Rows: Qa, Qb, Qc. Columns: x, y, z.
        labels = ["a", "b", "c"]
        comps = ["x", "y", "z"]
        for ii, lab in enumerate(labels):
            for jj, comp in enumerate(comps):
                row[f"Q_{lab}_{comp}_uC_cm2"] = r.Q[ii, jj]

        rows.append(row)

    return pd.DataFrame(rows)


def collect_cases_from_base_path(base_path: str | Path) -> list[tuple[str, Path]]:
    """
    For one HfO2 base path, collect:

        HfO2
        a-HZO
        b-HZO
        c-HZO

    Example
    -------
    base_path = xyz-000-011

    returns:
        [
            ("HfO2", xyz-000-011),
            ("a-HZO", xyz-000-011/a-xyz-000-011),
            ("b-HZO", xyz-000-011/b-xyz-000-011),
            ("c-HZO", xyz-000-011/c-xyz-000-011),
        ]
    """
    base_path = Path(base_path)
    base_name = base_path.name

    cases: list[tuple[str, Path]] = []

    # HfO2 itself
    if (base_path / "1-dipol").exists():
        cases.append(("HfO2", base_path))

    # HZO doped cases
    for dop in ["a", "b", "c"]:
        dop_dir = base_path / f"{dop}-{base_name}"
        if (dop_dir / "1-dipol").exists():
            cases.append((f"{dop}-HZO", dop_dir))

    return cases


def process_one_neb_path(
    base_path: str | Path,
    images: Iterable[int] = range(9),
    search_range: int = 2,
    metric: str = "norm",
    direction: Sequence[float] | None = None,
    output_csv: str | Path | None = None,
) -> pd.DataFrame:
    """
    Process one full NEB path directory.

    Example structure
    -----------------
    xyz-000-011/
    ├── 1-dipol/00 ... 08
    ├── a-xyz-000-011/1-dipol/00 ... 08
    ├── b-xyz-000-011/1-dipol/00 ... 08
    └── c-xyz-000-011/1-dipol/00 ... 08

    Returns
    -------
    df : pandas.DataFrame
        Contains raw polarization, corrected polarization, and branch.
    """
    base_path = Path(base_path)
    path_name = base_path.name

    all_dfs = []

    cases = collect_cases_from_base_path(base_path)

    if len(cases) == 0:
        raise RuntimeError(f"No valid 1-dipol cases found under {base_path}")

    for system, case_dir in cases:
        dipol_dir = case_dir / "1-dipol"

        print(f"[INFO] Processing {system}: {dipol_dir}")

        image_results = read_dipol_path(
            dipol_dir=dipol_dir,
            images=images,
            start_branch=(0, 0, 0),
            search_range=search_range,
            metric=metric,
            direction=direction,
        )

        df = images_to_dataframe(
            image_results=image_results,
            system=system,
            path_name=path_name,
            case_dir=case_dir,
        )

        all_dfs.append(df)

    result = pd.concat(all_dfs, ignore_index=True)

    if output_csv is not None:
        output_csv = Path(output_csv)
        result.to_csv(output_csv, index=False)
        print(f"[INFO] Wrote {output_csv}")

    return result



# ============================================================
# NEB energy and electric-field-corrected barrier utilities
# ============================================================

FIELD_COEFF = 6.241509074e-6
# eV / (Å^3 * μC cm^-2 * MV cm^-1)


def read_vasp_energy_from_outcar(
    outcar: str | Path,
    energy_type: str = "sigma0",
) -> float:
    """
    Read final energy from VASP OUTCAR.

    Parameters
    ----------
    outcar
        Path to OUTCAR.
    energy_type
        "sigma0":
            use energy(sigma->0), recommended for comparing static energies.
        "without_entropy":
            use energy without entropy.
        "toten":
            use free energy TOTEN.

    Returns
    -------
    energy : float
        Energy in eV.
    """
    outcar = Path(outcar)
    text = outcar.read_text(errors="ignore")

    # Example:
    # energy  without entropy=    -123.456  energy(sigma->0) =    -123.457
    sigma_pattern = re.compile(
        r"energy\s+without\s+entropy=\s*([-+0-9.Ee]+)"
        r"\s+energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)"
    )

    # Example:
    # free  energy   TOTEN  =      -123.456 eV
    toten_pattern = re.compile(
        r"free\s+energy\s+TOTEN\s+=\s*([-+0-9.Ee]+)\s+eV"
    )

    sigma_matches = sigma_pattern.findall(text)
    toten_matches = toten_pattern.findall(text)

    energy_type = energy_type.lower()

    if energy_type == "sigma0":
        if sigma_matches:
            return float(sigma_matches[-1][1])
        if toten_matches:
            return float(toten_matches[-1])
        raise RuntimeError(f"Cannot find sigma0 or TOTEN energy in {outcar}")

    elif energy_type == "without_entropy":
        if sigma_matches:
            return float(sigma_matches[-1][0])
        if toten_matches:
            return float(toten_matches[-1])
        raise RuntimeError(f"Cannot find without_entropy or TOTEN energy in {outcar}")

    elif energy_type == "toten":
        if toten_matches:
            return float(toten_matches[-1])
        raise RuntimeError(f"Cannot find TOTEN energy in {outcar}")

    else:
        raise ValueError(
            "energy_type should be one of: 'sigma0', 'without_entropy', 'toten'"
        )


def read_vasp_energy_from_oszicar(
    oszicar: str | Path,
    energy_type: str = "sigma0",
) -> float:
    """
    Fallback reader for OSZICAR.

    It tries to read E0 first, then F.

    Example OSZICAR final line:
    1 F= -.123456E+03 E0= -.123457E+03 d E = ...
    """
    oszicar = Path(oszicar)
    text = oszicar.read_text(errors="ignore").splitlines()

    e0_pattern = re.compile(r"E0=\s*([-+0-9.Ee]+)")
    f_pattern = re.compile(r"F=\s*([-+0-9.Ee]+)")

    e0_values = []
    f_values = []

    for line in text:
        e0_match = e0_pattern.search(line)
        f_match = f_pattern.search(line)

        if e0_match:
            e0_values.append(float(e0_match.group(1)))

        if f_match:
            f_values.append(float(f_match.group(1)))

    energy_type = energy_type.lower()

    if energy_type in {"sigma0", "without_entropy"}:
        if e0_values:
            return e0_values[-1]
        if f_values:
            return f_values[-1]

    elif energy_type == "toten":
        if f_values:
            return f_values[-1]
        if e0_values:
            return e0_values[-1]

    raise RuntimeError(f"Cannot find energy in OSZICAR: {oszicar}")


def read_single_neb_image_energy(
    image_dir: str | Path,
    energy_type: str = "sigma0",
) -> float:
    """
    Read one NEB image energy.

    It first tries OUTCAR, then OSZICAR.

    Parameters
    ----------
    image_dir
        Directory such as 00, 01, ..., 08.
    energy_type
        "sigma0", "without_entropy", or "toten".

    Returns
    -------
    energy : float
        Energy in eV.
    """
    image_dir = Path(image_dir)

    outcar = image_dir / "OUTCAR"
    oszicar = image_dir / "OSZICAR"

    if outcar.exists() and outcar.stat().st_size > 0:
        return read_vasp_energy_from_outcar(outcar, energy_type=energy_type)

    if oszicar.exists() and oszicar.stat().st_size > 0:
        return read_vasp_energy_from_oszicar(oszicar, energy_type=energy_type)

    raise FileNotFoundError(f"No OUTCAR or OSZICAR found in {image_dir}")


def read_neb_energy_path(
    case_dir: str | Path,
    images: Iterable[int] = range(9),
    energy_type: str = "sigma0",
) -> pd.DataFrame:
    """
    Read NEB image energies from one case directory.

    Example
    -------
    case_dir = xyz-000-011
    or
    case_dir = xyz-000-011/a-xyz-000-011

    It reads:
        case_dir/00
        case_dir/01
        ...
        case_dir/08

    Returns
    -------
    df_energy : pandas.DataFrame
    """
    case_dir = Path(case_dir)

    rows = []

    for i in images:
        image_dir = case_dir / f"{i:02d}"

        if not image_dir.exists():
            raise FileNotFoundError(f"NEB image directory not found: {image_dir}")

        energy = read_single_neb_image_energy(
            image_dir=image_dir,
            energy_type=energy_type,
        )

        rows.append(
            {
                "case_dir": str(case_dir),
                "image": i,
                "neb_image_dir": str(image_dir),
                "E_neb_eV": energy,
            }
        )

    df = pd.DataFrame(rows)

    e0 = float(df.loc[df["image"] == min(images), "E_neb_eV"].iloc[0])
    df["E_neb_rel_initial_eV"] = df["E_neb_eV"] - e0
    df["E_neb_rel_min_eV"] = df["E_neb_eV"] - df["E_neb_eV"].min()

    return df


def read_all_neb_energies_for_base_path(
    base_path: str | Path,
    images: Iterable[int] = range(9),
    energy_type: str = "sigma0",
) -> pd.DataFrame:
    """
    Read NEB energies for HfO2, a-HZO, b-HZO, c-HZO
    under one base path.

    Uses collect_cases_from_base_path from the previous code.
    """
    base_path = Path(base_path)
    path_name = base_path.name

    all_dfs = []

    cases = collect_cases_from_base_path(base_path)

    if len(cases) == 0:
        raise RuntimeError(f"No valid NEB cases found under {base_path}")

    for system, case_dir in cases:
        print(f"[INFO] Reading NEB energy for {system}: {case_dir}")

        df_energy = read_neb_energy_path(
            case_dir=case_dir,
            images=images,
            energy_type=energy_type,
        )

        df_energy.insert(0, "system", system)
        df_energy.insert(1, "path", path_name)

        all_dfs.append(df_energy)

    return pd.concat(all_dfs, ignore_index=True)


def merge_energy_and_polarization(
    df_pol: pd.DataFrame,
    df_energy: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge corrected polarization data with NEB energy data.

    Required columns in df_pol:
        system, path, case_dir, image,
        volume_A3,
        P_corr_x_uC_cm2, P_corr_y_uC_cm2, P_corr_z_uC_cm2

    Required columns in df_energy:
        system, path, case_dir, image, E_neb_eV
    """
    merge_cols = ["system", "path", "case_dir", "image"]

    df = pd.merge(
        df_pol,
        df_energy,
        on=merge_cols,
        how="inner",
        validate="one_to_one",
    )

    if len(df) == 0:
        raise RuntimeError(
            "Merged dataframe is empty. "
            "Please check system/path/case_dir/image columns."
        )

    return df


def make_field_vector(
    field_strength_MV_cm: float,
    direction: Sequence[float],
) -> np.ndarray:
    """
    Construct electric field vector in MV/cm.

    Parameters
    ----------
    field_strength_MV_cm
        Scalar field strength, unit MV/cm.
    direction
        Direction vector, for example:
            [1, 0, 0]
            [0, 1, 0]
            [0, 0, 1]
            [1, 1, 0]

    Returns
    -------
    E_vec : np.ndarray
        Electric field vector, unit MV/cm.
    """
    direction = np.array(direction, dtype=float)
    norm = np.linalg.norm(direction)

    if norm < 1e-12:
        raise ValueError("Electric-field direction vector is zero.")

    direction = direction / norm

    return field_strength_MV_cm * direction


def apply_electric_field_correction(
    df: pd.DataFrame,
    field_vector_MV_cm: Sequence[float],
    energy_col: str = "E_neb_eV",
    volume_col: str = "volume_A3",
    px_col: str = "P_corr_x_uC_cm2",
    py_col: str = "P_corr_y_uC_cm2",
    pz_col: str = "P_corr_z_uC_cm2",
) -> pd.DataFrame:
    """
    Apply electric-field correction to each NEB image.

    Formula
    -------
    F_i(E) = E_i - C * Ω_i * P_i · E

    where:
        C = 6.241509074e-6
        Ω: Å^3
        P: μC/cm^2
        E: MV/cm
        F: eV

    Parameters
    ----------
    df
        Merged energy + polarization dataframe.
    field_vector_MV_cm
        Electric field vector [Ex, Ey, Ez], unit MV/cm.

    Returns
    -------
    df_field
        Dataframe with additional columns:
            Efield_x_MV_cm
            Efield_y_MV_cm
            Efield_z_MV_cm
            PdotE
            field_correction_eV
            F_field_eV
            F_field_rel_initial_eV
            F_field_rel_min_eV
    """
    df_field = df.copy()

    E_vec = np.array(field_vector_MV_cm, dtype=float)

    if E_vec.shape != (3,):
        raise ValueError("field_vector_MV_cm should have shape (3,)")

    P = df_field[[px_col, py_col, pz_col]].to_numpy(dtype=float)

    pdotE = P @ E_vec

    field_corr = -FIELD_COEFF * df_field[volume_col].to_numpy(dtype=float) * pdotE

    df_field["Efield_x_MV_cm"] = E_vec[0]
    df_field["Efield_y_MV_cm"] = E_vec[1]
    df_field["Efield_z_MV_cm"] = E_vec[2]

    df_field["PdotE_uC_cm2_MV_cm"] = pdotE
    df_field["field_correction_eV"] = field_corr
    df_field["F_field_eV"] = df_field[energy_col].to_numpy(dtype=float) + field_corr

    # Relative energies are computed separately for each system/path/case.
    group_cols = ["system", "path", "case_dir"]

    df_field["F_field_rel_initial_eV"] = np.nan
    df_field["F_field_rel_final_eV"] = np.nan
    df_field["F_field_rel_min_eV"] = np.nan

    for _, idx in df_field.groupby(group_cols).groups.items():
        sub = df_field.loc[idx].sort_values("image")

        first_energy = float(sub.iloc[0]["F_field_eV"])
        final_energy = float(sub.iloc[-1]["F_field_eV"])
        min_energy = float(sub["F_field_eV"].min())

        df_field.loc[sub.index, "F_field_rel_initial_eV"] = (
            sub["F_field_eV"] - first_energy
        )

        df_field.loc[sub.index, "F_field_rel_final_eV"] = (
            sub["F_field_eV"] - final_energy
        )

        df_field.loc[sub.index, "F_field_rel_min_eV"] = (
            sub["F_field_eV"] - min_energy
        )

    return df_field


def calculate_neb_barrier_from_dataframe(
    df: pd.DataFrame,
    energy_col: str = "E_neb_eV",
    group_cols: Sequence[str] = ("system", "path", "case_dir"),
) -> pd.DataFrame:
    """
    Calculate NEB barrier for each system/path/case.

    For each group:
        forward barrier = max(E_i) - E_initial
        backward barrier = max(E_i) - E_final
        reaction energy = E_final - E_initial

    Parameters
    ----------
    df
        Dataframe containing image energies.
    energy_col
        Column name for energy, such as:
            E_neb_eV
            F_field_eV
    group_cols
        Grouping columns.

    Returns
    -------
    df_barrier : pandas.DataFrame
    """
    rows = []

    for group_key, sub in df.groupby(list(group_cols)):
        sub = sub.sort_values("image").copy()

        e_initial = float(sub.iloc[0][energy_col])
        e_final = float(sub.iloc[-1][energy_col])

        imax = sub[energy_col].idxmax()
        saddle_image = int(sub.loc[imax, "image"])
        e_saddle = float(sub.loc[imax, energy_col])

        forward_barrier = e_saddle - e_initial
        backward_barrier = e_saddle - e_final
        reaction_energy = e_final - e_initial

        if not isinstance(group_key, tuple):
            group_key = (group_key,)

        row = {
            col: val for col, val in zip(group_cols, group_key)
        }

        row.update(
            {
                "energy_col": energy_col,
                "initial_image": int(sub.iloc[0]["image"]),
                "final_image": int(sub.iloc[-1]["image"]),
                "saddle_image": saddle_image,
                "E_initial_eV": e_initial,
                "E_final_eV": e_final,
                "E_saddle_eV": e_saddle,
                "barrier_forward_eV": forward_barrier,
                "barrier_backward_eV": backward_barrier,
                "reaction_energy_eV": reaction_energy,
            }
        )

        rows.append(row)

    return pd.DataFrame(rows)


def process_one_neb_path_with_energy_and_polarization(
    base_path: str | Path,
    images: Iterable[int] = range(9),
    search_range: int = 2,
    unwrap_metric: str = "norm",
    unwrap_direction: Sequence[float] | None = None,
    energy_type: str = "sigma0",
    output_prefix: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full zero-field workflow:

    1. Read Berry-phase dipoles from 1-dipol/00...08
    2. Correct polarization branch
    3. Read NEB energies from 00...08
    4. Merge energy and polarization
    5. Calculate zero-field NEB barriers

    Returns
    -------
    df_merged
        Image-level dataframe.
    df_barrier
        Barrier summary dataframe.
    """
    base_path = Path(base_path)

    df_pol = process_one_neb_path(
        base_path=base_path,
        images=images,
        search_range=search_range,
        metric=unwrap_metric,
        direction=unwrap_direction,
        output_csv=None,
    )

    df_energy = read_all_neb_energies_for_base_path(
        base_path=base_path,
        images=images,
        energy_type=energy_type,
    )

    df_merged = merge_energy_and_polarization(
        df_pol=df_pol,
        df_energy=df_energy,
    )

    df_barrier = calculate_neb_barrier_from_dataframe(
        df=df_merged,
        energy_col="E_neb_eV",
    )

    if output_prefix is not None:
        output_prefix = Path(output_prefix)

        image_csv = output_prefix.with_name(output_prefix.name + "_images.csv")
        barrier_csv = output_prefix.with_name(output_prefix.name + "_barriers.csv")

        df_merged.to_csv(image_csv, index=False)
        df_barrier.to_csv(barrier_csv, index=False)

        print(f"[INFO] Wrote {image_csv}")
        print(f"[INFO] Wrote {barrier_csv}")

    return df_merged, df_barrier


def process_one_neb_path_under_field(
    base_path: str | Path,
    field_strength_MV_cm: float,
    field_direction: Sequence[float],
    images: Iterable[int] = range(9),
    search_range: int = 2,
    unwrap_metric: str = "norm",
    unwrap_direction: Sequence[float] | None = None,
    energy_type: str = "sigma0",
    output_prefix: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full electric-field-corrected workflow for one field value.

    Formula
    -------
    F_i(E) = E_i - C * Ω_i * P_i^corr · E

    Parameters
    ----------
    field_strength_MV_cm
        Field strength, unit MV/cm.
        Can be positive or negative.
    field_direction
        Field direction vector.

    Returns
    -------
    df_field
        Image-level dataframe with field correction.
    df_barrier_field
        Barrier summary under electric field.
    """
    df_merged, _ = process_one_neb_path_with_energy_and_polarization(
        base_path=base_path,
        images=images,
        search_range=search_range,
        unwrap_metric=unwrap_metric,
        unwrap_direction=unwrap_direction,
        energy_type=energy_type,
        output_prefix=None,
    )

    E_vec = make_field_vector(
        field_strength_MV_cm=field_strength_MV_cm,
        direction=field_direction,
    )

    df_field = apply_electric_field_correction(
        df=df_merged,
        field_vector_MV_cm=E_vec,
    )

    df_barrier_field = calculate_neb_barrier_from_dataframe(
        df=df_field,
        energy_col="F_field_eV",
    )

    df_barrier_field["field_strength_MV_cm"] = field_strength_MV_cm
    df_barrier_field["field_dir_x"] = E_vec[0]
    df_barrier_field["field_dir_y"] = E_vec[1]
    df_barrier_field["field_dir_z"] = E_vec[2]

    if output_prefix is not None:
        output_prefix = Path(output_prefix)

        image_csv = output_prefix.with_name(output_prefix.name + "_field_images.csv")
        barrier_csv = output_prefix.with_name(output_prefix.name + "_field_barriers.csv")

        df_field.to_csv(image_csv, index=False)
        df_barrier_field.to_csv(barrier_csv, index=False)

        print(f"[INFO] Wrote {image_csv}")
        print(f"[INFO] Wrote {barrier_csv}")

    return df_field, df_barrier_field


def scan_field_dependent_barriers(
    base_path: str | Path,
    field_values_MV_cm: Sequence[float],
    field_direction: Sequence[float],
    images: Iterable[int] = range(9),
    search_range: int = 2,
    unwrap_metric: str = "norm",
    unwrap_direction: Sequence[float] | None = None,
    energy_type: str = "sigma0",
    output_prefix: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Scan NEB barriers under a series of electric fields.

    This is useful for estimating coercive field by extrapolating
    barrier_forward_eV -> 0.

    Parameters
    ----------
    field_values_MV_cm
        List or array of field values, unit MV/cm.
        Example:
            np.linspace(-10, 10, 41)
    field_direction
        Electric field direction vector.

    Returns
    -------
    df_all_images
        Image-level dataframe for all fields.
    df_all_barriers
        Barrier dataframe for all fields.
    """
    base_path = Path(base_path)

    # Read zero-field energy and corrected polarization only once.
    df_merged, _ = process_one_neb_path_with_energy_and_polarization(
        base_path=base_path,
        images=images,
        search_range=search_range,
        unwrap_metric=unwrap_metric,
        unwrap_direction=unwrap_direction,
        energy_type=energy_type,
        output_prefix=None,
    )

    image_dfs = []
    barrier_dfs = []

    for field in field_values_MV_cm:
        E_vec = make_field_vector(
            field_strength_MV_cm=float(field),
            direction=field_direction,
        )

        df_field = apply_electric_field_correction(
            df=df_merged,
            field_vector_MV_cm=E_vec,
        )

        df_field["field_strength_MV_cm"] = float(field)

        df_barrier = calculate_neb_barrier_from_dataframe(
            df=df_field,
            energy_col="F_field_eV",
        )

        df_barrier["field_strength_MV_cm"] = float(field)
        df_barrier["Efield_x_MV_cm"] = E_vec[0]
        df_barrier["Efield_y_MV_cm"] = E_vec[1]
        df_barrier["Efield_z_MV_cm"] = E_vec[2]

        image_dfs.append(df_field)
        barrier_dfs.append(df_barrier)

    df_all_images = pd.concat(image_dfs, ignore_index=True)
    df_all_barriers = pd.concat(barrier_dfs, ignore_index=True)

    if output_prefix is not None:
        output_prefix = Path(output_prefix)

        image_csv = output_prefix.with_name(output_prefix.name + "_field_scan_images.csv")
        barrier_csv = output_prefix.with_name(output_prefix.name + "_field_scan_barriers.csv")

        df_all_images.to_csv(image_csv, index=False)
        df_all_barriers.to_csv(barrier_csv, index=False)

        print(f"[INFO] Wrote {image_csv}")
        print(f"[INFO] Wrote {barrier_csv}")

    return df_all_images, df_all_barriers






















# df = process_one_neb_path(
#     base_path="../Pca2_1/xyz-000-011",
#     images=range(9),
#     search_range=2,
#     output_csv="xyz-000-011_polarization_branch.csv",
# )

# print(df[[
#     "system",
#     "image",
#     "branch_a",
#     "branch_b",
#     "branch_c",
#     "P_raw_x_uC_cm2",
#     "P_raw_y_uC_cm2",
#     "P_raw_z_uC_cm2",
#     "P_corr_x_uC_cm2",
#     "P_corr_y_uC_cm2",
#     "P_corr_z_uC_cm2",
# ]])


# import numpy as np
# from polarization_neb import scan_field_dependent_barriers

field_values = np.linspace(0, 20, 20)

df_all_images, df_all_barriers = scan_field_dependent_barriers(
    base_path="../Pca2_1/xyz-000-011",
    field_values_MV_cm=field_values,
    field_direction=[1, 0, 0],
    images=range(9),
    search_range=3,
    energy_type="sigma0",
    output_prefix="xyz-000-011_field_scan_x",
)

print(df_all_barriers[[
    "system",
    "field_strength_MV_cm",
    "saddle_image",
    "barrier_forward_eV",
    "barrier_backward_eV",
    "reaction_energy_eV",
]])