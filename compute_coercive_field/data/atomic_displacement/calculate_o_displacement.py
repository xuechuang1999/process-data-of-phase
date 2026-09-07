#!/usr/bin/env python3
"""Calculate average O-atom displacement and path length for VASP NEB paths."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


BASE_PATH_RE = re.compile(r"^(?P<direction>[a-z]+)-(?P<start>[01]{3})-(?P<end>[01]{3})$")


@dataclass(frozen=True)
class AtomicMetrics:
    atom_index: int
    element: str
    net_displacement: float
    path_length: float
    max_single_step: float


@dataclass(frozen=True)
class Poscar:
    comment: str
    lattice: np.ndarray
    symbols: list[str]
    counts: list[int]
    frac: np.ndarray

    @property
    def natoms(self) -> int:
        return sum(self.counts)


@dataclass(frozen=True)
class PathMetrics:
    system: str
    base_path: str
    dopant_config: str
    folder_name: str
    relative_path: str
    image_sequence: str
    nimages: int
    metrics: list[AtomicMetrics]


def _is_int_list(words: Sequence[str]) -> bool:
    try:
        [int(word) for word in words]
    except (TypeError, ValueError):
        return False
    return True


def _effective_scale(raw_lattice: np.ndarray, scale: float, path: Path) -> float:
    if scale == 0.0:
        raise ValueError(f"{path}: POSCAR scale factor cannot be zero")
    if scale > 0.0:
        return scale

    raw_volume = abs(float(np.linalg.det(raw_lattice)))
    if raw_volume < 1.0e-14:
        raise ValueError(f"{path}: POSCAR lattice is singular")
    return (-scale / raw_volume) ** (1.0 / 3.0)


def read_poscar(path: Path) -> Poscar:
    """Read a VASP POSCAR/CONTCAR and normalize coordinates to Direct."""

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 8:
        raise ValueError(f"{path}: file is too short to be a POSCAR")

    scale = float(lines[1].split()[0])
    raw_lattice = np.array(
        [[float(value) for value in lines[index].split()[:3]] for index in range(2, 5)],
        dtype=float,
    )
    factor = _effective_scale(raw_lattice, scale, path)
    lattice = raw_lattice * factor

    cursor = 5
    symbols = lines[cursor].split()
    if _is_int_list(symbols):
        raise ValueError(f"{path}: element symbols are required")
    cursor += 1
    counts = [int(value) for value in lines[cursor].split()]
    cursor += 1
    if len(symbols) != len(counts):
        raise ValueError(f"{path}: numbers of symbols and counts differ")

    if lines[cursor].strip().lower().startswith("s"):
        cursor += 1

    coordinate_mode = lines[cursor].strip().lower()
    cursor += 1
    if coordinate_mode.startswith("d"):
        direct = True
    elif coordinate_mode.startswith(("c", "k")):
        direct = False
    else:
        raise ValueError(f"{path}: unknown coordinate mode {coordinate_mode!r}")

    natoms = sum(counts)
    if len(lines) < cursor + natoms:
        raise ValueError(f"{path}: expected {natoms} coordinate lines")

    coordinates = np.array(
        [
            [float(value) for value in lines[cursor + index].split()[:3]]
            for index in range(natoms)
        ],
        dtype=float,
    )
    if direct:
        frac = coordinates
    else:
        frac = (coordinates * factor) @ np.linalg.inv(lattice)

    return Poscar(
        comment=lines[0].rstrip(),
        lattice=lattice,
        symbols=symbols,
        counts=counts,
        frac=frac % 1.0,
    )


def image_dirs(path: Path) -> list[Path]:
    """Return sorted numeric NEB image directories, ignoring 1-dipol and similar."""

    candidates: list[tuple[int, Path]] = []
    for child in path.iterdir():
        if child.is_dir() and child.name.isdigit():
            candidates.append((int(child.name), child))
    return [item[1] for item in sorted(candidates)]


def structure_file(image_dir: Path) -> Path:
    """Return the relaxed structure file for one NEB image."""

    contcar = image_dir / "CONTCAR"
    if contcar.is_file():
        return contcar
    poscar = image_dir / "POSCAR"
    if poscar.is_file():
        return poscar
    raise FileNotFoundError(f"missing CONTCAR and POSCAR in {image_dir}")


def read_path_structures(path: Path) -> tuple[list[str], list[Poscar]]:
    images = image_dirs(path)
    if len(images) < 2:
        raise ValueError(f"{path}: fewer than two numeric NEB image directories")

    structures = [read_poscar(structure_file(image)) for image in images]
    first_symbols = structures[0].symbols
    if first_symbols is None:
        raise ValueError(f"{structure_file(images[0])}: POSCAR has no element symbols")

    expected_natoms = structures[0].natoms
    for image, structure in zip(images, structures):
        if structure.symbols is None:
            raise ValueError(f"{structure_file(image)}: POSCAR has no element symbols")
        if structure.natoms != expected_natoms:
            raise ValueError(
                f"{structure_file(image)}: atom count {structure.natoms} differs "
                f"from image 00 ({expected_natoms})"
            )

    elements = np.repeat(
        np.asarray(structures[0].symbols, dtype=object),
        np.asarray(structures[0].counts, dtype=int),
    ).tolist()
    return elements, structures


def minimum_image_cartesian(
    delta_frac: np.ndarray, lattice: np.ndarray
) -> np.ndarray:
    """Return the shortest Cartesian vector for nearby fractional coordinates.

    VASP CONTCAR positions are normally wrapped to [0, 1). Searching integer
    lattice translations around the nearest fractional vector also handles the
    slightly non-orthogonal lattices used by the Pmn2_1 phase.
    """

    nearest = np.rint(delta_frac)
    best_vector = None
    best_length = np.inf
    for shift_x in (-1, 0, 1):
        for shift_y in (-1, 0, 1):
            for shift_z in (-1, 0, 1):
                candidate_frac = delta_frac - (
                    nearest + np.array([shift_x, shift_y, shift_z], dtype=float)
                )
                candidate_cart = candidate_frac @ lattice
                length = float(np.linalg.norm(candidate_cart))
                if length < best_length:
                    best_length = length
                    best_vector = candidate_cart
    assert best_vector is not None
    return best_vector


def calculate_path(
    system: str,
    base_path_name: str,
    dopant_config: str,
    path: Path,
    root: Path,
) -> PathMetrics:
    elements, structures = read_path_structures(path)
    unwrapped_frac = [structures[0].frac.copy()]
    for previous, current in zip(structures[:-1], structures[1:]):
        mean_lattice = 0.5 * (previous.lattice + current.lattice)
        delta = current.frac - previous.frac
        delta_cart = np.apply_along_axis(
            minimum_image_cartesian, 1, delta, lattice=mean_lattice
        )
        step_frac = delta_cart @ np.linalg.inv(mean_lattice)
        unwrapped_frac.append(unwrapped_frac[-1] + step_frac)

    first_frac = np.asarray(unwrapped_frac[0], dtype=float)
    last_frac = np.asarray(unwrapped_frac[-1], dtype=float)
    net_lattice = 0.5 * (structures[0].lattice + structures[-1].lattice)
    net_cart = (last_frac - first_frac) @ net_lattice
    net_lengths = np.linalg.norm(net_cart, axis=1)

    step_lengths = np.zeros((len(structures) - 1, structures[0].natoms), dtype=float)
    for index, (previous, current) in enumerate(zip(structures[:-1], structures[1:])):
        mean_lattice = 0.5 * (previous.lattice + current.lattice)
        delta_frac = np.asarray(unwrapped_frac[index + 1]) - np.asarray(
            unwrapped_frac[index]
        )
        lengths = np.linalg.norm(delta_frac @ mean_lattice, axis=1)
        step_lengths[index] = lengths

    metrics: list[AtomicMetrics] = []
    for atom_index, element in enumerate(elements, start=1):
        if element.upper() != "O":
            continue
        metrics.append(
            AtomicMetrics(
                atom_index=atom_index,
                element=element,
                net_displacement=float(net_lengths[atom_index - 1]),
                path_length=float(step_lengths[:, atom_index - 1].sum()),
                max_single_step=float(step_lengths[:, atom_index - 1].max()),
            )
        )

    return PathMetrics(
        system=system,
        base_path=base_path_name,
        dopant_config=dopant_config,
        folder_name=path.name,
        relative_path=path.relative_to(root).as_posix(),
        image_sequence="|".join(image.name for image in image_dirs(path)),
        nimages=len(structures),
        metrics=metrics,
    )


def discover_cases(root: Path) -> tuple[list[Path], list[tuple[str, Path, Path]]]:
    """Find undoped base paths and their a/b/c doped partners."""

    bases: list[Path] = []
    for child in sorted(root.iterdir(), key=lambda item: item.name):
        if not child.is_dir() or BASE_PATH_RE.fullmatch(child.name) is None:
            continue
        if len(image_dirs(child)) >= 2:
            bases.append(child)

    pairs: list[tuple[str, Path, Path]] = []
    for base in bases:
        for dopant in ("a", "b", "c"):
            doped = base / f"{dopant}-{base.name}"
            if not doped.is_dir():
                continue
            if len(image_dirs(doped)) >= 2:
                pairs.append((dopant, base, doped))
    return bases, pairs


def atomic_rows(paths: list[tuple[str, str, str, str, Path, Path]]) -> list[dict]:
    rows: list[dict] = []
    for system, base_name, dopant, relative, path, root in paths:
        result = calculate_path(system, base_name, dopant, path, root)
        for metric in result.metrics:
            rows.append(
                {
                    "system": result.system,
                    "dopant_config": result.dopant_config,
                    "base_path": result.base_path,
                    "relative_path": result.relative_path,
                    "image_sequence": result.image_sequence,
                    "nimages": result.nimages,
                    "atom_index": metric.atom_index,
                    "element": metric.element,
                    "net_displacement_A": metric.net_displacement,
                    "path_length_A": metric.path_length,
                    "max_single_step_A": metric.max_single_step,
                }
            )
    return rows


def summary_frame(atomic: pd.DataFrame) -> pd.DataFrame:
    grouping = [
        "system",
        "dopant_config",
        "base_path",
        "relative_path",
        "image_sequence",
        "nimages",
    ]
    return (
        atomic.groupby(grouping, as_index=False, sort=False)
        .agg(
            n_O=("atom_index", "size"),
            avg_O_net_displacement_A=("net_displacement_A", "mean"),
            avg_O_path_length_A=("path_length_A", "mean"),
            std_O_net_displacement_A=("net_displacement_A", "std"),
            std_O_path_length_A=("path_length_A", "std"),
            max_O_net_displacement_A=("net_displacement_A", "max"),
            max_O_path_length_A=("path_length_A", "max"),
        )
        .sort_values(["system", "base_path", "dopant_config"], ignore_index=True)
    )


def comparison_frame(summary: pd.DataFrame) -> pd.DataFrame:
    undoped = summary[summary["dopant_config"] == "undoped"].copy()
    doped = summary[summary["dopant_config"].isin(("a", "b", "c"))].copy()

    rows: list[dict] = []
    for _, doped_row in doped.iterrows():
        matches = undoped[
            (undoped["system"] == doped_row["system"])
            & (undoped["base_path"] == doped_row["base_path"])
        ]
        if matches.empty:
            continue
        undoped_row = matches.iloc[0]
        row = {
            "system": doped_row["system"],
            "base_path": doped_row["base_path"],
            "dopant_config": doped_row["dopant_config"],
            "undoped_relative_path": undoped_row["relative_path"],
            "doped_relative_path": doped_row["relative_path"],
            "n_O_undoped": int(undoped_row["n_O"]),
            "n_O_doped": int(doped_row["n_O"]),
        }
        for metric in ("avg_O_net_displacement_A", "avg_O_path_length_A"):
            before = float(undoped_row[metric])
            after = float(doped_row[metric])
            prefix = metric.removeprefix("avg_O_").removesuffix("_A")
            row[f"undoped_{prefix}_A"] = before
            row[f"doped_{prefix}_A"] = after
            row[f"{prefix}_change_A"] = after - before
            row[f"{prefix}_change_percent"] = (
                (after - before) / before * 100.0 if before != 0.0 else np.nan
            )
        rows.append(row)

    columns = list(rows[0].keys()) if rows else []
    frame = pd.DataFrame(rows, columns=columns)
    if not frame.empty:
        frame = frame.sort_values(
            ["system", "base_path", "dopant_config"], ignore_index=True
        )
    return frame


def plot_summary(summary: pd.DataFrame, output_dir: Path) -> list[Path]:
    """Draw grouped bar charts for O net displacement and path length."""

    order = ["undoped", "a", "b", "c"]
    labels = ["HfO2", "a-HZO", "b-HZO", "c-HZO"]
    colors = ["#6E7B8B", "#D45B2D", "#2E7F73", "#7A5EA6"]

    grouped = (
        summary[summary["dopant_config"].isin(order)]
        .pivot_table(
            index=["system", "base_path"],
            columns="dopant_config",
            values=[
                "avg_O_net_displacement_A",
                "avg_O_path_length_A",
            ],
            aggfunc="first",
        )
        .reindex(order, axis=1, level=1)
    )

    group_names = [
        f"{system}\n{base}" for system, base in grouped.index
    ]
    x = np.arange(len(grouped), dtype=float)
    width = 0.19

    specs = [
        (
            "avg_O_net_displacement_A",
            "Average O net displacement",
            "o_net_displacement_bar.png",
        ),
        (
            "avg_O_path_length_A",
            "Average O path length",
            "o_path_length_bar.png",
        ),
    ]

    outputs: list[Path] = []
    for metric, title, filename in specs:
        fig, ax = plt.subplots(figsize=(max(7.0, 1.7 * len(grouped)), 5.0))
        for offset, dopant, label, color in zip(
            np.arange(4) - 1.5 * width, order, labels, colors
        ):
            values = grouped[(metric, dopant)].to_numpy(dtype=float)
            bars = ax.bar(
                x + offset,
                values,
                width,
                label=label,
                color=color,
                edgecolor="white",
                linewidth=0.7,
            )
            ax.bar_label(
                bars,
                fmt="%.3f",
                padding=2.5,
                fontsize=8.5,
                rotation=90,
            )

        ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
        ax.set_ylabel("Distance (Å)")
        ax.set_xticks(x, group_names)
        ax.tick_params(axis="x", labelsize=9)
        ax.grid(axis="y", color="#D8DDE3", linewidth=0.8, alpha=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        max_value = max(grouped[metric].max(axis=1).max(), 0.0)
        ax.set_ylim(0.0, max_value * 1.24)
        ax.legend(
            title="System",
            ncols=4,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.02),
        )
        fig.tight_layout()

        output = output_dir / filename
        fig.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(fig)
        outputs.append(output)
        print(f"[OK] Plot:             {output}")

    return outputs


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "roots",
        nargs="*",
        type=Path,
        default=[Path.cwd()],
        help="Directories containing xxx-nnn-nnn NEB paths (default: current directory)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir,
        help="Directory for output CSV files (default: data/atomic_displacement)",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Generate CSV files only",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cases: list[tuple[str, str, str, str, Path, Path]] = []
    for root in args.roots:
        system = root.resolve().name
        bases, doped_pairs = discover_cases(root)
        if not bases:
            print(f"[WARN] No complete base NEB path found under {root}")
            continue

        for base in bases:
            cases.append(
                    (system, base.name, "undoped", base.name, base, root)
            )
        for dopant, base, doped in doped_pairs:
            cases.append(
                (
                    system,
                    base.name,
                    dopant,
                    f"{dopant}-{base.name}",
                    doped,
                    root,
                )
            )

    if not cases:
        raise SystemExit("No complete NEB paths were found.")

    atomic = pd.DataFrame(atomic_rows(cases))
    summary = summary_frame(atomic)
    comparison = comparison_frame(summary)

    atomic_path = args.output_dir / "o_atomic_displacements.csv"
    summary_path = args.output_dir / "o_path_average_summary.csv"
    comparison_path = args.output_dir / "o_doping_comparison.csv"
    atomic.to_csv(atomic_path, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    comparison.to_csv(comparison_path, index=False, encoding="utf-8-sig")
    plots = [] if args.no_plot else plot_summary(summary, args.output_dir)

    print(f"[OK] Atomic O metrics: {atomic_path} ({len(atomic)} rows)")
    print(f"[OK] Path averages:    {summary_path} ({len(summary)} paths)")
    print(f"[OK] Doping compare:   {comparison_path} ({len(comparison)} pairs)")
    for plot in plots:
        print(f"[OK] Plot:             {plot}")


if __name__ == "__main__":
    main()
