#!/usr/bin/env python3
"""Plot per-O displacement and path length for one selected phase and NEB path."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator


CONFIG_ORDER = ["undoped", "a", "b", "c"]
CONFIG_LABELS = {
    "undoped": r"HfO$_{2}$",
    "a": r"a-HZO",
    "b": r"b-HZO",
    "c": r"c-HZO",
}
BAR_COLORS = {
    "undoped": "#9c9c9c",
    "a": "#ff0808",
    "b": "#1873b5",
    "c": "#ff8400",
}
METRICS = {
    "net_displacement_A": ("O net displacement", "o_atomic_net_displacement"),
    "path_length_A": ("O path length", "o_atomic_path_length"),
}
PHASES = {
    "pca2_1": r"$Pca2_{1}$",
    "pmn2_1": r"$Pmn2_{1}$",
}


rcParams.update(
    {
        "font.family": "Arial",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "axes.linewidth": 1.3,
        "lines.linewidth": 1.4,
        "lines.markersize": 5.0,
        "lines.markeredgewidth": 1.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def normalize_phase(phase: str) -> str:
    normalized = phase.strip().replace("_", "").replace("-", "").lower()
    matches = [name for name in PHASES if name.replace("_", "") == normalized]
    if not matches:
        choices = ", ".join(PHASES)
        raise ValueError(f"Unknown phase {phase!r}. Expected one of: {choices}")
    return matches[0]


def default_csv_path(phase: str) -> Path:
    return Path(__file__).resolve().parent / f"o_atomic_displacements-{phase}.csv"


def read_selected_atoms(csv_path: Path, phase: str, base_path: str) -> pd.DataFrame:
    if not csv_path.is_file():
        raise FileNotFoundError(f"Atomic-displacement CSV not found: {csv_path}")

    if f"-{phase}" not in csv_path.stem.lower():
        raise ValueError(
            f"CSV filename does not match phase {phase!r}: {csv_path.name}. "
            "Expected a phase-specific file such as "
            f"o_atomic_displacements-{phase}.csv"
        )

    atoms = pd.read_csv(csv_path, encoding="utf-8-sig")
    required = {
        "system",
        "dopant_config",
        "base_path",
        "atom_index",
        "element",
        *METRICS,
    }
    missing = sorted(required.difference(atoms.columns))
    if missing:
        raise ValueError(
            f"{csv_path} is missing required columns: {', '.join(missing)}"
        )

    selected = atoms.loc[atoms["base_path"] == base_path].copy()
    if selected.empty:
        available = sorted(atoms["base_path"].astype(str).unique())
        raise ValueError(
            f"Path {base_path!r} is not present in {csv_path}. "
            f"Available paths: {', '.join(available)}"
        )

    selected["dopant_config"] = selected["dopant_config"].astype(str).str.lower()
    selected = selected[selected["dopant_config"].isin(CONFIG_ORDER)].copy()
    present = set(selected["dopant_config"])
    absent = [config for config in CONFIG_ORDER if config not in present]
    if absent:
        raise ValueError(
            f"Path {base_path!r} in {csv_path} is missing configurations: "
            f"{', '.join(absent)}"
        )

    non_oxygen = selected.loc[selected["element"].astype(str).str.upper() != "O"]
    if not non_oxygen.empty:
        indices = sorted(non_oxygen["atom_index"].astype(int).unique().tolist())
        raise ValueError(f"Non-oxygen atoms are included at atom indices: {indices}")

    keys = ["dopant_config", "atom_index"]
    duplicates = selected.duplicated(keys, keep=False)
    if duplicates.any():
        examples = selected.loc[duplicates, keys].to_string(index=False)
        raise ValueError(
            f"{csv_path} has duplicate configuration/atom rows for {base_path}:\n"
            f"{examples}"
        )

    for metric in METRICS:
        selected[metric] = pd.to_numeric(selected[metric], errors="raise")

    selected["atom_index"] = selected["atom_index"].astype(int)
    return selected.sort_values(["atom_index", "dopant_config"]).reset_index(drop=True)


def plot_selected_atoms(
    atoms: pd.DataFrame,
    phase: str,
    base_path: str,
    output_dir: Path,
) -> list[Path]:
    atom_indices = sorted(atoms["atom_index"].unique().astype(int).tolist())
    x = np.arange(len(atom_indices), dtype=float)

    total_width = 0.72
    gap = 0.015
    bar_width = (total_width - gap * (len(CONFIG_ORDER) - 1)) / len(CONFIG_ORDER)
    offsets = (np.arange(len(CONFIG_ORDER)) - 1.5) * (bar_width + gap)

    outputs: list[Path] = []
    for metric, (metric_title, basename) in METRICS.items():
        fig, ax = plt.subplots(
            figsize=(max(8.5, 0.95 * len(atom_indices)), 3.2)
        )

        for offset, config in zip(offsets, CONFIG_ORDER):
            values = (
                atoms.loc[atoms["dopant_config"] == config]
                .set_index("atom_index")
                .reindex(atom_indices)[metric]
                .to_numpy(dtype=float)
            )
            bars = ax.bar(
                x + offset,
                values,
                bar_width,
                label=CONFIG_LABELS[config],
                color=BAR_COLORS[config],
                linewidth=0.5,
            )

            for bar in bars:
                height = float(bar.get_height())
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height * 1.012,
                    f"{height:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=5.5,
                    fontweight="bold",
                    rotation=90,
                )

        ax.set_title(
            f"{metric_title}, {PHASES[phase]}, {base_path}",
            fontsize=14,
            fontweight="bold",
            pad=12,
        )
        ax.set_ylabel(r"Distance ($\mathrm{\AA}$)", fontsize=12)
        ax.set_xlabel("O atom index", fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(atom_indices, fontsize=10, fontweight="bold")
        ax.tick_params(axis="y", labelsize=10)
        ax.tick_params(axis="x", width=2.0, length=4.0)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        largest_value = float(atoms[metric].max())
        ax.set_ylim(0.0, max(largest_value, 0.0) * 1.28)

        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        ax.legend(
            fontsize=10,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            frameon=False,
            ncols=4,
            borderaxespad=0.4,
            columnspacing=0.7,
            labelspacing=0.1,
            handletextpad=0.5,
            borderpad=0.1,
            handlelength=1.5,
        )

        fig.tight_layout()
        stem = f"{basename}-{phase}-{base_path}"
        png_path = output_dir / f"{stem}.png"
        pdf_path = output_dir / f"{stem}.pdf"
        fig.savefig(png_path, dpi=600)
        fig.savefig(pdf_path, transparent=True)
        plt.close(fig)
        outputs.extend([png_path, pdf_path])

    return outputs


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", help="Phase to plot: pca2_1 or pmn2_1")
    parser.add_argument("path", help="Exact NEB base path, for example xyz-000-011")
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help=f"Atomic-displacement CSV (default: o_atomic_displacements-<phase>.csv "
        f"beside this script)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for PNG/PDF files (default: the CSV's directory)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    phase = normalize_phase(args.phase)
    csv_path = args.csv.resolve() if args.csv is not None else default_csv_path(phase)
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else csv_path.parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = read_selected_atoms(csv_path, phase, args.path)
    charts = plot_selected_atoms(atoms, phase, args.path, output_dir)

    print(f"[OK] Read atomic CSV: {csv_path}")
    print(f"[OK] Selected rows:   {len(atoms)} O atoms x {len(CONFIG_ORDER)} systems")
    print(f"[OK] Charts:          {', '.join(path.name for path in charts)}")


if __name__ == "__main__":
    main()
