#!/usr/bin/env python3
"""Plot O displacement and path-length bar charts from existing summary CSV files."""

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
    "avg_O_net_displacement_A": (
        r"Average O net displacement",
        r"o_net_displacement_bar",
    ),
    "avg_O_path_length_A": (
        r"Average O path length",
        r"o_path_length_bar",
    ),
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
        "xtick.minor.size": 2.5,
        "ytick.major.size": 4,
        "ytick.minor.size": 2.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def read_summary(csv_path: Path) -> pd.DataFrame:
    """Read and validate one summary CSV produced by calculate_o_displacement.py."""

    if not csv_path.is_file():
        raise FileNotFoundError(f"Summary CSV not found: {csv_path}")

    summary = pd.read_csv(csv_path, encoding="utf-8-sig")
    required = {"system", "dopant_config", "base_path", *METRICS}
    missing = sorted(required.difference(summary.columns))
    if missing:
        raise ValueError(
            f"{csv_path} is missing required columns: {', '.join(missing)}"
        )

    summary["dopant_config"] = summary["dopant_config"].astype(str).str.lower()
    summary = summary[summary["dopant_config"].isin(CONFIG_ORDER)].copy()
    if summary.empty:
        raise ValueError(f"{csv_path} contains no undoped/a/b/c rows")

    keys = ["system", "base_path", "dopant_config"]
    duplicates = summary.duplicated(keys, keep=False)
    if duplicates.any():
        examples = summary.loc[duplicates, keys].to_string(index=False)
        raise ValueError(
            f"{csv_path} has duplicate path/configuration rows:\n{examples}"
        )

    for metric in METRICS:
        summary[metric] = pd.to_numeric(summary[metric], errors="raise")

    return summary


def make_wide_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Return one row per NEB path and one column per doping configuration."""

    group_order = summary[["system", "base_path"]].drop_duplicates()
    wide = summary.pivot(
        index=["system", "base_path"],
        columns="dopant_config",
        values=list(METRICS),
    ).reindex(CONFIG_ORDER, axis=1, level=1)
    wide = wide.reindex(pd.MultiIndex.from_frame(group_order))

    missing_bars = wide.isna().any(axis=1)
    if missing_bars.any():
        missing_groups = group_order.loc[missing_bars.values]
        details = "\n".join(
            f"  {system} / {base}"
            for system, base in missing_groups.itertuples(index=False)
        )
        raise ValueError(
            "These groups do not contain all four configurations "
            "(undoped, a, b, and c):\n" + details
        )

    return wide.astype(float)


def phase_label(system: str, csv_path: Path) -> str:
    """Use the CSV suffix as the phase label, normalized to publication notation."""

    suffix = csv_path.stem.removeprefix("o_path_average_summary")
    suffix = suffix.removeprefix("-").replace("_", "").lower()
    return {
        "pca21": r"$Pca2_{1}$",
        "pmn21": r"$Pmn2_{1}$",
    }.get(suffix, system or csv_path.stem)


def output_suffix(csv_path: Path) -> str:
    """Derive the filename suffix from a phase-specific summary CSV name."""

    suffix = csv_path.stem.removeprefix("o_path_average_summary")
    if not suffix:
        return ""
    return suffix if suffix.startswith("-") else f"-{suffix}"


def plot_one_summary(
    summary: pd.DataFrame,
    phase: str,
    csv_path: Path,
    output_dir: Path,
) -> list[Path]:
    """Draw separate grouped bar charts for displacement and path length."""

    wide = make_wide_summary(summary)
    path_names = [base for _, base in wide.index]
    x = np.arange(len(wide), dtype=float)

    # Keep bars within one cluster close together, while leaving clear space
    # between clusters for datasets with relatively few NEB paths.
    if len(wide) <= 8:
        total_width = 0.60
    elif len(wide) <= 12:
        total_width = 0.72
    else:
        total_width = 0.84
    gap = 0.01
    bar_width = (total_width - gap * (len(CONFIG_ORDER) - 1)) / len(CONFIG_ORDER)
    offsets = (np.arange(len(CONFIG_ORDER)) - 1.5) * (bar_width + gap)

    outputs: list[Path] = []
    for metric, (title, basename) in METRICS.items():
        fig, ax = plt.subplots(
            figsize=(max(9.5, 1.12 * len(wide)), 3.2)
        )

        bars_list = []
        for offset, config in zip(offsets, CONFIG_ORDER):
            values = wide[(metric, config)].to_numpy(dtype=float)
            bars = ax.bar(
                x + offset,
                values,
                bar_width,
                label=CONFIG_LABELS[config],
                color=BAR_COLORS[config],
                linewidth=0.5,
            )
            bars_list.append(bars)

            for bar in bars:
                height = float(bar.get_height())
                if height >= 0.0:
                    y_pos = height * 1.012
                    alignment, rotation = "bottom", 60
                else:
                    y_pos = height * 1.012
                    alignment, rotation = "top", -45
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0 + 0.01,
                    y_pos,
                    f"{height:.3f}",
                    ha="center",
                    va=alignment,
                    fontsize=5.5,
                    fontweight="bold",
                    rotation=rotation,
                )

        ax.set_title(f"{title}, {phase}", fontsize=14, fontweight="bold", pad=12)
        ax.set_ylabel(r"Distance ($\mathrm{\AA}$)", fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(path_names, rotation=0, ha="center", fontsize=9, fontweight="bold")
        ax.set_xlim(-0.62, len(wide) - 0.38)
        ax.tick_params(axis="y", labelsize=10)
        ax.tick_params(axis="x", width=2.0, length=4.0)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        largest_value = float(wide[metric].max(axis=1).max())
        ax.set_ylim(0.0, max(largest_value, 0.0) * 1.22)

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
        stem = f"{basename}{output_suffix(csv_path)}"
        png_path = output_dir / f"{stem}.png"
        pdf_path = output_dir / f"{stem}.pdf"
        fig.savefig(png_path, dpi=600)
        fig.savefig(pdf_path, transparent=True)
        plt.close(fig)
        outputs.extend([png_path, pdf_path])

    return outputs


def discover_summary_files(input_dir: Path) -> list[Path]:
    """Find phase-specific summary CSV files in the preferred plotting order."""

    files = {
        path.stem.removeprefix("o_path_average_summary-").lower(): path
        for path in input_dir.glob("o_path_average_summary-*.csv")
    }
    ordered = [files[name] for name in ("pca2_1", "pmn2_1") if name in files]
    ordered.extend(
        path for name, path in files.items() if name not in {"pca2_1", "pmn2_1"}
    )
    return ordered


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        nargs="*",
        type=Path,
        default=None,
        help="One or more o_path_average_summary CSV files "
        "(default: all files matching o_path_average_summary-*.csv beside this script)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for PNG/PDF files (default: each CSV's directory)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    default_dir = Path(__file__).resolve().parent
    csv_paths = args.csv or discover_summary_files(default_dir)
    if not csv_paths:
        raise SystemExit(
            f"No summary CSV files found. Expected o_path_average_summary-*.csv in {default_dir}"
        )

    for csv_path in csv_paths:
        csv_path = csv_path.resolve()
        output_dir = (
            args.output_dir.resolve()
            if args.output_dir is not None
            else csv_path.parent
        )
        output_dir.mkdir(parents=True, exist_ok=True)

        summary = read_summary(csv_path)
        phase = phase_label(str(summary["system"].iloc[0]), csv_path)
        charts = plot_one_summary(summary, phase, csv_path, output_dir)

        print(f"[OK] Read summary CSV: {csv_path} ({len(summary)} rows)")
        print(f"[OK] Bar charts:       {', '.join(path.name for path in charts)}")


if __name__ == "__main__":
    main()
