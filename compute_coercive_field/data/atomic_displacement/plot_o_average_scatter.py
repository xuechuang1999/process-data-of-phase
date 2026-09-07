#!/usr/bin/env python3
"""Plot average O displacement and path-length scatter charts from summary CSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator


CONFIG_ORDER = ["undoped", "a", "b", "c"]
CONFIG_CSV_LABELS = ["HfO2", "a-HZO", "b-HZO", "c-HZO"]
CONFIG_STYLES = {
    "undoped": {"marker": "s", "size": 30, "color": "dimgray", "label": r"HfO$_2$"},
    "a": {"marker": "o", "size": 30, "color": "#e74c3c", "label": "a-HZO"},
    "b": {"marker": "^", "size": 35, "color": "#2e86de", "label": "b-HZO"},
    "c": {"marker": "v", "size": 35, "color": "#f39c12", "label": "c-HZO"},
}
METRICS = {
    "avg_O_net_displacement_A": {
        "ylabel": r"Average O displacement ($\mathrm{\AA}$)",
        "data_title": "Average O displacement",
        "unit": "Å",
        "stem": "o_avg_net_displacement_scatter",
    },
    "avg_O_path_length_A": {
        "ylabel": r"Average O path length ($\mathrm{\AA}$)",
        "data_title": "Average O path length",
        "unit": "Å",
        "stem": "o_avg_path_length_scatter",
    },
}

# Publication x coordinates, including blank spaces between axis-permutation
# groups, matching the existing energy-barrier figures.
PATH_X = {
    "pca2_1": {
        "xyz-000-011": 1,
        "xyz-000-101": 2,
        "xyz-000-110": 3,
        "xyz-001-000": 4,
        "xyz-001-011": 5,
        "xyz-001-101": 6,
        "xyz-001-110": 7,
        "xzy-000-000": 9,
        "xzy-000-101": 10,
        "yxz-000-000": 12,
        "yxz-001-000": 13,
        "yzx-000-000": 15,
        "yzx-000-011": 16,
        "yzx-000-101": 17,
        "yzx-000-110": 18,
        "zxy-000-000": 20,
        "zxy-000-011": 21,
        "zxy-000-101": 22,
        "zxy-000-110": 23,
        "zyx-000-000": 25,
        "zyx-000-011": 26,
    },
    "pmn2_1": {
        "xyz-000-011": 1,
        "xyz-001-000": 2,
        "xyz-001-011": 3,
        "xyz-101-000": 4,
        "xyz-101-011": 5,
        "xzy-000-000": 7,
        "xzy-100-000": 8,
    },
}

# Keep phase labels in data-sparse corners for each phase. Pmn2_1 has high-value
# points at the right edge, so its label belongs in the upper-left corner.
PHASE_LABEL_LAYOUT = {
    "pmn2_1": (0.03, 0.95, "left"),
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
        "xtick.direction": "in",
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


def read_summary(csv_path: Path) -> pd.DataFrame:
    """Read and validate one average-displacement summary CSV."""

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
    """Return one row per NEB path and metric columns for the four systems."""

    group_order = summary[["system", "base_path"]].drop_duplicates()
    wide = summary.pivot(
        index=["system", "base_path"],
        columns="dopant_config",
        values=list(METRICS),
    ).reindex(CONFIG_ORDER, axis=1, level=1)
    wide = wide.reindex(pd.MultiIndex.from_frame(group_order))

    missing_points = wide.isna().any(axis=1)
    if missing_points.any():
        missing_groups = group_order.loc[missing_points.values]
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
    """Derive the publication phase label from the phase-specific CSV name."""

    suffix = csv_path.stem.removeprefix("o_path_average_summary")
    suffix = suffix.removeprefix("-").replace("_", "").lower()
    return {
        "pca21": r"$Pca2_{1}$",
        "pmn21": r"$Pmn2_{1}$",
    }.get(suffix, system or csv_path.stem)


def output_suffix(csv_path: Path) -> str:
    """Keep the phase suffix when naming output figures."""

    suffix = csv_path.stem.removeprefix("o_path_average_summary")
    if not suffix:
        return ""
    return suffix if suffix.startswith("-") else f"-{suffix}"


def phase_key(csv_path: Path) -> str:
    """Return the phase identifier used by PATH_X."""

    suffix = csv_path.stem.removeprefix("o_path_average_summary")
    return suffix.removeprefix("-").lower()
    return suffix


def write_plot_data(
    wide: pd.DataFrame,
    metric: str,
    path_names: list[str],
    x: np.ndarray,
    csv_path: Path,
    output_dir: Path,
) -> Path:
    """Write the exact wide-format values used by one scatter chart."""

    plot_data = pd.DataFrame(
        {
            "path": path_names,
            "x": x.astype(int),
            "HfO2": wide[(metric, "undoped")].to_numpy(dtype=float),
            "a-HZO": wide[(metric, "a")].to_numpy(dtype=float),
            "b-HZO": wide[(metric, "b")].to_numpy(dtype=float),
            "c-HZO": wide[(metric, "c")].to_numpy(dtype=float),
        }
    )
    data_stem = METRICS[metric]["stem"].removesuffix("_scatter")
    output_path = output_dir / f"{data_stem}_plot_data{output_suffix(csv_path)}.csv"
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["", "", *([METRICS[metric]["data_title"]] * len(CONFIG_CSV_LABELS))]
        )
        writer.writerow(
            ["", "", *([METRICS[metric]["unit"]] * len(CONFIG_CSV_LABELS))]
        )
        writer.writerow(["", "", *CONFIG_CSV_LABELS])
        for row in plot_data.to_numpy(dtype=object):
            writer.writerow(
                [row[0], int(row[1]), *(f"{value:.10g}" for value in row[2:])]
            )
    return output_path


def path_groups(path_names: list[str]) -> list[tuple[str, list[int]]]:
    """Group paths by their axis permutation, preserving CSV row order."""

    groups: list[tuple[str, list[int]]] = []
    for index, path_name in enumerate(path_names):
        group_name = path_name.split("-", 1)[0]
        if not groups or groups[-1][0] != group_name:
            groups.append((group_name, []))
        groups[-1][1].append(index)
    return groups


def plot_one_phase(
    wide: pd.DataFrame,
    phase: str,
    csv_path: Path,
    output_dir: Path,
    show: bool,
) -> list[Path]:
    """Draw one displacement and one path-length scatter chart for a phase."""

    path_names = [base for _, base in wide.index]
    mapping = PATH_X.get(phase_key(csv_path))
    if mapping is None:
        raise ValueError(
            f"No publication x-coordinate mapping is defined for {csv_path.name}. "
            "Supported phase suffixes: pca2_1 and pmn2_1."
        )
    expected_paths = list(mapping)
    if path_names != expected_paths:
        raise ValueError(
            f"The path order in {csv_path.name} does not match the publication "
            f"layout for {phase_key(csv_path)}.\n"
            f"Expected: {', '.join(expected_paths)}\n"
            f"Found:    {', '.join(path_names)}"
        )
    try:
        x = np.array([mapping[path_name] for path_name in path_names], dtype=float)
    except KeyError as error:
        raise ValueError(
            f"Path {error.args[0]!r} from {csv_path.name} is not in the "
            f"{phase_key(csv_path)} publication x-coordinate mapping."
        ) from error

    groups = path_groups(path_names)

    boundaries = [
        (x[indices[-1]] + x[groups[group_index + 1][1][0]]) / 2.0
        for group_index, (_, indices) in enumerate(groups[:-1])
    ]

    outputs: list[Path] = []

    for metric, metric_info in METRICS.items():
        fig, ax = plt.subplots(figsize=(4.2, 2.8))

        for config in CONFIG_ORDER:
            style = CONFIG_STYLES[config]
            ax.scatter(
                x,
                wide[(metric, config)].to_numpy(dtype=float),
                marker=style["marker"],
                s=style["size"],
                color=style["color"],
                label=style["label"],
            )

        for group_name, indices in groups:
            center = float(np.mean(x[indices]))
            axis_label = ",".join(list(group_name))
            ax.text(
                center,
                0.02,
                r"$\{" + axis_label + r"\}$",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontweight="normal",
                fontsize=8.5,
            )

        ax.set_ylabel(metric_info["ylabel"], fontsize=12)
        values = wide[metric].to_numpy(dtype=float).ravel()
        y_min = float(np.nanmin(values))
        y_max = float(np.nanmax(values))
        span = y_max - y_min
        pad = max(span * 0.10, abs(y_max) * 0.02, 1e-6)
        ax.set_ylim(y_min - pad, y_max + pad)
        ax.set_xlim(x.min() - 1.0, x.max() + 1.0)
        ax.set_xticks(boundaries)
        ax.set_xticklabels([])
        ax.tick_params(axis="y", which="both", left=True, right=False)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, 1.10),
            ncol=4,
            fontsize=10,
            handlelength=1.25,
            handleheight=1,
            borderaxespad=0.4,
            columnspacing=0.7,
            labelspacing=0.1,
            handletextpad=0.2,
            borderpad=0.1,
        )

        label_x, label_y, label_ha = PHASE_LABEL_LAYOUT.get(
            phase_key(csv_path), (0.97, 0.95, "right")
        )
        ax.text(
            label_x,
            label_y,
            phase,
            transform=ax.transAxes,
            ha=label_ha,
            va="top",
            fontsize=12,
            fontstyle="italic",
            fontweight="normal",
        )

        plt.tight_layout()
        stem = f"{metric_info['stem']}{output_suffix(csv_path)}"
        for extension, transparent in ((".png", True), (".pdf", True), (".svg", True)):
            output_path = output_dir / f"{stem}{extension}"
            fig.savefig(output_path, transparent=transparent)
            outputs.append(output_path)
        outputs.append(
            write_plot_data(wide, metric, path_names, x, csv_path, output_dir)
        )

        if show:
            plt.show()
        plt.close(fig)

    return outputs


def discover_summary_files(input_dir: Path) -> list[Path]:
    """Find the two phase summary CSVs in the preferred plotting order."""

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
        help="Directory for PNG/PDF/SVG files (default: each CSV's directory)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display each figure after saving it",
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
        wide = make_wide_summary(summary)
        phase = phase_label(str(summary["system"].iloc[0]), csv_path)
        charts = plot_one_phase(wide, phase, csv_path, output_dir, args.show)

        print(f"[OK] Read summary CSV: {csv_path} ({len(summary)} rows)")
        print(f"[OK] Scatter charts:  {', '.join(path.name for path in charts)}")


if __name__ == "__main__":
    main()
