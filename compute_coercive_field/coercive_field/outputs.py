"""CSV and plot output helpers."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import to_rgb
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .calculation import (
    barrier_at_field,
    corrected_relative_energies,
    switching_axis,
    target_delta_dipoles,
    target_polarization,
)
from .constants import EV_PER_EANG_MV_CM
from .models import PathResult


PLOT_STYLE = {
    "axes.unicode_minus": False,
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.serif": ["Arial", "Helvetica"],
    "text.latex.preamble": r"\usepackage{helvet} \renewcommand{\familydefault}{\sfdefault}",
    "pgf.rcfonts": False,
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
    "xtick.major.size": 3,
    "xtick.minor.size": 2.5,
    "ytick.major.size": 4,
    "ytick.minor.size": 2.5,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


SERIES_COLORS = ["#9a9c9e", "#ff0d0d", "#1f77b4", "#ff8000"]
SERIES_MARKERS = ["s", "o", "o", "o"]
FIGURE_SIZE = (3.8, 2.8)
GRID_FIGURE_SIZE = (7.6, 5.6)
LEGEND_KWARGS = {
    "frameon": False,
    "fontsize": 10,
    "handlelength": 2,
    "labelspacing": 0.2,
    "borderpad": 0.1,
    "handletextpad": 0.2,
}


def apply_plot_style() -> None:
    rcParams.update(PLOT_STYLE)


def display_label(label: str) -> str:
    return "HfO2" if label == "undoped-HfO2" else label


def sparse_integer_fields(x_min: float, x_max: float, max_markers: int = 8) -> np.ndarray:
    first = math.ceil(x_min)
    last = math.floor(x_max)
    if last < first:
        return np.array([], dtype=float)
    count = last - first + 1
    step = max(1, math.ceil(count / max_markers))
    return np.arange(first, last + 1, step, dtype=float)


def field_profile_fields(coercive_field: float) -> np.ndarray:
    if not math.isfinite(coercive_field) or coercive_field <= 0:
        return np.array([0.0])

    integer_fields = sparse_integer_fields(0.0, math.floor(coercive_field), max_markers=5)
    fields = list(integer_fields)
    if not any(abs(field - coercive_field) < 1e-6 for field in fields):
        fields.append(float(coercive_field))
    return np.array(sorted(fields), dtype=float)


def color_ramp(base_color: str, count: int) -> list[tuple[float, float, float]]:
    rgb = np.array(to_rgb(base_color), dtype=float)
    if count <= 1:
        return [tuple(rgb)]
    white = np.ones(3)
    blend_amounts = np.linspace(0.0, 0.65, count)
    return [tuple((1.0 - amount) * rgb + amount * white) for amount in blend_amounts]


def style_axis(ax: plt.Axes, xlabel: str | None = None, ylabel: str | None = None) -> None:
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=14, labelpad=4)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14, labelpad=6)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="major", length=4, top=False, right=False)
    ax.tick_params(which="minor", length=2.5, top=False, right=False)
    ax.tick_params(axis="x", which="major", width=1.2, length=4, labelsize=12)
    ax.tick_params(axis="y", which="major", width=1.2, length=4, labelsize=12)
    ax.tick_params(axis="x", which="minor", width=1.0, length=2.5)
    ax.tick_params(axis="y", which="minor", width=1.0, length=2.5)


def plot_curve(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    label: str,
    color: str,
    marker: str,
    *,
    smooth: bool,
    marker_x: np.ndarray | None = None,
    linestyle: str = "-",
    fill_marker: bool = True,
    alpha: float = 1.0,
    zorder: int = 3,
    clip_on: bool = False,
) -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    x_clean = x[mask]
    y_clean = y[mask]
    if len(x_clean) == 0:
        return

    mfc = color if fill_marker else "white"
    marker_positions = x_clean if marker_x is None else marker_x[np.isfinite(marker_x)]

    if smooth and len(x_clean) >= 4 and np.all(np.diff(x_clean) > 0):
        x_line = np.linspace(float(x_clean.min()), float(x_clean.max()), 300)
        y_line = PchipInterpolator(x_clean, y_clean)(x_line)
        markevery = [int(np.abs(x_line - xi).argmin()) for xi in marker_positions]
    else:
        x_line = x_clean
        y_line = y_clean
        markevery = [int(np.abs(x_line - xi).argmin()) for xi in marker_positions]

    ax.plot(
        x_line,
        y_line,
        label=label,
        color=color,
        linestyle=linestyle,
        marker=marker,
        markevery=markevery,
        markersize=6,
        linewidth=1.6,
        markerfacecolor=mfc,
        markeredgecolor=color,
        markeredgewidth=1.2,
        clip_on=clip_on,
        zorder=zorder,
        alpha=alpha,
    )


def vector_columns(prefix: str, vec: np.ndarray) -> dict[str, float]:
    return {
        f"{prefix}_x": float(vec[0]),
        f"{prefix}_y": float(vec[1]),
        f"{prefix}_z": float(vec[2]),
    }


def matrix_columns(prefix: str, mat: np.ndarray) -> dict[str, float]:
    rows: dict[str, float] = {}
    for i, row_name in enumerate(["a", "b", "c"]):
        for j, axis in enumerate(["x", "y", "z"]):
            rows[f"{prefix}_{row_name}_{axis}"] = float(mat[i, j])
    return rows


def write_outputs(results: list[PathResult], output_dir: Path, field_step: float, field_max: float | None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    finite_fields = [r.coercive_field_MV_cm for r in results if math.isfinite(r.coercive_field_MV_cm)]
    if field_max is None:
        field_max = max(finite_fields, default=20.0) * 1.25
        field_max = max(field_max, 5.0)
    field_grid = np.arange(0.0, field_max + 0.5 * field_step, field_step)

    raw_rows = []
    branch_rows = []
    projection_rows = []
    summary_rows = []
    barrier_rows = []
    profile_rows = []
    polarization_rows = []

    for result in results:
        first_branch = result.branch_images[0]
        last_branch = result.branch_images[-1]
        summary_rows.append(
            {
                "label": result.label,
                "path": str(result.path),
                "notes": "; ".join(result.notes),
                "branch_mode": result.branch_mode,
                "barrier_mode": result.barrier_mode,
                "continuity_threshold_uC_cm2": result.continuity_threshold_uC_cm2,
                "saddle_image": result.saddle_image,
                "zero_field_barrier_eV": result.saddle_energy_eV,
                "coercive_field_MV_cm": result.coercive_field_MV_cm,
                "limiting_image": result.limiting_image,
                "endpoint_rel_energy_eV": result.endpoint_rel_energy_eV,
                "endpoint_projected_dipole_eA": result.endpoint_projected_dipole_eA,
                **vector_columns("field_dir", result.field_direction),
                **vector_columns("branch_dir", result.branch_direction),
                **vector_columns("P_initial_uC_cm2", first_branch.branch_polarization_uC_cm2),
                **vector_columns("P_final_uC_cm2", last_branch.branch_polarization_uC_cm2),
                **vector_columns("target_P_initial_uC_cm2", result.target_initial_polarization_uC_cm2),
                **vector_columns("target_P_final_uC_cm2", result.target_final_polarization_uC_cm2),
            }
        )

        switch_axis = switching_axis(result)
        target_delta_mu_rows = target_delta_dipoles(result)
        initial_target_p = target_polarization(result, 0)
        for image_index, (raw, branch) in enumerate(zip(result.raw_images, result.branch_images)):
            target_p = target_polarization(result, image_index)
            delta_mu = target_delta_mu_rows[image_index]
            delta_p = target_p - initial_target_p
            branch_delta_mu = branch.branch_dipole_eA - first_branch.branch_dipole_eA
            branch_delta_p = branch.branch_polarization_uC_cm2 - first_branch.branch_polarization_uC_cm2
            delta_mu_parallel = float(delta_mu @ result.field_direction)
            delta_p_parallel = float(delta_p @ result.field_direction)
            slope = -EV_PER_EANG_MV_CM * delta_mu_parallel
            requirement = math.inf
            if raw.rel_energy_eV > 0 and delta_mu_parallel > 0:
                requirement = raw.rel_energy_eV / (EV_PER_EANG_MV_CM * delta_mu_parallel)

            raw_rows.append(
                {
                    "label": result.label,
                    "image": raw.image,
                    "sample_path": str(raw.sample_path),
                    "energy_path": str(raw.energy_path),
                    "structure_path": str(raw.structure_path),
                    "dipole_outcar_path": str(raw.dipole_outcar_path),
                    "energy_eV": raw.energy_eV,
                    "rel_energy_eV": raw.rel_energy_eV,
                    "volume_A3": raw.volume_A3,
                    **matrix_columns("lattice_A", raw.lattice_A),
                    **matrix_columns("polarization_quantum_uC_cm2", raw.quantum_vectors_uC_cm2),
                    **vector_columns("ion_dipole_eA", raw.ion_dipole_eA),
                    **vector_columns("electronic_dipole_eA", raw.electronic_dipole_eA),
                    **vector_columns("raw_dipole_eA", raw.raw_dipole_eA),
                    **vector_columns("raw_P_uC_cm2", raw.raw_polarization_uC_cm2),
                }
            )
            branch_rows.append(
                {
                    "label": result.label,
                    "image": raw.image,
                    "branch_n1": int(branch.branch_index[0]),
                    "branch_n2": int(branch.branch_index[1]),
                    "branch_n3": int(branch.branch_index[2]),
                    "branch_projection_uC_cm2": branch.branch_projection_uC_cm2,
                    "delta_branch_projection_from_prev_uC_cm2": branch.delta_branch_projection_from_prev_uC_cm2,
                    "max_abs_component_jump_from_prev_uC_cm2": branch.max_abs_component_jump_from_prev_uC_cm2,
                    "component_jump_ok_x": bool(
                        abs(branch.delta_branch_polarization_from_prev_uC_cm2[0])
                        <= result.continuity_threshold_uC_cm2
                    ),
                    "component_jump_ok_y": bool(
                        abs(branch.delta_branch_polarization_from_prev_uC_cm2[1])
                        <= result.continuity_threshold_uC_cm2
                    ),
                    "component_jump_ok_z": bool(
                        abs(branch.delta_branch_polarization_from_prev_uC_cm2[2])
                        <= result.continuity_threshold_uC_cm2
                    ),
                    "monotonic_along_branch": branch.monotonic_along_branch,
                    "component_continuity_ok": branch.component_continuity_ok,
                    **vector_columns("branch_P_uC_cm2", branch.branch_polarization_uC_cm2),
                    **vector_columns("branch_dipole_eA", branch.branch_dipole_eA),
                    **vector_columns(
                        "delta_branch_P_from_prev_uC_cm2",
                        branch.delta_branch_polarization_from_prev_uC_cm2,
                    ),
                }
            )
            projection_rows.append(
                {
                    "label": result.label,
                    "image": raw.image,
                    "rel_energy_eV": raw.rel_energy_eV,
                    "delta_mu_parallel_eA": delta_mu_parallel,
                    "delta_P_parallel_uC_cm2": delta_p_parallel,
                    "field_energy_slope_eV_per_MV_cm": slope,
                    "required_field_MV_cm": requirement,
                    "can_be_lowered_by_field": bool(delta_mu_parallel > 0),
                    **vector_columns("delta_mu_eA", delta_mu),
                    **vector_columns("delta_P_uC_cm2", delta_p),
                    **vector_columns("branch_delta_mu_eA", branch_delta_mu),
                    **vector_columns("branch_delta_P_uC_cm2", branch_delta_p),
                }
            )
            polarization_rows.append(
                {
                    "label": result.label,
                    "image": raw.image,
                    "P_parallel_field_uC_cm2": float(branch.branch_polarization_uC_cm2 @ result.field_direction),
                    "P_parallel_branch_uC_cm2": branch.branch_projection_uC_cm2,
                    "target_P_switching_uC_cm2": float(target_p @ switch_axis),
                    "target_P_parallel_field_uC_cm2": float(target_p @ result.field_direction),
                    "target_P_parallel_branch_uC_cm2": float(target_p @ result.branch_direction),
                    **vector_columns("target_P_uC_cm2", target_p),
                    **vector_columns("branch_P_uC_cm2", branch.branch_polarization_uC_cm2),
                }
            )

        for field in field_grid:
            barrier_rows.append(
                {
                    "label": result.label,
                    "field_MV_cm": float(field),
                    "barrier_eV": barrier_at_field(result, float(field)),
                }
            )
        profile_fields = [0.0]
        if math.isfinite(result.coercive_field_MV_cm):
            profile_fields.append(result.coercive_field_MV_cm)
        for field in profile_fields:
            corrected = corrected_relative_energies(result, field)
            for raw, value in zip(result.raw_images, corrected):
                profile_rows.append(
                    {
                        "label": result.label,
                        "field_MV_cm": float(field),
                        "image": raw.image,
                        "corrected_rel_energy_eV": float(value),
                    }
                )

    pd.DataFrame(summary_rows).to_csv(output_dir / "coercive_field_summary.csv", index=False)
    pd.DataFrame(raw_rows).to_csv(output_dir / "raw_vasp_data.csv", index=False)
    pd.DataFrame(branch_rows).to_csv(output_dir / "branch_correction.csv", index=False)
    pd.DataFrame(projection_rows).to_csv(output_dir / "field_projection.csv", index=False)
    pd.DataFrame(branch_rows).to_csv(output_dir / "polarization_branches.csv", index=False)
    pd.DataFrame(barrier_rows).to_csv(output_dir / "barrier_vs_field.csv", index=False)
    pd.DataFrame(profile_rows).to_csv(output_dir / "corrected_profiles.csv", index=False)
    pd.DataFrame(polarization_rows).to_csv(output_dir / "polarization_profiles.csv", index=False)

    plot_barriers(pd.DataFrame(barrier_rows), results, output_dir / "barrier_vs_field.png")
    plot_profiles(results, output_dir / "energy_profiles.png")
    plot_field_profile_families(results, output_dir / "field_profiles.png")
    plot_polarization_profiles(results, output_dir / "polarization_profiles.png")
    plot_branch_polarization_profiles(results, output_dir / "branch_polarization_profiles.png")


def plot_barriers(df: pd.DataFrame, results: list[PathResult], output: Path) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for idx, result in enumerate(results):
        sub = df[df["label"] == result.label]
        color = SERIES_COLORS[idx % len(SERIES_COLORS)]
        marker = SERIES_MARKERS[idx % len(SERIES_MARKERS)]
        label = display_label(result.label)
        fields = sub["field_MV_cm"].to_numpy(dtype=float)
        barriers = sub["barrier_eV"].to_numpy(dtype=float)
        x_min = float(fields.min())
        x_max = float(fields.max())
        integer_fields = sparse_integer_fields(x_min, x_max)
        plot_curve(
            ax,
            fields,
            barriers,
            label,
            color,
            marker,
            smooth=False,
            marker_x=integer_fields,
            clip_on=True,
        )
        if math.isfinite(result.coercive_field_MV_cm):
            ec_barrier = barrier_at_field(result, result.coercive_field_MV_cm)
            if abs(ec_barrier) < 1e-10:
                ec_barrier = 0.0
            ax.scatter(
                [result.coercive_field_MV_cm],
                [ec_barrier],
                color=color,
                marker=marker,
                s=25,
                linewidths=1.2,
                zorder=3,
            )
            # ax.axvline(result.coercive_field_MV_cm, color=color, linewidth=1, zorder=0, alpha=0.5, linestyle="--")

    ax.axhline(0.0, color="0.65", linewidth=0.8)
    ax.set_ylim(bottom=0.0)
    style_axis(ax, "Electric field (MV/cm)", "Switching barrier (eV/cell)")
    ax.legend(loc="upper right", **LEGEND_KWARGS)
    fig.tight_layout(pad=0.2)
    fig.savefig(output, transparent=False)
    plt.close(fig)


def plot_profiles(results: list[PathResult], output: Path) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=GRID_FIGURE_SIZE, sharex=True)
    axes_flat = axes.ravel()
    for panel_index, (ax, result) in enumerate(zip(axes_flat, results)):
        images = np.array([img.image for img in result.raw_images], dtype=float)
        color = SERIES_COLORS[panel_index % len(SERIES_COLORS)]
        marker = SERIES_MARKERS[panel_index % len(SERIES_MARKERS)]
        plot_curve(
            ax,
            images,
            corrected_relative_energies(result, 0.0),
            "0 MV/cm",
            color,
            marker,
            smooth=True,
            fill_marker=False,
        )
        if math.isfinite(result.coercive_field_MV_cm):
            corrected = corrected_relative_energies(result, result.coercive_field_MV_cm)
            plot_curve(
                ax,
                images,
                corrected,
                f"{result.coercive_field_MV_cm:.2f} MV/cm",
                color,
                marker,
                smooth=True,
            )
        ax.axhline(0.0, color="0.75", linewidth=0.8)
        ax.set_title(display_label(result.label), fontsize=12)
        ylabel = "Rel. free energy (eV/cell)" if panel_index % 2 == 0 else None
        style_axis(ax, ylabel=ylabel)
        ax.legend(loc="lower left", **LEGEND_KWARGS)
    for ax in axes_flat[-2:]:
        ax.set_xlabel("NEB image", fontsize=14, labelpad=4)
    for ax in axes_flat[len(results) :]:
        ax.axis("off")
    fig.tight_layout(pad=0.2)
    fig.savefig(output, transparent=False)
    plt.close(fig)


def plot_field_profile_families(results: list[PathResult], output: Path) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=GRID_FIGURE_SIZE, sharex=True)
    axes_flat = axes.ravel()
    for panel_index, (ax, result) in enumerate(zip(axes_flat, results)):
        images = np.array([img.image for img in result.raw_images], dtype=float)
        fields = field_profile_fields(result.coercive_field_MV_cm)
        colors = color_ramp(SERIES_COLORS[panel_index % len(SERIES_COLORS)], len(fields))
        marker = SERIES_MARKERS[panel_index % len(SERIES_MARKERS)]
        for idx, (field, color) in enumerate(zip(fields, colors)):
            corrected = corrected_relative_energies(result, float(field))
            label = "E=0" if field == 0 else f"{field:.2f} MV/cm"
            plot_curve(
                ax,
                images,
                corrected,
                label,
                color,
                marker,
                smooth=True,
                alpha=1.0,
                zorder=idx + 1,
            )
        ax.axhline(0.0, color="0.75", linewidth=0.8, zorder=0)
        ax.set_title(display_label(result.label), fontsize=12)
        ylabel = "Energy barrier (eV)" if panel_index % 2 == 0 else None
        style_axis(ax, ylabel=ylabel)
        ax.legend(loc="lower left", **LEGEND_KWARGS)
    for ax in axes_flat[-2:]:
        ax.set_xlabel("Coordinate path", fontsize=14, labelpad=4)
    for ax in axes_flat[len(results) :]:
        ax.axis("off")
    fig.tight_layout(pad=0.2)
    fig.savefig(output, transparent=False)
    plt.close(fig)


def plot_polarization_profiles(results: list[PathResult], output: Path) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=GRID_FIGURE_SIZE, sharex=True)
    axes_flat = axes.ravel()
    for panel_index, (ax, result) in enumerate(zip(axes_flat, results)):
        images = np.array([img.image for img in result.raw_images], dtype=float)
        axis = switching_axis(result)
        polarization = np.array(
            [target_polarization(result, idx) @ axis for idx in range(len(result.raw_images))],
            dtype=float,
        )
        color = SERIES_COLORS[panel_index % len(SERIES_COLORS)]
        marker = SERIES_MARKERS[panel_index % len(SERIES_MARKERS)]
        plot_curve(
            ax,
            images,
            polarization,
            display_label(result.label),
            color,
            marker,
            smooth=False,
        )
        ax.axhline(0.0, color="0.75", linewidth=0.8, zorder=0)
        ax.set_title(display_label(result.label), fontsize=12)
        ylabel = r"Polarization ($\mu$C/cm$^2$)" if panel_index % 2 == 0 else None
        style_axis(ax, ylabel=ylabel)
    for ax in axes_flat[-2:]:
        ax.set_xlabel("NEB image", fontsize=14, labelpad=4)
    for ax in axes_flat[len(results) :]:
        ax.axis("off")
    fig.tight_layout(pad=0.2)
    fig.savefig(output, transparent=False)
    plt.close(fig)


def plot_branch_polarization_profiles(results: list[PathResult], output: Path) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=GRID_FIGURE_SIZE, sharex=True)
    axes_flat = axes.ravel()
    for panel_index, (ax, result) in enumerate(zip(axes_flat, results)):
        images = np.array([img.image for img in result.raw_images], dtype=float)
        polarization = np.array(
            [img.branch_polarization_uC_cm2 @ result.field_direction for img in result.branch_images],
            dtype=float,
        )
        color = SERIES_COLORS[panel_index % len(SERIES_COLORS)]
        marker = SERIES_MARKERS[panel_index % len(SERIES_MARKERS)]
        plot_curve(
            ax,
            images,
            polarization,
            display_label(result.label),
            color,
            marker,
            smooth=True,
        )
        ax.axhline(0.0, color="0.75", linewidth=0.8, zorder=0)
        ax.set_title(display_label(result.label), fontsize=12)
        ylabel = r"Branch $P_\parallel$ ($\mu$C/cm$^2$)" if panel_index % 2 == 0 else None
        style_axis(ax, ylabel=ylabel)
    for ax in axes_flat[-2:]:
        ax.set_xlabel("NEB image", fontsize=14, labelpad=4)
    for ax in axes_flat[len(results) :]:
        ax.axis("off")
    fig.tight_layout(pad=0.2)
    fig.savefig(output, transparent=False)
    plt.close(fig)
