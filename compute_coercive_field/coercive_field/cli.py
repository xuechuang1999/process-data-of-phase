from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

from .calculation import analyze_path
from .outputs import write_outputs


def default_samples(root: Path) -> list[tuple[str, Path]]:

    import platform

    if platform.system() == "Windows":
        translation_path_root = str(root).split("\\")[-1]
    else:
        translation_path_root = str(root).split("/")[-1]

    return [
        ("HfO2", root),
        ("a-HZO", root / f"a-{translation_path_root}"),
        ("b-HZO", root / f"b-{translation_path_root}"),
        ("c-HZO", root / f"c-{translation_path_root}"),
    ]


def parse_sample_args(samples: Iterable[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    for item in samples:
        if "=" not in item:
            raise argparse.ArgumentTypeError("--sample must be LABEL=PATH/to/sample, such as --sample HfO2=path/to/HfO2")
        label, raw_path = item.split("=", 1)
        parsed.append((label, Path(raw_path)))
    return parsed


def build_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description="Compute coercive fields from VASP NEB and LCALCPOL outputs.")

    parser.add_argument("--root", type=Path, default=Path("Pmn2_1") / "xyz-101-000")

    parser.add_argument("--sample", action="append", default=[], help="Optional sample in LABEL=PATH form.")

    parser.add_argument("--energy-source", choices=["oszicar", "neb", "auto"], default="oszicar")

    parser.add_argument(
        "--initial-polarization",
        default="auto",
        help=(
            "Initial polarization in uC/cm^2. Use auto to infer the phase initial direction "
            "(Pca2_1: 001, Pmn2_1: 101 in lattice coordinates) and choose the matching VASP image 00 branch; "
            "or pass Px,Py,Pz."
        ),
    )

    parser.add_argument("--final-polarization", default="auto", help="Final polarization vector in uC/cm^2; use 'opposite-initial' or 'auto' to infer from the path.")

    parser.add_argument(
        "--continuity-threshold",
        default="auto",
        help="Max allowed component jump in uC/cm^2. Use auto for max(|polarization quantum|)/(2*image_count).",
    )

    parser.add_argument("--field-direction", default="0,0,-1", help="Specify the direction of the external field.")

    parser.add_argument("--field-coordinates", choices=["lattice", "cartesian"], default="lattice")

    parser.add_argument("--barrier-mode", choices=["to-saddle", "intermediate", "full-path"], default="full-path")

    parser.add_argument("--coercive-criterion", choices=["barrier", "downhill"], default="downhill")

    parser.add_argument("--field-step", type=float, default=0.1)

    parser.add_argument("--field-max", type=float, default=None)

    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    samples = parse_sample_args(args.sample) if args.sample else default_samples(args.root)

    results = [
        analyze_path(
            label=label,
            path=path,
            energy_source=args.energy_source,
            field_direction_arg=args.field_direction,
            field_coordinates=args.field_coordinates,
            barrier_mode=args.barrier_mode,
            coercive_criterion=args.coercive_criterion,
            initial_polarization=args.initial_polarization,
            final_polarization_arg=args.final_polarization,
            continuity_threshold_uC_cm2=args.continuity_threshold,
        )
        for label, path in samples
    ]
    write_outputs(results, args.output_dir, args.field_step, args.field_max)

    print("label, zero-field barrier (eV), coercive field (MV/cm), limiting constraint")
    for result in results:
        ec = "inf" if not math.isfinite(result.coercive_field_MV_cm) else f"{result.coercive_field_MV_cm:.4f}"
        limiting = result.limiting_image
        if result.coercive_criterion == "downhill" and result.limiting_image is not None:
            limiting = f"{result.limiting_image - 1}->{result.limiting_image}"
        print(f"{result.label}, {result.saddle_energy_eV:.6f}, {ec}, {limiting}")
        for note in result.notes:
            print(f"  note: {note}")
    print(f"wrote outputs to {args.output_dir.resolve()}")
