"""Command line entry point."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from .calculation import analyze_path
from .directions import parse_vector
from .outputs import write_outputs
from .final_polarization_direction import final_polarization_direction


def default_samples(root: Path) -> list[tuple[str, Path]]:
    """Return the default undoped/a/b/c sample paths under one NEB root.

    If root is ``Pca2_1/xyz-000-011``, doped samples are expected beside it as
    ``a-xyz-000-011``, ``b-xyz-000-011``, and ``c-xyz-000-011``.
    """

    import platform

    # The doped directory names reuse only the final path component of root.
    # Use platform-specific splitting because root may come from a raw CLI path.
    if platform.system() == "Windows":
        translation_path_root = str(root).split("\\")[-1]
    else:
        translation_path_root = str(root).split("/")[-1]

    # 默认要处理的样本路径，返回一个列表，每个元素是一个元组，包含样本标签和对应的路径
    return [
        ("undoped-HfO2", root),
        ("a-HZO", root / f"a-{translation_path_root}"),
        ("b-HZO", root / f"b-{translation_path_root}"),
        ("c-HZO", root / f"c-{translation_path_root}"),
    ]


# 根据输入的根目录，自动返回空间群和晶胞平移的表达式，方便后续自动解析终点极化方向
def parse_default_samples(root: Path) -> list[tuple[str, str]]:

    import platform
    if platform.system() == "Windows":
        space_group = str(root).split("\\")[0]
        translation_formula = str(root).split("\\")[1]
    else:
        space_group = str(root).split("/")[0]
        translation_formula = str(root).split("/")[1]

    return [space_group, translation_formula]

# 样本标签和路径的解析函数，通过--sample参数传入，传入方式为：--sample HfO2=path/to/HfO2 --sample a-HZO=path/to/a-HZO --sample b-HZO=path/to/b-HZO --sample c-HZO=path/to/c-HZO
def parse_sample_args(samples: Iterable[str]) -> list[tuple[str, Path]]:
    """Parse repeated ``--sample LABEL=PATH`` arguments."""

    parsed: list[tuple[str, Path]] = [] # initialize empty list, and specify the type for clarity
    for item in samples:
        if "=" not in item:
            raise argparse.ArgumentTypeError("--sample must be LABEL=PATH/to/sample, such as --sample HfO2=path/to/HfO2")
        label, raw_path = item.split("=", 1) # 使用split将字符串按照等号分割，并限制分割次数为1，得到样本标签和路径
        parsed.append((label, Path(raw_path)))
    return parsed

# 创建命令行解析器，定义所有的命令行参数和选项
def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface.
    创建命令行解析对象

    Most defaults mirror the README workflow: use OSZICAR energies, unfold
    Berry branches by three-component continuity, interpret directions in
    crystal-axis coordinates by default, apply a field along -c, and compute Ec
    with the downhill criterion over the full path.
    """

    parser = argparse.ArgumentParser(description="Compute coercive fields from VASP NEB and LCALCPOL outputs.") # create the configuration parser of all the command line arguments, create the parser object, and provide a description of the program, in the below lines are used to register the rules of parameters for the objects
    parser.add_argument("--root", type=Path, default=Path("Pca2_1") / "xyz-000-011") # 指定一个根目录，默认是 Pca2_1/xyz-000-011
    parser.add_argument("--sample", action="append", default=[], help="Optional sample in LABEL=PATH form.") # 指定传入的样本，格式为 --sample label=/path/to/sample，支持多次传入
    parser.add_argument("--energy-source", choices=["oszicar", "neb", "auto"], default="oszicar")#指定读取能量时使用的来源，支持 oszicar、neb、auto 三种方式，默认是 oszicar
    parser.add_argument("--initial-polarization", type=parse_vector, default=np.array([0.0, 0.0, 50.0])) # 指定初始极化方向，默认是 [0, 0, 50]，使用 parse_vector 函数解析输入的字符串为 numpy 数组
    parser.add_argument("--final-polarization", default="opposite-initial") # 指定最终极化方向，默认是 opposite-initial，表示与初始极化方向相反
    parser.add_argument(
        "--endpoint-direction",
        default="auto",
        help=(
            "Final-state direction used only to resolve equivalent Berry branches. "
            "Use auto, none, or x,y,z. Auto uses Pca2_1->{010} and Pmn2_1->[-110]."
        ),
    ) # 指定最终状态的方向，用于解决等效的 Berry 分支，默认是 auto，表示根据空间群自动选择方向
    parser.add_argument(
        "--endpoint-polarization",
        default="none",
        help="Explicit final branch polarization in uC/cm^2; overrides endpoint-direction magnitude.",
    ) # 指定最终分支的极化，单位是 uC/cm^2，默认是 none，表示不指定
    parser.add_argument("--branch-mode", choices=["center-zero", "endpoint-target"], default="center-zero") # 指定分支模式，支持 center-zero 和 endpoint-target 两种方式，默认是 center-zero
    parser.add_argument("--branch-direction", default="0,0,-1") # 指定分支方向，默认是 [0, 0, -1]
    parser.add_argument("--branch-coordinates", choices=["lattice-diagonal", "lattice", "cartesian"], default="lattice") # 指定分支坐标系，支持 lattice-diagonal、lattice、cartesian 三种方式，默认是 lattice
    parser.add_argument(
        "--branch-search-radius",
        type=int,
        default=2,
        help="Deprecated compatibility option; direct quantum-coordinate branch selection does not enumerate candidates.",
    ) # 指定分支搜索半径，默认是 2，表示在量子坐标系中选择分支时不枚举候选分支
    parser.add_argument("--continuity-threshold", type=float, default=80.0, help="Max allowed component jump in uC/cm^2.") # 指定连续性阈值，单位是 uC/cm^2，默认是 80.0，表示允许的最大分量跳跃
    parser.add_argument("--field-direction", default="0,0,-1", help="Specify the direction of the external field.") # 指定电场方向，默认是 [0, 0, -1]
    parser.add_argument("--field-coordinates", choices=["lattice-diagonal", "lattice", "cartesian"], default="lattice") # 指定电场坐标系，支持 lattice-diagonal、lattice、cartesian 三种方式，默认是 lattice
    parser.add_argument("--barrier-mode", choices=["to-saddle", "intermediate", "full-path"], default="full-path") # 指定barrier的计算方式，支持 to-saddle、intermediate、full-path 三种方式，默认是 full-path
    parser.add_argument("--coercive-criterion", choices=["barrier", "downhill"], default="downhill") # 指定矫顽场判据方式，支持 barrier 和 downhill 两种方式，默认是 downhill
    parser.add_argument("--field-step", type=float, default=0.1) # 指定电场步长，单位是 MV/cm，默认是 0.1
    parser.add_argument("--field-max", type=float, default=None) # 指定最大电场，单位是 MV/cm，默认是 None，表示不限制最大电场
    parser.add_argument("--output-dir", type=Path, default=Path("results")) # 指定输出目录，默认是 results
    return parser


def main() -> None:
    """Run the full command-line workflow."""

    parser = build_parser()
    args = parser.parse_args()
    samples = parse_sample_args(args.sample) if args.sample else default_samples(args.root)
    # Analyze every sample independently, then write combined CSV/figure output.
    results = [
        analyze_path(
            label=label,
            path=path,
            energy_source=args.energy_source,
            field_direction_arg=args.field_direction,
            field_coordinates=args.field_coordinates,
            branch_direction_arg=args.branch_direction,
            branch_coordinates=args.branch_coordinates,
            branch_mode=args.branch_mode,
            barrier_mode=args.barrier_mode,
            coercive_criterion=args.coercive_criterion,
            initial_polarization=args.initial_polarization,
            final_polarization_arg=args.final_polarization,
            endpoint_direction_arg=args.endpoint_direction,
            endpoint_polarization_arg=args.endpoint_polarization,
            continuity_threshold_uC_cm2=args.continuity_threshold,
            branch_search_radius=args.branch_search_radius,
        )
        for label, path in samples
    ]
    write_outputs(results, args.output_dir, args.field_step, args.field_max)

    # Keep terminal output intentionally compact; detailed per-image data is in
    # the CSV files.
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
