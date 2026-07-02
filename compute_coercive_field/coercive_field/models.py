"""Data models for the coercive-field workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class RawImage:
    """One NEB image after reading the raw VASP files.

    这个对象只保存“原始文件直接读到或直接换算得到”的量：
    OSZICAR/neb.dat 给零场能量，CONTCAR/POSCAR 给晶格和体积，
    1-dipol/OUTCAR 给 LCALCPOL 偶极矩。
    
    """

    label: str # 用于指定样本的标签名，例如 undoped-HfO2、a-HZO、b-HZO、c-HZO 等
    sample_path: Path # 指定相应标签的样本路径
    image: int # 该图像在 NEB 路径中的索引，从 0 开始，通过循环进行依次读取
    energy_path: Path # 能量信息读取的文件的路径
    structure_path: Path # 结构信息读取的文件的路径
    dipole_outcar_path: Path # 偶极矩信息读取的 OUTCAR 文件的路径
    energy_eV: float # 该图像的能量值，单位为 eV
    rel_energy_eV: float # 该图像的相对能量值，相对于路径中第一个图像的能量，单位为 eV
    lattice_A: np.ndarray # 3x3 的 numpy 数组，表示晶格向量，单位为 Angstrom
    volume_A3: float # 晶格体积，单位为 Angstrom^3
    ion_dipole_eA: np.ndarray # 3x1 的 numpy 数组，表示离子偶极矩，单位为 e*Angstrom
    electronic_dipole_eA: np.ndarray # 3x1 的 numpy 数组，表示电子偶极矩，单位为 e*Angstrom
    raw_dipole_eA: np.ndarray # 3x1 的 numpy 数组，表示原始偶极矩，即离子偶极矩和电子偶极矩的总和，单位为 e*Angstrom
    raw_polarization_uC_cm2: np.ndarray # 3x1 的 numpy 数组，表示原始极化，单位为 uC/cm^2，未经过平移时的极化值
    quantum_vectors_uC_cm2: np.ndarray # 3x3 的 numpy 数组，表示极化量子向量，单位为 uC/cm^2


@dataclass
class BranchImage:
    """One NEB image after choosing a continuous Berry polarization branch.

    VASP 的极化只在 polarization quantum 模意义下确定，所以这里记录
    每个 image 加了多少整数分支、修正后的极化/偶极矩，以及相邻
    image 的连续性诊断。
    """

    image: int
    branch_index: np.ndarray
    branch_polarization_uC_cm2: np.ndarray
    branch_dipole_eA: np.ndarray
    branch_projection_uC_cm2: float
    delta_branch_polarization_from_prev_uC_cm2: np.ndarray
    delta_branch_projection_from_prev_uC_cm2: float
    max_abs_component_jump_from_prev_uC_cm2: float
    monotonic_along_branch: bool
    component_continuity_ok: bool


@dataclass
class PathResult:
    """All data and final coercive-field result for one sample path.

    raw_images 和 branch_images 是逐 image 的数据；其余字段描述
    外场方向、分支方向、目标极化路径、矫顽场判据和最终结果。
    """

    label: str
    path: Path
    raw_images: list[RawImage]
    branch_images: list[BranchImage]
    field_direction: np.ndarray
    branch_direction: np.ndarray
    target_initial_polarization_uC_cm2: np.ndarray
    target_final_polarization_uC_cm2: np.ndarray
    branch_mode: str
    barrier_mode: str
    coercive_criterion: str
    continuity_threshold_uC_cm2: float
    saddle_image: int
    saddle_energy_eV: float
    coercive_field_MV_cm: float
    limiting_image: int | None
    endpoint_rel_energy_eV: float
    endpoint_projected_dipole_eA: float
    notes: list[str] = field(default_factory=list)
