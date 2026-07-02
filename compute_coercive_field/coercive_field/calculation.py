"""Coercive-field calculations.

本模块负责把已经读入并完成 Berry branch 修正的 VASP 数据转成外场
自由能曲线和矫顽场。核心数学链条如下：

1. 零场 NEB 相对能量:
       U_i = E_i - E_0

2. Berry branch 修正后的偶极矩变化:
       Delta_mu_i = mu_branch_i - mu_branch_0

3. 将偶极矩变化投影到外场方向:
       Delta_mu_parallel_i = Delta_mu_i dot e_hat

4. 外场修正后的相对自由能:
       E_vec = E * e_hat
       F_i(E_vec) = U_i - 0.01 * (Delta_mu_i dot E_vec)

   其中 E_vec 的单位是 MV/cm，Delta_mu 的单位是 e*Angstrom，
   0.01 = EV_PER_EANG_MV_CM，用来把二者乘积换算成 eV。

5. 矫顽场 Ec 不是靠扫描网格取近似值，而是对每个约束解析求解。
   每个约束都写成:
       Delta_E - 0.01 * E * Delta_mu_parallel <= 0

   若 Delta_E > 0 且 Delta_mu_parallel > 0:
       E_required = Delta_E / (0.01 * Delta_mu_parallel)

   最终:
       Ec = max(E_required_i)

   因为所有能垒/上坡约束必须同时被压低到不大于 0，所以取最大值。
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .branch import build_branch_images
from .constants import EV_PER_EANG_MV_CM, UC_CM2_PER_EANG_A3
from .directions import resolve_direction
from .models import PathResult
from .vasp_io import parse_structure, read_raw_images


def resolve_final_polarization(text: str, initial_polarization: np.ndarray) -> np.ndarray:
    """Resolve the final target polarization from CLI input."""

    if text == "opposite-initial":
        return -initial_polarization
    return np.array([float(part.strip()) for part in text.split(",")], dtype=float)


def infer_endpoint_direction(path: Path) -> np.ndarray | None:
    """Infer the non-180 endpoint direction from the phase name when possible.

    This direction is used only to resolve the final Berry branch in
    center-zero mode. It is not a target path for intermediate images.
    """

    text = str(path)
    if "Pca2_1" in text:
        return np.array([0.0, 1.0, 0.0])
    if "Pmn2_1" in text:
        direction = np.array([-1.0, 1.0, 0.0])
        return direction / np.linalg.norm(direction)
    return None


# 用来自动计算两个相邻的image之间可以接受的最大差值
def compute_maxium_difference_between_adjacent_images(quantum_polarization:np.ndarry, n_images:int) -> float:
    max_quantum_polarization = np.max(quantum_polarization)
    return max_quantum_polarization / (n_images - 1)


def resolve_endpoint_polarization(text: str) -> np.ndarray | None:
    """Resolve an explicit endpoint polarization vector from CLI input."""
    # 检测用户输入的终点极化参数，如果是 "auto" 或 "none"，则返回 None，否则将其解析为一个 numpy 数组，也是一个三维向量，表示终点极化的大小和方向
    if text in {"auto", "none"}:
        return None
    return np.array([float(part.strip()) for part in text.split(",")], dtype=float)


def resolve_endpoint_direction(text: str, path: Path) -> np.ndarray | None:
    """Resolve the endpoint direction used to choose among continuous branches."""

    if text == "none":
        return None
    if text == "auto":
        # 查看路径中是否包含 Pca2_1 或 Pmn2_1 等信息，如果有则返回对应的方向向量，否则返回 None
        return infer_endpoint_direction(path)
    direction = np.array([float(part.strip()) for part in text.split(",")], dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm == 0:
        raise ValueError("endpoint direction has zero norm")
    return direction / norm  #返回归一化之后的用户输入的方向向量，3x1的numpy数组


def target_polarization(result: PathResult, image_index: int) -> np.ndarray:
    """Return the diagnostic target switching polarization at one image.

    这个函数只生成“用户设想的线性极化路径”，用于画图和诊断，不参与
    外场能量修正。真正进入 F_i(E) 的是 LCALCPOL 数据经 Berry branch
    修正后的 branch_dipole。

    Formula:
        t_i = i / (N - 1)
        P_target_i = P_initial + t_i * (P_final - P_initial)

    Meaning of symbols:
        i: image_index, NEB image number.
        N: len(result.raw_images), total number of NEB images.
        P_initial: result.target_initial_polarization_uC_cm2, target initial
            polarization in uC/cm^2.
        P_final: result.target_final_polarization_uC_cm2, target final
            polarization in uC/cm^2.
        P_target_i: returned vector, target polarization at image i.
    """

    image_count = len(result.raw_images)
    # 反应坐标的线性进度 t_i。只有一个 image 时退化为 t_i = 0。
    progress = image_index / (image_count - 1) if image_count > 1 else 0.0
    # 线性插值得到诊断目标极化 P_target_i。
    return result.target_initial_polarization_uC_cm2 + progress * (
        result.target_final_polarization_uC_cm2 - result.target_initial_polarization_uC_cm2
    )


def switching_axis(result: PathResult) -> np.ndarray:
    """Return the axis used for plotting the target switching profile.

    Formula:
        e_switch = P_initial / |P_initial|

    Meaning of symbols:
        P_initial: result.target_initial_polarization_uC_cm2.
        e_switch: returned unit vector. If P_initial is zero, the code falls
            back to result.branch_direction.
    """

    axis = result.target_initial_polarization_uC_cm2
    norm = float(np.linalg.norm(axis))
    if norm == 0.0:
        return result.branch_direction
    return axis / norm


def target_delta_dipoles(result: PathResult) -> np.ndarray:
    """Convert diagnostic target polarization changes to target dipole changes.

    这里处理的是 target path 的诊断量，不用于 Ec。

    极化和偶极矩换算关系:
        P = UC_CM2_PER_EANG_A3 * mu / Omega

    因此:
        Delta_mu = Delta_P * Omega / UC_CM2_PER_EANG_A3

    这里使用每个 image 自己的 Omega_i，因为 NEB images 的弛豫晶胞体积
    可能略有不同。

    Formula:
        Delta_P_i = P_target_i - P_target_0
        Delta_mu_i = Delta_P_i * Omega_i / UC_CM2_PER_EANG_A3

    Meaning of symbols:
        Delta_P_i: target polarization change at image i, in uC/cm^2.
        P_target_i: target_polarization(result, i).
        Omega_i: raw.volume_A3, cell volume of image i in Angstrom^3.
        UC_CM2_PER_EANG_A3: conversion factor from
            e*Angstrom/Angstrom^3 to uC/cm^2.
        Delta_mu_i: returned target dipole change, in e*Angstrom.
    """

    initial_target = target_polarization(result, 0)
    rows = []
    for image_index, raw in enumerate(result.raw_images):
        # Delta_P_target_i = P_target_i - P_target_0。
        delta_p = target_polarization(result, image_index) - initial_target
        # Delta_mu_target_i = Delta_P_target_i * Omega_i / C_P。
        rows.append(delta_p * raw.volume_A3 / UC_CM2_PER_EANG_A3)
    return np.array(rows, dtype=float)


def branch_delta_dipoles(result: PathResult) -> np.ndarray:
    """Return Berry-branch-corrected dipole changes relative to image 00.

    这是外场能量修正真正使用的偶极矩变化。

    Formula:
        Delta_mu_i = mu_branch_i - mu_branch_0

    Meaning of symbols:
        mu_branch_i: result.branch_images[i].branch_dipole_eA, the dipole
            corresponding to the Berry-branch-corrected vector polarization.
        Delta_mu_i: returned vector, used directly for electric-field energy
            correction.
    """

    # 以 image 00 作为参考态，所有外场修正都只使用相对偶极矩变化。
    first = result.branch_images[0].branch_dipole_eA
    return np.array([image.branch_dipole_eA - first for image in result.branch_images], dtype=float)


def branch_delta_polarizations(result: PathResult) -> np.ndarray:
    """Return Berry-branch-corrected polarization changes relative to image 00.

    这个量主要写入 CSV 供检查；能量修正等价地使用 branch_delta_dipoles。

    Formula:
        Delta_P_i = P_branch_i - P_branch_0

    Meaning of symbols:
        P_branch_i: result.branch_images[i].branch_polarization_uC_cm2.
        Delta_P_i: branch-corrected vector polarization change in uC/cm^2.
    """

    # Delta_P_i = P_branch_i - P_branch_0。
    first = result.branch_images[0].branch_polarization_uC_cm2
    return np.array([image.branch_polarization_uC_cm2 - first for image in result.branch_images], dtype=float)


def corrected_relative_energies(result: PathResult, field_MV_cm: float) -> np.ndarray:
    """Apply the external-field correction to all relative NEB energies.

    本函数是外场能量修正的核心。

    代码采用显式矢量点乘写法:
        E_vec = E * e_hat
        F_i(E_vec) = U_i - 0.01 * (Delta_mu_i dot E_vec)

    单位:
        E_vec: MV/cm
        Delta_mu_i: e*Angstrom
        Delta_mu_i dot E_vec: e*Angstrom*MV/cm
        0.01 * (Delta_mu_i dot E_vec): eV
        F_i(E): eV

    这里的 U_i 已经是相对 image 00 的零场能量，所以 F_i(E) 也是相对
    image 00 的外场修正能量。image 00 的 Delta_mu_0 = 0，因此 F_0(E)=U_0。

    Equivalent vector form:
        F_i(E_vec) = U_i - Omega_i * Delta_P_i dot E_vec
        E_vec = E * e_hat

    Code formula:
        E_vec = field_MV_cm * e_hat
        field_work_i = Delta_mu_i dot E_vec
        F_i(E_vec) = U_i - EV_PER_EANG_MV_CM * field_work_i

    Meaning of symbols:
        U_i: img.rel_energy_eV, zero-field relative NEB energy in eV.
        E: field_MV_cm, external-field magnitude in MV/cm.
        e_hat: result.field_direction, Cartesian unit field direction.
        E_vec: Cartesian electric-field vector in MV/cm.
        Delta_mu_i: branch_delta_dipoles(result)[i], Berry-branch-corrected
            dipole change relative to image 00 in e*Angstrom.
        field_work_i: vector dot product Delta_mu_i dot E_vec.
        EV_PER_EANG_MV_CM: 0.01 eV per (e*Angstrom*MV/cm).
        F_i(E): returned corrected relative free energy in eV.
    """

    # U_i: zero-field relative NEB energies.
    rel = np.array([img.rel_energy_eV for img in result.raw_images])
    # E_vec = E * e_hat. Internally e_hat is Cartesian, after resolving the
    # user's crystallographic direction with the first image lattice.
    field_vector = field_MV_cm * result.field_direction
    # field_work_i = Delta_mu_i dot E_vec. This is the vector dot product in
    # e*Angstrom*MV/cm that enters -Delta_mu dot E.
    field_work = branch_delta_dipoles(result) @ field_vector
    # F_i(E_vec) = U_i - C_E * (Delta_mu_i dot E_vec).
    return rel - EV_PER_EANG_MV_CM * field_work


def barrier_indices(result: PathResult) -> range:
    """Return image indices included in the selected path criterion.

    这个函数只决定“哪些 image 或 step 参与判据”。它不计算能量。

    Meaning of modes:
        to-saddle: i = 1..saddle_image.
        intermediate: i = 1..N-2, excluding both endpoints.
        full-path: i = 1..N-1, excluding only image 00.
    """

    if result.barrier_mode == "to-saddle":
        return range(1, result.saddle_image + 1)
    if result.barrier_mode == "intermediate":
        return range(1, len(result.raw_images) - 1)
    if result.barrier_mode == "full-path":
        return range(1, len(result.raw_images))
    raise ValueError(f"unknown barrier mode: {result.barrier_mode}")


def path_value_at_field(result: PathResult, field_MV_cm: float) -> float:
    """Evaluate the active path criterion at one external field.

    这是给 barrier_vs_field.csv/png 使用的“路径判据函数”。它把一个电场
    E 映射成一个标量 path value。

    barrier 判据:
        B(E) = max_i F_i(E) - F_0(E)

    downhill 判据:
        D(E) = max_i [F_i(E) - F_{i-1}(E)]

    Ec 的定义是这个 path value 第一次降到 0 的电场。注意真正求 Ec
    用的是 coercive_requirements() 的解析式；本函数主要用于扫描和画图。

    Formulas:
        Barrier criterion:
            B(E) = max_i F_i(E) - F_0(E)

        Downhill criterion:
            D(E) = max_i [F_i(E) - F_{i-1}(E)]

    Meaning of symbols:
        i: image index included by barrier_indices(result).
        F_i(E): corrected_relative_energies(result, E)[i].
        B(E): whole-path barrier height in eV.
        D(E): largest adjacent uphill step in eV.
        Returned value: B(E) if criterion is "barrier", otherwise D(E).
    """

    corrected = corrected_relative_energies(result, field_MV_cm)
    if result.coercive_criterion == "barrier":
        indices = list(barrier_indices(result))
        if not indices:
            return 0.0
        # B(E) = max selected F_i(E) - F_0(E)。
        return float(np.max(corrected[indices]) - corrected[0])
    if result.coercive_criterion == "downhill":
        indices = list(barrier_indices(result))
        if not indices:
            return 0.0
        # D(E) = max selected adjacent uphill step。
        steps = [corrected[image] - corrected[image - 1] for image in indices]
        return float(np.max(steps))
    raise ValueError(f"unknown coercive criterion: {result.coercive_criterion}")


def barrier_at_field(result: PathResult, field_MV_cm: float) -> float:
    return path_value_at_field(result, field_MV_cm)


def coercive_requirements(
    rel_energies: np.ndarray,
    projected_delta: np.ndarray,
    candidate_indices: range,
    criterion: str,
) -> tuple[list[tuple[float, int]], list[int]]:
    """Compute the field required by each candidate constraint.

    这个函数把“求 Ec”转成一组线性不等式求解。

    每个约束的一般形式:
        Delta_E - 0.01 * E * Delta_mu_parallel <= 0

    其中 Delta_E 是零场下该约束两端的能量差，Delta_mu_parallel 是同一
    约束两端的投影偶极矩差。因为 F_i(E) 对 E 是线性的，所以每个约束
    都可以解析求出刚好变平的电场。

    若 Delta_E <= 0:
        该约束在零场下已经不是上坡，不限制 Ec。

    若 Delta_E > 0 且 Delta_mu_parallel > 0:
        外场能降低该上坡约束，所需电场为
        E_required = Delta_E / (0.01 * Delta_mu_parallel)

    若 Delta_E > 0 且 Delta_mu_parallel <= 0:
        当前外场方向不会降低该约束；电场越大也无法把它压到 0，
        因此记录到 impossible。

    Formulas:
        Barrier criterion for image i:
            Delta_E_i = U_i - U_0
            Delta_mu_parallel_i =
                (Delta_mu_i - Delta_mu_0) dot e_hat
            E_required_i =
                Delta_E_i / (0.01 * Delta_mu_parallel_i)

        Downhill criterion for step i-1 -> i:
            Delta_E_i = U_i - U_{i-1}
            Delta_mu_parallel_i =
                (Delta_mu_i - Delta_mu_{i-1}) dot e_hat
            E_required_i =
                Delta_E_i / (0.01 * Delta_mu_parallel_i)

    Meaning of code variables:
        rel_energies: array of U_i values in eV.
        projected_delta: array of Delta_mu_i dot e_hat in e*Angstrom.
        de: Delta_E_i, zero-field energy difference for the constraint.
        dm: Delta_mu_parallel_i, projected dipole difference for the constraint.
        requirements: finite E_required_i values in MV/cm.
        impossible: constraints with de > 0 and dm <= 0.
    """

    requirements: list[tuple[float, int]] = []
    impossible: list[int] = []
    if criterion == "barrier":
        for image in candidate_indices:
            # Barrier constraint compares image 00 -> image i:
            #   Delta_E = U_i - U_0
            #   Delta_mu_parallel = Delta_mu_i_parallel - Delta_mu_0_parallel
            # Since Delta_mu_0_parallel = 0, the subtraction is kept for clarity.
            de = float(rel_energies[image] - rel_energies[0])
            dm = float(projected_delta[image] - projected_delta[0])
            if de <= 0:
                # Already downhill/flat at zero field; it does not constrain Ec.
                continue
            if dm <= 0:
                # Positive uphill barrier but the selected field direction cannot lower it.
                impossible.append(image)
            else:
                # E_required = Delta_E / (0.01 * Delta_mu_parallel).
                requirements.append((de / (EV_PER_EANG_MV_CM * dm), image))
    elif criterion == "downhill":
        for image in candidate_indices:
            # Downhill constraint compares adjacent step image-1 -> image:
            #   Delta_E = U_i - U_{i-1}
            #   Delta_mu_parallel =
            #       Delta_mu_i_parallel - Delta_mu_{i-1}_parallel
            de = float(rel_energies[image] - rel_energies[image - 1])
            dm = float(projected_delta[image] - projected_delta[image - 1])
            if de <= 0:
                # This step is already not uphill at zero field.
                continue
            if dm <= 0:
                # The selected field direction cannot reduce this uphill step.
                impossible.append(image)
            else:
                # Field required to make this adjacent step exactly flat.
                requirements.append((de / (EV_PER_EANG_MV_CM * dm), image))
    else:
        raise ValueError(f"unknown coercive criterion: {criterion}")
    return requirements, impossible


def candidate_indices_for_mode(mode: str, image_count: int, saddle_image: int) -> range:
    """Map barrier-mode names to the image range checked by the criterion.

    这个函数必须和 barrier_indices() 保持一致。区别只是：
    barrier_indices() 用于扫描/画图；candidate_indices_for_mode() 用于
    解析求 Ec。
    """

    if mode == "to-saddle":
        return range(1, saddle_image + 1)
    if mode == "intermediate":
        return range(1, image_count - 1)
    if mode == "full-path":
        return range(1, image_count)
    raise ValueError(f"unknown barrier mode: {mode}")


def limiting_label(criterion: str, image: int) -> str:
    """Human-readable label for the limiting constraint."""

    if criterion == "downhill":
        return f"{image - 1}->{image}"
    return str(image)


def limiting_image_value(criterion: str, image: int) -> int:
    if criterion == "downhill":
        return image
    return image


def coercive_field_from_requirements(
    requirements: list[tuple[float, int]],
    impossible: list[int],
    criterion: str,
    mode: str,
    notes: list[str],
) -> tuple[float, int | None]:
    """Collapse per-constraint requirements into the final coercive field.

    coercive_requirements() 给出每个约束单独需要的电场。最终 Ec 必须让
    所有约束同时满足，所以取最大值。

    Formula:
        E_c = max_i E_required_i

    Meaning of symbols:
        E_required_i: field needed to flatten one barrier/downhill constraint.
        E_c: returned coercive field. The max is used because all constraints
            must be lowered at the same time.
        limiting image/step: image giving the maximum E_required_i.
    """

    if impossible:
        # 有任何无法被当前外场方向降低的正上坡约束，就不存在有限 Ec。
        notes.append(
            "no finite field for "
            + criterion
            + " constraints "
            + ",".join(limiting_label(criterion, image) for image in impossible)
            + f" in {mode} mode"
        )
        return math.inf, None
    if requirements:
        # 限制约束是 required field 最大的那一个。
        field, image = max(requirements, key=lambda item: item[0])
        return float(field), limiting_image_value(criterion, image)
    # 没有正上坡约束时，零场下路径已经满足所选判据，Ec = 0。
    return 0.0, None


def analyze_path(
    label: str,
    path: Path,
    energy_source: str,
    field_direction_arg: str,
    field_coordinates: str,
    branch_direction_arg: str,
    branch_coordinates: str,
    branch_mode: str,
    barrier_mode: str,
    coercive_criterion: str,
    initial_polarization: np.ndarray,
    final_polarization_arg: str,
    endpoint_direction_arg: str,
    endpoint_polarization_arg: str,
    continuity_threshold_uC_cm2: float,
    branch_search_radius: int,
) -> PathResult:
    """Analyze one sample path from raw VASP data through final Ec.

    这是单个 sample 的总控函数。它把输入文件中的 VASP 原始量转换成
    PathResult，后续 outputs.py 只负责把 PathResult 写出。

    High-level calculation flow:
        1. 读取并存储vasp原始数据
        2. 确定连续的branch
               P_branch_i = P_raw_i + n_i1*Q_a_i + n_i2*Q_b_i + n_i3*Q_c_i
           在center-zero模式下，计算流程为：
            首先将所有的原始极化值平移到同一条branch上，保证其连续性:
               Q_i^T c_i = P_branch_{i-1} - P_raw_i
               n_i = floor(c_i + 0.5)
               P_branch_i = P_raw_i + n_i @ Q_i
           在移动之后，将最中心的branch作为参考值，进行整体平移，计算每个image的branch修正后的极化值和偶极矩。
           对于非180度的终点极化路径，指定终点极化方向之后，将最中心的image到终点的四个image进行平移，使得最终极化方向沿着重点方向。
        3. 对于电场矫正，使用 
            Delta_mu_i = mu_branch_i - mu_branch_0
            其中 mu_branch_i 是经过 Berry branch 修正后的偶极矩，mu_branch_0 是第一个image的偶极矩。

        4. 对于电场下的矫正能量：
            E_vec = E * e_hat
            F_i(E_vec) = U_i - 0.01 * (Delta_mu_i dot E_vec)
            对于解析矫顽场Ec， 将 Delta_mu_i 投影到 e_hat 上，得到标量场强度 E 的未知系数。


    Main symbols:
        U_i: 零电场下image i的相对能量，单位 eV。
        Omega_i: image i的晶胞体积，单位 A^3。
        P_branch_i: 经过 Berry branch 修正后的image i的极化值。
        Delta_mu_i: 经过校正后的偶极矩变化（相对于 image 00），单位 e*Angstrom。
        e_hat: 外电场单位方向。
        E_c: 矫顽场 in MV/cm.
    """

    # 读取vasp计算出的原始数据，晶格参数，相对能量，晶胞体积，极化量子值。返回的raw_images是一个列表，每个元素都是一个RawImage对象。label用于指定名称，如“HfO2”，“a-HZO”、“b-HZO”、“c-HZO”; path用于指定每一个label对应的路径；energy用于指定能量从哪一个文件中读取。
    raw_images, notes = read_raw_images(label, path, energy_source)
    # 读取第一个图像的晶格参数--3x3的numpy数组，作为参考晶格，parse_structure函数返回的第二个值是晶胞体积，这里不需要，所以用 _ 忽略掉
    reference_lattice, _ = parse_structure(raw_images[0].structure_path)

    # 指定电场方向，默认是 [0, 0, -1]，
    # 这里的reverse-initial表示取初始极化方向的反方向来作为电场方向，
    # 但是在定义该参数的时候并没有说可以选择reverse-initial，所以需要后续修改！！
    if field_direction_arg == "reverse-initial":
        raw_field_direction = -initial_polarization
        # 将初始极化方向取反，得到外场方向，并归一化为单位向量， np.linalg.norm(raw_field_direction) 计算向量的欧几里得范数，即向量的长度。该等式为了得到单位向量。
        field_direction = raw_field_direction / np.linalg.norm(raw_field_direction)
    else:
        # 这里用的是第一个image的晶格参数作为参考晶格。field_direction_arg默认是0, 0, -1；field_coordinates默认是lattice，这里是为了考虑当
        field_direction = resolve_direction(field_direction_arg, field_coordinates, reference_lattice)

    # branch_direction的默认值为0,0,-1;
    # branch_coodinates默认是"lattice"，模式。
    if branch_direction_arg == "field":
        branch_direction = field_direction
    elif branch_direction_arg == "reverse-initial":
        raw_branch_direction = -initial_polarization
        branch_direction = raw_branch_direction / np.linalg.norm(raw_branch_direction)
    else:
        branch_direction = resolve_direction(branch_direction_arg, branch_coordinates, reference_lattice) # 这里使用branch的坐标参数依然使用第一个image的晶格参数作为参考晶格，仍有待商榷，个人认为是有问题的。

    # 将相对能赋值到一个列表中，并提取出能量值最高的image的索引值
    rel_energies = np.array([img.rel_energy_eV for img in raw_images])
    saddle_image = int(np.argmax(rel_energies))

    # final_polarization_arg为用户输入的最终极化方向-----需要想办法让它根据最后一个image自动推断出来，而不是用户输入。需要确定是最终极化方向还是最终极化值！目前写的final_polarization_arg的默认值为opposite-initial.   应该是最终极化值！！ 下面的endpoint_direction为最终极化方向。
    # initial_polarization为用户输入的初始极化方向-----可以由用户指定或者根据一个字典来自动判断，如Pca2_1: 001, Pmn2_1:101。同样需要确定是初始极化方向还是初始极化值！目前默认值为np.array([0.0, 0.0, 50.0])--为一个numpy数组   应该是最终极化值
    # resolve_final_polarization用于将用户输入的--final-polarization参数解析为一个numpy数组。这里之所以要传入initial_polarization，是因为对于final_polarization_arg为“opposite-initial”时，会返回-initial_polarization
    final_polarization = resolve_final_polarization(final_polarization_arg, initial_polarization)

    # ednpoint_direction_arg为用户输入的--endpoint-direction参数，默认值为auto。用于推断终点的极化方向。如果字段中含有Pca2_1，则。。。如果字段中含有Pmn2_1则。。。返回的是一个单位矢量
    endpoint_direction = resolve_endpoint_direction(endpoint_direction_arg, path)

    # 将用户输入的--endpoint-polarization参数解析为一个numpy数组，如果是auto或none，则返回None，否则返回一个三维向量
    endpoint_polarization = resolve_endpoint_polarization(endpoint_polarization_arg)
    if endpoint_direction_arg == "auto" and endpoint_direction is not None:
        notes.append(
            "auto endpoint branch direction "
            + ",".join(f"{value:.6g}" for value in endpoint_direction)
            + " used only to resolve final Berry branch"
        )

    # 为了后续计算矫顽场，进行build_branch_images()
    # 进行场修正之后的：
    #   P_branch_i = P_raw_i + n_i1*Q_a_i + n_i2*Q_b_i + n_i3*Q_c_i
    #
    # 默认的 center-zero 模式下，算法从 image 00 开始向最后一个 image 走。对于每个 image，它在极化量子基底中求解分支系数：
    #   Q_i^T c_i = P_branch_{i-1} - P_raw_i
    #   n_i = floor(c_i + 0.5)
    #   P_branch_i = P_raw_i + n_i @ Q_i
    # 创建一个连续的branch之后，通过将最中心的image的polarization尽可能接近于0，来实现全局的整数偏移。对于非180度的终点，终点方向/值可以用来解决最终分支的等效性。如何自动识别终点方向的值。
    branch_images = build_branch_images(
        raw_images=raw_images,
        branch_mode=branch_mode,
        branch_direction=branch_direction,
        center_image=len(raw_images) // 2,
        continuity_threshold_uC_cm2=continuity_threshold_uC_cm2,
        search_radius=branch_search_radius,
        initial_polarization=initial_polarization,
        final_polarization=final_polarization,
        endpoint_direction=endpoint_direction,
        endpoint_polarization=endpoint_polarization,
    )

    # Project branch-corrected dipole changes onto the external-field direction
    # only to get the linear coefficient used in the analytic Ec solve.
    # Math:
    #   Delta_mu_i = mu_branch_i - mu_branch_0
    #   projected_delta_i = Delta_mu_i dot e_hat
    #
    # The full corrected energy uses the vector dot product:
    #   F_i(E_vec) = U_i - 0.01 * (Delta_mu_i dot E_vec)
    # with E_vec = E * e_hat. Therefore projected_delta_i is just the
    # coefficient multiplying the unknown scalar field magnitude E.
    first_branch_dipole = branch_images[0].branch_dipole_eA
    projected_delta = np.array(
        [(image.branch_dipole_eA - first_branch_dipole) @ field_direction for image in branch_images],
        dtype=float,
    )
    # Select which images/steps become constraints in the chosen barrier mode.
    candidate_indices = candidate_indices_for_mode(barrier_mode, len(raw_images), saddle_image)
    # Solve every linear constraint for its required field.
    requirements, impossible = coercive_requirements(
        rel_energies,
        projected_delta,
        candidate_indices,
        coercive_criterion,
    )
    # Ec is the maximum finite requirement, unless any active constraint is impossible.
    coercive_field, limiting_image = coercive_field_from_requirements(
        requirements,
        impossible,
        coercive_criterion,
        barrier_mode,
        notes,
    )

    discontinuous = [img.image for img in branch_images if not img.component_continuity_ok]
    if discontinuous:
        notes.append(
            "branch component jump exceeds threshold at images "
            + ",".join(str(image) for image in discontinuous)
        )

    return PathResult(
        label=label,
        path=path,
        raw_images=raw_images,
        branch_images=branch_images,
        field_direction=field_direction,
        branch_direction=branch_direction,
        target_initial_polarization_uC_cm2=initial_polarization,
        target_final_polarization_uC_cm2=final_polarization,
        branch_mode=branch_mode,
        barrier_mode=barrier_mode,
        coercive_criterion=coercive_criterion,
        continuity_threshold_uC_cm2=continuity_threshold_uC_cm2,
        saddle_image=saddle_image,
        saddle_energy_eV=float(rel_energies[saddle_image]),
        coercive_field_MV_cm=coercive_field,
        limiting_image=limiting_image,
        endpoint_rel_energy_eV=float(rel_energies[-1]),
        endpoint_projected_dipole_eA=float(projected_delta[-1]),
        notes=notes,
    )
