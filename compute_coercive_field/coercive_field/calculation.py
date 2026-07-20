from __future__ import annotations

import itertools
import math
import re
from pathlib import Path

import numpy as np

from .branch import build_branch_images, continuous_center_zero_branch, direct_branch_to_reference
from .constants import EV_PER_EANG_MV_CM, UC_CM2_PER_EANG_A3
from .directions import resolve_direction, cartesian_to_lattice
from .final_polarization_direction import final_polarization_direction
from .models import PathResult, RawImage
from .vasp_io import parse_structure, read_raw_images


def parse_optional_vector(text: str) -> np.ndarray | None:
    text = text.strip()
    if text.lower() in {"auto", "vasp", "from-vasp"}:
        return None
    return np.array([float(part.strip()) for part in text.split(",")], dtype=float)


def parse_vector(text: str) -> np.ndarray:
    return np.array([float(part.strip()) for part in text.split(",")], dtype=float)


def phase_and_state_from_path(path: Path) -> tuple[str, str]:
    phase_names = {"pca21", "pca2_1", "pmn21", "pmn2_1", "fe1", "fe2"}
    parts = path.parts
    for index, part in enumerate(parts):
        if part.lower() not in phase_names:
            continue
        for candidate in parts[index + 1 :]:
            if re.search(r"[xyz]{3}-[01]{3}-[01]{3}", candidate):
                return part, candidate
        raise ValueError(f"could not find switching state after phase directory in {path}")
    raise ValueError(f"could not find phase directory in {path}")


def initial_lattice_direction_from_phase(phase: str) -> np.ndarray:
    phase = phase.lower()
    if phase in {"pca21", "pca2_1", "fe1"}:
        return np.array([0.0, 0.0, 1.0])
    if phase in {"pmn21", "pmn2_1", "fe2"}:
        return np.array([1.0, 0.0, 1.0])
    raise ValueError(f"unknown phase: {phase}")


def estimate_initial_magnitude(raw_images: list[RawImage], center_image: int) -> float:
    raw_polarizations = [image.raw_polarization_uC_cm2 for image in raw_images]
    quantum_vectors = [image.quantum_vectors_uC_cm2 for image in raw_images]
    branch_pols, _ = continuous_center_zero_branch(raw_polarizations, quantum_vectors, center_image)
    print(f"branch_pols: {branch_pols}")
    return float(np.linalg.norm(branch_pols[0]))


def calculate_initial_polarization(
    text: str,
    raw_images: list[RawImage],
    path: Path,
    reference_lattice: np.ndarray,
    center_image: int,
) -> tuple[np.ndarray, np.ndarray]:
    # 检测text中是否为"auto", "vasp", "from-vasp"其中之一，如果是则返回None，否则解析为向量
    explicit = parse_optional_vector(text)
    if explicit is not None:
        # 计算向量的模长
        norm = float(np.linalg.norm(explicit))
        # 得到单位矢量
        direction = explicit / norm if norm != 0.0 else np.zeros(3)
        return explicit, direction

    # phase_and_state_from_path函数返回的是相和末态表示符，phase为"pca21", "pca2_1", "pmn21", "pmn2_1", "fe1", "fe2"其中之一，_ 为xyz-101-000类似的表示符
    phase, _ = phase_and_state_from_path(path)
    # 根据相名获取初始极化方向的晶格坐标表示，pca2_1:[0.0, 0.0, 1.0]; pmn2_1:[1.0, 0.0, 1.0]
    initial_polarization_by_lattice_direction = initial_lattice_direction_from_phase(phase)

    # 将初态结构的极化朝向由晶格坐标转换为笛卡尔坐标，并归一化，得到初始化的方向。
    initial_polarization_by_cartesian_direction = resolve_direction(
        ",".join(str(value) for value in initial_polarization_by_lattice_direction),
        "lattice",
        reference_lattice,
    )

    # 另一个计算image0极化值的方法：1. 直接通过将晶格坐标表示的初始极化方向与初态结构的极化量子矢量相乘，
    # 得到初始极化近似值，然后通过direct_branch_to_reference函数将初始极化近似值调整为最接近的极化值
    # 因为极化量子是使用的一个3*3矩阵表示的，其是在笛卡尔坐标系下的表示
    initial_polarization_value_by_cartesian_direction, _ = direct_branch_to_reference(
        raw_images[0].raw_polarization_uC_cm2,
        raw_images[0].quantum_vectors_uC_cm2,
        initial_polarization_by_cartesian_direction @ raw_images[0].quantum_vectors_uC_cm2
        )
    print(f"initial_polarization_value_by_cartesian_direction----hhhhhh: {initial_polarization_value_by_cartesian_direction}")
    print(f"norm: {np.linalg.norm(initial_polarization_value_by_cartesian_direction)}")
    # # 转换为笛卡尔坐标表示的初始极化值：
    # initial_polarization_value_by_cartesian_direction = initial_polarization_value_by_lattice_direction @ reference_lattice


    # # 计算image0极化值的方法2: codex生成的代码
    # magnitude = estimate_initial_magnitude(raw_images, center_image)
    # target = magnitude * initial_polarization_by_cartesian_direction
    # initial_polarization_value_by_cartesian_direction, _ = direct_branch_to_reference(
    #     raw_images[0].raw_polarization_uC_cm2,
    #     raw_images[0].quantum_vectors_uC_cm2,
    #     target,
    # )
    return initial_polarization_value_by_cartesian_direction, initial_polarization_by_cartesian_direction


def resolve_final_polarization(text: str, initial_polarization: np.ndarray, path: Path, lattice_A: np.ndarray) -> np.ndarray:
    if text == "auto":
        phase, state = phase_and_state_from_path(path)
        _, rotation = final_polarization_direction(state, phase)

        # rotation为一个3*3的旋转矩阵，为晶轴坐标系下的旋转矩阵。
        # 将rotation矩阵转换为笛卡尔坐标系下的旋转矩阵
        # print(f"rotation: {rotation}")
        rotation = lattice_A.T @ rotation @ np.linalg.inv(lattice_A.T)
        final_polarization = rotation @ initial_polarization
        # print(f"rotation: {rotation}")
        # print(f"initial_polarization: {initial_polarization}")
        # print(f"path: {path}")
        # print(f"final_polarization: {final_polarization}")
        return final_polarization
    if text == "opposite-initial":
        return -initial_polarization
    return parse_vector(text)


def is_non_180_switch(initial_polarization: np.ndarray, final_polarization: np.ndarray) -> bool:
    initial_norm = float(np.linalg.norm(initial_polarization))
    final_norm = float(np.linalg.norm(final_polarization))
    if initial_norm == 0.0 or final_norm == 0.0:
        return False
    cosine = float((initial_polarization @ final_polarization) / (initial_norm * final_norm))
    return not np.isclose(cosine, -1.0, atol=1e-3)


def automatic_continuity_threshold(raw_images: list[RawImage]) -> float:
    max_quantum = 0.0
    for image in raw_images:
        for coeff in itertools.product((-1, 0, 1), repeat=3):
            if coeff == (0, 0, 0):
                continue
            # 计算沿着不同晶向的极化量子矢量的模长，并取最大值
            quantum = np.array(coeff, dtype=float) @ image.quantum_vectors_uC_cm2
            max_quantum = max(max_quantum, float(np.linalg.norm(quantum)))
    return max_quantum # 2 * max_quantum / len(raw_images)


def resolve_continuity_threshold(value: str | float, raw_images: list[RawImage]) -> tuple[float, bool]:
    if isinstance(value, str) and value.strip().lower() == "auto":
        return automatic_continuity_threshold(raw_images), True
    return float(value), False


def target_polarization(result: PathResult, image_index: int) -> np.ndarray:
    image_count = len(result.raw_images)
    progress = image_index / (image_count - 1) if image_count > 1 else 0.0
    return result.target_initial_polarization_uC_cm2 + progress * (
        result.target_final_polarization_uC_cm2 - result.target_initial_polarization_uC_cm2
    )


def switching_axis(result: PathResult) -> np.ndarray:
    axis = result.target_initial_polarization_uC_cm2
    norm = float(np.linalg.norm(axis))
    if norm == 0.0:
        return result.branch_direction
    return axis / norm


def target_delta_dipoles(result: PathResult) -> np.ndarray:
    initial_target = target_polarization(result, 0)
    rows = []
    for image_index, raw in enumerate(result.raw_images):
        delta_p = target_polarization(result, image_index) - initial_target
        rows.append(delta_p * raw.volume_A3 / UC_CM2_PER_EANG_A3)
    return np.array(rows, dtype=float)


def branch_delta_dipoles(result: PathResult) -> np.ndarray:
    first = result.branch_images[0].branch_dipole_eA
    return np.array([image.branch_dipole_eA - first for image in result.branch_images], dtype=float)


def branch_delta_polarizations(result: PathResult) -> np.ndarray:
    first = result.branch_images[0].branch_polarization_uC_cm2
    return np.array([image.branch_polarization_uC_cm2 - first for image in result.branch_images], dtype=float)


def corrected_relative_energies(result: PathResult, field_MV_cm: float) -> np.ndarray:
    rel = np.array([img.rel_energy_eV for img in result.raw_images])
    field_vector = field_MV_cm * result.field_direction
    field_work = branch_delta_dipoles(result) @ field_vector
    return rel - EV_PER_EANG_MV_CM * field_work


def barrier_indices(result: PathResult) -> range:
    return candidate_indices_for_mode(result.barrier_mode, len(result.raw_images), result.saddle_image)


def path_value_at_field(result: PathResult, field_MV_cm: float) -> float:
    corrected = corrected_relative_energies(result, field_MV_cm)
    indices = list(barrier_indices(result))
    if not indices:
        return 0.0
    if result.coercive_criterion == "barrier":
        return float(np.max(corrected[indices]) - corrected[0])
    if result.coercive_criterion == "downhill":
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
    requirements: list[tuple[float, int]] = []
    impossible: list[int] = []

    if criterion == "barrier":
        for image in candidate_indices:
            # 相对能量差值
            # 偶极矩差值
            delta_energy = float(rel_energies[image] - rel_energies[0])
            delta_dipole = float(projected_delta[image] - projected_delta[0])
            if delta_energy <= 0:
                continue
            if delta_dipole >= delta_energy:
                requirements.append((delta_energy / (delta_dipole), image))
            else:
                continue
        return requirements, impossible

    if criterion == "downhill":
        for image in candidate_indices:
            delta_energy = float(rel_energies[image] - rel_energies[image - 1])
            delta_dipole = float(projected_delta[image] - projected_delta[image - 1])
            if delta_energy <= 0:
                continue
            # if delta_dipole >= delta_energy:
            #     requirements.append((delta_energy / (delta_dipole), image))
            # else:
            #     continue

            requirements.append((delta_energy / (EV_PER_EANG_MV_CM * delta_dipole), image))
            print(f"image: {image}, delta_energy: {delta_energy}, delta_dipole: {delta_dipole}, requirement: {delta_energy / (EV_PER_EANG_MV_CM * delta_dipole)}")
        return requirements, impossible

    raise ValueError(f"unknown coercive criterion: {criterion}")


def candidate_indices_for_mode(mode: str, image_count: int, saddle_image: int) -> range:
    if mode == "to-saddle":
        return range(1, saddle_image + 1)
    if mode == "intermediate":
        return range(1, image_count - 1)
    if mode == "full-path":
        return range(1, image_count)
    raise ValueError(f"unknown barrier mode: {mode}")


def limiting_label(criterion: str, image: int) -> str:
    if criterion == "downhill":
        return f"{image - 1}->{image}"
    return str(image)


def limiting_image_value(criterion: str, image: int) -> int:
    return image


def coercive_field_from_requirements(
    requirements: list[tuple[float, int]],
    impossible: list[int],
    criterion: str,
    mode: str,
    notes: list[str],
) -> tuple[float, int | None]:
    if impossible:
        notes.append(
            "no finite field for "
            + criterion
            + " constraints "
            + ",".join(limiting_label(criterion, image) for image in impossible)
            + f" in {mode} mode"
        )
        return math.inf, None
    if requirements:
        field, image = max(requirements, key=lambda item: item[0])
        return float(field), limiting_image_value(criterion, image)
    return 0.0, None


# 计算F_i(E)=U_ks - \Omega \vec{P} \cdot \vec{E}，其中U_ks为vasp计算得到的能量，\Omega为体积，\vec{P}为极化矢量，\vec{E}为外加电场
# omega_i: 晶胞体积，单位为A^3
# P_i: 极化矢量，单位为uC/cm^2
# E_field_direction: 外加电场方向，单位为MV/cm
# E_field_magnitude: 外加电场强度，单位为MV/cm
# 返回值为F_i(E)，单位为eV
def omega_P_dot_E_energy(omega_i: float, P_i: np.ndarray, E_field_direction: np.ndarray, E_field_magnitude: float) -> float:
    print(f"omega_i: {omega_i}, P_i: {P_i}, E_field_direction: {E_field_direction}, E_field_magnitude: {E_field_magnitude}")
    omega_P_dot_E = omega_i * np.dot(P_i, E_field_direction * E_field_magnitude) / UC_CM2_PER_EANG_A3
    return omega_P_dot_E


def analyze_path(
    label: str,
    path: Path,
    energy_source: str,
    field_direction_arg: str,
    field_coordinates: str,
    barrier_mode: str,
    coercive_criterion: str,
    initial_polarization: str,
    final_polarization_arg: str,
    continuity_threshold_uC_cm2: str | float,
) -> PathResult:
    #=============================================================
    # 1. 从vasp的计算文件中提取原始数据。
    #    得到raw_images: list[RawImage]
    #=============================================================
    raw_images, notes = read_raw_images(label, path, energy_source)
    reference_lattice = np.array(raw_images[0].lattice_A)
    print(f"reference_lattice: {reference_lattice}")
    # print(f"raw_image[-1].image: {raw_images[-1].image}")

    center_image = len(raw_images) // 2

    #==============================================================
    # 2. 计算铁电极化连续时，相邻两个image的极化差值阈值
    # 默认的continuity_threshold为auto，continuity_threshold得到的为自动计算的阈值，而used_auto_threshold得到的为True
    #==============================================================
    continuity_threshold, used_auto_threshold = resolve_continuity_threshold(
        continuity_threshold_uC_cm2,
        raw_images,
    )
    # 如果是使用的自动计算的阈值，则打印输出计算的结果
    if used_auto_threshold:
        notes.append(
            f"auto branch continuity threshold {continuity_threshold:.6g} uC/cm^2 "
            "from 2* max(|integer-combined polarization quantum|)/(image_count)"
        )

    #==============================================================
    # 3. 推断初始极化
    # 首先calculate_initial_polarization函数中
    # 3.1. 检测initial_polarization是否为"auto", "vasp", "from-vasp"其中之一，如果是则返回None，否则解析为向量
    # 3.2. 如果返回了向量，则计算其模长和单位矢量，并返回
    # 3.3. 如果initial_polarization为"auto", "vasp", "from-vasp"其中之一，则解析path，检测path中是否包含相名（pca2_1, fe1）和末态标识符(xyz-001-000)
    # 3.4. 如果path中包含了相名，则会直接根据设定的相的初始极化方向，给出来初始极化方向，并根据传进去的reference_lattice将极化方向转换为笛卡尔坐标系下的单位矢量，至此，解决了初始极化方向的问题
    # 3.5. estimate_initial_magnitude()调用continuous_center_zero_branch()函数
    # initial_p 和initial_direction均使用的笛卡尔坐标系表示
    #==============================================================
    initial_p, initial_direction = calculate_initial_polarization(
        initial_polarization,
        raw_images,
        path,
        reference_lattice,
        center_image,
    )
    notes.append(
        "initial polarization "
        + ",".join(f"{value:.6g}" for value in initial_p)
        + " uC/cm^2 from direction "
        + ",".join(f"{value:.6g}" for value in initial_direction)
    )


    #==============================================================
    # 4. 计算末态极化值
    # final_polarization_arg为用户输入的--final-polarization对应的参数，可以为auto
    # final_polarization_arg为auto时，根据path中的相名和末态标识符，推断出最终极化方向，
    # 并将其转换为笛卡尔坐标系下的单位矢量
    #==============================================================
    final_p = resolve_final_polarization(final_polarization_arg, initial_p, path, raw_images[0].lattice_A)

    # # is_non_180_switch函数用于测试初始计划和最终极化方向是否为180度，如果为180度，则返回false
    # # 否则，返回true，表示为非180度的极化翻转，此时，endpoint_reference为最终极化方向final_p，否则为None
    # endpoint_reference = final_p if is_non_180_switch(initial_p, final_p) else None

    # endpoint_reference为最终极化方向final_p，直接替换上面的 `endpoint_reference = final_p if is_non_180_switch(initial_p, final_p) else None`
    endpoint_reference = final_p
    print(f"endpoint_reference: {endpoint_reference}")
    if endpoint_reference is not None:
        endpoint_reference, _ = direct_branch_to_reference(
            raw_images[-1].raw_polarization_uC_cm2,
            raw_images[-1].quantum_vectors_uC_cm2,
            endpoint_reference,
            )
        notes.append(
            "non-180 switch detected; final polarization "
            + ",".join(f"{value:.6g}" for value in final_p)
            + " uC/cm^2 used as endpoint branch reference"
        )

    #==============================================================
    # 5. 计算外加电场方向：笛卡尔坐标系下的方向
    # field_direction 默认为00-1方向，field-coordinates默认为lattice，
    # reference_lattice为image0的晶格矢量矩阵，
    # 转换为笛卡尔坐标系下的方向向量
    # field_direction得到笛卡尔坐标下的外加电场方向的单位矢量
    #==============================================================
    field_direction = resolve_direction(field_direction_arg, field_coordinates, reference_lattice)
    print(f"field_direction: {field_direction}")

    #==============================================================
    # 6. 构建branch
    #=============================================================
    branch_images = build_branch_images(
        initial_polarization=initial_p,
        raw_images=raw_images,
        projection_direction=field_direction,
        center_image=center_image,
        continuity_threshold_uC_cm2=continuity_threshold,
        endpoint_reference=endpoint_reference,
    )
    # print(f"The last image of the branch is {branch_images[-1].branch_polarization_uC_cm2} uC/cm^2")

    #=============================================================
    # 7. 计算外场下的影响
    #=============================================================
    rel_energies = np.array([img.rel_energy_eV for img in raw_images])
    #print(f"rel_energies: {rel_energies}")
    saddle_image = int(np.argmax(rel_energies))
    first_branch_dipole = branch_images[0].branch_dipole_eA

    # projected_delta为一个列表，表示每个分支图像的偶极矩在外加电场方向上的投影与第一个分支图像的偶极矩在外加电场方向上的投影的差值
    projected_delta = np.array(
        [(branch_image.branch_dipole_eA - first_branch_dipole) @ field_direction for branch_image in branch_images],
        dtype=float,
    )


    #=============================================================
    # 计算 \omega \vec{P} \cdot \vec{E} 能量项
    #=============================================================
    first_omega_P_E = omega_P_dot_E_energy(
        raw_images[0].volume_A3,
        branch_images[0].branch_polarization_uC_cm2,
        E_field_direction=field_direction,
        E_field_magnitude=1.0
        )
    delta_omega_P_E = np.array(
        [omega_P_dot_E_energy(raw_images[image].volume_A3, branch_images[image].branch_polarization_uC_cm2, E_field_direction=field_direction, E_field_magnitude=1.0)-first_omega_P_E for image in range(len(branch_images))],
        dtype=float,
    )
    print(f"delta_omega_P_E: {delta_omega_P_E}")
    # print("====="*6)
    # print(f"projected_delta: {projected_delta}")
    # print("====="*6)

    #=============================================================
    # 计算外场矫正下的能量
    # 所需参数：相对能，晶胞体积，branch 校正后的P，电场方向，电厂振幅
    #=============================================================
    # omega_P_E: list[np.ndarray | None] = [None] * len(raw_images)
    # field_magnitude = 1.0  # MV/cm，施加单位电场时
    # for image_index in range(len(branch_images)):
    #     omega_i = raw_images[image_index].volume_A3
    #     P_i = branch_images[image_index].branch_polarization_uC_cm2
    #     E_field_direction = field_direction
    #     E_field_magnitude = field_magnitude
    #     omega_P_E[image_index] = omega_P_dot_E_energy(omega_i, P_i, E_field_direction, E_field_magnitude)
    #     print(f"image_index: {image_index}, omega_P_E: {omega_P_E[image_index]} eV")


    discontinuous = [branch_image.image for branch_image in branch_images if not branch_image.component_continuity_ok]
    print(f"discontinuous: {discontinuous}")
    if discontinuous:
        print("========================================")
        notes.append(
            "branch component jump exceeds threshold at images "
            + ",".join(str(image) for image in discontinuous)
            + "; coercive field is not reliable until branch correction is fixed"
        )
        coercive_field = math.inf
        limiting_image = None
    else:
        # barrier_mode为用户输入--barrier-mode，默认为full-path
        #==============================================================
        # 根据输入的--barrier-mode参数，默认值为full-path，返回不同的范围值；
        # to-saddle: 1到saddle_image的范围
        # intermediate: 1到image_count-1的范围，意味着不考虑首尾两个端点的能量
        # full-path: 1到image_count的范围
        # 得到的是一个列表
        #==============================================================
        candidate_indices = candidate_indices_for_mode(barrier_mode, len(raw_images), saddle_image)
        # print("================================")
        # print(f"candidate_indices: {list(candidate_indices)}")

        #==============================================================
        # 
        #==============================================================
        requirements, impossible = coercive_requirements(
            rel_energies,
            delta_omega_P_E,
            candidate_indices,
            coercive_criterion,
        )
        print("================================")
        print(f"requirements: {requirements}")
        print(f"impossible: {impossible}")
        print("================================")

        #=============================================================
        # 根据requirements和impossible计算矫顽场和限制图像
        # requirements为一个列表，包含了当前image所需要的矫顽场大小
        # impossible为一个列表，包含了当前image无法满足矫顽场要求的图像索引
        # 
        #=============================================================
        coercive_field, limiting_image = coercive_field_from_requirements(
            requirements,
            impossible,
            coercive_criterion,
            barrier_mode,
            notes,
        )

    return PathResult(
        label=label,
        path=path,
        raw_images=raw_images,
        branch_images=branch_images,
        field_direction=field_direction,
        branch_direction=field_direction,
        target_initial_polarization_uC_cm2=initial_p,
        target_final_polarization_uC_cm2=final_p,
        branch_mode="center-zero",
        barrier_mode=barrier_mode,
        coercive_criterion=coercive_criterion,
        continuity_threshold_uC_cm2=continuity_threshold,
        saddle_image=saddle_image,
        saddle_energy_eV=float(rel_energies[saddle_image]),
        coercive_field_MV_cm=coercive_field,
        limiting_image=limiting_image,
        endpoint_rel_energy_eV=float(rel_energies[-1]),
        endpoint_projected_dipole_eA=float(projected_delta[-1]),
        notes=notes,
    )
