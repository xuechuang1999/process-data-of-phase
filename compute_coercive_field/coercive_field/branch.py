from __future__ import annotations

import numpy as np
from itertools import product

from .constants import UC_CM2_PER_EANG_A3
from .models import BranchImage, RawImage


def direct_branch_to_reference(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    reference: np.ndarray,
    branch_search: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    # 计算raw_polarization到reference的最短路径，并返回调整后的极化值和对应的index
    coeff = np.linalg.solve(quantum_vectors.T, reference - raw_polarization)
    nearest_index = np.floor(coeff + 0.5).astype(int)

    best_pol = None
    best_index = None
    best_distance = np.inf

    # 在nearest_index的周围搜索，找到最接近reference的极化值
    for offset in product(range(-branch_search, branch_search + 1), repeat=3):
        offset = np.array(offset, dtype=int)
        index = nearest_index + offset

        candidate = raw_polarization + index @ quantum_vectors
        distance = np.linalg.norm(candidate - reference)

        if distance < best_distance:
            best_distance = distance
            best_pol = candidate
            best_index = index

    # print("=========================================")
    # print(f"raw_polarization: {raw_polarization},\n reference: {reference},\n coeff: {coeff},\n index: {nearest_index}")
    # print(f"raw_polarization + best_index @ quantum_vectors: {raw_polarization + best_index @ quantum_vectors}")
    # print(f"index:{best_index}")
    # print(f"quantum_vectors:\n {quantum_vectors}")
    # print("=========================================")

    return best_pol, best_index


def continuous_center_zero_branch(
    raw_polarizations: list[np.ndarray],
    quantum_vectors: list[np.ndarray],
    center_image: int,
    branch_search: int = 1,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    image_count = len(raw_polarizations)
    branch_pols: list[np.ndarray | None] = [None] * image_count
    branch_indices: list[np.ndarray | None] = [None] * image_count

    # 将中心的image的极化值调整为0附近的值，并得到调整之后的polarization值和对应的index
    branch_pols[center_image], branch_indices[center_image] = direct_branch_to_reference(
        raw_polarizations[center_image],
        quantum_vectors[center_image],
        np.zeros(3),
        branch_search=branch_search,
    )

    # 将中心image左边的极化值调整为连续的值
    for image in range(center_image - 1, -1, -1):
        assert branch_pols[image + 1] is not None
        print(f"Now is process the {image}th image")
        branch_pols[image], branch_indices[image] = direct_branch_to_reference(
            raw_polarizations[image],
            quantum_vectors[image],
            branch_pols[image + 1],
            branch_search=branch_search,
        )

    # 将中心image右边的极化值调整为连续的值
    for image in range(center_image + 1, image_count):
        assert branch_pols[image - 1] is not None
        print(f"Now is process the {image}th image")
        branch_pols[image], branch_indices[image] = direct_branch_to_reference(
            raw_polarizations[image],
            quantum_vectors[image],
            branch_pols[image - 1],
            branch_search=branch_search,
        )

    return [pol for pol in branch_pols if pol is not None], [index for index in branch_indices if index is not None]


def shift_right_half_to_endpoint(
    raw_polarizations: list[np.ndarray],
    quantum_vectors: list[np.ndarray],
    branch_pols: list[np.ndarray],
    branch_indices: list[np.ndarray],
    center_image: int,
    endpoint_reference: np.ndarray,
    continuity_threshold_uC_cm2: float,
) -> None:
    desired_final_pol, desired_final_index = direct_branch_to_reference(
        raw_polarizations[-1],
        quantum_vectors[-1],
        endpoint_reference,
        branch_search=1,
    )
    current_score = float(np.linalg.norm(branch_pols[-1] - endpoint_reference))
    desired_score = float(np.linalg.norm(desired_final_pol - endpoint_reference))
    offset = desired_final_index - branch_indices[-1]

    if desired_score >= current_score or not np.any(offset != 0):
        return # 未改变数据

    # 创建一个浅拷贝的列表
    proposed_pols = list(branch_pols)
    proposed_indices = list(branch_indices)
    for image in range(center_image + 1, len(branch_pols)):
        proposed_pols[image] = proposed_pols[image] + offset @ quantum_vectors[image]
        proposed_indices[image] = proposed_indices[image] + offset

    for image in range(center_image + 1, len(branch_pols)):
        jump = float(np.max(np.abs(proposed_pols[image] - proposed_pols[image - 1])))
        if jump > continuity_threshold_uC_cm2:
            return # 未改变数据

    # 将
    for image in range(center_image + 1, len(branch_pols)):
        branch_pols[image] = proposed_pols[image]
        branch_indices[image] = proposed_indices[image]


def shift_left_half_to_initial(
    raw_polarizations: list[np.ndarray],
    quantum_vectors: list[np.ndarray],
    branch_pols: list[np.ndarray],
    branch_indices: list[np.ndarray],
    center_image: int,
    initial_polarization_reference: np.ndarray,
    continuity_threshold_uC_cm2: float,
) -> None:
    desired_initial_pol, desired_initial_index = direct_branch_to_reference(
        raw_polarizations[0],
        quantum_vectors[0],
        initial_polarization_reference,
        branch_search=1,
    )
    current_score = float(np.linalg.norm(branch_pols[0] - initial_polarization_reference))
    desired_score = float(np.linalg.norm(desired_initial_pol - initial_polarization_reference))
    offset = desired_initial_index - branch_indices[0]

    if desired_score >= current_score or not np.any(offset != 0):
        return # 未改变数据

    # 创建一个浅拷贝的列表
    proposed_pols = list(branch_pols)
    proposed_indices = list(branch_indices)
    for image in range(0, center_image):
        proposed_pols[image] = proposed_pols[image] + offset @ quantum_vectors[image]
        proposed_indices[image] = proposed_indices[image] + offset
    
    for image in range(1, center_image + 1):
        jump = float(np.max(np.abs(proposed_pols[image] - proposed_pols[image - 1])))
        if jump > continuity_threshold_uC_cm2:
            return # 未改变数据
    
    for image in range(0, center_image):
        branch_pols[image] = proposed_pols[image]
        branch_indices[image] = proposed_indices[image]


def build_branch_images(
    initial_polarization: np.ndarray,
    raw_images: list[RawImage],
    projection_direction: np.ndarray,
    center_image: int,
    continuity_threshold_uC_cm2: float,
    endpoint_reference: np.ndarray | None,
) -> list[BranchImage]:
    
    #==============================================================
    # 提取原始极化矢量和极化量子
    #==============================================================
    raw_polarizations = [image.raw_polarization_uC_cm2 for image in raw_images]
    quantum_vectors = [image.quantum_vectors_uC_cm2 for image in raw_images]

    #==============================================================
    # 构建连续的分支极化矢量和对应的index
    #==============================================================
    branch_pols, branch_indices = continuous_center_zero_branch(
        raw_polarizations,
        quantum_vectors,
        center_image,
    )

    #=============================================================
    # 存在endpoint_reference时，将右半部分的极化矢量调整为接近endpoint_reference
    #=============================================================
    if endpoint_reference is not None:
        shift_right_half_to_endpoint(
            raw_polarizations,
            quantum_vectors,
            branch_pols,
            branch_indices,
            center_image,
            endpoint_reference,
            continuity_threshold_uC_cm2,
        )

    if initial_polarization is not None:
        shift_left_half_to_initial(
            raw_polarizations,
            quantum_vectors,
            branch_pols,
            branch_indices,
            center_image,
            initial_polarization,
            continuity_threshold_uC_cm2,
        )

    # ==============================================================
    # 构建BranchImage对象列表
    # ==============================================================
    branch_images: list[BranchImage] = []
    previous_pol: np.ndarray | None = None
    previous_projection: float | None = None

    for raw_image, branch_pol, branch_index in zip(raw_images, branch_pols, branch_indices):
        print(f"pol: {branch_pol}")

        # 根据极化矢量计算分支的偶极矩，单位为eÅ
        branch_dipole = branch_pol * raw_image.volume_A3 / UC_CM2_PER_EANG_A3

        # 计算极化矢量在外加电场方向上的投影,
        # projection_direction为笛卡尔坐标系下 外加电场方向的单位矢量
        projection = float(branch_pol @ projection_direction)
        print(f"projection: {projection}")

        if previous_pol is None:
            delta = np.zeros(3)
            delta_projection = 0.0
            max_jump = 0.0
            monotonic = True
            continuous = True
        else:
            delta = branch_pol - previous_pol
            delta_projection = projection - float(previous_projection)
            max_jump = float(np.max(np.abs(delta)))
            monotonic = delta_projection >= -1e-8
            continuous = max_jump <= continuity_threshold_uC_cm2
            # print("continuous: ", continuous)

        branch_images.append(
            BranchImage(
                image=raw_image.image,
                branch_index=branch_index,
                branch_polarization_uC_cm2=branch_pol,
                branch_dipole_eA=branch_dipole,
                branch_projection_uC_cm2=projection,
                delta_branch_polarization_from_prev_uC_cm2=delta,
                delta_branch_projection_from_prev_uC_cm2=delta_projection,
                max_abs_component_jump_from_prev_uC_cm2=max_jump,
                monotonic_along_branch=monotonic,
                component_continuity_ok=continuous,
            )
        )
        # 刷新previous_pol和previous_projection为当前的pol和projection
        # 以便在for循环的下一步中使用
        previous_pol = branch_pol
        previous_projection = projection

    return branch_images
