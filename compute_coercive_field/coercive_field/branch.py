"""Berry-branch correction logic."""

from __future__ import annotations

import numpy as np

from .constants import UC_CM2_PER_EANG_A3
from .models import BranchImage, RawImage


def direct_branch_to_reference(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Choose the integer Berry branch directly from quantum coordinates.

    VASP 的 raw polarization 只确定到 polarization quantum 的整数倍。
    对一个参考极化 reference，先在 polarization quantum 基底中求解
    raw_polarization 到 reference 还需要平移多少个 quantum：

        coeff = solve(Q^T, reference - P_raw)

    这里 Q 的三行是 Q_a、Q_b、Q_c。对于非正交晶胞，不能直接用
    Cartesian 分量逐项相除；必须在 quantum vector 基底中求系数。

    然后用 floor(coeff + 0.5) 直接得到最接近 reference 的整数分支。
    这等价于“除以极化量子后取整数”，但参考点是 previous/target，
    而不是固定的 Cartesian 原点。这样可以保持 NEB 路径连续。

    如果要把单个极化值折回 [0,Q) 主区间，才使用 floor(coeff) 形式；
    本工作流需要的是相邻 image 连续，所以使用最近整数。

    Formula:
        Q^T c = P_ref - P_raw
        n = floor(c + 0.5)
        P_branch = P_raw + n1*Q_a + n2*Q_b + n3*Q_c

    Meaning of symbols:
        P_raw: raw_polarization, LCALCPOL raw polarization in uC/cm^2.
        Q_a, Q_b, Q_c: rows of quantum_vectors, polarization quantum vectors.
        P_ref: reference, previous image or target reference polarization.
        c: continuous branch coefficient in the quantum basis.
        n: returned integer Berry-branch translation.
    """

    coeff = np.linalg.solve(quantum_vectors.T, reference - raw_polarization)
    index = np.floor(coeff + 0.5).astype(int)
    return raw_polarization + index @ quantum_vectors, index


def fold_to_target(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """选择一个最接近目标极化的整数branch。"""

    return direct_branch_to_reference(raw_polarization, quantum_vectors, target)


def fold_near_previous(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    previous_polarization: np.ndarray,
    branch_direction: np.ndarray,
    side: str,
    search_radius: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Choose a branch that stays continuous with the previous image.

    This helper is retained for compatibility with older branch experiments.
    It now uses direct quantum-coordinate rounding instead of enumerating
    neighboring integer branches.

    Formula:
        Q^T c = P_previous - P_raw
        n = floor(c + 0.5)
        P_branch = P_raw + n @ Q

    Meaning of symbols:
        P_previous: previous_polarization, already selected branch at neighbor.
        Q: quantum_vectors, rows are Q_a,Q_b,Q_c.
        n: integer branch translation.
    """

    del branch_direction, side, search_radius
    return direct_branch_to_reference(raw_polarization, quantum_vectors, previous_polarization)


def fold_continuous_from_previous(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    previous_polarization: np.ndarray,
    search_radius: int,
) -> tuple[np.ndarray, np.ndarray]:
    """选择一个整数branch，使其与前一个image的极化值最连续。

    This is the branch correction used by the default center-zero workflow:
    no target polarization is supplied. We only use the polarization quantum
    vectors and choose integer translations that keep adjacent images on the
    same continuous branch.

    The integer branch is computed directly, not by brute-force enumeration.
    In the quantum-vector basis:

        Q^T c = P_previous - P_raw
        n = floor(c + 0.5)
        P_branch = P_raw + n @ Q

    The floor(c + 0.5) step chooses the nearest integer quantum translation
    to the previous corrected image, which keeps the path continuous.

    Meaning of symbols:
        P_raw: raw LCALCPOL polarization of the current image.
        Q_a,Q_b,Q_c: polarization quantum vectors of the current image.
        P_previous: already corrected polarization of the previous image.
        n1,n2,n3: integer branch translations from floor(c + 0.5).
    """

    # 当前代码不再使用 search_radius 参数，所以删除了该局部变量
    del search_radius

    return direct_branch_to_reference(raw_polarization, quantum_vectors, previous_polarization)


def fold_continuous_near_target(
    raw_polarization: np.ndarray,
    quantum_vectors: np.ndarray,
    previous_polarization: np.ndarray,
    target: np.ndarray,
    search_radius: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Choose a branch near the previous image, with a weak target preference.

    endpoint-target mode follows the path from image 00 to the end. The
    default continuity reference is the previous image. A very weak target
    pull is retained by moving the reference 1% from previous_polarization
    toward target before applying direct quantum-coordinate rounding.

    Formulas:
        P_ref = P_previous + 0.01*(P_target - P_previous)
        Q^T c = P_ref - P_raw
        n = floor(c + 0.5)
        P_branch = P_raw + n @ Q

    Meaning of symbols:
        P_previous: selected branch polarization of the previous image.
        P_target: weak target polarization for the current image.
        P_ref: direct-rounding reference, mostly previous with a weak target pull.
    """

    del search_radius
    weak_reference = previous_polarization + 0.01 * (target - previous_polarization)
    return direct_branch_to_reference(raw_polarization, quantum_vectors, weak_reference)


def endpoint_reference(
    initial_polarization: np.ndarray,
    endpoint_direction: np.ndarray | None,
    endpoint_polarization: np.ndarray | None,
) -> np.ndarray | None:
    """Build the final-image branch reference used to resolve endpoint degeneracy.

    This is not a target path through all images. It only tells the branch
    selector which physically equivalent final-image branch should be preferred
    when several continuous branches are possible.

    Formulas:
        if P_endpoint is given:
            P_ref = P_endpoint
        else:
            P_ref = |P_initial| * e_endpoint

    Meaning of symbols:
        P_initial: branch-corrected polarization of image 00.
        e_endpoint: unit vector for the expected non-180 endpoint direction.
        P_ref: final-image polarization reference in uC/cm^2.
    """

    if endpoint_polarization is not None:
        return endpoint_polarization
    if endpoint_direction is None:
        return None
    # 这里norm返回的是initial_polarization的模长，如果为0，则返回None，否则返回模长乘以endpoint_direction，表示最终的参考极化方向。但是这样计算是有问题的！！如果endpoint_direction为（1，0，1）会放大原来的极化值！！
    norm = float(np.linalg.norm(initial_polarization))
    if norm == 0:
        return None
    return norm * endpoint_direction


def refine_forward_branch_to_endpoint(
    raw_pols: list[np.ndarray],
    quanta: list[np.ndarray],
    branch_pols: list[np.ndarray | None],
    branch_indices: list[np.ndarray | None],
    center_image: int,
    endpoint_target: np.ndarray,
    search_radius: int,
    continuity_threshold_uC_cm2: float,
) -> None:
    """Select a continuous right-half branch that ends near the endpoint target.

    The default continuity-only branch can be ambiguous for non-180 switching:
    a closed-looking branch may be equivalent by continuity but end on a less
    physical final branch. This function no longer enumerates branch
    candidates. It computes the final-image branch directly from the endpoint
    target, converts that into one common integer offset for the right half of
    the path, and accepts the offset only if every adjacent component jump
    remains below the continuity threshold.

    Formulas:
        jump(i) = max_alpha |P_i - P_{i-1}|_alpha
        n_endpoint = floor(c_endpoint + 0.5)
        offset = n_endpoint - n_current_endpoint
        P_i' = P_i + offset @ Q_i

    Meaning of symbols:
        P_endpoint: endpoint_target, expected final branch polarization.
        P_N: candidate branch polarization at the final image.
        alpha: Cartesian component x, y, or z.
    """

    image_count = len(raw_pols)
    if center_image >= image_count - 1:
        raise ValueError("Center image must be before the final image for endpoint refinement.")

    assert branch_pols[center_image] is not None
    assert branch_indices[center_image] is not None
    assert branch_pols[-1] is not None
    assert branch_indices[-1] is not None

    # 通过fold_to_target函数计算出最终的极化值和最终的分支索引
    desired_final_pol, desired_final_index = fold_to_target(raw_pols[-1], quanta[-1], endpoint_target)
    # 计算当前极化值与目标极化值之间的距离
    current_endpoint_score = float(np.linalg.norm(branch_pols[-1] - endpoint_target))
    # 计算目标极化值与最终极化值之间的距离
    desired_endpoint_score = float(np.linalg.norm(desired_final_pol - endpoint_target))
    # 计算当前分支索引与目标分支索引之间的偏移量
    offset = desired_final_index - branch_indices[-1]

    # 如果目标极化值比当前极化值更接近终点极化值，并且偏移量不为零，则尝试将偏移量应用于右半部分的分支索引和分支极化值
    if desired_endpoint_score < current_endpoint_score and np.any(offset != 0):
        # 将branch_pols和branch_indices复制到proposed_pols和proposed_indices中，以便在不破坏原始数据的情况下进行修改
        proposed_pols = list(branch_pols)
        proposed_indices = list(branch_indices)
        for image in range(center_image + 1, image_count):
            assert proposed_pols[image] is not None and proposed_indices[image] is not None
            proposed_pols[image] = proposed_pols[image] + offset @ quanta[image]
            proposed_indices[image] = proposed_indices[image] + offset

        continuous = True
        for image in range(center_image + 1, image_count):
            assert proposed_pols[image] is not None and proposed_pols[image - 1] is not None
            delta = proposed_pols[image] - proposed_pols[image - 1]
            if float(np.max(np.abs(delta))) > continuity_threshold_uC_cm2:
                continuous = False
                break
        if continuous:
            for image in range(center_image + 1, image_count):
                branch_pols[image] = proposed_pols[image]
                branch_indices[image] = proposed_indices[image]
            return # 这里也没有返回值啊？？！！！

    # If the direct endpoint offset is not better, or it breaks continuity,
    # leave the continuity-selected branch unchanged.
    del search_radius


def build_branch_images(
    raw_images: list[RawImage],
    branch_mode: str,
    branch_direction: np.ndarray,
    center_image: int,
    continuity_threshold_uC_cm2: float,
    search_radius: int,
    initial_polarization: np.ndarray,
    final_polarization: np.ndarray,
    endpoint_direction: np.ndarray | None = None,
    endpoint_polarization: np.ndarray | None = None,
) -> list[BranchImage]:
    """Build the continuous Berry-branch data for an entire NEB path.

    对于每一个branch i 的极化值定义:
        P_branch_i = P_raw_i
            + n_i1*Q_a_i + n_i2*Q_b_i + n_i3*Q_c_i
        mu_branch_i = P_branch_i * Omega_i / UC_CM2_PER_EANG_A3
        mu_branch_i为 branch-corrected dipole，用于电场修正计算。？？待确定

    极化P的连续性诊断：
        1. 计算相邻两个image的极化差值：
        Delta_P_prev_i = P_branch_i - P_branch_{i-1}
        2. 计算相邻两个image的极化差值的最大分量：
        max_jump_i = max_alpha |Delta_P_prev_{i,alpha}|
        3. 判断是否满足连续性阈值：
        component_continuity_ok = max_jump_i <= continuity_threshold_uC_cm2

    默认的中心极化值最接近0的工作流：
        1. 通过仅使用相邻image的连续性和极化量子平移，将image连接成一个连续的分支。
        2. 对每个image应用相同的整数平移k=(k1,k2,k3)：
            P_branch_i = P_branch_i + k1*Q_a_i + k2*Q_b_i + k3*Q_c_i
        其中k矢量的选择是使得中间的image最接近0。


    符号意义：
        P_raw_i: raw_images[i].raw_polarization_uC_cm2. 为第i个image的原始极化值，单位为uC/cm^2。
        Q_a_i, Q_b_i, Q_c_i: raw_images[i].quantum_vectors_uC_cm2. 为第i个image的极化量子向量，单位为uC/cm^2。
        n_i1,n_i2,n_i3: 返回第i个image的branch indices。
        Omega_i: raw_images[i].volume_A3. 为第i个image的晶格体积，单位为Angstrom^3。
        mu_branch_i: branch_dipole_eA, 经过branch-correction的偶极矩，单位为e*Angstrom，用于电场修正计算。
    """

    # 提取每个image的原始极化值和极化量子
    raw_pols = [img.raw_polarization_uC_cm2 for img in raw_images]
    quanta = [img.quantum_vectors_uC_cm2 for img in raw_images]
    image_count = len(raw_images)

    branch_pols: list[np.ndarray | None] = [None] * image_count
    branch_indices: list[np.ndarray | None] = [None] * image_count

    if branch_mode == "center-zero":
        # 首先将每个image放在一个连续的branch上，不使用任何目标极化值；只添加极化量子的整数倍以最小化相邻跳跃。
        # 读取第一个image的极化值作为该branch的起点，以此为起点，将其余的image依次放在这个branch上。所以branch_pols[0]和branch_indices[0]分别为raw_pols[0]和全零的整数分支。
        branch_pols[0] = raw_pols[0]
        branch_indices[0] = np.zeros(3, dtype=int)
        for image in range(1, image_count):
            # 确保前一个image的branch_pols不为None，因为后面需要减去前一个image的极化值来计算两者差值
            assert branch_pols[image - 1] is not None
            branch_pols[image], branch_indices[image] = fold_continuous_from_previous(
                raw_pols[image],
                quanta[image],
                branch_pols[image - 1],
                search_radius,
            )

        # 判断中心的image的branch_pols是否为None，如果是None，则抛出异常，因为后续需要使用该值来计算全局平移
        assert branch_pols[center_image] is not None
        # 调用fold_to_target函数，使得中心的image尽可能地接近于0，将返回的全局平移的整数分支赋值到global_shift变量中
        _, global_shift = fold_to_target(branch_pols[center_image], quanta[center_image], np.zeros(3))
        # 将所有的image平移global_shift个极化量子。
        for image in range(image_count):
            assert branch_pols[image] is not None and branch_indices[image] is not None
            branch_pols[image] = branch_pols[image] + global_shift @ quanta[image]
            branch_indices[image] = branch_indices[image] + global_shift

        # Non-180 endpoints can have several branches that are all continuous
        # within the threshold. Use only the final endpoint direction/value to
        # choose among those alternatives; do not force a linear target path.
        # branch_pols[0] 不是None，
        assert branch_pols[0] is not None

        # 根据输入的初始极化值，终点方向和终点极化值，得到终点极化值：
        # 1. 如果终点极化值（endpoint_polarization）不为None，则直接返回终点极化值
        # 2. 如果终点方向不为None，则返回初始极化值的模长乘以终点方向向量。如何实现根据xxx-nnn-nnn的标识自动计算终点极化值。
        # 3. 如果branch_pols[0]的模长为0，则返回None，否则返回模长乘以终点方向向量。
        endpoint_target = endpoint_reference(branch_pols[0], endpoint_direction, endpoint_polarization)
        # 
        if endpoint_target is not None:
            refine_forward_branch_to_endpoint(
                raw_pols,
                quanta,
                branch_pols,
                branch_indices,
                center_image,
                endpoint_target,
                search_radius,
                continuity_threshold_uC_cm2,
            )
    elif branch_mode == "endpoint-target":
        # 已弃用，需要是已知initial_polarization和final_polarization,根据两者差值进行计算得到中间的polarization的目标值，随后调用fold_continuous_near_target函数进行计算得到每个image的branch_pols和branch_indices
        initial_target = initial_polarization
        branch_pols[0], branch_indices[0] = fold_to_target(raw_pols[0], quanta[0], initial_target)
        for image in range(1, image_count):
            assert branch_pols[image - 1] is not None
            progress = image / (image_count - 1) if image_count > 1 else 0.0
            target = initial_polarization + progress * (final_polarization - initial_polarization)
            branch_pols[image], branch_indices[image] = fold_continuous_near_target(
                raw_pols[image],
                quanta[image],
                branch_pols[image - 1],
                target,
                search_radius,
            )
    else:
        raise ValueError(f"unknown branch mode: {branch_mode}")

    branch_images: list[BranchImage] = []
    previous_pol: np.ndarray | None = None
    previous_projection: float | None = None
    for raw, pol, index in zip(raw_images, branch_pols, branch_indices):
        assert pol is not None and index is not None
        # 将经过平移后的polarization转换为相应的偶极矩，方便后续计算电场下的NEB能量路径
        branch_dipole = pol * raw.volume_A3 / UC_CM2_PER_EANG_A3
        # 计算当前image的polarization在branch_direction方向上的投影，得到的是
        projection = float(pol @ branch_direction)
        if previous_pol is None:
            delta = np.zeros(3)
            delta_projection = 0.0
            max_jump = 0.0
            monotonic = True
            continuous = True
        else:
            delta = pol - previous_pol
            delta_projection = projection - float(previous_projection)
            max_jump = float(np.max(np.abs(delta)))
            monotonic = delta_projection >= -1e-8
            continuous = max_jump <= continuity_threshold_uC_cm2
        branch_images.append(
            BranchImage(
                image=raw.image,
                branch_index=index,
                branch_polarization_uC_cm2=pol,
                branch_dipole_eA=branch_dipole,
                branch_projection_uC_cm2=projection,
                delta_branch_polarization_from_prev_uC_cm2=delta,
                delta_branch_projection_from_prev_uC_cm2=delta_projection,
                max_abs_component_jump_from_prev_uC_cm2=max_jump,
                monotonic_along_branch=monotonic,
                component_continuity_ok=continuous,
            )
        )
        previous_pol = pol
        previous_projection = projection
    return branch_images
