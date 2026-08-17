#!/usr/bin/env python3
# conda activate icet
import os
import argparse
import warnings
from typing import Dict, List, Tuple

from ase import Atoms
from ase.io import read, write
from ase.build import sort
from icet import ClusterSpace
from icet.tools.structure_generation import generate_sqs_from_supercells
from icet.input_output.logging_tools import set_log_config

# 配置日志
set_log_config(level="INFO")


def load_structure(filename: str) -> Atoms:
    """读取结构文件并清理不必要的残留信息（如占据率）。"""
    # 忽略 ASE 解析 CIF 时的警告
    warnings.filterwarnings("ignore", category=UserWarning, module="ase.io.cif")

    if filename.endswith('.cif'):
        atoms = read(filename, format='cif')
    else:
        atoms = read(filename)

    # 清理残留的占据率和对称位点标记，防止后续生成 CIF 时报错
    if 'occupancy' in atoms.info:
        del atoms.info['occupancy']
    if 'spacegroup_kinds' in atoms.arrays:
        del atoms.arrays['spacegroup_kinds']

    print(f"Loaded FILE: {filename}")
    print(f"Chemical symbol: {atoms.get_chemical_formula()}")
    return atoms

def cal_cut_off(atoms: Atoms, supercell_matrix: Tuple[int, int, int]) -> List[float]:
    supercell = atoms * supercell_matrix
    a, b, c = supercell.cell.cellpar()[:3]
    min_lattice = min(a, b, c)

    # PBC limit
    pbc_limit = min_lattice / 2.0

    # physical limit
    phys_limit_pair = 8.0
    phys_limit_triplet = 5.0
    phys_limit_quadtuplet = 4.0

    # get the min value between PBC limit and phys limit
    pair_cutoff = min(pbc_limit-0.1, phys_limit_pair)
    triplet_cutoff = min(pair_cutoff * 1.0, phys_limit_triplet)
    quadtuplet_cutoff = min(pair_cutoff * 1.0, phys_limit_quadtuplet)

    # get round
    # cut_off = round(min_lattice / 2 - 0.5)

    return [pair_cutoff, triplet_cutoff, quadtuplet_cutoff]



def setup_cluster_space(atoms: Atoms, cutoffs: List[float], substitutions: Dict[str, List[str]]) -> ClusterSpace:
    """根据替换规则设置 ClusterSpace。"""
    chemical_symbols = []
    for atom in atoms:
        # 如果当前元素在我们的替换字典中，则使用字典中的列表；否则保持原样
        if atom.symbol in substitutions:
            chemical_symbols.append(substitutions[atom.symbol])
        else:
            chemical_symbols.append([atom.symbol])

    cs = ClusterSpace(atoms, cutoffs, chemical_symbols)
    print("ClusterSpace successfully generated.")
    return cs


def build_sqs_structure(
    cs: ClusterSpace,
    atoms: Atoms,
    supercell_matrix: Tuple[int, int, int],
    target_concentrations: Dict[str, float],
    n_steps: int
) -> Atoms:
    """使用 ASE 扩胞并在固定的超胞形状内寻找最佳 SQS 结构。"""
    # 1. 生成指定的超胞
    supercell = atoms * supercell_matrix
    print(f"Supercell generated with {len(supercell)} atoms.")

    # 2. 生成 SQS
    sqs = generate_sqs_from_supercells(
        cluster_space=cs,
        supercells=[supercell],
        target_concentrations=target_concentrations,
        n_steps=n_steps,
        random_seed=42
    )

    print("Cluster vector of generated SQS:\n ", cs.get_cluster_vector(sqs))
    return sqs


def save_structure(sqs: Atoms, original_filename: str, order: List[str] = None):
    """排序原子顺序并将其保存为 VASP 和 CIF 格式。"""
    prefix = os.path.splitext(original_filename)[0]
    output_vasp = f"{prefix}_sqs_sorted.vasp"
    output_cif = f"{prefix}_sqs_sorted.cif"

    # 按化学元素排序
    if order is not None:
        element_order = {elem: i for i, elem in enumerate(order)}
        #sqs_sorted = sort(sqs, key=lambda x: element_order.get(x.symbol, 999))
        indices = sorted(range(len(sqs)), key=lambda i: element_order.get(sqs[i].symbol, 999))
        sqs_sorted = sqs[indices]
    else:
        sqs_sorted = sort(sqs)

    # 写入文件 (加入 vasp5=True 以显式输出元素符号)
    write(output_vasp, sqs_sorted, format="vasp", direct=True, vasp5=True)
    write(output_cif, sqs_sorted, format='cif')
    print(f"Has saved the SQS structure files as:\n - {output_vasp}\n - {output_cif}")


def print_cluster_vector_details(cs: ClusterSpace, sqs: Atoms):
    """打印 cluster vector 每一项对应的 orbit 信息。"""
    cv = cs.get_cluster_vector(sqs)

    print("\n===== Cluster vector details =====")
    print(f"{'idx':>4} {'order':>5} {'kind':>12} {'radius':>10} "
          f"{'mult':>6} {'sublattices':>12} {'mc_vector':>20} {'value':>14}")

    for orbit, value in zip(cs.as_list, cv):
        order = orbit["order"]

        if order == 0:
            kind = "zerolet"
        elif order == 1:
            kind = "singlet"
        elif order == 2:
            kind = "pair"
        elif order == 3:
            kind = "triplet"
        elif order == 4:
            kind = "quadruplet"
        else:
            kind = f"{order}-body"

        idx = orbit["index"]
        radius = orbit["radius"]
        mult = orbit["multiplicity"]
        sublattices = orbit.get("sublattices", ".")
        mc_vector = orbit.get("multi_component_vector", ".")

        print(f"{idx:4d} {order:5d} {kind:>12} {radius:10.4f} "
              f"{mult:6d} {str(sublattices):>12} {str(mc_vector):>20} {value:14.8f}")


def main():
    """主控函数，串联整个工作流。"""
    parser = argparse.ArgumentParser(description="Generate SQS structure from a given unit cell.")
    parser.add_argument("filename", help="The filename you want to read (e.g. POSCAR or cif)")
    parser.add_argument("--order", nargs="+",  help="Specify element order, e.g. --order Zr Hf O")
    parser.add_argument("--steps", type=int, default=10000, help="Set the steps number of MC")
    parser.add_argument("-sc_m", "--supercell_matrix", nargs=3, type=int, help="Create supercell along `a` `b` `c` direction")
    args = parser.parse_args()

    # ================= 配置参数区 =================
    # 将以前散落在代码各处的参数集中在这里，方便以后修改
    # CUTOFFS = [7.0, 5.0]
    SUBSTITUTIONS = {"Zr": ["Zr", "Hf"], "O": ["O"]}           # 替换规则：将 Zr 替换为 Zr 和 Hf
    TARGET_CONCENTRATIONS = {"A": {"Zr": 0.5, "Hf": 0.5}, "B": {"O": 1}}  # 目标浓度
    # SUPERCELL_MATRIX = (3, 3, 3)                   # 扩胞矩阵
    N_STEPS = 10000                                 # 蒙特卡洛步数
    # ==============================================

    # 1. 加载和清理结构
    atoms = load_structure(args.filename)

    # print(atoms)
    # cell_matrix = atoms.cell[:]
    # super_cell_matrix = cell_matrix * args.supercell_matrix
    # print(super_cell_matrix)

    # get the cutoff
    auto_cutoffs = cal_cut_off(atoms, args.supercell_matrix)
    print(f"The cutoff for pair, triplet, quadruplet are: {auto_cutoffs}")

    # 2. 设置 ClusterSpace
    # cs = setup_cluster_space(atoms, CUTOFFS, SUBSTITUTIONS)
    cs = setup_cluster_space(atoms, auto_cutoffs, SUBSTITUTIONS)
    print(cs)

    # 3. 扩胞并生成 SQS
    # sqs = build_sqs_structure(cs, atoms, SUPERCELL_MATRIX, TARGET_CONCENTRATIONS, args.steps)
    sqs = build_sqs_structure(cs, atoms, args.supercell_matrix, TARGET_CONCENTRATIONS, args.steps)
    
    print_cluster_vector_details(cs, sqs)

    # 4. 排序并保存
    save_structure(sqs, args.filename, args.order)


if __name__ == "__main__":
    main()
