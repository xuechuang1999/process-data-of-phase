"""Direction parsing and normalization."""

from __future__ import annotations

import argparse

import numpy as np


def parse_vector(text: str) -> np.ndarray:
    """Parse a command-line vector written as ``x,y,z``.
    将一个文本字符串解析为一个三维向量，文本格式为 "x,y,z"，其中 x、y、z 是浮点数，表示向量的三个分量。
    """

    parts = [float(part.strip()) for part in text.split(",")] # strip() 用于去除每个分量的前后空格，split(",") 用于将字符串按逗号分割为三个部分
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated numbers")
    arr = np.array(parts, dtype=float)
    if np.linalg.norm(arr) == 0:
        raise argparse.ArgumentTypeError("vector must be nonzero")
    return arr # 返回一个 numpy 数组，表示解析后的三维向量


def resolve_direction(text: str, coordinates: str, lattice: np.ndarray) -> np.ndarray:
    """Convert the user direction into a Cartesian unit vector.

    coordinates 有三种含义：
    - cartesian：用户输入已经是笛卡尔方向。
    - lattice：用户输入是晶轴(a,b,c)坐标，即完整晶格基矢的线性组合
      d1*a+d2*b+d3*c。例如 -1,0,0 表示沿 -a 方向。
    - lattice-diagonal：旧 Excel 兼容模式，只使用晶格矩阵对角元素，
      即按 a_x、b_y、c_z 三个长度近似归一化；非正交晶胞不应默认使用。

    Formulas:
        cartesian:
            r = d
        lattice:
            r = d1*a + d2*b + d3*c
        lattice-diagonal:
            r = (d1*a_x, d2*b_y, d3*c_z)
        e_hat = r / |r|

    Meaning of symbols:
        d = (d1,d2,d3): user vector from --field-direction or
            --branch-direction.
        a,b,c: lattice vectors from the first image structure.
        a_x,b_y,c_z: diagonal entries of the lattice matrix.
        r: unnormalized Cartesian direction.
        e_hat: returned Cartesian unit vector.
    """
    # text 为用户输入的--field-direction 或者 --branch-direction 参数的值，表示一个方向向量，格式为 "x,y,z"。对于电场默认值为 "0,0,-1"。
    # coordinates 为用户输入的坐标系类型，可能是 "cartesian"、 "lattice"、 "lattice-diagonal"---该项建议弃用，为了兼容旧版本的excel表格，但是其实计算过程中引入计算误差。对于电场默认为 "lattice"。
    # lattice 为第一个图像的晶格向量矩阵，是一个 3x3 的 numpy 数组，表示晶格的三个基矢。

    vector = parse_vector(text) # 将text解析为一个三维向量，返回一个 numpy 数组
    if coordinates == "cartesian":
        raw = vector
    elif coordinates == "lattice":
        raw = vector @ lattice # 这里使用矩阵乘法，将用户输入的晶轴坐标向量与晶格矩阵相乘，得到笛卡尔坐标下的方向向量，但是对于电场方向来说，这么做是否合适有待商榷。-----同时，这里的晶格矩阵是第一个image的晶格矩阵，也需要重新确定。以及没有进行归一化！！
    elif coordinates == "lattice-diagonal":
        raw = np.array(
            [
                vector[0] * lattice[0, 0],
                vector[1] * lattice[1, 1],
                vector[2] * lattice[2, 2],
            ],
            dtype=float,
        )
    else:
        raise ValueError(f"unknown coordinate system: {coordinates}")
    norm = float(np.linalg.norm(raw))
    if norm == 0:
        raise ValueError("direction has zero norm")
    return raw / norm
