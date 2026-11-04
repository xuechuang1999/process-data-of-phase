from __future__ import annotations

import argparse

import numpy as np


def parse_vector(text: str) -> np.ndarray:
    parts = [float(part.strip()) for part in text.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated numbers")
    arr = np.array(parts, dtype=float)
    if np.linalg.norm(arr) == 0:
        raise argparse.ArgumentTypeError("vector must be nonzero")
    return arr


def resolve_direction(text: str, coordinates: str, lattice: np.ndarray) -> np.ndarray:
    # text: 方向的字符串表示，例如 "1,0,0" 或 "0.5,0.5,0.5"
    # coordinates: "cartesian" 或 "lattice"，表示输入的方向是晶轴坐标还是笛卡尔坐标
    # lattice 晶格参数
    vector = parse_vector(text)
    if coordinates == "cartesian":
        raw = vector
    elif coordinates == "lattice":
        raw = vector @ lattice
    else:
        raise ValueError(f"unknown coordinate system: {coordinates}")
    norm = float(np.linalg.norm(raw))
    if norm == 0:
        raise ValueError("direction has zero norm")
    return raw / norm


def lattice_to_cartesian(vector_lattice: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    晶轴坐标 -> 笛卡尔坐标

    lattice 的每一行是一个晶格矢量 a, b, c。
    """
    return vector_lattice @ lattice


def cartesian_to_lattice(vector_cartesian: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    笛卡尔坐标 -> 晶轴坐标

    满足：
        vector_cartesian = vector_lattice @ lattice
    """
    return np.linalg.solve(lattice.T, vector_cartesian)