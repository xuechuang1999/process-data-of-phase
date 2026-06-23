"""Direction parsing and normalization."""

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
    vector = parse_vector(text)
    if coordinates == "cartesian":
        raw = vector
    elif coordinates == "lattice":
        raw = vector @ lattice
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

