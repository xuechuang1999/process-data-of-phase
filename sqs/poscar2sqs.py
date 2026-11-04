from __future__ import annotations
import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

@dataclass
class Poscar:
    comment: str
    scale: float
    lattice: List[List[float]]
    elements: List[str]
    counts: List[int]
    coord_type: str
    coords_frac: List[List[float]]
    site_elements: List[str]

def _mat_det3(m: List[List[float]]) -> float:
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)

def _mat_inv3(m: List[List[float]]) -> List[List[float]]:
    det = _mat_det3(m)
    if abs(det) < 1e-16:
        raise ValueError("Lattice matrix is singular (det ~ 0). Cannot invert.")
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    adj = [
        [(e*i - f*h), -(b*i - c*h), (b*f - c*e)],
        [-(d*i - f*g), (a*i - c*g), -(a*f - c*d)],
        [(d*h - e*g), -(a*h - b*g), (a*e - b*d)]
    ]
    inv = [[adj[r][c] / det for c in range(3)] for r in range(3)]
    return inv

def _mat_vec_mul(m: List[List[float]], v: List[float]) -> List[float]:
    return [
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]
    ]

def _vec_add(a: List[float], b: List[float]) -> List[float]:
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]