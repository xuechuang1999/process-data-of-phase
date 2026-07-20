from pathlib import Path

import numpy as np
import re


def string2matrix(string):
    prefix, suffix, translation = string.split("-")

    M = np.diag([-1 if i == "1" else 1 for i in suffix])

    P = np.zeros((3, 3), dtype=int)
    for index, i in enumerate([ord(i) - 119 for i in prefix]):
        P[index, i - 1] = 1

    translation = [0.5 if i == "1" else 0 for i in translation]
    return M @ P, translation


def extract_state_name(name):
    name = Path(name).name
    m = re.search(r"([xyz]{3}-[01]{3}-[01]{3})", name)
    if not m:
        raise ValueError(f"Cannot find state name in: {name}")
    return m.group(1)


def final_polarization_direction(folder_name, phase):
    state = extract_state_name(folder_name)
    rot, _ = string2matrix(state)

    phase = phase.lower()

    if phase in ["pca21", "pca2_1", "fe1"]:
        p0 = np.array([0, 0, 1], dtype=int)
    elif phase in ["pmn21", "pmn2_1", "fe2"]:
        p0 = np.array([1, 0, 1], dtype=int)
    else:
        raise ValueError(f"Unknown phase: {phase}")

    p_final = rot @ p0

    return p_final, rot
