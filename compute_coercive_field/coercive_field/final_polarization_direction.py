import re
import numpy as np
from pathlib import Path


def string2matrix(string):
    prefix, suffix, translation = string.split("-")

    M = np.diag([-1 if i == "1" else 1 for i in suffix])

    P = np.zeros((3, 3), dtype=int)
    for index, i in enumerate([ord(i) - 119 for i in prefix]):
        P[index, i - 1] = 1

    translation = [0.5 if i == "1" else 0 for i in translation]
    return M @ P, translation


def extract_state_name(name):
    """
    从 xyz-000-011 或 a-xyz-000-011 中提取 xyz-000-011
    """
    name = Path(name).name
    m = re.search(r"([xyz]{3}-[01]{3}-[01]{3})", name)
    if not m:
        raise ValueError(f"Cannot find state name in: {name}")
    return m.group(1)


def direction_to_label(v):
    """
    把 [0, 0, -1] 转成 [0 0 -1]
    把 [1, 0, 1] 转成 [1 0 1]
    """
    v = np.asarray(v, dtype=int)
    return "[" + " ".join(str(i) for i in v) + "]"


def final_polarization_direction(folder_name, phase):
    """
    phase:
        'pca21'  初始极化方向 001
        'pmn21'  初始极化方向 101
        'fe1'    初始极化方向 001, 也是  'pca21'
        'fe2'    初始极化方向 101, 也是  'pmn21'
    """
    state = extract_state_name(folder_name)
    rot, trans = string2matrix(state)

    phase = phase.lower()
    if phase in ["pca21", "pca2_1", "fe1"]:
        p0 = np.array([0, 0, 1], dtype=int)
    elif phase in ["pmn21", "pmn2_1", "fe2"]:
        p0 = np.array([1, 0, 1], dtype=int)
    else:
        raise ValueError(f"Unknown phase: {phase}")

    p_final = rot @ p0

    return {
        "folder": folder_name,
        "state": state,
        "rotation_matrix": rot,
        "translation": trans,
        "initial_direction": direction_to_label(p0),
        "final_direction_vector": p_final,
        "final_direction": direction_to_label(p_final),
    }

# # print(final_polarization_direction("xzy-001-000", "pca2_1"))
# for key, value in final_polarization_direction("xzy-001-000", "pmn2_1").items():
#     print(f"{key}: {value}")

#     # if key != "folder":
#     #     print(f"{key}: {value}")