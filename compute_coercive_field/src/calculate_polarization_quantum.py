#!/usr/bin/env python3
import sys
import numpy as np
from pathlib import Path

E_PER_A2_TO_UC_CM2 = 1602.176634  # e/Å^2 -> μC/cm^2
E_PER_A2_TO_C_M2 = 16.02176634    # e/Å^2 -> C/m^2


def read_lattice_from_poscar(path: str | Path) -> np.ndarray:
    path = Path(path)
    lines = path.read_text().splitlines()

    scale_values = [float(x) for x in lines[1].split()]
    lattice_raw = np.array(
        [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)],
        dtype=float,
    )

    if len(scale_values) == 1:
        scale = scale_values[0]

        # VASP: negative scale means target cell volume
        if scale < 0:
            target_volume = abs(scale)
            raw_volume = abs(np.linalg.det(lattice_raw))
            scale = (target_volume / raw_volume) ** (1.0 / 3.0)

        lattice = lattice_raw * scale

    elif len(scale_values) == 3:
        # Rare case: three scale factors
        lattice = lattice_raw * np.array(scale_values)

    else:
        raise ValueError("POSCAR/CONTCAR 第 2 行 scale 格式无法识别。")

    return lattice

def calculate_cell_volume(lattice: np.ndarray) -> float:
    """Calculate the volume of a unit cell given its lattice vectors.

    Args:
        lattice (np.ndarray): A 3x3 array representing the lattice vectors.

    Returns:
        float: The volume of the unit cell.
    """
    return abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))

def calculate_polarization_quantum(lattice: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    volume = abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))

    # rows: Qa, Qb, Qc; columns: x, y, z
    q_e_per_a2 = lattice / volume
    q_uc_cm2 = q_e_per_a2 * E_PER_A2_TO_UC_CM2
    q_c_m2 = q_e_per_a2 * E_PER_A2_TO_C_M2

    return q_e_per_a2, q_uc_cm2, q_c_m2


def main():
    if len(sys.argv) < 2:
        print("Usage: python calc_pol_quantum.py CONTCAR")
        sys.exit(1)

    contcar = sys.argv[1]
    lattice = read_lattice_from_poscar(contcar)

    a, b, c = lattice
    # volume = abs(np.dot(a, np.cross(b, c)))

    # # rows: Qa, Qb, Qc; columns: x, y, z
    # q_e_per_a2 = lattice / volume
    # q_uc_cm2 = q_e_per_a2 * E_PER_A2_TO_UC_CM2
    # q_c_m2 = q_e_per_a2 * E_PER_A2_TO_C_M2

    volume = calculate_cell_volume(lattice)

    q_e_per_a2, q_uc_cm2, q_c_m2 = calculate_polarization_quantum(lattice)

    np.set_printoptions(precision=8, suppress=False)

    print(f"File: {contcar}")
    print()
    print("Lattice vectors in Å:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print()
    print(f"Cell volume Ω = {volume:.10f} Å^3")
    print()
    print("Polarization quantum basis vectors:")
    print("Rows are Qa, Qb, Qc; columns are x, y, z.")
    print()
    print("In e/Å^2:")
    print(q_e_per_a2)
    print()
    print("In C/m^2:")
    print(q_c_m2)
    print()
    print("In μC/cm^2:")
    print(q_uc_cm2)
    print()
    print("Meaning:")
    print("P' = P + n1*Qa + n2*Qb + n3*Qc, where n1,n2,n3 are integers.")


if __name__ == "__main__":
    main()