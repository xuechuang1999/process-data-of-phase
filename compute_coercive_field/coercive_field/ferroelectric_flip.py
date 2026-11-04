#!/usr/bin/env python3
"""Generate the opposite-polarization configuration from a VASP POSCAR.

Two construction modes are available.

1. inversion (one input structure)

   s_minus = 2 * origin - s_plus   (modulo lattice translations)

   This constructs the inversion-related domain and reverses all components of
   the polarization vector.  It is the appropriate default for a periodic,
   chemically homogeneous bulk ferroelectric in zero field.

2. reference (polar structure plus a non-polar/midpoint reference)

   delta   = minimum_image(s_plus - s_ref)
   s_minus = s_ref - delta

   This reflects the distortion about the reference.  Atom order must be the
   same in the two POSCAR files.  The reference must be a physically meaningful
   midpoint between the two domains; an arbitrary high-symmetry parent is not
   automatically such a midpoint when non-polar modes are also present.

The script writes Direct coordinates and preserves the input atom ordering and
Selective-dynamics flags.  It intentionally has only one dependency: NumPy.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .vasp_io import Poscar, read_poscar, write_poscar


def minimum_image(delta_frac: np.ndarray) -> np.ndarray:
    """Map fractional displacement components to [-0.5, 0.5)."""
    return delta_frac - np.floor(delta_frac + 0.5)


def check_reference_compatibility(initial: Poscar, reference: Poscar) -> None:
    if initial.counts != reference.counts or initial.symbols != reference.symbols:
        raise ValueError(
            "Initial and reference POSCAR headers differ. They must contain the "
            "same species, counts, and atom order."
        )
    if initial.natoms != reference.natoms:
        raise ValueError("Initial and reference structures have different atom counts")


def invert_fractional_positions(initial: Poscar, origin: np.ndarray) -> np.ndarray:
    """Return the inversion-related domain about a fractional origin."""

    center = np.asarray(origin, dtype=float)
    if center.shape != (3,):
        raise ValueError(f"origin has shape {center.shape}, expected (3,)")
    return np.mod(2.0 * center - initial.frac, 1.0)


def reflect_about_reference(
    initial: Poscar,
    reference: Poscar,
    *,
    cell_tol: float = 1.0e-5,
    allow_cell_mismatch: bool = False,
) -> np.ndarray:
    """Reverse the initial distortion about an atom-matched reference."""

    check_reference_compatibility(initial, reference)
    mismatch = float(np.max(np.abs(initial.lattice - reference.lattice)))
    if mismatch > cell_tol and not allow_cell_mismatch:
        raise ValueError(
            f"Cell mismatch is {mismatch:.6g} A, larger than cell tolerance "
            f"{cell_tol:g} A. Use structures in the same cell, or explicitly "
            "allow a cell mismatch after checking the mapping."
        )

    displacement = minimum_image(initial.frac - reference.frac)
    return np.mod(reference.frac - displacement, 1.0)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Construct an opposite-polarization ferroelectric POSCAR."
    )
    parser.add_argument("input", type=Path, help="initial +P POSCAR/CONTCAR")
    parser.add_argument("output", type=Path, help="output -P POSCAR")
    sub = parser.add_subparsers(dest="mode", required=True)

    inv = sub.add_parser(
        "inversion",
        help="invert the entire structure; requires only the initial structure",
    )
    inv.add_argument(
        "--origin",
        nargs=3,
        type=float,
        default=(0.0, 0.0, 0.0),
        metavar=("SX", "SY", "SZ"),
        help="fractional inversion center (default: 0 0 0)",
    )

    ref = sub.add_parser(
        "reference",
        help="reflect atomic displacements about a supplied midpoint structure",
    )
    ref.add_argument("--reference", type=Path, required=True, help="midpoint POSCAR")
    ref.add_argument(
        "--cell-tol",
        type=float,
        default=1.0e-5,
        help="maximum absolute lattice mismatch in angstrom (default: 1e-5)",
    )
    ref.add_argument(
        "--allow-cell-mismatch",
        action="store_true",
        help="use fractional coordinates even if the two cells differ",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    initial = read_poscar(args.input)

    if args.mode == "inversion":
        origin = np.asarray(args.origin, dtype=float)
        flipped = invert_fractional_positions(initial, origin)
        method = f"inversion about fractional origin {origin.tolist()}"
    else:
        reference = read_poscar(args.reference)
        flipped = reflect_about_reference(
            initial,
            reference,
            cell_tol=args.cell_tol,
            allow_cell_mismatch=args.allow_cell_mismatch,
        )
        method = f"displacement reflection about {args.reference}"

    comment = f"{initial.comment} | opposite polarization: {method}"
    write_poscar(args.output, initial, flipped, comment)

    print(f"Wrote {args.output}")
    print(f"Atoms: {initial.natoms}")
    print(f"Method: {method}")
    print("Note: verify the result with Berry-phase polarization and relaxation.")


if __name__ == "__main__":
    main()
