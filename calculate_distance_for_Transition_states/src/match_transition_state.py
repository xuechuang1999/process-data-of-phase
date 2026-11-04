#!/usr/bin/env python3
"""Find the relaxed (or symmetry-equivalent) structure closest to a VASP TS.

The script scans the ``n_*`` directories below one or more phase roots.  For a
pristine structure it reads the structures already generated below
``reduced/<operator>``.  For ``a-n_*``, ``b-n_*`` and ``c-n_*`` structures it
applies the same operator to the doped relaxed structure in memory.  Crucially,
this transformation never reorders the doped atoms.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms
from ase.io import read

from reduced import align_cell, string2matrix


STATE_DIR_RE = re.compile(r"^\d+_.+")
DOPED_PREFIXES = ("a-", "b-", "c-")


@dataclass
class Candidate:
    phase: str
    state: str
    doping: str
    operator: str
    source: Path
    atoms: Atoms


def structure_file(directory: Path) -> Path | None:
    """Prefer the relaxed CONTCAR, falling back to POSCAR."""
    for name in ("CONTCAR", "POSCAR"):
        path = directory / name
        if path.is_file() and path.stat().st_size:
            return path
    return None


def transform_doped(atoms: Atoms, operator: str) -> Atoms:
    """Apply reduced.py's representative operation without sorting atoms."""
    rotation, translation = string2matrix(operator)
    # reduced.py represents signed axis permutations.  Using its exact formulas
    # keeps this result one-to-one with the corresponding pristine state.
    permutation = [ord(char) - 120 for char in operator[:3]]
    new_cell = np.asarray(atoms.cell)[:, permutation][permutation, :]
    for axis in range(3):
        if operator[axis + 4] == "1":
            new_cell[axis, :] *= -1
            new_cell[:, axis] *= -1
    new_cell = align_cell(new_cell)

    scaled = atoms.get_scaled_positions(wrap=True)
    new_scaled = np.transpose(
        [
            1.0 - scaled[:, permutation[axis]]
            if operator[axis + 4] == "1"
            else scaled[:, permutation[axis]]
            for axis in range(3)
        ]
    )
    # string2matrix is called as a validation too; directory names that are not
    # valid reduced.py operators should fail clearly rather than be misread.
    del rotation
    new_scaled = (new_scaled + np.asarray(translation)) % 1.0

    result = atoms.copy()
    result.set_cell(new_cell, scale_atoms=False)
    result.set_scaled_positions(new_scaled)
    return result


def operators_for(state_dir: Path) -> list[str]:
    reduced_dir = state_dir / "reduced"
    if not reduced_dir.is_dir():
        return []
    return sorted(
        child.name
        for child in reduced_dir.iterdir()
        if child.is_dir() and structure_file(child) is not None
    )


def collect_candidates(roots: Iterable[Path]) -> list[Candidate]:
    candidates: list[Candidate] = []
    for root in roots:
        phase = root.parent.name if root.name == "states" else root.name
        if not root.is_dir():
            raise FileNotFoundError(f"phase root does not exist: {root}")
        for state_dir in sorted(root.iterdir()):
            if not state_dir.is_dir() or not STATE_DIR_RE.fullmatch(state_dir.name):
                continue
            operators = operators_for(state_dir)
            for operator in operators:
                source = structure_file(state_dir / "reduced" / operator)
                if source:
                    candidates.append(
                        Candidate(
                            phase, state_dir.name, "pristine", operator,
                            source, read(source, format="vasp"),
                        )
                    )
            for prefix in DOPED_PREFIXES:
                doped_dir = state_dir / f"{prefix}{state_dir.name}"
                source = structure_file(doped_dir)
                if source is None:
                    continue
                doped = read(source, format="vasp")
                for operator in operators:
                    candidates.append(
                        Candidate(
                            phase, state_dir.name, prefix[0], operator, source,
                            transform_doped(doped, operator),
                        )
                    )
    return candidates


def _linear_sum(cost: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError as exc:
        raise RuntimeError(
            "scipy is required for unambiguous atom assignment; install scipy"
        ) from exc
    return linear_sum_assignment(cost)


def _assignment(
    reference: Atoms, candidate: Atoms, shift: np.ndarray
) -> tuple[float, float]:
    """Return Cartesian RMSD and maximum displacement for one origin shift."""
    ref_frac = reference.get_scaled_positions(wrap=True)
    cand_frac = candidate.get_scaled_positions(wrap=True)
    ref_symbols = np.asarray(reference.get_chemical_symbols())
    cand_symbols = np.asarray(candidate.get_chemical_symbols())
    # Measure both structures in the mean cell.  This includes positional
    # mismatch without allowing a slightly relaxed cell to dominate the score.
    metric_cell = 0.5 * (np.asarray(reference.cell) + np.asarray(candidate.cell))
    squared: list[float] = []
    for symbol in sorted(set(ref_symbols)):
        ri = np.flatnonzero(ref_symbols == symbol)
        ci = np.flatnonzero(cand_symbols == symbol)
        if len(ri) != len(ci):
            return float("inf"), float("inf")
        delta = ref_frac[ri, None, :] - (cand_frac[ci, :] + shift)[None, :, :]
        delta -= np.rint(delta)
        cart = delta @ metric_cell
        cost = np.einsum("ijk,ijk->ij", cart, cart)
        rows, cols = _linear_sum(cost)
        squared.extend(cost[rows, cols].tolist())
    values = np.asarray(squared)
    return float(np.sqrt(values.mean())), float(np.sqrt(values.max()))


def compare(reference: Atoms, candidate: Atoms) -> tuple[float, float, float]:
    """Compare structures under PBC, atom permutation, and origin translation."""
    if len(reference) != len(candidate):
        return float("inf"), float("inf"), float("inf")
    if sorted(reference.get_chemical_symbols()) != sorted(candidate.get_chemical_symbols()):
        return float("inf"), float("inf"), float("inf")

    ref_frac = reference.get_scaled_positions(wrap=True)
    cand_frac = candidate.get_scaled_positions(wrap=True)
    ref_symbols = np.asarray(reference.get_chemical_symbols())
    cand_symbols = np.asarray(candidate.get_chemical_symbols())
    # Any equivalent origin shift maps at least one atom onto an atom of the
    # same element.  Enumerating those anchor shifts is robust for small cells.
    rare_symbol = min(set(ref_symbols), key=lambda s: np.count_nonzero(ref_symbols == s))
    shifts = (
        ref_frac[ref_symbols == rare_symbol, None, :]
        - cand_frac[cand_symbols == rare_symbol][None, :, :]
    ).reshape(-1, 3)
    rmsd, max_displacement = min(
        (_assignment(reference, candidate, shift) for shift in shifts),
        key=lambda pair: pair[0],
    )
    ref_lengths = np.asarray(reference.cell.lengths())
    cand_lengths = np.asarray(candidate.cell.lengths())
    cell_strain = float(np.linalg.norm((cand_lengths - ref_lengths) / ref_lengths))
    return rmsd, max_displacement, cell_strain


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "transition_state", nargs="+", type=Path,
        help="one or more TS POSCAR/CONTCAR files",
    )
    parser.add_argument(
        "--root", action="append", required=True, type=Path,
        help="phase states root; repeat for PE and FE",
    )
    parser.add_argument("--top", type=int, default=20, help="number of matches")
    parser.add_argument("--csv", type=Path, help="write all compatible matches")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidates = collect_candidates(args.root)
    all_rows = []
    displayed = []
    for ts_path in args.transition_state:
        ts = read(ts_path, format="vasp")
        rows = []
        for item in candidates:
            rmsd, maximum, strain = compare(ts, item.atoms)
            if np.isfinite(rmsd):
                rows.append(
                    {
                        "transition_state": str(ts_path),
                        "rmsd_A": rmsd,
                        "max_displacement_A": maximum,
                        "cell_length_strain": strain,
                        "phase": item.phase,
                        "state": item.state,
                        "doping": item.doping,
                        "operator": item.operator,
                        "source": str(item.source),
                    }
                )
        rows.sort(key=lambda row: (row["rmsd_A"], row["cell_length_strain"]))
        all_rows.extend(rows)
        displayed.extend(rows[: args.top])
    if not all_rows:
        raise SystemExit("No composition-compatible candidate structures found.")

    fields = list(all_rows[0])
    print(",".join(fields))
    for row in displayed:
        print(",".join(str(row[field]) for field in fields))
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(all_rows)


if __name__ == "__main__":
    main()
