#!/usr/bin/env python3
"""Report which relaxed phase/state is closest to every VASP NEB image."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
from ase.io import read

from match_transition_state import Candidate, collect_candidates, compare, structure_file


NEB_DIR_RE = re.compile(r"^[xyz]{3}-[01]{3}-[01]{3}$")
IMAGE_DIR_RE = re.compile(r"^\d+$")
DOPED_PREFIXES = ("a-", "b-", "c-")

DEFAULT_NEB_ROOTS = (
    Path("/gs/home/xuechuang/work/ZrO2/hzo/phase-transition/fe1"),
    Path("/gs/home/xuechuang/work/ZrO2/hzo/phase-transition/fe2"),
)
DEFAULT_PHASE_ROOTS = (
    Path("/gs/home/xuechuang/work/ZrO2/hafnia/1-PE/states"),
    Path("/gs/home/xuechuang/work/ZrO2/hafnia/2-FE/states"),
)


@dataclass(frozen=True)
class NebJob:
    group: str
    path_operator: str
    doping: str
    directory: Path


def collect_neb_jobs(roots: Iterable[Path]) -> list[NebJob]:
    """Find pristine and a/b/c NEB calculations in deterministic order."""
    jobs: list[NebJob] = []
    for root in roots:
        if not root.is_dir():
            raise FileNotFoundError(f"NEB root does not exist: {root}")
        for neb_dir in sorted(root.iterdir()):
            if not neb_dir.is_dir() or not NEB_DIR_RE.fullmatch(neb_dir.name):
                continue
            jobs.append(NebJob(root.name, neb_dir.name, "pristine", neb_dir))
            for prefix in DOPED_PREFIXES:
                doped = neb_dir / f"{prefix}{neb_dir.name}"
                if doped.is_dir():
                    jobs.append(NebJob(root.name, neb_dir.name, prefix[0], doped))
    return jobs


def collect_images(job: NebJob) -> Iterator[tuple[str, Path]]:
    """Yield numeric VASP image directories, including zero-padded names."""
    image_dirs = sorted(
        (
            child
            for child in job.directory.iterdir()
            if child.is_dir() and IMAGE_DIR_RE.fullmatch(child.name)
        ),
        key=lambda path: (int(path.name), path.name),
    )
    for image_dir in image_dirs:
        source = structure_file(image_dir)
        if source is not None:
            yield image_dir.name, source


def match_job(
    job: NebJob, candidates: Iterable[Candidate], top: int
) -> list[dict[str, object]]:
    """Return up to ``top`` matches for each image in one NEB calculation."""
    compatible_doping = [
        item for item in candidates if item.doping == job.doping
    ]
    rows: list[dict[str, object]] = []
    for image, image_source in collect_images(job):
        atoms = read(image_source, format="vasp")
        ranked = []
        for item in compatible_doping:
            rmsd, maximum, strain = compare(atoms, item.atoms)
            if np.isfinite(rmsd):
                ranked.append((rmsd, strain, maximum, item))
        ranked.sort(key=lambda result: (result[0], result[1]))
        for rank, (rmsd, strain, maximum, item) in enumerate(ranked[:top], 1):
            rows.append(
                {
                    "neb_group": job.group,
                    "path_operator": job.path_operator,
                    "doping": job.doping,
                    "image": image,
                    "image_source": str(image_source),
                    "rank": rank,
                    "closest_phase": item.phase,
                    "closest_state": item.state,
                    "closest_operator": item.operator,
                    "rmsd_A": rmsd,
                    "max_displacement_A": maximum,
                    "cell_length_strain": strain,
                    "candidate_source": str(item.source),
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--neb-root", action="append", type=Path,
        help="NEB root; repeat for fe1 and fe2 (defaults to cluster paths)",
    )
    parser.add_argument(
        "--phase-root", action="append", type=Path,
        help="relaxed phase states root; repeat for PE and FE",
    )
    parser.add_argument(
        "--top", type=int, default=1,
        help="matches retained per image (default: 1)",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("neb_image_matches.csv"),
        help="output CSV path (default: neb_image_matches.csv)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.top < 1:
        raise SystemExit("--top must be at least 1")
    neb_roots = args.neb_root or list(DEFAULT_NEB_ROOTS)
    phase_roots = args.phase_root or list(DEFAULT_PHASE_ROOTS)

    jobs = collect_neb_jobs(neb_roots)
    if not jobs:
        raise SystemExit("No NEB directories matching xxx-nnn-nnn were found.")
    candidates = collect_candidates(phase_roots)
    if not candidates:
        raise SystemExit("No relaxed phase candidates were found.")

    rows = []
    for job in jobs:
        rows.extend(match_job(job, candidates, args.top))
    if not rows:
        raise SystemExit(
            "No composition-compatible phase match was found for any NEB image."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    # Also print the compact result to stdout; csv.writer correctly quotes paths.
    writer = csv.DictWriter(sys.stdout, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    print(
        f"\nMatched {len(rows)} row(s) from {len(jobs)} NEB calculation(s); "
        f"saved to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
