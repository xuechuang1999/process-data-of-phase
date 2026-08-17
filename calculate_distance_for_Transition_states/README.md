# Transition-state structure matching

`src/match_transition_state.py` compares a VASP transition-state structure
against relaxed PE/FE structures and all representatives made by
`src/reduced.py`.

For a doped directory (`a-n_*`, `b-n_*`, or `c-n_*`), the code applies the
parent structure's symmetry representative directly to the doped relaxed
structure **without sorting its atoms**.  The dopant-to-site correspondence is
therefore retained.  During comparison, atoms may be optimally permuted only
within the same chemical element, and distances use periodic boundary
conditions.

## Usage

Run from this repository (or put `src` on `PYTHONPATH`):

```bash
python src/match_transition_state.py /path/to/path1/TS/CONTCAR /path/to/path2/TS/CONTCAR \
  --root /gs/home/xuechuang/work/ZrO2/hafnia/1-PE/states \
  --root /gs/home/xuechuang/work/ZrO2/hafnia/2-FE/states \
  --top 20 \
  --csv matches.csv
```

One or several transition-state files can be supplied in the same run; the
reported top matches are ranked independently for each file.  The ranking uses
atomic Cartesian RMSD in ångström.  It also reports the
largest matched-atom displacement and relative cell-length strain.  Candidates
whose atom counts or compositions differ from the transition state are skipped.

Each `n_*` directory is expected to contain:

```text
n_name/
├── reduced/<operator>/{CONTCAR|POSCAR}
├── a-n_name/{CONTCAR|POSCAR}
├── b-n_name/{CONTCAR|POSCAR}
└── c-n_name/{CONTCAR|POSCAR}
```

`CONTCAR` is preferred; `POSCAR` is used as a fallback.  Python dependencies
are `numpy`, `ase`, `spglib`, and `scipy`.

## Match every image in NEB paths

`src/match_neb_images.py` scans `fe1` and `fe2` for path directories named
`xxx-nnn-nnn`, then scans their numeric image directories.  It also scans the
nested `a-xxx-nnn-nnn`, `b-xxx-nnn-nnn`, and `c-xxx-nnn-nnn` calculations.

The cluster paths described above are built-in defaults, so on the cluster the
short form is:

```bash
python src/match_neb_images.py --output neb_image_matches.csv
```

All paths can be overridden:

```bash
python src/match_neb_images.py \
  --neb-root /gs/home/xuechuang/work/ZrO2/hzo/phase-transition/fe1 \
  --neb-root /gs/home/xuechuang/work/ZrO2/hzo/phase-transition/fe2 \
  --phase-root /gs/home/xuechuang/work/ZrO2/hafnia/1-PE/states \
  --phase-root /gs/home/xuechuang/work/ZrO2/hafnia/2-FE/states \
  --top 1 \
  --output neb_image_matches.csv
```

Pristine NEB paths are compared only with pristine phase candidates.  Likewise,
`a`, `b`, and `c` paths are compared only with candidates carrying the same
doping label.  The CSV reports the NEB group, path operator, image number,
closest phase, closest `n_*` state, closest symmetry operator, RMSD, maximum
atomic displacement, cell strain, and both source files.
