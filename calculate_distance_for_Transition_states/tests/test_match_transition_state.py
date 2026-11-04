import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase import Atoms

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from match_transition_state import compare, transform_doped
from match_neb_images import collect_images, collect_neb_jobs


class MatchTests(unittest.TestCase):
    def test_compare_ignores_origin_and_same_element_order(self):
        reference = Atoms(
            ["Hf", "O", "O"],
            scaled_positions=[[0.1, 0.2, 0.3], [0.3, 0.4, 0.5], [0.7, 0.8, 0.9]],
            cell=[5, 5, 5],
            pbc=True,
        )
        candidate = Atoms(
            ["O", "Hf", "O"],
            scaled_positions=[[0.55, 0.65, 0.75], [0.35, 0.45, 0.55], [0.95, 0.05, 0.15]],
            cell=[5, 5, 5],
            pbc=True,
        )
        rmsd, maximum, strain = compare(reference, candidate)
        self.assertAlmostEqual(rmsd, 0.0, places=12)
        self.assertAlmostEqual(maximum, 0.0, places=12)
        self.assertAlmostEqual(strain, 0.0, places=12)

    def test_doped_transform_preserves_atom_sequence(self):
        atoms = Atoms(
            ["Hf", "Zr", "O"],
            scaled_positions=[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]],
            cell=np.diag([5.0, 6.0, 7.0]),
            pbc=True,
        )
        transformed = transform_doped(atoms, "xyz-000-000")
        self.assertEqual(transformed.get_chemical_symbols(), ["Hf", "Zr", "O"])
        np.testing.assert_allclose(
            transformed.get_scaled_positions(), atoms.get_scaled_positions(), atol=1e-12
        )

    def test_neb_discovery_and_numeric_image_sorting(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "fe1"
            neb = root / "xyz-010-101"
            doped = neb / "a-xyz-010-101"
            ignored = root / "not-a-neb-path"
            for directory in (neb / "10", neb / "02", doped, ignored):
                directory.mkdir(parents=True, exist_ok=True)
            atoms = Atoms(
                ["Hf", "O"],
                scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
                cell=[5, 5, 5],
                pbc=True,
            )
            atoms.write(neb / "10" / "POSCAR", format="vasp")
            atoms.write(neb / "02" / "CONTCAR", format="vasp")

            jobs = collect_neb_jobs([root])
            self.assertEqual(
                [(job.path_operator, job.doping) for job in jobs],
                [("xyz-010-101", "pristine"), ("xyz-010-101", "a")],
            )
            self.assertEqual(
                [image for image, _ in collect_images(jobs[0])], ["02", "10"]
            )


if __name__ == "__main__":
    unittest.main()
