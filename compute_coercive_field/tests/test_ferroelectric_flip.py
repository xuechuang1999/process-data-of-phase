from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from coercive_field.ferroelectric_flip import (
    invert_fractional_positions,
    minimum_image,
    reflect_about_reference,
)
from coercive_field.vasp_io import Poscar, parse_structure, read_poscar, write_poscar


class FerroelectricFlipTests(unittest.TestCase):
    def make_structure(self, frac: np.ndarray) -> Poscar:
        return Poscar(
            comment="test",
            lattice=np.diag([5.0, 5.1, 5.2]),
            symbols=["Hf", "O"],
            counts=[1, 2],
            frac=np.asarray(frac, dtype=float),
            selective=True,
            extras=[["T", "T", "T"], ["F", "F", "F"], ["T", "F", "T"]],
        )

    def test_inverting_twice_recovers_periodic_positions(self) -> None:
        initial = self.make_structure(
            np.array([[0.1, 0.2, 0.3], [0.9, 0.4, 0.7], [0.0, 0.8, 0.6]])
        )
        origin = np.array([0.25, 0.0, 0.5])

        first = invert_fractional_positions(initial, origin)
        intermediate = self.make_structure(first)
        second = invert_fractional_positions(intermediate, origin)

        np.testing.assert_allclose(minimum_image(second - initial.frac), 0.0, atol=1.0e-12)

    def test_reference_reflection_reverses_minimum_image_displacement(self) -> None:
        reference = self.make_structure(
            np.array([[0.98, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
        )
        initial = self.make_structure(
            np.array([[0.02, 0.2, 0.3], [0.35, 0.55, 0.6], [0.7, 0.8, 0.82]])
        )

        flipped = reflect_about_reference(initial, reference)

        initial_displacement = minimum_image(initial.frac - reference.frac)
        flipped_displacement = minimum_image(flipped - reference.frac)
        np.testing.assert_allclose(flipped_displacement, -initial_displacement, atol=1.0e-12)

    def test_shared_poscar_io_round_trip_preserves_metadata(self) -> None:
        initial = self.make_structure(
            np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "POSCAR"
            write_poscar(output, initial, initial.frac, "round trip")
            loaded = read_poscar(output)
            lattice, volume = parse_structure(output)

        self.assertEqual(loaded.symbols, initial.symbols)
        self.assertEqual(loaded.counts, initial.counts)
        self.assertEqual(loaded.extras, initial.extras)
        self.assertTrue(loaded.selective)
        np.testing.assert_allclose(loaded.frac, initial.frac)
        np.testing.assert_allclose(lattice, initial.lattice)
        self.assertAlmostEqual(volume, float(np.linalg.det(initial.lattice)))

    def test_parse_structure_does_not_require_complete_coordinates(self) -> None:
        header_only = """incomplete structure
1.0
5.0 0.0 0.0
0.0 5.1 0.0
0.0 0.0 5.2
Hf O
1 2
Direct
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            path.write_text(header_only, encoding="utf-8")
            lattice, volume = parse_structure(path)
            with self.assertRaises(ValueError):
                read_poscar(path)

        np.testing.assert_allclose(lattice, np.diag([5.0, 5.1, 5.2]))
        self.assertAlmostEqual(volume, 5.0 * 5.1 * 5.2)


if __name__ == "__main__":
    unittest.main()
