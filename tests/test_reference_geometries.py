import unittest

import numpy as np

from polyhedral_analysis.reference_geometries import (
    REFERENCE_GEOMETRIES,
    get_reference_geometry,
)


class TestReferenceGeometries(unittest.TestCase):

    def test_all_17_geometries_are_present(self):
        self.assertEqual(len(REFERENCE_GEOMETRIES), 17)

    def test_octahedron_coordinates(self):
        expected = np.array([
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        np.testing.assert_array_equal(
            REFERENCE_GEOMETRIES["Octahedron"], expected)

    def test_coordination_numbers(self):
        expected_sizes = {
            "Tetrahedron": 4,
            "Trigonal bipyramid": 5,
            "Square pyramid": 5,
            "Octahedron": 6,
            "Trigonal prism": 6,
            "Pentagonal bipyramid": 7,
            "Square-face capped trigonal prism": 7,
            "Face-capped octahedron": 7,
            "Cube": 8,
            "Square antiprism": 8,
            "Square-face bicapped trigonal prism": 8,
            "Triangular-face bicapped trigonal prism": 8,
            "Dodecahedron with triangular faces": 8,
            "Hexagonal bipyramid": 8,
            "Bicapped octahedron (opposed cap faces)": 8,
            "Bicapped octahedron (cap faces with one atom in common)": 8,
            "Bicapped octahedron (cap faces with one edge in common)": 8,
        }
        # This also implicitly checks that we haven't missed a geometry name
        self.assertEqual(set(expected_sizes.keys()),
                         set(REFERENCE_GEOMETRIES.keys()))
        for name, size in expected_sizes.items():
            self.assertEqual(len(REFERENCE_GEOMETRIES[name]), size,
                             msg=f"{name} should have {size} points")


class TestGetReferenceGeometry(unittest.TestCase):

    def test_returns_correct_array(self):
        result = get_reference_geometry("Octahedron")
        np.testing.assert_array_equal(
            result, REFERENCE_GEOMETRIES["Octahedron"])

    def test_returns_copy(self):
        result = get_reference_geometry("Octahedron")
        result[0, 0] = 999.0
        self.assertNotEqual(
            REFERENCE_GEOMETRIES["Octahedron"][0, 0], 999.0)

    def test_unknown_name_raises_key_error(self):
        with self.assertRaises(KeyError):
            get_reference_geometry("Nonexistent")


if __name__ == '__main__':
    unittest.main()
