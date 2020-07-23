import unittest
from unittest.mock import Mock, MagicMock, patch, mock_open
from polyhedral_analysis.atom import Atom
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from pymatgen import Site, Lattice
import numpy as np
import json

class TestAtomInit(unittest.TestCase):

    def test_atom_is_initialised(self):
        index = 123
        site = Mock(spec=Site)
        label = 'foo'
        atom = Atom(index=index, site=site, label=label)
        self.assertEqual(atom.index, index)
        self.assertEqual(atom.site, site)
        self.assertEqual(atom.label, label)
        self.assertEqual(atom.in_polyhedra, [])
        self.assertEqual(atom.neighbours, {})

    def test_atom_is_initialised_without_label(self):
        index = 123
        site = Mock(spec=Site)
        atom = Atom(index=index, site=site)
        self.assertEqual(atom.index, index)
        self.assertEqual(atom.site, site)
        self.assertEqual(atom.label, None)
        self.assertEqual(atom.in_polyhedra, [])
        self.assertEqual(atom.neighbours, {})

class TestAtom(unittest.TestCase):

    def setUp(self):
        index = 123
        site = Mock(spec=Site)
        label = 'foo'
        site.frac_coords = np.array([1.0, 2.0, 3.0])
        site.coords = np.array([10.0, 11.0, 12.0])
        site.lattice = Mock(spec=Lattice)
        self.site = site
        self.atom = Atom(index=index, site=site, label=label)

    def test_frac_coords(self):
        np.testing.assert_array_equal(
            self.atom.frac_coords, self.site.frac_coords)

    def test_coords(self):
        np.testing.assert_array_equal(self.atom.coords, self.site.coords)

    def test_lattice(self):
        self.assertEqual(self.atom.lattice, self.site.lattice)

    def test___lt__(self):
        index = 456
        site = MagicMock(spec=Site)
        label = 'bar'
        self.frac_coords = np.array([2.0, 3.0, 4.0])
        site.coords = np.array([20.0, 21.0, 22.0])
        site.lattice = Mock(spec=Lattice)
        other_atom = Atom(index=index, site=site, label=label)
        self.assertTrue(self.atom < other_atom)

    def test___lt__raises_TypeError_if_not_passed_Atom(self):
        with self.assertRaises(TypeError):
            self.atom < 'foo'    

    def test___eq___when_true(self):
        index = 123
        site = Mock(spec=Site)
        label = 'foo'
        site.frac_coords = np.array([1.0, 2.0, 3.0])
        site.coords = np.array([10.0, 11.0, 12.0])
        site.lattice = Mock(spec=Lattice)
        other_atom = Atom(index=index, site=site, label=label)
        self.assertEqual(self.atom, other_atom)

    def test___eq___when_false(self):
        index = 456
        site = Mock(spec=Site)
        label = 'foo'
        site.frac_coords = np.array([1.0, 2.0, 3.0])
        site.coords = np.array([10.0, 11.0, 12.0])
        site.lattice = Mock(spec=Lattice)
        other_atom = Atom(index=index, site=site, label=label)
        self.assertNotEqual(self.atom, other_atom)

    def test___eq__is_False_if_not_passed_Atom(self):
        self.assertNotEqual(self.atom, 'foo')

    def test___hash__(self):
        self.assertEqual(hash(self.atom), self.atom.index)

    def test_distance(self):
        index = 456
        site = Mock(spec=Site)
        label = 'foo'
        other_atom = Atom(index=index, site=site, label=label)
        self.atom.site.distance = Mock(return_value=2.3)
        self.assertEqual(self.atom.distance(other_atom), 2.3)

    def test_as_dict(self):
        self.atom.site.as_dict = Mock(return_value={'key': 'value'})
        mock_other_polyhedra = Mock(spec=CoordinationPolyhedron)
        mock_other_polyhedra.index = 15
        self.atom.in_polyhedra = [mock_other_polyhedra]
        expected_dict = {'index': 123,
                         'site': {'key': 'value'},
                         'label': 'foo', 'in_polyhedra': [15],
                         'neighbours': {},
                         '@module': 'polyhedral_analysis.atom',
                         '@class': 'Atom'}
        self.assertEqual(self.atom.as_dict(), expected_dict)

    def test_to_json(self):
        self.atom.site.as_dict = Mock(return_value={'key': 'value'})
        mock_other_polyhedra = Mock(spec=CoordinationPolyhedron)
        mock_other_polyhedra.index = 15
        self.atom.in_polyhedra = [mock_other_polyhedra]
        expected_dict = {'index': 123,
                         'site': {'key': 'value'},
                         'label': 'foo', 'in_polyhedra': [15],
                         'neighbours': {},
                         '@module': 'polyhedral_analysis.atom',
                         '@class': 'Atom'}
        self.assertEqual(self.atom.to(), json.dumps(expected_dict))

    def test_to_json_file(self):
        self.atom.site.as_dict = Mock(return_value={'key': 'value'})
        mock_other_polyhedra = Mock(spec=CoordinationPolyhedron)
        mock_other_polyhedra.index = 15
        self.atom.in_polyhedra = [mock_other_polyhedra]
        expected_dict = {'index': 123,
                         'site': {'key': 'value'},
                         'label': 'foo', 'in_polyhedra': [15],
                         'neighbours': {},
                         '@module': 'polyhedral_analysis.atom',
                         '@class': 'Atom'}
        with patch('polyhedral_analysis.atom.zopen', mock_open(), create=True) as m:
            self.atom.to(filename='filename')
        m.assert_called_once_with('filename', 'wt')
        m().write.assert_called_once_with(json.dumps(expected_dict))

    def test_to_json_file_with_fmt_equals_JSON(self):
        self.atom.site.as_dict = Mock(return_value={'key': 'value'})
        mock_other_polyhedra = Mock(spec=CoordinationPolyhedron)
        mock_other_polyhedra.index = 15
        self.atom.in_polyhedra = [mock_other_polyhedra]
        expected_dict = {'index': 123,
                         'site': {'key': 'value'},
                         'label': 'foo', 'in_polyhedra': [15],
                         'neighbours': {},
                         '@module': 'polyhedral_analysis.atom',
                         '@class': 'Atom'}
        with patch('polyhedral_analysis.atom.zopen', mock_open(), create=True) as m:
            self.atom.to(filename='filename', fmt='JSON')
        m.assert_called_once_with('filename', 'wt')
        m().write.assert_called_once_with(json.dumps(expected_dict))

    def test_to_json_file_with_unknown_fmt_raises_ValueError(self):
        self.atom.site.as_dict = Mock(return_value={'key': 'value'})
        mock_other_polyhedra = Mock(spec=CoordinationPolyhedron)
        mock_other_polyhedra.index = 15
        self.atom.in_polyhedra = [mock_other_polyhedra]
        with self.assertRaises(ValueError):
            self.atom.to(filename='filename', fmt='FOO')


if __name__ == '__main__':
    unittest.main()
