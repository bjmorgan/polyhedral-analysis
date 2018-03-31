import unittest
from unittest.mock import Mock, patch
from polyhedral_analysis.atom import Atom
from pymatgen import Site, Lattice
import numpy as np


class TestAtomInit( unittest.TestCase ):

    def test_atom_is_initialised( self ):
        index = 123
        site = Mock( spec=Site )
        label = 'foo'
        atom = Atom( index=index, site=site, label=label )
        self.assertEqual( atom.index, index )
        self.assertEqual( atom.site, site )
        self.assertEqual( atom.label, label )
        self.assertEqual( atom.in_polyhedra, [] )
        self.assertEqual( atom.neighbours, None )

    def test_atom_is_initialised_without_label( self ):
        index = 123
        site = Mock( spec=Site )
        atom = Atom( index=index, site=site )
        self.assertEqual( atom.index, index )
        self.assertEqual( atom.site, site )
        self.assertEqual( atom.label, None )
        self.assertEqual( atom.in_polyhedra, [] )
        self.assertEqual( atom.neighbours, None )

class TestAtom( unittest.TestCase ):

    def setUp( self ):
        index = 123
        site = Mock( spec=Site )
        label = 'foo'
        site.frac_coords = np.array( [ 1.0, 2.0, 3.0 ] )
        site.coords = np.array( [ 10.0, 11.0, 12.0 ] )
        site.lattice = Mock( spec=Lattice )
        self.site = site
        self.atom = Atom( index=index, site=site, label=label )
       
    def test_frac_coords( self ):
        np.testing.assert_array_equal( self.atom.frac_coords, self.site.frac_coords )

    def test_coords( self ):
        np.testing.assert_array_equal( self.atom.coords, self.site.coords )

    def test_lattice( self ):
        self.assertEqual( self.atom.lattice, self.site.lattice )

    def test___lt__( self ):
        index = 456
        site = Mock( spec=Site )
        label = 'bar'
        self.frac_coords = np.array( [ 2.0, 3.0, 4.0 ] )
        site.coords = np.array( [ 20.0, 21.0, 22.0 ] )
        site.lattice = Mock( spec=Lattice )
        other_atom = Atom( index=index, site=site, label=label )
        self.assertTrue( self.atom < other_atom )

    def test___eq___when_true( self ):
        index = 123
        site = Mock( spec=Site )
        label = 'foo' 
        site.frac_coords = np.array( [ 1.0, 2.0, 3.0 ] )
        site.coords = np.array( [ 10.0, 11.0, 12.0 ] )
        site.lattice = Mock( spec=Lattice )
        other_atom = Atom( index=index, site=site, label=label )
        self.assertEqual( self.atom, other_atom )

    def test___eq___when_false( self ):
        index = 456
        site = Mock( spec=Site )
        label = 'foo'
        site.frac_coords = np.array( [ 1.0, 2.0, 3.0 ] )
        site.coords = np.array( [ 10.0, 11.0, 12.0 ] )
        site.lattice = Mock( spec=Lattice )
        other_atom = Atom( index=index, site=site, label=label )
        self.assertNotEqual( self.atom, other_atom )

if __name__ == '__main__':
    unittest.main()
