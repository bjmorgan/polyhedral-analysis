import unittest
from unittest.mock import Mock, patch
from polyhedral_analysis.atom import Atom
from pymatgen import Site

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

if __name__ == '__main__':
    unittest.main()
