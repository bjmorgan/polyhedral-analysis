import unittest
from unittest.mock import Mock, patch
from pymatgen import Structure, Site
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.coordination_polyhedron import CoordinationPolyhedron
from polyhedral_analysis.atom import Atom

class TestConfigurationInit( unittest.TestCase ):

    def test_configuration_is_initialised( self ):
        mock_sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        mock_sites[0].species_string = 'A'
        mock_sites[1].species_string = 'B'
        mock_structure = Mock( spec=Structure )
        mock_structure.sites = mock_sites
        mock_polyhedron = Mock( spec=CoordinationPolyhedron )
        mock_atoms = [ Mock( spec=Atom ), Mock( spec=Atom ) ]
        mock_polyhedron.central_atom = mock_atoms[0]
        mock_polyhedron.vertices = [ mock_atoms[1] ]
        mock_recipe = Mock( spec=PolyhedraRecipe )
        mock_recipe.find_polyhedra = Mock( return_value=[ mock_polyhedron ] )
        with patch( 'polyhedral_analysis.configuration.Atom' ) as mock_Atom:
            mock_Atom.side_effect = mock_atoms
            configuration = Configuration( structure=mock_structure, recipes=[ mock_recipe ] )
        self.assertEqual( configuration.structure, mock_structure )
        self.assertEqual( configuration.config_number, None )
        self.assertEqual( configuration.central_atoms, [ mock_atoms[0] ] )
        self.assertEqual( configuration.coordination_atoms, [ mock_atoms[1] ] )
        self.assertEqual( configuration.atoms, mock_atoms )
        
if __name__ == '__main__':
    unittest.main()

