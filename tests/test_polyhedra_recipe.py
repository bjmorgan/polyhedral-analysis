import unittest
import numpy as np
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe

class TestPolyhedraRecipeInit( unittest.TestCase ):

    def test_polyhedra_recipe_init_raises_ValueError_for_invalid_method( self ):
        invalid_method = 'foo'
        with self.assertRaises( ValueError ):
            PolyhedraRecipe( method=invalid_method )

    def test_polyhedra_recipe_init_raises_ValueError_for_missing_central_atom_definition( self ):
        with self.assertRaises( ValueError ) as e:
           PolyhedraRecipe( method='distance cutoff', 
                            central_atom_list=None, 
                            central_atom_list_generator=None )
        self.assertEqual( str(e.exception), 'central_atom_list or central_atom_list_generator must be specified' )

    def test_polyhedra_recipe_init_raises_ValueError_for_missing_coordination_atom_definition( self ):
        with self.assertRaises( ValueError ) as e:
            PolyhedraRecipe( method='distance cutoff',
                             central_atom_list = [ 1, 2, 3 ],
                             coordination_atom_list=None,
                             coordination_atom_list_generator=None )
        self.assertEqual( str( e.exception ), 'coordination_atom_list or coordination_atom_list_generator must be specified' )  

if __name__ == '__main__':
    unittest.main()

