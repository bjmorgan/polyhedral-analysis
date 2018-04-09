import unittest
import numpy as np
from polyhedral_analysis.polyhedra_recipe import ( PolyhedraRecipe,
                                                   matching_sites )
from pymatgen import Structure, Lattice

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

class TestPolyhedraRecipeFunctions( unittest.TestCase ):

    def test_matching_sites( self ):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.0, 0.5, 0.5 ],
                             [ 0.5, 0.0, 0.5 ] ] )
        atom_list = [ 'Li' ] * len( coords )
        lattice = Lattice.from_parameters( a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        structure = Structure( lattice, atom_list, coords )    
        ref_coords = np.array( [ [ 0.1, 0.1, 0.1 ],
                                 [ 0.0, 0.4, 0.5 ] ] )
        ref_atom_list = [ 'Na' ] * len( ref_coords )
        ref_structure = Structure( lattice, ref_atom_list, ref_coords )
        matched_sites = matching_sites( structure, ref_structure )
        self.assertEqual( len( matched_sites ), 2 )
        self.assertEqual( matched_sites[0], [ structure[0], 0 ] )
        self.assertEqual( matched_sites[1], [ structure[2], 2 ] )

    def test_matching_sites_with_species( self ):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.0, 0.5, 0.5 ],
                             [ 0.5, 0.0, 0.5 ] ] )
        atom_list = [ 'Li', 'Mg', 'Mg', 'Mg' ]
        lattice = Lattice.from_parameters( a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        structure = Structure( lattice, atom_list, coords )
        ref_coords = np.array( [ [ 0.1, 0.1, 0.1 ],
                                 [ 0.0, 0.4, 0.5 ] ] )
        ref_atom_list = [ 'Na' ] * len( ref_coords )
        ref_structure = Structure( lattice, ref_atom_list, ref_coords )
        matched_sites = matching_sites( structure, ref_structure, species=['Li'] )
        self.assertEqual( len( matched_sites ), 1 )
        self.assertEqual( matched_sites[0], [ structure[0], 0 ] )


if __name__ == '__main__':
    unittest.main()

