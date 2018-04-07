from pymatgen.io.vasp import Xdatcar
from .configuration import Configuration
import re

def read_config_numbers_from_xdatcar( filename ):
    with open( filename ) as f:
        return [ int(s) for s in re.findall( 'Direct configuration=\s+(\d+)', f.read() ) ]

class Trajectory:

    def __init__( self, xdatcar, recipes, read_config_numbers=True, 
                  config_numbers=None, verbose=False ):
        """
        Class describing a single simulation trajectory.

        Args:
            xdatcar (str): filename for a VASP XDATCAR file.
            recipes (list(PolyhedraRecipe): List of `PolyhedraRecipe` recipes, where each recipe
                defines how to construct a set of `CoordinationPolyhedra` for each configuration.
            read_config_numbers (:opt:`bool`): Read configuration frame numbers from the XDATCAR.
            config_numbers (:opt:`list`): Optional list of integers to use as frame numbers for
                each configuration in the XDATCAR input. If this argument is set, it will override
                the `read_config_numbers` argument.
            verbose (:opt:`bool`): verbose output while parsing the XDATCAR input. Default is False.

        Returns:
            None

        Notes:
            if `read_config_numbers` is set to False, and `config_numbers` is not set,
            the XDATCAR configurations will be numbered 1, 2, 3, 4, ….
        """
        self.xdatcar = Xdatcar( xdatcar )
        if read_config_numbers and not config_numbers:
            self.config_numbers = read_config_numbers_from_xdatcar( xdatcar )
        elif config_numbers:
            self.config_numbers = config_numbers
        else:
            self.config_numbers = list( range( 1, len( self.xdatcar.structures ) + 1 ) ) 
        assert( len( self.config_numbers ) == len( self.xdatcar.structures ) )
        self.recipes = recipes
        # generate polyhedra configurations
        self.configurations = []
        for n, s in zip( self.config_numbers, self.structures ):
            if verbose:
                print( 'Reading configuration {}'.format( n ), flush=True )
            c = Configuration( structure=s, recipes=self.recipes, config_number=n )
            self.configurations.append( c )

    @property
    def structures( self ):
        return self.xdatcar.structures

    #TODO implement different ways of creating a Trajectory object, e.g. 
    #@class_method
    #def from_configurations( self, configurations ):
    #    TODO

    #@class_method
    #def from_xdatcar( self, xdatcar, … ):
    #    TODO
        
 
