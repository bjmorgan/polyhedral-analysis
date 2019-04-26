from pymatgen.io.vasp import Xdatcar
from .configuration import Configuration
import re
import copy

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
            read_config_numbers (:obj:`bool`, optional): Read configuration frame numbers from the XDATCAR.
            config_numbers (:obj:`list`, optional): Optional list of integers to use as frame numbers for
                each configuration in the XDATCAR input. If this argument is set, it will override
                the `read_config_numbers` argument.
            verbose (:obj:`bool`, optional): verbose output while parsing the XDATCAR input. Default is False.

        Returns:
            None

        Notes:
            if `read_config_numbers` is set to False, and `config_numbers` is not set,
            the XDATCAR configurations will be numbered 1, 2, 3, 4, ….
        """
        self.xdatcar = Xdatcar( xdatcar )
        if read_config_numbers and not config_numbers:
            config_numbers = read_config_numbers_from_xdatcar( xdatcar )
        elif config_numbers:
            config_numbers = config_numbers
        else:
            config_numbers = list( range( 1, len( self.xdatcar.structures ) + 1 ) ) 
        if len( config_numbers ) != len( self.xdatcar.structures ):
            raise ValueError( 'number of configuration numbers != number of structures' )
        self.recipes = recipes
        # generate polyhedra configurations
        self.configurations = []
        for n, s in zip( config_numbers, self.structures ):
            if verbose:
                print( 'Reading configuration {}'.format( n ), flush=True )
            c = Configuration( structure=s, recipes=self.recipes, config_number=n )
            self.configurations.append( c )

    def extend( self, other, offset=1 ):
        """
        Extend this Trajectory with the data from another Trajectory.

        Args:
            other (`Trajectory`): The trajectory to be appended.
            offset (int): Offset to apply to the frame numbers.

        Returns:
            None
        """
        extended_configurations = copy.deepcopy( other.configurations )
        for c in extended_configurations:
            c.config_number += offset + self.config_numbers[-1]
        self.configurations.extend( extended_configurations )
        #self.config_numbers.extend( [ n + offset + self.config_numbers[-1] for n in other.config_numbers ] ) 
        # TODO need to update individual configuration numbers

    @property
    def config_numbers( self ):
        return [ c.config_number for c in self.configurations ]

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
        
 
