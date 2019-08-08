from pymatgen.io.vasp import Xdatcar, Structure
from .configuration import Configuration
import re
import copy
from monty.io import zopen
import pickle
from .utils import flatten
import multiprocessing
from tqdm import tqdm, tqdm_notebook
from functools import partial

def read_config_numbers_from_xdatcar( filename ):
    with zopen( filename ) as f:
        file_string = f.read()
        try: # attempt to decode byte object
            file_string = file_string.decode()
        except AttributeError:
            pass
        return [ int(s) for s in re.findall( 'Direct configuration=\s+(\d+)', file_string ) ]

class Trajectory:

    def __init__( self, structures, configurations, config_numbers=None ):
        for s in structures:
            if not isinstance( s, Structure ):
                raise TypeError( "structures argument should contain pymatgen Structure objects" )
        for c in configurations:
            if not isinstance( c, Configuration ):
                raise TypeError( "configurations argument should contain Configuration objects" )
        self._structures = structures
        self._configurations = configurations
        if not config_numbers:
            config_numbers = list( range( 1, len( self.structures ) + 1 ) )
        if len( config_numbers ) != len( self.structures ):
            raise ValueError( 'number of configuration numbers != number of structures: {}, {}'.format( config_numbers, len( self.structures ) ) )
        for n, c in zip( config_numbers, self._configurations ):
            c.config_number = n

    @classmethod
    def _get_configuration( cls, args ):
         return Configuration( structure=args['structure'], 
                               recipes=args['recipes'], 
                               config_number=args['config_number'] )
    
    @classmethod
    def from_structures( cls, structures, recipes, config_numbers=None, verbose=False, 
                         ncores=None, progress=False ):
        """
        Generate a :obj:`Trajectory` object by applying one or more polyhedral recipes to
        a series of `pymatgen` :obj:`Structure` objects.,

        Args:
            structures (list(:obj:`pymatgen.Structure`)): A list of `pymatgen` 
                :obj:`Structure`` objects.
            recipes (list(PolyhedraRecipe): List of `PolyhedraRecipe` recipes, where each recipe
                defines how to construct a set of `CoordinationPolyhedra` for each configuration.
            config_numbers (:obj:`list`, optional): Optional list of integers to use as frame 
                numbers for each structure. If this argument is not set, the
                configurations will be numnbered 1, 2, 3 ….
            verbose (:obj:`bool`, optional): verbose output while parsing the input structures.
                Default is False.
            ncores (:obj:`int`, optional): TODO.
            progress (:obj:`str`, optional): Show a progress bar.
               Setting to ``True`` gives a simple progress bar.
               Setting to ``"notebook"`` gives a Jupyter notebook compatible progress bar.
               Default is ``False``.

        Returns:
            None
        """
        if progress == True:
            progress_bar = partial( tqdm, unit=' configurations' )
        elif progress == 'notebook':
            progress_bar = partial( tqdm_notebook, unit=' configurations' )
        else:
            progress_bar = lambda x: x
        if not config_numbers:
            config_numbers = list( range( 1, len( self.structures ) + 1 ) ) 
        if len( config_numbers ) != len( structures ):
            raise ValueError( 'number of configuration numbers != number of structures: {}, {}'.format( config_numbers, len( structures ) ) )
        # generate polyhedra configurations
        if ncores:
            with multiprocessing.Pool( ncores ) as p:
                args = [ { 'structure': s, 'recipes': recipes, 'config_number': n }
                    for n, s in zip( config_numbers, structures ) ]
                configurations = p.map( cls._get_configuration, progress_bar( args ) )
        else:
            configurations = []
            for n, s in progress_bar( zip( config_numbers, structures ) ):
                if verbose:
                    print( 'Reading configuration {}'.format( n ), flush=True )
                c = Configuration( structure=s, recipes=recipes, config_number=n )
                configurations.append( c )
        return cls( structures=structures, configurations=configurations, config_numbers=config_numbers )

    def extend( self, other, offset=0 ):
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

    def __add__( self, other ):
        """TODO"""
        summed_structures = self.structures + other.structures
        summed_configurations = self.configurations + other.configurations
        new_trajectory = Trajectory( summed_structures, summed_configurations )
        return new_trajectory

    @property
    def config_numbers( self ):
        return [ c.config_number for c in self.configurations ]

    @property
    def structures( self ):
        return self._structures

    @property
    def configurations( self ):
        return self._configurations

    #TODO implement different ways of creating a Trajectory object, e.g. 
    #@class_method
    #def from_configurations( self, configurations ):
    #    TODO

    @classmethod
    def from_xdatcar( cls, filename, recipes, read_config_numbers=True, config_numbers=None,
                      verbose=False, progress=None, ncores=None ):
        """
        Args:
            filename (str): Filename for a VASP XDATCAR file.
            recipes (list(:obj:`PolyhedraRecipe`): List of `PolyhedraRecipe` recipes, 
                where each recipe defines how to construct a set of `CoordinationPolyhedra` 
                for each configuration.
            read_config_numbers (:obj:`bool`, optional): Read configuration frame numbers from the XDATCAR.
            config_numbers (:obj:`list`, optional): Optional list of integers to use as frame numbers for
                each configuration in the XDATCAR input. If this argument is set, it will override
                the `read_config_numbers` argument.

                      config_numbers=None, verbose=False ):

        Returns:
            (:obj:`Trajectory`): a :obj:`Trajectory` object.

        Notes:
            If `read_config_numbers` is set to False, and `config_numbers` is not set,
            the XDATCAR configurations will be numbered 1, 2, 3, 4, ….
        """
        xdatcar = Xdatcar( filename )
        structures = xdatcar.structures
        if read_config_numbers and not config_numbers:
            config_numbers = read_config_numbers_from_xdatcar( filename )
        elif config_numbers:
            config_numbers = config_numbers
        else:
            config_numbers = list( range( 1, len( structures ) + 1 ) )
        if len( config_numbers ) != len( structures ):
            raise ValueError( 'number of configuration numbers != number of structures: {}, {}'.format( config_numbers, len( structures ) ) )
        return cls.from_structures( structures, recipes, config_numbers, verbose=verbose, progress=progress,
                                    ncores=ncores )

    @classmethod
    def from_xdatcars( cls, filenames, recipes, verbose=False, ncores=None, progress=None ):
        """
        TODO
        """
        if ncores and len( filenames ) > 1:
            with multiprocessing.Pool( ncores ) as p:
                xdatcars = p.map( _get_xdatcar, filenames )
        else:
            xdatcars = [ Xdatcar( f ) for f in filenames ]
        structures = flatten( [ x.structures for x in xdatcars ] )
        config_numbers = list( range( 1, len( structures ) + 1 ) )
        return cls.from_structures( structures, recipes, config_numbers, verbose, ncores=ncores,
            progress=progress ) 

def _get_xdatcar( filename ):
    """
    Internal method to support multiprocessing.
    """
    return Xdatcar( filename )

