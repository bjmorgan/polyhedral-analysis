import unittest
from unittest.mock import patch, Mock
import numpy as np
import matplotlib.pyplot as plt
from polyhedral_analysis.plotting import plot_orientation_distribution, _process_orientation_data

class TestOrientationDistribution(unittest.TestCase):
    def setUp(self):
        # Use more varied data for all tests
        self.orientations = [(0, 0), (45, 45), (90, 90), (135, 135), (180, 180)]
        self.orientations += [(theta, phi) for theta in range(0, 181, 30) for phi in range(-180, 181, 60)]

    def test_process_orientation_data_valid_input(self):
        result = _process_orientation_data(self.orientations)
        self.assertIn('phi_grid', result)
        self.assertIn('theta_grid', result)
        self.assertIn('z', result)
        self.assertEqual(result['phi_grid'].shape, (100, 100))
        self.assertEqual(result['theta_grid'].shape, (100, 100))
        self.assertEqual(result['z'].shape, (100, 100))

    def test_process_orientation_data_empty_input(self):
        with self.assertRaises(ValueError):
            _process_orientation_data([])

    def test_process_orientation_data_invalid_input(self):
        with self.assertRaises(ValueError):
            _process_orientation_data([(0, 0, 0)])

    def test_process_orientation_data_singular_matrix(self):
        with self.assertRaises(ValueError):
            _process_orientation_data([(0, 0), (0, 0), (0, 0)])

    @patch('matplotlib.pyplot.show')
    def test_plot_orientation_distribution(self, mock_show):
        fig = plot_orientation_distribution(self.orientations, title="Test Plot", fontsize=14)
        
        # Check if the figure is created
        self.assertIsInstance(fig, plt.Figure)
        
        # Check if the figure has one axis (the main plot)
        self.assertEqual(len(fig.axes), 2)  # Main plot and colorbar
        
        ax = fig.axes[0]  # Main plot axis
        
        # Check title
        self.assertEqual(ax.get_title(), "Test Plot")
        
        # Check labels
        self.assertEqual(ax.get_xlabel(), r'$\phi$ (degrees)')
        self.assertEqual(ax.get_ylabel(), r'$\theta$ (degrees)')
        
        # Check axis limits
        self.assertEqual(ax.get_xlim(), (-180, 180))
        self.assertEqual(ax.get_ylim(), (0, 180))
        
        # Check if contourf was called (by checking if there are collections)
        self.assertTrue(len(ax.collections) > 0)
        
        # Close the figure to free memory
        plt.close(fig)

    def test_plot_orientation_distribution_empty_input(self):
        with self.assertRaises(ValueError):
            plot_orientation_distribution([])

    def test_plot_orientation_distribution_invalid_input(self):
        with self.assertRaises(ValueError):
            plot_orientation_distribution([(0, 0, 0)])

if __name__ == '__main__':
    unittest.main()
	