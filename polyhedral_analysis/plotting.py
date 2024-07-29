import numpy as np
from scipy.stats import gaussian_kde # type: ignore
from cmcrameri import cm # type: ignore
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any

def _process_orientation_data(orientations: List[Tuple[float, float]]) -> Dict[str, Any]:
    """
    Process orientation data to prepare for plotting.

    Args:
        orientations (List[Tuple[float, float]]): List of (theta, phi) tuples.

    Returns:
        Dict[str, Any]: Dictionary containing processed data for plotting.

    Raises:
        ValueError: If orientations is empty, contains invalid data, or results in singular covariance matrix.
    """
    orientations_array = np.array(orientations)
    
    if orientations_array.size == 0:
        raise ValueError("The orientations list is empty.")
    
    if orientations_array.shape[1] != 2:
        raise ValueError("Each orientation should be a tuple of (theta, phi).")

    phi = np.linspace(-180, 180, 100)
    theta = np.linspace(0, 180, 100)
    phi_grid, theta_grid = np.meshgrid(phi, theta)

    positions = np.vstack([phi_grid.ravel(), theta_grid.ravel()])

    values = orientations_array.T[::-1]  # Swap and transpose
    
    try:
        kernel = gaussian_kde(values, bw_method=0.1)
        z = np.reshape(kernel(positions).T, phi_grid.shape)
    except np.linalg.LinAlgError:
        raise ValueError("The data resulted in a singular covariance matrix. Consider adding more varied data points.")

    return {
        'phi_grid': phi_grid,
        'theta_grid': theta_grid,
        'z': z
    }

def plot_orientation_distribution(orientations: List[Tuple[float, float]],
                                  title: Optional[str] = None,
                                  figsize: Tuple[int, int] = (10, 8),
                                  cmap = cm.lipari,
                                  fontsize: Optional[int] = None,
                                  plot: bool = False) -> plt.Figure:
    """
    Create a contour map of the probability distribution for (phi, theta) orientations.
    
    Args:
        orientations (List[Tuple[float, float]]): List of (theta, phi) tuples.
        title (Optional[str], optional): Title for the plot. Defaults to None.
        figsize (Tuple[int, int], optional): Figure size as (width, height). Defaults to (10, 8).
        cmap (matplotlib.colors.Colormap, optional): Colormap to use for the contour plot. Defaults to cm.lipari.
        fontsize (Optional[int], optional): Font size for labels and title. If None, uses matplotlib defaults.
        plot (bool, optional): Whether to display the plot immediately. Defaults to False.
    
    Returns:
        matplotlib.figure.Figure: The generated figure object.

    Raises:
        ValueError: If orientations is empty or contains invalid data.
    """
    processed_data = _process_orientation_data(orientations)

    fig, ax = plt.subplots(figsize=figsize)

    contour = ax.contourf(processed_data['phi_grid'], processed_data['theta_grid'], 
                          processed_data['z'], levels=20, cmap=cmap)
    if title:
        ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel(r'$\phi$ (degrees)', fontsize=fontsize)
    ax.set_ylabel(r'$\theta$ (degrees)', fontsize=fontsize)
    if fontsize:
        ax.tick_params(axis='both', which='major', labelsize=fontsize)

    cbar = fig.colorbar(contour)
    cbar.set_label('Probability Density', fontsize=fontsize)
    if fontsize:
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.ax.yaxis.get_offset_text().set_fontsize(fontsize)

    cbar.formatter.set_powerlimits((0, 0)) # type: ignore
    cbar.update_ticks()

    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 180)

    plt.tight_layout()
    
    if plot:
        plt.show()
    else:
        plt.close(fig)  # Close the figure to prevent it from being displayed
    
    return fig

# Example usage:
# orientations = [(theta1, phi1), (theta2, phi2), ...]
# fig = plot_orientation_distribution(orientations, title="Optional Title", fontsize=14)
# plt.show()