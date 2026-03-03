__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

import matplotlib.pyplot as plt

from ..plotting import curvature
from ..macrodata_refinement.utils import read

# ------ default style parameters ------
main_linewidth=1.3
main_linestyle='-'

E_zero_color='gray'
E_zero_linewidth=0.2
E_zero_linestyle='--'

k_color='gray'
k_linewidth=0.2
k_linestyle='-'
# --------------------------------------
c_plotter = curvature(main_linewidth, main_linestyle, E_zero_color, E_zero_linewidth, E_zero_linestyle, k_color, k_linewidth, k_linestyle)

def plot_curvature(
    path_read,
    axis=0,
    kticks=None,
    klabels=None,
    kbreaks=None,
    colors=['r', 'g', 'b'],
    curv_limit=None,
    ax=None,
    show=False,
    savefile=None):
    """
    A function to plot berry curvature along the path calculated with wannier90

    Parameters
    ----------
    path_read: str
        path where the pband file is, by default None
    axis: int
        axis along which berry curvature is calculated
    kticks: list
    
    klabels: list

    kbreaks: list

    colors: str
        color of each component of the curture. By default ['r', 'g', 'b']
    curv_limit: list
        interval between which curvature is plotted.
    ax: matplotlib.axes.Axes
        By default None
    show: bool
        wether figure is shown or not. By default False
    savefile: str
        path where output figure is saved. When None, it is not
        stored. By default 'None'
    """

    fig, ax = c_plotter.plot_curvature_path(path_read, axis, kticks, klabels, kbreaks, colors, curv_limit, ax)

    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile, bbox_inches='tight')
    
    return fig, ax