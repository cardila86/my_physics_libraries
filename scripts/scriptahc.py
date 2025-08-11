__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

import matplotlib.pyplot as plt

from ..plotting import ahc
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
ahc_plotter = ahc(main_linewidth, main_linestyle, E_zero_color, E_zero_linewidth, E_zero_linestyle, k_color, k_linewidth, k_linestyle)
r_ahc = read()


def plot_ahc(
    path_read=None,
    root='AHC',
    code='wannierberri',
    ahc_axis='',
    spin='both',
    iteration=[1],
    E_zero=0,
    norm_ahc=0,
    E_limit=None,
    ahc_limit=None,
    colors=[[255, 0, 0], [0, 255, 0], [0, 0, 255]],
    nbands=None,

    klabels= None,
    kticks=None,
    kbreaks=None,
    label=None,

    fig_orientation='vertical',
    scale='linear',
    ax=None,
    show=False,
    savefile=None):
    if code=='wannierberri':
        # ---------------- alert ----------------
        print('------------------------------------------')
        print('WARNING: wannierberri changes the UNITS of\nAHC depending on the version used to run\nthe calculation. Please check them.')
        print('------------------------------------------')
        # --------- read data ----------
        ahc_data = r_ahc._read_ahc_wannierberri(path_read, root, iteration)
         # --------- plot fig ----------
        fig, ax = ahc_plotter.plot_AHC_wannierberri(ahc_data, ahc_axis, E_limit, ahc_limit, E_zero, norm_ahc, colors, spin, iteration, fig_orientation, scale, ax)

    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile, bbox_inches='tight')

    return fig, ax