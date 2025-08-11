__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

import matplotlib.pyplot as plt

from ..plotting import dos
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
dos_plotter = dos(main_linewidth, main_linestyle, E_zero_color, E_zero_linewidth, E_zero_linestyle, k_color, k_linewidth, k_linestyle)
r_dos = read()



def plot_dos():
    pass

def plot_ldos(path_read,
    code='vaspkit',
    root=None,
    orbitals_tag='all',
    colors=None,
    E_limit=None,
    dos_limit=None,
    E_zero=0,
    E_vaspkit=False,
    fig_orientation='vertical',
    ax=None,
    show=False,
    savefile=None):
    """
    A function to plot local density of states

    Parameters
    ----------
    path_read: str
        path where the band file is, by default None
    code: str
        code used for generating the band file. Options are 
        'vasp', 'vaspkit', 'p4vasp', 'wannier90'. By default 'vaspkit'
    root: str
        root of the file where ldos is saved. By default None
    orbitals_tag: str
        orbitals to plot. By default 'all'
    colors: list
        list of colors to use when plotting the ldos. By default None
    E_limit: list
        Energy interval to plot. By default None
    dos_limit: list
        DOS interval to plot. By default None
    E_zero: float
        Set Fermi energy. By default 0
    E_vaspkit: bool
        Defines wheter Fermi was set at zero by vaspkit when 
        generating the files, or not. By default False
    fig_ortientation: str
        Orientation of the figure. By default 'vertical'
    
    ax: matplotlib.axes.Axes
        By default None
    show: bool
        wether figure is shown or not. By default False
    savefile: str
        path where output figure is saved. When None, it is not
        stored. By default 'None'
    """
    # --------- read data ---------
    path = f'{path_read}/{root}'
    if code=='vaspkit':
        E, orbitals_labels, orbitals, total = r_dos._read_LDoS_vaspkit(path, E_vaspkit)
    elif code=='p4vasp':
        E, orbitals, total = r_dos._read_LDoS_p4vasp(path, E_vaspkit)
    # ------- plotting data -------
    fig, ax = dos_plotter.plot_ldos_vaspkit(code, E, orbitals, total, orbitals_labels, orbitals_tag, colors, E_limit, dos_limit, E_zero, fig_orientation, ax)

    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile, bbox_inches='tight')

    return fig, ax