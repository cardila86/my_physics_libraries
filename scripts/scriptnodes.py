__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

from ..plotting import bands
from ..macrodata_refinement.utils import read, process


def plot_nodes():
    pass

def separate_nodes(
    path_read,
    kpath=None,
    Elimit=None,
    ktol=0.001,
    savefile=None):
    '''
    A function to filter the nodes calculated using wannierTools

    Parameters
    ----------
    path_read: str
        path and wannierTools' output file , by default None
    kpath: list
        list of kpoints in direct units along which the nodes are searched.
        By default None
    Elimit: list
        Energy interval. By default None
    ktol: float
        Tolerance for the k-path location
    savefile
        path and file name where filtered nodes are saved. By default None
    '''

    r_nodes = read()
    p_nodes = process()
    # -------------- reads info --------------
    kvec_car, gap, E, kvec_dir = r_nodes._read_nodes_wanniertools(path_read)
    nodes_separated, E_separeted, Egap_separeted, path_separated = p_nodes.separate_nodes(kvec_car, gap, E, kvec_dir,kpath, Elimit, ktol, savefile)

    return nodes_separated, E_separeted, Egap_separeted, path_separated