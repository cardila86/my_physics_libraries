__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

import matplotlib.pyplot as plt

from ..plotting import bands
from ..macrodata_refinement.utils import *
from ..macrodata_refinement.utils_vaspvis import *
from ..macrodata_refinement.utils_pyprocar import *

# ------ default style parameters ------
main_linewidth=1
main_linestyle='-'

background_linecolor='lightgray'

E_zero_color='gray'
E_zero_linewidth=0.2
E_zero_linestyle='--'

k_color='gray'
k_linewidth=0.2
k_linestyle='-'

mask=None
marker_size=1
cmap='hot_r'
vmin=None
vmax=None
marker='o'
ticks=None
# --------------------------------------
b_plotter = bands(main_linewidth, main_linestyle, E_zero_color, E_zero_linewidth, E_zero_linestyle, k_color, k_linewidth, k_linestyle, background_linecolor, marker_size)
r_bands = read()
# --------------------------------------

def plot_bands(
    path_read=None,
    code='vaspkit',
    E_zero=0,
    E_limit=None,
    color='k',
    nbands=None,
    spin=None,
    orbitals=None,
    atoms=None,
    root=None,

    klabels= None,
    kticks=None,
    kbreaks=None,
    label=None,
    cmap='bwr',
    vmin=None,
    vmax=None,

    mode='plain',
    fermi_path=None,
    E_vaspkit=False,
    cbar=True,
    legend=True,
    ax=None,
    show=False,
    savefile=None
    ):
    """
    A function to plot band structure

    Parameters
    ----------
    path_read: str
        path where the band file is, by default None
    code: str
        code used for generating the band file. Options are 
        'vasp', 'vaspkit', 'p4vasp', 'wannier90'. When choosing
        'p4vasp', it is supposed that input files contain pband.
        By default 'vaspkit'
    E_zero: float
        Set Fermi energy. By default 0
    E_limit: list
        Energy interval to plot. By default None
    color: str
        color of the band structure. When plotting pbbands, a list
        can be introduced to define the projections colors. By default 'k'
    nbands: float, list
        range of bands to plot. When float is introduced, it is shown
        the first 'nbands' bands. When list is introduced, it refers
        to the interval of bands. By default None
    spin: str
        spin to plot. Options are 'x', 'y', 'z', 'up', 'down'.
        By default None, which means no spin polarization plot.
    klabels: 
    kticks:
    kbreaks:
    label:
    fermi_path: str
        path to the file containing Fermi energy, when none is given,
        fermi energy is not read from any file, only E_zero parameter
        is taken into account. When given, both are considered. By
        default None
    
    E_vaspkit: bool
        Defines wheter Fermi was set at zero by vaspkit when 
        generating the files, or not. When True, FERMI_ENERGY
        file is read and shifted back to original. By default
        False
    
    ax: matplotlib.axes.Axes
        By default None
    show: bool
        wether figure is shown or not. By default False
    savefile: str
        path where output figure is saved. When None, it is not
        stored. By default 'None'
    """
    if code=='vaspkit' or  code=='wannier90':
        if spin is None and orbitals is None and atoms is None:
            # -------------- reads info ------------
            if code=='vaspkit':
                if klabels is None and kticks is None:
                    kpoints, E, klabels, kticks = r_bands._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=True, kticks_bool=True)
                elif klabels is None and kticks is not None:
                    kpoints, E, klabels = r_bands._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=True, kticks_bool=False)
                elif klabels is not None and kticks is None:
                    kpoints, E, kticks = r_bands._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=False, kticks_bool=True)
                else:
                    kpoints, E = r_bands._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=False, kticks_bool=False)   
            elif code=='wannier90':
                # path_read=path_read+'/'#+root
                if klabels is None and kticks is None:
                    kpoints, E, klabels, kticks = r_bands._read_bands_wannier90(path_read, klabels_bool=True, kticks_bool=True)
                elif klabels is None and kticks is not None:
                    kpoints, E, klabels = r_bands._read_bands_wannier90(path_read, klabels_bool=True, kticks_bool=False)
                elif klabels is not None and kticks is None:
                    kpoints, E, kticks = r_bands._read_bands_wannier90(path_read, klabels_bool=False, kticks_bool=True)
                else:
                    kpoints, E = r_bands._read_bands_wannier90(path_read, klabels_bool=False, kticks_bool=False)
            # -------------- plotting --------------
            fig, ax = b_plotter.plot_bands_processed(kpoints, E, E_zero, klabels, kticks, kbreaks, label, color, nbands, ax)
    elif code=='p4vasp':
        # -------------- reads info ------------
        if klabels is None and kticks is None:
            kpoints, bands, bands_projections, klabels, kticks = r_bands._read_pbands_p4vasp(path_read, root, klabels_bool=True, kticks_bool=True)
        elif klabels is None and kticks is not None:
            kpoints, bands, bands_projections, klabels = r_bands._read_pbands_p4vasp(path_read, root, klabels_bool=True, kticks_bool=False)
        elif klabels is not None and kticks is None:
            kpoints, bands, bands_projections, kticks = r_bands._read_pbands_p4vasp(path_read, root, klabels_bool=False, kticks_bool=True)
        else:
            kpoints, bands, bands_projections = r_bands._read_pbands_p4vasp(path_read, root, klabels_bool=False, kticks_bool=False)
        # -------------- plotting --------------
        if color is not list:
            color=[[255, 0, 0], [0, 255, 0], [0, 0, 255]]
        fig, ax = b_plotter.plot_pbands_p4vasp(kpoints, bands, bands_projections, E_zero, E_limit, klabels, kticks, kbreaks, color, label, nbands, ax)
    
    elif code=='vasp':
        if klabels is None and kticks is None:
            bands, kpoints, e_fermi_outcar, kticks, klabels, projection = r_bands._read_bands_vasp(path_read, fermi_path, spin, orbitals, atoms, klabels=True, kticks=True)
        elif klabels is not None and kticks is None:
            bands, kpoints, e_fermi_outcar, kticks, ___, projection = r_bands._read_bands_vasp(path_read, fermi_path, spin, orbitals, atoms, klabels=True, kticks=True)
        elif klabels is None and kticks is not None:
            bands, kpoints, e_fermi_outcar, ___, klabels, projection = r_bands._read_bands_vasp(path_read, fermi_path, spin, orbitals, atoms, klabels=True, kticks=True)
        elif klabels is not None and kticks is not None:
            bands, kpoints, e_fermi_outcar, ___, ____, projection = r_bands._read_bands_vasp(path_read, fermi_path, spin, orbitals, atoms, klabels=True, kticks=True)
            # bands, kpoints, e_fermi_outcar, kticks, klabels, orbitals_projections
        E_zero += e_fermi_outcar
        # --- check spinor without SOC---
        if projection is not None:
            if int(len(bands))==int(2*len(projection)):
                shape = list(projection.shape)
                projection_new = np.empty((shape[0]*2, shape[1], shape[2]))
                projection_new[0:shape[0]][:][:] = projection
                projection_new[shape[0]:][:][:] = projection
                projection = projection_new
        # ------------ plotting ------------
        if mode == 'plain':
            fig, ax = b_plotter.plot_bands_processed(kpoints,
                                bands,
                                E_zero,
                                klabels,
                                kticks,
                                kbreaks,
                                label,
                                color,
                                nbands,
                                ax)
        elif mode == 'parametric':
            fig, ax = b_plotter.plot_bands_parametric(kpoints,
                                bands,
                                projection,
                                E_zero,
                                klabels,
                                kticks,
                                kbreaks,
                                label,
                                color,
                                cmap,
                                nbands,
                                spin,
                                orbitals,
                                atoms,
                                vmin,
                                vmax,
                                cbar,
                                ax)
        elif mode == 'scatter':
            fig, ax = b_plotter.plot_bands_scatter(kpoints,
                                bands,
                                projection,
                                E_zero,
                                klabels,
                                kticks,
                                kbreaks,
                                label,
                                color,
                                cmap,
                                nbands,
                                spin,
                                orbitals,
                                atoms,
                                vmin,
                                vmax,
                                cbar,
                                ax)
        else:
            print("Invalid mode. Choose 'plain' or 'parametric'.")
            exit()
    #  --------------------- E_limit -----------------
    ax.set_ylabel(r'$E-E_{F} [eV]$')
    if E_limit is not None:
        ax.set_ylim(E_limit)
    if legend:
        lgnd = ax.legend(scatterpoints=1, fontsize=10, loc='upper right')
        # make all label markers the same size even when plotting different sizes
        for i in range(len(lgnd.legend_handles)):
            lgnd.legend_handles[i]._sizes =[20]
    if savefile is not None:
        plt.savefig(os.path.join(path_read, savefile), bbox_inches='tight')
    if show:
        plt.show()
    # --------------------------------------
    return fig, ax

def plot_surface(path_read=None,
    code='vaspk',
    E_zero=0,
    E_limit=None,
    color='k',
    nbands=None,
    spin=None,
    orbitals=None,
    atoms=None,
    kaxis=[0, 1],

    klabels= None,
    kticks=None,
    kbreaks=None,
    label=None,
    cmap='bwr',
    vmin=None,
    vmax=None,

    mode='plain',
    fermi_path=None,
    E_vaspkit=False,
    cbar=True,
    legend=True,
    ax=None,
    show=False,
    savefile=None
    ):
    """
    A function to plot energy surface

    Parameters
    ----------
    path_read: str
        path where the band file is, by default None
    code: str
        code used for generating the band file. Options are 
        'vasp', 'vaspkit', 'p4vasp', 'wannier90'. When choosing
        'p4vasp', it is supposed that input files contain pband.
        By default 'vaspkit'
    E_zero: float
        Set Fermi energy. By default 0
    E_limit: list
        Energy interval to plot. By default None
    color: str
        color of the band structure. When plotting pbbands, a list
        can be introduced to define the projections colors. By default 'k'
    nbands: float, list
        range of bands to plot. When float is introduced, it is shown
        the first 'nbands' bands. When list is introduced, it refers
        to the interval of bands. By default None
    spin: str
        spin to plot. Options are 'x', 'y', 'z', 'up', 'down'.
        By default None, which means no spin polarization plot.
    klabels: 
    kticks:
    kbreaks:
    label:
    fermi_path: str
        path to the file containing Fermi energy, when none is given,
        fermi energy is not read from any file, only E_zero parameter
        is taken into account. When given, both are considered. By
        default None
    
    E_vaspkit: bool
        Defines wheter Fermi was set at zero by vaspkit when 
        generating the files, or not. When True, FERMI_ENERGY
        file is read and shifted back to original. By default
        False
    
    ax: matplotlib.axes.Axes
        By default None
    show: bool
        wether figure is shown or not. By default False
    savefile: str
        path where output figure is saved. When None, it is not
        stored. By default 'None'
    """
    if code=='vasp':
        bands, kpoints, e_fermi_outcar, projection = r_bands._read_surface_vasp(path_read, fermi_path, spin, orbitals, atoms, klabels=True, kticks=True)

        kpoints = kpoints.T
        k1, k2 = kpoints[kaxis[0]], kpoints[kaxis[1]]
        kpoints = [k1, k2]
        
        E_zero += e_fermi_outcar
        #  --------------------- E_limit -----------------
        if E_limit is not None:
            for i in range(len(bands)):
                band = bands[i]
                mask = (band>=E_limit[0]) & (band<=E_limit[1])
                band[~mask] = np.nan
                bands[i] = band
        # --- check spinor without SOC---
        if projection is not None:
            if int(len(bands))==int(2*len(projection)):
                shape = list(projection.shape)
                projection_new = np.empty((shape[0]*2, shape[1], shape[2]))
                projection_new[0:shape[0]][:][:] = projection
                projection_new[shape[0]:][:][:] = projection
                projection = projection_new
        # ------------ plotting ------------
        fig, ax = b_plotter.plot_surface_scatter(
                            kpoints,
                            bands,
                            projection,
                            E_zero,
                            color,
                            nbands,
                            spin,
                            orbitals,
                            atoms,
                            vmin,
                            vmax,
                            cbar,
                            ax)
    #  -----------------------------------------------
    ax.set_zlabel(r'$E-E_{F} [eV]$')

    if savefile is not None:
        plt.savefig(os.path.join(path_read, savefile), bbox_inches='tight')
    if show:
        plt.show()
    # --------------------------------------
    return fig, ax

'''
# ---- chech for repaired PROCAR file ----
        if not os.path.exists(os.path.join(path_read, "PROCAR_repaired")):
            UtilsProcar.ProcarRepair(
                os.path.join(path_read, "PROCAR"),
                os.path.join(path_read, "PROCAR_repaired"))

        outcarparser = UtilsProcar()
        # --- read fermi from file ---
        e_fermi = 0
        if fermi_path is not None:
            e_fermi = outcarparser.FermiOutcar(os.path.join(fermi_path, "OUTCAR"))
        e_fermi += E_zero
        # ------- read procar --------
        recLat = None  # Will contain reciprocal vectors, if necessary
        outcarparser = UtilsProcar()
        recLat = outcarparser.RecLatOutcar(os.path.join(path_read, "OUTCAR"))

        file = ProcarParser()
        file.readFile(os.path.join(path_read, "PROCAR_repaired"), recLattice=recLat)
        data = ProcarSelect(file)
        if spin is not None:
            if spin=='x':
                data.selectIspin([1])
            elif spin=='y':
                data.selectIspin([2])
            elif spin=='z':
                data.selectIspin([3])
            elif spin=='up':
                data.selectIspin([0])
            elif spin=='down':
                data.selectIspin([1])
            else:
                raise ValueError("Invalid spin option. Choose 'x', 'y', 'z', 'up', or 'down'.")
        if orbitals is not None:
            pass
        else:
            orbitals = [-1]
        if atoms is not None:
            pass
        else:
            atoms = [-1]
        # ------- substract Fermi energy -----------
        data.bands = (data.bands.transpose() - np.array(e_fermi)).transpose()
        # ------- plot bands -----------
        plot = ProcarPlot(data.bands, data.spd, data.kpoints)
        if mode == "scatter":
            fig, ax =plot.scatterPlot(ax, mask=mask, size=size,
                            cmap=cmap, vmin=vmin,
                            vmax=vmax, marker=marker, ticks=ticks)
            # plt.colorbar()
            plt.ylabel(r"$E-E_F$ [eV]")

        elif mode == "plain":
            fig, ax =plot.plotBands(ax, size, marker=marker, ticks=ticks)
            plt.ylabel(r"$E-E_F$ [eV]")
        
        elif mode == "parametric":
            fig, ax =plot.parametricPlot(ax, cmap=cmap, vmin=vmin, vmax=vmax,
                                ticks=ticks)
            plt.ylabel(r"Energy [eV]")
            
        elif mode == "atomic":
            fig, ax =plot.atomicPlot(ax, cmap=cmap, vmin=vmin, vmax=vmax)
            plt.ylabel(r"Energy [eV]")
'''