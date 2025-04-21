from plottingClass import plottingTools
from plottingClass import styler
import time


path_read = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/1-opt/scf-soc/'

t0 = time.time()

plotter = plottingTools()

# data = plotter._read_LDoS(path_read)

fig, ax = plotter.plot_LDoS(
        path_read,
        roots=['PDOS_B_SOC.dat'],
        orbitals_tag=['s', 'py', 'pz', 'dz2', 'dxz'],
        # colors=['r', 'b'],
        E_limit=None,
        E_zero=0,
        E_vaspkit=False,
        ax=None,
        show=False,
        savefile=path_read+'PDOS_B_SOC.png')

print('Elapsed time: ', time.time()-t0)
