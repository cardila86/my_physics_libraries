from plotting_class import plottingTools
import matplotlib.pyplot as plt

#from plotting_class import styler
import time

t0 = time.time()

plotter = plottingTools()
fig, ax = plt.subplots()


path_read_vaspkit = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/NdPt2B/SG180/FM001/k131311/2-bands-SOC-2'
path_read_wannier90 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/NdPt2B/SG180/FM001/k11119/3-wannier-SOC'
ax = plotter.plot_bands(path_read_vaspkit,
                                    E_limit=[-1, 1],
                                    program='vaspkit',
                                    E_zero=8.1944,
                                    E_vaspkit=False,
                                    klabels=[r'$\Gamma$','M','K',r'$\Gamma$','A','L','H','A|L',
                                            'M|H', ''],
                                    kticks= None,
                                    kbreaks=[7, 8],
                                    label=None,
                                    color='k',
                                    ax=ax,
                                    show=False,
                                    savefile=None)
ax = plotter.plot_bands(path_read_wannier90,
                                    E_limit=[-1, 1],
                                    program='wannier90',
                                    E_zero=8.1944,
                                    E_vaspkit=False,
                                    root='wannier90',
                                    klabels=[r'$\Gamma$','M','K',r'$\Gamma$','A','L','H','A|L',
                                            'M|H', ''],
                                    kticks= None,
                                    kbreaks=[7, 8],
                                    label=None,
                                    color='r',
                                    ax=ax,
                                    show=True,
                                    savefile=None)

print('Elapsed time: ', time.time()-t0)
