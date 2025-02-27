from my_physics_libraries.plotting_class import plottingTools
import matplotlib.pyplot as plt
#from plotting_class import styler
import time

path_read = '/home/ficomaco/Documents/work/CarlosArdila/projects/kagome123/DFT-Model/p-orb/B-atom/2D/opt/scf-soc/bands'

t0 = time.time()

plotter = plottingTools()

fig, ax = plotter.plot_bands_projected_orbs(path_read,
                                        root='PBAND_B_SOC.dat',
                                        orbs=[0, [1, 2, 3]],
                                        E_limit=[-13, 15],
                                        E_zero=-2.7766,
                                        klabels=[r'$\Gamma$', 'M', 'K', r'$\Gamma$'],
                                        colors=[[214, 40, 40], [252, 191, 73], [0, 0, 0]],
                                        )


# ax = plotter.plot_bands(path_read,
#                                     E_limit=[-1, 1],
#                                     program='vaspkit',
#                                     E_zero=-2.8161,
#                                     E_vaspkit=False,
#                                     # klabels=[r'$\Gamma$','M','K',r'$\Gamma$','A','L','H','A|L',
#                                     #         'M|H', r'$K'],
#                                     kticks= None,
#                                     kbreaks=None,
#                                     label=None,
#                                     color='k',
#                                     show=False,
#                                     savefile=None)

print('Elapsed time: ', time.time()-t0)

plt.show()