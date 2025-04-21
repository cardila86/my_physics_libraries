from my_physics_libraries.plotting_class import plottingTools
import matplotlib.pyplot as plt
#from plotting_class import styler
import time

path_read = '/home/ficomaco/Documents/work/CarlosArdila/projects/kagome123/DFT-Model/d-orb/Sc-atom/2D/opt/scf-soc/bands'

t0 = time.time()

plotter = plottingTools()



plotter.plot_bands_projected_orbs_p4vasp(path_read,
                                        root='PBANDS_Sc_s.dat',
                                        E_limit=[-13, 15],
                                        E_zero=-2.7766,
                                        klabels=[r'$\Gamma$', 'M', 'K', r'$\Gamma$'],
                                        colors=[[255, 0, 0], [0, 255, 0], [0, 0, 255]],
                                        )




# fig, ax = plotter.plot_bands_projected_orbs(path_read,
#                                         root='PBANDS_Sc_P4vasp.dat',
#                                         orbs=[0, [1, 2, 3]],
#                                         E_limit=[-13, 15],
#                                         E_zero=-2.7766,
#                                         klabels=[r'$\Gamma$', 'M', 'K', r'$\Gamma$'],
#                                         colors=[[255, 190, 11], [58, 134, 255], [255, 0, 110]],
                                        # )


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
