from plotting_class import plottingTools
import matplotlib.pyplot as plt

#from plotting_class import styler
import time

path_read = '/home/carlos/Documents/scriptsCalculations/projects/kagome123/DFT/1-Model/d-orb/2D/opt/scf-soc/bands'

t0 = time.time()

plotter = plottingTools()
# fig, ax = plt.subplots()

# kpoints_0, E_0, orbs_projection_list = plotter._read_bands_vaspkit_projected(path_read, 'PBAND_Sc_SOC.dat', [0], fermi_vaspkit=False, klabels_bool=False, kticks_bool=False)
# kpoints_1, E_1 = plotter._read_bands_vaspkit(path_read, fermi_vaspkit=False, klabels_bool=False, kticks_bool=False)
# for i in range(len(E_0)):
#     if any(abs(E_0[i]-E_1[i]) > 0.1):
#         plt.plot(kpoints_0, E_0[i], 'r')
#         plt.plot(kpoints_1, E_1[i], 'b')
#         # print('Error: E_0[0]-E_1[0] > 0.1')
#         # dif = E_0[i]-E_1[i]
#         # mask = abs(dif) > 0.1
#         # print(dif)
#         # print(abs(E_0[0]-E_1[0]))
#         # print(i)
#         # exit()
# plt.show()
# exit()
fig, ax = plotter.plot_bands_projected_orbs(path_read,
                                        root='PBAND_Sc_SOC.dat',
                                        orbs=[0, [1, 2, 3], [4, 5, 6, 7, 8]],
                                        E_limit=[-1, 1],
                                        E_zero=-2.8161,
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

plt.show()