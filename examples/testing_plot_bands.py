from plotting_class import plottingTools
#from plotting_class import styler
import time

t0 = time.time()

plotter = plottingTools()

fig, ax = plotter.plot_bands_vaspkit(path_read='/home/carlos/Documents/scriptsCalculations/projects/RPt2B/NdPt2B/FM001/k131311/2-bands-SOC',
                                    E_limit=[-1, 1],
                                    E_zero=8.1944,
                                    E_vaspkit=False,
                                    klabels=[r'$\Gamma$','M','K',r'$\Gamma$','A','L','H','A|L',
                                            'M|H', r'$K|\Gamma$','M\'','K\'',r'$\Gamma$', 'A\'',
                                            'L\'', 'H\'', 'A\'|L\'','M\'|H\'','K\''],
                                    kticks= None,
                                    kbreaks=None,
                                    label=None,
                                    color='k',
                                    ax=None,
                                    show=True,
                                    savefile=None)


#bands = styler(fig, ax)
# bands.bandstructure_style(savefile='./testing/bands_updw/bands')

print('Elapsed time: ', time.time()-t0)
