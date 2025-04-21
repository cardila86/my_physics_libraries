from plotting_class import plottingTools
import time

path_read = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/NdPt2B/FM001/k11119/AHC/try3'   

t0 = time.time()

plotter = plottingTools()

fig, ax = plotter.plot_AHC_wannierberri(path_read,
                                        ahc_axis='all',
                                        E_limit=None,
                                        ahc_limit=None,
                                        E_zero=0,
                                        colors=['r', 'r', 'b', 'b', 'g', 'g'],
                                        spin='both',
                                        root='NdPt2B_FM001',
                                        iteration=[1],
                                        fig_orientation='vertical',
                                        scale='linear',
                                        ax=None,
                                        show=True,
                                        savefile=None)


print('Elapsed time: ', time.time()-t0)