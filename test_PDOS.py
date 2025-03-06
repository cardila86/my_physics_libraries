from my_physics_libraries.plotting_class import plottingTools
from my_physics_libraries.styling_class import styler

import numpy as np
import matplotlib.pyplot as plt

import time

path_read = '/home/ficomaco/Downloads'
root = 'exportDOSWCl3-orb-d.dat'

t0 = time.time()

plotter = plottingTools()



E, PDOS = plotter._read_PDOS_p4vasp_projected(path_read, root)

colors = ['gray', 
          'gray',
            [125/255, 45/255, 87/255]]
alphas = [0.5, 
          0.5,
            1]
linewidths = [0.5, 
              0.5,
                1]
labels = [None, None, 'd']

fig, axs = plt.subplots(1, 3)
for DOS, color, a, lw, orb in zip(PDOS, colors, alphas, linewidths, labels):
    axs[0].plot(E, DOS, color=color, alpha=a, linewidth=lw, label=orb)


style = styler(fig, axs[0])

style.set_labels(ylabel=r'DoS $[states/eV]$', xlabel='$E-E_{F} [eV]$', fontsize_ylabel=18)
# style.set_scale()
style.set_ticks(yticks_step=20, fontsize_yticks=14,
                fontsize_xticks=14)
style.set_size(figSize=[40, 8], figdpi=300)
style.set_lim(xlim=[-15, 15], ylim=[-80, 80])
style.set_legend(fontsize=20)

# plt.show()
plt.savefig(path_read+'/PDOS.png', bbox_inches='tight')