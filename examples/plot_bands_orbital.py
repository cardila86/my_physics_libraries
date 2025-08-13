import time
import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 1, figsize=[6, 4])

t0 = time.time()
import my_physics_libraries as mpl
from my_physics_libraries.styling import styler 

path_read='./scf-soc/bands'
fermi_path='./scf-soc'
# path_save='./spin_polarized_bands.png'
                  
for i, orbs in zip([0, 1, 2], [0, [1, 2, 3], [4, 5, 6, 7, 8]]):
   fig, axs[i] = mpl.plot_bands(path_read=path_read,
                  fermi_path=fermi_path,
                  code='vasp',
                  orbitals=[orbs],
                  mode='scatter',
                  E_limit=[-1, 1],
                  vmin=-1, vmax=1,
                  klabels=[r'$\Gamma$','M','K',r'$\Gamma$','A','L','H','A|L', 'M|H', 'K','','','', '','', '', '','',''],
                  kbreaks=[7, 8],
                  ax=axs[i],
                  )
   if i!=2:
      style = styler(fig, axs[i])
      style.set_ticks(erase_xticks=True)

print('running time: ', time.time()-t0)

# plt.savefig(path_save)
plt.show()