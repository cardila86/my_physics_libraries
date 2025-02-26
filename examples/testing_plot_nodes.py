from plotting_class import plottingTools
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import time

path_read='/home/carlos/Documents/scriptsCalculations/projects/RPt2B/NdPt2B/FM001/k131311/2-bands-SOC'

convergence = 1e-2
closestSite = 1
E_fermi = 8.1944


def center_bands(data, E_fermi):
    kpoints, bands = data
    for band in bands:
        band-=E_fermi

    return kpoints, bands

def getNodes_minAround(info1, Elimit=[-100, 100], pathSave=None):
    # This method takes the closest kpoint around a +- closestSite. It is designed for relative big convergence values
    kpoints, bands = info1

    knodes, enodes, dnodes = [], [], []
    ind = -closestSite-1
    # ==== Iteration over Info 1 ====
    for j in range(len(bands)-1):
        for i in range(len(kpoints)):
            e1 = bands[j][i]
            e2 = bands[j+1][i]
            k = kpoints[i]
            if abs(e2-e1)<convergence and (Elimit[0]<e1<Elimit[1] or Elimit[0]<e2<Elimit[1]):
                if abs(i-ind)<closestSite+1:
                    if abs(e2-e1)<dnodes[-1]:
                        knodes[-1] = k
                        enodes[-1] = e1
                        dnodes[-1] = abs(e2-e1)
                        ind = i

                else:
                #print('Nodo: ', k, e1, '. Delta E: ', abs(e2-e1))
                    knodes.append(k)
                    enodes.append(e1)
                    dnodes.append(abs(e2-e1))

                    ind = i

    if pathSave is not None:
        doc = np.column_stack((knodes, enodes))
        doc = np.column_stack((doc, dnodes))
        np.savetxt(pathSave + '/nodes.dat', tuple(doc))

    return knodes, enodes

def getNodes_local(info1, Elimit=[-100, 100], pathSave=None):
    # This method gets the local minimum. Better with small convergence values
    kpoints, bands = info1

    knodes, enodes, dnodes = [], [], []
    # ==== Iteration over Info 1 ====
    for i in range(len(kpoints)):
        k = kpoints[i]
        for j in range(len(bands)-1):
            e1 = bands[j][i]
            e2 = bands[j+1][i]

            if abs(e2-e1)<convergence and (Elimit[0]<e1<Elimit[1] or Elimit[0]<e2<Elimit[1]):
                print('Nodo: ', k, e1, '. Delta E: ', abs(e2-e1))
                knodes.append(k)
                enodes.append(e1)
                dnodes.append(abs(e2-e1))

    if pathSave is not None:
        doc = np.column_stack((knodes, enodes))
        doc = np.column_stack((doc, dnodes))
        np.savetxt(pathSave + '/nodes.dat', tuple(doc))

    return knodes, enodes


t0 = time.time()

plotter = plottingTools()

fig, ax = plotter.plot_bands_vaspkit(path_read,
                                    E_limit=[-1, 1],
                                    E_zero=E_fermi,
                                    E_vaspkit=False,
                                    klabels=[r'$\Gamma$','M','K',r'$\Gamma$','A','L','H','A|L',
                                            'M|H', r'$K|\Gamma$','M\'','K\'',r'$\Gamma$', 'A\'',
                                            'L\'', 'H\'', 'A\'|L\'','M\'|H\'','K\''],
                                    kticks= None,
                                    kbreaks=None,
                                    label=None,
                                    color='k',
                                    ax=None,
                                    show=False,
                                    savefile=None)

data = plotter._read_bands_vaspkit(path_read, fermi_vaspkit=False, klabels_bool=False, kticks_bool=False)
data = center_bands(data, E_fermi)
knodes, enodes = getNodes_minAround(data, Elimit=[-1, 1])
ax.scatter(knodes, enodes, c='b')

plt.show()

