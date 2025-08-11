__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

import numpy as np
import os

from ..macrodata_refinement.utils_pyprocar import *
from pymatgen.io.vasp.inputs import Kpoints

class read:
    def __init__(self):
        pass

    def _read_bands_vasp(self,
        path_read,
        fermi_path,
        spin,
        orbitals,
        atoms,
        klabels,
        kticks):
        '''
        atoms and orbitals are not implemented yet.
        
        '''
        # ---- check for repaired PROCAR file ----
        if not os.path.exists(os.path.join(path_read, "PROCAR_repaired")):
            utils_procar = UtilsProcar()
            utils_procar.ProcarRepair(os.path.join(path_read, "PROCAR"),
                os.path.join(path_read, "PROCAR_repaired"))

        outcarparser = UtilsProcar()
        # ---- check for processed band file -----
        if not os.path.exists(os.path.join(path_read, "bands.npy")):
            band_file=False
        else:
            band_file=True
        if not os.path.exists(os.path.join(path_read, "kpoints.txt")):
            kpoints_file=False
        else:
            kpoints_file=True
        if not os.path.exists(os.path.join(path_read, "klabels_kticks.txt")):
            klabels_kticks_file=False
        else:
            klabels_kticks_file=True
        if not os.path.exists(os.path.join(path_read, "spin_projections.npy")):
            spins_file=False
        else:
            spins_file=True
        if not os.path.exists(os.path.join(path_read, "orbitals_projections.npy")):
            orbitals_file=False
        else:
            orbitals_file=True
        if not os.path.exists(os.path.join(path_read, "atoms_projections.npy")):
            atoms_file=False
        else:
            atoms_file=True
        # --- read fermi from file ---
        e_fermi_outcar = 0
        if fermi_path is not None:
            e_fermi_outcar = outcarparser.FermiOutcar(os.path.join(fermi_path, "OUTCAR"))
        # ------- read procar --------
        recLat = None  # Will contain reciprocal vectors, if necessary
        outcarparser = UtilsProcar()
        recLat = outcarparser.RecLatOutcar(os.path.join(path_read, "OUTCAR"))

        if not band_file or not kpoints_file or (not spins_file and spin is not None) or (not orbitals_file and orbitals is not None) or (not atoms_file and atoms is not None) or not klabels_kticks_file:
            file = ProcarParser()
            file.readFile(os.path.join(path_read, "PROCAR_repaired"), recLattice=recLat)
            data = ProcarSelect(file)
            # -- bands --
            bands = data.bands.transpose()
            np.save(os.path.join(path_read, "bands.npy"), bands)
            # -- kpoints --
            kpoints = data.kpoints
            # -- klabels or kticks --
            kpoints_pathfile = Kpoints.from_file(os.path.join(path_read, "KPOINTS"))

            klabels = []
            labels = [word for l in kpoints_pathfile.labels for word in l.split()]
            for k in labels:
                if not k.isnumeric():
                    klabels.append(k)

            klabels_new = []
            for i in range(len(klabels)):
                if i==0:
                    klabels_new.append(klabels[i])
                    pair = 1
                else:
                    if klabels_new[-1] != klabels[i]:
                        if pair==1:
                            klabels_new.append(klabels[i])
                            pair = 0
                        else:
                            klabels_new[-1] = (klabels_new[-1]+'|'+klabels[i])
                            pair = 1
                    else:
                        pair = 1
            klabels = klabels_new

            kstep = len(kpoints) // (len(klabels)-1)
            kticks = [kpoints[(kstep*i)-1] for i in range(len(klabels))]
            kticks[0] = kpoints[0]
            # -- kpoints --
            tool = tools()
            distances = []
            for i in range(len(kticks)-1):
                distances.append(tool.distance(kticks[i], kticks[i+1]))
            # -- klabels or kticks --
            kpoints_new = np.array([])
            x = 0
            for i in range(len(distances)):
                l = np.linspace(x, x+distances[i], kstep)
                x += distances[i]
                kpoints_new = np.concatenate((kpoints_new, l))
            kpoints = kpoints_new

            kticks = [kpoints[kstep*(i)-1] for i in range(len(klabels))]
            kticks[0] = kpoints[0]
            with open(os.path.join(path_read, "klabels_kticks.txt"), "w") as f:
                for l, t in zip(kticks, klabels):
                    f.write(f"{l}\t{t}\n")
            np.savetxt(os.path.join(path_read, "kpoints.txt"), kpoints)
            # -- projections --
            if spin is not None:
                spin_projections = np.transpose(file.spd[:, :, :, -1, -1], axes=(1, 0, 2))  # Does not  filter orbs neither atoms
                np.save(os.path.join(path_read, "spin_projections.npy"), spin_projections)
        else:
            print('---------------------------------------------')
            print('======= reading processed input files =======')
            print('---------------------------------------------')
            bands = np.load(os.path.join(path_read, "bands.npy"))
            kpoints = np.loadtxt(os.path.join(path_read, "kpoints.txt"))
            kpoints = [float(k) for k in kpoints]
            k_labels_ticks = np.loadtxt(os.path.join(path_read, "klabels_kticks.txt"), dtype=str)
            k_labels_ticks = k_labels_ticks.transpose()
            kticks = k_labels_ticks[0]
            klabels = k_labels_ticks[1]
            kticks = [float(k) for k in kticks]
            if spin is not None:
                spin_projections = np.load(os.path.join(path_read, "spin_projections.npy"))            
            if atoms is not None:
                atoms_projections = np.load(os.path.join(path_read, "atoms_projections.npy"))       
            if orbitals is not None:
                orbitals_projections = np.load(os.path.join(path_read, "orbitals_projections.npy"))

        # ------- select data --------
        if spin is not None:
            if spin=='x':
                spin = 1
            elif spin=='y':
                spin = 2
            elif spin=='z':
                spin = 3
            elif spin=='up':
                spin = 0
            elif spin=='down':
                spin = 1
            else:
                raise ValueError("Invalid spin option. Choose 'x', 'y', 'z', 'up', or 'down'.")
            
            spin_projections = spin_projections[:, :, spin]
            return bands, kpoints, e_fermi_outcar, kticks, klabels, spin_projections
        elif orbitals is not None:
            pass
        elif atoms is not None:
            pass
        else:
            return bands, kpoints, e_fermi_outcar, kticks, klabels, None
        
    def _read_bands_vaspkit(self,
        path_read,
        fermi_vaspkit=False,
        klabels_bool=False,
        kticks_bool=False):
        '''
        Reads the band structure generated by vaspkit1.3.5. It
        works for both spin and spinless calculations.

        Parameters
        ----------
        path_read: str
            path where the band file is, by default None
    
        returns
        ----------
        kpoints, E, klabel (optional), kticks (optional)
        '''
        # ------ single band structure file for both spin ------
        path_band = path_read + '/REFORMATTED_BAND.dat'
        if os.path.isfile(path_band):
            data = np.loadtxt(path_band, skiprows=1)
            data = data.transpose()
            kpoints = data[0]
            E = data[1:]
        # ------- band structure files splitted by spin --------
        else:
            path_band_dw = path_read + '/REFORMATTED_BAND_DW.dat'
            path_band_up = path_read + '/REFORMATTED_BAND_UP.dat'
            data = np.loadtxt(path_band_dw, skiprows=1)
            data = data.transpose()
            kpoints = data[0]
            E = data[1:]
            data = np.loadtxt(path_band_up, skiprows=1)
            data = data.transpose()
            E = np.append(E, data[1:], axis=0)
        # ------- Adjust fermi level from FERMI_ENERGY file ----
        if fermi_vaspkit:
            path_Efermi = path_read + '/FERMI_ENERGY'
            Efermi = np.loadtxt(path_Efermi, skiprows=1)
            for i in range(len(E)):
                E[i] += Efermi
        # ------- kticks & klabels --------
        if klabels_bool or kticks_bool:
            klabel, kticks = [], []
            with open(path_read + '/KLABELS') as data_file:
                for data in data_file:
                    if len(data.split())==2:
                        label = data.split()[0]
                        tick = eval(data.split()[1])
                        klabel.append(label)
                        kticks.append(tick)
            kticks = np.array(kticks)
        # ------- return data --------
        if klabels_bool and kticks_bool:
            return kpoints, E, klabel, kticks
        elif klabels_bool and not kticks_bool:
            return kpoints, E, klabel
        elif not klabels_bool and kticks_bool:
            return kpoints, E, kticks
        else:
            return kpoints, E

    def _read_bands_wannier90(self,
        path_read,
        klabels_bool=False,
        kticks_bool=False):
        '''
        Reads the band structure generated by wannier90. It
        works for both spin and spinless calculations.

        Parameters
        ----------
        path_read: str
            path where the band file is, by default None
    
        returns
        ----------
        kpoints, E, klabel (optional), kticks (optional)
        '''
        # ------- Define paths --------
        path_band = path_read+'_band.dat'
        path_label_info = path_read+'_band.labelinfo.dat'
        # ------- band structure in 2 columns separated by empty rows --------
        data = np.loadtxt(path_band, skiprows=0)
        data = data.transpose()
        kpoints, bands = data[0], data[1]
        for i in range(len(kpoints)-1):
            if kpoints[i+1]<kpoints[i]:
                npoints = i+1
                kpoints = kpoints[:npoints]
                break
        E = [bands[npoints*i:npoints*(i+1)] for i in range(len(bands)//npoints)]
        # ------- kticks & klabels --------
        if klabels_bool or kticks_bool:
            knames, kticks = [], []
            with open(path_label_info) as dataFile:
                idx = -10  # Numero lejos del primer idx
                for data in dataFile:
                    name = data.split()[0]
                    tick = eval(data.split()[2])

                    # Mira si los ticks estan uno al lado de otro.
                    # En caso tal, los junta en uno solo
                    if abs(eval(data.split()[1])-idx)<2:
                        knames[-1] = knames[-1]+'|'+name
                    else:
                        knames.append(name)
                        kticks.append(tick)

                    idx = eval(data.split()[1])
        # ------- return data --------
        if klabels_bool and kticks_bool:
            return kpoints, E, knames, kticks
        elif klabels_bool and not kticks_bool:
            return kpoints, E, knames
        elif not klabels_bool and kticks_bool:
            return kpoints, E, kticks
        else:
            return kpoints, E

    def _read_pbands_p4vasp(self,
        path_read,
        klabels_bool,
        kticks_bool):
        '''
        Reads the projected band structure generated by p4vasp.

        Parameters
        ----------
        path_read: str
            path where the pband file is, by default None
    
        returns
        ----------
        kpoints, bands, bands_projection, klabel (optional), kticks (optional)
        '''
        # ----------- separate bands and projections -----------
        nlines, nkpoints, nbands = 0, 0, 0
        with open(path_read) as dataFile:
            for data in dataFile:
                if len(data.split())==3:

                    break
                elif len(data.split())==0:
                    if nkpoints==0:
                        nkpoints = nlines
                    nbands+=1
                    nlines+=1
                else:
                    nlines+=1
        kpoints, bands, bands_projection = [], [], []
        if os.path.isfile(path_read):
            # -------- bands, kpoints and projection --------
            data_projection = np.loadtxt(path_read, usecols=(0,1,2), skiprows=nlines)
            data_projection = data_projection.transpose()
            kpoints = data_projection[0][:nkpoints]*2*np.pi
            bands = data_projection[1]
            projection = data_projection[2]
            bands = [bands[nkpoints*j:nkpoints*(j+1)] for j in range(nbands)]
            # ------ checks whether the bands are duplicated ------
            spin = False
            for i in range(len(bands)-1):
                if bands[i][0]>bands[i+1][0]:
                    spin = True
            n = 0
            for i in range((len(projection)//nkpoints)//nbands):
                bands_projection.append([projection[nkpoints*(j+nbands*i):nkpoints*(j+1+nbands*i)] for j in range(nbands)])
            # ======= check that projection and bands match =======
            for i in range(len(bands)):
                assert len(kpoints) == len(bands[i]), f'ERROR: The number of kpoints and the lenght of each band do not match. (1).'
            for i in range(len(bands_projection)):
                for j in range(len(bands_projection[i])):
                    assert len(kpoints) == len(bands_projection[i][j]), 'ERROR: The number of kpoints in the bands and the projection file do not match. (2)'
            # =====================================================
        else:
            print(f'ERROR: path \'{path_read}\' not found')
        # ------- kticks & klabels --------
        if klabels_bool or kticks_bool:
            klabel, kticks = [], []
            with open(path_read + '/KLABELS') as data_file:
                for data in data_file:
                    if len(data.split())==2:
                        label = data.split()[0]
                        tick = eval(data.split()[1])
                        klabel.append(label)
                        kticks.append(tick)
            kticks = np.array(kticks)
        # ------- return data --------
        if klabels_bool and kticks_bool:
            return kpoints, bands, bands_projection, klabel, kticks
        elif klabels_bool and not kticks_bool:
            return kpoints, bands, bands_projection, klabel
        elif not klabels_bool and kticks_bool:
            return kpoints, bands, bands_projection, kticks
        else:
            return kpoints, bands, bands_projection

    def _read_ldos_p4vasp(self,
        path_read):
        '''
        Reads the projected DOS generated by p4vasp.

        Parameters
        ----------
        path_read: str
            path where the band file is, by default None
    
        returns
        ----------
        E, PDOS, DOS
        '''
        if os.path.isfile(path_read):
            # -------- bands witoout projections --------
            # data_bands = np.loadtxt(path_projection, usecols=(0,1))
            # data_bands = data_bands.transpose()
            # kpoints = data_bands[0]
            # bands = data_bands[1]
            # kpoints = kpoints[:nkpoints]*2*np.pi
            # bands = [bands[nkpoints*j:nkpoints*(j+1)] for j in range(len(bands)//nkpoints)]
    
            data_projection = np.loadtxt(path_read, usecols=(0,1))
            data_projection = data_projection.transpose()
            E = data_projection[0]
            PDOS = data_projection[1]
            npoints = 0
            for i in range(len(E)-1):
                if E[i+1]<E[i]:
                    npoints = i+1
                    break
            E = E[:npoints]
            PDOS = [PDOS[npoints*j:npoints*(j+1)] for j in range(len(PDOS)//npoints)]
            DOS  = PDOS[0]
            PDOS = PDOS[1:]
        
            return E, PDOS, DOS
        else:
            print(f'ERROR: File \'{path_read}\' not found.')

    def _read_ldos_vaspkit(self,
        path_read,
        fermi_vaspkit=False):
        '''
        Reads the local density of states generated by vaspkit1.3.5.

        Parameters
        ----------
        path_read: str
            path where the band file is, by default None
    
        returns
        ----------
        E, orbitals_labels, orbitals, total
        '''        
        if os.path.isfile(path_read):
            data = np.loadtxt(path_read)
            data = data.transpose()
            # IMPORTANTE: Parece que el usar este metodo de sort no incrementa
            # tanto el tiempo, pero podriamos probar usando el sort de numpy.
            data = list(zip(*sorted(zip(*data))))
            data = np.array(data)
            E = data[0]
            orbitals = data[1:-1]
            total = data[-1]
            # 
            with open(path_read) as f:
                orbitals_labels = f.readline()
                orbitals_labels = orbitals_labels.split()[1:-1]
        else:
            print(f'File {path_read} not found.')
            exit()
        # ------- Adjust fermi level from FERMI_ENERGY file ----
        if fermi_vaspkit:
            path_Efermi = path_read + '/FERMI_ENERGY'
            Efermi = np.loadtxt(path_Efermi, skiprows=1)
            E += Efermi
        
        
        return E, orbitals_labels, orbitals, total

    def _read_nodes_wanniertools(self,
        path_read):
        '''
        Reads the nodes calculated by wanniertools.

        Parameters
        ----------
        path_read: str
            path where the band file is, by default None
    
        returns
        ----------
        kvec_car, gap, E, kvec_dir
        '''        
        # ------------ reading of the nodes ----------------
        path_nodes = path_read + '/Nodes.dat'
        if os.path.isfile(path_nodes):
            data = np.loadtxt(path_nodes)
            data = data.transpose()
            
            kvec_car= np.array(data[0:3]).transpose()
            gap     = np.array(data[3])
            E       = np.array(data[4])
            kvec_dir= np.array(data[5:]).transpose()
        # ------- path does not match with any file --------
        else:
            print(f'ERROR: path \'{path_nodes}\' not found')
        
        return kvec_car, gap, E, kvec_dir

    def _read_curvature_wannier90(self,
        path_read,
        axis):
        '''
        Reads the nodes calculated by wanniertools.

        Parameters
        ----------
        path_read: str
            path where the band file is, by default None
        axis: int
            axis along which is projected the curvature. 0 is
            total curvature, 1 x, 2 y, and 3 z.
            
        returns
        ----------
        kpoints, curv
        '''        
        data = np.loadtxt(path_read + '/wannier90-curv.dat')
        kpoints=data[:,0]
        if axis == 0:
            curv=np.zeros(len(kpoints))
            for i in [1, 2, 3]:
                curv+=np.array(data[:,i])**2
            curv = np.sqrt(curv)
            curv = curv.tolist()
        else:
            curv=data[:,axis]

        return kpoints, curv

    def _read_ahc_wannierberri(self,
        path_read,
        root,
        iteration):
        # --------- read data ----------
        ahc_data = [np.loadtxt(path_read+'/'+root+f"-ahc_iter-{i:04d}.dat") for i in iteration]

        return ahc_data

class process:
    def __init__(self):
        pass

    def separate_nodes(self,
        kvec_car,
        gap,
        E,
        kvec_dir,
        kpath,
        Elimit,
        ktol,
        savefile):
         # --- separates nodes by location and energy ---
        nodes_separated = []
        E_separeted     = []
        Egap_separeted  = []
        path_separated  = []
        if kpath is not None:
            for i in range(len(kpath)-1):
                nodes_separated.append([])
                E_separeted.append([])
                Egap_separeted.append([])

                symm_point_1 = kpath[i]
                symm_point_2 = kpath[i+1]
                vec_path = np.array(symm_point_2) - np.array(symm_point_1)
                for j in range(len(kvec_dir)):
                    vec_a = np.array(kvec_dir[j])-np.array(symm_point_1)
                    distance = np.linalg.norm(np.cross(vec_path, vec_a))/np.linalg.norm(vec_path)
                    if distance<ktol:
                        if Elimit is not None:
                            if Elimit[0]<=E[j]<=Elimit[1]:
                                nodes_separated[i].append(kvec_dir[j])
                                E_separeted[i].append(E[j])
                                Egap_separeted[i].append(gap[j])
                                path_separated.append(i+1)
                        elif Elimit is None:
                            nodes_separated[i].append(kvec_dir[j])
                            E_separeted[i].append(E[j])
                            Egap_separeted[i].append(gap[j])
                            path_separated.append(i+1)
        else:
            if Elimit is not None:
                for i in range(len(E)):
                    if Elimit[0]<=E[i]<=Elimit[1]:
                        nodes_separated.append(kvec_dir[i])
                        E_separeted.append(E[i])
                        Egap_separeted.append(gap[i])
                        path_separated.append(0)
            else:
                nodes_separated = kvec_dir
                E_separeted = E
                Egap_separeted = gap
                path_separated = np.zeros(len(E))
        # --------------------------------------
        if savefile is not None:
            knodes_saveformat   = []
            E_saveformat        = []
            Egap_saveformat     = []
            if kpath is not None:
                for i in range(len(nodes_separated)):
                    for j in range(len(nodes_separated[i])):
                        knodes_saveformat.append(nodes_separated[i][j])
                        E_saveformat.append(E_separeted[i][j])
                        Egap_saveformat.append(Egap_separeted[i][j])
            else:
                knodes_saveformat   = nodes_separated
                E_saveformat        = E_separeted
                Egap_saveformat     = Egap_separeted
            knodes_saveformat = np.array(knodes_saveformat)
            knodes_saveformat = knodes_saveformat.transpose()
            if len(knodes_saveformat)==0:
                knodes_saveformat = [[], [], []]

            doc = np.column_stack([path_separated, knodes_saveformat[0]])
            doc = np.column_stack([doc, knodes_saveformat[1]])
            doc = np.column_stack([doc, knodes_saveformat[2]])
            doc = np.column_stack([doc, E_saveformat])
            doc = np.column_stack([doc, Egap_saveformat])
            np.savetxt(savefile, doc, header='path      k1           k2           k3           E           Egap', 
                        fmt=['%4d', ' % 2.8f', ' % 2.8f', ' % 2.8f', ' % 2.8f', ' % 2.8f'])
        # --------------------------------------

        return nodes_separated, E_separeted, Egap_separeted, path_separated

class tools:
    def __init__(self):
        pass
    
    def distance(self, p1, p2):
        if type(p1) is not np.ndarray:
            p1 = np.array(p1)
        if type(p2) is not np.ndarray:
            p2 = np.array(p2)

        squared_dist = np.sum((p1-p2)**2, axis=0)
        dist = np.sqrt(squared_dist)
        return dist