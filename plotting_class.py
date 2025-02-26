__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "January 06, 2025"

import matplotlib.pyplot as plt
import os
import numpy as np


class plottingTools:
    '''
    Class to plot physics graph.
    ==========
    Parameters
    ==========
        - property:
        - 
        -

    ======
    Return
    ======
        None.
    ======
    Tested compatibility
    ======
        - matplotlib 3.10.0
        - numpy 2.2.1
        - pyqt5
    '''
    def __init__(self, E_zero_color='gray', E_zero_linewidth=0.2, E_zero_linestyle='--',
                k_color='gray', k_linewidth=0.2, k_linestyle='-',
                main_linewidth=1.3, main_linestyle='-'):
        self.E_zero_color=E_zero_color
        self.E_zero_linewidth=E_zero_linewidth
        self.E_zero_linestyle=E_zero_linestyle
        self.k_color=k_color
        self.k_linewidth=k_linewidth
        self.k_linestyle=k_linestyle
        self.main_linewidth=main_linewidth
        self.main_linestyle=main_linestyle

    def _read_bands_vaspkit(self, path_read, fermi_vaspkit=False, klabels_bool=False, kticks_bool=False):
        '''
        Reads the band structure generated by vaspkit1.3.5. It
        works for both spin and spinless calculations.

        ==========
        Parameters
        ==========
            - property:
            - 
            -

        ======
        Return
        ======
            None.
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

    def _read_bands_wannier(self, path_read, klabels_bool=False, kticks_bool=False):
        '''
        Reads the band structure generated by wannier90.
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

    def __plot_bands_vaspkit_break(self,
                            path_read,
                            E_limit=[-1, 1],
                            E_zero=0,
                            E_vaspkit=False,
                            klabels= None,
                            kticks=None,
                            kbreaks=None,
                            label=None,
                            color='k',
                            ax=None,
                            show=False,
                            savefile=None):
        pass

    def plot_bands_vasp_orbs_atoms(self,   
                                   path_read,
                                   atoms=None,
                                   orbitals=None,
                                   spins=None,
                                   E_limit=[-1, 1],
                                   E_zero=0,
                                   klabels= None,
                                   kticks=None,
                                   cmap='jet',
                                   ax=None,
                                   show=False,
                                   savefile=None):
        '''
        Uses pyprocar 6.3.2
        '''
        
        import pyprocar
        

        bands = pyprocar.bandsplot(
                code='vasp', dirname=path_read, mode='parametric',
                fermi=E_zero, fermi_color = 'black', fermi_linestyle='dashed', fermi_linewidth=0.5,
                elimit=E_limit, cmap=cmap, linestyle=['solid', 'solid'], linewidth=[1, 1],
                kticks=kticks, knames=klabels,
                atoms=atoms, orbitals=orbitals, spins=spins,
                show=show, savefig=savefile, dpi=300, ax=ax)
        
        fig, ax = bands[0:2]
        return fig, ax    

    def plot_bands(self,
                           path_read,
                           program='vaspkit',
                           root='wannier90',
                           E_limit=[-1, 1],
                           E_zero=0,
                           E_vaspkit=False,
                           klabels= None,
                           kticks=None,
                           kbreaks=None,
                           label=None,
                           color='k',
                           nbands=None,
                           ax=None,
                           show=False,
                           savefile=None):
        # -------------- reads info --------------
        if program=='vaspkit':
            if klabels is None and kticks is None:
                kpoints, E, klabels, kticks = self._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=True, kticks_bool=True)
            elif klabels is None and kticks is not None:
                kpoints, E, klabels = self._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=True, kticks_bool=False)
            elif klabels is not None and kticks is None:
                kpoints, E, kticks = self._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=False, kticks_bool=True)
            else:
                kpoints, E = self._read_bands_vaspkit(path_read, fermi_vaspkit=E_vaspkit, klabels_bool=False, kticks_bool=False)
        elif program=='wannier90':
            path_read=path_read+'/'+root

            if klabels is None and kticks is None:
                kpoints, E, klabels, kticks = self._read_bands_wannier(path_read, klabels_bool=True, kticks_bool=True)
            elif klabels is None and kticks is not None:
                kpoints, E, klabels = self._read_bands_wannier(path_read, klabels_bool=True, kticks_bool=False)
            elif klabels is not None and kticks is None:
                kpoints, E, kticks = self._read_bands_wannier(path_read, klabels_bool=False, kticks_bool=True)
            else:
                kpoints, E = self._read_bands_wannier(path_read, klabels_bool=False, kticks_bool=False)
        # ------------ ax, fig objects -----------
        if ax is None and kbreaks is None:
            fig, ax = plt.subplots()
        elif ax is None and kbreaks is not None:
            num_ax = len(kbreaks)+1
            fig, ax = plt.subplots(num_ax)
        else:
            fig = None
            # IMPORTANTE: CREO QUE ES MEJOR GRAFICAR COMO SI FUERAN DISTINTOS EJES
            # referenceTicks = ax.get_xticks()
            # kpoints, kTicks = fixKpath(referenceTicks, kpoints, kTicks)
        # ------ plot klabels and kticks ------
        for ktick in kticks:
            ax.axvline(ktick, color=self.k_color, linewidth=self.k_linewidth, linestyle=self.k_linestyle)
        ax.set_xticks(kticks)
        ax.set_xticklabels(klabels)

        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # -------------- plot bands --------------
        if nbands is None:
            bands=E
        elif type(nbands) is int or type(nbands) is float:
            nbands = int(nbands)
            bands=E[0:nbands]
        elif type(nbands) is list:
            bands=[]
            for i in nbands:
                bands.append(E[int(i)])
        else:
            print('ERROR: nbands must be an integer, a float or a list of integers.\n'+
                  'a '+str(type(nbands))+' type was recieved. Please check inputs.') 
            exit()
        n = 0
        for band in bands:
            band-=E_zero
            if n==0 and label is not None:
                n+=1
                ax.plot(kpoints, band, c=color, linewidth=self.main_linewidth, linestyle=self.main_linestyle, label=label)
            else:
                ax.plot(kpoints, band, c=color, linewidth=self.main_linewidth, linestyle=self.main_linestyle)

        # ------------- set limits --------------
        ax.set_xlim([kticks[0], kticks[-1]])
        ax.set_ylim(E_limit)
        ax.set_ylabel(r'$E-E_{F} [eV]$')

        # ------ output -------
        if label is not None:
            ax.legend()
        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.
        if show:
            plt.show()
        if savefile is not None:
            plt.savefig(savefile, bbox_inches='tight')
        if fig is None:
            return ax
        else:
            return fig, ax

    def plot_pathIntensity(self,
                           path_read,
                           bands,
                           intensity,
                           E_limit=[-1, 1],
                           E_zero=0,
                           E_vaspkit=False,
                           klabels= None,
                           kticks=None,
                           kbreaks=None,
                           label=None,
                           color='k',
                           ax=None,
                           show=False,
                           savefile=None
                           ):
        '''
        IN DEVELOPEMENT
        '''
        # ------- some parameters -------
        kwidth = 0.2 
        # ======== read data ========
        bands = np.loadtxt(open(path_read + "/bandsWannierBerri.dat","r"))
        curv = np.loadtxt(open(path_read + '/curvatureWannierBerri.dat',"r"))
        kticksAndLabel = np.loadtxt(open(path_read + '/kticks.dat',"r"), str)

        kpoints = np.transpose(bands)[0]
        bands = np.transpose(bands)[1:]
        curv = np.transpose(curv)[1:]
        kticks = np.array(([float(i) for i in np.transpose(kticksAndLabel)[0]]))
        if klabels==[]:
            klabels = np.transpose(kticksAndLabel)[1]
        # ====== organize data ======
                           
    def plot_AHC_wannierberri(self,
                           path_read,
                           ahc_axis='all',
                           E_limit=None,
                           ahc_limit=None,
                           E_zero=0,
                           colors=None,
                           spin='both',
                           root='AHC',
                           iteration=[1],
                           fig_orientation='vertical',
                           scale='linear',
                           ax=None,
                           show=False,
                           savefile=None):
        # ---------------- alert ----------------
        print('------------------------------------------')
        print('WARNING: wannierberri changes the UNITS of\nAHC depending on the version used to run\nthe calculation. Please check them.')
        print('------------------------------------------')
        # ----------- fig, ax objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = None
        if colors is not None:
            ax.set_prop_cycle('color', colors)  # Colores a usar
        # --------- read data ----------
        ahc_data = [np.loadtxt(path_read+'/'+root+f"-ahc_iter-{i:04d}.dat") for i in iteration]
        # ----- ahc axis selection -----
        if ahc_axis=='all':
            ahc_label =fr'$\sigma [S/cm]$'
            if spin=='both':
                curv_label = [r'$\sigma_{yz}$ - down', r'$\sigma_{yz}$ - up', r'$\sigma_{xz}$ - down',
                              r'$\sigma_{xz}$ - up', r'$\sigma_{xy}$ - down', r'$\sigma_{xy}$ - up']
            if spin=='up':
                curv_label = [r'$\sigma_{yz}$ - up', r'$\sigma_{xz}$ - up', r'$\sigma_{xy}$ - up']
            if spin=='down':
                curv_label = [r'$\sigma_{yz}$ - down', r'$\sigma_{xz}$ - down', r'$\sigma_{xy}$ - down']
        else:    
            ahc_label=fr'$\sigma_{ahc_axis} [S/cm]$'
        # --------- orientation ---------
        if (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
            ax.set_xlabel(ahc_label)
            if fig_orientation=='vertical_rigth'or'vr':
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
            ax.set_ylabel(fr'$E-E_F$ [eV]')
            ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
            ax.set_xscale(scale)
        elif (fig_orientation=='horizontal'or'h'):
            ax.set_ylabel(ahc_label)
            ax.set_xlabel(fr'$E-E_F$ [eV]')
            ax.axvline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
            ax.set_yscale(scale)
        # --------- plotting ---------
        if ahc_axis=='all':
            n=0
            for x in range(3):
                for i in range(len(ahc_data)):
                    a=ahc_data[i]
                    E=a[:,0]
                    E-=E_zero
                    AHC_down=a[:,x+1]
                    AHC_up  =a[:,x+1+3]
                    # ----- curv -----
                    if (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
                        if spin=='both':
                            ax.plot(AHC_down, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                            n+=1
                            ax.plot(AHC_up, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='up':
                            ax.plot(AHC_up, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='down':
                            ax.plot(AHC_down, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)  
                    elif (fig_orientation=='horizontal'or'h'):
                        if spin=='both':
                            ax.plot(E, AHC_down, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                            n+=1
                            ax.plot(E, AHC_up, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='up':
                            ax.plot(E, AHC_up, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='down':
                            ax.plot(E, AHC_down, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)  
        else:
            for i in range(len(ahc_data)):
                a=ahc_data[i]
                E=a[:,0]
                E-=E_zero
                AHC_down=a[:,x+1]
                AHC_up  =a[:,x+1+3]
                # ----- curv -----
                if (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
                    if spin=='both':
                        ax.plot(AHC_down, E, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        ax.plot(AHC_up, E, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                    elif spin=='up':
                        ax.plot(AHC_up, E, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                    elif spin=='down':
                        ax.plot(AHC_down, E, linewidth=self.main_linewidth, linestyle=self.main_linestyle)  
                elif (fig_orientation=='horizontal'or'h'):
                    if spin=='both':
                        ax.plot(E, AHC_down, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        ax.plot(E, AHC_up, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                    elif spin=='up':
                        ax.plot(E, AHC_up, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                    elif spin=='down':
                        ax.plot(E, AHC_down, linewidth=self.main_linewidth, linestyle=self.main_linestyle)  
        # -------- legend ---------
        ax.legend(loc='best')
        # -------- set limits --------
        if ahc_limit is not None and (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
            ax.set_xlim(ahc_limit)
        elif ahc_limit is not None and (fig_orientation=='horizontal'or'h'):
            ax.set_ylim(ahc_limit)
        
        if E_limit is not None and (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
            ax.set_ylim(E_limit)
        elif E_limit is None and (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
            ax.set_ylim([min(E), max(E)])
        elif E_limit is not None and (fig_orientation=='horizontal'or'h'):
            ax.set_xlim(E_limit)
        elif E_limit is None and (fig_orientation=='horizontal'or'h'):
            ax.set_xlim([min(E), max(E)])

        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.
        if show:
            plt.show()
        if savefile is not None:
            plt.savefig(savefile, bbox_inches='tight')
        if fig is None:
            return ax
        else:
            return fig, ax    
    
    def plot_xy():
        # look script 'plot_curvature.py' or 'plot_AHC.py'
        pass

    def plot_grid():
        # For plotting Curvature on surface. 
        pass

    def _read_LDoS(self, path_read, fermi_vaspkit=False):
        '''
        Reads the local density of states generated by vaspkit1.3.5.
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
    
    def plot_LDoS(self,
                  path_read,
                  roots=None,
                  orbitals_tag='all',
                  colors=None,
                  E_limit=None,
                  E_zero=0,
                  E_vaspkit=False,
                  ax=None,
                  show=False,
                  savefile=None):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = None
        if colors is not None:
            ax.set_prop_cycle('color', colors)  # Colores a usar
        # --------- read data ---------
        E, orbitals_labels, orbitals, total = self._read_LDoS(path_read+roots[0], E_vaspkit)
        # ----- Sets zero energy level -----
        E -= E_zero
        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # --------- plot partial DOS ---------
        if orbitals_tag=='all':
            for orbs, label in zip(orbitals, orbitals_labels):
                ax.plot(E, orbs, label=label)
        else:
            for tag in orbitals_tag:
                index = orbitals_labels.index(tag)
                ax.plot(E, orbitals[index], label=tag)
        # --------- plot total DOS ---------
        ax.plot(E, total, color='gray', linewidth=1, linestyle='--')

        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('Density of States (DOS)')
        if E_limit is not None:
            ax.set_xlim(E_limit)
        else:
            ax.set_xlim([min(E), max(E)])

        ax.legend()
        # if tag_bool and not DOS_total:
        #     ax.set_title('Right-click to hide all\nMiddle-click to show all',
        #             loc='right')  # , size='medium')
        #     ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1),
        #             ncol=3)
        #     leg = interactive_legend()
        #     fig.subplots_adjust(right=0.55)
        
        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.
        if show:
            plt.show()
        if savefile is not None:
            plt.savefig(savefile, bbox_inches='tight')
        if fig is None:
            return ax
        else:
            return fig, ax
