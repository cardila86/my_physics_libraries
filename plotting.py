__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "August 05, 2025"

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os

class bands:
    def __init__(self,
        main_linewidth=1.3,
        main_linestyle='-',
        E_zero_color='gray',
        E_zero_linewidth=0.2,
        E_zero_linestyle='--',
        k_color='gray',
        k_linewidth=0.2,
        k_linestyle='-',
        background_linecolor='lightgray',
        marker_size=1):

        self.marker_size=marker_size

        self.main_linewidth=main_linewidth
        self.main_linestyle=main_linestyle

        self.background_linecolor=background_linecolor

        self.E_zero_color=E_zero_color
        self.E_zero_linewidth=E_zero_linewidth
        self.E_zero_linestyle=E_zero_linestyle

        self.k_color=k_color
        self.k_linewidth=k_linewidth
        self.k_linestyle=k_linestyle

    def plot_bands_scatter(self,
        kpoints,
        bands,
        parameter,
        E_zero,
        klabels,
        kticks,
        kbreaks,
        label,
        color,
        cmap,
        nbands,
        spin,
        orbitals,
        atoms,
        vmin,
        vmax,
        cbar_bool,
        ax):
        if type(kpoints)!=np.ndarray:
            kpoints = np.array(kpoints)
        if type(bands[0])!=np.ndarray:
            bands = [np.array(band) for band in bands]
        # ------ shift e fermi ------
        bands = [band-E_zero for band in bands]
        # ------------ ax, fig objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        # ------ plot klabels and kticks ------        
        for ktick in kticks:
            ax.axvline(ktick, color=self.k_color, linewidth=self.k_linewidth, linestyle=self.k_linestyle)
        ax.set_xticks(kticks)
        ax.set_xticklabels(klabels)

        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # ------------ discontinuities -----------
        if kbreaks is not None:
            for i in kbreaks:
                for j in range(len(kpoints)):
                    if j!=0 and abs(kpoints[j]-kticks[i])<abs(kpoints[j]-kpoints[j-1])/2: # Mira que este cerca a un ktick
                        for k in range(len(bands)):
                            bands[k][j] = np.nan
        # --- isolate nbands ---
        if nbands is None:
                pass
        elif type(nbands) is int or type(nbands) is float:
            bands=bands[0:int(nbands)]
            parameter=parameter[0:int(nbands)]
        elif type(nbands) is list:
            bands_new=[]
            parameter_new=[]
            for i in nbands:
                bands_new.append(bands[int(i)])
                parameter_new.append(parameter[int(i)])
        else:
            print('ERROR: nbands must be an integer, a float or a list of integers.\n'+
                'a '+str(type(nbands))+' type was recieved. Please check inputs.') 
            exit()
        # ======= distinguis between spin polarization and projected bands =======
        # --- spin polarized ---
        if spin is not None and orbitals is None and atoms is None:
            # -------- vmin and vmax --------
            if vmin is None:
                vmin = np.min(parameter)
            if vmax is None:
                vmax = np.max(parameter)

            norm = plt.Normalize(vmin, vmax)
            # ------ plotting --------
            # --- plot plain bands as background ---
            for band in bands:
                ax.plot(kpoints, band, c=self.background_linecolor, linewidth=self.main_linewidth/5, linestyle=self.main_linestyle)
            # --- scatter plot ---
            for i in range(len(bands)):
                ax.scatter(kpoints, bands[i], c=parameter[i], cmap=cmap, s=abs(parameter[i])*self.marker_size, norm=norm)
            if cbar_bool:
                    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label=None)
                    cbar.set_ticks([vmin,(vmax-vmin)/2+vmin,vmax])
        # --- projected bands ---
        elif orbitals is not None or atoms is not None:
            # organize parameters
            if orbitals is not None and atoms is None:
                pass
            elif orbitals is None and atoms is not None:
                parameter = parameter[:, :, :, 0]
            # colors and labels
            if type(color) != list:
                color_custom = False
            else:
                color_custom = True
                color = [np.array(i)/255 for i in color]
            color_list = [
            np.array([255, 0, 0])/255,
            np.array([0, 255, 0])/255,
            np.array([0, 0, 255])/255,
            np.array([255, 255, 0])/255,
            np.array([0, 255, 255])/255,
            np.array([255, 0, 255])/255,
            np.array([255, 255, 255])/255,
            np.array([0, 0, 0])/255,
            np.array([128, 128, 128])/255,
            np.array([139, 0, 0])/255,
            np.array([34, 139, 34])/255,
            np.array([0, 0, 128])/255,
            np.array([255, 165, 0])/255,
            np.array([128, 0, 128])/255,
            np.array([165, 42, 42])/255,
            np.array([0, 128, 128])/255,
            ]
            if color_custom:
                if len(orbitals)==1:
                    if type(orbitals[0])==list:
                        color_list.append(color[0])
                    else:
                        color_list[orbitals[0]] = color[0]
                else:
                    for i in range(len(orbitals)):
                        if type(orbitals[i])==list:
                            color_list.append(color[i])
                        else:
                            color_list[orbitals[i]] = color[i]

            color = color_list


            orbital_labels = [
            r"$s$",
            r"$p_{y}$",
            r"$p_{z}$",
            r"$p_{x}$",
            r"$d_{xy}$",
            r"$d_{yz}$",
            r"$d_{z^{2}}$",
            r"$d_{xz}$",
            r"$d_{x^{2}-y^{2}}$",
            r"$f_{y^{3}x^{2}}$",
            r"$f_{xyz}$",
            r"$f_{yz^{2}}$",
            r"$f_{z^{3}}$",
            r"$f_{xz^{2}}$",
            r"$f_{zx^{3}}$",
            r"$f_{x^{3}}$",
            ]
            # ----- check if there is a whole orbital -----
            for i in range(len(orbitals)):
                if type(orbitals[i])==list:
                    # -- orbitals --
                    if 1 in orbitals[i] and 2 in orbitals[i] and 3 in orbitals[i]:
                        orbitals[i] = len(orbital_labels)
                        orbital_labels.append('p')
                        if not color_custom:
                            color.append('b')
                    elif 4 in orbitals[i] and 5 in orbitals[i] and 6 in orbitals[i] and 7 in orbitals[i] and 8 in orbitals[i]:
                        orbitals[i] = len(orbital_labels)
                        orbital_labels.append('d')
                        if not color_custom:
                            color.append('y')
                    else:
                        orbitals[i] = len(orbital_labels)
                        orbital_labels.append('-')
                        if not color_custom:
                            color.append('g')
                else:
                    pass
            # ---------------------------------------------
            shape = parameter.shape
            # --- plot plain bands as background ---
            for band in bands:
                ax.plot(kpoints, band, c=self.background_linecolor, linewidth=self.main_linewidth/5, linestyle=self.main_linestyle)
            # -- plot scatter --
            kpoints = list(kpoints)*len(bands)
            for j in range(shape[-1]):
                s = abs(parameter[:, :, j])*self.marker_size
                ax.scatter(kpoints, bands, color=color[orbitals[j]], label=orbital_labels[orbitals[j]],s=s)
        # ------------- set limits --------------
        bool_klabels = [i=='' for i in klabels]
        klabels_filtered = [i for i in klabels if i!='']
        if any(bool_klabels):
            if len(klabels_filtered)!=2:
                ind = klabels.index('')
                ax.set_xlim([0,kticks[ind-1]])
            elif len(klabels_filtered)==2:
                ind1 = klabels.index(klabels_filtered[0])
                ind2 = klabels.index(klabels_filtered[1])
                ax.set_xlim([kticks[ind1],kticks[ind2]])
        else:
            ax.set_xlim([kticks[0], kticks[-1]])

        return fig, ax

    def plot_bands_parametric(self,
        kpoints,
        bands,
        parameter,
        E_zero,
        klabels,
        kticks,
        kbreaks,
        label,
        color,
        cmap,
        nbands,
        spin,
        orbitals,
        atoms,
        vmin,
        vmax,
        cbar_bool,
        ax):
        assert bands.shape == parameter.shape, "E and parameter must have the same shape."
        if type(kpoints)!=np.ndarray:
            kpoints = np.array(kpoints)
        if type(bands[0])!=np.ndarray:
            bands = [np.array(band) for band in bands]
        # ------ shift e fermi ------
        bands = [band-E_zero for band in bands]
        # ------------ ax, fig objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        # ------ plot klabels and kticks ------        
        for ktick in kticks:
            ax.axvline(ktick, color=self.k_color, linewidth=self.k_linewidth, linestyle=self.k_linestyle)
        ax.set_xticks(kticks)
        ax.set_xticklabels(klabels)

        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # ------------ discontinuities -----------
        if kbreaks is not None:
            for i in kbreaks:
                for j in range(len(kpoints)):
                    if j!=0 and abs(kpoints[j]-kticks[i])<abs(kpoints[j]-kpoints[j-1])/2: # Mira que este cerca a un ktick
                        for k in range(len(bands)):
                            bands[k][j] = np.nan
        # -------- vmin and vmax --------
        if vmin is None:
            vmin = np.min(parameter)
        if vmax is None:
            vmax = np.max(parameter)
        
        norm = plt.Normalize(vmin, vmax)
        # ------ plotting --------
        if nbands is None:
            pass
        elif type(nbands) is int or type(nbands) is float:
            bands=bands[0:int(nbands)]
            parameter=parameter[0:int(nbands)]
        elif type(nbands) is list:
            bands_new=[]
            parameter_new=[]
            for i in nbands:
                bands_new.append(bands[int(i)])
                parameter_new.append(parameter[int(i)])
        else:
            print('ERROR: nbands must be an integer, a float or a list of integers.\n'+
                  'a '+str(type(nbands))+' type was recieved. Please check inputs.') 
            exit()
        
        for i in range(len(bands)):
            points = np.array([kpoints, bands[i]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=cmap, norm=norm)

            lc.set_array(parameter[i])
            # lc.set_linewidth(bandwidth)

            line = ax.add_collection(lc)
            if i==0 and cbar_bool:
                cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                            ax=ax, orientation='vertical', label=None)
                # cbar.set_ticks([])

        # ------------- set limits --------------
        bool_klabels = [i=='' for i in klabels]
        klabels_filtered = [i for i in klabels if i!='']
        if any(bool_klabels):
            if len(klabels_filtered)!=2:
                ind = klabels.index('')
                ax.set_xlim([0,kticks[ind-1]])
            elif len(klabels_filtered)==2:
                ind1 = klabels.index(klabels_filtered[0])
                ind2 = klabels.index(klabels_filtered[1])
                ax.set_xlim([kticks[ind1],kticks[ind2]])
        else:
            ax.set_xlim([kticks[0], kticks[-1]])

        return fig, ax

    def plot_bands_processed(self,
        kpoints,
        E,
        E_zero,
        klabels,
        kticks,
        kbreaks,
        label,
        color,
        nbands,
        ax):
        # ------------ ax, fig objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        # ------ plot klabels and kticks ------        
        for ktick in kticks:
            ax.axvline(ktick, color=self.k_color, linewidth=self.k_linewidth, linestyle=self.k_linestyle)
        ax.set_xticks(kticks)
        ax.set_xticklabels(klabels)

        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # ------------ discontinuities -----------
        if kbreaks is not None:
            for i in kbreaks:
                for j in range(len(kpoints)):
                    if j!=0 and abs(kpoints[j]-kticks[i])<abs(kpoints[j]-kpoints[j-1])/2: # Mira que este cerca a un ktick
                        for k in range(len(E)):
                            E[k][j] = np.nan
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
        bool_klabels = [i=='' for i in klabels]
        klabels_filtered = [i for i in klabels if i!='']
        if any(bool_klabels):
            if len(klabels_filtered)!=2:
                ind = klabels.index('')
                ax.set_xlim([0,kticks[ind-1]])
            elif len(klabels_filtered)==2:
                ind1 = klabels.index(klabels_filtered[0])
                ind2 = klabels.index(klabels_filtered[1])
                ax.set_xlim([kticks[ind1],kticks[ind2]])
        else:
            ax.set_xlim([kticks[0], kticks[-1]])
        # ------ output -------
        if label is not None:
            ax.legend()
        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.

        return fig, ax

    def plot_pbands_p4vasp(self,
        kpoints,
        bands,
        bands_projections,
        E_zero,
        E_limit,
        klabels,
        kticks,
        kbreaks,
        colors,
        label,
        nbands,
        ax=None,
        show=False,
        savefile=None):
        # ------------ ax, fig objects -----------
        if ax is None and kbreaks is None:
            fig, ax = plt.subplots()
        elif ax is None and kbreaks is not None:
            num_ax = len(kbreaks)+1
            fig, ax = plt.subplots(num_ax)
        else:
            fig = ax.get_figure()
        # ------- Order Projection and asign colors -------
        # --- define colors ---
        for i in range(len(colors)):
            for j in range(len(colors[i])):
                colors[i][j] = colors[i][j]/255
            colors[i] = np.array(colors[i])
        # ---- checks that there are 3 list of projections ----
        print(f'There are {len(bands)} bands. The number of projections is: ', len(bands_projections))
        while len(bands_projections)<3:
            bands_projections.append(np.zeros(np.shape(bands_projections[0])))
        # la proyeccion es de la forma 3 proyecciones x nbandas x nkpoints. Para pasar
        # esto a color, la lista pasa a ser nbandas x 3 colores x nkpoints
        colors_projections = np.empty((len(bands), 3, len(kpoints))) 
        # i es la banda, j es la componente del color, k es el punto k...
        l = []
        for i in range(len(bands)):
            for j in range(3):
                for k in range(len(kpoints)):
                    # promedio ponderado
                    colors_projections[i][j][k] = 1 - (bands_projections[0][i][k]*(1-colors[0][j]) + bands_projections[1][i][k]*(1-colors[1][j]) + bands_projections[2][i][k]*(1-colors[2][j])) #/(colors[0][j]+colors[1][j]+colors[2][j])
                    # Checks for saturated values
                    if colors_projections[i][j][k]>1:
                        colors_projections[i][j][k] = 1
                    elif colors_projections[i][j][k]<0:
                        colors_projections[i][j][k] = 0                    
        # ----------- plot projected bands ------------
        for band, colors in zip(bands, colors_projections):
            # ------- asign colors -------
            band-=E_zero
            color_projection = np.column_stack([colors[0], colors[1], colors[2]])
            # ---- create segments ----
            points = np.array([kpoints, band]).transpose().reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, colors=color_projection, linewidths=self.main_linewidth)
 
            ax.add_collection(lc)
        # ------ plot klabels and kticks ------
        for ktick in kticks:
            ax.axvline(ktick, color=self.k_color, linewidth=self.k_linewidth, linestyle=self.k_linestyle)
        ax.set_xticks(kticks)
        ax.set_xticklabels(klabels)

        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # ---------- set limits and ylabel ----------
        bool_klabels = [i=='' for i in klabels]
        klabels_filtered = [i for i in klabels if i!='']
        if any(bool_klabels):
            if len(klabels_filtered)!=2:
                ind = klabels.index('')
                ax.set_xlim([0,kticks[ind-1]])
            elif len(klabels_filtered)==2:
                ind1 = klabels.index(klabels_filtered[0])
                ind2 = klabels.index(klabels_filtered[1])
                ax.set_xlim([kticks[ind1],kticks[ind2]])
        else:
            ax.set_xlim([kticks[0], kticks[-1]])
        ax.set_ylim(E_limit)
        ax.set_ylabel(r'$E-E_{F} [eV]$')
        # ------ output -------
        if label is not None:
            ax.legend()
        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.        
        return fig, ax

    def plot_pyprocar(self,
        path_read,
        atoms=None,
        orbitals=None,
        spins=None,
        mode='plain',
        E_limit=[-1, 1],
        E_zero=0,
        klabels= None,
        kticks=None,
        cmap='jet',
        cbar=False,
        ax=None,
        show=False,
        savefile=None):
        '''
        Uses pyprocar 6.1.3
        '''
        import pyprocar
        
        if not os.path.exists(path_read+'/PROCAR-repaired'):
            pyprocar.repair(path_read+'/PROCAR', path_read+'/PROCAR-repaired')
        
        bands = pyprocar.bandsplot(
                code='vasp', dirname=path_read, mode=mode,
                fermi=E_zero, fermi_color = 'black', fermi_linestyle='dashed', fermi_linewidth=0.5,
                cmap=cmap, linestyle=['solid', 'solid'], linewidth=[1, 1],
                spins=spins, atoms=atoms, orbitals=orbitals, 
                elimit=E_limit, knames=klabels, kticks=kticks,                
                show=show, savefig=savefile, dpi=300, ax=ax)
        
        fig, ax = bands.fig, bands.ax

        return fig, ax    

    def plot_surface_scatter(self,
        kpoints,
        bands,
        parameter,
        E_zero,
        color,
        nbands,
        spin,
        orbitals,
        atoms,
        vmin,
        vmax,
        cbar_bool,
        ax):
        if type(kpoints)!=np.ndarray:
            kpoints = np.array(kpoints)
        if type(bands[0])!=np.ndarray:
            bands = [np.array(band) for band in bands]
        # ------ shift e fermi ------
        bands = [band-E_zero for band in bands]
        # ------------ ax, fig objects -----------
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
        else:
            fig = ax.get_figure()
       # --- isolate nbands ---
        if nbands is None:
                pass
        elif type(nbands) is int or type(nbands) is float:
            bands=bands[0:int(nbands)]
            if parameter is not None:
                parameter=parameter[0:int(nbands)]
        elif type(nbands) is list:
            bands_new=[]
            parameter_new=[]
            for i in nbands:
                bands_new.append(bands[int(i)])
                if parameter is not None:
                    parameter_new.append(parameter[int(i)])
            
            if parameter is not None:
                parameter = parameter_new
            bands = bands_new
        else:
            print('ERROR: nbands must be an integer, a float or a list of integers.\n'+
                'a '+str(type(nbands))+' type was recieved. Please check inputs.') 
            exit()

        # ======= distinguis between spin polarization and projected bands =======
        # --- spin polarized ---
        if spin is not None and orbitals is None and atoms is None:
            # -------- vmin and vmax --------
            if vmin is None:
                vmin = np.min(parameter)
            if vmax is None:
                vmax = np.max(parameter)

            norm = plt.Normalize(vmin, vmax)
            # ------ plotting --------
            custom_cmap = LinearSegmentedColormap.from_list("my_cmap",[(0, 'blue'), (0.5, 'lightgray'), (1, 'red')],N=256)
            # --- scatter plot ---
            for i in range(len(bands)):
                ax.scatter(kpoints[0], kpoints[1], bands[i], c=parameter[i], cmap=custom_cmap, s=abs(parameter[i])*self.marker_size, norm=norm)
            if cbar_bool:
                    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=custom_cmap), ax=ax, orientation='vertical', label=None)
                    cbar.set_ticks([vmin,(vmax-vmin)/2+vmin,vmax])
        # --- projected bands ---
        elif orbitals is not None or atoms is not None:
            # organize parameters
            if orbitals is not None and atoms is None:
                pass
            elif orbitals is None and atoms is not None:
                parameter = parameter[:, :, :, 0]
            # colors and labels
            if type(color) != list:
                color_custom = False
            else:
                color_custom = True
                color = [np.array(i)/255 for i in color]
            color_list = [
            np.array([255, 0, 0])/255,
            np.array([0, 255, 0])/255,
            np.array([0, 0, 255])/255,
            np.array([255, 255, 0])/255,
            np.array([0, 255, 255])/255,
            np.array([255, 0, 255])/255,
            np.array([255, 255, 255])/255,
            np.array([0, 0, 0])/255,
            np.array([128, 128, 128])/255,
            np.array([139, 0, 0])/255,
            np.array([34, 139, 34])/255,
            np.array([0, 0, 128])/255,
            np.array([255, 165, 0])/255,
            np.array([128, 0, 128])/255,
            np.array([165, 42, 42])/255,
            np.array([0, 128, 128])/255,
            ]
            if color_custom:
                if len(orbitals)==1:
                    if type(orbitals[0])==list:
                        color_list.append(color[0])
                    else:
                        color_list[orbitals[0]] = color[0]
                else:
                    for i in range(len(orbitals)):
                        if type(orbitals[i])==list:
                            color_list.append(color[i])
                        else:
                            color_list[orbitals[i]] = color[i]

            color = color_list


            orbital_labels = [
            r"$s$",
            r"$p_{y}$",
            r"$p_{z}$",
            r"$p_{x}$",
            r"$d_{xy}$",
            r"$d_{yz}$",
            r"$d_{z^{2}}$",
            r"$d_{xz}$",
            r"$d_{x^{2}-y^{2}}$",
            r"$f_{y^{3}x^{2}}$",
            r"$f_{xyz}$",
            r"$f_{yz^{2}}$",
            r"$f_{z^{3}}$",
            r"$f_{xz^{2}}$",
            r"$f_{zx^{3}}$",
            r"$f_{x^{3}}$",
            ]
            # ----- check if there is a whole orbital -----
            for i in range(len(orbitals)):
                if type(orbitals[i])==list:
                    # -- orbitals --
                    if 1 in orbitals[i] and 2 in orbitals[i] and 3 in orbitals[i]:
                        orbitals[i] = len(orbital_labels)
                        orbital_labels.append('p')
                        if not color_custom:
                            color.append('b')
                    elif 4 in orbitals[i] and 5 in orbitals[i] and 6 in orbitals[i] and 7 in orbitals[i] and 8 in orbitals[i]:
                        orbitals[i] = len(orbital_labels)
                        orbital_labels.append('d')
                        if not color_custom:
                            color.append('y')
                    else:
                        orbitals[i] = len(orbital_labels)
                        orbital_labels.append('-')
                        if not color_custom:
                            color.append('g')
                else:
                    pass
            # ---------------------------------------------
            shape = parameter.shape
            # --- plot plain bands as background ---
            for band in bands:
                ax.scatter(kpoints[0], kpoints[1], band, c=self.background_linecolor)
            # -- plot scatter --
            kpoints = list(kpoints)*len(bands)
            for j in range(shape[-1]):
                s = abs(parameter[:, :, j])*self.marker_size
                ax.scatter(kpoints[0], kpoints[1], bands, color=color[orbitals[j]], label=orbital_labels[orbitals[j]],s=s)
        # ---
        elif spin is None and orbitals is None and atoms is None:
            for band in bands:
                ax.scatter(kpoints[0], kpoints[1], band, c=color, s=self.marker_size)
        return fig, ax



    def __discontinuities(self):
        pass

class dos:
    def __init__(self,
        main_linewidth=1.3,
        main_linestyle='-',
        E_zero_color='gray',
        E_zero_linewidth=0.2,
        E_zero_linestyle='--',
        k_color='gray',
        k_linewidth=0.2,
        k_linestyle='-'):

        self.main_linewidth=main_linewidth
        self.main_linestyle=main_linestyle

        self.E_zero_color=E_zero_color
        self.E_zero_linewidth=E_zero_linewidth
        self.E_zero_linestyle=E_zero_linestyle

        self.k_color=k_color
        self.k_linewidth=k_linewidth
        self.k_linestyle=k_linestyle

    # ------------ ax, fig objects -----------
    def plot_ldos_vaspkit(self,
        code,
        E,
        orbitals,
        total,
        orbitals_labels,        
        orbitals_tag,
        colors,
        E_limit,
        DOS_limit,
        E_zero,
        fig_ortientation,
        ax):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        if colors is not None:
            ax.set_prop_cycle('color', colors)  # Colores a usar
        # ------- figure orientation -------
        if fig_ortientation=='horizontal' or fig_ortientation=='Horizontal' or fig_ortientation=='h' or fig_ortientation=='H':
            axes_vertical = False
        elif fig_ortientation=='vertical' or fig_ortientation=='Vertical' or fig_ortientation=='v' or fig_ortientation=='V':
            axes_vertical = True
        # ----- Sets zero energy level -----
        E -= E_zero
        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # --------- plot partial DOS ---------
        if code=='vaspkit':
            if orbitals_tag=='all':
                for orbs, label in zip(orbitals, orbitals_labels):
                    if axes_vertical==False:
                        ax.plot(E, orbs, label=label, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                    elif axes_vertical==True:
                        ax.plot(orbs, E, label=label, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
            else:
                for tag in orbitals_tag:
                    index = orbitals_labels.index(tag)
                    if axes_vertical==False:
                        ax.plot(E, orbitals[index], label=tag, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                    elif axes_vertical==True:
                        ax.plot(orbitals[index], E, label=tag, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
        elif code=='p4vasp':
            for orb, tag_orb in zip(orbitals, orbitals_tag):
                if axes_vertical==False:
                    ax.plot(E, orb, label=tag_orb, linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                elif axes_vertical==True:
                    ax.plot(orb, E, label=tag_orb, linewidth=self.main_linewidth, linestyle=self.main_linestyle)            
        # --------- plot total DOS ---------
        if axes_vertical==False:
            ax.plot(E, total, color='gray', linewidth=self.main_linewidth*0.8, linestyle='--')
            ax.set_xlabel('Energy [eV]')
            ax.set_ylabel('DOS [states/eV]')
            if E_limit is not None:
                ax.set_xlim(E_limit)
            else:
                ax.set_xlim([min(E), max(E)])
            if DOS_limit is not None:
                ax.set_ylim(DOS_limit)
        elif axes_vertical==True:
            ax.plot(total, E, color='gray', linewidth=self.main_linewidth*0.8, linestyle='--')
            ax.set_ylabel('Energy [eV]')
            ax.set_xlabel('DOS [states/eV]')
            if E_limit is not None:
                ax.set_ylim(E_limit)
            else:
                ax.set_ylim([min(E), max(E)])
            if DOS_limit is not None:
                ax.set_xlim(DOS_limit)

        ax.legend()        
        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.

        return fig, ax

class nodes:
    def __init__(self,
        main_linewidth=1.3,
        main_linestyle='-',
        E_zero_color='gray',
        E_zero_linewidth=0.2,
        E_zero_linestyle='--',
        k_color='gray',
        k_linewidth=0.2,
        k_linestyle='-'):

        self.main_linewidth=main_linewidth
        self.main_linestyle=main_linestyle

        self.E_zero_color=E_zero_color
        self.E_zero_linewidth=E_zero_linewidth
        self.E_zero_linestyle=E_zero_linestyle

        self.k_color=k_color
        self.k_linewidth=k_linewidth
        self.k_linestyle=k_linestyle


    def plot_nodes_wanniertools(self,
        path_read,
        E_limit=[-1, 1],
        E_zero=0,
        path_k=None,
        klabels= None,
        kticks=None,
        kbreaks=None,
        label=None,
        color='k',
        nbands=None,
        ax=None,
        show=False,
        savefile=None):
        '''
        Plot nodes calculated by wannierTools and if wanted, their chirality.
        '''
        # ------------ ax, fig objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        # -------------- reads info --------------
        kvec_car, gap, E, kvec_dir = self._read_nodes_wannierTools(path_read)

class ahc:
    def __init__(self,
        main_linewidth=1.3,
        main_linestyle='-',
        E_zero_color='gray',
        E_zero_linewidth=0.2,
        E_zero_linestyle='--',
        k_color='gray',
        k_linewidth=0.2,
        k_linestyle='-'):

        self.main_linewidth=main_linewidth
        self.main_linestyle=main_linestyle

        self.E_zero_color=E_zero_color
        self.E_zero_linewidth=E_zero_linewidth
        self.E_zero_linestyle=E_zero_linestyle

        self.k_color=k_color
        self.k_linewidth=k_linewidth
        self.k_linestyle=k_linestyle


    def plot_AHC_wannierberri(self,
        ahc_data,
        ahc_axis,
        E_limit,
        ahc_limit,
        E_zero,
        norm_val,
        colors,
        spin,
        iteration,
        fig_orientation,
        scale,
        ax):
        # ----------- fig, ax objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        if colors is not None:
            ax.set_prop_cycle('color', colors)  # Colores a usar
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
                    if x==0:
                        E-=E_zero
                    AHC_down=a[:,x+1]
                    AHC_up  =a[:,x+1+3]
                    # ----- curv -----
                    if (fig_orientation=='vertical'or'v'or'vertical_rigth'or'vr'):
                        if spin=='both':
                            ax.plot(AHC_down/norm_val, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                            n+=1
                            ax.plot(AHC_up/norm_val, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='up':
                            ax.plot(AHC_up/norm_val, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='down':
                            ax.plot(AHC_down/norm_val, E, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)  
                    elif (fig_orientation=='horizontal'or'h'):
                        if spin=='both':
                            ax.plot(E, AHC_down/norm_val, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                            n+=1
                            ax.plot(E, AHC_up/norm_val, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='up':
                            ax.plot(E, AHC_up/norm_val, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
                        elif spin=='down':
                            ax.plot(E, AHC_down/norm_val, label=curv_label[n], linewidth=self.main_linewidth, linestyle=self.main_linestyle)  
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

        return fig, ax    

class curvature:
    def __init__(self,
        main_linewidth=1.3,
        main_linestyle='-',
        E_zero_color='gray',
        E_zero_linewidth=0.2,
        E_zero_linestyle='--',
        k_color='gray',
        k_linewidth=0.2,
        k_linestyle='-'):
        from .macrodata_refinement.utils import read
        self.r_c = read()

        self.main_linewidth=main_linewidth
        self.main_linestyle=main_linestyle

        self.E_zero_color=E_zero_color
        self.E_zero_linewidth=E_zero_linewidth
        self.E_zero_linestyle=E_zero_linestyle

        self.k_color=k_color
        self.k_linewidth=k_linewidth
        self.k_linestyle=k_linestyle


    def plot_curvature_path(self,
        path_read,
        axis,
        kticks,
        klabels,
        kbreaks,
        colors,
        curv_limit,
        ax):
        # ------------ ax, fig objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        # -------- axis name --------
        if int(axis) == 0:
            axisName = ''
        elif int(axis) == 1:
            axisName = '_x'
        elif int(axis) == 2:
            axisName = '_y'
        elif int(axis) == 3:
            axisName = '_z'
        # ---- E Zero ----
        ax.axhline(0, color=self.E_zero_color, linewidth=self.E_zero_linewidth, linestyle=self.E_zero_linestyle)
        # ---- kpoints ---
        if kticks is not None:
            for i in range(len(kticks)):
                if kbreaks is not None:
                    if i in kbreaks:
                        ax.axvline(kticks[i], color='k', linewidth=self.k_linewidth*1.8, linestyle=self.k_linestyle)    
                    else:
                        ax.axvline(kticks[i], color=self.k_color, linewidth=self.k_linewidth, linestyle=self.k_linestyle)    
        # ---- klabels ---
        if klabels is not None:
            ax.set_xticks(kticks,klabels)

        # ----- curv -----
        if axis==0:
            curvLabel = [r'$\Omega_x$', r'$\Omega_y$', r'$\Omega_z$']
            for i in [1, 2, 3]:
                # -------------- reads info --------------
                kpoints, curv = self.r_c._read_curvature_wannier90(path_read, i)
                # ----- plot -----
                ax.plot(kpoints, curv, color=colors[i-1], label=curvLabel[i-1], linewidth=self.main_linewidth, linestyle=self.main_linestyle)        
                # ----------------
        else:
            # -------------- reads info --------------
            kpoints, curv = self.r_c._read_curvature_wannier90(path_read, axis)    
            # ----- plot -----
            ax.plot(kpoints, curv, color=colors[0], linewidth=self.main_linewidth, linestyle=self.main_linestyle)
            # ----------------
        
        ax.legend(loc='best')

        # ------- Cut unnamed paths and axis limits-------
        boolKnames = [i=='' for i in klabels]
        knamesFiltered = [i for i in klabels if i!='']
        if any(boolKnames):
            if len(knamesFiltered)!=2:
                ind = klabels.index('')
                ax.set_xlim([0,kticks[ind-1]])
            elif len(knamesFiltered)==2:
                ind1 = klabels.index(knamesFiltered[0])
                ind2 = klabels.index(knamesFiltered[1])
                ax.set_xlim([kticks[ind1],kticks[ind2]])
        else:
            ax.set_xlim([kticks[0],kticks[-1]])

        if curv_limit is not None:
            ax.set_ylim(curv_limit)
        else:
            # ylimit = [min(curv)-0.025*(max(y)-min(y)),max(y)+0.025*(max(y)-min(y))]
            # ax.set_ylim(ylimit)
            pass
        # ----- set labels -----
        kBold = r'\mathbf{k}'
        ax.set_ylabel(fr'$-\Omega{axisName}({kBold})$  [ $\AA^2$ ]')

        ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.

        return fig, ax