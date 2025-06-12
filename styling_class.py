__author__ = "Carlos Ardila Gutierrez"
__maintainer__ = "Carlos Ardila Gutierrez"
__email__ = "carlos2248383@correo.uis.edu.co"
__date__ = "January 06, 2025"

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

class styler:
    '''
    Class to control style details of the graph.    
    '''
    def __init__(self, fig=None, ax=None):
        self.fig = fig
        self.ax = ax
        if ax is not None:
            self.ax.format_coord = lambda x, y: 'x={:g}, y={:g}'.format(x, y) # To show coords.

    def set_style(self, style='default'):
        '''
        Set the style of the graph. It uses the matplotlib styles.
        '''
        plt.style.use(style)

    def set_size(self, figSize=[10, 8], figdpi=300):
        '''
        Set the dimensions of the graph.
        '''
        # ------------- fig -------------
        figWidht, figHeight = figSize
        self.fig.set_size_inches(figWidht, figHeight)
        self.fig.set_dpi(figdpi)
        self.fig.subplots_adjust()
        # -------------------------------

    def set_title(self, maintitle=None, ax_title=None,
                    mtitle_fontsize=20, ax_title_fontsize=12):
        '''
        Set the titles of the graph.
        '''
        if maintitle is not None:
            self.fig.suptitle(maintitle, fontsize=mtitle_fontsize)
        if ax_title is not None:
            self.ax.set_title(ax_title, fontsize=ax_title_fontsize)

    def set_scale(self, scale_x='linear', scale_y='linear'):
        '''
        Set the scale used on each of the axes.
        '''
        self.ax.set_yscale(scale_x)
        self.ax.set_yscale(scale_y)

    def set_labels(self, xlabel=None, ylabel=None,
                    fontsize_xlabel=None, fontsize_ylabel=None):
        '''
        manage the labels of the graph.
        '''
        # ------------ Set the labels of the graph. ------------
        if xlabel is not None:
            self.ax.set_xlabel(xlabel)
        if ylabel is not None:
            self.ax.set_ylabel(ylabel)
        # ------------ Set the size of the text. ------------
        if fontsize_xlabel is not None:
            self.ax.xaxis.label.set_fontsize(fontsize_xlabel)
        if fontsize_ylabel is not None:
            self.ax.yaxis.label.set_fontsize(fontsize_ylabel)

    def set_ticks(self, xticks_step=None, yticks_step=None,
                  fontsize_xticks=None, fontsize_yticks=None,
                  erase_xticks=False, erase_yticks=False):
        # ------------ Control step of the ticks. ------------
        if xticks_step is not None:
            self.ax.xaxis.set_major_locator(ticker.MultipleLocator(xticks_step))
        if yticks_step is not None:
            self.ax.yaxis.set_major_locator(ticker.MultipleLocator(yticks_step))
        # ------------ set fontsize of the ticks. ------------
        if fontsize_xticks is not None:
            self.ax.tick_params(axis='x', labelsize=fontsize_xticks)
        if fontsize_yticks is not None:
            self.ax.tick_params(axis='y', labelsize=fontsize_yticks)
        # ------------------- erase ticks --------------------
        if erase_xticks:
            self.ax.set_xticks([])
        if erase_yticks:
            self.ax.set_yticks([])

    def set_lim(self, xlim=None, ylim=None):
        '''
        Set the limits of the graph.
        '''
        if xlim is not None:
            self.ax.set_xlim(xlim)
        if ylim is not None:
            self.ax.set_ylim(ylim)

    def set_frame_thickness(self, thickness=1):
        '''
        change te thickness of the frame of the graph.
        '''
        self.ax.spines['top'].set_linewidth(thickness)
        self.ax.spines['bottom'].set_linewidth(thickness)
        self.ax.spines['left'].set_linewidth(thickness)
        self.ax.spines['right'].set_linewidth(thickness)

    def set_legend(self, loc='best', fontsize=12):
        '''
        Set the legend of the graph.
        '''
        self.ax.legend(loc=loc, fontsize=fontsize)

    def return_fig_ax(self):
        return self.fig, self.ax

    def savefig(self, savefile, dpi=300, format='png'):
        '''
        Save the figure.
        '''
        self.fig.savefig(savefile+f'.{format}', dpi=dpi, format=format)

    def bandstructure_style(self, savefile=None):
        '''
        A default configuration for band structures figures.
        '''
        self.set_labels(ylabel=r'$E-E_{F} [eV]$', xlabel='', fontsize_ylabel=18)
        # self.set_scale()

        # self.ax.axhline(0, color='gray', linewidth=0.2)
        
        ylim = self.ax.get_ylim()
        self.set_ticks(yticks_step=min(abs(np.array(ylim)))//3, fontsize_yticks=14,
                        fontsize_xticks=14)

        self.set_size(figSize=[12, 8], figdpi=300)

        if savefile is not None:
            self.savefig(savefile)
        
        return self.fig, self.ax

    def style_2figs_h(self, fig, axs):
        self.fig = fig
        self.set_size(figSize=[12, 4], figdpi=300)
        # ------ left ------
        self.ax = axs[0]
        self.set_ticks(fontsize_xticks=12, fontsize_yticks=12)
        self.set_labels(fontsize_ylabel=12)
        fig, axs[0] = self.return_fig_ax()
        # ------ right ------
        self.ax = axs[1]
        self.set_ticks(fontsize_xticks=12, erase_yticks=True)
        self.set_labels(ylabel='', fontsize_ylabel=14)
        fig, axs[1] = self.return_fig_ax()
        # ----- space betweem subplots -----
        fig.subplots_adjust(wspace=0.045)

    def style_4figs_square(self, fig, axs):
        self.fig = fig
        self.set_size(figSize=[12, 8], figdpi=300)

        for i in range(2):
            # ------ left ------
            self.ax = axs[i][0]
            self.set_ticks(fontsize_xticks=12, fontsize_yticks=12)
            self.set_labels(fontsize_ylabel=12)
            fig, axs[i][0] = self.return_fig_ax()
            # ------ right ------
            self.ax = axs[i][1]
            self.set_ticks(fontsize_xticks=12, erase_yticks=True)
            self.set_labels(ylabel='', fontsize_ylabel=14)
            fig, axs[i][1] = self.return_fig_ax()
        for i in range(2):
            # ------- Top ----=--
            self.ax = axs[0][i]
            self.set_ticks(erase_xticks=True)
            fig, axs[0][i] = self.return_fig_ax()

        # ----- space betweem subplots -----
        fig.subplots_adjust(wspace=0.045, hspace=0.08)

    def __interactive_legend(ax=None):
        if ax is None:
            ax = plt.gca()
        if ax.legend_ is None:
            ax.legend()
        return InteractiveLegend(ax.get_legend())


class InteractiveLegend(object):
    '''
    Class taken from: tryoing to get the url...
    '''
    def __init__(self, legend):
        self.legend = legend
        self.fig = legend.axes.figure

        self.lookup_artist, self.lookup_handle = self._build_lookups(legend)
        self._setup_connections()

        self.update()

    def _setup_connections(self):
        for artist in self.legend.texts + self.legend.legendHandles:
            artist.set_picker(10) # 10 points tolerance

        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def _build_lookups(self, legend):
        labels = [t.get_text() for t in legend.texts]
        handles = legend.legendHandles
        label2handle = dict(zip(labels, handles))
        handle2text = dict(zip(handles, legend.texts))

        lookup_artist = {}
        lookup_handle = {}
        for artist in legend.axes.get_children():
            if artist.get_label() in labels:
                handle = label2handle[artist.get_label()]
                lookup_handle[artist] = handle
                lookup_artist[handle] = artist
                lookup_artist[handle2text[handle]] = artist

        lookup_handle.update(zip(handles, handles))
        lookup_handle.update(zip(legend.texts, handles))

        return lookup_artist, lookup_handle

    def on_pick(self, event):
        handle = event.artist
        if handle in self.lookup_artist:

            artist = self.lookup_artist[handle]
            artist.set_visible(not artist.get_visible())
            self.update()

    def on_click(self, event):
        if event.button == 3:
            visible = False
        elif event.button == 2:
            visible = True
        else:
            return

        for artist in self.lookup_artist.values():
            artist.set_visible(visible)
        self.update()

    def update(self):
        for artist in self.lookup_artist.values():
            handle = self.lookup_handle[artist]
            if artist.get_visible():
                handle.set_visible(True)
            else:
                handle.set_visible(False)
        self.fig.canvas.draw()

    def show(self):
        plt.show()

'''
def legend_thickness(self, thickness=1):
    plt.rc('legend',fontsize=20) # using a size in points
    plt.rc('legend',fontsize='medium')
''' 