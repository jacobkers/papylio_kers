# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:44:52 2018

@author: ivoseverins
"""
import wx
import wx.lib.mixins.inspection as wit
import PySide2
import matplotlib as mpl
# mpl.use('WXAgg')
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
# from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

import numpy as np
from pathlib2 import Path

from PySide2.QtWidgets import QMainWindow, QPushButton, QWidget, QVBoxLayout
from PySide2.QtGui import QKeySequence
from PySide2.QtCore import Qt

import sys
import time

import numpy as np

# from matplotlib.backends.qt_compat import QtWidgets
from PySide2 import QtWidgets
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure




class TracePlotWindow(QMainWindow):
    def __init__(self, dataset=None, plot_variables=['intensity', 'FRET'],
                 ylims=[(0, 35000), (0, 1)], colours=[('g', 'r'), ('b')], save_path=None):
        super().__init__()

        self.setWindowTitle("Traces")
        button = QPushButton("Press Me!")

        # Set the central widget of the Window.
        self.setCentralWidget(button)


        self.dataset = dataset
        self.plot_variables = plot_variables
        self.ylims = ylims
        self.colours = colours

        if save_path is None:
            self.save_path = save_path
        else:
            self.save_path = Path(save_path)

        self.canvas = TracePlotCanvas(self, width=5, height=4, dpi=100)

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(self.canvas)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

        self.molecule_index = 0

        self.show()


    @property
    def molecule_index(self):
        return self._molecule_index

    @molecule_index.setter
    def molecule_index(self, molecule_index):
        self._molecule_index = molecule_index
        self.molecule = self.dataset.isel(molecule=self._molecule_index)

    def next_molecule(self):
        if (self.molecule_index+1) < len(self.dataset.molecule):
            self.molecule_index += 1

    def previous_molecule(self):
        if self.molecule_index > 0:
            self.molecule_index -= 1

    def update_current_molecule(self):
        self.molecule_index = self.molecule_index

    @property
    def molecule(self):
        return self.canvas.molecule

    @molecule.setter
    def molecule(self, molecule):
        self.canvas.molecule = molecule

    def keyPressEvent(self, e):
        key = e.key()
        if key == Qt.Key_Right: # Right arrow
            self.next_molecule()
        elif key == Qt.Key_Left: # Left arrow
            self.previous_molecule()
        elif key == Qt.Key_Space: # Spacebar
            self.dataset.selected[dict(molecule=self.molecule_index)] = ~self.dataset.selected[dict(molecule=self.molecule_index)]
            self.update_current_molecule()
        elif key == Qt.Key_S: # S
            self.trace_panel.save()


class TracePlotCanvas(FigureCanvas):
    # Kader om plot als geselecteerd
    # Autosave function
    def __init__(self, parent=None, width=14, height=7, dpi=100):
        self.figure = mpl.figure.Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)  # , figsize=(2, 2))
        super().__init__(self.figure)
        self.parent = parent
        plot_variables = self.parent.plot_variables

        grid = self.figure.add_gridspec(len(plot_variables), 2, width_ratios=[10, 1]) #, height_ratios=(2, 7),
                         # left=0.1, right=0.9, bottom=0.1, top=0.9,
                         # wspace=0.05, hspace=0.05)

        self.plot_axes = {}
        self.histogram_axes = {}

        for i, plot_variable in enumerate(plot_variables):
            plot = self.figure.add_subplot(grid[i, 0])
            histogram = self.figure.add_subplot(grid[i, 1], sharey=plot)

            if i > 0:
                plot.sharex(self.plot_axes[plot_variables[0]])
                # histogram.sharex(self.histogram_axes[plot_variables[0]])

            plot.set_ylim(self.parent.ylims[i])
            plot.set_ylabel(plot_variable)

            if i == len(plot_variables)-1:
                plot.set_xlabel('Time')

            histogram.get_yaxis().set_visible(False)

            self.plot_axes[plot_variable] = plot
            self.histogram_axes[plot_variable] = histogram


        # self.intensity_plot = self.figure.add_subplot(grid[0, 0])
        # self.FRET_plot = self.figure.add_subplot(grid[1, 0], sharex=self.intensity_plot)
        # self.intensity_histogram = self.figure.add_subplot(grid[0, 1], sharey=self.intensity_plot)
        # self.FRET_histogram = self.figure.add_subplot(grid[1, 1], sharex=self.intensity_histogram, sharey=self.FRET_plot)

        # self.figure = plt.Figure(dpi=dpi, figsize=(2,2))
        #
        # self.axis = self.figure.gca()

        #self.figure, self.axes = mpl.figure.Figure().subplots(2,1)



        self._molecule = None

        self.plot_artists = {}
        self.histogram_artists = {}

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, molecule):
        self._molecule = molecule

        # g = molecule.intensity.sel(channel=0).values
        # r = molecule.intensity.sel(channel=1).values
        # e = molecule.FRET.values

        if not self.plot_artists:
            for i, plot_variable in enumerate(self.parent.plot_variables):
                self.plot_artists[plot_variable] = self.plot_axes[plot_variable].plot(molecule[plot_variable].T)
                if i==0:
                    self.title_artist = self.plot_axes[plot_variable].set_title('')
                for j, plot_artist in enumerate(self.plot_artists[plot_variable]):
                    plot_artist.set_color(self.parent.colours[i][j])
                #molecule.intensity.plot.line(x='frame', ax=self.plot_axes[plot_variable], color=self.parent.colours[i])
                self.histogram_artists[plot_variable] = self.histogram_axes[plot_variable].hist(molecule[plot_variable].T,
                                bins=50, orientation='horizontal', range=self.plot_axes[plot_variable].get_ylim(),
                                                                    color=self.parent.colours[i], alpha=0.5)[2]
                if not isinstance(self.histogram_artists[plot_variable], list):
                    self.histogram_artists[plot_variable] = [self.histogram_artists[plot_variable]]

            # self.artists += [self.intensity_plot.plot(g, c='g')]
            # self.artists += [self.intensity_plot.plot(r, c='r')]
            # self.artists += [self.FRET_plot.plot(e, c='b')]
            # self.artists += [[self.intensity_plot.set_title('test')]]
            # self.artists += [self.intensity_histogram.hist(g, bins=100, orientation='horizontal',
            #                                                range=self.intensity_plot.get_ylim(), color='g', alpha=0.5)[2]]
            # self.artists += [self.intensity_histogram.hist(r, bins=100, orientation='horizontal',
            #                                                range=self.intensity_plot.get_ylim(), color='r', alpha=0.5)[2]]
            # self.artists += [self.FRET_histogram.hist(e, bins=100, orientation='horizontal',
            #                                           range=self.FRET_plot.get_ylim(), color='b')[2]]

            #self.axes[1].plot(molecule.E(), animate=True)
            artists = [self.title_artist] + \
                      [a for b in self.plot_artists.values() for a in b] + \
                      [a for c in self.histogram_artists.values() for b in c for a in b]
            self.bm = BlitManager(self, artists)
            self.draw()

        # for axis in self.axes:
        #     axis.cla()

        for i, plot_variable in enumerate(self.parent.plot_variables):
            data = np.atleast_2d(molecule[plot_variable])
            if self._molecule.selected.item():
                selection_string = ' | Selected'
            else:
                selection_string = ''

            self.title_artist.set_text(f'File: {molecule.file.values} | Molecule: {molecule.molecule_in_file.values}' + selection_string)#| Sequence: {molecule.sequence_name.values}')
            for j in range(len(data)):
                self.plot_artists[plot_variable][j].set_ydata(data[j])
                n, _ = np.histogram(data[j], 50, range=self.plot_axes[plot_variable].get_ylim())
                for count, artist in zip(n, self.histogram_artists[plot_variable][j]):
                    artist.set_width(count)



        # self.artists[0][0].set_ydata(g)
        # self.artists[1][0].set_ydata(r)
        # self.artists[2][0].set_ydata(e)
        # self.artists[3][0].set_text(molecule.sequence_name.values)
        # n, _ = np.histogram(g, 100, range=self.intensity_plot.get_ylim())
        # for count, artist in zip(n, self.artists[4]):
        #     artist.set_width(count)
        # n, _ = np.histogram(r, 100, range=self.intensity_plot.get_ylim())
        # for count, artist in zip(n, self.artists[5]):
        #     artist.set_width(count)
        # n, _ = np.histogram(e, 100, range=self.FRET_plot.get_ylim())
        # for count, artist in zip(n, self.artists[6]):
        #     artist.set_width(count)
        #     #for count, rect in zip(n, bar_container.patches):
        # tell the blitting manager to do its thing
        self.bm.update()



        # self.axes[0].plot(molecule.intensity.T)
        # self.axes[1].plot(molecule.E())
        # self.canvas.draw()

        # empty_cell = (0, 0)
        #
        # donor_checkbox = wx.CheckBox(self, -1, label="Donor", name="Donor")
        # acceptor_checkbox = wx.CheckBox(self, -1, label="Acceptor", name="Acceptor")
        #
        # channel_sizer = wx.FlexGridSizer(2,2, gap=wx.Size(10,0))
        # channel_sizer.Add(wx.StaticText(self, label="Channels:"), 5, wx.EXPAND, 0)
        # channel_sizer.Add(donor_checkbox, 1, wx.EXPAND, 0)
        # channel_sizer.Add(empty_cell, 0, wx.EXPAND, 0)
        # channel_sizer.Add(acceptor_checkbox, 0, wx.EXPAND, 0)
        #
        # average_image_radio_button = wx.RadioButton(self, -1, label="Average image", name="Average image")
        # maximum_projection_radio_button = wx.RadioButton(self, -1, label="Maximum projection", name="Maximum projection")
        #
        # image_type_sizer = wx.FlexGridSizer(2,2, gap=wx.Size(10,0))
        # image_type_sizer.AddGrowableCol(0,1)
        # image_type_sizer.AddGrowableCol(1,3)
        # image_type_sizer.Add(wx.StaticText(self, label="Image type atestat:"), 0, wx.EXPAND, 0)
        # image_type_sizer.Add(average_image_radio_button, 0, wx.EXPAND, 0)
        # image_type_sizer.Add(empty_cell, 0, wx.EXPAND, 0)
        # image_type_sizer.Add(maximum_projection_radio_button, 0, wx.EXPAND, 0)
        # #
        # # image_type_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # # image_type_sizer.Add(wx.StaticText(self, label="Image type:"), 0, wx.EXPAND | wx.ALL, 10)
        # # image_type_sizer.Add(image_type_combobox, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 10)
        #
        # self.sizer = wx.BoxSizer(wx.VERTICAL)
        # self.sizer.Add(channel_sizer, 0, wx.EXPAND | wx.ALL, 10)
        # self.sizer.Add(image_type_sizer, 0, wx.EXPAND | wx.ALL, 10)
        # self.sizer.Add(wx.Slider(self, -1, value=0, minValue=0, maxValue=100,
        #                          name='Minimum'), 0, wx.EXPAND | wx.ALL, 10)
        # self.sizer.Add(wx.Slider(self, -1, value=0, minValue=0, maxValue=100,
        #                          name='Maximum'), 0, wx.EXPAND | wx.ALL, 10)
        # self.sizer.Add(wx.Button(self, -1, "Button 1"), 0, wx.EXPAND | wx.ALL, 10)
        # self.sizer.Add(wx.Button(self, -1, "Button 2"), 0, wx.EXPAND | wx.ALL, 10)
        # self.SetSizer(self.sizer)

    def save(self):
        if self.parent.save_path is not None:
            file_name = self.molecule.file.item().replace('\\' ,' - ')+f' - mol {self.molecule.molecule_in_file.item()}.png'
            file_path = self.parent.save_path.joinpath(file_name)
            self.figure.savefig(file_path, bbox_inches='tight')
        else:
            raise ValueError('No save_path set')








class TraceAnalysisFrame(QMainWindow):
    def __init__(self, parent=None, dataset=None, title='Traces', plot_variables=['intensity', 'FRET'],
                 ylims=[(0, 35000), (0, 1)], colours=[('g', 'r'), ('b')], save_path=None):
        wx.Frame.__init__(self, parent, title=title, size=(1400, 700))
        self.parent = parent
        self.dataset = dataset
        self.plot_variables = plot_variables
        self.ylims = ylims
        self.colours = colours

        if save_path is None:
            self.save_path = save_path
        else:
            self.save_path = Path(save_path)
        #self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.trace_panel = TraceAnalysisPanel(parent=self)
        # self.control_panel = ControlPanel(parent=self)
        self.Bind(wx.EVT_CHAR_HOOK, self.OnNavigationKey)

        self.molecule_index = 0

        self.Show()

    # @property
    # def molecules(self):
    #     return self._molecules
    #
    # @molecules.setter
    # def molecules(self, molecules):
    #     self._molecules = molecules


    @property
    def molecule_index(self):
        return self._molecule_index

    @molecule_index.setter
    def molecule_index(self, molecule_index):
        self._molecule_index = molecule_index
        self.molecule = self.dataset.isel(molecule=self.molecule_index)

    def next_molecule(self):
        if (self.molecule_index+1) < len(self.dataset.molecule):
            self.molecule_index += 1

    def previous_molecule(self):
        if self.molecule_index > 0:
            self.molecule_index -= 1

    def update_current_molecule(self):
        self.molecule_index = self.molecule_index

    @property
    def molecule(self):
        return self.panel.molecule

    @molecule.setter
    def molecule(self, molecule):
        self.trace_panel.molecule = molecule

    def OnNavigationKey(self, event):
        key_code = event.GetKeyCode()
        print(key_code)
        if key_code == 316: # Right arrow
            self.next_molecule()
        elif key_code == 314: # Left arrow
            self.previous_molecule()
        elif key_code == 32: # Spacebar
            self.dataset.selected[dict(molecule=self.molecule_index)] = ~self.dataset.selected[dict(molecule=1)]
            self.update_current_molecule()
        elif key_code == 83: # S
            self.trace_panel.save()

# class ControlPanel(wx.Panel):
# To file/molecule
# Selected button changes colour (spacebar)
# Y-axis limits
# Classification technique - change classification pannel based on specific technique

# class Threshold_classification_panel

# class HMM_classification_panel


class BlitManager:
    def __init__(self, canvas, animated_artists=()):
        """
        Parameters
        ----------
        canvas : FigureCanvasAgg
            The canvas to work with, this only works for sub-classes of the Agg
            canvas which have the `~FigureCanvasAgg.copy_from_bbox` and
            `~FigureCanvasAgg.restore_region` methods.

        animated_artists : Iterable[Artist]
            List of the artists to manage
        """
        self.canvas = canvas
        self._bg = None
        self._artists = []

        for a in animated_artists:
            self.add_artist(a)
        # grab the background on every draw
        self.cid = canvas.mpl_connect("draw_event", self.on_draw)

    def on_draw(self, event):
        """Callback to register with 'draw_event'."""
        cv = self.canvas
        if event is not None:
            if event.canvas != cv:
                raise RuntimeError
        self._bg = cv.copy_from_bbox(cv.figure.bbox)
        self._draw_animated()

    def add_artist(self, art):
        """
        Add an artist to be managed.

        Parameters
        ----------
        art : Artist

            The artist to be added.  Will be set to 'animated' (just
            to be safe).  *art* must be in the figure associated with
            the canvas this class is managing.

        """
        if art.figure != self.canvas.figure:
            raise RuntimeError
        art.set_animated(True)
        self._artists.append(art)

    def _draw_animated(self):
        """Draw all of the animated artists."""
        fig = self.canvas.figure
        for a in self._artists:
            fig.draw_artist(a)

    def update(self):
        """Update the screen with animated artists."""
        cv = self.canvas
        fig = cv.figure
        # paranoia in case we missed the draw event,
        if self._bg is None:
            self.on_draw(None)
        else:
            # restore the background
            cv.restore_region(self._bg)
            # draw all of the animated artists
            self._draw_animated()
            # update the GUI state
            cv.blit(fig.bbox)
        # let the GUI event loop process anything it has to do
        # cv.flush_events()

# class MainWindow(wx.Frame):
#    def __init__(self, parent, title):
#        wx.Frame.__init__(self, parent, title=title, size=(300, 700))
#        self.parent = parent
#        self.panel = TraceAnalysisPanel(parent=self)
#        # self.Bind(wx.EVT_CLOSE, self.OnClose)
#        self.Show()

if __name__ == "__main__":

    # # Check whether there is already a running QApplication (e.g., if running
    # # from an IDE).
    # qapp = QtWidgets.QApplication.instance()
    # if not qapp:
    #     qapp = QtWidgets.QApplication(sys.argv)
    #
    # app = ApplicationWindow()
    # app.show()
    # app.activateWindow()
    # app.raise_()
    # qapp.exec_()


    import trace_analysis as ta
    import os, sys
    mapping_path = Path(os.getcwd()).joinpath('trace_analysis').joinpath('mapping')
    sys.path.append(mapping_path)
    print(sys.path)
    from trace_analysis.experiment import Experiment
    exp = Experiment(r'D:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy')
    ds = exp.files[0].dataset


    from PySide2.QtWidgets import QApplication

    app = QApplication(sys.argv)
    frame = TracePlotWindow(ds)
        #, "Sample editor", plot_variables=['intensity', 'FRET'],  # 'classification'],
        #          ylims=[(0, 1000), (0, 1), (-1,2)], colours=[('g', 'r'), ('b'), ('k')])

    app.exec_()

    # # exp = ta.Experiment(r'D:\20200918 - Test data\Single-molecule data small')
    # #exp = ta.Experiment(r'P:\SURFdrive\Promotie\Data\Test data')
    # # exp = ta.Experiment(r'/Users/ivoseverins/SURFdrive/Promotie/Data/Test data')
    # # print(exp.files)
    # # m = exp.files[1].molecules[0]
    # # print(exp.files[2])
    # import xarray as xr
    # #file_paths = [p for p in exp.nc_file_paths if '561' in str(p)]
    # file_paths = [exp.nc_file_paths[0]]
    # with xr.open_mfdataset(file_paths, concat_dim='molecule', combine='nested') as ds:
    #     # ds_sel = ds.sel(molecule=ds.sequence_name=='HJ7_G')# .reset_index('molecule', drop=True) # HJ1_WT, HJ7_G116T
    #     app = wx.App(False)
    #     # app = wit.InspectableApp()
    #     frame = TraceAnalysisFrame(None, ds, "Sample editor", plot_variables=['intensity', 'FRET'], #'classification'],
    #              ylims=[(0, 1000), (0, 1), (-1,2)], colours=[('g', 'r'), ('b'), ('k')])
    #     # frame.molecules = exp.files[1].molecules
    #     print('test')
    #     import wx.lib.inspection
    #     wx.lib.inspection.InspectionTool().Show()
    #     app.MainLoop()





# Add time to existing .nc file
# for file in exp.files:
#     with xr.open_dataset(file.absoluteFilePath.with_suffix('.nc')) as ds:
#         i = ds.intensity.load()
#     test = i.assign_coords(time=file.movie.time)
#     test.to_netcdf(file.absoluteFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')


#
# from matplotlib import use
# use('TkAgg')
#
# import trace_analysis as ta
# exp = ta.Experiment(r'D:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy')
# #exp = ta.Experiment(r'J:\Ivo\20200221 - Magnetic tweezers setup (Old)\Data')
# # exp.files[-2].perform_mapping()
# # exp.files[-2].mapping.show_mapping_transformation()

# class B:
#     def __init__(self):
#         print('Badd')
#         super().__init__()
#
#
#
# class A:
#     def __init__(self):
#         print('A')
#
#
# def test(c):
#     return type(c.__name__, (c,B),{})
#
#
# @test
# class Bo(A):
#     def __init__(self):
#         print('Bo')
#         super().__init__()



# class B:
#     def __init__(self):
#         print('Badd')
#         super().__init__()
#
# # class PluginMetaClass(type):
# #     def __new__(cls, clsname, bases, attrs):
# #         bases_base = tuple(base for base in bases if not base.__name__ is clsname)
# #         attrs.pop('__qualname__')
# #         cls_base = type(clsname+'_base', bases_base, attrs)
# #         bases_main = tuple(base for base in bases if base.__name__ is clsname) + (cls_base,)
# #         return super().__new__(cls, clsname, bases_main, {})
# class PluginMetaClass(type):
#     def __new__(cls, clsname, bases_base, attrs):
#         # bases_base = tuple(base for base in bases if not base.__name__ is clsname)
#         attrs_base = attrs.copy()
#         attrs_base.pop('__qualname__')
#         #attrs_base.pop('__module__')
#         #attrs_base.pop('__classcell__')
#         cls_base = super().__new__(cls, clsname, bases_base, attrs_base)
#         #cls_base = type(clsname, bases_base, attrs)
#         added_bases = (B,)
#         bases_main = added_bases + (cls_base,)
#         test = super().__new__(cls, clsname+'main', bases_main,{})
#         print('test')
#         return test
#
# class A:
#     def __init__(self):
#         print('A')
#
# class Bo(A, metaclass=PluginMetaClass):
#     def __init__(self):
#         print('Bo')
#         super().__init__()
#



# exp = ta.Experiment(r'P:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy')
# # exp = ta.Experiment(r'D:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy')
# exp.files[-1].use_mapping_for_all_files()



# def add_class_to_class(base_class):
#     def add_class_to_class_decorator(added_class):
#         base_class.__bases__ += (added_class,)
#     return add_class_to_class_decorator
#
# @add_class_to_class(ta.File)
# class ExperimentPlugIn():
#     def test(self):
#         print(self.name)



# exp.files[0].find_coordinates()
#
#
# # #exp = ta.Experiment(r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20191209 - Single-molecule setup (TIR-I)')
# # exp.files[0].perform_mapping(transformation_type='nonlinear')
# #
# import matplotlib.pyplot as plt
# figure = plt.figure()
# #exp.files[0].show_average_image(figure=figure)
# plt.imshow(exp.files[0].movie.maximum_projection_image)
# exp.files[0].show_coordinates(figure=figure)
# #exp.files[0].mapping.show_mapping_transformation(figure=figure)



# exp.files[-1].use_mapping_for_all_files()

from trace_analysis.plotting import histogram
# exp.files[7].histogram(bins = 100, molecule_averaging=True, export=True)
# exp.histogram(bins = 100, molecule_averaging=True, export=True)
#
# import sys
# #sys.path.append(r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\fastqAnalysis')
# sys.path.append(r'D:\SURFdrive\Promotie\Code\Python\fastqAnalysis')
#
# from trace_analysis.traceAnalysisCode import Experiment
# from fastqAnalysis import FastqData
#
# from pathlib import Path # For efficient path manipulation
#
# path = Path(r'G:\Ivo\20190918 - Sequencer (MiSeq)\Analysis')
# #path = 'D:\\ivoseverins\\Desktop\\Sequencing data\\20180705\\'
# #path = 'C:\\Users\\Ivo Severins\\Desktop\\Sequencing data\\20180705\\'
# fileName = r'One_S1_L001_R1_001.fastq'
#
#
# data = FastqData(path.joinpath(fileName))
#
# data.selection(sequence = 'AA')
#
# data.matches_per_tile(sequence = 'TATCTGTATAATGAGAAATATGGAGTACAATTTTTTTTTTTTTTTTTTTT')









#import wx
#
#
#class OtherFrame(wx.Frame):
#    """
#    Class used for creating frames other than the main one
#    """
#
#    def __init__(self, title, parent=None):
#        wx.Frame.__init__(self, parent=parent, title=title)
#        self.Show()
#
#
#class MyPanel(wx.Panel):
#
#    def __init__(self, parent):
#        wx.Panel.__init__(self, parent)
#
#        btn = wx.Button(self, label='Create New Frame')
#        btn.Bind(wx.EVT_BUTTON, self.on_new_frame)
#        self.frame_number = 1
#
#    def on_new_frame(self, event):
#        title = 'SubFrame {}'.format(self.frame_number)
#        frame = OtherFrame(title=title)
#        self.frame_number += 1
#
#
#class MainFrame(wx.Frame):
#
#    def __init__(self):
#        wx.Frame.__init__(self, None, title='Main Frame', size=(800, 600))
#        panel = MyPanel(self)
#        self.Show()
#
#
#if __name__ == '__main__':
#    app = wx.App(False)
#    frame = MainFrame()
#    app.MainLoop()


# #!/usr/bin/env python
# import wx
# import wx.dataview
# import wx.lib.agw.aui as aui
# import os
#
# import wx.lib.agw.customtreectrl as CT
# #from traceAnalysisCode import Experiment
# import wx.lib.agw.hypertreelist as HTL
#
#
# import matplotlib as mpl
# from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
# from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
#
# from matplotlib import use
# use('WXAgg')
# from matplotlib import pyplot as plt
# #import matplotlib.pyplot as plt
#
#
#
# class MyFrame(wx.Frame):
#     """ We simply derive a new class of Frame. """
#     def __init__(self, parent, title):
#         wx.Frame.__init__(self, parent, title=title, size=(400,400))
#         tree_list = HTL.HyperTreeList(self)
#
#         tree_list.AddColumn("First column")
#
#         root = tree_list.AddRoot("Root")
#
#         parent = tree_list.AppendItem(root, "First child")
#         child = tree_list.AppendItem(parent, "First Grandchild")
#
#         tree_list.AppendItem(root, "Second child", ct_type=1)
#         self.Show(True)
#
# app = wx.App(False)
# frame = MyFrame(None, 'Small editor')
# app.MainLoop()