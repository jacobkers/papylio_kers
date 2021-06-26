# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:44:52 2018

@author: ivoseverins
"""
import wx
import wx.lib.mixins.inspection as wit
import matplotlib as mpl
mpl.use('WXAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar



class TraceAnalysisFrame(wx.Frame):
    def __init__(self, molecules, parent=None, title='Traces'):
        wx.Frame.__init__(self, parent, title=title, size=(1400,700))
        self.parent = parent
        self.panel = TraceAnalysisPanel(parent=self)
        self.molecules = molecules
        #self.Bind(wx.EVT_CLOSE, self.OnClose)

        self.Bind(wx.EVT_CHAR_HOOK, self.OnNavigationKey)

        self.Show()

    @property
    def molecules(self):
        return self._molecules

    @molecules.setter
    def molecules(self, molecules):
        self._molecules = molecules
        self.molecule_index = 0

    @property
    def molecule_index(self):
        return self._molecule_index

    @molecule_index.setter
    def molecule_index(self, molecule_index):
        self._molecule_index = molecule_index
        self.molecule = self.molecules[molecule_index]

    def next_molecule(self):
        if (self.molecule_index+1) < len(self.molecules):
            self.molecule_index += 1

    def previous_molecule(self):
        if self.molecule_index > 0:
            self.molecule_index -= 1

    @property
    def molecule(self):
        return self.panel.molecule

    @molecule.setter
    def molecule(self, molecule):
        self.panel.molecule = molecule

    def OnNavigationKey(self, event):
        key_code = event.GetKeyCode()
        print(key_code)
        if key_code == 316:
            self.next_molecule()
        elif key_code == 314:
            self.previous_molecule()


class TraceAnalysisPanel(wx.Panel):
    def __init__(self, parent, id=-1, dpi=None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, size=(1400,700), **kwargs)


        self.figure = mpl.figure.Figure(dpi=dpi)#, figsize=(2, 2))
        self.axes=self.figure.subplots(2,1)
        # self.figure = plt.Figure(dpi=dpi, figsize=(2,2))
        #
        # self.axis = self.figure.gca()

        #self.figure, self.axes = mpl.figure.Figure().subplots(2,1)

        self.canvas = FigureCanvas(self, -1, self.figure)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        self.SetSizer(sizer)

        self._molecule = None

        self.artists = []

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, molecule):
        self._molecule = molecule

        if not self.artists:
            self.artists += self.axes[0].plot(molecule.intensity[0, :].T, c='g')
            self.artists += self.axes[0].plot(molecule.intensity[1, :].T, c='r')
            self.artists += self.axes[1].plot(molecule.E().T, c='b')
            #self.axes[1].plot(molecule.E(), animate=True)
            self.bm = BlitManager(self.canvas, self.artists)
            self.canvas.draw()

        # for axis in self.axes:
        #     axis.cla()

        self.artists[0].set_ydata(molecule.intensity[0,:])
        self.artists[1].set_ydata(molecule.intensity[1,:])
        self.artists[2].set_ydata(molecule.E())

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

import trace_analysis as ta
#exp = ta.Experiment(r'D:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy')
exp = ta.Experiment(r'D:\20200918 - Test data\Single-molecule data small')
# print(exp.files)
# m = exp.files[1].molecules[0]
# print(exp.files[2])

app = wx.App(False)
# app = wit.InspectableApp()
frame = TraceAnalysisFrame(None, m, "Sample editor")
# frame.molecules = exp.files[1].molecules
print('test')
import wx.lib.inspection
wx.lib.inspection.InspectionTool().Show()
app.MainLoop()






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