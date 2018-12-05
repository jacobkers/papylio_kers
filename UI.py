# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:44:52 2018

@author: ivoseverins
"""

#!/usr/bin/env python
import wx
import os
from traceAnalysisCode import Experiment


        

class MainFrame(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title='Trace Analysis', size=(300,200))
        
        # Status bar
        self.CreateStatusBar()
        
        # File menu in Menu bar
        fileMenu = wx.Menu()
        fileMenuOpen = fileMenu.Append(wx.ID_OPEN, "&Open", "Open directory")
        fileMenuAbout = fileMenu.Append(wx.ID_ABOUT, "&About")
        fileMenuExit = fileMenu.Append(wx.ID_EXIT, "&Exit")
        
        self.Bind(wx.EVT_MENU, self.OnOpen, fileMenuOpen)
        self.Bind(wx.EVT_MENU, self.OnAbout, fileMenuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, fileMenuExit)

        # Menu bar        x
        menuBar = wx.MenuBar()
        menuBar.Append(fileMenu,"&File")
        self.SetMenuBar(menuBar)
        
        
        # TreeCtrl
        self.tree = wx.TreeCtrl(self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TR_HAS_BUTTONS | wx.TR_MULTIPLE)
        #self.Bind(wx.EVT_TREE_SEL_CHANGED, self.Test, self.tree)
        self.tree.Bind(wx.EVT_LEFT_DOWN, self.Test)
        #panel = TreePanel(self)
        
        
           
        self.Show(True)
    
    # File menu event handlers
    def OnOpen(self,event):
        self.experimentRoot = ''
        dlg = wx.DirDialog(self, "Choose a directory", self.experimentRoot)
        if dlg.ShowModal() == wx.ID_OK:
            self.experimentRoot = dlg.GetPath()
            print(self.experimentRoot)
            exp = Experiment(self.experimentRoot)
            
            self.root = self.tree.AddRoot(exp.name)
            
            for file in exp.files:
                self.tree.AppendItem(self.root, file.name)
            
            self.tree.ExpandAll()
    
    def OnAbout(self,event):
        dlg = wx.MessageDialog(self, 'Software for trace analysis', 'About Trace Analysis', wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def OnExit(self,event):
        self.Close(True) # Close program

    
    # TreeCtrl event handlers
    def Test(self, event):
         ht_item, ht_flags = self.tree.HitTest(event.GetPosition())
         print(ht_flags)
         
         if not self.tree.IsSelected(ht_item):
             self.tree.SelectItem(ht_item)
         else:
             self.tree.SelectItem(ht_item, select = False)
         
         

app = wx.App(False)
frame = MainFrame(None, 'MainFrame')
app.MainLoop()



del app


"""

class MyTree(wx.TreeCtrl):
    def __init__(self, parent, id, pos, size, style):
        wx.TreeCtrl.__init__(self, parent, id, pos, size, style)
        
class TreePanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        self.tree = MyTree(self, wx.ID_ANY, wx.DefaultPosition, (300,200), wx.TR_HAS_BUTTONS)
        
        self.root = self.tree.AddRoot('Something goes here')
        self.tree.SetPyData(self.root, ('key', 'value'))
        self.tree.AppendItem(self.root, 'Operating Systems')
        self.tree.Expand(self.root)
        


"""