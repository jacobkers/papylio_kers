# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:44:52 2018

@author: ivoseverins
"""

#!/usr/bin/env python
import wx
import wx.dataview
import os
from traceAnalysisCode import Experiment


        

class MainFrame(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title='Trace Analysis', size=(300,500))
        
        #self.panel1 = wx.Panel(self, wx.ID_ANY, size = (200,200))
        #self.panel2 = wx.Panel(self, wx.ID_ANY, size = (200,200))
        
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
        
        
        
        
        # TreeListCtrl
        self.tree = wx.dataview.TreeListCtrl(self, wx.ID_ANY, wx.DefaultPosition, wx.Size(200,300), 
                                             wx.dataview.TL_CHECKBOX | wx.dataview.TL_MULTIPLE)
        self.Bind(wx.dataview.EVT_TREELIST_ITEM_CHECKED, self.Test, self.tree)
        #self.Bind(wx.EVT_TREE_SEL_CHANGED, self.Test, self.tree)
        ##self.tree.Bind(wx.EVT_LEFT_DOWN, self.Test)
        #panel = TreePanel(self)
        
        #test = wx.Button(self, -1, 'Large button')
        #test = wx.Button(self, -1, 'Large button')
        
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(self.tree, 0,0,0)
        box.Add(wx.Button(self, -1, 'Small button'), 0, 0, 0)
        box.Add(wx.Button(self, -1, 'Large button'), 0, 0, 0)
        #box.Add(self.tree, 1,0,0)
        #box.Add(self.panel2, 1,0,0)
        self.SetSizerAndFit(box)
           
        self.Show(True)
        
        self.createTree(r'D:\ivoseverins\OneDrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\HJ A')
        
        print('done')
        
    # File menu event handlers
    def OnOpen(self,event):
        self.experimentRoot = ''
        dlg = wx.DirDialog(self, "Choose a directory", self.experimentRoot)
        if dlg.ShowModal() == wx.ID_OK:
            self.createTree(dlg.GetPath())
#            self.experimentRoot = dlg.GetPath()
#            print(self.experimentRoot)
#            exp = Experiment(self.experimentRoot)
#            
#            self.tree.AppendColumn('Files')
#            self.experimentRoot = self.tree.AppendItem(self.tree.GetRootItem(),exp.name)
#            
#            for file in exp.files:
#                self.tree.AppendItem(self.experimentRoot, file.name)
#            
#            self.tree.Expand(self.experimentRoot)
    
    # Temporary function to automate data import
    def createTree(self, experimentRoot):
        self.experimentRoot = experimentRoot
        print(self.experimentRoot)
        exp = Experiment(self.experimentRoot)
        
        self.tree.AppendColumn('Files')
        self.experimentRoot = self.tree.AppendItem(self.tree.GetRootItem(),exp.name)
        
        for file in exp.files:
            self.tree.AppendItem(self.experimentRoot, file.name, data = file)
        
        self.tree.Expand(self.experimentRoot)
    #### End of temporary function
    
    def OnAbout(self,event):
        dlg = wx.MessageDialog(self, 'Software for trace analysis', 'About Trace Analysis', wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def OnExit(self,event):
        self.Close(True) # Close program

    
    # TreeListCtrl event handlers
    def Test(self, event):
        item = event.GetItem()
        newItemCheckedState = bool(self.tree.GetCheckedState(item))
        file = self.tree.GetItemData(item)
        file.isSelected = newItemCheckedState
        
         

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