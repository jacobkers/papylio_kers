# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 01:02:04 2019

@author: iason
"""

# sized_dialog.py

import wx

class Example(wx.Frame):

   def __init__(self, parent, title):
      super(Example, self).__init__(parent, title = title)

      self.InitUI()
      self.Centre()
      self.Show()

   def InitUI(self):

      panel = wx.Panel(self)
      sizer = wx.GridBagSizer(0,0)

      text = wx.StaticText(panel, label = "Name:")
      sizer.Add(text, pos = (0, 0), flag = wx.ALL, border = 5)

      tc = wx.TextCtrl(panel)
      sizer.Add(tc, pos = (0, 1), span = (1, 2), flag = wx.EXPAND|wx.ALL, border = 5)

      text1 = wx.StaticText(panel, label = "address")
      sizer.Add(text1, pos = (1, 0), flag = wx.ALL, border = 5)

      tc1 = wx.TextCtrl(panel,style = wx.TE_MULTILINE)
      sizer.Add(tc1, pos = (1,1), span = (1, 3), flag = wx.EXPAND|wx.ALL, border = 5)

      text2 = wx.StaticText(panel,label = "age")
      sizer.Add(text2, pos = (2, 0), flag = wx.ALL, border = 5)

      tc2 = wx.TextCtrl(panel)
      sizer.Add(tc2, pos = (2,1), flag = wx.ALL, border = 5)

      text3 = wx.StaticText(panel,label = "Mob.No")
      sizer.Add(text3, pos = (2, 2), flag = wx.ALIGN_CENTER|wx.ALL, border = 5)

      tc3 = wx.TextCtrl(panel)
      sizer.Add(tc3, pos = (2,3),flag = wx.EXPAND|wx.ALL, border = 5)

      text4 = wx.StaticText(panel, label = "Description")
      sizer.Add(text4, pos = (3, 0), flag = wx.ALL, border = 5)

      tc4 = wx.TextCtrl(panel,style =  wx.TE_MULTILINE)
      sizer.Add(tc4, pos = (3,1), span = (1,3), flag = wx.EXPAND|wx.ALL, border = 5)
      sizer.AddGrowableRow(3)

      buttonOk = wx.Button(panel, label = "Ok")
      buttonClose = wx.Button(panel, label = "Close" )

      sizer.Add(buttonOk, pos = (4, 2),flag = wx.ALL, border = 5)
      sizer.Add(buttonClose, pos = (4, 3), flag = wx.ALL, border = 5)

      panel.SetSizerAndFit(sizer)

app = wx.App()
Example(None, title = 'GridBag Demo')
app.MainLoop()