# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:14:28 2019

@author: ikatechis
"""

from trace_analysis import Experiment
import os
from pathlib import Path, PureWindowsPath
from interactiveAnalysis import InteractivePlot

mainPath = 'G:/SM-data/20191101_dcas9_flow_DNA04_DNA20/'
mainPath += '#4.20_streptavidin_0.5nM_dcas9-crRNA-Cy5_10nM_DNA04-Cy3_G_flow'
mainPath = Path(mainPath)
exp = Experiment(mainPath)
file = exp.files[1]
i = InteractivePlot(file)
i.plot_initialize()
i.plot_molecule()

#plt.show()