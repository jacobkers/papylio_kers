# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 03:55:04 2019

@author: iason
"""
import numpy as np
import traceAnalysisCode as analysis
import pandas as pd
import os

#main_path = './traces'
#os.chdir(os.path.dirname(os.path.abspath(__file__)))

main_path = './iasonas/cas9_simulations/DNA11-20_30exposure'
exp = analysis.Experiment(main_path, 30)


for file in exp.files:
    file.autoThreshold('red', threshold=80)
