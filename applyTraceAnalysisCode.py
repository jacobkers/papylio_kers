# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:50:59 2018

@author: ivoseverins
"""



import os

# Change path to traceAnalysis directory (assuming this file is in that directory)
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from traceAnalysisCode import Experiment

# Define path to data, replace by your own directory
# (Now it is set to the twoColourExampleData folder in the traceAnalysis directory)
mainPath = './twoColourExampleData'


# Initialize an experiment
exp = Experiment(mainPath)

