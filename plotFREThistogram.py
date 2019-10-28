# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 11:09:03 2019

@author: pimam
"""

import os
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path, PureWindowsPath
import traceAnalysisCode as trace_ana

mainPath = PureWindowsPath('C:\\Users\\pimam\\Documents\\MEP\\tracesfiles')
mainPath = Path(mainPath)
exp = trace_ana.Experiment(mainPath)
file = exp.files[0]
nbins=100
plt.figure(1)
trace_ana.histogram(file.molecules, nbins, None)
plt.title(f'FRET histogram bins:{nbins}')
plt.xlabel('FRET')
plt.ylabel('count')
plt.savefig(f'{file.name}_FRET_hist.png', facecolor='white', dpi=300)