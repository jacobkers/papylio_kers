import os
import numpy as np
from matplotlib import pyplot as plt

from trace_analysis import Experiment
from trace_analysis import InteractivePlot


# Define path to data, replace by your own directory
mainPath = r'D:\pathToDataFolder'
mainPath=r'C:\MargreetDATA\werk\TUDelft\python\CMJ trace analysis\development_not_to_be_included_in_git\data Pim' 
mapping_file_index = 0# works fine for good and bad config vaues

mainPath=r'C:\MargreetDATA\werk\TUDelft\python\CMJ trace analysis\development_not_to_be_included_in_git\example_data_SHK\set 2 high FRET\01 50pM branch1-3 5nM MF34 50nM MF68' # gives error in polywarp cannot unpack non-iterable NoneType object
mapping_file_index = 0#above folder now runs with default config file (and run loop over try multiple times)

mainPath=r'C:\MargreetDATA\werk\TUDelft\python\CMJ trace analysis\development_not_to_be_included_in_git\20200227 FRET exchange 2 same targets in mol 100mM MgCl2\1. 100pM MF133 + 5nM MF34 + 500mM NaCl + 100mM MgCl2' # AttributeError: 'NoneType' object has no attribute 'make_average_image'
mapping_file_index = 2 # rename hel12.pks, otherwise error on nb of molecules
# Initialize an experiment
exp = Experiment(mainPath)

# Print files in experiment
print(exp.files)

# Perform mapping

mapping_file = exp.files[mapping_file_index]

# try 1 with current config
from trace_analysis.image_adapt.autoconfig  import autoconfig_AND_perform_mapping 
autoconfig_AND_perform_mapping(mapping_file_index, mainPath)

        
#mapping_file.use_mapping_for_all_files()
#
## Run for all files
#for file in exp.files:
#    if 0:
#        autoconfig(exp,file, opt='find_coordinates')
#    # Do not run for the mapping file
#    if file.is_mapping_file: continue
#
#    # Do not run if there are already molecules loaded, i.e. if pks or traces files are present
#    if file.molecules: continue
#
#    file.find_coordinates()
#    file.extract_traces()

### Show interactive plot
#file_index = 3
#molecules = exp.files[0].molecules
#int_plot = InteractivePlot(molecules)
#int_plot.plot_initialize()
#int_plot.plot_molecule()





