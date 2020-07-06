import os
import numpy as np
from matplotlib import pyplot as plt

from trace_analysis import Experiment
from trace_analysis import InteractivePlot


# Define path to data, replace by your own directory
mainPath = r'D:\pathToDataFolder'
mainPath=r'C:\MargreetDATA\werk\TUDelft\python\CMJ trace analysis\development_not_to_be_included_in_git\data Pim'
# Initialize an experiment
exp = Experiment(mainPath)

# Print files in experiment
print(exp.files)

# Perform mapping
mapping_file_index = -1
mapping_file = exp.files[mapping_file_index]

# now you need autoconfig:
if 1:
    from trace_analysis.image_adapt.autoconfig  import autoconfig 
    autoconfig(exp,mapping_file, opt='mapping')
    print('autoconfig done')

mapping_file.perform_mapping()
mapping_file.use_mapping_for_all_files()

# Run for all files
for file in exp.files:
    if 1:
        autoconfig(exp,file, opt='find_coordinates')
    # Do not run for the mapping file
    if file.is_mapping_file: continue

    # Do not run if there are already molecules loaded, i.e. if pks or traces files are present
    if file.molecules: continue

    file.find_coordinates()
    file.extract_traces()

### Show interactive plot
file_index = 3
molecules = exp.files[0].molecules
int_plot = InteractivePlot(molecules)
int_plot.plot_initialize()
int_plot.plot_molecule()





