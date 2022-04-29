import os
import numpy as np
from trace_analysis import Experiment # Import before matplotlib

from matplotlib import pyplot as plt
from sys import platform

# if platform == "darwin":
#     from matplotlib import use
#     use('WXAgg')

# from trace_analysis import InteractivePlot
from trace_analysis.image_adapt.autoconfig import autoconfig_AND_perform_mapping, autoconfig


# Define path to data, replace by your own directory
# mainPath = r'D:\pathToDataFolder'
mainPath = None # If mainPath is None, a window will appear to select a directory

# Initialize an experiment
exp = Experiment(mainPath)

# Print files in experiment
exp.print_files()

# Perform mapping
perform_new_mapping = False
if perform_new_mapping:
    mapping_file_index = int(input('\nEnter file number for mapping... '))
    mapping_file = exp.files[mapping_file_index]

    # autoconfig_AND_perform_mapping(mapping_file_index, mainPath)

    mapping_file.perform_mapping()
    figure_mapping = plt.figure(101)
    mapping_file.show_image(figure=figure_mapping)
    mapping_file.mapping.show_mapping_transformation(figure=figure_mapping, show_source=True)
    plt.show(block=False)
    plt.pause(0.1)

# Run for specific files
file_indices_of_interest = range(2, 110)

selected_files = [file for file_index, file in enumerate(exp.files) if file_index in file_indices_of_interest]

for file in selected_files:
    # autoconfig(exp, file, opt='find_coordinates')

    # Do not run for the mapping file
    if file.is_mapping_file:
        continue

    # Do not run if there are already molecules loaded, i.e. if pks or traces files are present
    # if file.molecules:
    #     continue

    print(f'\nWorking on {file}')
    file.find_coordinates()
    file.extract_traces()


# show movie
for file in selected_files:
    figure_handle_movie = plt.figure()
    # file.show_image(figure=figure_handle_movie)
    file.show_image(figure=figure_handle_movie, vmin=120, vmax=200)
    if file.is_mapping_file:
        file.mapping.show_mapping_transformation(figure=figure_handle_movie)
    else:
        if file.coordinates.size > 0:
            file.show_coordinates(figure=figure_handle_movie)
plt.show()
plt.pause(0.1)

# Show interactive plot
file_index = 3
molecules = exp.files[0].molecules
int_plot = InteractivePlot(molecules)
int_plot.plot_initialize()
int_plot.plot_molecule()
