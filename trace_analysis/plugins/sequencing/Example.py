import numpy as np
from pathlib import Path # For efficient path manipulation
import matplotlib.pyplot as plt
from git import Repo

import trace_analysis as ta
#from fastqAnalysis import FastqData
from trace_analysis.mapping.geometricHashing import SequencingDataMapping
from trace_analysis.coordinate_transformations import transform
from trace_analysis.plotting import scatter_coordinates

# This part makes an analysis.txt file in the same directory as the this .py file. (please do not commit this)
# This file contains the git commit used and all differences to this commit.
# writePath = Path(Path(__file__).parent)
# writeFile = writePath.joinpath('git_repository_version_and_differences.txt')
# repo = Repo(Path(ta.__file__).parent.parent)
#
# with writeFile.open("a") as f:
#     f.write('------------------------------------------\n\n')
#     f.write(f'Trace_analysis version: {repo.head.object.hexsha} \n\n')
#     t = repo.head.commit.tree
#     f.write(repo.git.diff(t))
#     f.write('\n\n------------------------------------------\n\n')


# =================
# FLUORESCENCE DATA
# =================

# Import experiment
# IMPORTANT: uses only a single colour
experiment_path = r'Path_to_experiment_folder'
exp = ta.Experiment(experiment_path, colours=['r'])

# Find coordinates for the files that have no molecules yet
# (remove the if statement if you want all files to be reanalyzed)
for file in exp.files:
    #if len(file.molecules)==0:
        file.find_coordinates()
        file.extract_traces()

# Make a boxplot of the number of molecules per file
exp.boxplot_number_of_molecules()

# Select only molecules with a average intensity larger than a threshold
files = exp.files
for file in files:
    file.molecules = [molecule for molecule in file.molecules if (np.mean(molecule.intensity)>8500)]
# Select only files with at least two molecules
files = [file for file in files if (len(file.molecules) > 2)]

# Make a boxplot of the number of molecules per file
exp.boxplot_number_of_molecules()

# Show coordinates of a file on top of the image
exp.files[0].show_coordinates()
plt.imshow(exp.files[0].movie.average_image, vmin=200, vmax=400)

# ===============
# SEQUENCING DATA
# ===============

# Import sequencing data into experiment
fastqFilePath = Path(r'Path_to_fastq_file')
exp.import_sequencing_data(fastqFilePath)

# ============================
# MAPPING BY GEOMETRIC HASHING
# ============================

# Select only the sequences that are used for mapping and generate hash table
mapping_sequence = 'TATCTGTATAATGAGAAATATGGAGTACAATTTTTTTTTTTTTTTTTTTT'
number_of_allowed_mismatches = 0
# exp.select_sequencing_data_for_mapping(mapping_sequence, number_of_allowed_mismatches)
exp.generate_mapping_hashtable(mapping_sequence, number_of_allowed_mismatches,
                               imaged_surface=2, maximum_distance_tile=3000, tuple_size=4)

# For each file in experiment find a match
# NOTE: 'scale': [-1,1] means a reflection with respect to the y axis.
for file in exp.files:
    file.find_sequences(maximum_distance_file=1000, tuple_size=4, initial_transformation={'scale': [-1,1]},
                        hash_table_distance_threshold=0.01,
                        alpha=0.1, test_radius=10, K_threshold=10e9,
                        nearest_neighbour_match_distance_threshold=25)

# Make and export a sequencing match plot for files that have a sequencing match
for file in exp.files:
    if file.sequencing_match:
        file.plot_sequencing_match()


# # Calculate rotation and magnification for all found matches
# rotation = np.array([file.sequence_match.rotation for file in exp.files if file.sequence_match])
# magnification = np.array([file.sequence_match.magnification for file in exp.files if file.sequence_match])
#
# # Make a heatmap scatter plot
# from scipy.stats import gaussian_kde
# x = rotation
# y = magnification
#
# xy = np.vstack([x,y])
# z = gaussian_kde(xy,0.01)(xy)
# fig, ax = plt.subplots()
# ax.scatter(x, y, c=z, s=100, edgecolor='')
# plt.title('Matches')
# plt.xlabel('Rotation (deg)')
# plt.ylabel('Magnfication')
# plt.show()
#
# # Calculate mean rotation and magnification and their standard deviations
# mean_rotation = np.mean(rotation)
# std_rotation = np.std(rotation)
# mean_magnification = np.mean(magnification)
# std_magnification = np.std(magnification)
# print(mean_rotation,std_rotation)
# print(mean_magnification,std_magnification)

