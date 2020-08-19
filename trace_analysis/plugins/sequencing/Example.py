import numpy as np
from pathlib import Path # For efficient path manipulation
import matplotlib.pyplot as plt
from git import Repo
import time

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
start = time.time()
exp.generate_mapping_hashtable(mapping_sequence, number_of_allowed_mismatches,
                               imaged_surface=2, maximum_distance_tile=3000, tuple_size=4)
end = time.time()
print(f'Hashtable generation: {end - start} seconds')

# Or just use the coordinates from the loc files
# path_to_loc_folder = Path(r'Path_to_loc_folder')
# tile_coordinate_files = ['2101.loc', '2102.loc', '2103.loc', '2104.loc']
# tile_coordinate_sets = [np.loadtxt(path_to_loc_folder.joinpath(f)) for f in tile_coordinate_files]
# start = time.time()
# exp.generate_mapping_hashtable_from_coordinate_set(tile_coordinate_sets, maximum_distance_tile=3000, tuple_size=4)
# end = time.time()
# print(f'Hashtable generation: {end - start} seconds')

# For each file in experiment find a match
# NOTE: 'scale': [-1,1] means a reflection with respect to the y axis.

start = time.time()
for file in exp.files:
    file.sequencing_match = None
    file.find_sequences(maximum_distance_file=1000, tuple_size=4, initial_transformation={'scale': [-1,1]},
                        hash_table_distance_threshold=0.01,
                        alpha=0.1, test_radius=10, K_threshold=10e6, # original K_threshold = 10e9
                        magnification_range=[3.3,3.4], rotation_range=[-1,1])
    print(file)
end = time.time()
print(f'Matching: {end - start} seconds')

matched_files = [file for file in exp.files if file.sequencing_match]
number_of_matches = len(matched_files)

# Improve mapping by performing a linear least-squares fit on all nearest neighbours within the distance threshold
for file in matched_files:
    nearest_neighbour_match_distance_threshold = 25
    file.sequencing_match.nearest_neighbour_match(nearest_neighbour_match_distance_threshold)

# Import sequencing data and locations into the files
for file in matched_files:
    file.get_all_sequences_from_sequencing_data()
    # Export pks file with all molecules that have a linked sequence
    file.export_pks_file()


# Make and export a sequencing match plot for files that have a sequencing match
for file in matched_files: #exp.files # files_correct: #exp.files:
    if file.sequencing_match:
        file.plot_sequencing_match()
        # plt.close()

# =================================
# ASSESSING OVERALL MAPPING QUALITY
# =================================

# Calculate rotation and magnification for all found matches
rotation = np.array([file.sequencing_match.rotation for file in exp.files if file.sequencing_match])
magnification = np.array([file.sequencing_match.magnification for file in exp.files if file.sequencing_match])
mean_magnification = np.mean(magnification, axis=1)

# Make a heatmap scatter plot
from scipy.stats import gaussian_kde
x = rotation
y = mean_magnification

xy = np.vstack([x,y])
z = gaussian_kde(xy,0.01)(xy)
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=100, edgecolor='')
plt.title('Matches')
plt.xlabel('Rotation (deg)')
plt.ylabel('Magnfication (average over two directions)')
plt.show()

# Make a heatmap scatter plot
from scipy.stats import gaussian_kde
selection = ((rotation>179) | (rotation<-179)) & (mean_magnification>3.3) & (mean_magnification<3.4)
x = rotation[selection]
y = mean_magnification[selection]
xy = np.vstack([x,y])
z = gaussian_kde(xy,0.01)(xy)
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=100, edgecolor='')
plt.title('Matches')
plt.xlabel('Rotation (deg)')
plt.ylabel('Magnfication (average over two directions)')
plt.show()


# Calculate mean rotation and magnification and their standard deviations
mean_rotation = np.mean(rotation)
std_rotation = np.std(rotation)
mean_magnification = np.mean(magnification)
std_magnification = np.std(magnification)
print(mean_rotation,std_rotation)
print(mean_magnification,std_magnification)


files_correct = [file for i, file in enumerate(matched_files) if selection[i]]
files_incorrect = [file for i, file in enumerate(matched_files) if not selection[i]]

# Calculate mean hash table distances
hash_table_distances_correct = [file.sequencing_match.hash_table_distance for file in files_correct]
hash_table_distances_incorrect = [file.sequencing_match.hash_table_distance for file in files_incorrect]

print(np.mean(hash_table_distances_correct))
print(np.mean(hash_table_distances_incorrect))
print(np.std(hash_table_distances_correct))
print(np.std(hash_table_distances_incorrect))

# Number source tuples checked in hash table
tuples_checked_in_hash_table_correct = [file.sequencing_match.hash_table_distances_checked for file in files_correct]
tuples_checked_in_hash_table_incorrect = [file.sequencing_match.hash_table_distances_checked for file in files_incorrect]

print(np.mean(tuples_checked_in_hash_table_correct))
print(np.mean(tuples_checked_in_hash_table_incorrect))
print(np.std(tuples_checked_in_hash_table_correct))
print(np.std(tuples_checked_in_hash_table_incorrect))

# Number of source tuples that were checked with the destination point set
tuples_checked_with_destination_correct = [file.sequencing_match.tuples_checked for file in files_correct]
tuples_checked_with_destination_incorrect = [file.sequencing_match.tuples_checked for file in files_incorrect]

print(np.mean(tuples_checked_with_destination_correct))
print(np.mean(tuples_checked_with_destination_incorrect))
print(np.std(tuples_checked_with_destination_correct))
print(np.std(tuples_checked_with_destination_incorrect))

# Number of points in file
number_of_molecules_correct = [file.number_of_molecules for file in files_correct]
number_of_molecules_incorrect = [file.number_of_molecules for file in files_incorrect]

print(np.mean(number_of_molecules_correct))
print(np.mean(number_of_molecules_incorrect))
print(np.std(number_of_molecules_correct))
print(np.std(number_of_molecules_incorrect))