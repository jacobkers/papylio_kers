import numpy as np
from pathlib import Path # For efficient path manipulation
import matplotlib.pyplot as plt
from git import Repo
import time
from skimage.transform import AffineTransform

import trace_analysis as ta
from trace_analysis.plugins.sequencing.fastqAnalysis import FastqData

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


# -----------------------------------
# Load experiment
# -----------------------------------
experiment_path = r'D:\20200918 - Test data\Single-molecule data small'
exp = ta.Experiment(experiment_path)

# In case movies are rotated
for file in exp.files:
    try:
        file.movie.rot90 = 1
    except AttributeError:
        pass

# -----------------------------------
# Define files
# -----------------------------------
files_green_laser = [file for file in exp.files if ('561' in file.name) and ('640' not in file.name)]
files_red_laser = [file for file in exp.files if '642' in file.name]

file_green = files_green_laser[0]
file_red = files_red_laser[0]

mapping_file = [file for file in exp.files if 'DUAL' in file.name][0]
mapping_file.use_mapping_for_all_files()

# -----------------------------------
# Load sequencing data
# -----------------------------------
# Possibly we can make something in the future that recognizes the R1, R2 and I1 tags in the filenames.
fastq_files = sorted(Path(r'D:\20200918 - Test data\Sequencing data').glob('*.fastq'))

FastqDataObjects = [FastqData(file_path) for file_path in fastq_files]

for fastqData in FastqDataObjects:
    print(len(fastqData))

R1=sum(FastqDataObjects[1::3])
R2=sum(FastqDataObjects[2::3])
I1=sum(FastqDataObjects[0::3])

sequencing_data = {
    'Read1': R1,
    'Read2': R2,
    'Read1_HJ': R1,
    'Read2_HJ': R2,
    'Index1': I1,
    'IndexL1': R1,
    'IndexL2': R2
}

exp.sequencing_data = R1

# -----------------------------------
# Select sequencing data for mapping
# -----------------------------------
mapping_sequence = 'ACTGACTGTAACAACAACAACAATAACAACAACAACAATAACAACAACAACAAT'
number_of_allowed_mismatches = 2
exp.select_sequencing_data_for_mapping(mapping_sequence, number_of_allowed_mismatches)

# -----------------------------------
# Find high-intensity coordinates
# -----------------------------------
configuration = exp.configuration['find_coordinates'].copy()
configuration['peak_finding']['minimum_intensity_difference'] = 1500

for file in files_red_laser:
    file.find_coordinates(configuration=configuration)


# -----------------------------------
# Geometric hashing 2.0
# -----------------------------------
# # This part is likely not yet working well
# start = time.time()
# exp.generate_mapping_hashtable(mapping_sequence, number_of_allowed_mismatches,
#                                imaged_surface=2, maximum_distance_tile=3000, tuple_size=4)
# end = time.time()
# print(f'Hashtable generation: {end - start} seconds')
#
# # Or just use the coordinates from the loc files
# # path_to_loc_folder = Path(r'Path_to_loc_folder')
# # tile_coordinate_files = ['2101.loc', '2102.loc', '2103.loc', '2104.loc']
# # tile_coordinate_sets = [np.loadtxt(path_to_loc_folder.joinpath(f)) for f in tile_coordinate_files]
# # start = time.time()
# # exp.generate_mapping_hashtable_from_coordinate_set(tile_coordinate_sets, maximum_distance_tile=3000, tuple_size=4)
# # end = time.time()
# # print(f'Hashtable generation: {end - start} seconds')
#
# # For each file in experiment find a match
# # NOTE: 'scale': [-1,1] means a reflection with respect to the y axis.
#
# start = time.time()
# for file in exp.files:
#     file.sequencing_match = None
#     file.find_sequences(maximum_distance_file=1000, tuple_size=4, initial_transformation={'scale': [-1,1]},
#                         hash_table_distance_threshold=0.01,
#                         alpha=0.1, test_radius=10, K_threshold=10e6, # original K_threshold = 10e9
#                         magnification_range=[3.3,3.4], rotation_range=[-1,1])
#     print(file)
# end = time.time()
# print(f'Matching: {end - start} seconds')
#
# matched_files = [file for file in exp.files if file.sequencing_match]
# number_of_matches = len(matched_files)
#
# # Improve mapping by performing a linear least-squares fit on all nearest neighbours within the distance threshold
# for file in matched_files:
#     nearest_neighbour_match_distance_threshold = 25
#     file.sequencing_match.nearest_neighbour_match(nearest_neighbour_match_distance_threshold)
#
# # Make and export a sequencing match plot for files that have a sequencing match
# for file in matched_files: #exp.files # files_correct: #exp.files:
#     if file.sequencing_match:
#         file.plot_sequencing_match()
#         # plt.close()

# -----------------------------------
# ASSESSING OVERALL MAPPING QUALITY
# -----------------------------------
#
# # Calculate rotation and magnification for all found matches
# rotation = np.array([file.sequencing_match.rotation for file in exp.files if file.sequencing_match])
# magnification = np.array([file.sequencing_match.magnification for file in exp.files if file.sequencing_match])
# mean_magnification = np.mean(magnification, axis=1)
#
# # Make a heatmap scatter plot
# from scipy.stats import gaussian_kde
# x = rotation
# y = mean_magnification
#
# xy = np.vstack([x,y])
# z = gaussian_kde(xy,0.01)(xy)
# fig, ax = plt.subplots()
# ax.scatter(x, y, c=z, s=100, edgecolor='')
# plt.title('Matches')
# plt.xlabel('Rotation (deg)')
# plt.ylabel('Magnfication (average over two directions)')
# plt.show()
#
# # Make a heatmap scatter plot
# from scipy.stats import gaussian_kde
# selection = ((rotation>179) | (rotation<-179)) & (mean_magnification>3.3) & (mean_magnification<3.4)
# x = rotation[selection]
# y = mean_magnification[selection]
# xy = np.vstack([x,y])
# z = gaussian_kde(xy,0.01)(xy)
# fig, ax = plt.subplots()
# ax.scatter(x, y, c=z, s=100, edgecolor='')
# plt.title('Matches')
# plt.xlabel('Rotation (deg)')
# plt.ylabel('Magnfication (average over two directions)')
# plt.show()


# -----------------------------------
# Geometric hashing 3.0
# -----------------------------------

# Define initial file transformation
initial_magnification = np.array([3.71, -3.71])
initial_rotation = 0.4

initial_file_transformation = AffineTransform(matrix=None, scale=initial_magnification,
                                                rotation=initial_rotation/360*np.pi*2,
                                                shear=None, translation=None)

# TODO: Add timer to generate_mapping_hashtable and find_sequences methods, by making a decorator function. [IS: 10-08-2020]
start = time.time()
exp.generate_mapping_hashtable3(mapping_sequence, number_of_allowed_mismatches,
                               imaged_surface=1, initial_file_transformation=initial_file_transformation)
end = time.time()
print(f'Hashtable generation: {end - start} seconds')

matched_files = []
start = time.time()
for file in files_red_laser:
    if (file.number_of_molecules>4):
        print(file)
        file.sequencing_match = None
        # file.find_sequences3(distance=15, alpha=0.8, sigma=5, K_threshold=10**7)
        file.find_sequences3(distance=10, alpha=0.9, sigma=5, K_threshold=10**9, channel=1)
        if file.sequencing_match:
            matched_files.append(file)
end = time.time()
print(f'Matching: {end - start} seconds')

for file in matched_files:
    print(file.name + ' __ ' + str(file.sequencing_match.tile))

# Improve mapping by performing a linear least-squares fit on all nearest neighbours within the distance threshold
for file in matched_files:
    print(file)
    # try:
    nearest_neighbour_match_distance_threshold = 25
    file.sequencing_match.nearest_neighbour_match(nearest_neighbour_match_distance_threshold, transformation_type='linear')
    # except:
    #     print(file)

exp.show_sequencing_matches()

for file in matched_files:
    file.plot_sequencing_match()
    file.export_sequencing_match()


# -----------------------------------
# Find sequencing matches using stage coordinates
# -----------------------------------
# Find mapping between stage coordinates and sequencing coordinates
exp.map_sequencing_and_stage_coordinates()
exp.show_stage_to_sequencing_mappings()

# Find sequencing matches using stage coordinates
for file in files_red_laser:
    file.find_sequences_using_stage_coordinates(channel=1, show=True, save=False)

exp.show_sequencing_matches()

# Improve mapping by performing a linear least-squares fit on all nearest neighbours within the distance threshold
for file in files_red_laser:
    print(file)
    file.sequencing_match.nearest_neighbour_match(25, 'linear')
    file.sequencing_match.nearest_neighbour_match(10, 'linear')
    # file.sequencing_match.nearest_neighbour_match(10, 'polynomial', order=4)

file_red.sequencing_match.show_mapping_transformation(inverse=True, crop=True)


# -----------------------------------
# Classify traces in sequencing data of file
# -----------------------------------
# Doing this on all sequencing data is still slow.

sequences_dict = {
    'CalSeq': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCAACAATGCCTAGCCGATCCGTAATGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
    'MapSeq': 'ACTGACTGTAACAACAACAACAATAACAACAACAACAATAACAACAACAACAAT',
    'HJ1': 'NNNNNNNNNNNNNNNNNNNNCCCACCGCTCTTCTCAACTGGGTTTTCCCAGTTGAGAGCTTGCTAGGGTTTTCCCTAGCAAGCCGCTGCTACGGTTTTCCGTAGCAGCGAGAGCGGTGGG',
    'HJ3': 'NNNNNNNNNNNNNNNNNNNNCCCACCGCTCAACTCAACTGGGTTTTCCCAGTTGAGTCCTTGCTAGGGTTTTCCCTAGCAAGGGGCTGCTACGGTTTTCCGTAGCAGCCTGAGCGGTGGG',
    'HJ7': 'NNNNNNNNNNNNNNNNNNNNCCCACCGCTCGGCTCAACTGGGTTTTCCCAGTTGAGCGCTTGCTAGGGTTTTCCCTAGCAAGCCGCTGCTACGGTTTTCCGTAGCAGCGCGAGCGGTGGG',
    'HJ7_G116N': 'NNNNNNNNNNNNNNNNNNNNCCCACCGCTCGGCTCAACTGGGTTTTCCCAGTTGAGCNCTTGCTAGGGTTTTCCCTAGCAAGCCGCTGCTACGGTTTTCCGTAGCAGCGCGAGCGGTGGG',
    'HJ_general': 'NNNNNNNNNNNNNNNNNNNNCCCACCGCTCNNCTCAACTGGGTTTTCCCAGTTGAGNNCTTGCTAGGGTTTTCCCTAGCAAGNNGCTGCTACGGTTTTCCGTAGCAGCNNGAGCGGTGGG',
    'HJ1_variable_WT':      'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTNNNNNNNNNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNNNNNNNNGANNNNNNNNNN',
    'HJ3_variable_WT':      'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNTCNNNNNNNNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNNNNNNNNNNNNNNCTNNNNNNNNNN',
    'HJ7_variable_WT':      'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNNNNNNNNNNNNNNCGNNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNNNNNNNNGCNNNNNNNNNN',
    'HJ7_variable_G116C':   'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNNNNNNNNGCNNNNNNNNNN',
    'HJ7_variable_G116A':   'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNNNNNNNNNNNNNNCANNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNNNNNNNNGCNNNNNNNNNN',
    'HJ7_variable_G116T':   'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNNNNNNNNNNNNNNCTNNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNNNNNNNNGCNNNNNNNNNN',
    'HJ1_indexL1': 'GTCA' + 'N' * (120 - 4),
    'HJ3_indexL1': 'AAGG' + 'N' * (120 - 4),
    'HJ7_indexL1': 'TACC' + 'N' * (120 - 4),
    'HJ7_G116N_indexL1': 'AGTC' + 'N' * (120 - 4)
}

criteria_dict = {'CalSeq': 'CalSeq<5',
                 'MapSeq': 'MapSeq<5',
                 'HJ1_WT': '(HJ1_indexL1==0) & (HJ_general<10) & (HJ1_variable_WT==0)',
                 'HJ3_WT': '(HJ3_indexL1==0) & (HJ_general<10) & (HJ3_variable_WT==0)',
                 'HJ7_WT': '(HJ7_indexL1==0) & (HJ_general<10) & (HJ7_variable_WT==0)',
                 'HJ7_G': '(HJ7_G116N_indexL1==0) & (HJ_general<10) & (HJ7_variable_WT==0)',
                 'HJ7_G116C': '(HJ7_G116N_indexL1==0) & (HJ_general<10) & (HJ7_variable_G116C==0)',
                 'HJ7_G116A': '(HJ7_G116N_indexL1==0) & (HJ_general<10) & (HJ7_variable_G116A==0)',
                 'HJ7_G116T': '(HJ7_G116N_indexL1==0) & (HJ_general<10) & (HJ7_variable_G116T==0)',
                 }

for file_r in files_red_laser:
    file_r.get_sequencing_data_for_file() # Only necessary since classification on exp.sequencing_data is slow
    file_r.sequencing_data.classify(sequences_dict, criteria_dict)


for file_g, file_r in zip(files_green_laser, files_red_laser):
    file_g.sequencing_data = file_r.sequencing_data
    file_g.sequencing_match = file_r.sequencing_match

# -----------------------------------
# Find new coordinates green file
# -----------------------------------
configuration = exp.configuration['find_coordinates'].copy()
configuration['peak_finding']['minimum_intensity_difference'] = 300
configuration['peak_finding']['maximum_intensity_difference'] = 100000
configuration['peak_finding']['filter_neighbourhood_size'] = 5
configuration['channels'] = ['donor','acceptor']
configuration['method'] = 'average_channels'
for file_g in files_green_laser:
    file_g.find_coordinates(configuration=configuration) # Using new configuration that finds all visible spots, not only the mapping sequence

file_green.show_coordinates_in_image()

# -----------------------------------
# Optimize sequencing_match using all visible coordinates
# -----------------------------------
seq_names = ['MapSeq','HJ']
for file_g in files_green_laser:
    file_g.optimize_sequencing_match_using_visible_sequences(visible_sequence_names=seq_names, distance_threshold=15,
                                                             transformation_type='polynomial', order=4)
file_green.plot_sequencing_match()
file_green.sequencing_match.show_mapping_transformation(inverse=True)

# -----------------------------------
# Determine to what sequence each coordinate corresponds
# -----------------------------------
for file_g in files_green_laser:
    file_g.get_sequencing_data_for_file()
    file_g.sequencing_data.classify(sequences_dict, criteria_dict)
    file_g.determine_sequences_at_current_coordinates(visible_sequence_names=seq_names, distance_threshold=2.5)
    file_g.export_pks_file()

# -----------------------------------
# Extract traces
# -----------------------------------
for file_g in files_green_laser:
    file_g.extract_traces()

for molecule in file_green.molecules:
    molecule.isSelected = False
    if (molecule.sequence_index is not None) and (molecule.sequence_name is not 'MapSeq'):
        molecule.isSelected = True
