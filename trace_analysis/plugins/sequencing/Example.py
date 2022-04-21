import numpy as np
from pathlib import Path # For efficient path manipulation
import matplotlib.pyplot as plt
from git import Repo
import time
import xarray as xr
from skimage.transform import AffineTransform

import trace_analysis as ta
from trace_analysis.plugins.sequencing.fastqAnalysis import FastqData
from trace_analysis.mapping.mapping import Mapping2

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



#####################################
# SEQUENCING DATA PROCESSING
#####################################

# TODO: Make a python function to convert a complete bam/sam file to nc dataset with written out sequences

# -----------------------------------
# Open existing dataset
# -----------------------------------
# sequencing_dataset_path = r'G:\Ivo\20211005 - Objective-type TIRF (BN)\Analysis\HJ_general_mapped.nc'
# sequencing_dataset = xr.open_dataset(sequencing_dataset_path)
# sequencing_dataset = sequencing_dataset.set_index(sequence=['tile', 'x', 'y'])
#
# sequencing_dataset['read_aligned'] = sequencing_dataset.read_aligned.astype(str)
# sequencing_dataset['quality_aligned'] = sequencing_dataset.quality_aligned.astype(str)

#####################################
# SINGLE-MOLECULE DATA PROCESSING
#####################################

# experiment_path = r'O:\Ivo\20211005 - Objective-type TIRF (BN)'
experiment_path = r'H:\Desktop\20211105 - Test'
# experiment_path = r'C:\Users\Ivo Severins\Desktop\20211005 - Test'
# experiment_path = r'C:\Users\severins\Desktop\20211005 - Test'
experiment_path = r'N:\tnw\BN\CMJ\Shared\Ivo\PhD_data\20211005 - Objective-type TIRF (BN)'
exp = ta.Experiment(experiment_path)

files_channel_mapping = exp.files[exp.files.name.regex('Mapping')]
files_green_laser = exp.files[exp.files.name.regex('TIRF 561')]
files_red_laser = exp.files[exp.files.name.regex('TIRF 642')]

# -----------------------------------
# Channel mapping
# -----------------------------------
channel_mapping_file = files_channel_mapping[0]
channel_mapping_file.perform_mapping()

channel_mapping_file.show_average_image()
channel_mapping_file.mapping.show_mapping_transformation(figure=plt.gcf(), show_source=True)

# -----------------------------------
# Find coordinates for sequence mapping
# -----------------------------------
exp.import_config_file()
configuration = exp.configuration['find_coordinates'].copy()
# configuration['peak_finding']['minimum_intensity_difference'] = 4000
configuration['peak_finding']['minimum_times_background'] = 11 # First try was 7
configuration['channels'] = ['acceptor']
configuration['method'] = 'by_channel'

files_red_laser.find_coordinates(configuration=configuration)

# -----------------------------------
# Find coordinates, extract traces and determine kinetics
# -----------------------------------

exp.import_config_file()
configuration = exp.configuration['find_coordinates'].copy()
# configuration['peak_finding']['minimum_intensity_difference'] = 4000
configuration['peak_finding']['minimum_times_background'] = 1.2
configuration['channels'] = ['donor', 'acceptor']
configuration['method'] = 'sum_channels'

files_green_laser.find_coordinates(configuration=configuration)
files_green_laser.show_coordinates_in_image()

files_green_laser.extract_traces2()


#####################################
# SEQUENCING AND SINGLE-MOLECULE MAPPING
#####################################

mapping_sequence_name = 'MapSeq'
# TODO: Select sequence from sequencing dataset

filepath_sequencing_data_for_mapping = r'G:\Ivo\20211011 - Sequencer (MiSeq)\Analysis\sequencing_data_MapSeq.csv'
exp.import_sequencing_data_for_mapping(filepath_sequencing_data_for_mapping, surface=0)
exp.generate_tile_mappings(files_red_laser)

# -----------------------------------
# Finding rotation and scale with geometric hashing
# -----------------------------------

# TODO: Geometric hashing example here
# TODO: Put geometric hashing in Mapping2



# -----------------------------------
# Finding translation with cross correlation
# -----------------------------------

# Previously found transformation
# Transformation based on sm pixel to sequencing MiSeq mapping with
# 'rotation': 0.006500218506032994, 'scale': [3.697153992993506, -3.697153992993506] using pixel size 0.125 µm
# exp.tile_mappings.transformation = AffineTransform(scale=[29.57723194394805, -29.57723194394805], rotation=0.006500218506032994)
exp.tile_mappings.transformation = AffineTransform(scale=[29.51, -29.51], rotation=0.003)


exp.tile_mappings.cross_correlation(peak_detection='auto', gaussian_width=7, divider=20, plot=False)

# bounds = ((0.98, 1.2), (-0.005, 0.005), (-250, 250), (-250, 250))
# exp.tile_mappings.kernel_correlation(bounds, sigma=25, crop='source',
#                                      strategy='best1bin', maxiter=1000, popsize=50, tol=0.01,
#                                      mutation=0.25, recombination=0.7, seed=None, callback=None, disp=False,
#                                      polish=False, init='sobol', atol=0, updating='immediate', workers=1,
#                                      constraints=())

exp.tile_mappings.save()
exp.tile_mappings.show_mapping_transformation(crop='source', save=True)

# -----------------------------------
# Determine translations for all tiles
# -----------------------------------
exp.tile_mappings.scatter_parameters('translation', 'translation', 'x', 'y', save=True)  # To determine correct mapping indices
exp.tile_mappings.estimate_translations(indices=np.arange(11), save=True)
exp.tile_mappings.scatter_parameters('translation', 'translation', 'x', 'y', save=True) # Rename original file before saving
exp.tile_mappings.save()
exp.tile_mappings.show_mapping_transformation(crop='source', save=True)

# -----------------------------------
# Matching statistics
# -----------------------------------

# TODO: Show how to extract matching statics, i.e. numbers/percentages of points matched.

# -----------------------------------
# Import sequencing data
# -----------------------------------

filepath_sequencing_data = r'G:\Ivo\20211011 - Sequencer (MiSeq)\Analysis\sequencing_data_HJ_general.csv'
exp.import_sequencing_data(filepath_sequencing_data, surface=0)

# -----------------------------------
# Obtain file sequencing data and generate sequencing match
# -----------------------------------
files = files_green_laser[0:14*30]

# TODO: Make this automatically save the sequencing data
# TODO: Automatic import and export when setting or getting sequencing data
files_green_laser.get_sequencing_data(margin=5)
files_green_laser.generate_sequencing_match(overlapping_points_threshold=25)

# -----------------------------------
# Finetune the sequencing matches
# -----------------------------------

sequencing_matches = exp.sequencing_matches(files_green_laser)
#3087 & 3088 & 3434
sequencing_matches.cross_correlation(divider=1/10, gaussian_width=7, crop=True, plot=False)


sequencing_matches.transformation = AffineTransform()
sequencing_matches.transformation_inverse = AffineTransform()
bounds = ((0.97, 1.03), (-0.05, 0.05), (-5, 5), (-5, 5))
sequencing_matches.kernel_correlation(bounds, sigma=0.125, crop=True,
                                         strategy='best1bin', maxiter=1000, popsize=50, tol=0.01,
                                         mutation=0.25, recombination=0.7, seed=None, callback=None, disp=False,
                                         polish=True, init='sobol', atol=0, updating='immediate', workers=1,
                                         constraints=())

bounds = ((0.99, 1.01), (-0.01, 0.01), (-1, 1), (-1, 1))
sequencing_matches.kernel_correlation(bounds, sigma=0.125, crop=True,
                                         strategy='best1bin', maxiter=1000, popsize=50, tol=0.001,
                                         mutation=0.25, recombination=0.7, seed=None, callback=None, disp=False,
                                         polish=True, init='sobol', atol=0, updating='immediate', workers=1,
                                         constraints=())

# -----------------------------------
# Find pairs and insert sequencing data into file dataset
# -----------------------------------

sequencing_matches.find_distance_threshold(maximum_radius=2)
sequencing_matches.determine_matched_pairs()
plt.figure()
plt.hist(np.hstack(sequencing_matches.pair_distances()), bins=100)

sequencing_matches.destination_distance_threshold = 0.2  # 0.506
sequencing_matches.determine_matched_pairs()
sequencing_matches.save()

sequencing_matches.show_mapping_transformation()

files_green_laser.insert_sequencing_data_into_file_dataset()








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
