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
writePath = Path(Path(__file__).parent)
writeFile = writePath.joinpath('git_repository_version_and_differences.txt')
repo = Repo(Path(ta.__file__).parent.parent)

with writeFile.open("a") as f:
    f.write('------------------------------------------\n\n')
    f.write(f'Trace_analysis version: {repo.head.object.hexsha} \n\n')
    t = repo.head.commit.tree
    f.write(repo.git.diff(t))
    f.write('\n\n------------------------------------------\n\n')


# =================
# FLUORESCENCE DATA
# =================

# Import experiment
# IMPORTANT: uses only a single colour
exp = ta.Experiment(r'Path_to_fluorescence_data', colours=['r'])

# Find coordinates for the files that have no molecules yet
# (remove the if statement if you want all files to be reanalyzed)
for file in exp.files:
    if len(file.molecules)==0:
        file.find_coordinates()

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

# Select only the sequences that are used for mapping
mapping_sequence = 'TATCTGTATAATGAGAAATATGGAGTACAATTTTTTTTTTTTTTTTTTTT'
number_of_allowed_mismatches = 0
exp.select_sequencing_data_for_mapping(mapping_sequence, number_of_allowed_mismatches)

# Set the configuration for the geometric hashing algorithm
# The magnification would be of the order imagewidth in pixels * magnification

# ============================
# MAPPING BY GEOMETRIC HASHING
# ============================

# Translation, rotation and magnification/scaling
configuration = {
                    'initial_image_transformation': {'reflection': 0, 'rotation': np.pi, 'magnification': 3.336},
                    'mapping_configuration': {'mode': 'similarity',
                                              'nBins': 250,
                                              'rotationRange': np.array([-180, 180]) / 360 * 2 * np.pi,
                                              'magnificationRange': np.array([3000, 6000])
                                              },
                    'bases_findMatch': 500
                }

# Mapping translation only
configuration = {
                    'initial_image_transformation': {'reflection': 0, 'rotation': np.pi, 'magnification': 3.336},
                    'mapping_configuration': {'mode': 'translation',
                                              'nBins': 500,
                                              },
                    'bases_findMatch': 500
                }

# Map all files to a specific tile
# Prism-type TIRF: 2... tiles
# Objective-type TIRF: 1... tiles
tile = 2102
exp.map_files_to_sequencing_data(tile, configuration=configuration)

# In case you want to remove the sequencing_match from the files, which is added again by the next line.
# for file in exp.files:
#     if hasattr(file,'sequence_match'):
#         delattr(file, 'sequence_match')

# Transfer the matches between files and tile that were found with mapping to the File objects
# Only if the match point count is larger or equal to a threshold
exp.sequencing_mapping_to_files(minimal_number_of_matching_points=4)

# Improve the mapping by performing an additional linear least squares fit (This may not work yet for all molecules)
# (If non-linear mapping is ready, we can also use that as it likely provides a better fit)
exp.sequencing_mapping_improvement()

# Make and export a sequencing match plot for files that have a sequencing match
for file in exp.files:
    if hasattr(file,'sequence_match'):
        file.plot_sequencing_match()


# Calculate rotation and magnification for all found matches
rotation = np.array([match.rotation for match in exp.seqmap.matches])
magnification = np.array([match.magnification[0] for match in exp.seqmap.matches])

# Make a heatmap scatter plot
from scipy.stats import gaussian_kde
x = rotation
y = magnification

xy = np.vstack([x,y])
z = gaussian_kde(xy,0.01)(xy)
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=100, edgecolor='')
plt.title('Matches')
plt.xlabel('Rotation (deg)')
plt.ylabel('Magnfication')
plt.show()


# Select matches that form the peak in heat map scatter plot
selection = np.all(np.array([rotation<-179,rotation>-181,magnification>3.25,magnification<3.4]),axis=0)
matches = [match for i, match in enumerate(exp.seqmap.matches) if selection[i]]

selected_rotations = [match.rotation for match in matches]
selected_magnifications = [match.magnification[0] for match in matches]

# Calculate mean rotation and magnification and their standard deviations
mean_rotation = np.mean(selected_rotations)
std_rotation = np.std(selected_rotations)
mean_magnification = np.mean(selected_magnifications)
std_magnification = np.std(selected_magnifications)
print(mean_rotation,std_rotation)
print(mean_magnification,std_magnification)

# Show that the selection only contains correct matches
for file in exp.files:
    if hasattr(file,'sequence_match'):
        if (file.sequence_match.rotation <-179) & (file.sequence_match.rotation > -181) & \
            (file.sequence_match.magnification[0]>3.25) & (file.sequence_match.magnification[0]<3.4):
            #file.plot_sequencing_match()
            print(file)

# Put the matched sequencing data into each File and give each molecule its sequence
# This also generates a fastq file for each file
for file in exp.files:
    if hasattr(file,'sequence_match'):
        file.get_all_sequences_from_sequencing_data()
        # Plotting of the coordinates on top of the average image
        # fig = plt.figure()
        # file.show_coordinates(figure=fig)
        # file.show_average_image(figure=fig)

