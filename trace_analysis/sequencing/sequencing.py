import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from trace_analysis.mapping.geometricHashing import SequencingDataMapping
from trace_analysis.sequencing.fastqAnalysis import FastqData
from trace_analysis.mapping.icp import icp, nearest_neighbor_pair

def import_sequences(experiment, fastq_file):
    experiment.fastq_file = Path(fastq_file)
    experiment.sequencing_data = FastqData(experiment.fastq_file)

def map_sequences_to_molecules(files, sequencing_data, mapping_sequence, tile, write_path, match_threshold = 5):
    sequencing_data.matches_per_tile(sequence=mapping_sequence)
    Nmatch = sequencing_data.number_of_matches(mapping_sequence)
    sequencing_data_for_mapping = sequencing_data.select(Nmatch == len(mapping_sequence), copyData=True)
    sequencing_data_for_mapping.show_tiles()
    sequencing_data_for_mapping.export_positions_per_tile()

    tile = sequencing_data_for_mapping.get_tile_object(tile=tile)

    seqmap = SequencingDataMapping(tile, files, 'similarity',
                                   write_path, nBins=250,
                                   rotationRange=np.array([-180, 180]) / 360 * 2 * np.pi,
                                   magnificationRange=np.array([3000, 6000]))
    #                               rotationRange = np.array([-5,5])/360*2*np.pi+np.pi/2,
    #                               magnificationRange = np.array([2750,3000]))
    # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': -35/180*np.pi, 'magnification':1/125 }
    # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': np.pi, 'magnification': 2*3.336}
    seqmap.initial_image_transformation = {'reflection': 0, 'rotation': np.pi, 'magnification': 3.336}
    # seqmap.initial_image_transformation = {}
    seqmap.bases_findMatch = 500
    seqmap.histogram_matches(export=True)
    seqmap.give_matches_to_files()

    matches = [match for match in seqmap.matches if match.count >= match_threshold]
    # matches = [match for match in seqmap.matches if match.percentMatch > 0.69]
    for match in matches:
        match.nearest_neighbour_match(distance_threshold=25)
        seqmap.plot_match(match)

    return seqmap, matches

# Probably it would be better to move the function below somewhere else [IS 12-02-2020]
def within_bounds(coordinates, bounds, margin=0):
    bounds = np.sort(bounds)
    criteria = np.array([(coordinates[:, 0] > (bounds[0, 0] + margin)),
                         (coordinates[:, 0] < (bounds[0, 1] - margin)),
                         (coordinates[:, 1] > (bounds[1, 0] + margin)),
                         (coordinates[:, 1] < (bounds[1, 1] - margin))
                         ])

    return criteria.all(axis=0)

def give_molecules_closest_sequence(files):
    for file in files:
        sequencing_data = file.experiment.sequencing_data
        indices_within_tile = sequencing_data.selection(tile=int(file.sequence_match.tile.name))
        sequencing_data = sequencing_data.select(indices_within_tile, copyData=True)

        #raise Warning("Not implemented for the donor channel yet")
        boundaries = file.sequence_match.transform_coordinates(file.movie.channel_boundaries('a').T).T

        indices_within_bounds = within_bounds(sequencing_data.coordinates, boundaries)
        sequencing_data.select(indices_within_bounds, copyData=False)

        tile_coordinates = sequencing_data.coordinates
        image_coordinates = file.sequence_match.transform_coordinates(file.coordinates)

        from trace_analysis.plotting import scatter_coordinates
        plt.figure()
        scatter_coordinates([tile_coordinates, image_coordinates])


        distances, image_indices, tile_indices = nearest_neighbor_pair(image_coordinates, tile_coordinates)

        distance_threshold = 25
        image_indices = image_indices[distances < distance_threshold]
        tile_indices = tile_indices[distances < distance_threshold]

        print(*[sequence.tostring() for sequence in sequencing_data.sequence[tile_indices]], sep="\n")

        plt.figure()
        scatter_coordinates([tile_coordinates, file.sequence_match.destination, image_coordinates])
        from trace_analysis.plotting import show_point_connections
        show_point_connections(image_coordinates[image_indices], tile_coordinates[tile_indices])

        for molecule in file.molecules: molecule.sequence = np.array([], dtype=bytes)
        for i in np.arange(len(image_indices)):
            file.molecules[image_indices[i]].sequence = sequencing_data.sequence[tile_indices[i]]
            file.molecules[image_indices[i]].sequencing_data = sequencing_data.select(i, copyData=True)