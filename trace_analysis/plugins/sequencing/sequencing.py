import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from skimage.transform import AffineTransform
# from trace_analysis.experiment import Experiment
# from trace_analysis.file import File
from trace_analysis.mapping.geometricHashing import SequencingDataMapping
from trace_analysis.mapping.mapping import Mapping2
from .fastqAnalysis import FastqData
from .geometricHashing2 import geometric_hash, find_match_after_hashing
from .plotting import plot_sequencing_match
from trace_analysis.mapping.icp import icp, nearest_neighbor_pair


class Experiment:
    def import_sequencing_data(self, fastq_file_path):
        self.fastq_file_path = Path(fastq_file_path)
        self.sequencing_data = FastqData(self.fastq_file_path)

    def generate_mapping_hashtable(self, mapping_sequence, number_of_allowed_mismatches,
                                   imaged_surface=None, maximum_distance_tile=None, tuple_size=None):

        self.select_sequencing_data_for_mapping(mapping_sequence, number_of_allowed_mismatches)

        if imaged_surface in ['top', 1]:
            self.sequencing_data_for_mapping = self.sequencing_data_for_mapping[self.sequencing_data_for_mapping.tile < 2000]
        elif imaged_surface in ['bottom', 2]:
            self.sequencing_data_for_mapping = self.sequencing_data_for_mapping[self.sequencing_data_for_mapping.tile > 2000]

        tile_coordinate_sets = [tile.coordinates for tile in self.sequencing_data_for_mapping.tiles]
        # TODO: get maximum_distance_tile and tuple_size from configuration
        self.geometric_hash_data = geometric_hash(tile_coordinate_sets, maximum_distance_tile, tuple_size)

    def generate_mapping_hashtable_from_coordinate_set(self, tile_coordinate_sets, maximum_distance_tile, tuple_size):
        self.geometric_hash_data = geometric_hash(tile_coordinate_sets, maximum_distance_tile, tuple_size)

    def select_sequencing_data_for_mapping(self, mapping_sequence, number_of_allowed_mismatches):
        self.mapping_sequence = mapping_sequence
        self.sequencing_data.matches_per_tile(sequence=mapping_sequence)

        number_of_matches = self.sequencing_data.number_of_matches(mapping_sequence)
        # sequencing_data_for_mapping = sequencing_data.select(Nmatch == len(mapping_sequence), copyData=True)
        self.sequencing_data_for_mapping = self.sequencing_data[number_of_matches >=
                                                                (len(mapping_sequence) - number_of_allowed_mismatches)]
        self.sequencing_data_for_mapping.show_tiles()
        self.sequencing_data_for_mapping.export_positions_per_tile()

    # def map_files_to_sequencing_data(self, tile, configuration=None):
    #     sequencing_mapping_path = self.mainPath.joinpath('sequencing_mapping')
    #     sequencing_mapping_path.mkdir(exist_ok=True)
    #
    #     presets = {'TIR-I single-molecule':
    #                    {
    #                        'initial_image_transformation': {'reflection': 0, 'rotation': np.pi, 'magnification': 3.336},
    #                        'mapping_configuration': {'mode': 'similarity',
    #                                                 'nBins': 250,
    #                                                 'rotationRange': np.array([-180, 180]) / 360 * 2 * np.pi,
    #                                                 'magnificationRange': np.array([3000, 6000])
    #                                                 },
    #                        'bases_findMatch': 500
    #                     }
    #                }
    #     # rotationRange = np.array([-5,5])/360*2*np.pi+np.pi/2,
    #     # magnificationRange = np.array([2750,3000])
    #     # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': -35/180*np.pi, 'magnification':1/125 }
    #     # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': np.pi, 'magnification': 2*3.336}
    #     # seqmap.initial_image_transformation = {}
    #
    #     if type(configuration) is str:
    #         configuration = presets[configuration]
    #
    #     tile = self.sequencing_data_for_mapping.get_tile_object(tile=tile)
    #     self.seqmap = SequencingDataMapping(tile, self.files, sequencing_mapping_path, **configuration['mapping_configuration'])
    #     self.seqmap.initial_image_transformation = configuration['initial_image_transformation']
    #     self.seqmap.bases_findMatch = configuration['bases_findMatch']
    #     self.seqmap.histogram_matches(export=True)

    def sequencing_mapping_to_files(self, minimal_number_of_matching_points):
        self.seqmap.give_matches_to_files(match_threshold=minimal_number_of_matching_points)

    # TODO: Put this improvement in file
    def sequencing_mapping_improvement(self):
        for match in self.seqmap.matches:
            match.nearest_neighbour_match(distance_threshold=25)


# def map_sequences_to_molecules(files, sequencing_data_for_mapping, mapping_sequence, tile, write_path, match_threshold = 5):
#
#

# tile = sequencing_data_for_mapping.get_tile_object(tile=tile)
#
# seqmap = SequencingDataMapping(tile, files, 'similarity',
#                                write_path, nBins=250,
#                                rotationRange=np.array([-180, 180]) / 360 * 2 * np.pi,
#                                magnificationRange=np.array([3000, 6000]))
# #                               rotationRange = np.array([-5,5])/360*2*np.pi+np.pi/2,
# #                               magnificationRange = np.array([2750,3000]))
# # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': -35/180*np.pi, 'magnification':1/125 }
# # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': np.pi, 'magnification': 2*3.336}
# seqmap.initial_image_transformation = {'reflection': 0, 'rotation': np.pi, 'magnification': 3.336}
# # seqmap.initial_image_transformation = {}
# seqmap.bases_findMatch = 500
# seqmap.histogram_matches(export=True)
# seqmap.give_matches_to_files()

# matches = [match for match in seqmap.matches if match.count >= match_threshold]
# # matches = [match for match in seqmap.matches if match.percentMatch > 0.69]
# for match in matches:
#     match.nearest_neighbour_match(distance_threshold=25)
#     seqmap.plot_match(match)

# return seqmap, matches

# Probably it would be better to move the function below somewhere else [IS 12-02-2020]
def within_bounds(coordinates, bounds, margin=0):
    bounds = np.sort(bounds)
    criteria = np.array([(coordinates[:, 0] > (bounds[0, 0] + margin)),
                         (coordinates[:, 0] < (bounds[1, 0] - margin)),
                         (coordinates[:, 1] > (bounds[0, 1] + margin)),
                         (coordinates[:, 1] < (bounds[1, 1] - margin))
                         ])

    return criteria.all(axis=0)


class File:
    # def map_file_sequences_to_molecules(self, mapping_sequence, tile, match_threshold=5):
    #     sequencing_mapping_path = self.experiment.mainPath.joinpath('sequencing_mapping')
    #     sequencing_mapping_path.mkdir(exist_ok=True)
    #     # Probably we will not have to export seqmap, at least not to file.
    #     self.seqmap, self.matches = map_sequences_to_molecules([self], self.experiment.sequencing_data, mapping_sequence,
    #                                                             tile, sequencing_mapping_path, match_threshold)
    #

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.sequencing_data = None
        self.sequencing_match = None

        self.importFunctions['.fastq'] = self.import_sequencing_data

        if self.experiment.import_all is True:
            self.findAndAddExtensions()

    @property
    def sequences(self):
        if len(self.molecules) > 0:
            return np.vstack([molecule.sequence for molecule in self.molecules])
        else:
            return np.array([])

    @sequences.setter
    def sequences(self, sequences):
        for i, molecule in enumerate(self.molecules):
            molecule.sequence = sequences[i]

    @property
    def sequencing_tile(self):
        if self.sequencing_data:
            if len(self.sequencing_data.tiles)==1:
                return self.sequencing_data.tiles[0]
            else:
                raise ValueError('File contains sequences from multiple tiles')
        else:
            return None

    # @property
    # def sequencing_data(self):
    #     return self._sequencing_data
    #
    # @sequencing_data.setter
    # def sequencing_data(self, sequencing_data):
    #     # self.molecules = []
    #     # self.number_of_molecules = len(sequencing_data)
    #     # for i, molecule in enumerate(self.molecules):
    #     #     molecule.sequencing_data = sequencing_data[i]

    def import_sequencing_data(self):
        self.sequencing_data = FastqData(self.absoluteFilePath.with_suffix('.fastq'))
        self.sequences = self.sequencing_data.sequence

    def find_sequences(self, maximum_distance_file, tuple_size, initial_transformation={},
                       hash_table_distance_threshold=0.01,
                       alpha=0.1, test_radius=10, K_threshold=10e9,
                       magnification_range=None, rotation_range=None,
                       nearest_neighbour_match_distance_threshold=25):
        # TODO: Make geometric hashing reflection invariant
        initial_transformation = AffineTransform(**initial_transformation)

        # TODO: make the following line more general and remove bounds dependence in geometric hashing
        coordinate_vertices_file = initial_transformation(self.movie.channel_vertices('a'))

        #self.geometric_hash_data = geometric_hash(initial_transform(self.coordinates), maximum_distance_file, tuple_size)

        #match.destination_index = destination_index
        match = find_match_after_hashing(initial_transformation(self.coordinates), maximum_distance_file, tuple_size, coordinate_vertices_file,
                                         *self.experiment.geometric_hash_data,
                                         hash_table_distance_threshold, alpha, test_radius, K_threshold,
                                         magnification_range, rotation_range)
        if match:
            # self.sequencing_tile = self.experiment.sequencing_data_for_mapping.tiles[match.destination_index]
            match.source = self.coordinates
            match.initial_transformation = initial_transformation
            match.transformation = match.transformation @ initial_transformation.params
            # TODO: Base this on some better criteria
            #match.nearest_neighbour_match(nearest_neighbour_match_distance_threshold)
            self.sequencing_match = match
            #self.get_all_sequences_from_sequencing_data()

    def get_all_sequences_from_sequencing_data(self):
        # raise Warning('Only works on acceptor channel for now')
        coordinate_bounds_file = self.movie.channel_boundaries('a')
        coordinate_bounds_tile = self.sequencing_match.transform_coordinates(coordinate_bounds_file)

        tile = self.experiment.sequencing_data_for_mapping.tiles[self.sequencing_match.destination_index]

        self.sequencing_data = \
            self.experiment.sequencing_data.get_selection(tile=tile.number,
                                                          x=coordinate_bounds_tile[:, 0],
                                                          y=coordinate_bounds_tile[:, 1])

        self.molecules = []
        self.coordinates = self.sequencing_match.transform_coordinates(self.sequencing_data.coordinates,
                                                                     inverse=True)
        self.sequencing_data.export_fastq(self.relativeFilePath)

    def plot_sequencing_match(self):
        #plot_sequencing_match(self.sequencing_match)
        #name = f'Tile: {self.sequencing_tile.name}, File: {str(self.relativeFilePath)}'
        name = self.name + '_sequencing_mapping'
        print(name)
        plot_sequencing_match(self.sequencing_match, self.absoluteFilePath.parent, name, unit='um')

    # This is probably not the way to go
    # def give_molecules_closest_sequence(self):
    #     sequencing_data = self.experiment.sequencing_data
    #     indices_within_tile = sequencing_data.selection(tile=int(self.sequence_match.tile.name))
    #     sequencing_data = sequencing_data.select(indices_within_tile, copyData=True)
    #
    #     #raise Warning("Not implemented for the donor channel yet")
    #     boundaries = self.sequence_match.transform_coordinates(self.movie.channel_boundaries('a'))
    #
    #     indices_within_bounds = within_bounds(sequencing_data.coordinates, boundaries)
    #     sequencing_data.select(indices_within_bounds, copyData=False)
    #
    #     tile_coordinates = sequencing_data.coordinates
    #     image_coordinates = self.sequence_match.transform_coordinates(self.coordinates)
    #
    #     from trace_analysis.plotting import scatter_coordinates
    #     plt.figure()
    #     scatter_coordinates([tile_coordinates, image_coordinates])
    #
    #
    #     distances, image_indices, tile_indices = nearest_neighbor_pair(image_coordinates, tile_coordinates)
    #
    #     distance_threshold = 25
    #     image_indices = image_indices[distances < distance_threshold]
    #     tile_indices = tile_indices[distances < distance_threshold]
    #
    #     print(*[sequence.tostring() for sequence in sequencing_data.sequence[tile_indices]], sep="\n")
    #
    #     plt.figure()
    #     scatter_coordinates([tile_coordinates, self.sequence_match.destination, image_coordinates])
    #     from trace_analysis.plotting import show_point_connections
    #     show_point_connections(image_coordinates[image_indices], tile_coordinates[tile_indices])
    #
    #     for molecule in self.molecules: molecule.sequence = np.array([], dtype=bytes)
    #     for i in np.arange(len(image_indices)):
    #         self.molecules[image_indices[i]].sequence = sequencing_data.sequence[tile_indices[i]]
    #         self.molecules[image_indices[i]].sequencing_data = sequencing_data.select(i, copyData=True)


class Molecule:
    @property
    def sequence(self):
        return self.file.sequencing_data.sequence[self.index, :].tostring()

    @property
    def sequencing_data(self):
        return self.file.sequencing_data[self.index]
