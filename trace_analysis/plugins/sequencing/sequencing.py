import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as pth
import math
from pathlib import Path
from skimage.transform import AffineTransform
import pandas as pd
import xarray as xr
# from trace_analysis.experiment import Experiment
# from trace_analysis.file import File
from trace_analysis.mapping.geometricHashing import SequencingDataMapping
from trace_analysis.mapping.mapping import Mapping2, crop_coordinates_indices
from trace_analysis.mapping.polynomial import PolynomialTransform
from .fastqAnalysis import FastqData
from .geometricHashing2 import geometric_hash, find_match_after_hashing
from .geometricHashing3 import GeometricHashTable
from .plotting import plot_sequencing_match, plot_matched_files_in_tile
from trace_analysis.mapping.icp import icp, nearest_neighbor_pair


class Experiment:
    def import_sequencing_data(self, fastq_file_path):
        self.fastq_file_path = Path(fastq_file_path)
        self.sequencing_data = FastqData(self.fastq_file_path)

    def generate_mapping_hashtable(self, mapping_sequence, number_of_allowed_mismatches,
                                   imaged_surface=None, maximum_distance_tile=None, tuple_size=None):

        # TODO: Add timer to generate_mapping_hashtable and find_sequences methods, by making a decorator function. [IS: 10-08-2020]

        self.select_sequencing_data_for_mapping(mapping_sequence, number_of_allowed_mismatches)

        if imaged_surface in ['top', 1]:
            self.sequencing_data_for_mapping = self.sequencing_data_for_mapping[self.sequencing_data_for_mapping.tile < 2000]
        elif imaged_surface in ['bottom', 2]:
            self.sequencing_data_for_mapping = self.sequencing_data_for_mapping[self.sequencing_data_for_mapping.tile > 2000]

        tile_coordinate_sets = [tile.coordinates for tile in self.sequencing_data_for_mapping.tiles]
        # TODO: get maximum_distance_tile and tuple_size from configuration
        self.geometric_hash_data = geometric_hash(tile_coordinate_sets, maximum_distance_tile, tuple_size)

    def generate_mapping_hashtable3(self, mapping_sequence, number_of_allowed_mismatches,
                                    imaged_surface=None, initial_file_transformation=None, maximum_distance_tile=None, tuple_size=None):

        # TODO: Add timer to generate_mapping_hashtable and find_sequences methods, by making a decorator function. [IS: 10-08-2020]

        self.select_sequencing_data_for_mapping(mapping_sequence, number_of_allowed_mismatches)

        if imaged_surface in ['top', 1]:
            self.sequencing_data_for_mapping = self.sequencing_data_for_mapping[self.sequencing_data_for_mapping.tile < 2000]
        elif imaged_surface in ['bottom', 2]:
            self.sequencing_data_for_mapping = self.sequencing_data_for_mapping[self.sequencing_data_for_mapping.tile > 2000]

        tile_coordinate_sets = [tile.coordinates for tile in self.sequencing_data_for_mapping.tiles]

        # initial_magnification = np.array([3.67058194, -3.67058194])
        # initial_rotation = 0.6285672733195177  # degrees
        # initial_file_transformation = AffineTransform(matrix=None, scale=initial_magnification,
        #                                                 rotation=initial_rotation/360*np.pi*2,
        #                                                 shear=None, translation=None)
        self.geometric_hashtable = GeometricHashTable(tile_coordinate_sets,
                                                      source_vertices=self.files[0].movie.channels[1].vertices,
                                                      initial_source_transformation=initial_file_transformation)

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
    #     sequencing_mapping_path = self.main_path.joinpath('sequencing_mapping')
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

    # def sequencing_mapping_to_files(self, minimal_number_of_matching_points):
    #     self.seqmap.give_matches_to_files(match_threshold=minimal_number_of_matching_points)

    def map_sequencing_and_stage_coordinates(self):
        tiles = self.sequencing_data.tiles
        self.stage_to_sequencing_mappings = []

        for tile in tiles:
            matched_files_in_tile = [file for file in self.files if file.sequencing_match and
                                     file.sequencing_match.tile == tile.number]

            if len(matched_files_in_tile) == 0:
                self.stage_to_sequencing_mappings.append(None)
                continue

            stage_coordinates = np.array([file.movie.stage_coordinates[0] for file in matched_files_in_tile])
            # stage_coordinates_in_pixels = np.array([file.movie.stage_coordinates_in_pixels[0] for file in matched_files_in_tile])
            sequencing_coordinates = np.vstack(
                [file.sequencing_match.transform_coordinates(np.array([[0, 0]])) for file in matched_files_in_tile])

            stage_to_tile_mapping = \
                Mapping2(stage_coordinates, sequencing_coordinates, transformation_type='linear',
                         source_name='Stage coordinates', destination_name='Sequencing_coordinates',
                         source_unit='μm', destination_unit='FASTQ',
                         name=f'Stage to sequencing coordinates - tile {tile.number}')
            stage_to_tile_mapping.direct_match('linear')
            stage_to_tile_mapping.tile = tile.number

            stage_to_tile_mapping.save(self.main_path.joinpath(f'stage_to_tile_{tile.number}.mapping'))
            self.stage_to_sequencing_mappings.append(stage_to_tile_mapping)

    def show_stage_to_sequencing_mappings(self):
        for mapping in self.stage_to_sequencing_mappings:
            if mapping is not None:
                mapping.show_mapping_transformation(source_colour='forestgreen', destination_colour='k',
                                                    save_path=self.main_path)

    def tile_boundaries_in_stage_coordinates(self):
        # TODO: Use real tile boundaries here
        boundaries = []
        for mapping in self.stage_to_sequencing_mappings:
            if mapping:
                boundaries.append(mapping.transform_coordinates([[0, 0], [30000, 30000]], inverse=True).T)
            else:
                boundaries.append(np.zeros((2,2)))
        return np.sort(boundaries)

    @property
    def files_with_sequencing_match(self):
        return [file for file in self.files if file.sequencing_match is not None]

    @property
    def sequencing_matches(self):
        return [file.sequencing_match for file in self.files if file.sequencing_match is not None]

    def show_sequencing_matches(self, show_file_coordinates=False):
        plot_matched_files_in_tile(self.files_with_sequencing_match, show_file_coordinates=show_file_coordinates,
                                   save=True)

    def sequencing_match_info_per_file(self, distance_threshold=25):
        columns = pd.MultiIndex.from_product([['File coordinates', 'Sequencing coordinates'],
                                              ['Matched', 'Total', 'Fraction']])
        df = pd.DataFrame(columns=columns)

        for file in self.files_with_sequencing_match:
            df.loc[file.name, ('File coordinates', 'Total')] = file.sequencing_match.source_cropped.shape[0]
            df.loc[file.name, ('Sequencing coordinates', 'Total')] = file.sequencing_match.destination_cropped.shape[0]
            df.loc[file.name, (['File coordinates','Sequencing coordinates'], 'Matched')] = \
                file.sequencing_match.number_of_matched_points(distance_threshold)
            df.loc[file.name, (['File coordinates','Sequencing coordinates'], 'Fraction')] = \
                file.sequencing_match.fraction_of_points_matched(distance_threshold)

        return df

    def sequencing_match_info_mean(self, distance_threshold=25):
        df = self.sequencing_match_info_per_file(distance_threshold)
        df = pd.DataFrame([df.mean(axis=0), df.std(axis=0)], index=['Mean', 'Std']).T
        return df.apply(std_string, axis=1)

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
    #     sequencing_mapping_path = self.experiment.main_path.joinpath('sequencing_mapping')
    #     sequencing_mapping_path.mkdir(exist_ok=True)
    #     # Probably we will not have to export seqmap, at least not to file.
    #     self.seqmap, self.matches = map_sequences_to_molecules([self], self.experiment.sequencing_data, mapping_sequence,
    #                                                             tile, sequencing_mapping_path, match_threshold)
    #

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.sequencing_data = None
        self.sequencing_match = None
        self.sequencing_match_old = None

        self.importFunctions['.fastq'] = self.import_sequencing_data
        self.importFunctions['_sequencing_match.mapping'] = self.import_sequencing_match

        if self.experiment.import_all is True:
            self.findAndAddExtensions()

    @property
    def sequences(self):
        if len(self.molecules) > 0:
            return np.vstack([molecule.sequence for molecule in self.molecules])
        else:
            return np.array([])

    # @sequences.setter
    # def sequences(self, sequences):
    #     for i, molecule in enumerate(self.molecules):
    #         molecule.sequence = sequences[i]

    @property
    def sequence_indices(self):
        if len(self.molecules) > 0:
            return np.array([molecule.sequence_index for molecule in self.molecules])
        else:
            return np.array([])

    @sequence_indices.setter
    def sequence_indices(self, sequence_indices):
        for i, molecule in enumerate(self.molecules):
            molecule.sequence_index = sequence_indices[i]

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
        # self.sequences = self.sequencing_data.sequence

    def import_sequencing_match(self):
        self.sequencing_match = Mapping2(load=self.absoluteFilePath.with_name(self.name+'_sequencing_match.mapping'))

    def export_sequencing_match(self):
        self.sequencing_match.save(self.absoluteFilePath.with_name(self.name+'_sequencing_match.mapping'))
        if self.sequencing_match_old is not None:
            self.sequencing_match_old.save(self.absoluteFilePath.with_name(self.name + '_sequencing_match_old.mapping'))

    def find_sequences(self, maximum_distance_file, tuple_size, initial_transformation={},
                       hash_table_distance_threshold=0.01,
                       alpha=0.1, test_radius=10, K_threshold=10e9,
                       magnification_range=None, rotation_range=None,
                       channel=0,
                       nearest_neighbour_match_distance_threshold=25):
        # TODO: Make geometric hashing reflection invariant
        initial_transformation = AffineTransform(**initial_transformation)

        # TODO: make the following line more general and remove bounds dependence in geometric hashing
        source_vertices = self.movie.channel_vertices(channel)
        coordinate_vertices_file = initial_transformation(source_vertices)

        #self.geometric_hash_data = geometric_hash(initial_transform(self.coordinates), maximum_distance_file, tuple_size)

        coordinates = self.coordinates_from_channel(channel)

        #match.destination_index = destination_index
        match = find_match_after_hashing(initial_transformation(coordinates), maximum_distance_file, tuple_size, coordinate_vertices_file,
                                         *self.experiment.geometric_hash_data,
                                         hash_table_distance_threshold, alpha, test_radius, K_threshold,
                                         magnification_range, rotation_range)
        if match:
            match.tile = self.experiment.sequencing_data_for_mapping.tiles[match.destination_index].number
            match.channel = channel
            match.source = self.coordinates
            match.initial_transformation = initial_transformation
            match.transformation = match.transformation @ initial_transformation.params
            match.source_vertices = source_vertices
            match.calculate_inverse_transformation()
            # TODO: Base this on some better criteria
            #match.nearest_neighbour_match(nearest_neighbour_match_distance_threshold)
            self.sequencing_match = match
            self.export_sequencing_match()
            #self.get_all_sequences_from_sequencing_data()

    def find_sequences3(self, distance=15, alpha=0.9, sigma=10, K_threshold=10e2, channel=0,
                        nearest_neighbour_match_distance_threshold=25):

        coordinates = self.coordinates_from_channel(channel)

        match = self.experiment.geometric_hashtable.query(coordinates, distance, alpha, sigma, K_threshold)

        if match:
            # match.destination_index = 0
            match.tile = self.experiment.sequencing_data_for_mapping.tiles[match.destination_index].number
            match.channel = channel

            # TODO: Base this on some better criteria
            #match.nearest_neighbour_match(nearest_neighbour_match_distance_threshold)
            self.sequencing_match = match
            self.export_sequencing_match()
            #self.get_all_sequences_from_sequencing_data()

    def find_sequences_using_stage_coordinates(self, channel=0, show=False, save=True):
        boundaries = self.experiment.tile_boundaries_in_stage_coordinates()
        tile_selection = np.array([self.movie.stage_coordinates > boundaries[:, :, 0],
                                   self.movie.stage_coordinates < boundaries[:, :, 1]])
        tile_index = np.where(np.all(tile_selection, axis=(0, 2)))[0][0]

        source_vertices = self.movie.channels[channel].vertices
        source_vertices_in_stage_coordinates = \
            self.movie.stage_coordinates + np.flip(source_vertices * self.movie.pixel_size, axis=1)


        source_vertices_in_sequencing_coordinates = \
            self.experiment.stage_to_sequencing_mappings[tile_index]\
                .transform_coordinates(source_vertices_in_stage_coordinates)

        match = Mapping2(source_vertices, source_vertices_in_sequencing_coordinates, transformation_type='linear')
        match.direct_match()
        match.source_vertices = source_vertices
        match.destination_vertices = np.array([[1720, 1330], [29720, 1330], [29720, 29330], [1720, 29330]]) # Tile coordinate bounds MiSeq
        match.tile = self.experiment.sequencing_data_for_mapping.tiles[tile_index].number
        match.channel = channel

        match.source = self.coordinates_from_channel(channel)
        match.destination = self.experiment.sequencing_data_for_mapping.tiles[tile_index].coordinates

        self.sequencing_match = match
        if show:
            self.plot_sequencing_match()
        if save:
            self.export_sequencing_match()


    def get_sequencing_data_for_file(self): #, margin=10):

        # sequencing_coordinates_in_image = \
        #     self.sequencing_match.transform_coordinates(self.experiment.sequencing_data.coordinates, inverse=True)
        # coordinate_selection = \
        #     pth.Path(self.sequencing_match.source_vertices).contains_points(sequencing_coordinates_in_image)
        #
        # self.sequencing_data = \
        #     self.experiment.sequencing_data.get_selection(tile=self.sequencing_match.tile,
        #                                                   boolean_selection=coordinate_selection)

        def interpolate(coordinates):
            coordinates = np.vstack([coordinates, coordinates[0]])
            return np.linspace(coordinates, np.roll(coordinates, -1, axis=0)).reshape(-1, 2, order='F')

        source_vertices_interpolated = interpolate(self.sequencing_match.source_vertices)
        destination_vertices_interpolated = self.sequencing_match.transform_coordinates(source_vertices_interpolated)
        self.sequencing_data = \
            self.experiment.sequencing_data.get_selection(tile=self.sequencing_match.tile,
                                                          coordinates_within_vertices=destination_vertices_interpolated)

        self.sequencing_data.export_fastq(self.relativeFilePath)

    def optimize_sequencing_match_using_visible_sequences(self, visible_sequence_names='', distance_threshold=25,
                                                          show=False, save=True, **kwargs):
        self.sequencing_match_old = copy.deepcopy(self.sequencing_match)

        sequencing_data = self.sequencing_data.get_selection(in_name=visible_sequence_names)

        self.sequencing_match.source = self.coordinates_from_channel(self.sequencing_match.channel)
        self.sequencing_match.destination = sequencing_data.coordinates
        self.sequencing_match.nearest_neighbour_match(distance_threshold, **kwargs)

        if show:
            self.plot_sequencing_match()
        if save:
            self.export_sequencing_match()

    def determine_sequences_at_current_coordinates(self, visible_sequence_names='', distance_threshold=None):
        selection = self.sequencing_data.selection(in_name=visible_sequence_names)

        sequence_coordinates_in_file = \
            self.sequencing_match.transform_coordinates(self.sequencing_data.coordinates[selection], inverse=True)

        distances, molecule_indices, sequence_indices = \
            nearest_neighbor_pair(self.coordinates_from_channel(self.sequencing_match.channel),
                                  sequence_coordinates_in_file,
                                  distance_threshold=distance_threshold)

        cum_selection = np.cumsum(selection)-1

        # for molecule_index, sequence_index, distance in zip(molecule_indices, sequence_indices, distances):
        #     self.molecules[molecule_index].sequence_index = np.argmax(cum_selection == sequence_index)
        #     self.molecules[molecule_index].distance_to_sequence = distance

        sequence_indices2 = np.array([np.argmax(cum_selection == sequence_index)
                                     for molecule_index, sequence_index, distance
                                     in zip(molecule_indices, sequence_indices, distances)])

        #molecule_index = self.molecule.to_index()[molecule_indices]
        selected_sequencing_data = self.sequencing_data[sequence_indices2]
        molecule = self.molecule_in_file
        sequencing_dataset = xr.Dataset(
            {
                'sequence':             (('molecule', 'nucleotide'), selected_sequencing_data.sequence),
                'quality':              (('molecule', 'nucleotide'), selected_sequencing_data.quality),
                'distance_to_sequence': ('molecule', distances)
            },
            coords=
            {
                'sequence_name':    ('molecule', selected_sequencing_data.name),
                'molecule':         ('molecule', molecule_indices),
                'sequence_in_file': ('molecule', sequence_indices2)
            }
        )

        sequencing_dataset = sequencing_dataset.reindex_like(molecule.set_index(molecule='molecule_in_file'),
                                        fill_value={'sequence_name': '',
                                                    'sequence_in_file': np.array(np.nan).astype(pd.UInt16Dtype)})
        sequencing_dataset = sequencing_dataset.reset_index('molecule').rename(molecule_='molecule_in_file')

        # Engine netcdf4 has some locking problems.
        sequencing_dataset.to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

    def use_sequences_as_molecules(self):
        self.molecules = []
        file_coordinates = self.sequencing_match.transform_coordinates(self.sequencing_data.coordinates, inverse=True)
        self.set_coordinates_of_channel(file_coordinates, channel=self.sequencing_match.channel)

    def plot_sequencing_match(self):
        #plot_sequencing_match(self.sequencing_match)
        #name = f'Tile: {self.sequencing_tile.name}, File: {str(self.relativeFilePath)}'
        filename = self.name + '_sequencing_mapping'
        print(filename)

        title = f'File: {self.relativePath.joinpath(self.name)} - Tile: {self.sequencing_match.tile}'

        # TODO: Put this somewhere, where it makes sense.
        def MiSeq_pixels_to_um(pixels):
            return 958 / 2800 * (pixels - 1000) / 10

        # TODO: Put this somewhere, where it makes sense.
        # TODO: Make this more general, get this value from tif metadata
        def Fluo_pixels_to_um(pixels):
            return pixels * 0.125

        plot_sequencing_match(self.sequencing_match, self.absoluteFilePath.parent, title, filename,
                              'um', MiSeq_pixels_to_um, Fluo_pixels_to_um)

    # from trace_analysis.plotting import show_point_connections
    # show_point_connections(self.coordinates[molecule_indices], sequence_coordinates_in_file[sequence_indices])


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
    slots = 'sequence_index'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.sequence_index = None

    @property
    def sequence_name(self):
        if self.sequence_index is not None:
            return self.file.sequencing_data.name[self.sequence_index]

    @property
    def sequence(self):
        if self.sequence_index is not None:
            return self.file.sequencing_data.sequence[self.sequence_index, :].tobytes().decode('UTF-8')

    @property
    def sequencing_data(self):
        if self.sequence_index is not None:
            return self.file.sequencing_data[self.sequence_index]


def std_string(value_and_uncertainty):
    value, uncertainty = value_and_uncertainty
    exponent = math.floor(math.log10(uncertainty))
    decimals = -exponent if exponent < 0 else 0
    units = exponent if exponent > 0 else 1
    return f'{np.round(value,-exponent):{units}.{decimals}f}±{np.round(uncertainty,-exponent):{units}.{decimals}f}'

def mean_and_std_string(values):
    return std_string(np.mean(values), np.std(values))