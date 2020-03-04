import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
# from trace_analysis.experiment import Experiment
# from trace_analysis.file import File
from trace_analysis.mapping.geometricHashing import SequencingDataMapping
from .fastqAnalysis import FastqData
from trace_analysis.mapping.icp import icp, nearest_neighbor_pair


class Experiment:
    def import_sequencing_data(self, fastq_file_path):
        self.fastq_file_path = Path(fastq_file_path)
        self.sequencing_data = FastqData(self.fastq_file_path)

    def select_sequencing_data_for_mapping(self, mapping_sequence, number_of_allowed_mismatches):
        self.mapping_sequence = mapping_sequence
        self.sequencing_data.matches_per_tile(sequence=mapping_sequence)

        number_of_matches = self.sequencing_data.number_of_matches(mapping_sequence)
        # sequencing_data_for_mapping = sequencing_data.select(Nmatch == len(mapping_sequence), copyData=True)
        self.sequencing_data_for_mapping = self.sequencing_data[number_of_matches >=
                                                                (len(mapping_sequence)-number_of_allowed_mismatches)]
        self.sequencing_data_for_mapping.show_tiles()
        self.sequencing_data_for_mapping.export_positions_per_tile()

    def map_files_to_sequencing_data(self, tile, configuration=None):
        sequencing_mapping_path = self.mainPath.joinpath('sequencing_mapping')
        sequencing_mapping_path.mkdir(exist_ok=True)

        presets = {'TIR-I single-molecule':
                       {
                           'initial_image_transformation': {'reflection': 0, 'rotation': np.pi, 'magnification': 3.336},
                           'mapping_configuration': {'mode': 'similarity',
                                                    'nBins': 250,
                                                    'rotationRange': np.array([-180, 180]) / 360 * 2 * np.pi,
                                                    'magnificationRange': np.array([3000, 6000])
                                                    },
                           'bases_findMatch': 500
                        }
                   }
        # rotationRange = np.array([-5,5])/360*2*np.pi+np.pi/2,
        # magnificationRange = np.array([2750,3000])
        # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': -35/180*np.pi, 'magnification':1/125 }
        # seqmap.initial_image_transformation = {'reflection': 0, 'rotation': np.pi, 'magnification': 2*3.336}
        # seqmap.initial_image_transformation = {}

        if type(configuration) is str:
            configuration = presets[configuration]

        tile = self.sequencing_data_for_mapping.get_tile_object(tile=tile)
        self.seqmap = SequencingDataMapping(tile, self.files, sequencing_mapping_path, *configuration['mapping_configuration'])
        self.seqmap.initial_image_transformation = configuration['initial_image_transformation']
        self.seqmap.bases_findMatch = configuration['bases_findMatch']
        self.seqmap.histogram_matches(export=True)

    def sequencing_mapping_to_files(self, minimal_number_of_matching_points):
        self.seqmap.give_matches_to_files(match_threshold=minimal_number_of_matching_points)

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
                         (coordinates[:, 0] < (bounds[0, 1] - margin)),
                         (coordinates[:, 1] > (bounds[1, 0] + margin)),
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

    def plot_sequencing_match(self):
        self.experiment.seqmap.plot_match(self.sequence_match)

    def get_all_sequences_from_sequencing_data(self):
        #raise Warning('Only works on acceptor channel for now')
        coordinate_bounds_file = self.movie.channel_boundaries('a')
        coordinate_bounds_tile = self.sequence_match.transform_coordinates(coordinate_bounds_file)

        self.sequencing_data = \
            self.experiment.sequencing_data.get_selection(tile=int(self.sequence_match.tile.name),
                                                          x=coordinate_bounds_tile[0],
                                                          y=coordinate_bounds_tile[1])
        self.molecules = []
        self.coordinates = self.sequence_match.transform_coordinates(self.sequencing_data.coordinates,
                                                                     inverse=True)



    # This is probably not the way to go
    def give_molecules_closest_sequence(self):
        sequencing_data = self.experiment.sequencing_data
        indices_within_tile = sequencing_data.selection(tile=int(self.sequence_match.tile.name))
        sequencing_data = sequencing_data.select(indices_within_tile, copyData=True)

        #raise Warning("Not implemented for the donor channel yet")
        boundaries = self.sequence_match.transform_coordinates(self.movie.channel_boundaries('a').T).T

        indices_within_bounds = within_bounds(sequencing_data.coordinates, boundaries)
        sequencing_data.select(indices_within_bounds, copyData=False)

        tile_coordinates = sequencing_data.coordinates
        image_coordinates = self.sequence_match.transform_coordinates(self.coordinates)

        from trace_analysis.plotting import scatter_coordinates
        plt.figure()
        scatter_coordinates([tile_coordinates, image_coordinates])


        distances, image_indices, tile_indices = nearest_neighbor_pair(image_coordinates, tile_coordinates)

        distance_threshold = 25
        image_indices = image_indices[distances < distance_threshold]
        tile_indices = tile_indices[distances < distance_threshold]

        print(*[sequence.tostring() for sequence in sequencing_data.sequence[tile_indices]], sep="\n")

        plt.figure()
        scatter_coordinates([tile_coordinates, self.sequence_match.destination, image_coordinates])
        from trace_analysis.plotting import show_point_connections
        show_point_connections(image_coordinates[image_indices], tile_coordinates[tile_indices])

        for molecule in self.molecules: molecule.sequence = np.array([], dtype=bytes)
        for i in np.arange(len(image_indices)):
            self.molecules[image_indices[i]].sequence = sequencing_data.sequence[tile_indices[i]]
            self.molecules[image_indices[i]].sequencing_data = sequencing_data.select(i, copyData=True)

class Molecule:
    @property
    def sequence(self):
        return self.file.sequencing_data.sequence[self.file.experiment.files.index(self.file),:]

    @property
    def sequencing_data(self):
        return self.file.sequencing_data[self.file.experiment.files.index(self.file)]
