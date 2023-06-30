from trace_analysis.collection import Collection
from trace_analysis.file import File
from trace_analysis.netcdf_operations import merge_datasets, reorder_datasets_using_sequence_subset

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# TODO: Make sure that the collection can only contain File objects
class FileCollection(Collection):
    def __getattr__(self, item):
        attrs = super(FileCollection, self).__getattr__(item)
        if callable(attrs):  # Test if attrs is a function, if this is the case add another function to concatenate possible xarray objects.
            def f2(*args, **kwargs):
                output = attrs(*args, **kwargs)
                if output is not None and isinstance(output[0], xr.DataArray):
                    output = xr.concat(output, dim=output[0].dims[0])
                return output
            return f2

        elif isinstance(attrs[0], xr.DataArray):
            xr.concat(attrs, dim=attrs[0].dims[0])
        elif isinstance(attrs[0], xr.Dataset):
            if 'molecule' in attrs[0].dims:
                attrs = xr.concat(attrs, dim='molecule')
        elif item == 'dwells':
            #TODO is this still necessary or is this done by the xr.DataArray elif?
            attrs = xr.concat(attrs, dim='dwell')
        return attrs

    @property
    def experiment(self):
        return self[0].experiment

    def select(self, search_string, variable='relativeFilePath'):
        # TODO: Make this accept keyword arguments where the key is the variable and the item is the search_string.
        # TODO: Make this accept mulitple keyword arguments.
        return self[getattr(self, variable).str.regex(search_string)]

    def show_histogram(self, *args, **kwargs):
        figure, axis = File.show_histogram(self.serial, *args, **kwargs)
        axis.set_title('')
        return figure, axis

    def histogram_FRET_intensity_total(self, selected=False, frame_range=None, average=True, axis=None, **hist2d_kwargs):
        FRET_values = self.get_FRET(selected=selected, frame_range=frame_range, average=average).values.flatten()
        total_intensity_values = self.get_intensity_total(selected=selected, frame_range=frame_range, average=average).values.flatten()

        if axis is None:
            figure, axis = plt.subplots()
        axis.hist2d(FRET_values, total_intensity_values, range=((-0.05, 1.05), None), **hist2d_kwargs)
        axis.set_xlabel('FRET')
        axis.set_ylabel('Total intensity (a.u.)')

        return axis

    @property
    def cycle_time(self):
        return self[0].cycle_time

    def analyze_dwells_combined(self, *args, **kwargs):
        return File.analyze_dwells(self.serial, *args, **kwargs)

    def print(self):
        for i, file in enumerate(self):
            print(f"{i:3d}.  {file.relativeFilePath}")

    def merge_datasets(self, filepath_out=None, init_file_index=0, with_selected_only=False, with_sequence_only=False):
        #TODO: remove sequencing part, or move to the sequencing plugin
        if filepath_out is None:
            filepath_out = self[0].absoluteFilePath.parent / 'merged_dataset.nc'
        filepaths_in = self.serial.absoluteFilePath.with_suffix('.nc')
        merge_datasets(filepaths_in, filepath_out, concat_dim='molecule', init_file=filepaths_in[init_file_index],
                       with_selected_only=with_selected_only, with_sequence_only=with_sequence_only)

    def reorder_datasets_using_sequence_subset(self, folderpath_out):
        files = self[self.has_sequencing_match]
        filepaths_in = files.serial.absoluteFilePath.with_suffix('.nc')
        reorder_datasets_using_sequence_subset(filepaths_in, folderpath_out, concat_dim='molecule')