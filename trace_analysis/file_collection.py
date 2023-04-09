from trace_analysis.collection import Collection
from trace_analysis.file import File

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
                    output = xr.concat(output, dim='molecule')
                return output
            return f2

        elif isinstance(attrs[0], xr.DataArray):
            attrs = xr.concat(attrs, dim='molecule')
        elif item == 'dwells':
            attrs = xr.concat(attrs, dim='dwell')
        return attrs

    @property
    def experiment(self):
        return self[0].experiment

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

    def analyze_dwells(self,  *args, **kwargs):
        return File.analyze_dwells(self.serial, *args, **kwargs)