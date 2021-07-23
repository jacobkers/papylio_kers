import numpy as np
import pandas as pd
import xarray as xr
import tifffile

from trace_analysis.movie.movie import Movie


class TifMovie(Movie):
    def __init__(self, arg, *args, **kwargs):
        super().__init__(arg, *args, **kwargs)
        
        self.writepath = self.filepath.parent
        self.name = self.filepath.with_suffix('').name
        
        self.threshold = {  'view':             (0,200),
                            'point-selection':  (45,25)
                            }

        self.read_header()
        self.create_frame_info()  # Possibly move to Movie later on

    def _read_header(self):
        with tifffile.TiffFile(self.filepath) as tif:
            tif_tags = {}
            for tag in tif.pages[0].tags.values():
                name, value = tag.name, tag.value
                tif_tags[name] = value
            self.width = tif_tags['ImageWidth']
            self.height = tif_tags['ImageLength']
            self.number_of_frames = len(tif.pages)
            self.data_type = np.dtype(f"uint{tif_tags['BitsPerSample']}")

            self.datetime = pd.to_datetime([page.tags['DateTime'].value for page in tif.pages])
            self.time = xr.DataArray((self.datetime-self.datetime[0]).total_seconds(), dims='frame',
                                     coords={}, attrs={'units': 's'})

            try:
                if tif.metaseries_metadata:
                    pixel_size_x = tif.metaseries_metadata['PlaneInfo']['spatial-calibration-x']
                    pixel_size_y = tif.metaseries_metadata['PlaneInfo']['spatial-calibration-y']
                    self.pixel_size = np.array([pixel_size_x, pixel_size_y])
                    self.pixel_size_unit = tif.metaseries_metadata['PlaneInfo']['spatial-calibration-units']
                    stage_position_x = tif.metaseries_metadata['PlaneInfo']['stage-position-x']
                    stage_position_y = tif.metaseries_metadata['PlaneInfo']['stage-position-y']
                    self.stage_coordinates = np.array([[stage_position_x, stage_position_y]])
                    self.stage_coordinates_in_pixels = self.stage_coordinates / self.pixel_size

            except AttributeError:
                pass

    def _read_frame(self, frame_number):
        with tifffile.TiffFile(self.filepath) as tif:
            tifpage = tif.pages
            if self.number_of_frames == 1:
                # return -1,0,0,0
                im = tifpage[0].asarray()
            elif (self.number_of_frames - 1) >= frame_number:
                im = tifpage[frame_number].asarray()
            else:
                im = tifpage[self.number_of_frames - 1].asarray()
                print('pageNb out of range, printed image {0} instead'.format(self.number_of_frames))
        return im


if __name__ == "__main__":
    movie = TifMovie(r'.\Example_data\tif\movie.tif')

