#1 general settings box ---------------------------------
from papylio.gui.common_layouts import get_button_value


def pass_buttons_to_config_for_find_coordinates(self):
#semi-temporal function to line up classic config with buttons in gui.
# Later to be replaced by direct kwargs outputfor 'find_coordinates'
    #1.1 channels:
    self.parent.experiment.configuration['find_coordinates']['channels'][0] = get_button_value(
        self.button_general_channel)
    #1.2 illuminations
    self.parent.experiment.configuration['find_coordinates']['illumination'] = get_button_value(
        self.button_general_illumination)
    #1.3 projection_type:
    self.parent.experiment.configuration['find_coordinates']['projection_type']=get_button_value(
        self.button_general_projection_type)
    #1.4 method:
    self.parent.experiment.configuration['find_coordinates']['method'] = get_button_value(
        self.button_general_method)

    #2. projection image box--------------------------------------
    #2.1 type
    self.parent.experiment.configuration['find_coordinates']['projection_image']['projection_type'] = get_button_value(
    self.button_projection_image_type)
    #2.2 frame range
    self.parent.experiment.configuration['find_coordinates']['projection_image']['frame_range'] = get_button_value(
    self.button_projection_image_frame_range)

    #2.3 illumination
    self.parent.experiment.configuration['find_coordinates']['projection_image']['illumination'] = get_button_value(
    self.button_projection_image_illumination)

    #3. peakfinding (from dynamic form building, see extraction def)

    #4. coordinate optimization box-------------------------------
    #4.1 margin
    self.parent.experiment.configuration['find_coordinates']['coordinate_optimization'][''] = get_button_value(
    self.button_coordinates_within_margin)
    #4.2 Gauss fit width
    self.parent.experiment.configuration['find_coordinates']['coordinate_optimization']['']  = get_button_value(
    self.button_coordinates_after_gaussian_fit_width)

def pass_buttons_to_config_for_extraction(self):
    # semi-temporal function to line up classic config with buttons in gui.
    # Later to be replaced by direct kwargs outputfor 'find_coordinates'
    #5. extract_traces box-------------------------------
    #5.1 channel:
    self.parent.experiment.configuration['extraction']['channel'] = get_button_value(
    self.button_extract_channel)
    #5.2 mask
    self.parent.experiment.configuration['extraction']['mask_size'] = get_button_value(
    self.button_extract_mask_size)
    #5.3 neighbourhood_size: 11
    self.parent.experiment.configuration['extraction']['neighbourhood_size'] = get_button_value(
    self.button_extract_neighbourhood_size)
    #5.4 background
    self.parent.experiment.configuration['extraction']['subtract_background'] = get_button_value(
    self.button_extract_subtract_background)
    #5.5 illumination correction
    self.parent.experiment.configuration['extraction']['correct_illumination'] = get_button_value(
    self.button_extract_correct_illumination)
