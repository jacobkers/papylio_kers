#1 general settings box ---------------------------------


def pass_buttons_to_config_for_find_coordinates(self)
    #1.1 channels:
    #1.2 illuminations
    self.parent.experiment.configuration['find_coordinates']['channels'][0] = get_button_value(self.button_general_channel)
    self.parent.experiment.configuration['find_coordinates']['illumination'] = get_button_value(self.button_general_illumination)

    #1.3 projection_type:
    self.button_general_projection_type
    #1.4 method:
    self.button_general_method

    #2. projection image box--------------------------------------

    #2.1 type
    self.button_projection_image_type
    #2.2 frame range
    self.button_projection_image_frame_range

    #2.3 illumination
    self.button_projection_image_illumination

    #3. peakfinding (from dynamic form building, see extraction def)

    #4. coordinate optimization box-------------------------------
    #4.1 margin
    self.button_coordinates_within_margin
    #4.2 Gauss fit width
    self.button_coordinates_after_gaussian_fit_width

    def pass_buttons_to_config_for_extraction(self)
    #5. extract_traces box-------------------------------
    #5.1 channel:
    self.parent.experiment.configuration['extraction']['channel']
    self.button_extract_channel
    #5.2 mask
    self.parent.experiment.configuration['extraction']['mask_size']
    self.button_extract_mask_size
    #5.3 neighbourhood_size: 11
    self.parent.experiment.configuration['extraction']['neighbourhood_size']
    self.button_extract_neighbourhood_size
    #5.4 background
    self.parent.experiment.configuration['extraction']['subtract_background']
    self.button_extract_subtract_background
    #5.5 illumination correction
    self.parent.experiment.configuration['extraction']['correct_illumination']
    self.button_extract_correct_illumination
    #in config:




