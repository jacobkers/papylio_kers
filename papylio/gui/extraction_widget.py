
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, \
    QComboBox, QLineEdit, QSpinBox, QFormLayout
from PySide2.QtCore import Qt

import matplotlib.pyplot as plt


from papylio import File
from papylio.gui.common_layouts import (Expander, ImageCanvas, HelpDialog,
                                        build_control_layouts,make_push_button,
                                        build_form,build_parameters_input, Group_Box)

from papylio.peak_finding import (find_peaks_absolute_threshold,
                                  find_peaks_adaptive_threshold,
                                  find_peaks_local_maximum,
                                  find_peaks_local_maximum_auto,
                                  find_peaks_relative_local_maximum)

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)


class ExtractionWidget(QWidget):
    def __init__(self, parent=None):
        super(ExtractionWidget, self).__init__(parent)
        self.parent = parent
        self.methods_spot_detection = {}
        self.method_forms_spot_detection = {}  # method_name -> (widget, inputs)

        self.parent.model.itemChanged.connect(self.onItemChange)



        #build a flexible form for the spot_detection channel:-------------------------
        group_spot_detection = Group_Box("Spot Detection")
        form_spot_detection = QFormLayout(group_spot_detection)

        # --- Method selector spot_detection---
        self.method_selector_spot_detection = QComboBox()
        self.method_selector_spot_detection.setToolTip("Choose peak_find_method")
        self.method_selector_spot_detection.currentTextChanged.connect(self._update_method_panel_spot_detection)
        form_spot_detection.addRow("Method:", self.method_selector_spot_detection)
        # --- Dynamic options container spot_detection ---
        self.stack_spot_detection = QWidget()
        self.stack_spot_detection_layout = QVBoxLayout(self.stack_spot_detection)
        self.stack_spot_detection_layout.setContentsMargins(0, 0, 0, 0)
        form_spot_detection.addRow("Options:", self.stack_spot_detection)

        #this one collects the various settings frames in horizontal fashion
        advanced_layout = QHBoxLayout()
        advanced_layout.setAlignment(Qt.AlignLeft)

        #general settings box ---------------------------------
        frame_general = Group_Box("General")
        advanced_general_layout = QFormLayout(frame_general)
        self.button_extract_chan_combobox = QComboBox() #channel
        self.button_extract_chan_combobox.addItems(['donor', 'acceptor'])
        advanced_general_layout.addRow("channels", self.button_extract_chan_combobox)
        self.button_extract_illumination_gen = QLineEdit()
        self.button_extract_illumination_gen.setPlaceholderText("0")
        advanced_general_layout.addRow("illumination:", self.button_extract_illumination_gen)

        #projection box--------------------------------------
        frame_projection = Group_Box("Projection")
        advanced_projection_layout = QFormLayout(frame_projection)
        self.button_projection_channel_combobox = QComboBox()
        button_project_chan_options = ['by_channel', 'average_channels', 'sum_channels']
        self.button_projection_channel_combobox.addItems(button_project_chan_options)
        advanced_projection_layout.addRow("channels:", self.button_projection_channel_combobox)

        self.button_projection_combobox = QComboBox()
        self.button_projection_combobox.addItems(['average', 'maximum'])
        advanced_projection_layout.addRow("projection type:", self.button_projection_combobox)

        self.button_projection_frame_range = QLineEdit()
        self.button_projection_frame_range.setPlaceholderText("[0, 20]")
        advanced_projection_layout.addRow("frame range:", self.button_projection_frame_range)

        self.button_illumination_proj = QLineEdit()
        self.button_illumination_proj.setPlaceholderText("0")
        advanced_projection_layout.addRow("illumination:", self.button_illumination_proj)

        # coordinate optimization box-------------------------------
        frame_coord_opt = Group_Box("Coordinate Optimization")
        advanced_coord_opt_layout = QFormLayout(frame_coord_opt)
        # more buttons to be added here:
        self.button_coordinates_within_margin = QLineEdit()
        self.button_coordinates_within_margin.setPlaceholderText("10")
        advanced_coord_opt_layout.addRow("within_margin:",
                                         self.button_coordinates_within_margin)
        self.button_coordinates_after_gaussian_fit_width = QLineEdit()
        self.button_coordinates_after_gaussian_fit_width.setPlaceholderText("3")
        advanced_coord_opt_layout.addRow("gaussian_fit:" ,self.button_coordinates_after_gaussian_fit_width)


        # extract_traces box-------------------------------
        frame_extract_traces = Group_Box("Extract Traces")
        advanced_extract_traces_layout = QFormLayout(frame_extract_traces)
        #channel:
        self.button_xtr_channel = QComboBox()
        self.button_xtr_channel.addItems(['all'])
        advanced_extract_traces_layout.addRow("channel:", self.button_xtr_channel)
        #mask
        self.button_mask_size = QLineEdit()
        self.button_mask_size.setPlaceholderText("11")
        self.button_mask_size.setToolTip("float number or presets: TIR-T, TIR-V, TIR-S 1.5x 2x2, TIR-S 1x 2x2, BN-TIRF")
        advanced_extract_traces_layout.addRow("mask size:", self.button_mask_size)
        #neighbourhood_size: 11
        self.button_neighbourhood_size= QSpinBox()
        self.button_neighbourhood_size.setValue(11)
        advanced_extract_traces_layout.addRow("neighbourhood_size:", self.button_neighbourhood_size)

        self.button_subtract_background = QComboBox()
        self.button_subtract_background.addItems(['False', 'True'])
        advanced_extract_traces_layout.addRow("subtract_background:", self.button_subtract_background)

        self.button_correct_illumination = QComboBox()
        self.button_correct_illumination.addItems(['False', 'True'])
        advanced_extract_traces_layout.addRow("subtract_background:", self.button_correct_illumination)


        #add general frame to tab layout
        advanced_layout.addWidget(frame_general)
        advanced_layout.addWidget(frame_projection)
        advanced_layout.addWidget(group_spot_detection)
        advanced_layout.addWidget(frame_coord_opt)
        advanced_layout.addWidget(frame_extract_traces)

        #build panel layout:
        extraction_advanced = Expander("Advanced")
        extraction_advanced.setContentLayout(advanced_layout)


        #main action:
        extraction_controls = build_control_layouts(
            [make_push_button('Find coordinates', self.find_coordinates,"extraction selected file(s)"),
             make_push_button('Extract traces', self.extract_traces, "extraction selected file(s)"),
             make_push_button('Help', self.show_help, None)])

        #collect:
        extraction_controls_layout = QVBoxLayout()
        extraction_controls_layout.setAlignment(Qt.AlignRight)
        extraction_controls_layout.addWidget(extraction_advanced)
        extraction_controls_layout.addWidget(extraction_controls)

        #pack in widget:
        self.extraction_controls = QWidget()
        self.extraction_controls.setLayout(extraction_controls_layout)
        self.extraction_controls.setMinimumWidth(150)

        #add all to tab:
        extraction_tab_layout = QHBoxLayout()
        extraction_tab_layout.addWidget(self.extraction_controls)

        self.setLayout(extraction_tab_layout)

        #collect peak finding methods for building flexible GUI forms
        self.register_method('local-maximum-auto', find_peaks_local_maximum_auto)  #default  in GUI
        self.register_method('absolute-threshold', find_peaks_absolute_threshold)
        self.register_method('adaptive-threshold', find_peaks_adaptive_threshold)
        self.register_method('local-maximum', find_peaks_local_maximum)
        self.register_method('relative-local-maximum', find_peaks_relative_local_maximum)


    def update_plots(self):
        selected_files = self.parent.experiment.selectedFiles + [None]

    def onItemChange(self, item):
        if isinstance(item.data(), File):
            file = item.data()
            file.isSelected = (True if item.checkState() == Qt.Checked else False)
            print(f'{file}: {file.isSelected}')

        else:
            self.update = False
            for i in range(item.rowCount()):
                item.child(i).setCheckState(item.checkState())
            self.update = True

        if self.update:
            self.update_plots()

    # -------------------------------------------------------------------------
    # Register methods dynamically and create their forms
    # -------------------------------------------------------------------------
    def register_method(self, name, func):
        """Register a peak finding method, introspect arguments,
        and build forms for spot_detection and acceptor channels"""

        #spot_detection-------------------------:
        form_widget_spot_detection, inputs_spot_detection = build_form(func)
        self.methods_spot_detection[name] = func
        self.method_forms_spot_detection[name] = (form_widget_spot_detection, inputs_spot_detection)
        self.method_selector_spot_detection.addItem(name) #add options to appropriate selector box
        if self.method_selector_spot_detection.count() == 1: # First registered method becomes default
            self._update_method_panel_spot_detection(name)



    def _update_method_panel_spot_detection(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_spot_detection_layout.count())):
            widget = self.stack_spot_detection_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms_spot_detection:
            form_widget, _ = self.method_forms_spot_detection[name]
            self.stack_spot_detection_layout.addWidget(form_widget)

    def perform_extraction(self, t):

        selected_files = self.parent.experiment.selectedFiles

        #get methods and corresponding parameters
        #spot_detection:
        method_name_spot_detection = self.method_selector_spot_detection.currentText()
        _, inputs_spot_detection = self.method_forms_spot_detection[method_name_spot_detection]
        method_name_acceptor = self.method_selector_acceptor.currentText()
        _, inputs_acceptor = self.method_forms_acceptor[method_name_acceptor]

        # Collect args for peak finding
        spot_detection_kwargs=build_parameters_input(method_name_spot_detection, inputs_spot_detection)
        acceptor_kwargs = build_parameters_input(method_name_acceptor, inputs_acceptor)

        #jk-read in a default configuration and allocate GUI values:
        panel_config= self.parent.experiment.configuration['extraction']
        panel_config['method']= self.button_extraction_method_combobox.currentText()
        panel_config['distance_threshold']= self.button_extraction_dist_treshold.text
        panel_config['peak_finding']['spot_detection'] = spot_detection_kwargs
        panel_config['peak_finding']['acceptor'] = acceptor_kwargs

        if selected_files:
            selected_files.serial.perform_extraction(**panel_config)
            self.update_plots()


    #TODO: added from main, to be edited
    def find_coordinates(self):
        selected_files = self.parent.experiment.selectedFiles
        if selected_files:
            selected_files.movie.determine_spatial_background_correction(use_existing=True)
            selected_files.find_coordinates()

            def find_coordinates(self):
                selected_files = self.experiment.selectedFiles
                if selected_files:
                    selected_files.movie.determine_spatial_background_correction(use_existing=True)
                    selected_files.find_coordinates()
                    self.image_canvas.refresh()
                    self.update_plots()

            def extract_traces(self):
                selected_files = self.experiment.selectedFiles
                if selected_files:
                    selected_files.extract_traces()
                    # self.image_canvas.refresh()
                    self.update_plots()
            self.update_plots()

    def extract_traces(self):
        selected_files = self.parent.experiment.selectedFiles
        if selected_files:
            selected_files.extract_traces()
            # self.image_canvas.refresh()
            self.update_plots()


    def show_help(self):
        help_text = """
            <html>
              <body style="font-family: sans-serif; font-size: 10pt;">

                <h2>Extraction</h2>

                <p>
                  Find spot coordinates and extract traces.
                </p>

                <ul>
                  <li>run as-is using the defaults</li>
                  <li>use 'advanced' to change settings  2</li>
                  <li>press 'Find coordinates' and 'Extract traces' 3
                </ul>

                <p>
                  For more help, see
                  <a href="https://papylio.readthedocs.io/en/stable/user_guide/molecule_localization.html">
                    Molecule Localization in Papylio
                  </a>.
                </p>

                <h3>Example</h3>

                <p>
                  Settings examples can be found in the documentation.
                </p>

              </body>
            </html>
            """

        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()


