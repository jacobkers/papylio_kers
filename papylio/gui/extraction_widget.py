
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, \
    QComboBox, QLineEdit, QLabel, QFormLayout, QFrame
from PySide2.QtCore import Qt

import matplotlib.pyplot as plt


from papylio import File
from papylio.gui.common_layouts import (Expander, ImageCanvas, HelpDialog,
                                        build_control_layouts,make_push_button,
                                        build_form,build_parameters_input)

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
        frame_spot_detection = QFrame()
        frame_spot_detection.setFrameShape(QFrame.StyledPanel)
        frame_spot_detection.setFrameShadow(QFrame.Plain)
        frame_spot_detection.setLineWidth(2)

        form_spot_detection = QFormLayout(frame_spot_detection)
        spot_detection_label = QLabel("spot_detection")
        spot_detection_label.setToolTip("Tune spot_detection")
        form_spot_detection.addRow(spot_detection_label)
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


        #non-dynamic extraction settings
        #-------------------------------------------------------

        # advanced extraction grid layout:


        #this one collects the various settings frames in horizontal fashion
        advanced_layout = QHBoxLayout()
        advanced_layout.setAlignment(Qt.AlignLeft)

        #first frame: general settings---------------------------------
        frame_general = QFrame()
        frame_general.setFrameShape(QFrame.StyledPanel)
        frame_general.setFrameShadow(QFrame.Plain)
        frame_general.setLineWidth(2)
        advanced_general_layout = QFormLayout(frame_general)
        general_label = QLabel("general ")
        advanced_general_layout.addRow(general_label)
        self.button_extract_chan_combobox = QComboBox() #channel
        button_extract_chan_options = ['donor', 'acceptor']
        self.button_extract_chan_combobox.addItems(button_extract_chan_options)
        advanced_general_layout.addRow("channels", self.button_extract_chan_combobox)
        self.button_extract_illumination_entry = QLineEdit()
        self.button_extract_illumination_entry.setPlaceholderText("0")
        advanced_general_layout.addRow("illumination:", self.button_extract_illumination_entry)

        frame_projection = QFrame()
        frame_projection.setFrameShape(QFrame.StyledPanel)
        frame_projection.setFrameShadow(QFrame.Plain)
        frame_projection.setLineWidth(2)
        advanced_projection_layout = QFormLayout(frame_projection)
        projection_label = QLabel("projection")
        advanced_projection_layout.addRow(projection_label)
        self.button_extract_method_combobox = QComboBox()
        button_extract_method_options = ['by_channel', 'average_channels', 'sum_channels']
        self.button_extract_method_combobox.addItems(button_extract_method_options)
        advanced_projection_layout.addRow("method:", self.button_extract_method_combobox)

        #add general frame to tab layout
        advanced_layout.addWidget(frame_general)
        advanced_layout.addWidget(frame_projection)
        advanced_layout.addWidget(frame_spot_detection)


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
        extraction_controls_layout.setAlignment(Qt.AlignTop)
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


