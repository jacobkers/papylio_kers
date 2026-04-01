
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
        self.methods_donor = {}
        self.method_forms_donor = {}  # method_name -> (widget, inputs)

        self.parent.model.itemChanged.connect(self.onItemChange)


        #build a flexible form for the donor channel:-------------------------
        frame_donor = QFrame()
        frame_donor.setFrameShape(QFrame.StyledPanel)
        frame_donor.setFrameShadow(QFrame.Plain)
        frame_donor.setLineWidth(2)

        form_donor = QFormLayout(frame_donor)
        donor_label = QLabel("DONOR spot detection")
        donor_label.setToolTip("Tune the peak detection for the donor channel")
        form_donor.addRow(donor_label)
        # --- Method selector donor---
        self.method_selector_donor = QComboBox()
        self.method_selector_donor.setToolTip("Choose peak_find_method")
        self.method_selector_donor.currentTextChanged.connect(self._update_method_panel_donor)
        form_donor.addRow("Method:", self.method_selector_donor)
        # --- Dynamic options container donor ---
        self.stack_donor = QWidget()
        self.stack_donor_layout = QVBoxLayout(self.stack_donor)
        self.stack_donor_layout.setContentsMargins(0, 0, 0, 0)
        form_donor.addRow("Options:", self.stack_donor)


        #genral extraction settings
        #-------------------------------------------------------
        #         #what we want to inspect:
        #         """"
        #         find_coordinates:
        #         channels:
        #         - donor
        #         illumination: 0
        #         projection_type: average
        #         method: by_channel
        #         projection_image:
        #         projection_type: average
        #         frame_range:
        #         - 0
        #         - 20
        #         illumination: 0
        #
        #     sliding_window:
        #     use_sliding_window: false
        #     frame_increment: 20
        #     minimal_point_separation: 2
        #
        # peak_finding [inspect]:
        #     method: local - maximum - auto
        #     coordinate_optimization:
        #     coordinates_within_margin:
        #     margin: 10
        # coordinates_after_gaussian_fit:
        #     gaussian_width: 3
        #     background:
        #     method: ROI_minimum
        #     frames_for_background:
        #     first_frame: 0
        #     last_frame: 9
        #         """
        #main extraction buttons definition:
        #basic extraction grid layout:
        #donor_acceptor block
        donor_acceptor_layout = QHBoxLayout()
        donor_acceptor_layout.addWidget(frame_donor)

        self.button_extraction_label = QLabel("method:")
        self.button_extraction_method_combobox = QComboBox()
        self.button_extraction_options = ['icp', 'nn']
        self.button_extraction_method_combobox.addItems(self.button_extraction_options)
        self.button_extraction_dist_treshold = QLineEdit()
        self.button_extraction_dist_treshold.setPlaceholderText("distance_treshold")
        self.button_transformation = QComboBox()
        self.button_transformation_options = ['polynomial', 'linear', 'nonlinear']
        self.button_transformation.addItems(self.button_transformation_options)
        self.button_initial_translation = QComboBox()
        self.button_initial_translation_options = ['width/2' , [1024,0]]
        self.button_initial_translation.addItems(self.button_initial_translation_options)


        # advanced extraction grid layout:
        frame = QFrame()
        frame.setFrameShape(QFrame.StyledPanel)
        frame.setFrameShadow(QFrame.Plain)
        frame.setLineWidth(2)

        extraction_advanced_layout = QVBoxLayout()
        extraction_advanced_layout.setAlignment(Qt.AlignLeft)

        extraction_advanced2_layout = QFormLayout(frame)
        extraction_advanced2_layout.addRow("Method", self.button_extraction_method_combobox)
        extraction_advanced2_layout.addRow("Distance treshold:",self.button_extraction_dist_treshold)
        extraction_advanced2_layout.addRow("Transformation_type:", self.button_transformation)
        extraction_advanced2_layout.addRow("Initial_translation:", self.button_initial_translation)

        extraction_advanced_layout.addWidget(frame)
        extraction_advanced_layout.addLayout(extraction_advanced2_layout)


        #build panel layout:
        extraction_advanced = Expander("Advanced")
        extraction_advanced.setContentLayout(extraction_advanced_layout)

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
        and build forms for donor and acceptor channels"""

        #donor-------------------------:
        form_widget_donor, inputs_donor = build_form(func)
        self.methods_donor[name] = func
        self.method_forms_donor[name] = (form_widget_donor, inputs_donor)
        self.method_selector_donor.addItem(name) #add options to appropriate selector box
        if self.method_selector_donor.count() == 1: # First registered method becomes default
            self._update_method_panel_donor(name)



    def _update_method_panel_donor(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_donor_layout.count())):
            widget = self.stack_donor_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms_donor:
            form_widget, _ = self.method_forms_donor[name]
            self.stack_donor_layout.addWidget(form_widget)

    def perform_extraction(self, t):

        selected_files = self.parent.experiment.selectedFiles

        #get methods and corresponding parameters
        #donor:
        method_name_donor = self.method_selector_donor.currentText()
        _, inputs_donor = self.method_forms_donor[method_name_donor]
        method_name_acceptor = self.method_selector_acceptor.currentText()
        _, inputs_acceptor = self.method_forms_acceptor[method_name_acceptor]

        # Collect args for peak finding
        donor_kwargs=build_parameters_input(method_name_donor, inputs_donor)
        acceptor_kwargs = build_parameters_input(method_name_acceptor, inputs_acceptor)

        #jk-read in a default configuration and allocate GUI values:
        panel_config= self.parent.experiment.configuration['extraction']
        panel_config['method']= self.button_extraction_method_combobox.currentText()
        panel_config['distance_threshold']= self.button_extraction_dist_treshold.text
        panel_config['peak_finding']['donor'] = donor_kwargs
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


