
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, \
    QComboBox, QLineEdit, QSpinBox, QFormLayout, QButtonGroup, QRadioButton, QCheckBox, QMessageBox
from PySide2.QtCore import Qt, Signal
from papylio import File
from papylio.gui.common_layouts import (Expander, HelpDialog,Group_Box,
                                        build_control_layouts,make_push_button,
                                        build_form,build_parameters_input, get_button_value)

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

#for method registry:
from scipy.ndimage import gaussian_filter,median_filter,minimum_filter

"""
There are three types of background subtraction that can be performed. 
A temporal background subtraction that corrects variations in background over time 
(but not in x and y); 

a spatial background correction that corrects variations over the x and y 
(but not over time); 

and a general background correction, which is a single value for the entire movie. 

Any of these can be applied. 
If more than one are used, they have to be determined in the order 
1) temporal, 
2) spatial, 
3) general,
 
... as this is also the order in which they are applied. 
Determining the corrections in a different order, 
will delete the corrections that should have been determined later.
(from Papylio docs, IS)

"""
class MovieCorrectionsWidget(QWidget):
    def __init__(self, parent=None):
        super(MovieCorrectionsWidget, self).__init__(parent)

        # self.parent = parent
        self.methods_spatial_background = {}
        self.method_forms_spatial_background = {}
        # movie corrections box 1: shading ---------------------------------
        frame_shading_correction = Group_Box(title="Shading", highlight=False)
        frame_shading_correction.setToolTip('"see help')
        shading_correction_layout = QFormLayout(frame_shading_correction)
        # method:
        self.method_shading = QComboBox()
        self.method_shading.setToolTip("Choose method")
        self.method_shading.addItems(['any', 'any'])
        #frame range:
        self.button_frame_range = QLineEdit()
        self.button_frame_range.setText("[0, 20]")
        # skipbox:
        self.skip_shading_checkbox = QCheckBox()
        # fill box:
        shading_correction_layout.addRow("method:", self.method_shading)
        # shading_correction_layout.addRow("frame range:", self.button_frame_range)
        shading_correction_layout.addRow("skip", self.skip_shading_checkbox)

        #movie corrections box 2: temporal ---------------------------------
        frame_temporal_correction = Group_Box(title="Temporal", highlight=False)
        frame_temporal_correction.setToolTip('"see help')
        temporal_correction_layout = QFormLayout(frame_temporal_correction)
        #method:
        self.method_temporal = QComboBox()
        self.method_temporal.setToolTip("Choose method")
        self.method_temporal.addItems(['BaSiC', 'any'])
        # #frame range:
        # self.button_frame_range = QLineEdit()
        # self.button_frame_range.setText("[0, 20]")
        #skipbox:
        self.skip_temporal_checkbox = QCheckBox()
        # fill box:
        temporal_correction_layout.addRow("method:", self.method_temporal)
        #temporal_correction_layout.addRow("frame range:", self.button_frame_range)
        temporal_correction_layout.addRow("skip", self.skip_temporal_checkbox)



        # dynamic form building box 3: spatial
        frame_spatial_correction = Group_Box(title="Spatial", highlight=False)
        frame_spatial_correction.setToolTip('"see help')
        form_spatial_background = QFormLayout(frame_spatial_correction)
        # --- Method selector spatial background---
        self.method_selector_spatial_background = QComboBox()
        self.method_selector_spatial_background.setToolTip("Choose peak_find_method")
        self.method_selector_spatial_background.currentTextChanged.connect(self._update_method_panel_spatial_background)
        form_spatial_background.addRow("Method:", self.method_selector_spatial_background)
        # --- Dynamic options container spatial background ---
        self.stack_spatial_background = QWidget()
        self.stack_spatial_background_layout = QVBoxLayout(self.stack_spatial_background)
        self.stack_spatial_background_layout.setContentsMargins(0, 0, 0, 0)
        form_spatial_background.addRow("Options:", self.stack_spatial_background)
        # frame range:
        self.button_spatial_frame_range = QLineEdit()
        self.button_spatial_frame_range.setText("[0, 20]")
        form_spatial_background.addRow("frame range", self.button_spatial_frame_range)
        #skipbox:
        self.skip_spatial_checkbox = QCheckBox()
        # fill box:
        form_spatial_background.addRow("skip", self.skip_spatial_checkbox)

        # movie corrections box 4: general ---------------------------------
        frame_general_correction = Group_Box(title="General", highlight=False)
        frame_general_correction.setToolTip('"see help')
        general_correction_layout = QFormLayout(frame_general_correction)
        # method:
        self.method_general = QComboBox()
        self.method_general.setToolTip("Choose method")
        self.method_general.addItems(['BaSiC', 'any'])
        self.skip_general_checkbox = QCheckBox()
        # fill box:
        general_correction_layout.addRow("method:", self.method_general)
        general_correction_layout.addRow("skip", self.skip_general_checkbox)

        this_tab_layout = QHBoxLayout()
        this_tab_layout.setAlignment(Qt.AlignLeft)
        this_tab_layout.addWidget(frame_shading_correction)
        this_tab_layout.addWidget(frame_temporal_correction)
        this_tab_layout.addWidget(frame_spatial_correction)
        this_tab_layout.addWidget(frame_general_correction)


        # main action:
        start_help_button = build_control_layouts([
            make_push_button('Apply', self.apply_corrections, "Apply correction(s)"),
            make_push_button('Help', self.show_main_help, None)])



        start_tab_layout = QVBoxLayout()
        start_tab_layout.addLayout(this_tab_layout)
        start_tab_layout.addStretch()
        start_tab_layout.addWidget(start_help_button)
        self.setLayout(start_tab_layout)

        self.file = None
        self.experiment = None
        # TODO: building
        # collect spatial filter methods for building flexible GUI forms
        #'median_filter', 'gaussian_filter', 'minimum_filter'
        self.register_method_for_spatial_background('median_filter', median_filter)
        self.register_method_for_spatial_background('minimum_filter', minimum_filter)
        self.register_method_for_spatial_background('gaussian_filter', gaussian_filter)

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        self._file = file

    @property
    def experiment(self):
        return self._experiment

    @experiment.setter
    def experiment(self, experiment):
        self._experiment = experiment
        if experiment is not None:
            self.update_button_settings()

    def apply_corrections(self):
        # TODO: bring in panel settings and decide on shading approach
        file=self.file
        # if 1: #not skip
        #    #do shading correction (note: is this GUI-handy? Needs specification of files...
        #    files_darkfield_correction[0].use_for_darkfield_correction()
        #    # For green illumination
        #    exp.determine_flatfield_and_darkfield_corrections(files_green_laser[::10], method='BaSiC',
        #                                                      illumination_index=0, frame_index=2,
        #                                                      estimate_darkfield=False, l_s=5, l_d=5)
        #    # For red illumination
        #    exp.determine_flatfield_and_darkfield_corrections(files_red_laser_before[::10], method='BaSiC',
        #                                                      illumination_index=1,
        #                                                      frame_index=2, estimate_darkfield=False, l_s=5, l_d=5)
        #
           # Uses the second frame of each file to determine the flatfield correction. l_s and l_d are parameters for the 'BaSiC' algorithm.

        if not self.skip_temporal_checkbox.isChecked(): #not skip
            mth_t=get_button_value(self.method_temporal)
            file.movie.determine_temporal_background_correction(method=mth_t)
        if not self.skip_spatial_checkbox.isChecked():
            # spatial background (Gui_box 3):
            method_name_spatial_background = self.method_selector_spatial_background.currentText()
            _, inputs_spatial_background = self.method_forms_spatial_background[method_name_spatial_background]
            # kwargs for peak finding
            frs = get_button_value(self.button_spatial_frame_range)
            spatial_background_kwargs = build_parameters_input(method_name_spatial_background, inputs_spatial_background)
            file.movie.determine_spatial_background_correction(frame_range=frs, **spatial_background_kwargs)
        if not self.skip_general_checkbox.isChecked():
            file.movie.determine_general_background_correction(method='fit_background_peak')

    def register_method_for_spatial_background(self, name, func):
        #TODO: rewrite this correctly!!
        """Register a peak finding method, introspect arguments,
        and build forms for spot_detection"""
        #spatial background:
        skip_inputs=['input', 'output','footprint', 'origin', 'cval']
        form_widget_spatial_background, inputs_spatial_background = build_form(func,skip_inputs)
        self.methods_spatial_background[name] = func
        self.method_forms_spatial_background[name] = (form_widget_spatial_background, inputs_spatial_background)
        self.method_selector_spatial_background.addItem(name) #add options to appropriate selector box
        if self.method_selector_spatial_background.count() == 1: # First registered method becomes default
            self._update_method_panel_spatial_background(name)

    def _update_method_panel_spatial_background(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_spatial_background_layout.count())):
            widget = self.stack_spatial_background_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms_spatial_background:
            form_widget, _ = self.method_forms_spatial_background[name]
            self.stack_spatial_background_layout.addWidget(form_widget)

    def show_main_help(self):
        help_text = """
                <html>
                  <body style="font-family: sans-serif; font-size: 10pt;">

                    <h2>Background corrections</h2>

                    <p>
                    Optionally, one can perform shading (illumination) correction.
                    Next, there are three types of background subtraction, to be performed in this order: 
                    </p>
                    
                    <p>
                    <ol>
                        <li> Shading correction 
                        <li> Temporal background subtraction 
                        <li> Spatial background correction
                        <li> Single-value background correction
                    </ol>    
                    </p>
                    
                    <h3>Temporal background subtraction</h3>
                    <p> 
                    text on temporal background subtraction
                    </p>
                    
                    <ul>
                      <li>Box 1 steps</li>
                      <li> nnnnn </li>
                    </ul>
                    
                    <p>
                      For background, see
                      <a href="https://papylio.readthedocs.io/en/stable/user_guide/background_subtraction.html">
                        Background Subtraction</a>
                      </a>.
                    </p>
                  </body>
                </html>
                """
        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()
