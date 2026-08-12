
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, \
    QComboBox, QLineEdit, QSpinBox, QFormLayout, QButtonGroup, QRadioButton, QCheckBox, QMessageBox
from PySide2.QtCore import Qt, Signal
from papylio import File
from papylio.gui.common_layouts import (Expander, HelpDialog,Group_Box,
                                        build_control_layouts,make_push_button,
                                        build_form,build_parameters_input, get_button_value)
from papylio.movie.movie import Channel
#for registry:
from papylio.peak_finding import (find_peaks_absolute_threshold,
                                  find_peaks_adaptive_threshold,
                                  find_peaks_local_maximum,
                                  find_peaks_local_maximum_auto,
                                  find_peaks_relative_local_maximum)

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
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


        # movie corrections box 3: spatial---------------------------------
        frame_spatial_correction = Group_Box(title="Spatial", highlight=False)
        frame_spatial_correction.setToolTip('"see help')
        spatial_correction_layout = QFormLayout(frame_spatial_correction)
        # method:
        self.method_spatial = QComboBox()
        self.method_spatial.setToolTip("choose Filter")
        self.method_spatial.addItems(['gaussian', 'median', 'mean'])
        #frame range:
        self.button_frame_range = QLineEdit()
        self.button_frame_range.setText("[0, 20]")
        #skipbox:
        self.skip_spatial_checkbox = QCheckBox()
        # fill box:
        spatial_correction_layout.addRow("method:", self.method_spatial)
        # spatial_correction_layout.addRow("frame range:", self.button_frame_range)
        spatial_correction_layout.addRow("skip", self.skip_spatial_checkbox)

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

        advanced_layout = QHBoxLayout()
        advanced_layout.setAlignment(Qt.AlignLeft)
        advanced_layout.addWidget(frame_shading_correction)
        advanced_layout.addWidget(frame_temporal_correction)
        advanced_layout.addWidget(frame_spatial_correction)
        advanced_layout.addWidget(frame_general_correction)


        # main action:
        start_help_button = build_control_layouts([
            make_push_button('Apply', self.show_main_help, "Apply correction(s)"),
            make_push_button('Help', self.show_main_help, None)])



        start_tab_layout = QVBoxLayout()
        start_tab_layout.addLayout(advanced_layout)
        start_tab_layout.addStretch()
        start_tab_layout.addWidget(start_help_button)
        self.setLayout(start_tab_layout)

        self.file = None
        self.experiment = None

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        self._file = file
        # if file is None:
        #     self.setDisabled(True)
        # else:
        #     self.setDisabled(False)
        #     self.update_button_settings()

    @property
    def experiment(self):
        return self._experiment

    @experiment.setter
    def experiment(self, experiment):
        self._experiment = experiment
        if experiment is not None:
            self.update_button_settings()

    def update_button_settings(self):
        self.button_movie_rotation.blockSignals(True)
        self.button_movie_rotation.setValue(self.experiment.configuration['movie']['rot90'])
        self.button_movie_rotation.blockSignals(False)

    def pass_buttons_to_config_for_setup(self, value):
    #line up 'classic config' with buttons in gui.
        self.experiment.configuration['movie']['rot90'] = value
        self.experiment.configuration.save()
        self.experiment.files.movie.rot90 = value

        QMessageBox.warning(self, "Rotation change", "User alert: please first unselect all files, then delete all previously generated projection images [*ave*.tif] and datafiles [.nc]. See also 'Help'")
        # self.file.movie.rot90 = self.button_movie_rotation.value()

    def on_channel_selection_changed(self, button):
        id = self.channel_selector.id(button)
        if id == 1:
            for file in self.experiment.files:
                file.movie.channels = [Channel(file.movie, 'green', 'g', other_names=['donor', 'd'])]
                file.movie.channel_arrangement = [[[0, ]]]
        if id == 2:
            for file in self.experiment.files:
                file.movie.channels = [Channel(file.movie, 'green', 'g', other_names=['donor', 'd']),
                                       Channel(file.movie, 'red', 'r', other_names=['acceptor', 'a'])]
                file.movie.channel_arrangement = [[[0, 1]]]

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
