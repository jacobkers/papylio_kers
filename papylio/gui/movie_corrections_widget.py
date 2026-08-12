
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, \
    QComboBox, QLineEdit, QSpinBox, QFormLayout, QButtonGroup, QRadioButton, QLabel, QMessageBox
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

        # 1 movie corrections box 1: temporal ---------------------------------
        frame_temporal_correction = Group_Box(title="Background treatment", highlight=False)
        frame_temporal_correction.setToolTip('"see help')


        temporal_correction_layout = QFormLayout(frame_temporal_correction)

        # todo change buttons etc
        # 1.1 channels:
        self.button_movie_rotation = QSpinBox()
        self.button_movie_rotation.setRange(-3, 3)
        self.button_movie_rotation.setValue(0)
        self.button_movie_rotation.valueChanged.connect(self.pass_buttons_to_config_for_setup)
        self.button_movie_rotation.setToolTip('If changed, remove all analyzed data and projection images')
        temporal_correction_layout.addRow("Rotation", self.button_movie_rotation)


        self.number_of_channels_selector_1 = QRadioButton("1")
        self.number_of_channels_selector_1.setToolTip('This needs to be manually set every time on reopening the gui')

        self.number_of_channels_selector_2 = QRadioButton("2")
        self.number_of_channels_selector_2.setChecked(True)
        self.number_of_channels_selector_2.setToolTip('This needs to be manually set every time on reopening the gui')

        self.channel_selector = QButtonGroup()
        self.channel_selector.addButton(self.number_of_channels_selector_1, 1)
        self.channel_selector.addButton(self.number_of_channels_selector_2, 2)
        # self.channel_selector.setToolTip("Choose number of channels")
        self.channel_selector.buttonClicked.connect(self.on_channel_selection_changed)

        channel_layout = QHBoxLayout()
        channel_layout.setContentsMargins(0, 0, 0, 0)
        channel_layout.addWidget(self.number_of_channels_selector_1)
        channel_layout.addWidget(self.number_of_channels_selector_2)

        temporal_correction_layout.addRow("Number of channels", channel_layout)


        advanced_layout = QHBoxLayout()
        advanced_layout.setAlignment(Qt.AlignLeft)
        advanced_layout.addWidget(frame_temporal_correction)

        start_help_button = build_control_layouts([
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
                    There are three types of background subtraction, to be performed in this order: 
                    </p>
                    
                    <p>
                    <ol>
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
