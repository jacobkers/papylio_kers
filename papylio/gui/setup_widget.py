
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, \
    QComboBox, QLineEdit, QSpinBox, QFormLayout
from PySide2.QtCore import Qt, Signal
from papylio import File
from papylio.gui.common_layouts import (Expander, ImageCanvas, HelpDialog,Group_Box,
                                        build_control_layouts,make_push_button,
                                        build_form,build_parameters_input, get_button_value)
#for registry:
from papylio.peak_finding import (find_peaks_absolute_threshold,
                                  find_peaks_adaptive_threshold,
                                  find_peaks_local_maximum,
                                  find_peaks_local_maximum_auto,
                                  find_peaks_relative_local_maximum)

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class SetUpWidget(QWidget):
    def __init__(self, parent=None):
        super(SetUpWidget, self).__init__(parent)

        # 1 movie settings box ---------------------------------
        frame_movie = Group_Box(title="Movie", highlight=True)
        frame_movie.setToolTip('setting these requires reloading: press refresh')
        movie_setup_layout = QFormLayout(frame_movie)
        # 1.1 channels:
        self.button_movie_rotation = QSpinBox()
        self.button_movie_rotation.setValue(1)
        self.button_movie_rotation.setToolTip('press refresh to have effect')
        movie_setup_layout.addRow("movie rotation", self.button_movie_rotation)

        advanced_layout = QHBoxLayout()
        advanced_layout.setAlignment(Qt.AlignRight)
        advanced_layout.addWidget(frame_movie)

        # build panel layout:
        setup_advanced = Expander("Advanced")
        setup_advanced.setContentLayout(advanced_layout)

        self.parent = parent
        start_help_button = build_control_layouts([
            make_push_button('Read Me', self.show_main_help, None)])

        start_tab_layout = QVBoxLayout()
        #TODO: development: keep invisible as long as it doesn't function:
        #start_tab_layout.addWidget(setup_advanced)
        start_tab_layout.addWidget(start_help_button)
        self.setLayout(start_tab_layout)

    def pass_buttons_to_config_for_setup(self):
    #line up 'classic config' with buttons in gui.
        self.parent.experiment.configuration['movie']['rot90']=get_button_value(self.button_movie_rotation)

    def show_main_help(self):
        help_text = """
                <html>
                  <body style="font-family: sans-serif; font-size: 10pt;">

                    <h2>Welcome</h2>

                    <p>
                      This gui is based on the Papylio framework
                    </p>

                    <ul>
                      <li>Select and view movies and variables in the top panel</li>
                      <li>Walk the pipeline via the tabs in the bottom panel</li>
                    </ul>

                    <p>
                      For background, see
                      <a href="https://papylio.readthedocs.io/en/stable/user_guide/index.html">
                        Papylio documentation
                      </a>.
                    </p>

                    <h3>Tips</h3>

                    <p>
                      <ul>
                        <li>Hover over buttons for help notes.</li>
                        <li>Find more detailed info under the 'Help' buttons per tab</li>
                    </ul>

                    </p>
                    
                    <h3>Known Issues</h3>

                    <p>
                    For known issues in Papylio, see
                      <a href="https://github.com/Chirlmin-Joo-lab/papylio/issues">
                        Papylio Issues
                      </a>.

                    </p>

                  </body>
                </html>
                """
        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()
