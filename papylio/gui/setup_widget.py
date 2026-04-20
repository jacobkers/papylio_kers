
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
    def __init__(self, top_tabs, parent=None):
        super(ExtractionWidget, self).__init__(parent)
        self.top_tabs = top_tabs
        self.parent = parent
        self.methods_spot_detection = {}
        self.method_forms_spot_detection = {}  # method_name -> (widget, inputs)
        self.parent.model.itemChanged.connect(self.onItemChange)

        #this one collects the various settings frames in horizontal fashion
        advanced_layout = QHBoxLayout()
        advanced_layout.setAlignment(Qt.AlignLeft)

        #1 general settings box ---------------------------------
        frame_general = Group_Box("General")
        advanced_general_layout = QFormLayout(frame_general)
        #1.1 channels:
        self.button_general_channel = QComboBox()
        self.button_general_channel.addItems(['donor', 'acceptor'])
        advanced_general_layout.addRow("channels", self.button_general_channel)
        #1.2 illuminations
        self.button_general_illumination = QSpinBox()
        self.button_general_illumination.setValue(0)
        advanced_general_layout.addRow("illumination:", self.button_general_illumination)
        #1.3 projection_type:
        self.button_general_projection_type = QComboBox()
        self.button_general_projection_type.addItems(['average', 'maximum'])
        advanced_general_layout.addRow("projection:", self.button_general_projection_type)
        #1.4 method:
        self.button_general_method = QComboBox()
        self.button_general_method.addItems(['by_channel', 'average_channels', 'sum_channels'])
        advanced_general_layout.addRow("channels:", self.button_general_method)

        #2. projection image box--------------------------------------
        #note that it is a subsection 'proj_image'
        frame_projection = Group_Box("Projection Image")
        advanced_projection_layout = QFormLayout(frame_projection)
        #2.1 type
        self.button_projection_image_type = QComboBox()
        self.button_projection_image_type.addItems(['average', 'maximum'])
        advanced_projection_layout.addRow("projection type:", self.button_projection_image_type)
        #2.2 frame range
        self.button_projection_image_frame_range = QLineEdit()
        self.button_projection_image_frame_range.setText("[0, 20]")
        advanced_projection_layout.addRow("frame range:", self.button_projection_image_frame_range)
        #2.3 illumination
        self.button_projection_image_illumination= QSpinBox()
        self.button_projection_image_illumination.setValue(0)
        advanced_projection_layout.addRow("illumination:", self.button_projection_image_illumination)

        #3. peakfinding (dynamic form building)
        # build a flexible form for the spot_detection channel:-------------------------
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

        #4. coordinate optimization box-------------------------------
        frame_coord_opt = Group_Box("Coordinate Optimization")
        advanced_coord_opt_layout = QFormLayout(frame_coord_opt)
        #4.1 margin
        self.button_coordinates_within_margin = QSpinBox()
        self.button_coordinates_within_margin.setValue(10)
        advanced_coord_opt_layout.addRow("within_margin:",
                                         self.button_coordinates_within_margin)
        #4.2 Gauss fit width
        self.button_coordinates_after_gaussian_fit_width = QSpinBox()
        self.button_coordinates_after_gaussian_fit_width.setValue(3)
        advanced_coord_opt_layout.addRow("gaussian_fit:" ,self.button_coordinates_after_gaussian_fit_width)


        #5. extract_traces box-------------------------------
        frame_extract_traces = Group_Box("Extract Traces")
        advanced_extract_traces_layout = QFormLayout(frame_extract_traces)
        #5.1 obsoletes/not used:

        #background_correction = None,
        #alpha_correction = None,
        #gamma_correction = None

        #5.2 mask
        self.button_extract_mask_size = QLineEdit()
        self.button_extract_mask_size.setPlaceholderText("11")
        self.button_extract_mask_size.setText("11")
        self.button_extract_mask_size.setToolTip("float number or presets: TIR-T, TIR-V, TIR-S 1.5x 2x2, TIR-S 1x 2x2, BN-TIRF")
        advanced_extract_traces_layout.addRow("mask size:", self.button_extract_mask_size)
        #5.3 neighbourhood_size: 11
        self.button_extract_neighbourhood_size= QSpinBox()
        self.button_extract_neighbourhood_size.setValue(11)
        advanced_extract_traces_layout.addRow("neighbourhood_size:", self.button_extract_neighbourhood_size)

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
        self.extract_button = make_push_button(
            'Extract traces',
            self.extract_traces,
            "extraction selected file(s)"
        )
        self.extract_button.clicked.connect(self.on_extract) #for tab switching

        self.find_coordinates_button= make_push_button(
            'Find coordinates',
            self.find_coordinates,"extraction selected file(s)"
        )
        self.find_coordinates_button.clicked.connect(self.on_find_coordinates)  # for tab switching


        extraction_controls = build_control_layouts([
            self.find_coordinates_button,
            self.extract_button,
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
