
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QComboBox, QLineEdit, QLabel
from PySide2.QtCore import Qt

import matplotlib.pyplot as plt

from matplotlib.figure import Figure

from papylio.gui.advanced_widget import Expander, ImageCanvas
import matchpoint as mp

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)


class MappingWidget(QWidget):
    def __init__(self, parent=None):
        super(MappingWidget, self).__init__(parent)
        self.parent = parent

        # imagery
        self.fig1 = Figure(figsize=(4, 3))
        self.map_overlay_image_canvas = FigureCanvas(self.fig1)
        map_overlay_image_layout = QVBoxLayout()
        map_overlay_image_layout.addWidget(self.map_overlay_image_canvas)

        self.map_overlay_image = QWidget()
        self.map_overlay_image.setLayout(map_overlay_image_layout)
        self.map_overlay_image.setToolTip("overlay image shows result of mapping")

        # imagery
        self.map_image_canvas = ImageCanvas(self, width=4, height=3, dpi=100)
        self.map_image_canvas.setToolTip("shows raw image (1st of multiple selection)")
        map_image_toolbar = NavigationToolbar(self.map_image_canvas, self)
        map_image_layout = QVBoxLayout()
        map_image_layout.addWidget(map_image_toolbar)
        map_image_layout.addWidget(self.map_image_canvas)

        # Create a placeholder widget to hold our toolbar and canvas.
        self.map_image = QWidget()
        self.map_image.setLayout(map_image_layout)


        #map settings
        #mapping buttons definition:
        self.button_map_label = QLabel("method:")
        self.button_map_method_combobox = QComboBox()
        self.button_map_options = ['icp', 'nn', 'foo']
        self.button_map_method_combobox.addItems(self.button_map_options)


        button_map_dist_treshold = QLineEdit()
        button_map_dist_treshold.setPlaceholderText("dist_treshold")
        button_map_margin = QLineEdit()
        button_map_margin.setPlaceholderText("edge_margin")

        #donor mapping:
        button_map_donor_label = QLabel("donor_pks")
        button_map_donor_method_combobox = QComboBox()
        button_map_donor_options = ['local-maximum-auto', 'other']
        button_map_donor_method_combobox.addItems(button_map_donor_options)
        button_map_donor_fract_diff = QLineEdit()
        button_map_donor_fract_diff.setPlaceholderText("fract._diff")
        button_map_donor_ns_min = QLineEdit()
        button_map_donor_ns_min.setPlaceholderText("nbh_min")
        button_map_donor_ns_max = QLineEdit()
        button_map_donor_ns_max.setPlaceholderText("nbh_max")
        # acceptor mapping:
        button_map_acceptor_label = QLabel("acceptor_pks")
        button_map_acceptor_method_combobox = QComboBox()
        button_map_acceptor_options = ['local-maximum-auto', 'other']
        button_map_acceptor_method_combobox.addItems(button_map_acceptor_options)
        button_map_acceptor_fract_diff = QLineEdit()
        button_map_acceptor_fract_diff.setPlaceholderText("fraction_diff")
        button_map_acceptor_ns_min = QLineEdit()
        button_map_acceptor_ns_min.setPlaceholderText("filt_nbh_min")
        button_map_acceptor_ns_max = QLineEdit()
        button_map_acceptor_ns_max.setPlaceholderText("filt_nbh_max")
        #mapping grid layout:
        map_buttons_layout = QGridLayout()
        map_buttons_layout.setAlignment(Qt.AlignLeft)
        map_buttons_layout.addWidget(self.button_map_label, 0, 0)
        map_buttons_layout.addWidget(self.button_map_method_combobox, 0, 1)


        #add below buttons under 'advanced':---------------------------------------------------
        advanced_layout = QGridLayout()
        advanced_layout.addWidget(button_map_dist_treshold, 1, 0)
        advanced_layout.addWidget(button_map_margin, 1, 1)

        advanced_layout.addWidget(button_map_donor_label, 2, 0)
        advanced_layout.addWidget(button_map_donor_method_combobox, 2, 1)
        advanced_layout.addWidget(button_map_donor_fract_diff, 2, 2)
        advanced_layout.addWidget(button_map_donor_ns_min, 2, 3)
        advanced_layout.addWidget(button_map_donor_ns_max, 2, 4)

        advanced_layout.addWidget(button_map_acceptor_label, 3, 0)
        advanced_layout.addWidget(button_map_acceptor_method_combobox, 3, 1)
        advanced_layout.addWidget(button_map_acceptor_fract_diff, 3, 2)
        advanced_layout.addWidget(button_map_acceptor_ns_min, 3, 3)
        advanced_layout.addWidget(button_map_acceptor_ns_max, 3, 4)


        #build panel layout:
        map_buttons = QWidget()
        map_buttons.setLayout(map_buttons_layout)

        map_advanced = Expander("Advanced")
        map_advanced.setContentLayout(advanced_layout)

        #extra:
        perform_mapping_button = QPushButton('Map it')
        perform_mapping_button.setToolTip("Map selected file(s) using above settings")
        perform_mapping_button.clicked.connect(self.perform_mapping)

        #collect:
        map_controls_layout = QGridLayout()
        map_controls_layout.setAlignment(Qt.AlignTop)
        map_controls_layout.addWidget(map_buttons,0,0)
        map_controls_layout.addWidget(map_advanced,1,0)  # test!
        map_controls_layout.addWidget(perform_mapping_button,3,0)


        #pack in widget:
        self.map_controls = QWidget()
        self.map_controls.setLayout(map_controls_layout)
        self.map_controls.setMinimumWidth(150)

        #add all to tab:
        mapping_tab_layout = QHBoxLayout()
        mapping_tab_layout.addWidget(self.map_controls)
        mapping_tab_layout.addWidget(self.map_image)
        mapping_tab_layout.addWidget(self.map_overlay_image)


        self.setLayout(mapping_tab_layout)

    def perform_mapping(self, t):
        print(t)
        self.fig1.clear()
        ax1 = self.fig1.add_subplot(111)

        selected_files = self.parent.experiment.selectedFiles

        #jk-read in a default config and change one value:
        panel_config= self.parent.experiment.configuration['mapping']
        panel_method = self.button_map_method_combobox.currentText()
        panel_config['method']=panel_method

        plot_file = selected_files[0]
        plot_file.mapping.show_mapping_transformation(axis=ax1)

        self.map_overlay_image_canvas.draw_idle()

        if selected_files:
            selected_files.serial.perform_mapping(**panel_config)
            self.map_image_canvas.refresh()




