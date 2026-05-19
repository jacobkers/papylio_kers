#import sys
#import json
from PySide2.QtWidgets import QHBoxLayout,  \
    QPushButton, QTabWidget, QComboBox, QFormLayout, QWidget, QLabel, QVBoxLayout, QSpinBox
from PySide2.QtCore import Qt

#from papylio import File
from papylio.gui.common_layouts import (HelpDialog, Group_Box, get_button_value,
                                        build_control_layouts, make_push_button)
from papylio.plotting import wysiwyg_export
from matplotlib.figure import Figure
#import numpy as np

from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class HistogramWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent

        # imagery
        self.fig_histogram = Figure(figsize=(14, 3))
        self.histogram_canvas = FigureCanvas(self.fig_histogram)


        # main control buttons
        histogram_controls=build_control_layouts(
            [make_push_button('Export',self.export_histogram, "export plot contents"),
            ])

        #tune panel on the right

        histogram_panel_layout=QVBoxLayout()
        histogram_panel_layout.setAlignment(Qt.AlignRight)
        histogram_panel_layout.addWidget(histogram_controls)

        #tab layout
        histogram_tab_layout = QHBoxLayout()
        histogram_tab_layout.addWidget(self.histogram_canvas)
        histogram_tab_layout.addLayout(histogram_panel_layout)



        #other
        other_graph_layout = QHBoxLayout()
        #other_graph_layout.addWidget(box2)

        tabs = QTabWidget()
        tabs.setTabPosition(QTabWidget.North)
        tabs.setMovable(False)
        tabs.setDocumentMode(True)

        tab1 = QWidget(self)
        tab1.setLayout(histogram_tab_layout)
        tabs.addTab(tab1, 'Histogram')
        tab2 = QWidget(self)
        tab2.setLayout(other_graph_layout)
        tabs.addTab(tab2, 'Other')


        self.kinetics_widget = QWidget()
        kinetics_layout = QHBoxLayout()
        kinetics_layout.addWidget(tabs)
        #self.kinetics_widget.setLayout(kinetics_layout)
        self.setLayout(kinetics_layout)

        self.file = None

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        self._file = file
        if file is None:
            self.setDisabled(True)
        else:
            self.setDisabled(False)

    def export_histogram(self):
        fpath=self.parent.experiment.analysis_path / 'Histograms'
        wysiwyg_export(self.fig_histogram, filepath=fpath,  filename="histogram_export", filetype="csv")

    def refresh_histogram(self):
        self.fig_histogram.clear()
        axes = self.fig_histogram.subplots(1, 2, sharex=True)

        selected_files = self.parent.experiment.selectedFiles

        self.histogram_canvas.draw_idle()

        if selected_files:
            #TODO: here, perform the histogram function (possibly both 1D and 2D)
            test_file=selected_files[0]
            bins = [i * 0.01 for i in range(200)]
            test_file.show_histogram(variable='FRET', axis=axes[0], bins=bins)
            #---UNDER CONSTRUCTION-------
            #selected_files.serial.determine_dwells_from_classification(variable=variable, selected=True, inactivate_start_and_end_states=True)
            #selected_files.serial.analyze_dwells(method=method, number_of_exponentials=[1, N_exponents])
            #selected_files.serial.plot_histogram_analysis(plot_range=(0, N_exponents), axes=axes, log=False)



