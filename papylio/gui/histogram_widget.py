import sys
import json
from PySide2.QtWidgets import QHBoxLayout,  \
    QPushButton, QTabWidget, QComboBox, QSizePolicy
from PySide2.QtWidgets import QWidget, QLabel, QVBoxLayout
from PySide2.QtCore import Qt

from papylio import File
from papylio.gui.common_layouts import (HelpDialog,
                                        build_control_layouts, make_push_button)
from papylio.plotting import wysiwyg_export
from matplotlib.figure import Figure
import numpy as np
import matplotlib as mpl
from matplotlib.backends.backend_qtagg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class HistogramWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent

        self.fig_histogram = Figure(figsize=(14, 3))
        self.histogram_canvas = HistogramCanvas(self.fig_histogram)

        # Create toolbar, passing canvas as first parameter, parent (self, the MainWindow) as second.
        #histogram_layout = QVBoxLayout()
        #histogram_layout.addWidget(self.histogram_canvas)

        # main control buttons
        histogram_controls = build_control_layouts(
            [make_push_button('Export', self.export_histogram, "export plot contents"),
             ])

        # tune panel on the right

        histogram_panel_layout = QVBoxLayout()
        histogram_panel_layout.setAlignment(Qt.AlignRight)
        histogram_panel_layout.addWidget(histogram_controls)

        # tab layout
        histogram_tab_layout = QHBoxLayout()
        histogram_tab_layout.addWidget(self.histogram_canvas)
        #histogram_tab_layout.addLayout(histogram_panel_layout)

        # other
        other_graph_layout = QHBoxLayout()
        # other_graph_layout.addWidget(box2)

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
        # self.kinetics_widget.setLayout(kinetics_layout)
        self.setLayout(kinetics_layout)

        self.file = None

        # Create a placeholder widget to hold our toolbar and canvas.
        #self.setLayout(histogram_layout)

    @property
    def file(self):
        return self.image_canvas.file

    @file.setter
    def file(self, file):
        self.histogram_canvas.file = file
        if file is None:
            self.setDisabled(True)
        else:
            self.setDisabled(False)

    def export_histogram(self):
        fpath=self.parent.experiment.analysis_path / 'Histograms'
        wysiwyg_export(self.fig_histogram, filepath=fpath,  filename="histogram_export", filetype="csv")

class HistogramCanvas(FigureCanvas):
    def __init__(self, parent=None, width=14, height=7, dpi=100):
        self.fig_histogram= mpl.figure.Figure(figsize=(width, height), dpi=dpi,
                                        constrained_layout=True)  # , figsize=(2, 2))
        super().__init__(self.fig_histogram)
        self.parent = parent
        self._file = None

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        if file is not None and file is not self._file:
            self._file = file
            self.refresh()
        elif file is None:
            self._file = None
            self.figure.clf()
            self.draw()

    def refresh(self):
        #self.figure.clf()
        #TODO here the histogram should be built and shown
        self.fig_histogram.clear()
        axes = self.fig_histogram.subplots(1, 2, sharex=True)
        #selected_files = self.parent.experiment.selectedFiles
        #self.fig_histogram.clf()

        #if selected_files:
        # TODO: here, perform the histogram function (possibly both 1D and 2D)
        #test_file = selected_files[0]
        bins = [i * 0.01 for i in range(200)]
        self.file.show_histogram(variable='FRET', axis=axes[0], bins=bins)

        self.file.histogram_2D_intensity_per_channel(
            ax=axes[1],
            bins=50,
            show_marginal=False
        )

        self.draw()