import sys
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QTabWidget, QTableWidget, QComboBox, QLineEdit
from PySide2.QtGui import QStandardItem, QStandardItemModel
from PySide2.QtCore import Qt
import numpy as np

from matplotlib.backends.backend_qtagg import (
    NavigationToolbar2QT as NavigationToolbar)


class AlignmentWidget(QWidget):
    def __init__(self, parent=None, *args, **kwargs):
        super(AlignmentWidget, self).__init__(parent)

        #containing an image display (including a toolbar)  and action buttons.
        # a. Create toolbar, passing canvas as first parameter, parent (self, the MainWindow) as second.
        #image_toolbar = NavigationToolbar(self.image_canvas, self)
        #image_layout = QVBoxLayout()
        #image_layout.addWidget(image_toolbar)
        #image_layout.addWidget(self.image_canvas)

        # b. Create a placeholder widget to hold our toolbar and the canvas as defined above.
        #toolbar is the column of control buttons to the right
        self.image = QWidget()
        #self.image.setLayout(image_layout)

        controls_layout = QGridLayout()
        controls_layout.setAlignment(Qt.AlignTop)

        # c. The  main action buttons are defined here (for now, the are grouped together:
        # map:
        perform_mapping_button = QPushButton('go: (err)')
        perform_mapping_button.setToolTip("select bead slide in tree pane and press; afterward close pop-ups")
        if 0: perform_mapping_button.clicked.connect(self.perform_mapping)
        controls_layout.addWidget(perform_mapping_button, 1, 0, 1, 2)

        self.controls = QWidget()
        self.controls.setLayout(controls_layout)
        self.controls.setMinimumWidth(200)

        movie_layout = QHBoxLayout()
        movie_layout.addWidget(self.image)
        movie_layout.addWidget(self.controls)

        self.setLayout(movie_layout)

        self._experiment = None

    @property
    def experiment(self):
        return self._experiment

    @experiment.setter
    def experiment(self, experiment):
        self._experiment = experiment