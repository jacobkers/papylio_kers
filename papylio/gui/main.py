import sys
import platform

from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QTabWidget, QTableWidget, QComboBox, QLineEdit
from PySide2.QtGui import QStandardItem, QStandardItemModel, QIcon
from PySide2.QtCore import Qt
import matplotlib as mpl
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt

import papylio as pp

from papylio import Experiment, File
from papylio.trace_plot import TracePlotWindow
from papylio.gui.selection_widget import SelectionWidget

class MainWindow(QMainWindow):
    def __init__(self, main_path=None):
        super().__init__()

        system = platform.system()
        if system == "Windows":
            extension = 'ico'
        else:  # macOS or Linux
            extension = "png"
        self.setWindowIcon(QIcon("icon."+extension))
        self.setWindowTitle("Papylio v" + pp.__version__ )

        #I. set up the left section of the GUI,
        # displaying a tree of selectable files in th dataset-directory:
        self.tree = QTreeView(self)
        self.model = QStandardItemModel()
        self.root = self.model.invisibleRootItem()
        self.model.setHorizontalHeaderLabels(['Name', 'Count'])
        self.tree.header().setDefaultSectionSize(180)
        self.tree.setModel(self.model)

        self.tree.setFocusPolicy(Qt.NoFocus)
        self.tree.setFixedWidth(256)
        self.update = True
        self.model.itemChanged.connect(self.onItemChange)
        self.image_canvas = ImageCanvas(self, width=5, height=4, dpi=100)

        #II. Build the 'extraction' tab, containing an image display (including a toolbar)  and action buttons.
        # a. Create toolbar, passing canvas as first parameter, parent (self, the MainWindow) as second.
        image_toolbar = NavigationToolbar(self.image_canvas, self)
        image_layout = QVBoxLayout()
        image_layout.addWidget(image_toolbar)
        image_layout.addWidget(self.image_canvas)

        # b. Create a placeholder widget to hold our toolbar and the canvas as defined above.
        #toolbar is the column of control buttons to the right
        self.image = QWidget()
        self.image.setLayout(image_layout)
        controls_layout = QGridLayout()
        controls_layout.setAlignment(Qt.AlignTop)

        # c. The  main action buttons are defined here:
        perform_mapping_button = QPushButton('Perform mapping')
        perform_mapping_button.clicked.connect(self.perform_mapping)
        controls_layout.addWidget(perform_mapping_button, 1, 0, 1, 2)
        find_molecules_button = QPushButton('Find coordinates')
        find_molecules_button.clicked.connect(self.find_coordinates)
        controls_layout.addWidget(find_molecules_button, 2, 0, 1, 2)
        extract_traces_button = QPushButton('Extract traces')
        extract_traces_button.clicked.connect(self.extract_traces)
        controls_layout.addWidget(extract_traces_button, 3, 0, 1, 2)

        self.controls = QWidget()
        self.controls.setLayout(controls_layout)
        self.controls.setMinimumWidth(200)

        extraction_layout = QHBoxLayout()
        extraction_layout.addWidget(self.image)
        extraction_layout.addWidget(self.controls)

        tabs = QTabWidget()
        tabs.setTabPosition(QTabWidget.North)
        tabs.setMovable(False)
        tabs.setDocumentMode(True)

        tab1 = QWidget(self)
        tab1.setLayout(extraction_layout)

        # III. Build a tab for trace evaluation by eye:
        tabs.addTab(tab1, 'Movie')
        self.traces = TracePlotWindow(parent=self, width=4, height=3, show=False)

        tabs.addTab(self.traces, 'Traces')
        self.selection = SelectionWidget()

        # IV. Build a tab for trace evaluation by tresholds and more:
        tabs.addTab(self.selection, 'Selection (beta)')
        tabs.currentChanged.connect(self.setTabFocus)

        # V. build an 'experiment' pane with tree and refresh button in vertical order:
        refresh_button = QPushButton('Refresh')
        refresh_button.clicked.connect(self.refresh)
        experiment_layout = QVBoxLayout()
        experiment_layout.addWidget(refresh_button)
        experiment_layout.addWidget(self.tree)

        # VI. now, assemble the various panes
        layout = QHBoxLayout()
        layout.addLayout(experiment_layout)
        layout.addWidget(tabs)

        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

        self.show()

        # self.experiment = Experiment(
        #     r'D:\SURFdrive\Promotie\Code\Python\papylio\twoColourExampleData\20141017 - Holliday junction - Copy')
        self.experiment = pp.Experiment(main_path, main_window=self)
        self.addExperiment(self.experiment)
        self.traces.save_path = self.experiment.analysis_path.joinpath('Trace_plots')

    def keyPressEvent(self, e):
        self.traces.keyPressEvent(e)

    def setTabFocus(self, e):
        if e == 0:
            self.image.setFocus()
        if e == 1:
            self.traces.setFocus()

    def midChange(self, input):
        input = int(input)
        self.experiment.configuration['find_coordinates']['peak_finding']['minimum_intensity_difference'] = input
        self.experiment.configuration.save()

    def perform_mapping(self, t):
        print(t)
        selected_files = self.experiment.selectedFiles
        if selected_files:
            selected_files.serial.perform_mapping()
            self.image_canvas.refresh()
            plt.show()

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

    def update_plots(self):
        selected_files = self.experiment.selectedFiles + [None]
        self.image_canvas.file = selected_files[0]
        if selected_files[0] is not None:
            self.traces.dataset = selected_files[0].dataset
            self.selection.file = selected_files[0]
        else:
            self.traces.dataset = None
            self.selection.file = None

    def addExperiment(self, experiment):
        #jk_note: when uncommenting a path here, code basically auto-loads:
        # experiment = Experiment(r'C:\Users\myname\personalpaths\Papylio example dataset')
        self.root.appendRow([
        QStandardItem(experiment.name),
        QStandardItem(0),
        ])
        experimentNode = self.root.child(self.root.rowCount() - 1)
        for file in experiment.files:
            print('addfile'+file.name)
            self.addFile(file, experimentNode)

        self.tree.expandAll()

        print('add')

    def addFile(self, file, experimentNode):
        folders = file.relativePath.parts

        parentItem = experimentNode
        parentItem.setCheckable(True)
        for folder in folders:

            # Get the folderItems and folder names for the current folderItem
            nodeItems = [parentItem.child(i) for i in range(parentItem.rowCount())]# if item.type == 'folder']
            nodeItemNames = [item.text() for item in nodeItems]

            if folder not in nodeItemNames:
                # Add new item for the folder and set parentItem to this item
                parentItem.appendRow([
                    QStandardItem(folder),
                    QStandardItem(0),
                ])
                parentItem = parentItem.child(parentItem.rowCount() - 1)
                parentItem.setCheckable(True)
            else:
                # Set parent item to the found folderItem
                parentItem = nodeItems[nodeItemNames.index(folder)]

        parentItem.appendRow([
            QStandardItem(file.name),
            QStandardItem(0),
        ])
        item = parentItem.child(parentItem.rowCount() - 1)
        item.setCheckable(True)
        if file.isSelected:
            item.setCheckState(Qt.Checked)
        else:
            item.setCheckState(Qt.Unchecked)
        item.setData(file)
        #self.FileItems.append(item)

        # self.insertDataIntoColumns(item)

        return item

    def refresh(self):
        self.root.removeRows(0, 1)
        self.experiment = Experiment(self.experiment.main_path)
        self.addExperiment(self.experiment)


class ImageCanvas(FigureCanvas):
    def __init__(self, parent=None, width=14, height=7, dpi=100):
        self.figure = mpl.figure.Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)  # , figsize=(2, 2))
        super().__init__(self.figure)
        self.parent = parent

        # self.axis = self.figure.gca()

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
        self.figure.clf()
        self._file.movie.determine_spatial_background_correction(use_existing=True)
        self._file.show_coordinates_in_image(figure=self.figure)
        self.draw()


if __name__ == '__main__':
    from multiprocessing import Process, freeze_support
    freeze_support()

    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec_()
