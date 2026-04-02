import sys
import PySide2
import platform

import sys
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QTabWidget, QTableWidget, QComboBox, QLineEdit, QLabel
from PySide2.QtGui import QStandardItem, QStandardItemModel, QIcon
from PySide2.QtCore import Qt
import matplotlib as mpl

from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

import matplotlib.pyplot as plt
import numpy as np

import papylio as pp
from papylio import Experiment, File
from papylio.trace_plot import TracePlotWindow
from papylio.gui.selection_widget import SelectionWidget
from papylio.gui.classification_widget import ClassificationWidget
from papylio.gui.mapping_widget import MappingWidget
from papylio.gui.extraction_widget import ExtractionWidget
from papylio.gui.kinetics_widget import KineticsWidget
from papylio.gui.common_layouts import ImageCanvas, HelpDialog

class MainWindow(QMainWindow):
    # def __init__(self):
    #     super().__init__()
    #     # model = QFileSystemModel()
    #     # model.setRootPath(QDir.currentPath())
    #
    #
    #     self.model = TreeModel()
    #
    #     self.tree = QTreeView()
    #     self.tree.setModel(self.model)
    #
    #      #experiment = Experiment(r'C:\Users\ivoseverins\surfdrive\Promotie\Code\Python\papylio\twoColourExampleData\20141017 - Holliday junction - Copy')
    #     #self.model.addExperiment(experiment)
    #
    #     self.setCentralWidget(self.tree)

    def __init__(self, main_path=None):
        super().__init__()

        system = platform.system()
        if system == "Windows":
            extension = 'ico'
        else:  # macOS or Linux
            extension = "png"
        self.setWindowIcon(QIcon("icon." + extension))
        self.setWindowTitle("Papylio v" + pp.__version__)

        self.tree = QTreeView(self)
        layout = QVBoxLayout()
        layout.addWidget(self.tree)
        self.model = QStandardItemModel()
        self.root = self.model.invisibleRootItem()
        self.model.setHorizontalHeaderLabels(['Name', 'Count'])
        self.tree.header().setDefaultSectionSize(180)
        self.tree.setModel(self.model)

        self.tree.setFocusPolicy(Qt.NoFocus)
        self.tree.setFixedWidth(256)
        self.update = True

        self.model.itemChanged.connect(self.onItemChange)

        # imagery (currently goes to extraction tab)
        self.image_canvas = ImageCanvas(self, width=4, height=4, dpi=100)

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        image_toolbar = NavigationToolbar(self.image_canvas, self)
        image_layout = QVBoxLayout()
        image_layout.addWidget(image_toolbar)
        image_layout.addWidget(self.image_canvas)

        # Create a placeholder widget to hold our toolbar and canvas.
        self.image = QWidget()
        self.image.setLayout(image_layout)

        # main buttons:
        main_help_button = QPushButton('Read me')
        main_help_button.clicked.connect(self.show_main_help)

        start_tab_layout = QVBoxLayout()
        start_tab_layout.addWidget(main_help_button)

        tabs = QTabWidget()
        tabs.setTabPosition(QTabWidget.North)
        tabs.setMovable(True)
        tabs.setDocumentMode(True)

        tab0 = QWidget(self)
        tab0.setLayout(start_tab_layout)
        tabs.addTab(tab0, 'Start')
        self.mapping = MappingWidget(parent=self)
        tabs.addTab(self.mapping, 'Mapping')
        self.extraction = ExtractionWidget(parent=self)
        tabs.addTab(self.extraction, 'Extraction')
        self.selection = SelectionWidget(parent=self)
        tabs.addTab(self.selection, 'Selection')
        self.classification = ClassificationWidget(parent=self)
        self.classification.classificationChanged.connect(self.update_plots)
        tabs.addTab(self.classification, 'Classification')
        self.kinetics = KineticsWidget(parent=self)
        tabs.addTab(self.kinetics, 'Kinetics')
        tabs.currentChanged.connect(self.setTabFocus)

        # refresh & tree
        experiment_layout = QVBoxLayout()
        refresh_button = QPushButton('Refresh')
        refresh_button.clicked.connect(self.refresh)
        experiment_layout.addWidget(refresh_button)
        experiment_layout.addWidget(self.tree)

        # right side has a viewing pane (top) and
        # viewing panes
        top_tabs = QTabWidget()
        top_tabs.setTabPosition(QTabWidget.North)
        top_tabs.setMovable(False)
        top_tabs.setDocumentMode(True)
        self.traces = TracePlotWindow(parent=self, width=4, height=5, show=False)
        top_tabs.addTab(self.traces, 'Traces')
        top_tabs.addTab(self.image, 'Frame')

        top_layout = QVBoxLayout()
        top_layout.addWidget(top_tabs)
        # ... a bottom pane (pipeline)
        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(tabs)

        # build main panel
        # full left is file tree
        left_layout = QVBoxLayout()
        left_layout.addLayout(experiment_layout)

        right_layout = QVBoxLayout()
        right_layout.addLayout(top_layout)
        right_layout.addLayout(bottom_layout)

        super_layout = QHBoxLayout()
        super_layout.addLayout(left_layout)
        super_layout.addLayout(right_layout)

        widget = QWidget()
        widget.setLayout(super_layout)

        self.setCentralWidget(widget)
        self.show()
        self.showMaximized()

        self.experiment = Experiment(
            r'C:\Users\jkerssemakers\OneDrive - Delft University of Technology\Documents\GitHub\Papylio example dataset')
        #self.experiment = pp.Experiment(main_path, main_window=self)
        self.addExperiment(self.experiment)
        self.traces.save_path = self.experiment.analysis_path.joinpath('Trace_plots')

    def keyPressEvent(self, e):
        self.traces.keyPressEvent(e)

    def setTabFocus(self, e):
        if e == 0:
            self.image.setFocus()
        if e == 1:
            self.traces.setFocus()

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
            self.classification.file = selected_files[0]
        else:
            self.traces.dataset = None
            self.selection.file = None
            self.classification.file = None

    def addExperiment(self, experiment):

        # experiment = Experiment(r'D:\SURFdrive\Promotie\Code\Python\papylio\twoColourExampleData\20141017 - Holliday junction - Copy')
        # experiment = Experiment(r'C:\Users\ivoseverins\surfdrive\Promotie\Code\Python\papylio\twoColourExampleData\20141017 - Holliday junction - Copy')
        self.root.appendRow([
            QStandardItem(experiment.name),
            QStandardItem(0),
        ])
        experimentNode = self.root.child(self.root.rowCount() - 1)
        for file in experiment.files:
            print('addfile' + file.name)
            self.addFile(file, experimentNode)

        self.tree.expandAll()

        print('add')

    def addFile(self, file, experimentNode):
        folders = file.relativePath.parts

        parentItem = experimentNode
        parentItem.setCheckable(True)
        for folder in folders:

            # Get the folderItems and folder names for the current folderItem
            nodeItems = [parentItem.child(i) for i in range(parentItem.rowCount())]  # if item.type == 'folder']
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
        # self.FileItems.append(item)

        # self.insertDataIntoColumns(item)

        return item

    def refresh(self):
        self.root.removeRows(0, 1)
        self.experiment = Experiment(self.experiment.main_path)
        self.addExperiment(self.experiment)

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
                      For background, see the
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

                  </body>
                </html>
                """
        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()


if __name__ == '__main__':
    from multiprocessing import Process, freeze_support

    freeze_support()

    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec_()


