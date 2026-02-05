import sys
import json
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QTabWidget, QTableWidget, QComboBox, QLineEdit
from PySide2.QtGui import QStandardItem, QStandardItemModel
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QWidget, QLabel, QVBoxLayout
from PySide2.QtGui import QPixmap
from PySide2.QtCore import Qt
from papylio.gui.common_layouts import ImageCanvas,Expander,HelpDialog

import numpy as np

class KineticsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        #build fake graphs:------------------------
        graphs_path = r'C:\Users\jkerssemakers\OneDrive - Delft University of Technology\ChJ_recent\ChJ25_Jacob\2026_01_17 GUI stubs'
        #graphs_path = r'M:\tnw\bn\alg\Shared\Jacob\ChJ_lab\2026_01_17 GUI stubs'
        graph1_path = graphs_path + r'\Renee_on_PapylioTest_01.jpeg'
        graph2_path = graphs_path + r'\Renee_on_PapylioTest_02.jpeg'
        graph3_path = graphs_path + r'\Renee_on_PapylioTest_03.jpeg'
        graph4_path = graphs_path + r'\dwell_time analysis.jpg'
        box1 = QLabel()
        pixmap = QPixmap(graph4_path)
        box1.setPixmap(pixmap)
        box1.setAlignment(Qt.AlignCenter)
        box1.setPixmap(
            pixmap.scaled(
                box1.size(),
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )
        )

        box2 = QLabel()
        pixmap = QPixmap(graph2_path)
        box2.setPixmap(pixmap)
        box2.setAlignment(Qt.AlignCenter)
        box2.setPixmap(
            pixmap.scaled(
                box2.size(),
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )
        )




        #main
        dwell_action_button = QPushButton('Get Dwells')
        dwell_action_button.clicked.connect(self.perform_mapping)
        dwell_help_button = QPushButton('Help!')
        dwell_help_button.clicked.connect(self.show_dwell_help)

        dwell_controls_layout = QVBoxLayout()
        dwell_controls_layout.addWidget(dwell_action_button)
        dwell_controls_layout.addWidget(dwell_help_button)

        dwell_controls = QWidget()
        dwell_controls.setLayout(dwell_controls_layout)

        dwell_times_tab_layout = QHBoxLayout()
        dwell_times_tab_layout.addWidget(dwell_controls)
        dwell_times_tab_layout.addWidget(box1)

        #other
        other_graph_layout = QHBoxLayout()
        other_graph_layout.addWidget(box2)


        tabs = QTabWidget()
        tabs.setTabPosition(QTabWidget.North)
        tabs.setMovable(False)
        tabs.setDocumentMode(True)

        tab1 = QWidget(self)
        tab1.setLayout(dwell_times_tab_layout)
        tabs.addTab(tab1, 'Dwell Time')
        tab2 = QWidget(self)
        tab2.setLayout(other_graph_layout)
        tabs.addTab(tab2, 'Other')


        self.kinetics_widget = QWidget()
        kinetics_layout = QHBoxLayout()
        kinetics_layout.addWidget(tabs)
        #self.kinetics_widget.setLayout(kinetics_layout)
        self.setLayout(kinetics_layout)

    def show_dwell_help(self):
        help_text = help_text = """\
                Dwell Times

                -----
                • description
                • ref to docs
                
                example 
                -------------------------------------
                
                settings examples
                ----------------
                """

        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()


