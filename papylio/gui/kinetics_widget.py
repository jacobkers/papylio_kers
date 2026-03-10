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

from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

class KineticsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent

        # imagery
        self.fig_kinetics = Figure(figsize=(5, 3))
        self.dwell_kinetics_canvas = FigureCanvas(self.fig_kinetics)

        dwell_variable_label=QLabel("variable:")
        dwell_variable_combobox = QComboBox()
        dwell_variable_items = ['FRET', 'intensity', '.....']
        dwell_variable_combobox.addItems( dwell_variable_items)
        dwell_options_layout=QHBoxLayout()
        dwell_options_layout.addWidget(dwell_variable_label)
        dwell_options_layout.addWidget(dwell_variable_combobox)

        dwell_options = QWidget()
        dwell_options.setLayout(dwell_options_layout)

        # main control buttons
        dwell_controls=build_control_layouts(
            [make_push_button("Get Dwells", self.perform_dwell_times_sequence,"fit and show dwell times "),
            make_push_button('Export',self.export_kinetics, "export plot contents"),
            make_push_button('Help',self.show_dwell_help, None)])
        #tab layout
        dwell_times_tab_layout = QHBoxLayout()
        dwell_times_tab_layout.addWidget(self.dwell_kinetics_canvas)
        dwell_times_tab_layout.addWidget(dwell_controls)

        #other
        other_graph_layout = QHBoxLayout()
        #other_graph_layout.addWidget(box2)

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

    def export_kinetics(self):
        fpath=self.parent.experiment.analysis_path / 'Dwell time analysis'
        wysiwyg_export(self.fig_kinetics, filepath=fpath,  filename="kinetics_export", filetype="csv")

    def perform_dwell_times_sequence(self):
        self.fig_kinetics.clear()
        axes = self.fig_kinetics.subplots(1, 2, sharex=True)

        selected_files = self.parent.experiment.selectedFiles

        self.dwell_kinetics_canvas.draw_idle()
        #here it would be good to pass button choices: variable, method
        if selected_files:
            selected_files.serial.determine_dwells_from_classification(variable='FRET', selected=True, inactivate_start_and_end_states=True)
            selected_files.serial.analyze_dwells(method='histogram_fit', number_of_exponentials=[1, 2])
            selected_files.serial.plot_dwell_analysis(plot_range=(0, 2), axes=axes, log=False)

    def show_dwell_help(self):

        help_text = """
                <html>
                  <body style="font-family: sans-serif; font-size: 10pt;">
                
                    <h2>Dwell Times</h2>
                
                    <p>
                      Choose which variable and method should be used for the dwell time analysis.
                    </p>
                
                    <ul>
                      <li>Select the variable for dwell time extraction</li>
                      <li>Select the analysis method</li>
                      <li>'Export' saves as-seen data from the graph panels to .csv and .png
                    </ul>
                
                    <p>
                      For more help, see the
                      <a href="https://papylio.readthedocs.io/en/stable/user_guide/dwell_time_analysis/index.html">
                        dwell time analysis documentation
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


