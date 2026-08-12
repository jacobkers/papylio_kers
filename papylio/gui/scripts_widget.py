#import sys
#import json
from PySide2.QtWidgets import QHBoxLayout,  \
    QPushButton, QTabWidget, QComboBox, QFormLayout, QWidget, QLabel, QVBoxLayout, QSpinBox
from PySide2.QtCore import Qt

#from papylio import File
from papylio.gui.common_layouts import (HelpDialog, Group_Box, get_button_value,
                                        build_control_layouts, make_push_button)

#import numpy as np

from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

from PySide2.QtWidgets import QWidget, QVBoxLayout, QPlainTextEdit, QPushButton

from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager

import sys
import io

class ScriptBox(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.context = context   # objects available to script

        console = make_console({
            "project": self.experiment
        })

        layout = QVBoxLayout(self)
        layout.addWidget(console)



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



    def show_script_help(self):

        help_text = """
                <html>
                  <body style="font-family: sans-serif; font-size: 10pt;">
                
                    <h2>Custom Scripts</h2>
                
                    <p>
                      Paste and execute Papylio scripts
                    </p>
                
                    <ul>
                      <li>aaa</li>
                      <li>bbbb</li>
                      <li>cccc
                    </ul>
                
                    <p>
                      For more help, see the
                      <a href="https://papylio.readthedocs.io/en/stable/SPARXS/1_single_molecule_data_analysis.html">
                        Data Analysis with Papylio
                      </a>.
                    </p>
                
                  </body>
                </html>
                """
        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()

def make_console(context):

    kernel_manager = QtInProcessKernelManager()
    kernel_manager.start_kernel()

    kernel = kernel_manager.kernel
    kernel.shell.push(context)

    kernel_client = kernel_manager.client()
    kernel_client.start_channels()

    console = RichJupyterWidget()
    console.kernel_manager = kernel_manager
    console.kernel_client = kernel_client

    return console

