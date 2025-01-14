from PySide2.QtWidgets import QApplication

from papilio.gui.main import MainWindow
import sys

from multiprocessing import Process, freeze_support
freeze_support()

app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec_()
