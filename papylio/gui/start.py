from PySide2.QtWidgets import QApplication

from papylio.gui.main import MainWindow
import sys

from multiprocessing import Process, freeze_support


def start_gui():
    freeze_support()

    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec_()


if __name__ == '__main__':
    start_gui()