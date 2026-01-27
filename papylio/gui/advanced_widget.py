from PySide2.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QToolButton,
    QLabel, QSizePolicy, QGridLayout
)
from PySide2.QtCore import Qt
import sys


class Expander(QWidget):
    def __init__(self, title, parent=None):
        super().__init__(parent)

        self.toggle_button = QToolButton(text=title, checkable=True, checked=False)
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(Qt.RightArrow)
        self.toggle_button.setStyleSheet("QToolButton { border: none; }")

        self.content = QWidget()
        self.content.setVisible(False)

        layout = QGridLayout(self)
        layout.setSpacing(3)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setContentsMargins(16, 4, 0, 4)
        layout.addWidget(self.toggle_button)
        layout.addWidget(self.content)

        self.toggle_button.toggled.connect(self.on_toggled)

    def on_toggled(self, checked):
        self.content.setVisible(checked)
        self.toggle_button.setArrowType(
            Qt.DownArrow if checked else Qt.RightArrow
        )

    def setContentLayout(self, content_layout):
        self.content.setLayout(content_layout)

