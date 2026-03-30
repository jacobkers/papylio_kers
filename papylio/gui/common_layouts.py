#this widget contains layouts that are used in multiple tabs.
#They are stored here to avoid circular imports -jk

import sys
from PySide2.QtWidgets import  QApplication, QSizePolicy
import matplotlib as mpl
import inspect


from matplotlib.backends.backend_qtagg import (
    FigureCanvas)

from PySide2.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,QToolButton, QTextBrowser,
    QLabel, QSizePolicy, QGridLayout, QDialog, QPushButton, QTextEdit, QDialog, QLineEdit
)


from PySide2.QtCore import Qt
import sys

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

class HelpDialog(QDialog):
    def __init__(self, parent=None, help_text=""):
        super().__init__(parent)
        self.setWindowTitle("Help")
        self.resize(900, 800)

        layout = QVBoxLayout(self)

        self.text = QTextBrowser()
        self.text.setOpenExternalLinks(True)
        self.text.setHtml(help_text)

        layout.addWidget(self.text)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)

        layout.addWidget(self.text)
        layout.addWidget(close_button)

def make_push_button(text, method, tooltip):
    #build a standardized push button
    btn = QPushButton(text)
    btn.setFixedHeight(40)  # small height
    btn.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    btn.setStyleSheet("QPushButton { padding: 2px 6px; font-size: 16px; }")
    btn.clicked.connect(method)
    if tooltip is not None:
        btn.setToolTip(tooltip)
    return btn

def build_control_layouts(button_list):
#controls are a  limited set (up to 5) of main action buttons
    controls_layout = QHBoxLayout()
    for b in button_list:
        controls_layout.addWidget(b, alignment=Qt.AlignBottom)
    # Optional: add stretch at the end to push buttons to the left
    controls_layout.addStretch()
    controls = QWidget()
    controls.setLayout(controls_layout)
    return controls


# -------------------------------------------------------------------------
# Register methods dynamically and create their forms
# (copied from classification widget with suffix _gen
# -------------------------------------------------------------------------
def register_method_gen(self, name, func):
    """Register a classification method, introspect arguments, and build a form."""

    self.methods[name] = func

    # --- build the form for the function ---
    form_widget = QWidget()
    form = QFormLayout(form_widget)
    inputs = {}

    sig = inspect.signature(func)
    for param_name, param in sig.parameters.items():
        if param_name in ['image','traces', 'classification', 'selection']:
            continue

        default = param.default if param.default is not inspect.Parameter.empty else None
        annotation = param.annotation

        # Pick appropriate input type
        if annotation == int or isinstance(default, int):
            widget = QSpinBox()
            widget.setRange(-1_000_000, 1_000_000)
            if default is not None:
                widget.setValue(default)
        elif annotation == float or isinstance(default, float):
            widget = QDoubleSpinBox()
            widget.setRange(-1e9, 1e9)
            widget.setDecimals(6)
            if default is not None:
                widget.setValue(default)
        else:
            widget = QLineEdit()
            if default not in (None, inspect.Parameter.empty):
                widget.setText(str(default))

        form.addRow(f"{param_name}:", widget)
        inputs[param_name] = widget

    self.method_forms[name] = (form_widget, inputs)
    self.method_selector.addItem(name)

    # First registered method becomes default
    if self.method_selector.count() == 1:
        self._update_method_panel_gen(name)


def _update_method_panel_gen(self, name):
    # Clear the old form
    for i in reversed(range(self.stack_layout.count())):
        widget = self.stack_layout.itemAt(i).widget()
        if widget:
            widget.setParent(None)

    # Add new form
    if name in self.method_forms:
        form_widget, _ = self.method_forms[name]
        self.stack_layout.addWidget(form_widget)


#below a series of functions to link GUI elements to Papylio methods in generic fashion
#Note to self: not sure if we are going to use this
def get_parameters(func):
    # for auto_generating GUI panels:
    #build a list of parameters for a given function.
    #NOTE: for nested settings such as in a dict, this will not go deeper
    sig = inspect.signature(func)
    params = {}

    for name, p in sig.parameters.items():
        if p.default is not inspect._empty:
            params[name] = p.default
        else:
            params[name] = None  # or some placeholder

    return params


def create_fields(func, layout):
    # for auto_generating GUI panels:
    #expand a of parameter fields
    # using the default values of this functions parameters
    param_defaults = get_parameters(func)
    widgets = {}

    for name, default in param_defaults.items():
        field = QLineEdit()
        if default is not None:
            field.setText(str(default))

        layout.addWidget(field)
        widgets[name] = field

    return widgets

def convert_value(text, default):
    # for auto_generating GUI panels:
    # restore type of the default value
    if default is None:
        return text  # fallback
    try:
        return type(default)(text)
    except Exception:
        return text

def call_with_widgets(func, widgets, defaults):
    # for auto_generating GUI panels:
    # activate the function
    kwargs = {}

    for name, widget in widgets.items():
        text = widget.text()
        default = defaults[name]
        kwargs[name] = convert_value(text, default)

    return func(**kwargs)

def bind_function_to_gui(func, layout, button):
    # for auto_generating GUI panels:
    # link a specified button to calling a function with
    # auto-generated parameters fields
    defaults = get_parameters(func)
    widgets = create_fields(func, layout)

    def on_click():
        result = call_with_widgets(func, widgets, defaults)
        print(result)

    button.clicked.connect(on_click)

    return widgets