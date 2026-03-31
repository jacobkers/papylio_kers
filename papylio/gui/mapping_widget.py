
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QComboBox, QLineEdit, QLabel, QFormLayout, QDoubleSpinBox,QSpinBox
from PySide2.QtCore import Qt

import matplotlib.pyplot as plt

from matplotlib.figure import Figure

from papylio import File
from papylio.gui.common_layouts import (Expander, ImageCanvas, HelpDialog,
                                        build_control_layouts,make_push_button)

from papylio.peak_finding import (find_peaks_absolute_threshold,
                                  find_peaks_adaptive_threshold,
                                  find_peaks_local_maximum,
                                  find_peaks_local_maximum_auto,
                                  find_peaks_relative_local_maximum)

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

import inspect

class MappingWidget(QWidget):
    def __init__(self, parent=None):
        super(MappingWidget, self).__init__(parent)
        self.parent = parent
        self.methods = {}
        self.method_forms = {}  # method_name -> (widget, inputs)

        self.parent.model.itemChanged.connect(self.onItemChange)

        # imagery
        # #--------------------------------------------
        self.fig1 = Figure(figsize=(4, 3))
        self.map_overlay_image_canvas = FigureCanvas(self.fig1)
        map_overlay_image_layout = QVBoxLayout()
        map_overlay_image_layout.addWidget(self.map_overlay_image_canvas)
        self.map_overlay_image = QWidget()
        self.map_overlay_image.setLayout(map_overlay_image_layout)
        self.map_overlay_image.setToolTip("overlay image shows result of mapping")
        self.map_image_canvas = ImageCanvas(self, width=4, height=3, dpi=100)
        self.map_image_canvas.setToolTip("shows raw image (1st of multiple selection)")
        map_image_toolbar = NavigationToolbar(self.map_image_canvas, self)
        map_image_layout = QVBoxLayout()

        # Create a placeholder widget to hold our toolbar and canvas.
        self.map_image = QWidget()
        self.map_image.setLayout(map_image_layout)


        #build a flexible form for the donor channel:-------------------------
        form_donor = QFormLayout()
        # --- Classification name input ---
        donor_label = QLabel("Spot detection donor")
        donor_label.setToolTip("Tune the peak detection for the donor channel")
        form_donor.addRow(donor_label)

        # --- Method selector donor---
        self.method_selector_donor = QComboBox()
        self.method_selector_donor.setToolTip("Choose peak_find_method")
        self.method_selector_donor.currentTextChanged.connect(self._update_method_panel_donor)
        form_donor.addRow("Method:", self.method_selector_donor)

        # --- Dynamic options container -donor ---
        self.stack_donor = QWidget()
        self.stack_donor_layout = QVBoxLayout(self.stack_donor)
        self.stack_donor_layout.setContentsMargins(0, 0, 0, 0)
        form_donor.addRow("Options:", self.stack_donor)

        # build a flexible form for the acceptor channel:-------------------------
        form_acceptor = QFormLayout()
        # --- Classification name input ---
        acceptor_label = QLabel("Spot detection acceptor")
        acceptor_label.setToolTip("Tune the peak detection for the donor channel")
        form_acceptor.addRow(acceptor_label)

        # --- Method selector acceptor---
        self.method_selector_acceptor = QComboBox()
        self.method_selector_acceptor.setToolTip("Choose peak_find_method")
        self.method_selector_acceptor.currentTextChanged.connect(self._update_method_panel_acceptor)
        form_acceptor.addRow("Method:", self.method_selector_acceptor)

        # --- Dynamic options container -acceptor ---
        self.stack_acceptor = QWidget()
        self.stack_acceptor_layout = QVBoxLayout(self.stack_acceptor)
        self.stack_acceptor_layout.setContentsMargins(0, 0, 0, 0)
        form_acceptor.addRow("Options:", self.stack_acceptor)


        #map settings
        #-------------------------------------------------------
        #main mapping buttons definition:
        self.button_map_label = QLabel("method:")
        self.button_map_method_combobox = QComboBox()
        self.button_map_options = ['icp', 'nn']
        self.button_map_method_combobox.addItems(self.button_map_options)
        self.button_map_dist_treshold = QLineEdit()
        self.button_map_dist_treshold.setPlaceholderText("dist_treshold")
        self.button_map_margin = QLineEdit()
        self.button_map_margin.setPlaceholderText("edge_margin")

        # add below buttons under 'advanced':---------------------------------------------------

        #basic mapping grid layout:
        #donor_acceptor block
        donor_acceptor_layout = QHBoxLayout()
        donor_acceptor_layout.addLayout(form_donor)
        donor_acceptor_layout.addLayout(form_acceptor)

        map_basics_layout = QVBoxLayout()
        map_basics_layout.addLayout(donor_acceptor_layout)

        map_basics_layout.setAlignment(Qt.AlignLeft)

        map_basics_layout.addWidget(self.button_map_label )
        map_basics_layout.addWidget(self.button_map_method_combobox)

        # advanced mapping grid layout:
        map_advanced_layout = QGridLayout()
        map_advanced_layout.addWidget(self.button_map_margin, 1, 1)
        map_advanced_layout.addWidget(button_map_donor_label, 2, 0)
        map_advanced_layout.addWidget(button_map_donor_method_combobox, 2, 1)
        map_advanced_layout.addWidget(button_map_donor_fract_diff, 2, 2)
        map_advanced_layout.addWidget(button_map_donor_ns_min, 2, 3)
        map_advanced_layout.addWidget(button_map_donor_ns_max, 2, 4)
        map_advanced_layout.addWidget(button_map_acceptor_label, 3, 0)
        map_advanced_layout.addWidget(button_map_acceptor_method_combobox, 3, 1)
        map_advanced_layout.addWidget(button_map_acceptor_fract_diff, 3, 2)
        map_advanced_layout.addWidget(button_map_acceptor_ns_min, 3, 3)
        map_advanced_layout.addWidget(button_map_acceptor_ns_max, 3, 4)

        #build panel layout:
        map_basics = QWidget()
        map_basics.setLayout(map_basics_layout)
        map_advanced = Expander("Advanced")
        map_advanced.setContentLayout(map_advanced_layout)


        #main action:
        map_controls = build_control_layouts(
            [make_push_button('Map', self.perform_mapping,"Map selected file(s)"),
             make_push_button('Help', self.show_help, None)])


        #collect:
        map_controls_layout = QVBoxLayout()
        map_controls_layout.setAlignment(Qt.AlignTop)
        map_controls_layout.addWidget(map_basics)
        map_controls_layout.addWidget(map_advanced)
        map_controls_layout.addWidget(map_controls)

        #pack in widget:
        self.map_controls = QWidget()
        self.map_controls.setLayout(map_controls_layout)
        self.map_controls.setMinimumWidth(150)

        #add all to tab:
        mapping_tab_layout = QHBoxLayout()
        mapping_tab_layout.addWidget(self.map_controls)
        mapping_tab_layout.addWidget(self.map_image)
        mapping_tab_layout.addWidget(self.map_overlay_image)

        self.setLayout(mapping_tab_layout)
        #collect peak finding methods for building flexible GUI forms
        self.register_method('absolute_treshold', find_peaks_absolute_threshold)
        self.register_method('adaptive_treshold', find_peaks_adaptive_threshold)
        self.register_method('local_maximum', find_peaks_local_maximum)
        self.register_method('auto-local_maximum', find_peaks_local_maximum_auto)
        self.register_method('relative local_maximum', find_peaks_relative_local_maximum)


    def update_plots(self):
        selected_files = self.parent.experiment.selectedFiles + [None]
        self.map_image_canvas.file = selected_files[0]
         #if selected_files[0] is not None:
         #     self.traces.dataset = selected_files[0].dataset
         #      self.selection.file = selected_files[0]
         #else:
         #     self.traces.dataset = None
         #      self.selection.file = None

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

    # -------------------------------------------------------------------------
    # Register methods dynamically and create their forms
    # -------------------------------------------------------------------------
    def register_method(self, name, func):
        """Register a classification method, introspect arguments, and build a form."""

        self.methods[name] = func

        # --- build the form for the function ---
        form_widget = QWidget()
        form = QFormLayout(form_widget)
        inputs = {}

        sig = inspect.signature(func)
        for param_name, param in sig.parameters.items():
            if param_name in ['image']:
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

        self.method_selector_donor.addItem(name)
        self.method_selector_acceptor.addItem(name)

        # First registered method becomes default
        if self.method_selector_donor.count() == 1:
            self._update_method_panel_donor(name)
        if self.method_selector_acceptor.count() == 1:
            self._update_method_panel_acceptor(name)

    def _update_method_panel_donor(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_donor_layout.count())):
            widget = self.stack_donor_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms:
            form_widget, _ = self.method_forms[name]
            self.stack_donor_layout.addWidget(form_widget)

    def _update_method_panel_acceptor(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_acceptor_layout.count())):
            widget = self.stack_acceptor_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms:
            form_widget, _ = self.method_forms[name]
            self.stack_acceptor_layout.addWidget(form_widget)

    def perform_mapping(self, t):
        print(t)
        self.fig1.clear()
        ax1 = self.fig1.add_subplot(111)

        selected_files = self.parent.experiment.selectedFiles

        #jk need to do a mapping here
        donor_kwargs = 1
        acceptor_kwargs =1

        #jk-read in a default configuration and change one value:
        panel_config= self.parent.experiment.configuration['mapping']
        panel_config['method']= self.button_map_method_combobox.currentText()
        panel_config['distance_threshold']= self.button_map_dist_treshold.text
        panel_config['coordinates_within_margin'] = self.button_map_margin.text
        panel_config['peak_finding']['donor'] = donor_kwargs
        panel_config['peak_finding']['acceptor'] = acceptor_kwargs

        plot_file = selected_files[0]
        plot_file.mapping.show_mapping_transformation(axis=ax1)
        self.map_overlay_image_canvas.draw_idle()
        if selected_files:
            selected_files.serial.perform_mapping(**panel_config)
            self.update_plots()
            self.map_image_canvas.refresh()

    def show_help(self):
        help_text = help_text = """\
               <html>
                  <body style="font-family: sans-serif; font-size: 10pt;">
                
                    <h2>Mapping</h2>
                
                    <p>
                      Mapping links coordinates in the used color channels
                    </p>
                
                    <ul>
                      <li>adjust main settings</li>
                      <li>adjust 'advanced' settings</li>
                      <li>description: [options], first value is default
                    </ul>
                
                    <p>
                    
                    <h3>Basic settings</h3>
                    
                    <h4>mapping</h4>
                    mapping: choose how pairs are found between source and target
                    <ul>
                      <li><code>method</code> [<code>icp</code> | <code>nn</code>]
                        <ul>
                          <li><code>icp</code> = <code>iterative_closest_point</code>>
                          <li><code>nn</code> = <code>nearest_neighbor</code>: two-way nearest neighbor within a distance threshold</li>
                          <li><code>direct</code> = <code>direct_match</code>: map points 1:1 
                          <li><code>distance_threshold[3]</code>: beyond this, no 'nn' pairs are considered
                          <li><code>transformation_type['polynomial',...]</code>: method to perform link [EDIT]
                          <ul>
                             <li><code>'linear' | 'affine'</code> : skimage.transform.AffineTransform
                             <li><code>'similarity'</code> :  skimage.transform.SimilarityTransform
                             <li><code> 'nonlinear'</code> :  IDL polywarp transform
                             <li><code> 'polynomial'</code> : skimage.transform.PolynomialTransform
                          </ul>
                          <li><code>initial_translation ['width/2'| n_pixels]/code>: initial shift for chosen mapping method
                       </ul>
                   </ul>  
                     
                <h3>Advanced settings</h3> 
                peak_finding: choose how spots are detected per donor or acceptor channel           
                <ul>                    
                    <li><code>donor</code>:
                        <ul>
                            <li><code> method[local-maximum-auto] </code>
                            <li><code> filter_neighbourhood_size_min[10]</code> [tbd]
                            <li><code> filter_neighbourhood_size_max[5]</code> [tbd]
                        </ul>
                    <li><code>acceptor</code>:
                        <ul>
                        <li><code> method: local-maximum-auto</code>
                        <li><code> filter_neighbourhood_size_min[10]</code>
                        <li><code> filter_neighbourhood_size_max[5]</code>
                        </ul>
                    <li><code>coordinate_optimization/code>
                        <ul>
                            <li> coordinates_after_gaussian_fit: gaussian_width: 5
                            <li> coordinates_within_margin: margin: 10
                        </ul>
    
                <h3>More help</h3>
                <p>
                See the 
                    <a href="https://papylio.readthedocs.io/en/stable/user_guide/channel_mapping.html">
                mapping documentation
                    </a>.
                </p>
                   

                
                </body>
                </html>
                """

        self.help_dialog = HelpDialog(self, help_text)
        # dialog.exec_()  # modal
        self.help_dialog.show()

