
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QComboBox, QLineEdit, QLabel, QFormLayout, QDoubleSpinBox,QSpinBox,QFrame
from PySide2.QtCore import Qt

import matplotlib.pyplot as plt

from matplotlib.figure import Figure

from papylio import File
from papylio.gui.common_layouts import (Expander, ImageCanvas, HelpDialog,
                                        build_control_layouts,make_push_button,
                                        build_form,build_parameters_input)

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
        self.methods_donor = {}
        self.method_forms_donor = {}  # method_name -> (widget, inputs)
        self.methods_acceptor = {}
        self.method_forms_acceptor = {}  #

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
        frame_donor = QFrame()
        frame_donor.setFrameShape(QFrame.StyledPanel)
        frame_donor.setFrameShadow(QFrame.Plain)
        frame_donor.setLineWidth(2)

        form_donor = QFormLayout(frame_donor)
        donor_label = QLabel("DONOR spot detection")
        donor_label.setToolTip("Tune the peak detection for the donor channel")
        form_donor.addRow(donor_label)
        # --- Method selector donor---
        self.method_selector_donor = QComboBox()
        self.method_selector_donor.setToolTip("Choose peak_find_method")
        self.method_selector_donor.currentTextChanged.connect(self._update_method_panel_donor)
        form_donor.addRow("Method:", self.method_selector_donor)
        # --- Dynamic options container donor ---
        self.stack_donor = QWidget()
        self.stack_donor_layout = QVBoxLayout(self.stack_donor)
        self.stack_donor_layout.setContentsMargins(0, 0, 0, 0)
        form_donor.addRow("Options:", self.stack_donor)




        # build a flexible form for the acceptor channel:-------------------------
        frame_acceptor = QFrame()
        frame_acceptor.setFrameShape(QFrame.StyledPanel)
        frame_acceptor.setFrameShadow(QFrame.Plain)
        frame_acceptor.setLineWidth(2)

        form_acceptor = QFormLayout(frame_acceptor)
        acceptor_label = QLabel("ACCEPTOR spot detection ")
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

        #genral map settings
        #-------------------------------------------------------
        #main mapping buttons definition:


        #basic mapping grid layout:
        #donor_acceptor block
        donor_acceptor_layout = QHBoxLayout()
        donor_acceptor_layout.addWidget(frame_donor)
        donor_acceptor_layout.addWidget(frame_acceptor)






        self.button_map_label = QLabel("method:")
        self.button_map_method_combobox = QComboBox()
        self.button_map_options = ['icp', 'nn']
        self.button_map_method_combobox.addItems(self.button_map_options)
        self.button_map_dist_treshold = QLineEdit()
        self.button_map_dist_treshold.setPlaceholderText("distance_treshold")
        self.button_transformation = QComboBox()
        self.button_transformation_options = ['polynomial', 'linear', 'nonlinear']
        self.button_transformation.addItems(self.button_transformation_options)
        self.button_initial_translation = QComboBox()
        self.button_initial_translation_options = ['width/2' , [1024,0]]
        self.button_initial_translation.addItems(self.button_initial_translation_options)



        # advanced mapping grid layout:
        frame = QFrame()
        frame.setFrameShape(QFrame.StyledPanel)
        frame.setFrameShadow(QFrame.Plain)
        frame.setLineWidth(2)

        map_advanced_layout = QVBoxLayout()
        map_advanced_layout.setAlignment(Qt.AlignLeft)

        map_advanced2_layout = QFormLayout(frame)
        map_advanced2_layout.addRow("Method", self.button_map_method_combobox)
        map_advanced2_layout.addRow("Distance treshold:",self.button_map_dist_treshold)
        map_advanced2_layout.addRow("Transformation_type:", self.button_transformation)
        map_advanced2_layout.addRow("Initial_translation:", self.button_initial_translation)

        map_advanced_layout.addWidget(frame)
        #map_advanced_layout.addLayout(map_advanced2_layout)
        map_advanced_layout.addLayout(donor_acceptor_layout)


        #build panel layout:
        map_advanced = Expander("Advanced")
        map_advanced.setContentLayout(map_advanced_layout)

        #main action:
        map_controls = build_control_layouts(
            [make_push_button('Map', self.perform_mapping,"Map selected file(s)"),
             make_push_button('Help', self.show_help, None)])


        #collect:
        map_controls_layout = QVBoxLayout()
        map_controls_layout.setAlignment(Qt.AlignTop)
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
        self.register_method('local-maximum-auto', find_peaks_local_maximum_auto)  #default  in GUI
        self.register_method('absolute-threshold', find_peaks_absolute_threshold)
        self.register_method('adaptive-threshold', find_peaks_adaptive_threshold)
        self.register_method('local-maximum', find_peaks_local_maximum)
        self.register_method('relative-local-maximum', find_peaks_relative_local_maximum)




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
        """Register a peak finding method, introspect arguments,
        and build forms for donor and acceptor channels"""

        #donor-------------------------:
        form_widget_donor, inputs_donor = build_form(func)
        self.methods_donor[name] = func
        self.method_forms_donor[name] = (form_widget_donor, inputs_donor)
        self.method_selector_donor.addItem(name) #add options to appropriate selector box
        if self.method_selector_donor.count() == 1: # First registered method becomes default
            self._update_method_panel_donor(name)

        # acceptor-------------------------:
        form_widget_acceptor, inputs_acceptor = build_form(func)
        self.methods_acceptor[name] = func
        self.method_forms_acceptor[name] = (form_widget_acceptor, inputs_acceptor)
        self.method_selector_acceptor.addItem(name)  # add options to appropriate selector box
        if self.method_selector_acceptor.count() == 1:  # First registered method becomes default
            self._update_method_panel_acceptor(name)


    def _update_method_panel_donor(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_donor_layout.count())):
            widget = self.stack_donor_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms_donor:
            form_widget, _ = self.method_forms_donor[name]
            self.stack_donor_layout.addWidget(form_widget)

    def _update_method_panel_acceptor(self, name):
        # Clear the old form
        for i in reversed(range(self.stack_acceptor_layout.count())):
            widget = self.stack_acceptor_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        # Add new form
        if name in self.method_forms_acceptor:
            form_widget, _ = self.method_forms_acceptor[name]
            self.stack_acceptor_layout.addWidget(form_widget)

    def perform_mapping(self, t):
        print(t)
        self.fig1.clear()
        ax1 = self.fig1.add_subplot(111)

        selected_files = self.parent.experiment.selectedFiles

        #get methods and corresponding parameters
        #donor:
        method_name_donor = self.method_selector_donor.currentText()
        _, inputs_donor = self.method_forms_donor[method_name_donor]
        method_name_acceptor = self.method_selector_acceptor.currentText()
        _, inputs_acceptor = self.method_forms_acceptor[method_name_acceptor]

        # Collect args for peak finding
        donor_kwargs=build_parameters_input(method_name_donor, inputs_donor)
        acceptor_kwargs = build_parameters_input(method_name_acceptor, inputs_acceptor)

        #jk-read in a default configuration and allocate GUI values:
        panel_config= self.parent.experiment.configuration['mapping']
        panel_config['method']= self.button_map_method_combobox.currentText()
        panel_config['distance_threshold']= self.button_map_dist_treshold.text
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


