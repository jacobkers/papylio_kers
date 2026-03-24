
from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QTreeView, QApplication, QMainWindow, \
    QPushButton, QComboBox, QLineEdit, QLabel
from PySide2.QtCore import Qt

import matplotlib.pyplot as plt

from matplotlib.figure import Figure

from papylio import File
from papylio.gui.common_layouts import (Expander, ImageCanvas, HelpDialog,
                                        build_control_layouts,make_push_button, bind_function_to_gui)

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)


class MappingWidget(QWidget):
    def __init__(self, parent=None):
        super(MappingWidget, self).__init__(parent)
        self.parent = parent

        self.parent.model.itemChanged.connect(self.onItemChange)

        # imagery
        self.fig1 = Figure(figsize=(4, 3))
        self.map_overlay_image_canvas = FigureCanvas(self.fig1)
        map_overlay_image_layout = QVBoxLayout()
        map_overlay_image_layout.addWidget(self.map_overlay_image_canvas)

        self.map_overlay_image = QWidget()
        self.map_overlay_image.setLayout(map_overlay_image_layout)
        self.map_overlay_image.setToolTip("overlay image shows result of mapping")

        # imagery
        self.map_image_canvas = ImageCanvas(self, width=4, height=3, dpi=100)
        self.map_image_canvas.setToolTip("shows raw image (1st of multiple selection)")

        map_image_toolbar = NavigationToolbar(self.map_image_canvas, self)
        map_image_layout = QVBoxLayout()
        #map_image_layout.addWidget(map_image_toolbar)
        #map_image_layout.addWidget(self.map_image_canvas)

        # Create a placeholder widget to hold our toolbar and canvas.
        self.map_image = QWidget()
        self.map_image.setLayout(map_image_layout)


        #map settings
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
        #donor mapping:
        button_map_donor_label = QLabel("donor_pks")
        button_map_donor_method_combobox = QComboBox()
        button_map_donor_options = ['local-maximum-auto', 'other']
        button_map_donor_method_combobox.addItems(button_map_donor_options)
        button_map_donor_fract_diff = QLineEdit()
        button_map_donor_fract_diff.setPlaceholderText("fract._diff")
        button_map_donor_ns_min = QLineEdit()
        button_map_donor_ns_min.setPlaceholderText("nbh_min")
        button_map_donor_ns_max = QLineEdit()
        button_map_donor_ns_max.setPlaceholderText("nbh_max")
        # acceptor mapping:
        button_map_acceptor_label = QLabel("acceptor_pks")
        button_map_acceptor_method_combobox = QComboBox()
        button_map_acceptor_options = ['local-maximum-auto', 'other']
        button_map_acceptor_method_combobox.addItems(button_map_acceptor_options)
        button_map_acceptor_fract_diff = QLineEdit()
        button_map_acceptor_fract_diff.setPlaceholderText("fraction_diff")
        button_map_acceptor_ns_min = QLineEdit()
        button_map_acceptor_ns_min.setPlaceholderText("filt_nbh_min")
        button_map_acceptor_ns_max = QLineEdit()
        button_map_acceptor_ns_max.setPlaceholderText("filt_nbh_max")

        #basic mapping grid layout:
        map_basics_layout = QGridLayout()
        map_basics_layout.setAlignment(Qt.AlignLeft)
        map_basics_layout.addWidget(self.button_map_label, 0, 0)
        map_basics_layout.addWidget(self.button_map_method_combobox, 0, 1)
        map_basics_layout.addWidget(self.button_map_dist_treshold, 0, 2)

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
        mapping_tab_layout.addWidget(self.map_image)
        mapping_tab_layout.addWidget(self.map_overlay_image)
        mapping_tab_layout.addWidget(self.map_controls)


        self.setLayout(mapping_tab_layout)

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

    def perform_mapping(self, t):
        print(t)
        self.fig1.clear()
        ax1 = self.fig1.add_subplot(111)

        selected_files = self.parent.experiment.selectedFiles

        #jk-read in a default config and change one value:
        panel_config= self.parent.experiment.configuration['mapping']

        panel_config['method']= self.button_map_method_combobox.currentText()
        panel_config['distance_threshold']= self.button_map_dist_treshold.text
        panel_config['coordinates_within_margin'] = self.button_map_margin.text

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

