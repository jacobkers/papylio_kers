Getting started
===============
A regular analysis consists of mapping the emission channels, localizing molecules and
extracting the intensity traces.

First import the library and create an Experiment object.

.. code-block:: python

    import trace_analysis as ta
    exp = ta.Experiment()
    exp.print_files()

Calling "Experiment" without any arguments will open a dialog box where you can select the main folder containing
the single-molecule movies.

Creating the experiment object will either load an exiting configuration (.config) file that is present in the selected path,
or it will create an default config file which can then be adapted.

Emission channel mapping
------------------------
As a first step perform the mapping. To that end find the index of the mapping file (in the code below it is set to 0),
run the mapping and show the mapping result.

.. code-block:: python

    mapping_file_index = 0
    mapping_file = exp.files[mapping_file_index]

    mapping_file.perform_mapping()

    figure, axis = mapping_file.show_image()
    mapping_file.mapping.show(axis=axis, show_source=True)

Set the appropriate setting under the mapping heading in the config file.
Especially the "minimum_intensity_difference" will need to be adjusted separately for the different emission channels.
An initial guess for the value can be determined from the difference between background and molecule intensity in a plotted image.
Make sure that most of the spots in the images are detected.
In addition, it is important that the detected spots and the matched points are homogeneously spread over the field of view.

The mapping is saved in a .mapping file and is automatically loaded when importing an experiment.
Note that if multiple mapping files are present, the first mapping file in exp.files is used.
So it is good practice to remove the .mapping files that should not be used.

Background subtraction
----------------------


Molecule localization
---------------------
Once the emission channel mapping is obtained, the molecule coordinates can be determined.
As a first step determine the files you want to analyze.

.. code-block:: python
    # All files
    file_selection = exp.files

    # Selection of files (files 2 to 10)
    file_selection = exp.files[2:10]

To accurately localize molecules in images, likely the standard settings need to be adjusted.
Therefore, find the optimal settings by trial and error using one or several files.
Settings for this step can be found under the heading "find_coordinates" in the configuration file.

.. code-block:: python

    test_file = file_selection[0]
    test_file.find_coordinates()
    test_file.show_coordinates_in_image()

After adjusting the settings, the coordinates for the selection of files can be extracted.

.. code-block:: python

    file_selection.find_coordinates()

The found coordinates are stored in a .nc file with the same name as the movie file.
Note that each time "find_coordinates" is run the .nc file is overwritten.

Trace extraction
----------------
Now the intensity traces can be extracted at the found positions.

.. code-block:: python

    file_selection.extract_traces()

The intensity and corresponding FRET traces are added to the existing .nc file.

Trace visualization
-------------------
Traces can be visualized for a specific file using the "show_trace" method. This will open a window showing the traces.
The y-limits of the plots can be adjusted using the ylims keyword argument. In addition the colors of the plots can be changed.

.. code-block:: python

    test_file = file_selection[0]
    test_file.show_traces(plot_variables=['intensity', 'FRET'],
                         ylims=[(0, 35000), (0, 1)],
                         colours=[('green', 'red'), ('blue')])

You can go backward and forward through the traces by clicking the left and right arrows.
Clicking 's' will save the current plot in the "Trace_plots" directory in the main experiment folder.



