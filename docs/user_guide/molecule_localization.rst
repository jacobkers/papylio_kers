Molecule localization
=====================

Once the emission channel mapping is obtained, the molecule coordinates can be determined.
As a first step determine the files you want to analyze.

.. code-block:: python

    # All files
    file_selection = exp.files

    # Selection of files (files 2 to 10)
    file_selection = exp.files[2:10]

To accurately localize molecules in images, it is likely that the standard settings need to be adjusted.
Therefore, find the optimal settings by trial and error using one or several files.
Settings for this step can be found under the heading ``find_coordinates`` in the config file.

.. code-block:: python

    test_file = file_selection[0]
    test_file.find_coordinates()
    test_file.show_coordinates_in_image()

After adjusting the settings, the coordinates for the selection of files can be extracted.

.. code-block:: python

    file_selection.find_coordinates()

The found coordinates are stored in an .nc file with the same name as the movie file.
Note that each time "find_coordinates" is run the .nc file is overwritten.
