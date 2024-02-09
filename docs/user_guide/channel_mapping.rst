Channel mapping
===============

As a first step perform the mapping. To that end, find the index of the mapping file (in the code below it is set to 0),
run the mapping and show the mapping result.

.. code-block:: python

    mapping_file_index = 0
    mapping_file = exp.files[mapping_file_index]

    mapping_file.perform_mapping()

    figure, axis = mapping_file.show_image()
    mapping_file.mapping.show(axis=axis, show_source=True)

Set the appropriate setting under the mapping heading in the config file.
Especially the ``minimum_intensity_difference`` will need to be adjusted separately for the different emission channels.
An initial guess for the value can be determined from the difference between background and molecule intensity in a plotted image.
Make sure that most of the spots in the images are detected.
In addition, it is important that the detected spots and the matched points are homogeneously spread over the field of view.

The mapping is saved in a .mapping file and is automatically loaded when importing an experiment.
Note that if multiple mapping files are present, the first mapping file in ``exp.files`` is used.
So it is good practice to remove the .mapping files that should not be used.