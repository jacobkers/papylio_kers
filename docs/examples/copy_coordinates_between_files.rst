Copy coordinates between files
==============================
Copying coordinates between files can for example be useful if you obtain two movies with separate illuminations. You may want to find the coordinates using one illumination and extract the intensity traces from the movie with the other illumination.

You can do that easily do this by setting the coordinates of one file to the coordinates of the other file.

.. code-block:: python

    import trace_analysis as ta

    exp = ta.Experiment(r'path_to_data')

    exp.files[1].find_coordinates()
    exp.files[0].coordinates = exp.files[1].coordinates