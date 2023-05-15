Save dataset with selected molecules
====================================
If you have selected specific molecules, for example while going through the time traces,
then you can save a netCDF dataset with only selected molecules by running:

.. code-block:: python

    exp.files[i].save_dataset_selected()

The saved nc file will have the same name as the original file, but with '_selected' appended to it.