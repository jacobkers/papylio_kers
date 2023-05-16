Merge file datasets
===================
To merge datasets from multiple files into a single nc file, you can use the following function:

.. code-block:: python

    exp.files[0:10].merge_datasets(filepath_out='path/file.nc', with_selected_only=False)

This merges the datasets of the first 10 files and saves them in the filepath_out.