Trace selection
===============

To select or filter molecules, selection variables can be set.
A selection variable is a boolean array specifying whether a molecule is included or not.

For example you can perform a selection based on the 5-frame rolling average of the total intensity,
this can eliminate molecules with multiple fluorescent labels.
To this end first calculate the 5-frame average of the rolling intensity using xarray.
Then determine the selection by finding the traces that are lower than a threshold for all frames.
Add this selection to the file dataset using the `set_variable` attribute of the `File` class.
Any name can be used that starts with `selection_`.

.. code-block::

    intensity_total_maximum_threshold = 20000
    intensity_total_rolling = file.intensity_total.rolling(frame=5, center=True).mean().dropna('frame')
    selection = (intensity_total_rolling < intensity_total_maximum_threshold).all('frame')
    file.set_variable(selection, name='selection_intensity_total')

Another selection criterion may be requiring that there is acceptor signal at the end of the trace.

.. code-block::

    intensity_red_threshold = 5000
    intensity_red_end = file.intensity.sel(channel=1, frame=slice(-5,None)).mean('frame')
    selection = intensity_red_threshold < intensity_red_end
    file.set_variable(selection, name='selection_intensity_red_end')

The selections in a file can be obtained using the `file.selections` property.

Now an overall selection can be made using

.. code-block::

    file.apply_selections()

This will combine all selections variables, where a molecule will be selected if all of the selections are True.
If you want to set the selection based on one or multiple selection variables, the names of the selection variables can
be passed to `apply_selections`:

.. code-block::

    file.apply_selections(['selection_intensity_total', 'selection_intensity_red_end'])

This can be useful for inspecting other combinations of selection criteria, while keeping
the individual selection variables untouched.