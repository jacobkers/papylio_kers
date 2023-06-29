Dwell time analysis
===================

Dwell times can be extracted from the overall classification using the ``determine_dwells_from_classification`` method:.

.. code-block::

    determine_dwells_from_classification(variable='FRET', selected=True)

The method will create an additional .nc file named ``<filename>_dwells.nc``.
From which the contents can be obtained using:

.. code-block::

    file.dwells

For each dwell the following information is provided:

* ``file``:  the file from which it originates
* ``molecule_in_file``: the molecule index in the original file
* ``state``: the state of the dwell
* ``frame_count``: the dwell duration in frames
* ``duration``: the dwell duration in time
* ``mean_<variable_name>``: the mean of the used variable

The dwells with positive states that have a negative state neighbor are inactivated (i.e. set to state ``-128``) by default.

The dwell times can be analyzed using

.. code-block::

    file.analyze_dwells(plot=True, state_names={0: 'Low FRET state', 1: 'High FRET state'}, logy=False):

The ``analyze_dwells`` method will construct a dwell time histogram for each of the given states and will fit
a single decaying exponential function

.. math::

   count = Ae^{-kt}


The method will create a new file named ``<filename>_dwell_analysis.nc`` containing the fit values, errors and covariance.
In addition, if ``plot`` is set to `True`, the plot of the histogram and fit will be saved.
The ``logy`` variable can be set to ``True`` for a logarithmic y-axis.