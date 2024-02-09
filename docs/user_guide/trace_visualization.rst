Trace visualization
===================

Traces can be visualized for a specific file using the :py:meth:`.show_traces` method. This will open a window showing the traces.
The y-limits of the plots can be adjusted using the ``ylims`` keyword argument. In addition the colors of the plots can be changed.

.. code-block:: python

    test_file = file_selection[0]
    test_file.show_traces(plot_variables=['intensity', 'FRET'],
                         ylims=[(0, 35000), (0, 1)],
                         colours=[('green', 'red'), ('blue')])

You can go backward and forward through the traces by pressing the left and right arrows, respectively.
Pressing ``s`` will save the current plot in the ``Trace_plots`` directory in the main experiment folder.



