Trace classification
====================

For trace classification, classification variables can be set.
A classification variable is a boolean or integer array specifying the state of each molecule at each frame.

Examples of simple classifications are minimum and maximum threshold.
First determine the classification.
To make the classification less dependent on noise we can, for example, apply a rolling median of the classification.
Then add the classification to the file dataset using the ``set_variable`` attribute of the ``File`` class.
Any name can be used that starts with ``classification_``.

.. code-block::

    window_size = 11
    intensity_total_min = 3500

    classification = file.intensity_total > intensity_total_min
    classification = classification.astype(int).rolling(frame=window_size, center=True, min_periods=1).median().astype(bool)
    file.set_variable(classification, name='classification_donor_active', dims=('molecule','frame'))

    classification = file.intensity_total < intensity_total_max
    classification = classification.astype(int).rolling(frame=window_size, center=True, min_periods=1).median().astype(bool)
    file.set_variable(classification, name='classification_single_dye', dims=('molecule','frame'))

Now an overall classification can be made using the ``file.apply_classifications()`` method,
which combines multiple classifications.
This can be useful for inspecting other combinations of classifications, while keeping
the individual classification variables untouched.

The ``apply_classifications`` method uses keyword arguments where the keys are the classification variable names
and the values are the state index or indices to use for the specific classification.

.. code-block::

    file.apply_classifications(classification_donor_active=-1, classification_single_dye=-2)

The default state is always 0, which is overwritten if a state is assigned.
If boolean classifications are passed
then for a positive state index the ``True`` values will be used for state assignment,
and for a negative state index the ``False`` values will be used for state assignment.

In the example, the state is zero,
except if the ``classification_donor_active`` is ``False``, then the state -1 is assigned; or
if the ``classification_single_dye`` is ``False``, then the state -2 is assigned.

The order of giving the keyword arguments is important, as the classifications are applied in order given.
In the example, if both ``classification_donor_active`` and ``classification_single_dye`` are ``False``,
then the final state will be ``-2``.
In addition, both the ``False`` and ``True`` state
s can be assigned with specific states.
The example below would be equivalent to the example above.

.. code-block::

    file.apply_classifications(classification_donor_active=[-1,0], classification_single_dye=[-2,0])

In line with this, if an integer classification is given then the passed states
will be assigned to the states present in the classification in ascending order.
Passing ``None`` excludes the state from being used.

.. code-block::

    file.apply_classifications(classification=[None, 0, 1])

Here the lowest state in ``classification`` is not applied,
the second state is set to ``0`` and the third state is set to ``1``.