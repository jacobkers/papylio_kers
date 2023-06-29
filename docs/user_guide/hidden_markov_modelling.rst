Hidden Markov modelling
=======================

Hidden Markov modelling can be used to classify traces.

If desired, first apply general selections to exclude traces from being analyzed. (See also: Trace selection)
And additionally apply classifications to exclude negatively classified parts of the traces
from being used in the fitting process of the hidden markov model. (See also: :doc:`/user_guide/trace_classification`)

.. code-block::

    file.apply_selections(['selection_intensity_total'])
    file.apply_classifications(classification_donor_active=-1, classification_single_dye=-2)

Hidden Markov modelling can then be applied.
Currently only a single variable can be used and only a two-state model is supported.
Both a one-state and a two-state model are fit and based on the Bayesian information criterion the best model is selected.

.. code-block::

    file.classify_hmm(variable='FRET')

The ``classify_hmm`` method will add several parameters to the file dataset.

* ``number_of_states``: the number of states used in the model for each trace.
* ``state_mean``: the mean value, e.g. FRET, for the each state.
* ``state_standard_deviation``: the standard deviation for the each state.
* ``transition_probability``: the transition probability from each state to each state
* ``start_probability``: the probability to start in a states
* ``end_probability``: the probability to end in a state
* ``transition_rate``: the transition rate from each state to each state

It will also add several selections:

* ``selection_complex_rates`: if ``False``, then calculating the transition rate from the transition probability resulted in a complex value.
* ``selection_lower_rate_limit`: if ``False``, then the calculated transition rate is lower than can be measured in the time window.

Finally it will add the found classification as the ``classification_hmm`` variable in the file dataset.
Any region not used for fitting using input classifications will be set to the value ``-1``.
For the region that are fitted, the state will be ``0`` in case a one-state model is used and
``0`` or ``1`` in case a two-state model is used.

Selections and classifications have to be applied to the general classification:

.. code-block::

    file.apply_selections(['selection_intensity_total',
                           'selection_complex_rates', 'selection_lower_rate_limit'])

    file.apply_classifications(classification_donor_active=-1, classification_single_dye=-2,
                               classification_hmm=[None, 0, 1])
