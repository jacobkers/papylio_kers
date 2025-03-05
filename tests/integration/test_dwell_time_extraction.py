import numpy as np
import pytest

from papylio.analysis.dwell_time_extraction import dwell_times_from_classification

@pytest.fixture
def dwells():
    classification = np.array([-1, -1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0]*10)
    classification = np.repeat(classification[None, :], 5, axis=0)
    traces = np.random.random(classification.shape)

    return dwell_times_from_classification(classification, traces=traces, cycle_time=0.1)

