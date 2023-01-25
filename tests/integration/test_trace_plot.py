import pytest
import numpy as np
import xarray as xr
import sys

@pytest.fixture
def dataset(shared_datadir):
    return xr.open_dataset(shared_datadir / 'BN_TIRF_output_test_file' / 'TIRF 561 0001.nc' )

def test_trace_plot(dataset):
    from trace_analysis.trace_plot import TracePlotWindow

    from PySide2.QtWidgets import QApplication

    app = QApplication(sys.argv)
    frame = TracePlotWindow(dataset)
        #, "Sample editor", plot_variables=['intensity', 'FRET'],  # 'classification'],
        #          ylims=[(0, 1000), (0, 1), (-1,2)], colours=[('g', 'r'), ('b'), ('k')])

    app.exec_()

