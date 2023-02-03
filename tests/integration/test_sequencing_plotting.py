import numpy as np
import xarray as xr
from trace_analysis.plugins.sequencing.plotting import double_mutations, plot_double_mutations
import matplotlib.pyplot as plt

general_sequence = 'ACTGACTG'

# Generate test dataset
dm = double_mutations(general_sequence)

da = xr.DataArray(np.arange(len(dm)), dims=('sequence',), coords={'sequence': dm})
da.name = 'Test dataset'

plot_double_mutations(general_sequence, da)

plot_double_mutations(general_sequence, da, da_annotation=da)