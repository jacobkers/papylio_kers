import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))

from trace_analysis.experiment import Experiment
from trace_analysis.file import File
from trace_analysis.molecule import Molecule
from trace_analysis.image_adapt.movie import Movie
from trace_analysis.image_adapt.sifx_file import SifxFile
from trace_analysis.analysis.interactiveAnalysis import InteractivePlot


# import trace_analysis.sequencing.sequencing
#
# Experiment = type(Experiment.__name__, (Experiment,) + tuple(Experiment.plugins), {})


