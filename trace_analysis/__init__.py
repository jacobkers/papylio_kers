import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))

from trace_analysis.traceAnalysisCode import Experiment, File, Molecule
from trace_analysis.image_adapt.movie import Movie
from trace_analysis.image_adapt.sifx_file import SifxFile