import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1]))

from trace_analysis.experiment import Experiment
from trace_analysis.file import File
# from trace_analysis.molecule import Molecule
from trace_analysis.movie.movie import Movie
from trace_analysis.analysis.interactiveAnalysis import InteractivePlot

# import trace_analysis.sequencing.sequencing
#
# Experiment = type(Experiment.__name__, (Experiment,) + tuple(Experiment.plugins), {})

try:
    from ._version import version as __version__
except ImportError:
    try:
        import setuptools_scm
        __version__ = setuptools_scm.get_version(version_scheme="post-release")
    except ImportError:
        __version__ = '0.0.0'