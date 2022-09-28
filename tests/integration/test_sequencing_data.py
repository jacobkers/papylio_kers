import pytest
from trace_analysis.plugins.sequencing.sequencing_data import parse_sam


@pytest.fixture
def sam_filepath(shared_datadir):
    return shared_datadir / 'sequencing' / 'aligned.sam'


def test_parse_sam(sam_filepath):
    parse_sam(sam_filepath, remove_duplicates=True, add_aligned_sequence=True, extract_sequence_subset=False,
              chunksize=10, write_csv=False, write_nc=True, write_filepath=None)
