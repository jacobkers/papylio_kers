import pytest
from trace_analysis.plugins.sequencing.sequencing_data import parse_sam, fastq_generator, add_sequence_data_to_dataset


@pytest.fixture
def read1_fastq_filepath(shared_datadir):
    return shared_datadir / 'sequencing' / 'read1.fastq'

@pytest.fixture
def index1_fastq_filepath(shared_datadir):
    return shared_datadir / 'sequencing' / 'index1.fastq'

@pytest.fixture
def sam_filepath(shared_datadir):
    return shared_datadir / 'sequencing' / 'aligned.sam'


def test_parse_sam(sam_filepath):
    parse_sam(sam_filepath, remove_duplicates=True, add_aligned_sequence=True, extract_sequence_subset=False,
              chunksize=10, write_csv=False, write_nc=True, write_filepath=None)

def test_fastq_generator(read1_fastq_filepath):
    fqg = fastq_generator(read1_fastq_filepath)
    next(fqg)

def test_add_sequence_data_to_dataset(sam_filepath, index1_fastq_filepath):
    test_parse_sam(sam_filepath)
    nc_filepath = sam_filepath.with_suffix('.nc')
    add_sequence_data_to_dataset(nc_filepath, index1_fastq_filepath, 'index1')

