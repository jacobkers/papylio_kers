Sequence alignment
==================

The output of the sequencing run is a fastq file with sequences and other additional information such as quality scores
and coordinates. All these sequences have to be aligned to the input sequence(s) because there might be some insertions,
deletions or mismatches. If multiple samples are present, this step is also to determine which sequence matches which sample.
The output of sequence alignment is a sam file, which is needed for the step where the single-molecule and sequencing data are linked.

We suggest using Bowtie 2 for the sequence alignment. All files needed for sequence alignment can be found in the sequence_analysis folder.
If this is the first time, run the Sequencing_setup file to create the right environment. More detailed instructions on
for example how to run Linux on Windows can be found in the Analysis.sh file in the sequence_analysis folder.

In short, these are the steps for sequence alignment:

#. Copy everything from the sequence_analysis folder to the folder with your sequencing data (not on a network location)
#. Adjust the sequences in the reference file
#. Open a Linux terminal
#. Navigate to the folder with your sequencing data
#. Run the lines below:

.. code-block::

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate sequence_analysis

    zcat *R1_001.fastq.gz > Read1.fastq
    # zcat *I1_001.fastq.gz > Index1.fastq # In case the index is sequenced.

    bowtie2-build Reference.fasta Reference
    bowtie2 -x Reference -U Read1.fastq -S Alignment.sam --local --np 0 --very-sensitive-local --n-ceil L,0,1 --threads 4 --score-min G,20,4 --norc

Some of the Bowtie 2 options:

--local  In this mode, Bowtie 2 does not require that the entire read aligns from one end to the other. Rather, some characters may be omitted ("soft clipped") from the ends in order to achieve the greatest possible alignment score.
--np  Sets penalty for positions where the read, reference, or both, contain an ambiguous character such as N.
--very-sensitive-local  Preset with: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
--n-ceil  Sets a function governing the maximum number of ambiguous characters (e.g. N) allowed in a read as a function of read length. For instance, L,0,1 gives f(x) = 0 + 1 * x, where x is the read length.
--threads  By default bowtie2-build is using only one thread. Increasing the number of threads will speed up the index building considerably in most cases.
--score-min  Sets a function governing the minimum alignment score needed for an alignment to be considered "valid" (i.e. good enough to report). Here G,20,8 gives f(x) = 20 + 8 * log(x), where x is the read length.
--norc  With this setting, Bowtie 2 will not attempt to align unpaired reads against the reverse-complement (Crick) reference strand.
-L  Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default 20 in --very-sensitive-local mode. For very short reads with a stretch of N's in the middle, lower this setting.


Check the manual of Bowtie2 for an explanation of all the options: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml