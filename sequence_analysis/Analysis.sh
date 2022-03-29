# CONDA_BASE=$(conda info --base)
# source $CONDA_BASE/etc/profile.d/conda.sh
# CONDA=$(conda info --base)"/condabin/conda"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sequence_analysis

# bowtie2-build Reference.fasta Reference

# bowtie2 -x Reference -U *R1_001.fastq -S Alignment.sam --local --np 0 --very-sensitive-local --n-ceil L,0,1 --threads 4 --score-min G,20,4 --norc

cat *R1_001.fastq > Read1.fastq

bwa index Reference.fasta
bwa mem Reference.fasta Read1.fastq -k 5 -T 15 -Y -L 10 -t 12 > Alignment.sam

# samtools view -o Alignment.bam Alignment.sam  # use the -c option to just count alignment records
# samtools sort Alignment.bam -o Alignment.sorted.bam
# samtools index Alignment.sorted.bam

while IFS="" read -r line || [ -n "$line" ]; do
  if  [[ $line == \>* ]];
  then
    name=${line#">"}
    name=$(tr -dc '[[:print:]]' <<< "$name")
    read sequence
    sequence=$(tr -dc '[[:print:]]' <<< "$sequence")
    echo Working on $name
    samtools view Alignment.sam | awk -v refname="$name" -v refseq="$sequence" -f ./Align_string.awk > sequencing_data_$name.csv
  fi
done <./Reference.fasta

# To extract specific positions
# gawk -v pos_string="31 32 56 57 82 83 108 109" -f ./Extract_section.awk -i inplace sequencing_data_HJ_general.csv

read -p 'Press Enter to continue...' var

# samtools mpileup -f Reference.fasta Alignment.sorted.bam -l Position\ list.BED -o test.txt --output-extra QNAME -O

#files=$(find . -name "*.fastq")
#files=$(echo $files | tr "\n" "," )

# bedtools intersect -a Alignment.sorted.bam -b Position\ list.BED -wa > test.bam

# awk 'NR%4==2 {print substr($0, 0, 4)}' Empty_S1_L001_R1_001.fastq | head

# samtools view Alignment.sam | awk '$3=="HJ_general" {print substr($10, 0, 4)}' > indexL1.txt

# samtools view Alignment.sam | awk '$3=="HJ_general" { split("31 32 56 57 82 83 108 109", a, " "); print a; for (i in a) {$3 = $3 substr($10, a[i], 1); print a[i]; print substr($10, a[i], 1)}; print $3}' | head

# samtools view Alignment.sorted.bam | awk 'BEGIN {print "Tile\tx\ty"} $3=="MapSeq" {split($1,a,":"); print a[5] "\t" a[6] "\t" a[7]}' > mapping_coordinates.txt

# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example
