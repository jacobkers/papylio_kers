# When running Linux for the first time on Windows 10:
# 1. In the search field of the start menu search for "Turn Windows features on and off"
# 2. Scroll down to "Windows Subsystem for Linux", check the box and click the OK button
# 3. In settings go to "Update & Security", select "For developers" and turn on the developer mode
# 4. In the Microsoft store search for and install "Ubuntu 22.04.2 LTS"
# 5. Restart the computer
# 6. Run the code in the sequencing_setup file in Ubuntu

# Steps to follow:
# 1. Copy everything from the sequence_analysis folder to the folder with your sequencing data (not on network location)
# 2. Unzip the fastq file and rename it to Read1
# 3. Adjust the sequences in the reference file
# 4. Open Ubuntu
# 5. Navigate to the folder with your sequencing data 
# 6. Run the lines below

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sequence_analysis

bowtie2-build Reference.fasta Reference
bowtie2 -x Reference -U Read1.fastq -S Alignment.sam --local --np 0 --very-sensitive-local --n-ceil L,0,1 --threads 4 --score-min G,20,4 --norc
# Check the manual for all the options: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml