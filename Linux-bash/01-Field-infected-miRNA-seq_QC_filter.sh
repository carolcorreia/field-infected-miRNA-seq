#############################################################
#  MicroRNA-seq of Blood Serum from Field-infected Animals  #
# - Linux bioinformatics workflow for data pre-processing - #
#############################################################

# Author: Carolina N. Correia 
# Last updated on: 28/07/2020

################################
# Download and files check sum #
################################

# Create and enter working directory:
mkdir $HOME/storage/miRNAseqValidation
cd !$

# Download data from MSU server into Stampede server

# Change directory permissions to read and execute only:
chmod -R 555 $HOME/storage/miRNAseqValidation

# Create md5sum working folder:
mkdir -p $HOME/scratch/miRNAseqValidation/md5check
cd !$

# Create bash script to perform md5sum check:
for file in `find $HOME/storage/miRNAseqValidation -name md5.txt`; \
do echo "cd `dirname $file` && md5sum -c `basename $file` >> \
$HOME/scratch/miRNAseqValidation/md5check/md5_UCD.txt" >> md5sum.sh; \
done

# Run script on Stampede:
chmod 755 ./md5sum.sh
nohup ./md5sum.sh &

# Check that all files passed the check:
grep -c 'OK' md5_UCD.txt

# In 2020, miRNA-seq files were moved to Rodeo server and their
# integrity re-checked.

for file in `find /home/workspace/ccorreia/miRNASeq_field/fastq -name md5.txt`; \
do echo "cd `dirname $file` && md5sum -c `basename $file` >> \
/home/workspace/ccorreia/miRNASeq_field/md5_UCD.txt" >> md5sum.sh; \
done

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering
cd !$

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/workspace/ccorreia/miRNASeq_field/fastq \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering $file" \
>> fastqc.sh; done

# Split and run all scripts on Stampede:
split -d -l 24 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Deleted all the HTML files:
rm -r *.html

# Check if all the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering/tmp; \
done

for file in \
`find /home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> summary_pre-filtering.txt; \
done

for file in \
`find /home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Check sequence quality:
grep 'Per base sequence quality' summary_pre-filtering.txt >> seq_quality.txt
wc -l seq_quality.txt # 48 lines
grep PASS seq_quality.txt | wc -l # 48 lines
grep WARN seq_quality.txt | wc -l # 0 line

# Check if sequence contain adapters:
grep 'Adapter Content' summary_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt # 48 lines
grep PASS adapter_content.txt | wc -l # 0 lines
grep WARN adapter_content.txt | wc -l # 0 lines
grep FAIL adapter_content.txt | wc -l # 48 lines

# Check sequence length:
grep 'Sequence length' basic_stats_pre-filtering.txt >> seq_length.txt
wc -l seq_length.txt # 48 lines
less seq_length.txt # 50 for all

# Check per base sequence quality
grep 'Per base sequence quality' summary_pre-filtering.txt >> base_quality.txt
wc -l base_quality.txt # 48 lines
grep PASS base_quality.txt | wc -l # 48 lines
grep WARN base_quality.txt | wc -l # 0 lines

# Transfer compressed folders to personal laptop via SCP
# and check HTML reports:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/quality_check/pre-filtering/tmp .

# Remove temporary folder and its files:
rm -r tmp

#############################################
# Trimming of adapter sequence within reads #
#############################################

# Required software is cutadapt (version 2.10).
# Consult manual for details:
# https://cutadapt.readthedocs.io/en/stable/guide.html
# Sequence used for adapter trimming comes from :
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-14.pdf

# Create and enter working directory:
mkdir /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed
cd !$

# Create bash script to trim the Illumina RNA 3’ Adapter (RA3) of each
# FASTQ file while keeping the sequencing lane information:
for file in `find /home/workspace/ccorreia/miRNASeq_field/fastq \
-name *fastq.gz`; \
do outfile=`basename $file | perl -p -e 's/.fastq.gz//'`; \
echo "cutadapt -a TGGAATTCTCGGGTGCCAAGG --overlap 10 --match-read-wildcards \
--cores=20 --minimum-length 17 --discard-untrimmed \
-o ./${outfile}_trim.fastq.gz $file" \
>> cutadapt.sh; \
done

# Split and run scripts on Stampede:
split -d -l 24 cutadapt.sh cutadapt.sh.
for script in `ls cutadapt.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Finished' cutadapt.sh.00.nohup
grep -c 'Finished' cutadapt.sh.01.nohup

# Generate a master file containing cutadapt trimming stats results:
for file in `ls /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/*.nohup`; \
do grep -oP "E\d+_\w+_L00(5|6)\_R1_001.fastq.gz" $file \
>> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/filename.txt; \
grep "Total reads processed\:" $file | \
perl -p -e 's/\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)/$1/' \
>> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/processed.txt; \
grep "Reads with adapters\:" $file | \
perl -p -e 's/\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)\s.*/$1/' \
>> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/reads_with_adapters.txt; \
grep "Reads that were too short\:" $file | \
perl -p -e 's/\w*\s\w*\s\w*\s\w*\s\w*\:\s*(\d*.\d*.\d*)\s.*/$1/' \
>> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/short.txt; \
grep "Reads written (passing filters)\:" $file | \
perl -p -e 's/\w*\s\w*\s.\w*\s\w*.\:\s*(\d*.\d*.\d*)\s.*/$1/' \
>> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/trimmed.txt; \
paste filename.txt processed.txt reads_with_adapters.txt short.txt trimmed.txt \
> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/trimmed_stats.txt; \
done

echo -e "Sample\tTotal reads processed\tReads with \
adapters\tReads that were too short\tReads written (passing filters)\t" \
| cat - /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/trimmed_stats.txt \
> /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/trimming_stats.txt

rm -r filename.txt processed.txt reads_with_adapters.txt \
short.txt trimmed.txt trimmed_stats.txt

# Transfer trimming stats file to personal laptop via SCP:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/fastq_trimmed/trimming_stats.txt .

###############################################
# FastQC quality check of trimmed FASTQ files #
###############################################

# Required software is FastQC v0.11.8, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter working directory:
mkdir /home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering
cd !$

# Create bash script to perform FastQC quality check on all trim.fastq.gz files:
for file in `find /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed \
-name *trim.fastq.gz`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering $file" \
>> fastqc.sh; done

# Split and run all scripts on Stampede:
split -d -l 24 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Deleted all the HTML files:
rm -r *.html

# Check if all the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering/tmp; \
done

for file in \
`find /home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering/tmp \
-name summary.txt`; do more $file >> summary_post-filtering.txt; \
done

for file in \
`find /home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Check sequence quality:
grep 'Per base sequence quality' summary_post-filtering.txt >> seq_quality.txt
wc -l seq_quality.txt # 48 lines
grep PASS seq_quality.txt | wc -l # 48 lines
grep WARN seq_quality.txt | wc -l # 0 line

# Check if sequence contain adapters:
grep 'Adapter Content' summary_post-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt # 48 lines
grep PASS adapter_content.txt | wc -l # 48 lines
grep WARN adapter_content.txt | wc -l # 0 lines
grep FAIL adapter_content.txt | wc -l # 0 lines

# Check sequence length:
grep 'Sequence length' basic_stats_post-filtering.txt >> seq_length.txt
wc -l seq_length.txt # 48 lines
less seq_length.txt # 17-41 for all

# Check per base sequence quality
grep 'Per base sequence quality' summary_post-filtering.txt >> base_quality.txt
wc -l base_quality.txt # 48 lines
grep PASS base_quality.txt | wc -l # 48 lines
grep WARN base_quality.txt | wc -l # 0 lines

# Check per base N content
grep 'Per base N content' summary_post-filtering.txt >> N_content.txt
wc -l N_content.txt # 48 lines
grep PASS N_content.txt | wc -l # 48 lines
grep WARN N_content.txt | wc -l # 0 lines

# Transfer compressed folders to personal laptop via SCP
# and check HTML reports:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/quality_check/post_filtering/tmp .

# Remove temporary folder and its files:
rm -r tmp
