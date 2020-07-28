#############################################################################
#         MicroRNA-seq of Blood Serum from Field-infected Animals           #
# - Linux bioinformatics workflow for known and novel miRNAs, and isomiRs - #
#############################################################################

# Author: Carolina N. Correia 
# Last updated on: 28/07/2020

############################################
# Genome assembly preparation (Btau_5.0.1) #
############################################

# Create and enter the reference genome directory:
mkdir -p /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/source_file
cd !$

# Download the reference genome Btau_5.0.1 from NCBI into Rodeo:
nohup wget \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/205/GCA_000003205.6_Btau_5.0.1/GCA_000003205.6_Btau_5.0.1_genomic.fna.gz &

# Uncompress the reference genome Btau_5.0.1 from NCBI:
gunzip GCA_000003205.6_Btau_5.0.1_genomic.fna.gz

# Modify headers in the reference genome Btau_5.0.1 from NCBI, so that:
# 1) Sequence headers contain only >chrNUMBER (>chr1, >ch2, ...), in order to match
# the gff3 file from miRBase. Unplaced contig headers will not be changed.
perl -p -i -e \
's/^>(.*)(Bos taurus breed Hereford chromosome )(.{1,2})(\,.*)$/>chr$3/' \
GCA_000003205.6_Btau_5.0.1_genomic.fna

perl -p -i -e 's/^>(.*)(Bos taurus chromosome )(.{1,2})(\,.*)$/>chr$3/' \
GCA_000003205.6_Btau_5.0.1_genomic.fna

# Check genome fasta file:
grep '>' GCA_000003205.6_Btau_5.0.1_genomic.fna
grep '>chr' GCA_000003205.6_Btau_5.0.1_genomic.fna

# Create and enter the reference genome annotation directory:
mkdir /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/annotation_file
cd !$

# Download the miRNA annotation file from miRBase version 22.1 2018 (based on
# alternate genome Btau_5.0.1, GCA_000003205.6 from NCBI):
nohup wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/bta.gff3 &

##########################################################
# Preparation of Bos taurus miRNA sequences from miRBase #
##########################################################

# Required software is miRDeep2 v.0.1.3, consult manual/tutorial for details:
# https://github.com/rajewsky-lab/mirdeep2

# Go to working directory:
mkdir /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta
cd !$

# Download the various FASTA files for mature, high confidence mature,
# precursor(hairpin), and high confidence precursor miRNA sequences
# from miRBase (version 22.1, 2018):
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin_high_conf.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature_high_conf.fa.gz

# Uncompress the miRNA FASTA files from miRBase:
for file in `ls \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/*.gz`; \
do gunzip $file; \
done

# Extract Bos taurus mature and precursor sequences from miRBase fasta files:
extract_miRNAs.pl mature.fa bta > bta_mature.fa
extract_miRNAs.pl mature_high_conf.fa bta > bta_mature_high_conf.fa
extract_miRNAs.pl mature.fa chi,oar,ssc,eca,cel,hsa,mmu,rno > other_mature.fa
extract_miRNAs.pl hairpin.fa bta > bta_hairpin.fa
extract_miRNAs.pl hairpin_high_conf.fa bta > bta_hairpin_high_conf.fa

#######################################
# Index reference genome using Bowtie #
#######################################

# Required software is Bowtie 1.2.3, consult manual for details:
# http://bowtie-bio.sourceforge.net/manual.shtml

# Create and enter the Index reference genome directory:
mkdir -p /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/bowtie1.2.3
cd !$

# Index the reference genome UMD3.1 using Bowtie:
nohup bowtie-build --threads 20 \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/source_file/GCA_000003205.6_Btau_5.0.1_genomic.fna \
Btau_5.0.1_index &

#############################################################
# Preprocessing of miRNA-seq data using miRDeep2: mapper.pl #
#############################################################

# Required software is miRDeep2 v.2.0.0.8, consult manual/tutorial for details:
# https://www.mdc-berlin.de/8551903/en/

# Create and enter the miRDeep2 directory for mapping work:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/mapper
cd !$

# Create symbolic links to FASTQ files:
for file in \
`ls $HOME/scratch/miRNAseqValidation/fastq_sequence/tmp/*_trim.fastq`; \
do ln -s \
$file $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/`basename $file`; \
done

# Run mapper.pl in one FASTQ file to see if it's working well:
mapper.pl E10_trim.fastq -e -h -m -o 3 -l 17 -r 50 -q -v -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0/Btau_UMD3.1_multi \
-s test_collapsed.fa -t test.arf

###### RUN BOWTIW WITH BEST STRATA TO AVOID SEQ BIAS
###### <s> A comma-separated list of files containing unpaired reads
###### to be aligned, or, if -c is specified, the unpaired read sequences
###### themselves. E.g., this might be lane1.fq,lane2.fq,lane3.fq,lane4.fq, or,
###### if -c is specified, this might be GGTCATCCT,ACGGGTCGT. Reads may be a
###### mix of different lengths. If - is specified, Bowtie gets the reads from
###### the “standard in” filehandle.
# Create bash script to map miRNA reads to the reference genome:
for file in `ls *_trim.fastq`; \
do outfile=`basename $file | perl -p -e 's/_trim\.fastq//'`; \
echo "mapper.pl $file -e -h -m -o 3 -l 17 -r 50 -q -v -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0/Btau_UMD3.1_multi \
-s ${outfile}_collapsed.fa -t ${outfile}.arf" \
>> mapper.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 12 mapper.sh mapper.sh.
for script in `ls mapper.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# After mapping is finished,
# generate a master file containing mapping statistics:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/mapper.sh.0*.nohup`; \
do grep -oP "log_\d+" $file >> ./log_id.txt; \
grep "total:" $file >> ./totals.txt; \
paste log_id.txt totals.txt > ./stats.txt; \
done

# Using awk, keep only columns of interest from stats.txt:
awk '{print $1, $3, $4, $5, $6, $7}' stats.txt > stats2.txt

# Adding header to stats2.txt and save mapper_stats.txt:
echo -e "Log_Id Input_reads Total_Mapped_Reads Total_Unmapped_Reads \
Percentage_Mapped_Reads Percentage_Unmapped_Reads" | cat - ./stats2.txt > \
./mapper_stats.txt

# Delete unnecessary files:
rm -r log_id.txt totals.txt stats.txt stats2.txt

# Get filenames from log ids:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/mapper_logs/mapper.log_*`; \
do basename $file | perl -p -e 's/mapper\.(log_\d+)_\d+/$1/' >> log.txt; \
grep -oP "E\d*_collapsed.fa" $file | perl -p -e 's/(\d+_\w*\d+)_collapsed\.fa/$1/' \
>> sample.txt; \
paste log.txt sample.txt > ./id_sample.txt; \
done

# Delete unnecessary files:
rm -r log.txt sample.txt

# Transfer mapper_stats.txt and id_sample.txt to laptop via SCP.

################################################################
# Quantification of known miRNAs using miRDeep2: quantifier.pl #
################################################################

# Create and enter the working directory:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier
cd !$

# Run quantifier.pl in one FASTA file to see if it's working well:
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10_collapsed.fa -t bta

# Create a shell script to quantify the mapped mature miRNAs:
# [You will need the mature and precursor (hairpin) miRNA FASTA files for
# Bos taurus sequences only. Please refer to the 
# BioValidation-miRNA-seq_QC_filter.sh script]
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/$outfile; \
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $file -t bta" >> quantifier.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 quantifier.sh quantifier.sh.
for script in `ls quantifier.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence
cd !$
# Create a shell script to quantify the mapped high confidence mature and
# precursor (hairpin) miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/$outfile; \
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
-r $file -t bta" >> quantifier_high_conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 quantifier_high_conf.sh quantifier_high_conf.sh.
for script in `ls quantifier_high_conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read counts files from regular and high confidence miRNAs for
# easier transfer into laptop using WinSCP:
mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep2/regular
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep2/regular/${outfile}_expressed.csv; \
done

mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep2/high_conf
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep2/high_conf/${outfile}_expressed_high_conf.csv; \
done

# After transferring into laptop, delete duplicate count folders:
cd $HOME/scratch/miRNAseqValidation/Counts
rm -r mirdeep2/

########################################################################
# Identification of known and novel miRNAs using miRDeep2: miRDeep2.pl #
######################################################################## 

# Create and enter working directory:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep
cd !$

# Copy and modify the required fasta files since miRdeep2 software
# requires no space in headers and no characters other than acgtunACGTUN in
# the sequences:
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa \
./Btau_UMD3.1_multi.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
./bta_mature-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
./high_conf_mature.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/other_mature-miRNA.fa \
./other_mature-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
./bta_hairpin-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
./high_conf_hairpin.fa

perl -p -i -e 's/^(>.*?)\s.*$/$1/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*$)/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_mature.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_hairpin-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) (.*?) (.*?) (.*?) .*$/$1_$2/' \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_hairpin.fa
sed -e '/^[^>]/s/[^ATGCatgc]/N/g' high_conf_hairpin.fa > high_conf_hairpin_clean.fa

# Run mirdeep.pl in one FASTA file to see if it's working well:
miRDeep2.pl $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10_collapsed.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10.arf \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_hairpin-miRNA.fa -t Cow

# Create bash script for identification and quantification of known and
# novel miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/$outfile; \
miRDeep2.pl $file \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mapper/${outfile}.arf \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/bta_hairpin-miRNA.fa -t Cow" \
>> miRdeep2.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 miRdeep2.sh miRdeep2.sh.
for script in `ls miRdeep2.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence
cd !$

# Create bash script for identification and quantification of
# high confidence miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo \
"mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence/$outfile; \
miRDeep2.pl $file \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mapper/${outfile}.arf \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_mature.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/other_mature-miRNA.fa \
$HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_conf_hairpin_clean.fa \
-t Cow" \
>> miRdeep2_high-conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 miRdeep2_high-conf.sh miRdeep2_high-conf.sh.
for script in `ls miRdeep2_high-conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read counts files from regular and high confidence miRNAs for
# transfering into laptop using WinSCP:
mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep.pl/regular
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*mirdeep\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep.pl/regular/${outfile}_exp_mirdeep.csv; \
done

mkdir -p $HOME/scratch/miRNAseqValidation/Counts/mirdeep.pl/high_conf
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/high_confidence/E* \
-name result*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*mirdeep\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
$HOME/scratch/miRNAseqValidation/Counts/mirdeep.pl/high_conf/${outfile}_novel_hc.csv; \
done

# Collet data from predicted novel miRNAs into one folder:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/novel_miRNAs
cd !$
for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/E* \
-name result_*.csv`; \
do outfile=`basename $file | perl -p -e 's/result_.+(\.csv)/novel_miRNAs$1/'`; \
sample=`dirname $file | perl -p -e 's/.+mirdeep\/(E\d+).*/$1/'`; \
cp $file ./$sample\_$outfile; \
done

# Clean data to keep only info about novel miRNAs:
for file in `ls $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/novel_miRNAs/E*_novel_miRNAs.csv`; \
do outfile=`basename $file | perl -p -e 's/(E\d+)_novel_miRNAs.csv/$1_clean.csv/'`; \
sed '/provisional id/,$!d' $file | sed '/mature miRBase miRNAs/,$d' > $outfile; \
done

# Transfer cleaned files to laptop with SCP.

##############################
# isomiR count summarisation #
##############################

# Process miRDeep2 output files to collect all isoforms read counts files:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/isomiR_counts
cd !$

for file in `find $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep \
-name "miRBase.mrd"`; \
do outfile=`echo $file | perl -p -e 's/.+mirdeep\/(E\\d+).+/$1_isomiR.txt/'`; \
echo "perl $HOME/storage/Scripts/Get_isomiR_count.pl -mrd $file \
-output $HOME/scratch/miRNAseqValidation/mirdeep2/mirdeep/isomiR_counts/$outfile" \
>> isomiR_summarisation.sh; \
done

# Run the script on Stampede
chmod 755 isomiR_summarisation.sh
nohup ./isomiR_summarisation.sh &

# Transfer isomiR count files to laptop with SCP.

########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent sense genes analyses were performed using the R statistical
# and the edgeR package. Please refer to file:
# BioValidation-miRNAseq_edgeR_pipeline.R

