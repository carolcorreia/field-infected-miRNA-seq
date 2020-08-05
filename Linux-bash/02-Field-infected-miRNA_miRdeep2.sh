#############################################################################
#         MicroRNA-seq of Blood Serum from Field-infected Animals           #
# - Linux bioinformatics workflow for known and novel miRNAs, and isomiRs - #
#############################################################################

# Author: Carolina N. Correia 
# Last updated on: 05/08/2020

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
# the gff3 file from miRBase. Unplaced contig headers will only contain their ID.
perl -p -i -e \
's/^>(.*)(Bos taurus breed Hereford chromosome )(.{1,2})(\,.*)$/>chr$3/' \
GCA_000003205.6_Btau_5.0.1_genomic.fna

perl -p -i -e 's/^>(.*)(Bos taurus chromosome )(.{1,2})(\,.*)$/>chr$3/' \
GCA_000003205.6_Btau_5.0.1_genomic.fna

perl -p -i -e 's/^>(.*)(Bos taurus breed Hereford UnplacedContig)(.*)$/>$1/' \
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

#############################
# Index genome using Bowtie #
#############################

# Required software is Bowtie 1.2.3, consult manual for details:
# http://bowtie-bio.sourceforge.net/manual.shtml

# Create and enter the Index reference genome directory:
mkdir -p /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/bowtie1.2.3
cd !$

# Index the reference genome UMD3.1 using Bowtie:
nohup bowtie-build --threads 20 \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/source_file/GCA_000003205.6_Btau_5.0.1_genomic.fna \
Btau_5.0.1_index &

########################################################################
# Merge and uncompress miRNA-seq FASTQ files to be used with miRDeep2  #
########################################################################

# Create and enter temporary working directory:
mkdir /home/workspace/ccorreia/miRNASeq_field/fastq_merged
cd !$

# Create bash script to uncompress and merge trim.fast.gz from lanes
# 005 and 006 for each library:
for file in `find /home/workspace/ccorreia/miRNASeq_field/fastq_trimmed \
-name *_L005_R1_001_trim.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_L005_)/_L006_/'`; \
sample=`basename $file | perl -p -e 's/(E\d+_)\w+_L\d+_R\d_\d*_trim.fastq.gz/$1/'`; \
echo "zcat $file $file2 > ${sample}trim.fastq" \
>> uncompress_merge.sh; \
done

# Run script:
chmod 755 uncompress_merge.sh
nohup ./uncompress_merge.sh &

#############################
# Map miRNA reads to genome #
#############################

# Required software is miRDeep2 v.0.1.3, consult manual/tutorial for details:
# https://www.mdc-berlin.de/8551903/en/

# Create and enter the miRDeep2 directory for mapping work:
mkdir -p /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper
cd !$

# Create bash script to map miRNA reads to the reference genome:
for file in `find /home/workspace/ccorreia/miRNASeq_field/fastq_merged \
-name *_trim.fastq`; \
do outfile=`basename $file | perl -p -e 's/_trim\.fastq//'`; \
echo "mapper.pl $file -e -h -i -l 18 -m -q -r 50 -v -o 20 \
-p /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/bowtie1.2.3/Btau_5.0.1_index \
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
`ls /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/mapper.sh.0*.nohup`; \
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
`ls /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/mapper_logs/mapper.log_*`; \
do basename $file | perl -p -e 's/mapper\.(log_\d+)_\d+/$1/' >> log.txt; \
grep -oP "E\d*_collapsed.fa" $file | perl -p -e 's/(\d+_\w*\d+)_collapsed\.fa/$1/' \
>> sample.txt; \
paste log.txt sample.txt > ./id_sample.txt; \
done

# Delete unnecessary files:
rm -r log.txt sample.txt

# Transfer mapper_stats.txt and id_sample.txt to laptop via SCP:
scp \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/mapper_stats.txt .

scp \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/id_sample.txt .

################################
# Quantify known mature miRNAs #
################################

# Create and enter the working directory:
mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier
cd !$

# Create a shell script to quantify miRNAs that mapped to the genome:
for file in \
`ls /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/$outfile; \
cd /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/$outfile; \
quantifier.pl -W -t bta \
-p /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_hairpin.fa \
-m /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_mature.fa \
-r $file" >> quantifier.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 12 quantifier.sh quantifier.sh.
for script in `ls quantifier.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all mature miRNA read counts:
mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/mature_counts
cd !$

for file in \
`find /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/(.*)\/.*$/$1/'`; \
cp $file \
/home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/mature_counts/${outfile}_expressed.csv; \
done

# Transfer count files to laptop via SCP:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/mature_counts .

################################################
# Quantify known high confidence mature miRNAs #
################################################

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence
cd !$

# Create a shell script to quantify high confidence miRNAs that mapped to the genome:
for file in \
`ls /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence/$outfile; \
cd /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence/$outfile; \
quantifier.pl -W -t bta \
-p /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_mature_high_conf.fa \
-m /home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_hairpin_high_conf.fa \
-r $file" >> quantifier_high_conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 12 quantifier_high_conf.sh quantifier_high_conf.sh.
for script in `ls quantifier_high_conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all high confidence mature miRNA read counts:
mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence/mat_high_conf_counts
cd !$

for file in `find /home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence/E* \
-name miRNAs_expressed_all_samples*.csv`; \
do outfile=`echo $file | perl -p -e 's/^.*quantifier\/.*\/(.*)\/.*$/$1/'`; \
cp $file \
/home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence/mat_high_conf_counts/${outfile}_expressed_high_conf.csv; \
done

# Transfer count files to laptop via SCP:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/miRNASeq_field/mirdeep2/quantifier/high_confidence/mat_high_conf_counts .

###################################################
# Identification of known and novel mature miRNAs #
###################################################

# Create and enter working directory:
mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mirdeep2.pl
cd !$

# Remove all spaces for miRDeep2.pl script:
sed -i 's/[[:blank:]]//g' \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/source_file/GCA_000003205.6_Btau_5.0.1_genomic.fna

# Run mirdeep.pl in one FASTA file to see if it's working well:
miRDeep2.pl /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/E10_collapsed.fa \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/source_file/GCA_000003205.6_Btau_5.0.1_genomic.fna \
/home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/E10.arf \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_mature.fa \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/other_mature.fa \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_hairpin.fa \
-t Cow

# Create bash script for identification and quantification of known and
# novel mature miRNAs:
for file in \
`ls /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mirdeep2.pl/$outfile; \
cd /home/workspace/ccorreia/miRNASeq_field/mirdeep2/mirdeep2.pl/$outfile; \
miRDeep2.pl $file \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/source_file/GCA_000003205.6_Btau_5.0.1_genomic.fna \
/home/workspace/ccorreia/miRNASeq_field/mirdeep2/mapper/${outfile}.arf \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_mature.fa \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/other_mature.fa \
/home/workspace/genomes/bostaurus/Btau_5.0.1_NCBI/miRBase_fasta/bta_hairpin.fa \
-t Cow" \
>> miRdeep2.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 12 miRdeep2.sh miRdeep2.sh.
for script in `ls miRdeep2.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

###################################################################
# Identification of known and novel mature high confidence miRNAs #
###################################################################

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

