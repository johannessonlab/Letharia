#!/bin/bash

# SNAP HMM training for Letharia lupina
#############################################################################
# Pipeline to create an HMM training file for Letharia lupina, based on 
# http://gmod.org/wiki/MAKER_Tutorial#Training_ab_initio_Gene_Predictors
# https://biowize.wordpress.com/2012/06/01/training-the-snap-ab-initio-gene-predictor/
# Check also the 00README of SNAP and Campbell et al. (2014)

# The script first calls MAKER with the Trinity transcripts of Letharia lupina
# (pure culture). MAKER creates a gff file that is then transformed into a ZFF
# file for SNAP. The ZFF file is processed and then used to create an HMM
# file. Then this HMM is used for a second round of training.

# The script was based on Pa_SNAPhmm.sh

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2016/12/07
# ++++++++++++++++++++++++++++++++++++++++++++++

#SBATCH -A XXXXXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH --mail-user something@email.com
#SBATCH --mail-type=ALL

GENOME=$1

if [[ -z "$GENOME" ]] ; then
        echo 'Must specify a genome fasta file to be annotated'
        echo "Usage: $0 genome.fa [ sampleID ]"
        exit 1
fi

# A name for the sample
basegenome=$(basename ${GENOME%.fa*})
sampleID=${2:-$basegenome}

# Load UPPMAX modules
# ***************** 
module load bioinfo-tools
module load maker/2.31.8 &> /dev/null
# ***************** 

# --- Variables ---
pathMAKER="/pica/sw/apps/bioinfo/maker/2.31.8/milou/bin"
zffcleaner="path2scripts/zffcleaner.py"
Transcripts="L.lupinapure_Trinity.fasta"
# -----------------

# Prepare the local copies of transcripts and proteins with better names
cp $Transcripts transcripts.fa
sed -i -E 's/(>TRINITY_[0-9A-Z_a-z]+)(.*)/\1/' transcripts.fa

# ==================== SNAP ROUND 1 =======================
# ---------------------------------------------------------
echo "*** Training SNAP round 1 ***"
mkdir -p SNAP1 
cd SNAP1

### PREPARE MAKER CONFIGURATION FILES ###
# Produce a sample configuration file
maker -CTL

# Rename them to match the sample
mv maker_bopts.ctl $sampleID'_bopts.ctl'
mv maker_exe.ctl $sampleID'_exe.ctl'
mv maker_opts.ctl $sampleID'_opts.ctl'

# Use the latest versions in UPPMAX
sed -i 's|blast/2.2.31+|blast/2.4.0+|g' $sampleID'_exe.ctl'

# *** Settings for Letharia data ***
# Set the current genome
sed -i "s|^genome=|genome=$GENOME|g" $sampleID'_opts.ctl'

# #-----EST Evidence
# ESTs or assembled mRNA-seq in fasta format
sed -i "s|^est=|est=../transcripts.fa|g" $sampleID'_opts.ctl' # Add that to the opts file

# #-----Gene Prediction
# infer gene predictions directly from ESTs, 1 = yes, 0 = no
sed -i "s|^est2genome=0|est2genome=1|g" $sampleID'_opts.ctl' # Add that to the opts file

echo "*** Running MAKER ***"
maker $sampleID'_opts.ctl' $sampleID'_bopts.ctl' $sampleID'_exe.ctl'

# ---------------------------------------------------------

# Collect all gff files into a single directory
mkdir -p MAKERgffs 
cd MAKERgffs
gff3_merge -d ../*.maker.output/*_index.log 
cd ..

echo "*** Retraining SNAP ***"

# Transform the gff files into the non-standard ZFF file for SNAP
maker2zff $(ls -1 MAKERgffs/*.gff) # The output are genome.ann and genome.dna

# In order for a gene model to be considered suitable for training, it has to
# pass several quality filters imposed by the maker2zff script. By default, a
# gene model must have half of its splice sites confirmed by an EST/mRNA-seq
# alignment; half of its exons must overlap an EST/mRNA-seq alignment; and its
# annotation edit distance must be less than 0.5. All of these criteria can be
# modified on the command line. (see Campbell et al. 2014)

# The basic steps for training SNAP are first to filter the input gene models,
# then capture genomic sequence immediately surrounding each model locus, and
# finally uses those captured segments to produce the HMM.

# Get some stats about the GFF file
fathom genome.ann genome.dna -gene-stats &> $basegenome"_zff.stats" # It includes the stderr
# Verify that the genes have no obvious errors
fathom genome.ann genome.dna -validate &> $basegenome"_zff.errors"
tail -n1 $basegenome"_zff.errors" # Tell me the summary

# The SNAP manual recommends:
# "You may find some errors and warnings. Check these out in some kind of genome
# browser and remove those that are real errors."
# I decided to remove them instead, or it would take way too long.

$zffcleaner genome.ann $basegenome"_zff.errors" # output is genome_clean.ann

echo "Preparing data for forge ..."
# Next, break up the sequences into fragments with one gene per sequence with
# the following command:

fathom -categorize 1000 genome_clean.ann genome.dna

# There will be up to 1000 bp on either side of the genes. You will find
# several new files.

    # alt.ann, alt.dna (genes with alternative splicing)
    # err.ann, err.dna (genes that have errors)
    # olp.ann, olp.dna (genes that overlap other genes)
    # wrn.ann, wrn.dna (genes with warnings)
    # uni.ann, uni.dna (single gene per sequence)

# Convert the uni genes to plus stranded with the command:

fathom -export 1000 -plus uni.ann uni.dna

# You will find 4 new files:

#     export.aa   proteins corresponding to each gene
#     export.ann  gene structure on the plus strand
#     export.dna  DNA of the plus strand
#     export.tx   transcripts for each gene

# The parameter estimation program, forge, creates a lot of files.
echo "Starting parameter estimation ..."
mkdir params
cd params
forge ../export.ann ../export.dna
cd ..

# Last is to build an HMM.
echo "Building the first HMM file ..."
hmm-assembler.pl Letharia_$sampleID params > $basegenome.hmm

cd .. # Leave SNAP1

# ==================== SNAP ROUND 2 =======================
# ---------------------------------------------------------
# We need to run MAKER again with the new HMM file we just built for SNAP.
echo "*** Training SNAP round 2 ***"
mkdir -p SNAP2 
cd SNAP2

# *** Prepare opts file ***
cp ../SNAP1/$sampleID'_opts.ctl' $sampleID'_opts2.ctl'

# #-----Gene Prediction
# infer gene predictions directly from ESTs, 1 = yes, 0 = no
sed -i "s|^est2genome=1|est2genome=0|g" $sampleID'_opts2.ctl' # Add that to the opts file
# SNAP HMM file
sed -i "s|^snaphmm=|snaphmm=../SNAP1/$basegenome.hmm|g" $sampleID'_opts2.ctl'

echo "*** Running MAKER again ***"
maker $sampleID'_opts2.ctl' ../SNAP1/$sampleID'_bopts.ctl' ../SNAP1/$sampleID'_exe.ctl' # With new options

# Collect all gff files into a single directory
mkdir -p MAKERgffs 
cd MAKERgffs
gff3_merge -d ../*.maker.output/*_index.log 
# cp $basegenome*/$basegenome*datastore/*/*/$sampleID*/$sampleID*.gff MAKERgffs
cd ..

echo "*** Retraining SNAP ***"
# Transform the gff files into a the non-standard ZFF file for SNAP
maker2zff $(ls -1 MAKERgffs/*.gff) # The output are genome.ann and genome.dna

# Get some stats about the GFF file
fathom genome.ann genome.dna -gene-stats &> $basegenome"_zff.stats" # It includes the stderr
# Verify that the genes have no obvious errors
fathom genome.ann genome.dna -validate &> $basegenome"_zff.errors"
tail -n1 $basegenome"_zff.errors"

$zffcleaner genome.ann $basegenome"_zff.errors" # output is genome_clean.ann

echo "Preparing data for forge ..."
fathom -categorize 1000 genome_clean.ann genome.dna

# Convert the uni genes to plus stranded
fathom -export 1000 -plus uni.ann uni.dna

echo "Starting parameter estimation ..."
mkdir params
cd params
forge ../export.ann ../export.dna
cd ..

# The new HMM
echo "Building the first HMM file ..."
hmm-assembler.pl Letharia_$sampleID params > $basegenome.hmm

cd .. # Leave SNAP2

# *** CLEANING ***
rm transcripts.fa

# Remove folder that contains FASTA indexes and BLAST database files created from the input EST, protein, and repeat databases.
# rm -r SNAP*/$basegenome".maker.output"/mpi_blastdb
