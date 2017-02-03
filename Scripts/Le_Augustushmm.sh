#!/bin/bash

#SBATCH -A XXXXXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH --mail-user something@email.com
#SBATCH --mail-type=ALL

# AUGUSTUS HMM training for Letharia lupina
#############################################################################
# Pipeline to create an HMM training file for Letharia lupina (pure culture), based on 
# http://bioinf.uni-greifswald.de/augustus/binaries/README.autoAug

# See also Pa_Augustushmm.sh

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2016/08/16
# ++++++++++++++++++++++++++++++++++++++++++++++

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
module load blat/34 # 35 might be too new? (see /sw/apps/bioinfo/Scipio/Scipio-1.4-install-README.md) 
module load augustus/3.2.2  &> /dev/null # It will tell me about 'module help augustus/3.2.2'
module load Scipio/1.4  &> /dev/null # It will tell me about 'module help augustus/3.2.2'
module load BioPerl/1.6.924_Perl5.18.4 # for scipio
module load BioPerl/1.6.1 # The only one that works with scipio
# ***************** 

# --- Variables ---
pathAUGUSTUS="somepath/Augustus/augustus-3.2.2"
prots="somepath/SNAP/L.lupinapure.scaffolds_ed.SNAP_sample1000.fa"
# -----------------

# ==================== autoAug.pl =======================
# ---------------------------------------------------------
export PATH=$PATH:/proj/b2015200/Pa_scripts/ThirdParty # To make pslCDnaFilter available to autoAug.pl (this scripts seems to belong to The UCSC Genome Browser, see http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
export PATH=$PATH:/sw/apps/bioinfo/Scipio/1.4/milou/ # To make pslCDnaFilter available to scipio.pl 

ln -sf $pathAUGUSTUS/scripts/autoAug.pl
ln -sf $pathAUGUSTUS/scripts/autoAugPred.pl
ln -sf $pathAUGUSTUS/scripts/helpMod.pm
cp $pathAUGUSTUS/scripts/autoAugTrain.pl .

# Comment the check for errors in scipio call
sed -i 's/system("$cmdString")==0 or die("Program/system("$cmdString")==0; # or die("Program/' autoAugTrain.pl
# Correct the faulty call to yaml2gff.1.4.pl
sed -i -e 's;cat scipio.yaml | ;;' -e 's;> scipio.gff; scipio.yaml > scipio.gff;' autoAugTrain.pl

# **** Make a local copy of the AUGUSTUS_CONFIG directory to write my new species ****
echo OLD AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH
oldAUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH
# Copy the whole thing
source $AUGUSTUS_CONFIG_COPY
# To be able to erase it later
chmod -R 777 augustus_config # $ rm -r augustus_config

# Make a link to the actual binaries so autoAug.pl can use them
ln -s $oldAUGUSTUS_CONFIG_PATH/../bin
ln -s $oldAUGUSTUS_CONFIG_PATH/../scripts
ln -s $oldAUGUSTUS_CONFIG_PATH/../src

# Trick Augustus to use the local installation
AUGUSTUS_CONFIG_PATH=$PWD/augustus_config

# With SNAP Letharia's proteins
perl autoAug.pl -g $GENOME --trainingset=$prots --species=$sampleID -v -v -v --singleCPU --maxIntronLen=1000 --useexisting #--optrounds=0

# -t traingenes.gff is the set of training gene structures in GFF format (e.g.
# as produced by PASA). Alternatively you can use -t traingenes.gb if you have
# a training set in Genbank format. This has the advantage that the training
# genes can be based on a different DNA sequence than genome.fa
# --trainingset=genesfile             genesfile contains training genes in Genbank, GFF or protein FASTA format
# --cdna=cdna.fa                      a fasta file with cDNA sequences (ESTs, mRNA)
# --useexisting                       use and change the present config and parameter files if they exist for 'species'
# --maxIntronLen=n                    maximal length of an intron as used by PASA and BLAT, not by AUGUSTUS (default 100000) # I don't think introns that long are common in fungi

# For this pipeline you only need a training set with coding regions, even
# if you want AUGUSTUS to predict the UTRs. It will construct a training set of
# UTRs from EST alignments.

# In a first test run we recommend to use the option --optrounds=0 to skip the
# optimization, which usually takes most time. 

##################### About autoAug.pl #####################

# The purpose of this pipeline script is to run the AUGUSTUS training process and gene
# prediction algorithm automatically on a given eukaryotic genome with
# available cDNA evidence. The complete process contains the following steps.
# These are encapsulated in this automatic training pipeline.

# -- Construct an initial training set of genes (protein coding regions thereof), e.g. using PASA
# -- Train AUGUSTUS to predict coding regions (i.e. without Untranslated Regions (UTR)) using the initial training set. 
# -- Predict coding regions of genes in your genome (both ab initio and using hints from cDNA)
# -- Train the UTR model of AUGUSTUS using EST-supported genes and EST alignments.
# -- Predict genes including UTRs in your genome using hints from cDNA.

# The first step is optional and can be replaced by preparing your own training set of genes. For a schema of the
# pipeline look at the poster in the doc folder. 
