#!/bin/bash

#SBATCH -A XXXXXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-user something@email.com
#SBATCH --mail-type=ALL

# Make ab initio gene predictions
#############################################################################
# Script to run a number of ab initio gene predictors on a multifasta file
# with previously produced training files.

# - SNAP
# - Augustus
# - GeneMark
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2016/12/23
# ++++++++++++++++++++++++++++++++++++++++++++++

GENOME=$1
VERSION="1.4"

if [[ -z "$GENOME" ]] ; then
        echo 'Must specify a genome fasta file to be annotated'
        echo "Usage: $0 genome.fa [ sampleID ]"
        echo "Version $VERSION"
        exit 1
fi

echo "$0 Version $VERSION"

# A name for the sample
basegenome=$(basename ${GENOME%.fa*})
sampleID=${2:-$basegenome}

# Load UPPMAX modules
# ***************** 
module load bioinfo-tools
module load snap/2013-02-16
module load GeneMark/4.32-es 
module load augustus/3.2.2 &> /dev/null
module load biopython/1.68
# ***************** 

# **** The number of cores used ****
# Set number of cores based on default sbatch command or by user
if [[ $SLURM_CPUS_ON_NODE ]] ; then # If the variable exists
	CORES=$SLURM_CPUS_ON_NODE
elif [[ condition ]]; then
	CORES=1
fi

# --- Variables ---
pathToGeneMark="/pica/sw/apps/bioinfo/GeneMark/4.32-es/milou"
# Scripts
SNAPgff2_gff3="/path2scripts/SNAPgff2_gff3.py"
GMgff2EVM="/path2scripts/GMgff2EVM.py"
subsetfastaIDbio="/path2scripts/subsetfastaIDbio.py"
augustus_GTF_to_EVM_GFF3="/somepath/EVM/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl"

# Training files for Letharia
snapHMM="L.lupinapure.scaffolds_ed.hmm"
localAugustusPATH="somepath2augustusconfigdirectory"
augustusSp="lupina"
GeneMarkMod="Letharia_GeneMark.mod"
# -----------------

# --------------------------------------
echo "### SNAP ###"
# --------------------------------------
snap $snapHMM $GENOME -gff -aa $basegenome.SNAP.faa > $basegenome.SNAP.gff2

echo "Editing the gff2 SNAP file ..."
$SNAPgff2_gff3 $basegenome.SNAP.gff2 -n -p $sampleID > $basegenome.SNAP.gff3

# Clean
rm $basegenome.SNAP.gff2

echo

# --------------------------------------
echo "### GeneMark-ES ###"
# --------------------------------------
# For an incomprehensible reason, GeneMark outputs gff3 file with the
# sequences names only as "seq" or set by the parameter "-s". Therefore,
# multifasta sequences are a disaster because all the features of all the
# contigs get merged. So I had to do the looping to create a final decent gff3
# file. Quite nasty, I know.

# gmhmme3 -m $GeneMarkMod -o $basegenome.GeneMark.gff -p -n -f gff3 $GENOME

# Loop trough all the fasta sequences
for seq in $(grep '>' $GENOME); do 
	# Get the name of each sequence in the fasta file
	nomseq=$(echo $seq | cut -d'>' -f2); 
	
	# Extract the individual sequence from the genome multifasta file
	$subsetfastaIDbio $GENOME $nomseq

	echo "Running gmhmm3 for $nomseq ..."
	gmhmme3 -m $GeneMarkMod -o $nomseq.temp.GeneMark.gff -f gff3 -s $nomseq "$basegenome"_"$nomseq".fa

	# Erase the fasta sequence
	rm "$basegenome"_"$nomseq".fa
done

echo "Merging gmhmm3 gff3 files ..."

printf "##gff-version 3\n# Eukaryotic GeneMark.hmm version 3.51\n" > $basegenome.GeneMark.gff3 

for file in $(ls *temp.GeneMark.gff); do 
	cat $file | grep -v "#" >> $basegenome.GeneMark.gff
done

# Fix it to make it compatible with EVM
$GMgff2EVM $basegenome.GeneMark.gff > $basegenome.GeneMark.gff3

rm *.temp.GeneMark.gff $basegenome.GeneMark.gff

# -p write protein translation
# -n write mucleotide sequence
# -b <output file> output statistics of predicted introns
# -d <file name> provide input for GeneMark.hmm plus
# -s <string> sequence tag in GFF output format
# -f <format> output prediction in [lst|gff3|gtf] format; default [lst]

# Notes I found through trial-error: 
# - the options above can only be used once at a time, it seems.
# - GeneMark doesn't like long names, I found (no more than ~60 characters?).
#   It causes the output gff file to be completely empty, even without header.
# - GeneMark will most likely print only the header if there are tracks of N's
#   that partitions the input sequence.
#   See http://gmod.827538.n3.nabble.com/GeneMark-ES-problem-td3096383.html
#   "Also look for runs of NNNNNNNN in your contigs.  GeneMark auto-splits
#   contigs when there are 40 or more N’s in a row.  The result is that a
#   contig you thought was large is split into many smaller contigs. If none
#   of these smaller contigs are above GeneMark’s default threshold (20,000
#   bp). Then GeneMark will not produce any results.  This is a big gotcha for
#   some assemblies."

echo

# --------------------------------------
echo "### Augustus ###"
# --------------------------------------
# Ab initio, no hints
augustus --species=$augustusSp $GENOME --AUGUSTUS_CONFIG_PATH=$localAugustusPATH/augustus_config > $basegenome.augustus.gtf

# --genemodel=partial, --genemodel=intronless, --genemodel=complete, --genemodel=atleastone or --genemodel=exactlyone
#   partial      : allow prediction of incomplete genes at the sequence boundaries (default)
#   intronless   : only predict single-exon genes like in prokaryotes and some eukaryotes
#   complete     : only predict complete genes
#   atleastone   : predict at least one complete gene
#   exactlyone   : predict exactly one complete gene

# --alternatives-from-evidence=true/false
#   report alternative transcripts when they are suggested by hints

# The posterior probabilities are estimated using a sampling algorithm. The parameter --sample==n adjusts the number 
# of sampling iterations. The higher 'n' is the more accurate is the estimation but it usually isn't important 
# that the posterior probability is very accurate. Every 30 sample iterations take about the same time as one run without 
# sampling, e.g. --sample=60 takes about 3 times as much time as --sample=0 (which was standard up to version 1.6).
# The default is
# --sample=100
# If you do not need the posterior probabilities or alternative transcripts, say
# --sample=0

## *** Get the protein sequences ***
echo "Getting Augustus proteins"
getAnnoFasta.pl $basegenome.augustus.gtf
mv $basegenome.augustus.aa $basegenome.augustus.faa # Rename to match SNAP

## *** Get decent gff ***
echo "Getting decent gff file ..."
$augustus_GTF_to_EVM_GFF3 $basegenome.augustus.gtf > $basegenome.augustus.gff

# Start the format that is more decent
printf "##gff-version 3\n# Augustus 3.2.2\n" > $basegenome.augustus.gff3

# Remove empty lines and ad missing ';' at the end of each line
cat $basegenome.augustus.gff | grep -v -e '^$' | perl -pe 's/(.*)/$1;/' >> $basegenome.augustus.gff3

# Remove the formats that are not useful
rm $basegenome.augustus.gtf $basegenome.augustus.gff
