#!/bin/bash

#SBATCH -A XXXXXXXX
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH --mail-user something@email.com
#SBATCH --mail-type=ALL

# ==== MAPPING RNAseq Illumina reads to reference genome ===
# Script that maps RNAseq paired-end reads to a given reference genome using
# STAR and SAMtools. It uses Doug's script mergepileupcolumns to produce a
# file of coverage per site. The output files include a sorted BAM file.

# $ wget https://raw.githubusercontent.com/douglasgscofield/bioinfo/master/scripts/mergePileupColumns
# $ chmod a+x mergePileupColumns

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2016
# +++++++++++++++++++++++++++++++++++++++++++++++++

star_path="STAR"

# Fancy option
#---------------------------
while getopts ":f:a:b:s:h" varname; do
	case $varname in
	    f)
	    	# echo "-f was triggered, Parameter: $OPTARG"
			if [[ ! -e "$OPTARG" || ! -s "$OPTARG" ]] ; 			# -e file exists, -s file is not zero size
				then 
					echo "ERROR: Could not find reference genome fasta file $OPTARG or file is size zero" >&2
					exit 1
				else
					genome=$OPTARG 								# To make things clearer
			fi 
	    	;;
	    a)
			if [[ ! -e "$OPTARG" || ! -s "$OPTARG" ]] ; 			# -e file exists, -s file is not zero size
				then 
					echo "ERROR: Could not find read 1 file $OPTARG or file is size zero" >&2
					exit 1
				else
					reads1=$OPTARG 								# To make things clearer
			fi 
	    	;;
	    b)
			if [[ ! -e "$OPTARG" || ! -s "$OPTARG" ]] ; 			# -e file exists, -s file is not zero size
				then 
					echo "ERROR: Could not find read 2 file $OPTARG or file is size zero" >&2
					exit 1
				else
					reads2=$OPTARG 								# To make things clearer
			fi 
	    	;;
	    s)
	    	genomename=$OPTARG 										# Assign a different name to the output files
	    	;;
	    h)
	    	echo "Usage: bash `basename $0` -f path/to/reference_genome_file.fa -a path/to/reads1.fastq.gz -b path/to/reads2.fastq.gz [ -s output_string ] [ -h for help ]"
			exit 1
	    	;;
	    \?)					 										# If it didn't recognize the flag, it assigned to '?'
	    	echo "ERROR: Invalid option -$OPTARG" >&2 					# >&2 -- Print into stderr
	    	exit 1
	    	;;
	    :)
	    	echo "ERROR: Option -$OPTARG requires an argument." >&2
	    	exit 1
	    	;;
	esac
done
#---------------------------

VERSION="2.1"

# ---------------------------
# Check input
# ---------------------------
if [[ -z "$genome" || -z "$reads1" || -z "$reads2" ]] ; then 								# -z string, True if the length of string is zero.
	echo "Usage: bash `basename $0` -f path/to/reference_genome_file.fa -a path/to/reads1.fastq.gz -b path/to/reads2.fastq.gz [ -s output_string ] [ -h for help ]"
	echo "$0 Version $VERSION"
	exit 1
fi

# Assign a name to output
if [[ -z "$genomename" ]]; then
	genomename=$(basename $genome | cut -d'.' -f1)
fi
# ---------------------------

# *****************
module load bioinfo-tools
module load star/2.5.1b
module load samtools/1.3
module load cufflinks/2.2.1
# *****************

echo "$0 Version $VERSION"

# **** The number of cores used ****
# Set number of cores based on default sbatch command or by user
if [[ $SLURM_CPUS_ON_NODE ]] ; then # If the variable exists
	CORES=$SLURM_CPUS_ON_NODE
elif [[ condition ]]; then
	CORES=1
fi
echo "Using $CORES core(s) and $SLURM_NNODES node(s)"

RAM=$(( $CORES * 8000000000 * $SLURM_NNODES)) # Calculate RAM limit based on No. of cores

echo "Using $RAM bytes of RAM"

set -e # The script will terminate after the first line that fails (returns nonzero exit code)
date
# ---------------------------
echo 'Generating STAR index from genome fasta sequence...'
# ---------------------------
# Make a new directory to save the indexing
mkdir -p $genomename'_'GenomeDir

star --runMode genomeGenerate --genomeDir ./$genomename'_'GenomeDir --genomeFastaFiles $genome --runThreadN $CORES --limitGenomeGenerateRAM $RAM --genomeLoad NoSharedMemory --genomeSAindexNbases 3 # Needed for L. lucida's scaffold
# genomeLoad=NoSharedMemory, shared memory is not used. This option is recommended if the shared memory is not configured properly on your server.
# --genomeSAIndexNbases 4 or 5 for small genomes (formula log2(numBases)/2 -1)

echo
# ---------------------------
echo 'Mapping...'
# Unstranded RNA-seq data
star --genomeDir ./$genomename'_'GenomeDir --readFilesIn $reads1 $reads2 --runThreadN $CORES --alignIntronMax 1000 --readFilesCommand zcat --outFileNamePrefix $genomename --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS 
# star --genomeDir ./$genomename'_'GenomeDir --readFilesIn $reads1 $reads2 --runThreadN $CORES --alignIntronMax 1000 --readFilesCommand zcat --outFileNamePrefix $genomename --outSAMstrandField intronMotif # --outSAMattributes ALL 
# outSAMattributes To make it compatible with SAMtools
# outSAMstrandField intronMotif to make it compatible with Cufflinks/Cuffdiff

# For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments
# with XS strand attribute, which STAR will generate with --outSAMstrandField
# intronMotif option. As required, the XS strand attribute will be generated for
# all alignments that contain splice junctions. The spliced alignments that have
# undefined strand (i.e. containing only non-canonical unannotated junctions)
# will be suppressed.

# If you have stranded RNA-seq data, you do not need to use any specific STAR
# options. Instead, you need to run Cufflinks with the library option --library-
# type options. For example, cufflinks ... --library-type fr-firststrand should
# be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq.
# This option has to be used only for Cufflinks runs and not for STAR runs.

# ---------------------------

echo
# ---------------------------
echo 'Indexing the reference for SAMtools...'
#---------------------------
samtools faidx $genome
echo "Indexed!"

# echo
# ---------------------------
echo 'Transforming SAM files to BAM file...'				# Very Slow
# ---------------------------
# samtools view -Shu $genomename'_Aligned.out.sam' | samtools sort - $genomename'_Aligned.out'.sorted   # Using samtools view -u gives uncompressed bam, and these do not have the EOF marker. --> bug of SAMtools if reports an error
samtools view -Shu $genomename'_Aligned.out.sam' | samtools sort -l 5 -O bam -T $genomename'_Aligned.out'.sorted -@ $CORES > $genomename'_Aligned.out'.sorted.bam # Using samtools view -u gives uncompressed bam, and these do not have the EOF marker. --> bug of SAMtools if reports an error

# Get some statistics of the result
samtools flagstat $genomename'_Aligned.out'.sorted.bam > $genomename'_'flagstat.txt

# ---------------------------
echo 'Indexing the sorted BAM...'
#---------------------------
samtools index $genomename'_Aligned.out'.sorted.bam
# Useful for visualizing with IGV

echo
# ---------------------------
echo 'Getting coverage from BAM file...'
# ---------------------------
# Following https://github.com/douglasgscofield/bioinfo/tree/master/scripts#mergepileupcolumns
# and https://github.com/douglasgscofield/bioinfo/tree/master/scripts#windowwig
samtools mpileup -AB -d1000000 -q0 -Q0 -f $genome $genomename'_Aligned.out'.sorted.bam | /proj/b2015200/tools/mergePileupColumns | cut -f1,2,4  > $genomename'_coverage.txt'  #  produce three columns: scaffold, position, coverage

echo
# ---------------------------
echo 'Running Cufflinks...'
# ---------------------------
# # Nonstranded RNAseq
# cufflinks $genomename'_Aligned.out'.sorted.bam --num-threads $CORES 2>&1  # Redirect the stderr to stdout

# Stranded RNAseq 
cufflinks $genomename'_Aligned.out'.sorted.bam --library-type fr-firststrand --num-threads $CORES 2>&1  # Redirect the stderr to stdout

# fr-unstranded: Standard Illumina Reads from the left-most end of the
# fragment (in transcript coordinates) map to the transcript strand, and the
# right-most end maps to the opposite strand.

# fr-firststrand:dUTP, NSR, NNSR Same as above except we enforce the rule that
# the right-most end of the fragment (in transcript coordinates) is the first
# sequenced (or only sequenced for single-end reads). Equivalently, it is
# assumed that only the strand generated during first strand synthesis is
# sequenced.

# fr-secondstrand: Ligation, Standard SOLiD Same as above except we enforce
# the rule that the left-most end of the fragment (in transcript coordinates)
# is the first sequenced (or only sequenced for single-end reads).
# Equivalently, it is assumed that only the strand generated during second
# strand synthesis is sequenced.

echo
# ---------------------------
echo 'Transforming Cufflinks GTF to GFF3...'
# ---------------------------
gffread -E transcripts.gtf -o- > $genomename'_'RNAcuff.gff # add | tail -n +3 to take out the headers

echo 
# ---------------------------
echo 'Making Cufflinks GFF3 a valid GFF3...'
# ---------------------------
/proj/b2013277/scripts/gffread2EVM.py $genomename'_'RNAcuff.gff -n > $genomename'_'RNAcuff.gff3 && rm $genomename'_'RNAcuff.gff

echo "DONE!!"
date
