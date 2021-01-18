# -*- snakemake -*-

### LichenPloidy: Small pipeline to infer the ploidy of a genome or metagenome
#############################################################################

# This pipeline is designed to calculate the minor allele frequency (MAF)
# distribution of a given sample based on the read counts of biallelic SNPs.
# In this case, the objective is to get MAF distributions of all *Letharia*
# samples using the *L. lupina* pure culture as a reference.

# In addition, I also explore the distribution of alleles along the contigs
# and I try to estimate a rough distance with the other lineages.

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/06/15 - 2020/07/10
# +++++++++++++++++++++++++++++++++++++++++++++++++

from Bio import SeqIO
import gzip # For unzipping files

# -------------------------------------------------
Illumina = config["Illumina"] # path to the folder with all the Illumina data

# Samples
samples = config["samples"]
reference = config["reference"] # full path of the fasta file used as reference
refsampleid = config["refsampleid"]
gffmat = config["gffmat"] # Annotation file of the scaffold containing the MAT idiomorph

# Repeat related stuff
TElibrary = config["TElibrary"] # Repeat library from RepeatModeler
TElibraryid = config["TElibraryid"]

# Scripts
PloidyMAF = config["PloidyMAF"]
HomoeologsDistance = config["HomoeologsDistance"]


# Other parameters
minlen = config["minlen"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: referencegenome, neededscripts, renameconsensi, indexbwa, snpsvcfnomiss, remove_small_scf, plotMAF, HomoeologsDistance, bedtoolsTEs
# ----------

rule all:
	input:
		f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_MAF-noTEs.pdf",
		f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_lupina_unMAF-noTEs.pdf",
		f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_lupina_LOH-noTEs.png",
		f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_lupina_mat-noTEs.png",
		f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_MAF-noTEs_distance.pdf",
		f"RepeatMasker/{refsampleid}.fa.out.gff"

# ------- PREPARE ALL DATA --------
rule referencegenome:
	""" Prepare a copy for the reference genome to avoid putting extra files near the original file """
	input:
		reference
	output:
		f"data/{refsampleid}.fa"
	shell:
		"""
		ln -s {input} {output}
		"""

rule neededscripts:
	""" Get the scripts needed for the rest of the pipeline """
	# About ancient() see https://bitbucket.org/snakemake/snakemake/pull-requests/119/ancient-flag-on-input-files/diff
	output:
		"scripts/renameRMDLconsensi.pl",
		"scripts/totalcovergff.py"
	shell:
		"cd scripts; "
		"wget https://raw.githubusercontent.com/genomicrocosm/physaliaTEcourse/master/Practical2_Computational_annotation/renameRMDLconsensi.pl"

rule renameconsensi:
	input:
		rawlib = TElibrary,
		renamer = "scripts/renameRMDLconsensi.pl"
	output:
		"results/{TElibraryid}_RM.lib"
	shell:
		""" 
		perl {input.renamer} {input.rawlib} {wildcards.TElibraryid} {output}
		"""

# ------- RepeatMasking -------

rule repeatmasker:
	""" Use RepeatMasker to find regions that should be filtered out """
	input:
		assembly = f"data/{refsampleid}.fa",
		TElib = expand("results/{TElibraryid}_RM.lib", TElibraryid = TElibraryid) # For some reason it won't take the f-string
	output:
		f"RepeatMasker/{refsampleid}.fa.out.gff"
	params:
		time = "6:00:00",
		threads = 16,
	shell:
		""" 
		RepeatMasker -pa {params.threads} -a -xsmall -gccalc -gff -excln -lib {input.TElib} -dir RepeatMasker {input.assembly}
		"""

# ---------------------------------

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = f"data/{refsampleid}.fa"
	output:
		index = f"data/{refsampleid}.fa.bwt"
	version: "1"
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	""" Map Illumina reads with BWA """
	input:
		genome = f"data/{refsampleid}.fa",
		index = f"data/{refsampleid}.fa.bwt",
		read1 = Illumina + "/{sample}.1.fq.gz",
		read2 = Illumina + "/{sample}.2.fq.gz",
	output:
		bwaoutput = temp("mapping/{sample}-to-" + f"{refsampleid}.bam.sorted"),
		log = "logs/bwa_mem/{sample}-to-" + f"{refsampleid}.log"
	params:
		time = "3:30:00",
		threads = 10, 
		refbase = refsampleid,
		rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina",
	version: "1"
	shell:
		"""
		(bwa mem {input.genome} {input.read1} {input.read2} -t {params.threads} -R '{params.rg}' -M | samtools view -Su - | samtools sort -l 5 -O bam -T {wildcards.sample}'-to-'{params.refbase} -@ {params.threads} > {output.bwaoutput}) 2> {output.log}
		# -l 5 following Doug
		"""

rule markduplicates:
	""" Mark duplicates in BAM """
	# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
	input:
		bwaoutput = "mapping/{sample}-to-" + f"{refsampleid}.bam.sorted"
	output:
		mdoutput = "mapping/{sample}-to-" + f"{refsampleid}.sorted.debup.bam",
		mdmetrics = "mapping/{sample}-to-" + f"{refsampleid}.sorted.metrics.txt"
	params:
		time = "1:30:00",
		threads = 10,
	version: "1"
	shell:
		"""
		# Using normal Picard
		picard MarkDuplicates I={input.bwaoutput} O={output.mdoutput} M={output.mdmetrics} ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR="temp"

		"""	
		# # VALIDATION_STRINGENCY=ValidationStringency
		# #                               Validation stringency for all SAM files read by this program.  Setting stringency to
		# #                               SILENT can improve performance when processing a BAM file in which variable-length data
		# #                               (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This
		# #                               option can be set to 'null' to clear the default value. Possible values: STRICT,
		# #                               LENIENT, SILENT
		# # CREATE_INDEX=Boolean          Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
		# #                               false. This option can be set to 'null' to clear the default value. Possible values:
		# #                               true, false
		# # TMP_DIR (File)  Default value: null. This option may be specified 0 or more times.

# -------------- VarScan -----------------

rule VarScan:
	""" Call variants using VarScan """
	input:
		ref = f"data/{refsampleid}.fa",
		bams = expand("mapping/{sample}-to-{refsampleid}.sorted.debup.bam", sample = samples, refsampleid = refsampleid)	
	output:
		"variants/Lichens-snps.vcf.gz"
	params:
		time = "10:00:00",
		threads = 1,
	shell:
		"samtools mpileup -f {input.ref} {input.bams} | varscan mpileup2snp --p-value 0.1 --min-var-freq 0.005 | bgzip > {output}"
		# "samtools mpileup -f {input.ref} {input.bams} | varscan mpileup2snp --p-value 0.05 | bgzip > {output}"

	# -q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]

	# --min-coverage	Minimum read depth at a position to make a call [8]
	# --min-reads2	Minimum supporting reads at a position to call variants [2]
	# --min-avg-qual	Minimum base quality at a position to count a read [15]
	# --min-var-freq	Minimum variant allele frequency threshold [0.01]
	# --min-freq-for-hom	Minimum frequency to call homozygote [0.75]
	# --p-value	Default p-value threshold for calling variants [99e-02]
	# --variants	Report only variant (SNP/indel) positions [0]

# -------------- filtering -----------------

rule snpsvcfnomiss:
	""" Filter out sites with missing data """
	input:
		vcf = "variants/Lichens-snps.vcf.gz"
	output:
		filteredvcf = "variants/Lichens-snps-miss1.vcf.gz"
	shell:
		"zless {input.vcf} | grep -v ':-:-:-:-' | bgzip > {output.filteredvcf}; "

rule remove_small_scf:
	""" Remove small scaffolds from VarScan vcf and output a more canonical vcf """
	input:
		ref = f"data/{refsampleid}.fa",
		vcf = "variants/Lichens-snps-miss1.vcf.gz"
	output:
		vcf = f"variants/Lichens-snps-miss1-{int(minlen/1000)}kp.vcf"
	run:	
		bigscf = []
		totalseqs, filtercount = [0, 0]

		# Read fasta file
		for seq_record in SeqIO.parse(input.ref, "fasta"):
			totalseqs += 1
			if len(seq_record) >= minlen:
				bigscf.append(seq_record.id)
				filtercount += 1

		print("Number of sequences in input fasta file: " + str(totalseqs))
		print("Number of sequences larger than " + str(minlen) + " bps: " + str(filtercount))

		vcfopen = gzip.open(input.vcf, 'rt') # 't' is text mode, to interpret the tabs and new lines
		outputvcf = open(output.vcf, "w")

		for line in vcfopen:
			if "Chrom" in line:
				outputvcf.write("##fileformat=VCFv4.X\n")

				# Add explanations
				outputvcf.write(f"##FORMAT=<ID=Cons,Type=String,Description='Consensus in IUPAC confidence'>\n")
				outputvcf.write(f"##FORMAT=<ID=Cov,Number=R,Type=Integer,Description='Coverage'>\n")
				outputvcf.write(f"##FORMAT=<ID=Reads1,Number=R,Type=Integer,Description='Allele 1'>\n")
				outputvcf.write(f"##FORMAT=<ID=Reads2,Number=R,Type=Integer,Description='Allele 2'>\n")
				outputvcf.write(f"##FORMAT=<ID=Freq,Type=String,Description='Frequency of Allele 2'>\n")
				outputvcf.write(f"##FORMAT=<ID=P-value,Number=R,Type=Float,Description='P-value'>\n")

				# Get the names of the samples
				individuals = '\t'.join(sorted(samples))

				# Print a new header
				outputvcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{individuals}\n")

			else:
				line = re.sub(r"[\s]", '\t', line) # Replace white spaces from VarScan with tabs (wtf, why are they there anyway?)
				columns = line.rstrip("\n").split('\t')
				CHROM, POS, REF, ALT, QUAL, FILTER, SamplesRef, SamplesHet, SamplesHom, SamplesNC = columns[0:10]

				FORMAT = "Cons:Cov:Reads1:Reads2:Freq:P-value"
				indivs = '\t'.join(columns[10:])
				indivs = indivs.rstrip("\t")

				if CHROM in bigscf:
					newline = f"{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t.\t{FORMAT}\t{indivs}\n"
					outputvcf.write(newline)

		vcfopen.close()
		outputvcf.close()

# -------------- Removing TEs -----------------

rule bedtoolsTEs:
	""" Filter vcf with the repeats from RepeatMasker """
	input:
		vcf = f"variants/Lichens-snps-miss1-{int(minlen/1000)}kp.vcf",
		gfffile = f"RepeatMasker/{refsampleid}.fa.out.gff"
	output:
		filteredvcf = f"variants/Lichens-snps-miss1-{int(minlen/1000)}kp-noTEs.vcf",
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.filteredvcf}
		
		# Filter out the repeats
		bedtools intersect -a {input.vcf} -b {input.gfffile} -v >> {output.filteredvcf}
		"""

# -------------- Plotting -----------------

rule plotMAF:
	""" Plot the Minor Allele Frequency distributions """
	input:
		vcf = f"variants/Lichens-snps-miss1-{int(minlen/1000)}kp-noTEs.vcf",
		gff = gffmat
	output:
		maf = f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_MAF-noTEs.pdf",
		cov = f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_cov-noTEs.pdf", # coverage
		unmaf = f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_lupina_unMAF-noTEs.pdf", 
		loh = f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_lupina_LOH-noTEs.png", 
		mat = f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_lupina_mat-noTEs.png", 
	conda: 
		"envs/plotmaf.yaml"
	script:
		PloidyMAF	

rule HomoeologsDistance:
	""" Plot the Minor Allele Frequency distributions """
	input:
		vcf = f"variants/Lichens-snps-miss1-{int(minlen/1000)}kp-noTEs.vcf",
	output:
		distance = f"results/Lichens-snps-miss1-{int(minlen/1000)}kp_MAF-noTEs_distance.pdf",
	conda: 
		"envs/plotmaf.yaml"
	script:
		HomoeologsDistance	
