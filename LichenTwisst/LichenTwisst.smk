### LichenDistances: Trying to guess the origin of the homeologs within the triploid L. lupina metagenome
#############################################################################

# In this pipeline I attempt to phase the triploid L. lupina metagenome using
# the pure-cultre L. lupina genome as reference. The idea is to take the
# alternative allele along the genome and make little SNP windows. Inspired by
# the basic logic of Twisst (Martin et al. 2017 Genetics) for plotting, I make
# little trees from the SNP windows, and classify the found topologies based
# on who is the sister of the lupinas.

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2021/03/08
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

from ete3 import Tree

# -------------------------------------------------
# Data
vcffile = config["vcffile"]
snpsMGdf = config["snpsMGdf"]
gffmat = config["gffmat"] # Annotation file of the scaffold containing the MAT idiomorph

# Parameters
matctg = config["matctg"]

# Environments
renvironment = config["renvironment"]

# Scripts
Trees4TwisstMAT = config["Trees4TwisstMAT"]
Trees4TwisstAll = config["Trees4TwisstAll"]
LethariaTwisstPlotMAT = config["LethariaTwisstPlotMAT"]
LethariaTwisstPlotAll = config["LethariaTwisstPlotAll"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: Trees4Twisst_mat, topos2tableMAT, TwisstPlotMAT, topos2tableAll
# ----------


rule all:
	input:
		"results/LupinaMG_MAF_topologies.pdf",
		"results/dummy.txt"

## --- Contig with the mating type locus ---

rule Trees4Twisst_mat:
	""" Make trees in windows from the vcf file """
	input:
		vcf = vcffile
	output:
		windowdata = "data/Letharia_windows_MAT.txt",
		trees = f"trees/{matctg}.tre"
	params:
		matctg = matctg
	conda: 
		renvironment
	script:
		Trees4TwisstMAT

rule topos2tableMAT:
	""" Append information on the topology type to an input data frame of sequence data windows """
	input:
		windowdata = "data/Letharia_windows_MAT.txt",
		trees = f"trees/{matctg}.tre"
	output:
		windowdata = "data/Letharia_windows_MAT_topos.txt"
	run:
		# Read the window data file
		tabs = [line.rstrip("\n") for line in open(input.windowdata, 'r')] 			# Read tab file into a list and remove the new line
		with open(input.trees, "rt") as tf: alltrees = [Tree(line) for line in tf.readlines()] # Read tab file into a list of trees

		# Define the interesting clades
		MGlupina = ["L.lupina_MG", "L.lupina_culture"]
		MGvulpina = ["L.lupina_MG", "L.vulpina"]
		MGrugosa = ["L.lupina_MG", "L.rugosa"]
		lupinavulpina = ["L.lupina_culture", "L.vulpina"]
		lupinarugosa = ["L.lupina_culture", "L.rugosa"]

		# Open output file
		ofile = open(output.windowdata, 'w')

		# Header
		ofile.write(f"{tabs[0]}\tsisMG\tsislupina\n")

		for i in range(1, len(tabs)):
			t = alltrees[i - 1] # The first line of the windows data file is the header, so go back one for the trees

			# Root with L. columbiana (= L. lucida)
			t.set_outgroup("L.columbiana")

			## What category of tree is it?
			# The sister of the MG lupina
			if t.check_monophyly(MGlupina, target_attr="name")[0]:
				MGsister = "lupina"
			elif t.check_monophyly(MGvulpina, target_attr="name")[0]:
				MGsister = "vulpina"
			elif t.check_monophyly(MGrugosa, target_attr="name")[0]:
				MGsister = "rugosa"
			else:
				MGsister = "other"

			# The sister of the pure culture lupina
			if t.check_monophyly(MGlupina, target_attr="name")[0]:
				lupinasister = "MG"
			elif t.check_monophyly(lupinavulpina, target_attr="name")[0]:
				lupinasister = "vulpina"
			elif t.check_monophyly(lupinarugosa, target_attr="name")[0]:
				lupinasister = "rugosa"
			else:
				lupinasister = "other"

			# Write the new line
			ofile.write(f"{tabs[i]}\t{MGsister}\t{lupinasister}\n")
			
rule TwisstPlotMAT:
	""" Plot a Twisst-ish plot of the MAT contig""" 
	input:
		windowdata = "data/Letharia_windows_MAT_topos.txt",
		snpsMGdf = snpsMGdf,
		gff = gffmat
	output:
		mat = "results/LupinaMG_MAF_topologies.pdf",
	params:
		matctg = matctg
	conda: 
		renvironment
	script:
		LethariaTwisstPlotMAT

## --- Other contigs ---

rule Trees4TwisstAll:
	""" Make trees in windows from the vcf file """
	input:
		vcf = vcffile
	output:
		windowdata = "data/Letharia_windows_all.txt",
	params:
		time = "24:00:00",
		threads = 3,
	conda: 
		renvironment
	script:
		Trees4TwisstAll

rule topos2tableAll:
	""" Append information on the topology type to an input data frame of sequence data windows """
	input:
		windowdata = "data/Letharia_windows_all.txt",
		# trees = f"trees/{matctg}.tre"
	output:
		windowdata = "data/Letharia_windows_all_topos.txt"
	run:
		# Read the window data file
		tabs = [line.rstrip("\n") for line in open(input.windowdata, 'r')] 			# Read tab file into a list and remove the new line

		# Open output file
		ofile = open(output.windowdata, 'w')

		# Header
		ofile.write(f"{tabs[0]}\tsisMG\tsislupina\n")

		# Get a list of all the contigs
		allcontigs = set([tab.split("\t")[0] for tab in tabs[1:]]) # Ignore the header to [1:]

		# Make dictionary of contigs to keep track of how many windows are there per contig
		ctgsdic = dict((key, 0) for key in allcontigs)

		# Define the interesting clades
		MGlupina = ["L.lupina_MG", "L.lupina_culture"]
		MGvulpina = ["L.lupina_MG", "L.vulpina"]
		MGrugosa = ["L.lupina_MG", "L.rugosa"]
		lupinavulpina = ["L.lupina_culture", "L.vulpina"]
		lupinarugosa = ["L.lupina_culture", "L.rugosa"]


		for ctg in allcontigs:
			# Read its corresponding tree file
			with open(f"trees/{ctg}.tre", "rt") as tf: alltrees = [Tree(line) for line in tf.readlines()] # Read tab file into a list of trees

			for tab in tabs:
				if ctg in tab:
					t = alltrees[ctgsdic[ctg]]

					# Root with L. columbiana (= L. lucida)
					t.set_outgroup("L.columbiana")

					## What category of tree is it?
					# The sister of the MG lupina
					if t.check_monophyly(MGlupina, target_attr="name")[0]:
						MGsister = "lupina"
					elif t.check_monophyly(MGvulpina, target_attr="name")[0]:
						MGsister = "vulpina"
					elif t.check_monophyly(MGrugosa, target_attr="name")[0]:
						MGsister = "rugosa"
					else:
						MGsister = "other"

					# The sister of the pure culture lupina
					if t.check_monophyly(MGlupina, target_attr="name")[0]:
						lupinasister = "MG"
					elif t.check_monophyly(lupinavulpina, target_attr="name")[0]:
						lupinasister = "vulpina"
					elif t.check_monophyly(lupinarugosa, target_attr="name")[0]:
						lupinasister = "rugosa"
					else:
						lupinasister = "other"

					# Write the new line
					ofile.write(f"{tab}\t{MGsister}\t{lupinasister}\n")

					ctgsdic[ctg] +=1 # add it at the end to keep base 0

rule TwisstPlotAll:
	""" Plot a Twisst-ish plot of all big contigs""" 
	input:
		windowdata = "data/Letharia_windows_all_topos.txt",
		snpsMGdf = snpsMGdf,
	output:
		dummy = "results/dummy.txt",
	params:
		time = "3:00:00",
		threads = 2,
	conda: 
		renvironment
	script:
		LethariaTwisstPlotAll
