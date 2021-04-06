#!/usr/bin/env Rscript

### Quantify differences between Letharia species and homeologs within L. lupina triploids
#############################################################################
# https://www.molecularecologist.com/2016/02/26/quick-and-dirty-tree-building-in-r/
# https://rdrr.io/cran/phangorn/f/vignettes/Networx.Rmd
# https://www.biostars.org/p/179953/
# https://wiki.duke.edu/display/AnthroTree/2.2+Simple+Parsimony+Analysis+in+R

# Part of the LichenTwisst.smk pipeline
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-06-18
# Version 1
# =======================================
# https://github.com/knausb/vcfR
# browseVignettes(package="vcfR")
# http://t-redactyl.io/blog/2016/02/creating-plots-in-r-using-ggplot2-part-7-histograms.html
# https://stackoverflow.com/questions/4725339/percentage-on-y-lab-in-a-faceted-ggplot-barchart

library(ggplot2)
library(vcfR)
library(tidyr) # For gather
library(plyr) # For revalue and count (it has to be before dplyr to avoid conflicts)
library(dplyr)
library(phangorn)

# ============================
# Data
# ============================
vcf <- read.vcfR(snakemake@input$vcf, verbose = FALSE)

outputfile <- snakemake@output$windowdata

# ============================
# Extensive filtering to chose the SNPs to calculate distance
# ============================
cat("Filtering ...\n")
### ----- Get the fixed part of the vcf file  ----
vcffix <- getFIX(vcf) %>% data.frame %>% select(-c("ID", "FILTER", "QUAL")) # Extract the fixed part of the vcf (and remove some irrelevant columns)

# Rename POS to SITE to avoid confusion later
names(vcffix)[2] <- "SITE"
# Make it numeric
vcffix$SITE <- as.numeric(as.character(vcffix$SITE))

# Modify column of CHROM + POS to match vcfdf
vcffix2 <- vcffix %>% unite(POS, CHROM, SITE, sep = "_", remove = FALSE)

### ----- The variable part of the vcf file ----

# Get the total coverage of the site
coverage <- extract.gt(vcf, element="Cov", as.numeric = TRUE) %>% data.frame() %>% tibble::rownames_to_column("POS") %>% gather("species", "cov", -POS)
# Get depth of the first allele
reads1 <- extract.gt(vcf, element="Reads1", as.numeric = TRUE) %>% data.frame() %>% tibble::rownames_to_column("POS") %>% gather("species", "allele1", -POS)
reads2 <- extract.gt(vcf, element="Reads2", as.numeric = TRUE) %>% data.frame() %>% tibble::rownames_to_column("POS") %>% gather("species", "allele2", -POS)

## Merge dataframe and filter out sites with more than 2 alleles
vcfdf <- cbind(coverage, allele1 = reads1$allele1, allele2 = reads2$allele2) %>% filter((allele1 + allele2) == cov)
# Warning: Notice that done like this some sites will be missing for some
# species, but because I have so much data I think I can afford a few
# differences in the total number of sites

## Calculate the Minor Allele Frequency per site 
vcfdf <- mutate(vcfdf, maf = pmin(allele1, allele2)/cov) 

# Rename the pure culture so it's pretty 
vcfdf$species <- as.factor(vcfdf$species)
levels(vcfdf$species)[levels(vcfdf$species) == "L.lupinapure_postQC"] <- "L.lupina_culture"

### ----- Put the data together ----
cat("Putting data together ...\n")
vcfdf_gt <- merge(vcfdf, vcffix2, by = "POS") # This takes a bit of time

## The allele1 is the same as the reference allele, so let's split the MAF by alelles (REF vs. ALT)
# Remove bad quality variants based on coverage, and sites with indels
vcfdf_gtf <- vcfdf_gt %>% filter(cov < 300 & cov > 50) %>% filter(!grepl('-|,', ALT))

# Fix the levels to remove the filtered-out ones
vcfdf_gtf$ALT <- factor(vcfdf_gtf$ALT)

##---
cat("Genotyping ...\n")
## In the case of the metagenome
# Let's assume that the sites with alternative alleles with maf < 0.1 are sequencing errors. 
# The inverse, reference alleles with maf < 0.1 are likely not mistakes, because then 0.9 reads 
# support a difference from the reference. However, alternative alleles with maf = 0 are correct!

vcfdf_gtf_lupina <- vcfdf_gtf %>% filter(species == "L.lupina") %>% 
  mutate(refcov = allele1/cov, altcov = allele2/cov) %>% 
  filter(altcov == 0 || altcov > 0.1)

## For the haploid samples we filter out the things with low frequency
# Now get the rest of the samples
vcfdf_gtf_notlupi <- vcfdf_gtf %>% filter(species != "L.lupina") %>% 
  mutate(refcov = allele1/cov, altcov = allele2/cov) %>% 
  filter(maf < 0.05)

# Put back together the surviving sites
vcfdf_gtf_clean <- rbind(vcfdf_gtf_lupina, vcfdf_gtf_notlupi)

## Add extra columns and modify the data frame
vcfdf_gtf_alleles <- vcfdf_gtf_clean %>%
  select(-c("allele1", "allele2")) %>%                                                         # Remove the raw counts because I won't need them anymore
  gather("allele", "cov_allele", -c(POS, species, cov, maf, CHROM, SITE, REF, ALT) )           # Make a long format with the allele and its frequency

### Make a matrix with the genotypes
# First the haploids
# Retain only the major allele
genotypeguess <- vcfdf_gtf_alleles %>% filter(species != "L.lupina") %>% 
  filter(cov_allele > 0.95) %>% 
  mutate(genotype = ifelse(allele == "refcov", REF, ALT))

# Now the Metagenome
# The logic here is to break into two genotypes, one like the reference and the alternative
genotypeguess_MG <- vcfdf_gtf_alleles %>% filter(species == "L.lupina", allele == "altcov") %>% 
  mutate(genotype = ifelse(cov_allele < 0.2, REF,ALT)) # This is the key to genotyping. If something looks fixed for reference, keep reference. Otherwise, above 0.2, call the alternative

# Put them together
genotypeguess_all <- rbind(genotypeguess, genotypeguess_MG)
# For some reason it gets the genotype in numbers, not in the base :/
genotypeguess_all$genotype <- revalue(as.factor(genotypeguess_all$genotype), c("1"="A", "2"="C", "3" = "G", "4" = "T"))

# Retain only sites that have 5 occurrences: the 4 haploid samples and the alternative genotype of the MG lupina
goodsites <- genotypeguess_all %>% plyr::count("POS") %>% filter(freq == 5) %>% .$POS # I have to say exactly where "count" came from to avoid conflict with dplyr

# Remove the sites with missing data
genotypeguess_all_nomiss <- genotypeguess_all %>% filter(POS %in% goodsites)

# # How many sites survived? 
# genotypeguess_all_nomiss$POS %>% unique %>% length()
# 949316

# Sort by POS (and by species), it takes a while but it's CRUCIAL to retain homology in the final alignment!!!!
genotypeguess_all_nomiss <- genotypeguess_all_nomiss[order(genotypeguess_all_nomiss$POS, genotypeguess_all_nomiss$species), ]

# Reset the row numbers (optional)
rownames(genotypeguess_all_nomiss) <- NULL

### --- Make alignment ---
cat("Making alignments and trees ...\n")
# Change the name of the genotyping of "L.lupina" to make it clear that it's the (alternative + fixed) genotypes in the metagenome
levels(genotypeguess_all_nomiss$species)[levels(genotypeguess_all_nomiss$species) == "L.lupina"] <- "L.lupina_MG"

# Get the species names
species <- genotypeguess_all_nomiss %>% .$species %>% levels

### --- Read it back ---

# A function to create windows of SNPs from a given contig
winSNPsperwin <- function(genotypes, path2alignments, windowsize = 50){
  species <- genotypes %>% .$species %>% levels
  genotypes <- genotypes[order(genotypes$SITE), ]
  
  # Make results folders
  if (!dir.exists(paste0(path2alignments, "/fastas")) ){
    dir.create(paste0(path2alignments, "/fastas"))
    dir.create(paste0(path2alignments, "/trees"))  
  }

  # Make results folders
  if (!dir.exists(paste0(path2alignments, "/trees")) ){
    dir.create(paste0(path2alignments, "/trees"))  
  }
  
  # Make report to keep track of sites
  report <- data.frame()
  
  for (chromi in unique(genotypes$CHROM)){
    print(chromi)
    
    # Get the sites of this chromosome
    localgenotypes <- genotypes %>% filter(CHROM == chromi)
    
    localsites <- localgenotypes %>% .$SITE %>% unique()
    
    # Define the name of the output
    slicetreename <- paste0(path2alignments, "/trees/", chromi, ".tre")

    # What is the length of the full alignment for this chromosome? (all variable sites in this contig)
    alignlen <- localgenotypes %>% filter(species == species[1]) %>% nrow()

    starts <- seq(1, alignlen, by = windowsize)
    n <- length(starts)    # Find the length of the vector "starts"

    # Loop through each window and make haplotypes of them
    for (i in 1:n){ # Under this design, the last window has less sites
      # print(i)
      startwin <- starts[i]

      if (starts[i] == tail(starts,1)){ # Last window
        endwin <- alignlen
      } else {
        endwin <- starts[i+1] - 1
      }
      
      slicefasta <- paste0(path2alignments, "/fastas/", chromi, "_w", startwin, "-", endwin, ".fa")
      # print(slicefasta)
      
      # What are the sites in this section?
      slicesites <- localsites[startwin:endwin]
      
      ## --- Write the fasta for this window ---
      fileConn <- file(slicefasta) # Open the file
      alignment <- "" # start it empty
      for (sp in species) {
        print(sp)
        thissp <- localgenotypes %>% filter(SITE %in% slicesites, species == sp)
        haplotype <- paste(thissp$genotype, collapse = '') # Make the sequence
        alignment <- c(alignment, paste0(">", sp, "\n", haplotype))
      }
      # See ?ape::read.FASTA() of another way I could have printed the fasta file..
      writeLines(alignment[2:6], fileConn) # Ignore the empty first element
      close(fileConn) # Close the file
      
      ## --- Now make a NJ tree ---
      lichenal <- ape::read.dna(slicefasta, format = "fasta")
      lichen_phyDat <- phyDat(lichenal, type = "DNA", levels = NULL)
      
      # Compute a distance matrix of pairwise differences from DNA
      lichendist <- dist.hamming(lichen_phyDat) # Uncorrected distances
      
      # Compute NJ tree
      slicetree <- NJ(lichendist)
      write.tree(slicetree, file = slicetreename, append = TRUE)
      
      ## --- Window data file in Twisst style, sort of ---
      physicalstart <- slicesites[1]
      physicalend <- slicesites %>% tail(1)
      midwinpoint = (physicalend - physicalstart)/2 + physicalstart
      # put together
      slicereport <- data.frame(scaffold = chromi, start = physicalstart, end = physicalend, mid = midwinpoint, sites = nrow(thissp), startwin = startwin, endwin = endwin)
      print(slicereport)
      report <- rbind(report, slicereport)
    }
  }
  return(report)
}

# Get the base working directory
workingdir <- dirname(outputfile) %>% dirname 

# Get the lenght of the contigs
bigchrs <- genotypeguess_all_nomiss %>% separate(POS, into = c(NA, NA, NA, "len", NA, NA, NA), sep ="_")
genotypeguess_all_nomiss <- cbind(genotypeguess_all_nomiss, len = as.numeric(bigchrs$len)) # Append the new column

# Make the trees for this contig
window_data <- winSNPsperwin(genotypeguess_all_nomiss, workingdir) 

# Output the window data for later
write.table(window_data, file = outputfile, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done!\n")

### ---------------
