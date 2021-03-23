#!/usr/bin/env Rscript

### Inferring ploidy based on the Minor Allele Frequency distribution
#############################################################################

# By plotting the Minor allele frequency distributions (MAF) within the
# Illumina reads of a given sample we can infer the ploidy level (see
# Ament-Velasquez et al. 2016 Mol Ecol). 

# Part of the LichenPloidy.smk pipeline
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-06-18 - 2021-03-23
# Version 2
# =======================================
# https://github.com/knausb/vcfR
# browseVignettes(package="vcfR")
# http://t-redactyl.io/blog/2016/02/creating-plots-in-r-using-ggplot2-part-7-histograms.html
# https://stackoverflow.com/questions/4725339/percentage-on-y-lab-in-a-faceted-ggplot-barchart


library(ggplot2)
library(dplyr)
library(vcfR)
library(tidyr) # For gather

# ============================
# Data
# ============================
cat("Reading data ...\n")
vcf <- read.vcfR(snakemake@input$vcf, verbose = FALSE)
gff_file <- snakemake@input$gff

### ----- 

# Get the total coverage of the site
coverage <- extract.gt(vcf, element="Cov", as.numeric = TRUE) %>% data.frame() %>% tibble::rownames_to_column("POS") %>% gather("species", "cov", -POS)
# Get depth of the first allele
reads1 <- extract.gt(vcf, element="Reads1", as.numeric = TRUE) %>% data.frame() %>% tibble::rownames_to_column("POS") %>% gather("species", "allele1", -POS)
reads2 <- extract.gt(vcf, element="Reads2", as.numeric = TRUE) %>% data.frame() %>% tibble::rownames_to_column("POS") %>% gather("species", "allele2", -POS)

## Merge database and filter out sites with more than 2 alleles
vcfdf <- cbind(coverage, allele1 = reads1$allele1, allele2 = reads2$allele2) %>% filter((allele1 + allele2) == cov)
# Warning: Notice that done like this some sites will be missing for some
# species, but because I have so much data I think I can afford a few
# differences in the total number of sites

## Calculate the Minor Allele Frequency per site and remove homozygous sites
vcfdf <- mutate(vcfdf, maf = pmin(allele1, allele2)/cov) %>% filter(maf > 0)

# Rename the pure culture so it's pretty (LETHARIA SPECIFIC)
vcfdf$species <- as.factor(vcfdf$species)
levels(vcfdf$species)[levels(vcfdf$species) == "L.lupinapure_postQC"] <- "L.lupina_culture"

# ============================
# Plot MAF while excluding sites with bad coverage 
# ============================
cat("Plotting MAF distribution ...\n")
# The filter for maf >= 0.01 is done because the pure culture has more coverage and thus can get smaller MAF than the metagenomes.
# The tricky part in geom_histogram is there to rescale the histogram counts so that the bar areas sum 1, but per panel
# https://stackoverflow.com/questions/4725339/percentage-on-y-lab-in-a-faceted-ggplot-barchart
mafplot <- ggplot(vcfdf %>% filter(cov < 300 & cov > 50) %>% filter(maf >=0.01), aes(x = maf, fill = species)) + 
  # geom_histogram(fill = "darkgoldenrod3") + ylab("Count") +  ## Plot the raw counts
  geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
  facet_grid(. ~ species) + xlab("Minor allele frequency") + ylab("Frequency of sites") +
  theme_light() +
  coord_cartesian(ylim=c(0, 0.25)) + # Remove the loooong values to improve clarity 
  scale_fill_manual(values = c("L.columbiana" = "#ecbdecff", 
                               "L.lupina" = "#c7e1a8ff", 
                               "L.lupina_culture" = "#c7e1a8ff",
                               "L.rugosa" = "#f9d5bdff", 
                               "L.vulpina" = "#b7d2fdff")) + 
  theme(strip.text.x = element_text(size = 13),
        axis.title=element_text(size=13)) + 
  guides(fill=FALSE) # Remove legend

ggsave(plot = mafplot, snakemake@output$maf, width = 10, height = 2.5)

# ============================
## Plot the coverage (remove the super extreme ones)
# ============================
cat("Plotting coverage distribution ...\n")
# vignette('sequence_coverage')
covplot <- ggplot(vcfdf %>% filter(cov < 500), aes(x=species, y=cov)) + 
  geom_violin(fill="#C0C0C0", adjust=1.0, scale = "count", trim=TRUE) +
  theme_bw() + ylab("Depth of coverage") +
  theme(axis.title.x = element_blank()) #+ 
# stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") # It requites package Hmisc and for some reason I can't install it

ggsave(plot = covplot, snakemake@output$cov, width = 6.5, height = 2.5)

# NOTES:
# I tried removing SNPs that are indels but that didn't make any effect, so it's not worth it.

# ============================
## Plot the allele frequencies distribution unfolded
# ============================
cat("Plotting unfolded MAF distribution ...\n")
#### ---- Get the fixed part of the vcf file
vcffix <- getFIX(vcf) %>% data.frame # Extract the fixed part of the vcf

## Massage the dataframe a bit more
names(vcffix)[2] <- "SITE" # Change the name of column so I can preserve it later on
# Modify colum of CHROM + POS to match vcfdf (and remove some irrelevant columns)
vcffix2 <- vcffix %>% unite(POS, CHROM, SITE, sep = "_", remove = FALSE) %>% select(-c("ID", "FILTER", "QUAL"))
vcffix2$SITE <- as.numeric(as.character(vcffix2$SITE)) # Make SITE numeric 

# Give the fixed info to the lupina metagenome
vcfdf_lupinaMG <- merge(vcfdf %>% filter(species == "L.lupina"), vcffix2, by = "POS")

## The allele1 is the same as the reference allele, so let's split the MAF by alelles (REF vs. ALT)
# Remove bad quality variants based on coverage
vcfdf_lupinaMG <- vcfdf_lupinaMG %>% filter(cov >= 50 & cov < 150) %>% filter(!grepl('-|,', ALT))

vcfdf_lupinaMG_alleles <- vcfdf_lupinaMG %>% mutate(REF = allele1/cov, ALT = allele2/cov) %>% # Make new columns with the frequencies of each allele
  select(-c("allele1", "allele2")) %>%                                                        # Remove the raw counts because I won't need them anymore
  gather("allele", "cov_allele", -c(POS, species, cov, maf, CHROM, SITE) )                    # Make a long format with the allele and its frequency

# Plot the unfolded distribution
unfoldedMAF <- ggplot(vcfdf_lupinaMG_alleles, aes(x = cov_allele, fill = allele)) + 
  # geom_histogram(alpha = 0.5, position="identity") + 
  geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, position="identity") + # Re-scale to sum 1 across all sites
  xlab("Minor allele frequency") + ylab("Frequency of sites") +
  theme_light() +
  # facet_grid(. ~ allele) +
  scale_fill_manual(values = c("darkgoldenrod3", "chartreuse4"), labels = c("other", "lupina")) # Lupina is green in the phylogenies

ggsave(plot = unfoldedMAF, snakemake@output$unmaf, width = 5, height = 3)

# ============================
## Plot the allele frequencies along contigs in windows
# ============================
cat("Plotting allele frequencies along contigs ...\n")
# Recover the length of the scaffold
# vcfdf_lupinaMG_alleles %>% head %>% separate(POS, into = c("NODE", "n", "length", "actual_len", "cov2", "x", "SITE"), sep ="_")
vcfdf_lupinaMG_df <- vcfdf_lupinaMG_alleles %>% separate(POS, into = c(NA, NA, NA, "len", NA, NA, "SITE"), sep ="_") # This makes the SITE column as character rather than factors
# Make SITE numeric
vcfdf_lupinaMG_df$SITE <- as.numeric(vcfdf_lupinaMG_df$SITE)

# Write this data frame to be used later
LlupinaMG_SNPs <- snakemake@output$SNPs
write.table(vcfdf_lupinaMG_df, file = LlupinaMG_SNPs, sep = "\t", quote = FALSE, row.names = FALSE)

### ---------------
# A function to create windows of MAF per windows using a non-overlapping window size
winmafperwin <- function(dfmaf, chrlines, windowsize = 50000, minsnps = 0){
  ## Notice that chrlines is a data frame with two columns, chr and positions, like such:
  #             chr positions
  # 1 chromosome_1   8813524
  # 2 chromosome_2   5165621
  
  # Make empty data frame
  winmaf <- data.frame(CHROM = c(), species = c(), cov_allele = c())
  # Do this per sample
  for (sample in dfmaf$species %>% unique()){
    print(sample)
    samplemaf <- dfmaf %>% filter(species == sample) # Get the data for only that sample
    
    # Calculate the window cov_allele per chromosome
    for (chromi in chrlines$CHROM){
      chrsamplemaf <- samplemaf %>% filter(CHROM == chromi)
      
      # What's the length of the chromosome?
      chrlen <- chrlines %>% filter(CHROM == chromi) %>% .$len
      
      starts <- seq(1, chrlen-windowsize, by = windowsize)
      n <- length(starts)    # Find the length of the vector "starts"
      
      # Loop through each window and get the MAF of SNPs in it
      for (i in 2:n-1){
        currentwin <- chrsamplemaf %>% filter(SITE >= starts[i] & SITE < starts[i+1])
        chunk <- currentwin$cov_allele %>% median()
        nsnps_win <- currentwin$cov_allele %>% length()
        
        if (nsnps_win >= minsnps) {
          winmaf <- rbind(winmaf, data.frame(CHROM = chromi, species = sample, SITE = starts[i], cov_allele = chunk))
        } else {
          winmaf <- rbind(winmaf, data.frame(CHROM = chromi, species = sample, SITE = starts[i], cov_allele = NA))
        }
        
      }
    }
  }
  return(winmaf)
}
### ---------------

# Filter the data to have only the big scaffods, the alternative allele, and remove obvious sequencing errors 
vcfdf_lupinaMG_df_ALT <- vcfdf_lupinaMG_df %>% 
  filter(allele == "ALT") %>%
  filter(as.numeric(len) >= 1000000) %>%
  filter(maf >= 0.01) %>% 
  sample_n(size = 50000) # Use slice_sample(n = 50000) for higher versions of dyplir
   
# Make a data frame with the contigs and their lengths
relevantchr <- data.frame(CHROM = vcfdf_lupinaMG_df %>% filter(as.numeric(len) >= 1000000) %>% .$CHROM %>% unique()) # Just the big contigs
chrlines <- data.frame(CHROM = relevantchr$CHROM, len = relevantchr %>% separate(CHROM, into = c(NA, NA, NA, "len", NA, NA), sep ="_") %>% .$len %>% as.numeric)

## Calculate the windows
winmaf_lupina <- winmafperwin(vcfdf_lupinaMG_df %>% filter(allele == "ALT"), windowsize = 7500, chrlines = chrlines)

## Make a data frame with arrows so I can highlight some particular areas
arrowsloh <- data.frame(SITE = c(950000, 1550000, 1327000, 90000, 870000),
                        cov_allele = c(0.6, 0.6, 0.3, 0.6, 0.6),
                        y2 = c(0.3, 0.3, 0.6, 0.3, 0.3),
                        species = "L.lupina",
                        CHROM = c("NODE_1_length_1806572_cov_94.6139", "NODE_1_length_1806572_cov_94.6139", "NODE_2_length_1652252_cov_91.7152", "NODE_5_length_1381321_cov_89.3474", "NODE_7_length_1185015_cov_91.9694") )

## Finally plot the alternative allele frequency along the pure culture lupina reference
vcfdf_lupinaMG_LOH <- ggplot(winmaf_lupina, aes(x = SITE, y = cov_allele)) + 
  facet_grid(CHROM ~ ., scales = "free") +
  xlab(expression(paste("Scaffold position (bp) in the ", italic("L. lupina"), " pure culture assembly"))) + 
  ylab(expression(paste('"Other" allele frequency in the ', italic("L. lupina"), " metagenome reads"))) +
  scale_y_continuous(limits = c(0, 1), breaks=c(0, 0.33, 0.66, 1)) + # Change the ticks to match the expectations
  geom_point(data = vcfdf_lupinaMG_df_ALT, aes(x = as.numeric(SITE), y = cov_allele), alpha = 0.1, size = 0.5, colour = "darkgoldenrod2") + # Plot the raw data too
  geom_line(size = 0.7, alpha = 0.8, colour = "darkgoldenrod3") +
  geom_segment( aes(x = SITE, y = cov_allele, xend = SITE, yend = y2), data = arrowsloh, arrow = arrow(length = unit(0.08, "npc")) ) + # Put arrows
  theme_bw() +
  theme(strip.text.y = element_text(size = 4)) # Change the size of the scaffolds' name

ggsave(plot = vcfdf_lupinaMG_LOH, snakemake@output$loh, width = 10, height = 10)

# ============================
## Plot the allele frequencies along the mating type scaffold in windows
# ============================
cat("Plotting allele frequencies along the MAT contig ...\n")
# Genes in the idiomorph
gff2genenames <- function(gff){
  # version 2
  genenames <- c()
  for(i in 1:nrow(gff)){
    attrib <- gff[i,9] %>% as.character() %>% strsplit(.,";") %>% .[[1]]
    name <- attrib[pmatch("Name=", attrib)] %>% strsplit(.,"=") %>% .[[1]] %>% .[2]
    genenames <- c(genenames, name)
  }
  gff <- cbind(gff, genenames) %>% select(V1, V3, V4, V5, genenames)
  return(gff)
}

gff <- read.table(gff_file) %>% filter(V3 == "gene") %>% data.frame %>% gff2genenames() %>% filter(!genenames %in% c("APN2", "SLA2"))

# Make a data frame with the gff
bottomy <- -0.003 # A lower position for the genes
genes <- data.frame(allele = "REF", x1 = gff$V4, y1 = bottomy, x2 = gff$V5, y2 = bottomy)

# Make a data frame with the contigs and their lengths
relevantchr_mat <- data.frame(CHROM = "NODE_87_length_133277_cov_84.8955")
chrlines_mat <- data.frame(CHROM = relevantchr_mat$CHROM, len = relevantchr_mat %>% separate(CHROM, into = c(NA, NA, NA, "len", NA, NA), sep ="_") %>% .$len %>% as.numeric)

## Calculate the windows
vcfdf_lupinaMG_df_mat <- vcfdf_lupinaMG_df %>% 
  filter(CHROM == "NODE_87_length_133277_cov_84.8955") %>%
  filter(maf >= 0.01) # Remove obvious sequencing errors

## I want to have both alleles in the same plot
winmaf_lupina_mat_alt <- winmafperwin(vcfdf_lupinaMG_df_mat %>% filter(allele == "ALT"), 
                                      windowsize = 2500, chrlines = chrlines_mat)

winmaf_lupina_mat_ref <- winmafperwin(vcfdf_lupinaMG_df_mat %>% filter(allele == "REF"), 
                                      windowsize = 2500, chrlines = chrlines_mat)
# Put them in the same data.frame
winmaf_lupina_mat <- rbind(cbind(winmaf_lupina_mat_alt, allele = "ALT"), cbind(winmaf_lupina_mat_ref, allele = "REF"))

vcfdf_lupinaMG_mat <- ggplot(winmaf_lupina_mat, aes(x = SITE, y = cov_allele, colour = allele)) +
  geom_line(size = 0.7, alpha = 0.8) +
  theme_bw() + xlab("Scaffold position") + ylab("Allele frequency in reads") +
  geom_segment(data = genes, aes(x = x1, y = y1, xend = x2, yend = y2), size = 1.5, colour = "black") + # draw the genes
  annotate("text", x = genes$x1[1] + 2000 + (genes[length(genes$x1),2] - genes$x1[1])/2, y = bottomy - 0.04, label = "MAT idiomorph") +
  geom_point(data = vcfdf_lupinaMG_df_mat, aes(x = SITE, y = cov_allele), alpha = 0.5, size = 0.5) +
  scale_y_continuous(breaks=c(0, 0.33, 0.66, 1)) + # Change the ticks to match the expectations
  scale_colour_manual(values = c("darkgoldenrod3", "chartreuse4"), labels = c("other", "lupina")) # Lupina is green in the phylogenies

ggsave(plot = vcfdf_lupinaMG_mat, snakemake@output$mat, width = 10, height = 3)

cat("Done!\n")