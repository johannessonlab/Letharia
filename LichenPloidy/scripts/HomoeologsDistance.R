#!/usr/bin/env Rscript

### Quantify differences between Letharia species and homeologs within L. lupina triploids
#############################################################################

# Part of the LichenPloidy.smk pipeline
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

# ============================
# Data
# ============================
vcf <- read.vcfR(snakemake@input$vcf, verbose = FALSE)

# ============================
# Extensive filtering to chose the SNPs to calculate distance
# ============================
cat("Filtering ...\n")
### ----- Get the fixed part of the vcf file  ----
vcffix <- getFIX(vcf) %>% data.frame %>% select(-c("ID", "FILTER", "QUAL")) # Extract the fixed part of the vcf (and remove some irrelevant columns)

# Rename POS to SITE to avoid confusion later
names(vcffix)[2] <- "SITE"
# MAke it numeric
vcffix$SITE <- as.numeric(as.character(vcffix$SITE))

# Modify colum of CHROM + POS to match vcfdf
vcffix2 <- vcffix %>% unite(POS, CHROM, SITE, sep = "_", remove = FALSE)

### ----- The variable part of the vcf file ----

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

## Calculate the Minor Allele Frequency per site 
vcfdf <- mutate(vcfdf, maf = pmin(allele1, allele2)/cov) 

# Rename the pure culture so it's pretty 
vcfdf$species <- as.factor(vcfdf$species)
levels(vcfdf$species)[levels(vcfdf$species) == "L.lupinapure_postQC"] <- "L.lupina_culture"

### ----- Put the data together ----
vcfdf_gt <- merge(vcfdf, vcffix2, by = "POS") # This takes a bit of time

## The allele1 is the same as the reference allele, so let's split the MAF by alelles (REF vs. ALT)
# Remove bad quality variants based on coverage, and sites with indels
vcfdf_gtf <- vcfdf_gt %>% filter(cov < 300 & cov > 50) %>% filter(!grepl('-|,', ALT))
# Fix the levels to remove the filtered-out ones
vcfdf_gtf$ALT <- factor(vcfdf_gtf$ALT, levels = vcfdf_gtf$ALT %>% unique %>% rev %>% droplevels)

## Add extra columns and modify the data frame
vcfdf_gtf_alleles <- vcfdf_gtf %>% mutate(lupina = allele1/cov, alternative = allele2/cov) %>% # Make new columns with the frequencies of each allele
  select(-c("allele1", "allele2")) %>%                                                         # Remove the raw counts because I won't need them anymore
  gather("allele", "cov_allele", -c(POS, species, cov, maf, CHROM, SITE, REF, ALT) )           # Make a long format with the allele and its frequency


# ============================
# Distance at sites that are polymorphic (heterozygous) in the lupina metagenome
# ============================
cat("Extracting heterozygous sites ...\n")
# Remove obvious sequencing mistakes in the lupina MG
vcfdf_gtf_alleles_lupina <- vcfdf_gtf_alleles %>% filter(species == "L.lupina", maf >= 0.20) 

# Now get those sites from the other species
vcfdf_gtf_alleles_spp <- vcfdf_gtf_alleles %>% filter(POS %in% vcfdf_gtf_alleles_lupina$POS)

## We think the rest are haploid, so filter out the things with low frequency
vcfdf_gtf_alleles_notlupi <- vcfdf_gtf_alleles_spp %>% 
  filter(species != "L.lupina") %>% 
  filter(cov_allele >= 0.9)

# Put them back together
vcfspp <- rbind(vcfdf_gtf_alleles_notlupi, vcfdf_gtf_alleles_lupina)

### Make a matrix with the genotypes
# First the haploids
prematrix <- vcfspp %>% filter(species != "L.lupina") %>% mutate(genotype = ifelse(allele == "lupina", REF, ALT))
prematrix$genotype <- revalue(as.factor(prematrix$genotype), c("1"="A", "2"="C", "3" = "G", "4" = "T")) # For some reason it gets the genotype in numbers, not in the base :/

## Now the triploid lupina, but only the allele that is different from the pure culture:
# prematrix <- rbind(prematrix, vcfspp %>% filter(species == "L.lupina", allele != "lupina") %>% mutate(genotype = ifelse(allele == "lupina", REF, ALT)) ) # If I do it like this it stays numeric
prematrix <- rbind(prematrix, vcfspp %>% filter(species == "L.lupina", allele != "lupina") %>% mutate(genotype = ALT))

## How many sites did I get per species?
# prematrix %>% select(species, genotype) %>% count("species")

# ============================
# Distance at sites that are fixed for the alternative allele
# ============================
cat("Extracting sites fixed for the alternative allele in the lupina MG ...\n")
# Get sites in the lupina MG that are fixed for the alternative allele
vcfdf_gtf_alleles_altfix <- vcfdf_gtf_alleles %>% filter(species == "L.lupina", allele == "alternative", cov_allele >= 0.9)

# Now get those sites from the other species
vcfdf_gtf_alleles_spp_altfix <- vcfdf_gtf_alleles %>% filter(POS %in% vcfdf_gtf_alleles_altfix$POS) %>% filter(cov_allele >= 0.9)

### Make a matrix with the genotypes
prematrix_altfix <- vcfdf_gtf_alleles_spp_altfix %>% mutate(genotype = ifelse(allele == "lupina", REF, ALT))
prematrix_altfix$genotype <- revalue(as.factor(prematrix_altfix$genotype), c("1"="A", "2"="C", "3" = "G", "4" = "T")) # For some reason it gets the genotype in numbers, not in the base :/

# ============================
# Distance at sites that are fixed for the lupina pure culture allele
# ============================
cat("Extracting sites fixed for the lupina culture allele in the lupina MG ...\n")
# Get sites in the lupina MG that are fixed for the pure culture allele
vcfdf_gtf_alleles_reffix <- vcfdf_gtf_alleles %>% filter(species == "L.lupina", allele == "lupina", cov_allele >= 0.9)

# Now get those sites from the other species and filter out sites that are probably sequencing errors
vcfdf_gtf_alleles_spp_reffix <- vcfdf_gtf_alleles %>% filter(POS %in% vcfdf_gtf_alleles_reffix$POS) %>% filter(cov_allele >= 0.9)

### Make a matrix with the genotypes
prematrix_reffix <- vcfdf_gtf_alleles_spp_reffix %>% mutate(genotype = ifelse(allele == "lupina", REF, ALT))
prematrix_reffix$genotype <- revalue(as.factor(prematrix_reffix$genotype), c("1"="A", "2"="C", "3" = "G", "4" = "T")) # For some reason it gets the genotype in numbers, not in the base :/

# In this case only the L.lupina_ref column is meaningfull, L.lupina_alt is not really the alternative because it's the area fixed for the pure culture so they are the same thing

# ============================
# Put them together
# ============================
cat("Calculating distances ...\n")
## Calculate the proportion of sites that are not shared between species
## Each species input is a vector of sites, and they correpond to each other (species1[1] is homologous to species2[1])
snpdistance <- function(species1, species2){
  varsites <- (species1 != species2) %>% as.numeric() %>% sum
  invarsites <- (species1 == species2) %>% as.numeric() %>% sum
  
  distance <- varsites/(varsites + invarsites)
  return(distance)
}

## From the prematrix produce a data frame with distances against the other species
prematrix2distance <- function(prematrix){
  # Remove all sites that have missing data (there are 5 taxa)
  goodsites <- prematrix %>% plyr::count("POS") %>% filter(freq == 5) %>% .$POS # I have to say exactly where "count" came from to avoid conflict with dplyr
  prematrix <- prematrix %>% filter(POS %in% goodsites)
  
  ## --- Transform from long to wide format and extract the alternative allele (different from the pure culture) from the lupina MG
  ## http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
  # The arguments to spread():
  # - data: Data object
  # - key: Name of column containing the new column names
  # - value: Name of column containing values
  matrixl <- spread(prematrix %>% select(POS, species, genotype), species, genotype)
  # To avoid confusion, rename the lupina MG column
  names(matrixl)[3] <- "L.lupina_alt" 
  
  ## Put together a data frame of distances in a wide format
  # cat("Calculating distances ...")
  lethariadiv <- data.frame(Species = c("L.columbiana", "L.rugosa", "L.vulpina"),
                            L.lupina_ref = c(snpdistance(matrixl$L.lupina_culture, matrixl$L.columbiana), 
                                             snpdistance(matrixl$L.lupina_culture, matrixl$L.rugosa),
                                             snpdistance(matrixl$L.lupina_culture, matrixl$L.vulpina)),
                            L.lupina_alt = c(snpdistance(matrixl$L.lupina_alt, matrixl$L.columbiana),
                                             snpdistance(matrixl$L.lupina_alt, matrixl$L.rugosa),
                                             snpdistance(matrixl$L.lupina_alt, matrixl$L.vulpina)) )
  
  return(lethariadiv)
}

lethariadiv_all <- rbind(cbind(region = "polymorphic", prematrix2distance(prematrix)),
                         cbind(region = "fixed_alternative", prematrix2distance(prematrix_altfix)),
                         cbind(region = "fixed_culture", prematrix2distance(prematrix_reffix)) )

# Replace the alternative lupina for the parts that are fixed for the pure culture with NAs
lethariadiv_all <- lethariadiv_all %>% mutate(L.lupina_alt = ifelse(region == "fixed_culture", NA, L.lupina_alt))

cat("Plotting ...\n")
# Make a long format that ggplot likes
lethariadiv_all_long <- lethariadiv_all %>% gather(key = allele, value = Distance, -region, -Species )

distanceplot <- ggplot(lethariadiv_all_long, aes(x = Species, y = Distance, colour = allele, shape = allele)) + 
  geom_point(size = 3.5) + 
  ylab("Proportion of sites that have different alleles") +
  facet_grid(. ~ region) +
  theme_bw() + 
  ylim(0, 1) +
  scale_color_manual(values = c("darkgoldenrod3", "chartreuse4"), labels = c("other", "culture")) + # Lupina is green in the phylogenies
  scale_shape_manual(values = c(19, 15), labels = c("other", "culture"))

# ggsave(plot = distanceplot, "/Users/Lorena/Dropbox/PhD_UU/Analyses/Letharia/3_LethariaPloidy/results/Lichens-snps-miss1-100kp_MAF-noTEs_distance.pdf", width = 8.3, height = 3.2)
ggsave(plot = distanceplot, snakemake@output$distance, width = 8.3, height = 3.2)

cat("Done!")
