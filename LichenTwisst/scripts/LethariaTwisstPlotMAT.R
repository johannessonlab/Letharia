#!/usr/bin/env Rscript

### Plot topologies along the Letharia contigs
#############################################################################
# Part of the LichenPloidy.smk pipeline
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2021-03-08
# Version 1
# =======================================

library(ggplot2)
library(dplyr)
library(tidyr) # For separate

############################## input files ######################################
#coordinates file for each window plus topologies
window_data_file <- snakemake@input$windowdata

# The lupina MG SNPs
LlupinaMG_SNPs <- snakemake@input$snpsMGdf

# The annotation of the MAT locus
gff_file <- snakemake@input$gff

# The contig containing the MAT locus
matctg <- snakemake@params$matctg 

################################# Plot contig ##################################

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
        nsnps_win <- currentwin$cov_allele %>% length() #NEW
        
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

# Gene annotation
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

# ============================
## Plot the allele frequencies along the mating type scaffold in windows
# ============================
cat("Processing allele frequencies along the MAT contig ...\n")

# Clean a bit further
window_data <- read.table(window_data_file, header = TRUE) %>% 
  #filter(scaffold == matctg) %>% 
  mutate(sizewin = end - start) %>% 
  filter(sizewin <= 10000) # Windows of 50 SNPs that are larger than that probably have too much missing data

# Read the SNPs data set
vcfdf_lupinaMG_df <- read.table(LlupinaMG_SNPs, header = TRUE)

# Make a data frame with the contigs and their lengths
relevantchr <- data.frame(CHROM = matctg) # MAT contig

chrlines <- data.frame(CHROM = relevantchr$CHROM, len = relevantchr %>% separate(CHROM, into = c(NA, NA, NA, "len", NA, NA), sep ="_") %>% .$len %>% as.numeric)

# Filter for the contig of interest
vcfdf_lupinaMG_df_mat <- vcfdf_lupinaMG_df %>% 
  filter(CHROM %in% relevantchr$CHROM) %>%
  filter(maf >= 0.01) # Remove obvious sequencing errors

## Calculate the windows
## I want to have both alleles in the same plot
winmaf_lupina_mat_alt <- winmafperwin(vcfdf_lupinaMG_df_mat %>% filter(allele == "ALT"), 
                                      windowsize = 2500, chrlines = chrlines)

winmaf_lupina_mat_ref <- winmafperwin(vcfdf_lupinaMG_df_mat %>% filter(allele == "REF"), 
                                      windowsize = 2500, chrlines = chrlines)
# Put them in the same data.frame
winmaf_lupina_mat <- rbind(cbind(winmaf_lupina_mat_alt, allele = "ALT"), cbind(winmaf_lupina_mat_ref, allele = "REF"))

# Genes in the idiomorph
gff <- read.table(gff_file) %>% filter(V3 == "gene") %>% data.frame %>% gff2genenames() %>% filter(!genenames %in% c("APN2", "SLA2"))

# Make a data frame with the gff
bottomy <- -0.03 # A lower position for the genes
genes <- data.frame(allele = "REF", x1 = gff$V4, y1 = bottomy, x2 = gff$V5, y2 = bottomy)

toposbottomy <- -0.1
toposbottomy2 <- toposbottomy - 0.09
# ----
cat("Plotting ...\n")

# Finally plot
vcfdf_lupinaMG_mat <- ggplot(winmaf_lupina_mat, aes(x = SITE, y = cov_allele, colour = allele)) +
  geom_line(size = 0.7, alpha = 0.8) +
  theme_bw() + xlab("Scaffold position") + ylab("Allele frequency in reads") +
  geom_segment(data = genes, aes(x = x1, y = y1, xend = x2, yend = y2), size = 1.5, colour = "black") + # draw the genes
  annotate("text", x = genes$x1[1] + 2000 + (genes[length(genes$x1),2] - genes$x1[1])/2, y = bottomy + 0.08, label = "MAT idiomorph", size = 2.5) +
  geom_point(data = vcfdf_lupinaMG_df_mat, aes(x = SITE, y = cov_allele), alpha = 0.5, size = 0.5) +
  scale_y_continuous(breaks=c(0, 0.33, 0.66, 1)) + # Change the ticks to match the expectations +
  scale_colour_manual(values = c("darkgoldenrod3", "chartreuse4"), labels = c("other", "lupina")) # Lupina is green in the phylogenies


## Plot the segments with colors representing the sister species of the two lupinas in a given window
# I had to do it in this awkward way because sometimes the data frames don't 
# have those categories, and if I use the topology in aes then it appears in the legend, 
# mixed with the colors of the SNP frequencies.I couldn't figure out how to prevent that.
if ("vulpina" %in% window_data$sisMG) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sisMG == "vulpina"), aes(x = start, y = toposbottomy, xend = end, yend = toposbottomy), size = 5, colour = "#0075DC") #blue
}
if ("lupina" %in% window_data$sisMG) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sisMG == "lupina"), aes(x = start, y = toposbottomy, xend = end, yend = toposbottomy), size = 5, colour = "chartreuse4")
}
if ("rugosa" %in% window_data$sisMG) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sisMG == "rugosa"), aes(x = start, y = toposbottomy, xend = end, yend = toposbottomy), size = 5, colour = "#FF5005") #carmin
}   
if ("other" %in% window_data$sisMG) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sisMG == "other"), aes(x = start, y = toposbottomy, xend = end, yend = toposbottomy), size = 5, colour = "gray")
}  

if ("vulpina" %in% window_data$sislupina) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sislupina == "vulpina"), aes(x = start, y = toposbottomy2, xend = end, yend = toposbottomy2), size = 5, colour = "#0075DC") #blue
}
if ("MG" %in% window_data$sislupina) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sislupina == "MG"), aes(x = start, y = toposbottomy2, xend = end, yend = toposbottomy2), size = 5, colour = "darkgoldenrod3")
}
if ("rugosa" %in% window_data$sislupina) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sislupina == "rugosa"), aes(x = start, y = toposbottomy2, xend = end, yend = toposbottomy2), size = 5, colour = "#FF5005") #carmin
}
if ("other" %in% window_data$sislupina) {
  vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + geom_segment(data = cbind(window_data, allele = "REF") %>% filter(sislupina == "other"), aes(x = start, y = toposbottomy2, xend = end, yend = toposbottomy2), size = 5, colour = "gray")
}   

vcfdf_lupinaMG_mat <- vcfdf_lupinaMG_mat + annotate("text", x = 9200, y = toposbottomy+0.01, label = 'Sister of "Other"', size = 3.3) +
  annotate("text", x = 9200, y = toposbottomy2+0.01, label = "Sister of pure culture", size = 3.3) 

# Finally, save the plot
ggsave(plot = vcfdf_lupinaMG_mat, snakemake@output$mat, width = 10, height = 3.1)

