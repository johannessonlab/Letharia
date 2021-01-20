#!/usr/bin/env Rscript

### Exploring the repeat content of Letharia lupina
#############################################################################

# This script produces figures S5, S6 and S7. It uses the output of RepeatMasker 
# with the repeat library from McKenzie et al. 2020. It also uses tables of coverage
# produced from these same RepeatMasker files with my script totalcovergff.py, like so:

# $ python totalcovergff.py RepeatMasker/Lecol_assembly_V1.0.fa.out.gff -f McKenzieData/Lecol_assembly_V1.0.fasta > RepeatContent_Lecol_assembly_V1.0.txt

# Here McKenzieData/Lecol_assembly_V1.0.fasta is the corresponding genome assembly.
# The script is available here: https://github.com/SLAment/Genomics/blob/master/GenomeAnnotation/totalcovergff.py

# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-06-18 - 2020-07-13
# Version 1
# =======================================

library(ggplot2)
library(dplyr)

setwd("/path/to/repo") # set working directory

# ============================
# Functions
# ============================
readgff <- function(gff_file){
  gff <- read.table(gff_file) # a gtf actually
  names(gff) <- c("sequence", "source", "feature", "start", "end", "identity", "strand", "phase", "Target", "Motif", "start_cons", "end_cons")
  # Remove some irrelevant columns
  gff <- select(gff, -c("source", "feature", "phase", "Target"))
  
  # Add length of alignment
  gff <- gff %>% mutate(len_ref = end - start + 1) %>%
    mutate(Motif= gsub("Motif:", "", Motif)) %>% 
    filter(!grepl(")n", Motif)) %>% filter(!grepl("-rich", Motif))
  
  # Make it a factor, not a character
  gff$sequence <- factor(gff$sequence)
  gff$Motif <- factor(gff$Motif)
  
  return(gff)
}

read_totalcovergff <- function(tetable){
  repeatcontent <-  read.table(tetable, header = TRUE) %>% filter(Contig != "Total")
  ## Levels ordered by size of contig
  # https://rstudio-pubs-static.s3.amazonaws.com/7433_4537ea5073dc4162950abb715f513469.html
  repeatcontent$Contig <- factor(repeatcontent$Contig, levels = repeatcontent$Contig[order(repeatcontent$Len_ctg, decreasing = TRUE)])
  return(repeatcontent)
}

coverTE <- function(gff){
  # Let's make a new simplified data.frame with the coverage per contig and per Motif
  coverage <- data.frame()
  for (seq in levels(gff$sequence)) { # For every contig
    print(seq)
    # For each family
    for (motif in levels(gff$Motif)) {
      # How much does it cover 
      thisctggff <- gff %>% filter(sequence == seq, Motif == motif)
      totallen <- thisctggff$len_ref %>% sum
      coverage <- rbind(coverage, data.frame(sequence = seq, Motif = motif, cov = totallen))
    }
  }
  # This works nicely because even 0s are kept.
  return(coverage)
}

# ============================
# Data
# ============================
gff_file <- "data/Lelup_v1.1.fa.out.gff"
gff_file_pure <- "data/lupinapure.fa.out.gff"
gff_file_col <- "data/Lecol_assembly_V1.0.fa.out.gff"

# Read gffs
gff <- readgff(gff_file)
gff_pure <- readgff(gff_file_pure)
gff_col <- readgff(gff_file_col)

# ------

# Tables with repeat content summarized
repeatcontent_lelup <- read_totalcovergff("data/RepeatContent_Lelup_v1.1.txt")
repeatcontent_pure <- read_totalcovergff("data/RepeatContent_lupinapure.txt")
repeatcontent_col <- read_totalcovergff("data/RepeatContent_Lecol_assembly_V1.0.txt")


# ============================
# Plot Repeat content of Lelup
# ============================

Lelup_TEcov <- ggplot(repeatcontent_lelup, aes(x = Contig, y = Len_ctg, fill = coverage)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab(expression(paste("Scaffold in the ", italic("L. lupina"), " McKenzie assembly"))) + 
  ylab("Size of contig (bp)") + labs(fill='Repeat\ncontent') +
  scale_fill_continuous(type = "viridis", limits=c(0, 100))

ggsave(plot = Lelup_TEcov, "results/FigS5_RepeatCov_Lelup.tiff", width = 7, height = 3)

Lpure_TEcov <- ggplot(repeatcontent_pure %>% filter(Len_ctg >= 150000), aes(x = Contig, y = Len_ctg, fill = coverage)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5)) +
  xlab(expression(paste("Scaffold in the ", italic("L. lupina"), " pure culture assembly (larger than 150 kb)"))) + 
  ylab("Size of contig (bp)") + labs(fill='Repeat\ncontent') +
  scale_fill_continuous(type = "viridis", limits=c(0, 100))

ggsave(plot = Lpure_TEcov, "results/RepeatCov_Lpure.pdf", width = 10, height = 5)

Lcol_TEcov <- ggplot(repeatcontent_col %>% filter(Len_ctg >= 150000), aes(x = Contig, y = Len_ctg, fill = coverage)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5)) +
  xlab(expression(paste("Scaffold in the ", italic("L. lupina"), " McKenzie assembly"))) + 
  ylab("Size of contig (bp)") + labs(fill='Repeat\ncontent') +
  scale_fill_continuous(type = "viridis", limits=c(0, 100))

ggsave(plot = Lcol_TEcov, "results/RepeatCov_Lecol.pdf", width = 7, height = 3)


# ============================
# Play
# ============================

## Make a table with counts of occurance per Motif in the gff
motifcounts <- gff$Motif %>% table() %>% data.frame
names(motifcounts) <- c("Motif", "Freq")
# Reorder by abundance
motifcounts <- motifcounts[order(motifcounts$Freq, decreasing = TRUE),]

# ## How many Motifs are there in each contig
# gff %>% filter(sequence == "contig_1") %>% .$Motif %>% table() %>% data.frame %>% filter(Freq > 0) %>% dim
# gff %>% filter(sequence == "contig_2") %>% .$Motif %>% table() %>% data.frame %>% filter(Freq > 0) %>% dim
# gff %>% filter(sequence == "contig_3") %>% .$Motif %>% table() %>% data.frame %>% filter(Freq > 0) %>% dim
# gff %>% filter(sequence == "contig_4") %>% .$Motif %>% table() %>% data.frame %>% filter(Freq > 0) %>% dim

# Ok, so it's not like there are more things in contig_1

## Plot the abundance per contig but only with the abundant Motifs

# Remove the least frequent repeats
commonmotifs <- motifcounts %>% filter(Freq >= 200)

ggplot(gff %>% filter(Motif %in% commonmotifs$Motif), aes(x = len_ref, fill = Motif)) +
  geom_histogram() + 
  facet_wrap(~sequence, nrow = 5)

# Create a new data frame with total abundances for each motif (this takes a bit of time)
coverage <- coverTE(gff)

## Levels ordered by size of contig
coverage$sequence <- factor(coverage$sequence, levels = levels(repeatcontent_lelup$Contig))

## All contigs
ggplot(coverage %>% filter(Motif %in% commonmotifs$Motif), aes(x = Motif, y = cov, colour = Motif)) +
  geom_point() + 
  facet_wrap(~sequence, nrow = 6)

## Plot families that are in contig_1
coverage_ctg1 <- coverage %>% filter(sequence == "contig_1", cov != 0)
ggplot(coverage %>% filter(Motif %in% coverage_ctg1$Motif) %>% filter(Motif %in% commonmotifs$Motif), aes(x = Motif, y = cov, colour = Motif)) +
  geom_point() + 
  facet_wrap(~sequence, nrow = 6)

# So it seems that the TEs in contig_1 are definetely present in some of the other alleles. 

## What is the most abundant TE in contig_1?
coverage_ctg1 <- coverage_ctg1[order(coverage_ctg1$cov, decreasing = TRUE),]
# coverage_ctg1 %>% head
# sequence               Motif    cov
# 200 contig_1 lelup-families4-970 284193
# 176 contig_1  lelup-families3-45 242492
# 169 contig_1  lelup-families2-10 220853
# 186 contig_1 lelup-families4-322 201266
# 197 contig_1 lelup-families4-693 194941
# 183 contig_1  lelup-families4-27 160850

# Ok, so there are several, not like one dominant guy. How abundant are those on the other contigs?

topTEcontig1 <- ggplot(coverage %>% filter(Motif %in% coverage_ctg1$Motif[1:5]), aes(x = sequence, y = cov, colour = Motif, shape = Motif)) +
  geom_point(size = 2) + 
  xlab(expression(paste("Scaffold in the ", italic("L. lupina"), " McKenzie assembly"))) + 
  ylab("Repeat coverage (bp)") + labs(colour='Top repeats\nin contig_1', shape = 'Top repeats\nin contig_1') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(plot = topTEcontig1, "results/FigS7_TopRepeatsContig_1.tiff", width = 9, height = 4)

## Are these same contigs in the pure culture?

# First let's get the big contigs
bigcontigs <- repeatcontent_pure %>% filter(Len_ctg >= 1000000) %>% .$Contig %>% factor
gff_pure_filter <- gff_pure %>% filter(sequence %in% bigcontigs)
# # Fix the levels to reflect the order of alignment to the reference
gff_pure_filter$sequence <- factor(gff_pure_filter$sequence, levels = gff_pure_filter$sequence %>% unique %>% rev %>% droplevels)

# Create a new data frame with total abundances for each motif
coverage_pure <- coverTE(gff_pure_filter)

## Levels ordered by size of contig
coverage_pure$sequence <- factor(coverage_pure$sequence, levels = levels(repeatcontent_pure$Contig))


topTEcontig1_pure <- ggplot(coverage_pure %>% filter(Motif %in% coverage_ctg1$Motif[1:5]), aes(x = sequence, y = cov, colour = Motif, shape = Motif)) +
  geom_point(size = 2) + 
  xlab(expression(paste("Scaffold in the ", italic("L. lupina"), " pure-culture assembly"))) + 
  ylab("Repeat coverage (bp)") + labs(colour='Top repeats\nin contig_1', shape = 'Top repeats\nin contig_1') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))

ggsave(plot = topTEcontig1_pure, "results/FigS6_TopRepeatsContig_1_pure.tiff", width = 9, height = 4)

## What about the L. columbiana assembly?

# First let's get the big contigs
bigcontigs_col <- repeatcontent_col %>% filter(Len_ctg >= 150000) %>% .$Contig %>% factor
gff_col_filter <- gff_col %>% filter(sequence %in% bigcontigs_col)
# # Fix the levels to reflect the order of alignment to the reference
gff_col_filter$sequence <- factor(gff_col_filter$sequence, levels = gff_col_filter$sequence %>% unique %>% rev %>% droplevels)

# Create a new data frame with total abundances for each motif
coverage_col <- coverTE(gff_col_filter)

## Levels ordered by size of contig
coverage_col$sequence <- factor(coverage_col$sequence, levels = levels(repeatcontent_col$Contig))


topTEcontig1_col <- ggplot(coverage_col %>% filter(Motif %in% coverage_ctg1$Motif[1:5]), aes(x = sequence, y = cov, colour = Motif, shape = Motif)) +
  geom_point(size = 2) + 
  xlab(expression(paste("Scaffold in the ", italic("L. lupina"), " pure-culture assembly"))) + 
  ylab("Repeat coverage (bp)") + labs(colour='Top repeats\nin contig_1', shape = 'Top repeats\nin contig_1') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))

ggsave(plot = topTEcontig1_col, "results/TopRepeatsContig_1_lecol.pdf", width = 9, height = 4)

# Done!
