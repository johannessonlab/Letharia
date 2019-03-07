#!/usr/bin/env Rscript

# Alignment of Letharia MAT vs the outgroup genera
# =======================================
# NUCmer was run like such:
# $ nucmer -b 200 -c 20 -p MelaneliaVsRugosa_nc /Users/Lorena/Dropbox/PhD_UU/Analyses/Letharia/FinalMATannotation/rugosa/MAT_rugosa.fa Melanelia_stygia_MAT.fas --maxmatch
# $ show-coords -r -B MelaneliaVsRugosa_nc.delta > MelaneliaVsRugosa_nc.tab
# 
# $ nucmer -b 200 -c 20 -p HypogymniaVsColumbiana_nc /Users/Lorena/Dropbox/PhD_UU/Analyses/Letharia/FinalMATannotation/columbiana/MAT_columbiana.fa  Hypogymnia_subphysodes_MAT.fas --maxmatch
# $ show-coords -r -B HypogymniaVsColumbiana_nc.delta > HypogymniaVsColumbiana_nc.tab
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-02-08
# =======================================
Version <- 1.0
# ============================
# Load the necessary libraries
# ============================
library(ggplot2, quietly = TRUE)
library(grid)
library(dplyr, warn.conflicts = FALSE)

# MAT1-1 (Melanelia stygia)
fileMAT1 <- "MelaneliaVsRugosa_nc.tab"
coordsMAT1 <- read.delim(fileMAT1, header = FALSE)
names(coordsMAT1) <- c("QUERY","DATE","LEN_Q","PROGRAM","REF_FILE","REF","S2","E2","S1","E1","IDY","SYM","LEN_2","V14","V15","V16","V17","SENSE","LEN_R","V20","V21")

# MAT1-2 (Hypogymnia subphysodes)
fileMAT2 <- "HypogymniaVsColumbiana_nc.tab"
coordsMAT2 <- read.delim(fileMAT2, header = FALSE)
names(coordsMAT2) <- c("QUERY","DATE","LEN_Q","PROGRAM","REF_FILE","REF","S2","E2","S1","E1","IDY","SYM","LEN_2","V14","V15","V16","V17","SENSE","LEN_R","V20","V21")

# MAT1-2 (Cetraria islandica)
fileMAT2cet <- "CetrariaVsColumbiana_nc.tab"
coordsMAT2cet <- read.delim(fileMAT2cet, header = FALSE)
names(coordsMAT2cet) <- c("QUERY","DATE","LEN_Q","PROGRAM","REF_FILE","REF","S2","E2","S1","E1","IDY","SYM","LEN_2","V14","V15","V16","V17","SENSE","LEN_R","V20","V21")

# ============================
# Processing
# ============================
# Prepare a function to produce the DotPlots from a given data set

columbiana <- data.frame(Gene= c("APN2", "lorf", "MAT1-2-1", "MAT1-2-L", "pseudoalpha", "SLA2"),
                         S1 = c(351, 3042, 3567, 5198, 6647, 8787),
                         end = c(2631, 3416, 4801, 5969, 6901, 12717),
                         LEN_Q = 0
                         )

rugosa <- data.frame(Gene= c("APN2", "lorf", "MAT1-1-L", "MAT1-1-1", "SLA2"),
                         S1 = c(351, 3066, 4075, 5755, 8811),
                         end = c(2650, 3507, 5427, 7021, 12743),
                         LEN_Q = 0
)

# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# ============================
# Plotting
# ============================

lineplotter <- function(data, multiref = FALSE, linesize = 1, forwardcol = "coral2", revcol = "cyan3"){ 
  # As it is, it can only take one contig alignment
  df <- data.frame()
  
  # Make an empty plot, with panels of height equal to lenght of queries
  query <- as.character(data$QUERY[1]) # Name of the query
  reference <- as.character(data$REF[1])
  
  p <- ggplot(data, aes(x=S1, y=1)) + geom_blank() + ylab(query) + xlab(reference) + #+ ylim(0, max(data$LEN_Q)) # Make axis the length of that query sequence
    geom_hline(yintercept=1) + geom_hline(yintercept=2) 
    
  # Plot the alignments 
  for (i in 1:dim(data)[1]){
    currentrow <- data[i,]
    
    if (currentrow$S2 > currentrow$E2){ # Reverse
      colorsito <- revcol
      currentdf <- data.frame(S1 = c(currentrow$S1, currentrow$E1, currentrow$S1, currentrow$E1), y = c(1, 2, 2, 1), col = colorsito) # Reverse order of coordinates
    }
    else{ # Forward
      colorsito <- forwardcol
      currentdf <- data.frame(S1 = c(currentrow$S1, currentrow$S1, currentrow$E1, currentrow$E1), y = c(1, 2, 2, 1), col = colorsito)
    }
    
    # Record the data for the polygons
    df <- rbind(df, currentdf)
    
    # Plot the alignment
    p <- p + geom_segment(data = currentrow, aes(x = S1, xend = E1, y = 2, yend = 2), colour = colorsito, size = linesize) +
      geom_point(data = currentrow, aes(x=S1, y = 2), colour = colorsito, size = linesize) + geom_point(data = currentrow, aes(x=E1, y = 2), colour = colorsito, size = linesize) 
  }
  
  # Make a polygons
  p <- p + geom_polygon(data = df, aes(x = S1, y = y, fill = alpha(col, 0.5))) + scale_fill_identity() # scale_fill_identity() equates the group to the color

  return(p)
}

hypogymnia <- lineplotter(coordsMAT2, linesize = 1.5, forwardcol = "cyan3") + theme_bw() + # White background
  theme(strip.text.y = element_blank(), 
        panel.border=element_blank(),  # Remove the gray label on the side and border of plot
        legend.position="none", 
        axis.title.y=element_blank(), # Remove y axis text
        axis.text.y=element_blank(),
        axis.ticks.y =element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
        ) + 
  geom_segment(data = columbiana, aes(x = S1, xend = end, yend = 1, colour = Gene), size = 5) +
  labs( x = expression(italic("Letharia columbiana")~" MAT1-2 (bp)")) +
  ggtitle(expression(italic("Hypogymnia subphysodes")~" MAT1-2")) +
  scale_color_manual(values=c("black", "#666666ff", "#5f8dd3ff", "#006680ff", "pink", "black")) + 
  annotate(geom="text", x=1500, y=1, label="APN2", color="white", size = 2) +
  annotate(geom="text", x=3230, y=1, label="lorf", color="white", size = 2) +
  annotate(geom="text", x=4200, y=1, label="MAT1-2-1", color="black", size = 2) +
  annotate(geom="text", x=5580, y=1, label="MAT1-2-L", color="white", size = 2) +
  annotate(geom="text", x=6780, y=0.95, label="*", color="black", size = 6) +
  annotate(geom="text", x=10800, y=1, label="SLA2", color="white", size = 2) 

cetraria <- lineplotter(coordsMAT2cet, linesize = 1.5, forwardcol = "cyan3") + theme_bw() + # White background
  theme(strip.text.y = element_blank(), 
        panel.border=element_blank(),  # Remove the gray label on the side and border of plot
        legend.position="none", 
        axis.title.y=element_blank(), # Remove y axis text
        axis.text.y=element_blank(),
        axis.ticks.y =element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
  ) + 
  geom_segment(data = columbiana, aes(x = S1, xend = end, yend = 1, colour = Gene), size = 5) +
  labs( x = expression(italic("Letharia columbiana")~" MAT1-2 (bp)")) +
  ggtitle(expression(italic("Cetraria islandica")~" MAT1-2")) +
  scale_color_manual(values=c("black", "#666666ff", "#5f8dd3ff", "#006680ff", "pink", "black")) + 
  annotate(geom="text", x=1500, y=1, label="APN2", color="white", size = 2) +
  annotate(geom="text", x=3230, y=1, label="lorf", color="white", size = 2) +
  annotate(geom="text", x=4200, y=1, label="MAT1-2-1", color="black", size = 2) +
  annotate(geom="text", x=5580, y=1, label="MAT1-2-L", color="white", size = 2) +
  annotate(geom="text", x=6780, y=0.95, label="*", color="black", size = 6) +
  annotate(geom="text", x=10800, y=1, label="SLA2", color="white", size = 2) 

melanelia <- lineplotter(coordsMAT1, linesize = 1.5, forwardcol = "cyan3") + theme_bw() + # White background
  theme(strip.text.y = element_blank(), 
        panel.border=element_blank(),  # Remove the gray label on the side and border of plot
        legend.position="none", 
        axis.title.y=element_blank(), # Remove y axis text
        axis.text.y=element_blank(),
        axis.ticks.y =element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
  ) + 
  geom_segment(data = rugosa, aes(x = S1, xend = end, yend = 1, colour = Gene), size = 5) +
  labs( x = expression(italic("Letharia rugosa")~" MAT1-1 (bp)")) +
  ggtitle(expression(italic("Melanelia stygia")~" MAT1-1")) +
  scale_color_manual(values=c("black", "#666666ff", "#501644ff", "#800000ff", "black")) + 
  annotate(geom="text", x=1500, y=1, label="APN2", color="white", size = 2) +
  annotate(geom="text", x=3300, y=1, label="lorf", color="white", size = 2) +
  annotate(geom="text", x=4750, y=1, label="MAT1-1-L", color="white", size = 2) +
  annotate(geom="text", x=6400, y=1, label="MAT1-1-1", color="white", size = 2) +
  annotate(geom="text", x=10800, y=1, label="SLA2", color="white", size = 2) 

pdf("SyntenyLetharia.pdf", width = 8, height = 6)
multiplot(hypogymnia, cetraria, melanelia, cols = 1) 
dev.off()
