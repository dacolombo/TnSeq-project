################# Dependencies ###################
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")

library("GenomicRanges")
library("circular")
library("ggplot2")


############## Data preprocessing ################
# Import data
# 
setwd('/home/daniele/github/TnSeq-project')

fitness_data <- na.omit(read.table("data/Tn_seq_fitness_data_Opijnen_et_al_2009.txt",
                                   header=TRUE, stringsAsFactors = FALSE, sep ="\t"))

gene_coordinates <- read.table("data/GCF_000006885.1_ASM688v1_genomic_olt.txt",
                               header=FALSE, stringsAsFactors=FALSE, sep = "\t",
                               col.names=c('locus','start','end'))

# Merge data
merged_fitness <- merge(fitness_data, gene_coordinates, by='locus')

# Convert data into genomic ranges
fitness_gr <- GRanges(
    seqnames = Rle("chr", nrow(merged_fitness)),
    ranges = IRanges(start=merged_fitness$start, end=merged_fitness$end, names=merged_fitness$locus),
    avg_fitness = merged_fitness$average_fitness)


# Add gene categories based on each gene fitness
fitness_breakpoints <- c(0,0.96,1.04,Inf)
fitness_categories <- c('Disadvantageous','Neutral','Advantageous')
fitness_gr$gene_category <- cut(fitness_gr$avg_fitness, breaks=fitness_breakpoints, labels=fitness_categories)


# Compute circular coordinates
genome_length <- max(end(fitness_gr))
fitness_gr$middle_coordinate <- (end(fitness_gr)-start(fitness_gr))/2 + start(fitness_gr)
fitness_gr$radians_coordinate <- fitness_gr$middle_coordinate*2*pi/genome_length
fitness_gr$degrees_coordinate <- fitness_gr$middle_coordinate*360/genome_length



compute_window_statistics <- function(gr, window_size=100000) {
  # Define vector containing window start coordinates
  windows_start <- seq(from=0, to=max(end(gr)), by=window_size)
  
  # Build GRanges object containing the windows
  window_gr <- GRanges(seqnames=Rle("chr", length(windows_start)),
                       ranges=IRanges(start=windows_start, width=window_size))
  
  # Retrieve genome length (used to convert coordinates in radians)
  genome_length <- max(end(gr))
  
  # Build GRangesList object: at each loop of the lapply we consider one window
  # and we compute the statistics for the group of genes overlapping with the
  # window.
  # The result would be a list with each element being a GRanges object
  # containing the genes overlapping with a specific window.
  
  window_statistics <- lapply(seq(1:length(window_gr)), function(X){
    # Select window
    window <- window_gr[X]
    
    # Subset by selected window
    window_overlapping_gr <- subsetByOverlaps(fitness_gr, window)
    
    # Compute window-related statistics
    window_middle_coord <- round((end(window)-start(window))/2 + start(window))
    window_middle_coord_radians <- window_middle_coord*2*pi/genome_length
    window_avg_fitness <- mean(window_overlapping_gr$avg_fitness)
    window_statistics <- data.frame(coordinate=window_middle_coord,
                                    radians=window_middle_coord_radians,
                                     fitness=window_avg_fitness)
     
    return(window_statistics)
    }
  )
  
  # Return dataframe
  return(do.call(rbind.data.frame, window_statistics))
}

window_statistics <- compute_window_statistics(fitness_gr)
window_statistics





############## Linear Model ################
# Fit linear model using radians
linear_model <- lm(formula=fitness ~ sin(radians)+cos(radians), data=window_statistics)
summary(linear_model)

# Predict window avg fitness
predicted_fitness <- predict(linear_model, window_statistics)
plot(window_statistics$coordinate, predicted_fitness)
plot(window_statistics$coordinate, window_statistics$fitness)

# Correct fitness by removing smile effect
corrected_fitness <- window_statistics$fitness - predicted_fitness + 1
plot(window_statistics$coordinate, corrected_fitness)


# Use the linear model on all the data
all_df <- data.frame(coordinate=fitness_gr$middle_coordinate,
                     radians=fitness_gr$radians_coordinate,
                     fitness=fitness_gr$avg_fitness)
corrected_fitness_all <- all_df$fitness - predict(linear_model, all_df) + 1
plot(all_df$coordinate, corrected_fitness_all)




############## Plots ################

# Plots
scatter.smooth(sinpi(df_linear_regression$radians), df_linear_regression$fitness)

circular::plot.circular(rep(df_linear_regression$radians, df_linear_regression$fitness^100), bins=24, stack=TRUE)

circular::rose.diag(rep(df_linear_regression$radians, df_linear_regression$fitness^100))

ggplot(df_linear_regression, aes(x=radians, y=fitness)) +
  geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  ylim(-1,1.5) +
  coord_polar(start = 0)





# Old stuff
smoothScatter(x=(geneCoord$V2[fData$average_fitness!=0]+geneCoord$V3[fData$average_fitness!=0])/2,
              y=fData$average_fitness[fData$average_fitness!=0],
              ylab = "avg gene fitness",xlab="genomic coordinates",main="fitness vs genome location")
