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

fitness_data <- na.omit(read.table("data/Tn_seq_fitness_data_Opijnen_et_al_2009.txt", header=TRUE, stringsAsFactors = FALSE, sep ="\t"))

gene_coordinates <- read.table("data/GCF_000006885.1_ASM688v1_genomic_olt.txt", header=FALSE, stringsAsFactors=FALSE, sep = "\t", col.names=c('locus','start','end'))

# Merge data
merged_fitness <- merge(fitness_data, gene_coordinates, by='locus')

# Convert data into genomic ranges
fitness_gr <- GRanges(
    seqnames = Rle("chr", nrow(merged_fitness)),
    ranges = IRanges(start=merged_fitness$start, end=merged_fitness$end, names=merged_fitness$locus),
    avg_fitness = merged_fitness$average_fitness)

# Build window of 100000 nucleotides
window_size <- 100000
windows_start <- seq(from=0, to=max(end(fitness_gr)), by=window_size)
window_gr <- GRanges(seqnames = Rle("chr", length(windows_start)),
                     ranges = IRanges(start=windows_start, width=window_size))

# Compute average fitness in each window of 100000 nucleotides
window_separated_gr <- lapply(seq(1:length(window_gr)), function(X){
                                  window_ranges <- subsetByOverlaps(fitness_gr, window_gr[X])
                                  window_ranges$window_avg_fitness <- mean(window_ranges$avg_fitness)
                                  return(window_ranges)
                     }
)

window_separated_gr <- GRangesList(window_separated_gr)



# Update fitness_gr with new window_avg_fitness metadata column,
# taking unique ranges since one range could overlap two windows
fitness_gr <- unique(unlist(window_separated_gr))


# Add gene categories based on each gene fitness
fitness_breakpoints <- c(0,0.96,1.04,Inf)
fitness_categories <- c('Disadvantageous','Neutral','Advantageous')
fitness_gr$gene_category <- cut(fitness_gr$avg_fitness, breaks=fitness_breakpoints, labels=fitness_categories)


############## Data circularization ################

# Compute middle coordinate of each gene in radians
# Proportion:
# given l = length of the circular genome
#       p = position to convert in radians
# The following proportion must be solved for x:
# l : 2*pi = p : x
# Obtaining the following formula:
# x = p * 2*pi/l
genome_length <- max(end(fitness_gr))
middle_coordinates <- (end(fitness_gr)-start(fitness_gr))/2 + start(fitness_gr)
fitness_gr$radians_coordinate <- middle_coordinates*2*pi/genome_length
fitness_gr$degrees_coordinate <- middle_coordinates*360/genome_length



# Extract X (middle position of window) and Y (avg fitness) values for linear regression model
windows_start <- min(start(window_separated_gr))
windows_end <- max(end(window_separated_gr))
window_middle_position <- round(((windows_end - windows_start)/2) + windows_start)
window_middle_position_radians <- window_middle_position*2*pi/genome_length
window_middle_position_degrees <- window_middle_position*360/genome_length

df_linear_regression <- data.frame(radians=window_middle_position_radians,
                                   degrees=window_middle_position_degrees,
                                   fitness=unique(fitness_gr$window_avg_fitness))


############## Linear Model ################
# Fit linear model using radians
linear_model <- lm(formula=fitness ~ sin(radians)+cos(radians), data=df_linear_regression)
summary(linear_model)

# Predict window avg fitness
predicted_fitness <- predict(linear_model, df_linear_regression)
plot(window_middle_position, predicted_fitness)
plot(window_middle_position, df_linear_regression$fitness)

# Correct fitness by removing smile effect
corrected_fitness <- df_linear_regression$fitness - predicted_fitness + 1
plot(window_middle_position, corrected_fitness)


# Use the linear model on all the data
df_linear_regression_all <- data.frame(radians=fitness_gr$radians_coordinate,
                                       fitness=fitness_gr$avg_fitness)
corrected_fitness_all <- fitness_gr$avg_fitness - predict(linear_model, df_linear_regression_all) + 1
plot(middle_coordinates, corrected_fitness_all)




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
