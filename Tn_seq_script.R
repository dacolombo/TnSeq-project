################# Dependencies ###################
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")

library("GenomicRanges")


############## Data preprocessing ################
# Import data
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


# Old stuff
smoothScatter(x=(geneCoord$V2[fData$average_fitness!=0]+geneCoord$V3[fData$average_fitness!=0])/2,
              y=fData$average_fitness[fData$average_fitness!=0],
              ylab = "avg gene fitness",xlab="genomic coordinates",main="fitness vs genome location")
