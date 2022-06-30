################# Dependencies ###################
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")

library("GenomicRanges")
#library("circular")
library("cowplot")
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
fitness_gr$middle_coordinate <- (end(fitness_gr)+start(fitness_gr))/2
fitness_gr$radians_coordinate <- fitness_gr$middle_coordinate*2*pi/genome_length
fitness_gr$degrees_coordinate <- fitness_gr$middle_coordinate*360/genome_length



############## Compute windows ################

compute_middle_coord <- function(gr, gnome_length) {
  start <- min(start(gr))
  end <- max(end(gr))
  middle_coord <- round((end+start)/2)
  middle_coord_radians <- middle_coord*2*pi/genome_length
  return(c(middle_coord, middle_coord_radians))
}

compute_avg_fitness <- function(gr) {
  avg_fitness <- mean(gr$avg_fitness)
  return(avg_fitness)
}

compute_window_statistics <- function(gr, window_size=100000) {
  # Define vector containing window start coordinates
  windows_start <- seq(from=min(start(gr)), to=max(end(gr)), by=window_size)
  
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
  
  window_statistics <- lapply(X=seq(1:length(window_gr)), FUN=function(X){
    # Select window
    window <- window_gr[X]
    
    # Subset by selected window
    window_overlapping_gr <- subsetByOverlaps(gr, window)
    
    # Compute window-related statistics
    window_coordinates <- compute_middle_coord(window, genome_length)
    window_avg_fitness <- compute_avg_fitness(window_overlapping_gr)
    return(c(window_coordinates, window_avg_fitness))
    }
  )
  
  # Return dataframe
  df <- as.data.frame(do.call(rbind, window_statistics))
  colnames(df) <- c('coordinate','radians','fitness')
  # fitness_breakpoints <- c(0,0.96,1.04,Inf)
  # fitness_categories <- c('Disadvantageous','Neutral','Advantageous')
  # df$category_before <- cut(df$fitness, breaks=fitness_breakpoints, labels=fitness_categories)
  return(df)
}

window_statistics <- compute_window_statistics(fitness_gr)
window_statistics



############## Compute groups ################

compute_group_statistics <- function(gr, group_size=100) {
  # Divide GRanges object into groups of 200 genes
  number_of_genes <- length(gr)
  number_of_groups <- ceiling(number_of_genes/group_size)
  gr_groups <- split(gr, rep(1:number_of_groups, each=group_size, length.out=number_of_genes))

  # Retrieve genome length (used to convert coordinates in radians)
  genome_length <- max(end(gr))

  # Compute statistics for each group of genes
  group_statistics <- lapply(seq(1:number_of_groups), function(X){
    group <- unlist(gr_groups[X])
    group_coordinates <- compute_middle_coord(group, genome_length)
    group_avg_fitness <- compute_avg_fitness(group)
    return(c(group_coordinates, group_avg_fitness))
  })

  # Return dataframe
  df <- as.data.frame(do.call(rbind, group_statistics))
  colnames(df) <- c('coordinate','radians','fitness')
  return(df)
}

group_statistics <- compute_group_statistics(fitness_gr)
group_statistics



############## Linear Model: coordinates ################
# Fit linear model using half of the genome
first_half_gr <- subsetByOverlaps(fitness_gr, GRanges('chr', IRanges(0, genome_length/2)))
first_half_df_windows <- compute_window_statistics(first_half_gr, window_size=100000)

second_half_gr <- subsetByOverlaps(fitness_gr, GRanges('chr',IRanges(genome_length/2, genome_length)))
second_half_df_windows <- compute_window_statistics(second_half_gr, window_size=100000)

# Fit model on the first half
lm_linear_coord <- lm(formula=fitness ~ coordinate, data=first_half_df_windows)
summary(lm_linear_coord)

# Predict first half (train)
first_half_pred_windows <- predict(lm_linear_coord, first_half_df_windows)
mean((first_half_pred_windows-first_half_df_windows$fitness)^2)

# Predict second half (test)
a <- coef(lm_linear_coord)['(Intercept)']
b <- coef(lm_linear_coord)['coordinate']
second_half_pred_windows <- a - b*(second_half_df_windows$coordinate - genome_length)
mean((second_half_pred_windows - second_half_df_windows$fitness)^2)

# Correct fitness values based on linear model
corrected_first_half_windows <- first_half_df_windows$fitness - first_half_pred_windows + 1
corrected_second_half_windows <- second_half_df_windows$fitness - second_half_pred_windows + 1

# Show correction on windows
window_statistics_linear <- rbind(first_half_df_windows, second_half_df_windows)
correction_plots_windows(window_statistics_linear, x='coordinate')


# Use the linear model on all the data
first_half_df_all <- data.frame(coordinate=first_half_gr$middle_coordinate,
                                fitness=first_half_gr$avg_fitness)


first_half_pred_all <- predict(lm_linear_coord, first_half_df_all)
first_half_df_all$corrected_fitness <- first_half_df_all$fitness - first_half_pred_all + 1



second_half_df_all <- data.frame(coordinate=second_half_gr$middle_coordinate,
                                 fitness=second_half_gr$avg_fitness)
second_half_pred_all <- a - b*(second_half_df_all$coordinate - genome_length)
second_half_df_all$corrected_fitness <- second_half_df_all$fitness - second_half_pred_all + 1

all_df_linear <- rbind(first_half_df_all, second_half_df_all)
all_df_linear$category_before <- fitness_gr$gene_category
all_df_linear$category_after <- cut(all_df_linear$corrected_fitness, breaks=fitness_breakpoints, labels=fitness_categories)
summary(all_df_linear$category_before)
summary(all_df_linear$category_after)


# Plot before and after correction
correction_plots_windows <- function(df, x) {
  
  plot1 <- ggplot(df, aes(x=.data[[x]], y=fitness)) +
    geom_point() +
    geom_smooth() +
    ylim(0.9, 1.1) +
    theme(plot.margin=unit(c(3,0,0,1), "lines"))
  
  plot2 <- ggplot(df, aes(x=.data[[x]], y=corrected_fitness)) +
    geom_point() +
    geom_smooth() +
    ylim(0.9, 1.1) +
    theme(plot.margin=unit(c(3,0,0,1), "lines"))
  
  plot_grid(plot1, plot2, ncol=2, label_x=0,
            labels=c('Before correction', 'After correction'))
}

correction_plots <- function(df, x) {
  # First scatterplot, before correction
  plot1 <- ggplot(df, aes(x=.data[[x]], y=fitness)) + 
    geom_point(aes(color=category_before)) +
    geom_smooth(color='black') +
    ylim(0.8, 1.2) +
    theme(plot.margin = unit(c(3,0,0,1), "lines"))
  
  # Second scatterplot, after correction
  plot2 <- ggplot(df, aes(x=.data[[x]], y=corrected_fitness)) + 
    geom_point(aes(color=category_after)) +
    geom_smooth(color='black') +
    ylim(0.8, 1.2) +
    theme(plot.margin = unit(c(3,0,0,1), "lines"))
  
  # First barplot, before correction
  plot3 <- ggplot(df, aes(x=category_before, fill=category_before)) +
    geom_bar(stat="count") +
    ylim(0,1300) +
    theme(plot.margin = unit(c(3,0,0,1), "lines"))
  
  # Second barplot, after correction
  plot4 <- ggplot(df, aes(x=category_after, fill=category_after)) +
    geom_bar(stat="count") +
    ylim(0,1300) +
    theme(plot.margin = unit(c(3,0,0,1), "lines"))
  
  # Show plots
  plot_grid(plot1, plot2, plot3, plot4, ncol=2, nrow=2,
            labels=c('Before correction', 'After correction'), label_x=0)
}


correction_plots(all_df_linear,x='coordinate')




############# Linear Model: radians ################
# Fit linear model using radians
window_statistics_radians <- window_statistics
lm_radians <- lm(formula=fitness ~ sin(radians) + cos(radians), data=window_statistics_radians)
summary(lm_radians)

# Correct windows fitness by removing smile effect
predicted_fitness <- predict(lm_radians, window_statistics_radians)
window_statistics_radians$corrected_fitness <- window_statistics_radians$fitness - predicted_fitness + 1
correction_plots_windows(window_statistics_radians, x='radians')
correction_plots_windows(window_statistics_radians, x='coordinate')


# Use the linear model on all the data
all_df_radians <- data.frame(radians=fitness_gr$radians_coordinate,
                             fitness=fitness_gr$avg_fitness,
                             category_before=factor(fitness_gr$gene_category))
all_df_radians$corrected_fitness <- all_df_radians$fitness - predict(lm_radians, all_df_radians) + 1
all_df_radians$category_after <- cut(all_df_radians$corrected_fitness, breaks=fitness_breakpoints, labels=fitness_categories)
summary(all_df_radians$category_before)
summary(all_df_radians$category_after)

# Plot before and after
correction_plots(all_df_radians,x='radians')




############## Plots ################

# Plots
scatter.smooth(sinpi(df_linear_regression$radians), df_linear_regression$fitness)

circular::plot.circular(rep(df_linear_regression$radians, df_linear_regression$fitness^100), bins=24, stack=TRUE)

circular::rose.diag(rep(df_linear_regression$radians, df_linear_regression$fitness^100))

ggplot(df_linear_regression, aes(x=radians, y=fitness)) +
  geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  ylim(-1,1.5) +
  coord_polar(start = 0)

window_plot_df <- compute_window_statistics(fitness_gr, 5000)
window_plot_df$gene_category <- cut(window_plot_df$fitness, breaks=fitness_breakpoints, labels=fitness_categories)
window_plot_df$corrected_fitness <- window_plot_df$fitness - predict(lm_radians, window_plot_df) + 1
window_plot_df$corrected_gene_category <- cut(window_plot_df$corrected_fitness, breaks=fitness_breakpoints, labels=fitness_categories)


# PLOT SHOWING SMILE EFFECT WITH CIRCULAR COORDINATES #
library(scales)
pi_scales <- math_format(.x * pi, format = function(x) x / pi)


# Plot before correction
ggplot(window_plot_df) + 
  geom_bar(aes(x=radians, y=fitness, fill=gene_category), stat="identity") + 
  ylim(-0.5,1.3) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  coord_polar(start=0) +
  scale_x_continuous(labels=pi_scales, breaks=seq(0, 2*pi, pi/2))
  

# Plot after correction
ggplot(window_plot_df) + 
  geom_bar(aes(x=radians, y=corrected_fitness, fill=corrected_gene_category), stat="identity") + 
  ylim(-0.5,1.3) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  coord_polar(start=0) +
  scale_x_continuous(labels=pi_scales, breaks=seq(0, 2*pi, pi/2))





# Old stuff
smoothScatter(x=(geneCoord$V2[fData$average_fitness!=0]+geneCoord$V3[fData$average_fitness!=0])/2,
              y=fData$average_fitness[fData$average_fitness!=0],
              ylab = "avg gene fitness",xlab="genomic coordinates",main="fitness vs genome location")
