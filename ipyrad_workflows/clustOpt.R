# FINAL R SCRIPT FOR ALL CLUSTOPT ANALYSES - Chris McDaniels 2025-09-11
# We will be using the downloaded files to look at 5 of 7 metrics as described in the paper "An empirical pipeline for choosing the optimal clustering threshold in RADseq studies" by McCartney-Melstad et al. More details can be seen here: https://doi.org/10.1111/1755-0998.13029 We will be using code and modified code from their associated GitHub repository: https://github.com/atcg/clustOpt
# Some of the modifications were made with the help of ChatGPT.

# Files you will need (with different species abbreviation):
  # hyci_vcfListFile.txt (created during the missingness step)
  # ALL_Samples_Sequenced_20250821.xlsx (download)
  # ALL_Samples_Sequenced_20250821.csv (you need to convert this from the xlsx version)
  # For each clustOpt threshold:
    # hyci_clustOpt##_missHM012.missingness (generated during the missingness step)
    # hyci_clustOpt##_s5_consensus.txt
    # hyci_clustOpt##.vcf
    # hyci_clustOpt##_stats.txt

# Make sure to clear environment for each analysis/metric since some variable names might be repetitive

# The Total Number of SNPs section assumes you have clustOpt runs 85-97. Change this if needed in the indicated line of code.

################################################################################
# MISSINGNESS
################################################################################

library(SNPRelate)
library(pheatmap)
library(grImport2)
library(dendextend)

# enter your species abbreviation !!!!
prefix <- "hyci"

#empty variables that will be filled as the code runs
clusters <- c()
PCCs <- c()

### Make sure that there is exactly one flag--the VCF file input list
args = paste0(prefix, "_vcfListFile.txt") #add vcf List File here
if (length(args)!=1) {
  stop("Exactly one argument must be supplied (<vcf file list>).n", call.=FALSE)
}

### Get and print a random seed for SNPRelate clustering
randomSeed <- sample(1:100000000, 1)
set.seed(randomSeed)
paste("Random seed for snpgdsHCluster = ", randomSeed, sep="")

### Core functions
makeMatrix <- function(vcfFilename, printStats = 0) {
  intermediateFileName <- paste(gsub(".vcf","", vcfFilename, fixed=T), '_missHM012.missingness', sep="")
  missInter <- read.csv(intermediateFileName, sep="\t", header=T)
  
  interMatrix <- matrix(nrow=length(unique(missInter$Sample1)), ncol=length(unique(missInter$Sample1)))
  interMatrix[lower.tri(interMatrix, diag=T)] <- (missInter$oneMiss / (missInter$oneMiss + missInter$neitherMiss))
  colnames(interMatrix) <- unique(missInter$Sample1)
  rownames(interMatrix) <- unique(missInter$Sample2)
  interMatrix[upper.tri(interMatrix, diag=F)] <- t(interMatrix)[upper.tri(interMatrix, diag=F)]
  
  if (printStats == 1) {
    minMissingness <- min(interMatrix, na.rm=T)
    maxMissingness <- max(interMatrix, na.rm=T)
    cat(paste("Minimum pairwise missingness for ", vcfFilename, ": ", minMissingness, sep=""), "\n")
    cat(paste("Maximum pairwise missingness for ", vcfFilename, ": ", maxMissingness, sep=""), "\n")
  }
  return(interMatrix)
}

makeVectorHeatmap <- function(missingnessMatrix, dendrogramz, name) {
  # Make the heat map SVG
  heatmapFile <- paste(name, ".heatmap.svg", sep = "")
  svg(heatmapFile, height=4, width=4.44)
  heatmapReordered <- missingnessMatrix[labels(dendrogramz), labels(dendrogramz)]
  pheatmap(heatmapReordered, breaks=mat_breaks, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=F, border_color = F)
  dev.off()
  
  # Make the dendrogram SVGs
  dendro1 <- paste(name, ".vDendro.svg", sep="")
  dendro2 <- paste(name, ".hDendro.svg", sep="")
  svg(dendro1, height=4.33, width=.6)
  par(mar=c(0,0,0,0), mgp=c(0,0,0))
  plot(rev(dendrogramz), type="rectangle", horiz = T, leaflab = "none", xaxt="n", yaxt="n")
  dev.off()
  
  svg(dendro2, height=.6, width=4.33)
  par(mar=c(0,0,0,0), mgp=c(0,0,0))
  plot(dendrogramz, type="rectangle", horiz = F, leaflab = "none", xaxt="n", yaxt="n")
  dev.off()
  
  dendroSVG <- readPicture(file = dendro1)
  dendroBsvg <- readPicture(file = dendro2)
  heatmapSVG <- readPicture(heatmapFile)
  outputname <- paste(name, ".pdf", sep="")
  pdf(outputname, width=10,height=10)
  plot(1, bty="n", type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 9), xaxt="n", yaxt="n")
  grid.picture(heatmapSVG, x = 5.4, y=4.85, hjust="centre", vjust="centre", width=6.66, height=6, default.units="in")
  grid.picture(dendroSVG, x=2.11, y=4.85, height=6.26, width=.75, default.units="in", distort=T)
  grid.picture(dendroBsvg, x=5.11, y=7.84, height=.75, width=6.25, default.units="in", distort=T)
  
  text(x=5.11, y=9, cex=1, labels = name)
  dev.off()
  
  
  
  # Delete intermediate files
  unlink(x = c(heatmapFile, dendro1, dendro2))
}


### Put the VCF filenames into a vector
vcfFiles <- scan(args[1], what = character())

### We need to define the quantile breaks across the collection of VCF files, so that
### they all use the same informative color scale for missingness. 
### Read in the first matrix,then rbind all subsequent matrices to the first matrix. 
### Then we'll get the breaks from superMatrix
firstFileName <- paste(gsub(".vcf","", vcfFiles[1], fixed=T), '_missHM012.missingness', sep="")
missFirst <- read.csv(firstFileName, sep="\t", header=T)

superMatrix <- matrix(nrow=length(unique(missFirst$Sample1)), ncol=length(unique(missFirst$Sample1)))
superMatrix[lower.tri(superMatrix, diag=T)] <- (missFirst$oneMiss / (missFirst$oneMiss + missFirst$neitherMiss))
colnames(superMatrix) <- unique(missFirst$Sample1)
rownames(superMatrix) <- unique(missFirst$Sample2)
superMatrix[upper.tri(superMatrix, diag=F)] <- t(superMatrix)[upper.tri(superMatrix, diag=F)]

if (length(vcfFiles) > 1) {
  for (i in 2:length(vcfFiles)) {
    matrixForBinding <- makeMatrix(vcfFiles[i], printStats = 0)
    
    # Make sure that the full and intermediate matrices line up
    if (!all.equal(colnames(matrixForBinding), colnames(superMatrix))) {
      stop("Col names don't match for ", args[1], " and ", args[i], ": Exiting now", call.=TRUE)
    }
    
    # Append the intermediate matrix to the full matrix
    superMatrix <- rbind(superMatrix, matrixForBinding)
  }
} else {
  # Here we have only a single VCF file, which is perfectly fine. We actually don't 
  # need to do anything in this case.
}

minMissingness <- min(superMatrix, na.rm=T)
maxMissingness <- max(superMatrix, na.rm=T)
paste("Minimum pairwise missingness across all VCFs (lower legend boundary): ", minMissingness, sep="")
paste("Maximum pairwise missingness across all VCFs (upper legend boundary): ", maxMissingness, sep="")


### Infer the quantile breaks from superMatrix:
quantile_breaks <- function(xs, n = 101) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = T)
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(superMatrix, n = 101)


for (i in 1:length(vcfFiles)) {
  # Just print out the minimum and maximum missingness values here
  missingMatrix <- makeMatrix(vcfFiles[i], printStats = 1)
}

### Now make the heat maps for each threshold 
for (i in 1:length(vcfFiles)) {
  # Make and load the GDS file:
  gdsOut <- paste(gsub(".vcf","", vcfFiles[i], fixed=T), '_missHM012.gds', sep="")
  
  snpgdsVCF2GDS(vcfFiles[i], gdsOut)
  
  gdsInter <- snpgdsOpen(gdsOut)
  # Compute the dendrogram:
  ibsInter <- snpgdsHCluster(snpgdsIBS(gdsInter, num.thread=2, autosome.only=FALSE)) # added the autosome.only
  
  # Make the matrix
  missingMatrix <- makeMatrix(vcfFiles[i], printStats = 0)
  heatmapFile = paste(gsub(".vcf","", vcfFiles[i], fixed=T), '.heatmap', sep="")
  
  # Print out the correlation coefficient of missingness as a function of relatedness:
  missGenCorrelation <- cor(c(missingMatrix[lower.tri(missingMatrix, diag=F)]), c(ibsInter$dist[lower.tri(ibsInter$dist, diag = F)]), method="pearson")
  cat("Correlation (PCC) between missingness and genetic similarity for ", vcfFiles[i], ": ", missGenCorrelation, "\n", sep="")
  
  # --- NEW: extract cluster number from filename ---
  # Example filename: "hyci_clustOpt85.vcf"
  cluster_num <- as.integer(sub(".*clustOpt(\\d+)\\.vcf", "\\1", vcfFiles[i]))
  
  # Store the values
  clusters <- c(clusters, cluster_num)
  PCCs <- c(PCCs, missGenCorrelation)
  
  # Plot the heatmaps
  #makeVectorHeatmap(missingMatrix, ibsInter$dendrogram, heatmapFile)
  #snpgdsClose(gdsInter)
  #file.remove(gdsOut)
}

# After finishing all files, write CSV
results_df <- data.frame(clust = clusters, PCC = PCCs)
write.csv(results_df, paste0(prefix, "_clustOpt_miss.csv"), row.names = FALSE)


# PLOT OUTPUTS -----------------------------------------------------------------

library(ggplot2)

miss_csv <- read.csv(paste0(prefix, "_clustOpt_miss.csv"))
miss_plot <- ggplot(miss_csv, aes(x = clust, y = PCC)) + geom_point(size=2) + labs(title="Pearson's correlation coefficient between pairwise \ngenetic dissimilarity and data missingness at \ndifferent clustering thresholds", x = "Clustering threshold (% similarity)", y = "PCC between genetic distance and missingness") + scale_x_continuous(labels = as.character(miss_csv$clust), breaks = miss_csv$clust)
miss_plot

ggsave(miss_plot, filename=paste0(prefix, "_miss.png"), bg="transparent", height=6, width=6, units="in")




################################################################################
# Slope of Genetic Divergence versus Geographic Distance
################################################################################

# CREATE CSV FILE OF SAMPLES AND GEOGRAPHIC COORDINATES ------------------------
library(readr)

# enter your species abbreviation !!!!
prefix <- "hyci"

# Set file paths
vcf_file <- paste0(prefix, "_clustOpt85.vcf") # this can be any clustOpt number you used, it does not matter
master_csv <- "ALL_Samples_Sequenced_20250821.csv"
output_csv <- paste0(prefix, "_clustOpt_IBD_latlong.csv")

# Extract samples from the vcf file
vcf_line11 <- readLines(vcf_file, n = 11)[11] # Read line 11 of the VCF file which contains sample names
vcf_fields <- strsplit(vcf_line11, "\t")[[1]] # Split the line by tabs
vcf_samples_trimmed <- vcf_fields[10:length(vcf_fields)] # Extract sample names (everything after the 9th column)
vcf_samples_clean <- sub("\\.trimmed$", "", vcf_samples_trimmed) # Remove ".trimmed" from sample names

# Read in excel sheet data
master_data <- read_csv(master_csv, col_select = c(UWBC_sampleName, decimalLatitude, decimalLongitude))

# Match VCF sample names (cleaned) to the master CSV, keeping VCF order
matched_data <- master_data[match(vcf_samples_clean, master_data$UWBC_sampleName), ]

# Add the VCF-ordered sample names as a new column to preserve order
matched_data$sample <- vcf_samples_trimmed  # overwrite just in case
colnames(matched_data) <- c("UWBC_sampleName", "lat", "long", "sample")

# Check for unmatched samples
if (any(is.na(matched_data$lat)) || any(is.na(matched_data$long))) {
  warning("Some samples from the VCF were not found in the master CSV!")
}

# Write the output CSV
write.csv(matched_data, output_csv, row.names = FALSE)


# FIND THE IBD SLOPES ----------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")

library(geosphere)
library(SNPRelate)

#empty variables that will be filled as the code runs
clusters <- c()
slopeList <- c()

### Make sure that there are exactly two flags–the VCF file input list [1] and the latlong file [2]
args = c(paste0(prefix, "_vcfListFile.txt"), paste0(prefix, "_clustOpt_IBD_latlong.csv"))
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied (<vcf file> and <latLongFile>).n", call.=FALSE)
}

# Pull in the VCF filenames
vcfFiles <- scan(args[1], what = character(), quiet=T)

# Set up the geographic distance matrix:
latLongFile <- read.csv(args[2], sep=",", header=T) #changes the sep
geoDistMatrix <- matrix(nrow=length(latLongFile$sample), ncol=length(latLongFile$sample))
rownames(geoDistMatrix) <- latLongFile$sample
colnames(geoDistMatrix) <- latLongFile$sample

# Calculate the pairwise geographic distance between each pair of points, in km:
for (i in 1:length(latLongFile$sample)) {
  for (j in 1:length(latLongFile$sample)) {
    geoDistMatrix[i,j] <- distm(c(latLongFile$long[i], latLongFile$lat[i]), c(latLongFile$long[j], latLongFile$lat[j]), fun = distHaversine) / 1000
  }
}

for (i in 1:length(vcfFiles)) {
  gdsOut <- paste(gsub(".vcf","", vcfFiles[i], fixed=T), '.pcaVar.gds', sep="")
  snpgdsVCF2GDS(vcfFiles[i], gdsOut, verbose=F)
  genofile <- snpgdsOpen(gdsOut)
  ibsMat <- snpgdsIBS(genofile, num.thread=2, verbose=F, autosome.only=FALSE) #added the autosome.only)
  
  #    cat(head(rownames(geoDistMatrix)))
  #    cat(head(ibsMat$sample.id))
  #    cat(match(rownames(geoDistMatrix), ibsMat$sample.id), "\n")
  
  # Rearrange IBS matrix to the same order as geoDistMatrix
  ibsMatRearranged <- ibsMat$ibs[match(rownames(geoDistMatrix), ibsMat$sample.id), match(rownames(geoDistMatrix), ibsMat$sample.id)]
  
  slopeToReport <- lm(100-ibsMatRearranged[lower.tri(ibsMatRearranged, diag=F)] ~ geoDistMatrix[lower.tri(geoDistMatrix, diag=F)])$coefficients[2]
  cat("IBD slope for ", vcfFiles[i], ": ", 100 * 100 * slopeToReport, "% increased SNP divergence per 100km\n", sep="")
  
  # --- NEW: extract cluster number from filename ---
  # Example filename: "hyci_clustOpt85.vcf"
  cluster_num <- as.integer(sub(".*clustOpt(\\d+)\\.vcf", "\\1", vcfFiles[i]))
  
  # Store the values
  clusters <- c(clusters, cluster_num)
  slopeList <- c(slopeList, (100 * 100 * slopeToReport))
  
  # -----------------------
  
  snpgdsClose(genofile)    
  file.remove(gdsOut)
}


# After finishing all files, write CSV
results_df <- data.frame(clust = clusters, pct_SNP_divergence_per_100km = slopeList)
write.csv(results_df, paste0(prefix, "_clustOpt_IBD.csv"), row.names = FALSE)

# PLOT OUTPUTS -----------------------------------------------------------------

library(ggplot2)

IBD_csv <- read.csv(paste0(prefix, "_clustOpt_IBD.csv"))
IBD_plot <- ggplot(IBD_csv, aes(x = clust, y = pct_SNP_divergence_per_100km)) + geom_point(size=2) + labs(title="Slopes of genetic isolation by geographic \ndistance across samples", x = "Clustering threshold (% similarity)", y = "% increased SNP divergence per 100 km")  + scale_x_continuous(labels = as.character(IBD_csv$clust), breaks = IBD_csv$clust)
IBD_plot

ggsave(IBD_plot, filename=paste0(prefix, "_IBD.png"), bg="transparent", height=6, width=6, units="in")



################################################################################
# Total Number of SNPs
################################################################################

# enter your species abbreviation !!!!
prefix <- "hyci"

# Define the range of cluster thresholds you have files for; CHANGE THIS IF NEEDED
clust_values <- 85:97

# Initialize a vector to store SNP counts
snp_counts <- numeric(length(clust_values))

# Loop through each threshold, read the corresponding stats file, and extract the SNP count
for (i in seq_along(clust_values)) {
  clust <- clust_values[i]
  
  # Construct the filename (change prefix if needed)
  filename <- paste0(prefix, "_clustOpt", clust, "_stats.txt")
  
  # Read all lines from the stats file
  lines <- readLines(filename)
  
  # Find the line containing "snps matrix size"
  snp_line <- lines[grepl("snps matrix size", lines)]
  
  # Extract the SNP count using a regex pattern
  # Example line: snps matrix size: (51, 312453), 65.92% missing sites.
  snp_count <- as.numeric(sub(".*\\(\\d+,\\s*(\\d+)\\).*", "\\1", snp_line))
  
  # Store the SNP count
  snp_counts[i] <- snp_count
}

# Create the dataframe
clustOpt_df <- data.frame(
  clust = clust_values,
  SNPs = snp_counts
)

# Write to CSV
write.csv(clustOpt_df, file = paste0(prefix, "_clustOpt_SNPs.csv"), row.names = FALSE)


# PLOT OUTPUTS -----------------------------------------------------------------
library(scales)
library(ggplot2)

SNPs_csv <- read.csv(paste0(prefix, "_clustOpt_SNPs.csv"))
SNPs_plot <- ggplot(SNPs_csv, aes(x = clust, y = SNPs)) + geom_point(size=2) + labs(title="Total SNPs recovered across different clustering \nthresholds", x = "Clustering threshold (% similarity)", y = "Total SNPs")  + scale_x_continuous(labels = as.character(SNPs_csv$clust), breaks = SNPs_csv$clust) + scale_y_continuous(labels = label_comma())
SNPs_plot

ggsave(SNPs_plot, filename=paste0(prefix, "_SNPs.png"), bg="transparent", height=6, width=6, units="in")


################################################################################
# Fraction of Variance Explained by the Main Principal Components of Genetic Variation
################################################################################

# Load SNPRelate
library(SNPRelate)

# enter your species abbreviation !!!!
prefix <- "hyci"

# Create empty variables
clusters <- c()
varList <- c()

### Make sure that there is exactly one flag–the VCF file input list
args = c(paste0(prefix, "_vcfListFile.txt"),"8") #  8 is the number of PCS
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied (<vcf file list> and <numPCsToSum>).n", call.=FALSE)
}

### Put the VCF filenames into a vector
vcfFiles <- scan(args[1], what = character(), quiet=T)

for (i in 1:length(vcfFiles)) {
  gdsOut <- paste(gsub(".vcf","", vcfFiles[i], fixed=T), '.pcaVar.gds', sep="") # makes gds file name
  snpgdsVCF2GDS(vcfFiles[i], gdsOut, verbose=F) #makes the gds file
  genofile <- snpgdsOpen(gdsOut) # opens gds file
  pcaz <- snpgdsPCA(genofile, num.thread=2, verbose=F, autosome.only=FALSE) #added the autosome.only
  varExplained <- sum(pcaz$varprop[1:args[2]])
  cat("Variance explained by first ", args[2], " PCs for ", vcfFiles[i], ": ", varExplained, "\n", sep="")
  
  # --- NEW: extract cluster number from filename ---
  # Example filename: "hyci_clustOpt85.vcf"
  cluster_num <- as.integer(sub(".*clustOpt(\\d+)\\.vcf", "\\1", vcfFiles[i]))
  
  # Store the values
  clusters <- c(clusters, cluster_num)
  varList <- c(varList, varExplained)
  
  # ------------
  snpgdsClose(genofile)
  file.remove(gdsOut)
}

# After finishing all files, write CSV
results_df <- data.frame(clust = clusters, variance = varList)
write.csv(results_df, paste0(prefix, "_clustOpt_PCA.csv"), row.names = FALSE)


# PLOT OUTPUTS -----------------------------------------------------------------

library(ggplot2)

PCA_csv <- read.csv(paste0(prefix, "_clustOpt_PCA.csv"))
PCA_plot <- ggplot(PCA_csv, aes(x = clust, y = variance)) + geom_point(size=2) + labs(title="Cumulative variance of all biallelic SNPs recovered \n by ipyrad explained by the first 8 principal \ncomponents across different clustering thresholds.", x = "Clustering threshold (% similarity)", y = "Cumulative variance in first n PCs") + scale_x_continuous(labels = as.character(PCA_csv$clust), breaks = PCA_csv$clust)
PCA_plot

ggsave(PCA_plot, filename=paste0(prefix, "_PCA.png"), bg="transparent", height=6, width=6, units="in")

################################################################################
# Heterozygosity
################################################################################

library(dplyr)

# enter your species abbreviation !!!!
prefix <- "hyci"

# Directory containing your files (change if needed)
file_dir <- "."  # current working directory

# List all files matching the pattern
files <- list.files(path = file_dir, pattern = paste0("^", prefix, "_clustOpt\\d+_s5_consensus\\.txt$"), full.names = TRUE)

# Initialize an empty data frame to store results
combined_data <- data.frame()

for (file in files) {
  # Extract clustering threshold from filename
  clust_val <- as.integer(sub(".*clustOpt(\\d+)_s5_consensus\\.txt", "\\1", basename(file)))
  
  # Read the file with whitespace separator, allow multiple spaces/tabs between columns
  df <- read.table(file, header = TRUE, strip.white = TRUE, row.names = NULL)
  
  # Extract sample name (1st col) and heterozygosity (last col)
  sample_col <- df[[1]]
  het_col <- df[[ncol(df)]]
  
  # Create new data frame with desired columns + clust
  temp_df <- data.frame(sample = sample_col,
                        heterozygosity = het_col,
                        clust = clust_val)
  
  # Append to combined data
  combined_data <- bind_rows(combined_data, temp_df)
}

# Write combined data to CSV
write.csv(combined_data, paste0(prefix, "_clustOpt_het.csv"), row.names = FALSE)


# PLOT OUTPUTS -----------------------------------------------------------------

library(ggplot2)

het_csv <- read.csv(paste0(prefix, "_clustOpt_het.csv"))
het_plot <- ggplot(het_csv, aes(x = clust, y = heterozygosity, group = clust)) + geom_boxplot() + scale_x_continuous(labels = as.character(het_csv$clust), breaks = het_csv$clust) + labs(title="Heterozygosity", x = "Clustering threshold (% similarity)", y = "Heterozygosity")
het_plot

ggsave(het_plot, filename=paste0(prefix, "_het.png"), bg="transparent", height=6, width=6, units="in")
