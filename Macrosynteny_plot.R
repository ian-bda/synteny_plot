library(circlize)
library(dplyr)
library(readr)

# Set the directory containing your CSV files
setwd('xxxx')

# List all CSV files
file_list <- list.files(pattern = "*.csv")

# Define expected chromosomes for each column
expected_chrom1 <- c("20", "20_2", "16", "17", "14")
expected_chrom2 <- c("16", "19", "14", "20_2", "20")

# Read and standardize column types while preserving original chromosome names
synteny_data <- lapply(file_list, function(file) {
  read_csv(file, col_types = cols(
    gene1 = col_character(),
    chrom1 = col_character(),
    start1 = col_double(),
    start2 = col_double(),
    gene2 = col_character(),
    chrom2 = col_character(),
    .default = col_character()
  )) %>%
    mutate(
      chrom1 = ifelse(chrom1 == "20-2", "20_2", chrom1),
      chrom2 = ifelse(chrom2 == "20-2", "20_2", chrom2)
    )
}) %>%
  bind_rows()

# Check and report any unexpected or missing chromosomes
unexpected_in_chrom1 <- setdiff(unique(synteny_data$chrom1), expected_chrom1)
unexpected_in_chrom2 <- setdiff(unique(synteny_data$chrom2), expected_chrom2)
missing_in_chrom1 <- setdiff(expected_chrom1, unique(synteny_data$chrom1))
missing_in_chrom2 <- setdiff(expected_chrom2, unique(synteny_data$chrom2))

if(length(unexpected_in_chrom1) > 0) {
  warning("Unexpected chromosomes in chrom1: ", paste(unexpected_in_chrom1, collapse = ", "))
}
if(length(unexpected_in_chrom2) > 0) {
  warning("Unexpected chromosomes in chrom2: ", paste(unexpected_in_chrom2, collapse = ", "))
}
if(length(missing_in_chrom1) > 0) {
  warning("Missing chromosomes in chrom1: ", paste(missing_in_chrom1, collapse = ", "))
}
if(length(missing_in_chrom2) > 0) {
  warning("Missing chromosomes in chrom2: ", paste(missing_in_chrom2, collapse = ", "))
}

# Display unique chromosome names to verify
print("Unique values in chrom1:")
print(sort(unique(synteny_data$chrom1)))
print("Unique values in chrom2:")
print(sort(unique(synteny_data$chrom2)))

# plot ------------------------------------------------------------------

# Convert chromosome numbers to factors and scale the positions
synteny_data$chrom1 <- factor(synteny_data$chrom1)
synteny_data$chrom2 <- factor(synteny_data$chrom2)
synteny_data$start1 <- synteny_data$start1 / 1e5
synteny_data$start2 <- synteny_data$start2 / 1e5

# Get unique chromosomes and their ranges
chroms <- unique(c(as.character(synteny_data$chrom1), as.character(synteny_data$chrom2)))
chrom_ranges <- data.frame(
  chrom = chroms,
  start = 0,
  end = sapply(chroms, function(chr) {
    max(c(synteny_data$start1[synteny_data$chrom1 == chr],
          synteny_data$start2[synteny_data$chrom2 == chr]))
  })
)

# Calculate the average end value of other chromosomes
avg_end <- mean(chrom_ranges$end[chrom_ranges$chrom != "17"])

# Set a minimum size for chromosome 17 (e.g., half of the average)
min_size_17 <- avg_end / 2

# Adjust the end value for chromosome 17 if it's below the minimum size
chrom_ranges$end[chrom_ranges$chrom == "17"] <- max(chrom_ranges$end[chrom_ranges$chrom == "17"], min_size_17)

# Normalize chromosome lengths to the range 0-100
chrom_ranges$length <- chrom_ranges$end - chrom_ranges$start
total_length <- sum(chrom_ranges$length)
chrom_ranges$prop <- chrom_ranges$length / total_length
chrom_ranges$end_norm <- cumsum(chrom_ranges$prop) * 100
chrom_ranges$start_norm <- c(0, chrom_ranges$end_norm[-nrow(chrom_ranges)])

# Create the xlim matrix for circos.initialize
xlim_matrix <- cbind(chrom_ranges$start_norm, chrom_ranges$end_norm)
rownames(xlim_matrix) <- chrom_ranges$chrom


# Initialize the circular layout
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = 5)

# Modified circos.initialize to use normalized chromosome lengths
circos.initialize(factors = chrom_ranges$chrom, xlim = xlim_matrix)  # Initialize with 0-100 range

# Add base ideogram track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05)

# Add chromosome labels outside the sections
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(10),  # Increased distance using mm_y(10)
              CELL_META$sector.index, cex = 0.8, niceFacing = TRUE) # Adjusted text position
}, bg.border = NA)

# Add links between syntenic regions, force color to black
for(i in 1:nrow(synteny_data)) {
  # Calculate normalized positions for links
  
  chr1 <- as.character(synteny_data$chrom1[i])
  chr2 <- as.character(synteny_data$chrom2[i])
  
  start1_norm <- (synteny_data$start1[i] / max(synteny_data$start1[synteny_data$chrom1 == chr1])) * (xlim_matrix[chr1,2] - xlim_matrix[chr1,1]) + xlim_matrix[chr1,1]
  start2_norm <- (synteny_data$start2[i] / max(synteny_data$start2[synteny_data$chrom2 == chr2])) * (xlim_matrix[chr2,2] - xlim_matrix[chr2,1]) + xlim_matrix[chr2,1]
  
  circos.link(chr1, start1_norm,
              chr2, start2_norm,
              col = rgb(0.5, 0.5, 0.5, 0.1))  # Light grey with high transparency
  
}

# Add ticks and labels - Corrected
for (chr in chroms) {
  # Calculate tick positions based on the normalized range of each chromosome
  start_pos <- xlim_matrix[chr, 1]
  end_pos <- xlim_matrix[chr, 2]
  tick_positions <- seq(start_pos, end_pos, length.out = 5)  # 5 ticks: 0%, 25%, 50%, 75%, 100%
  
  # Generate labels (0, 25, 50, 75, 100)
  labels <- c("0", "25", "50", "75", "100")
  
  circos.axis(sector.index = chr, 
              major.at = tick_positions, 
              labels = labels, 
              labels.cex = 0.6,
              labels.facing = "outside")
}
