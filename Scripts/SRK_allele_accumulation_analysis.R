############################################################
# SRK allele accumulation analysis - FIXED TO MATCH CORRECT METHOD
# Uses same approach as your working R code
############################################################

cat("Starting SRK allele accumulation analysis\n")

############################################################
# 1. Load data
############################################################

geno <- read.table(
  "SRK_individual_genotypes.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_info <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors = FALSE
)

############################################################
# 2. Merge and filter ingroup - FIXED individual matching
############################################################

# Filter metadata to ingroup only
sample_info <- sample_info[sample_info$Ingroup == 1, ]

# Filter genotype data to match ingroup individuals
# FIXED: Use exact matching like your working code
geno <- geno[geno$Individual %in% sample_info$SampleID, ]

cat("Individuals after ingroup filtering:", nrow(geno), "\n")

# Add population info using correct matching
geno$Population <- sample_info$Pop[match(geno$Individual, sample_info$SampleID)]

############################################################
# 3. CORRECT allele processing - columns 2 onwards
############################################################

# FIXED: Alleles are columns 2 onwards (like your working code)
allele_start_col <- 2
allele_end_col <- ncol(geno) - 1  # Exclude the Population column we just added

cat("Allele columns: from", allele_start_col, "to", allele_end_col, "\n")
cat("Total allele columns:", allele_end_col - allele_start_col + 1, "\n")

# Extract allele matrix - CORRECT method
allele_matrix <- as.matrix(geno[, allele_start_col:allele_end_col])

# Convert to numeric and handle NAs
allele_matrix <- apply(allele_matrix, 2, as.numeric)
allele_matrix[is.na(allele_matrix)] <- 0

# Add row names
rownames(allele_matrix) <- geno$Individual

cat("Genotype matrix dimensions:", nrow(allele_matrix), "x", ncol(allele_matrix), "\n")

# Safety check
cat("Genotype value summary:\n")
print(summary(as.vector(allele_matrix)))

############################################################
# 4. CORRECT allele accumulation function
############################################################

allele_accumulation <- function(mat, pop_name = "Unknown", n_perm = 1000, seed = 123) {

  set.seed(seed)

  n_ind <- nrow(mat)
  
  # CORRECT allele count - same method as your working code
  allele_counts <- colSums(mat)
  present_alleles <- allele_counts > 0
  true_alleles <- sum(present_alleles)
  
  cat("Population", pop_name, "- Individuals:", n_ind, "- Alleles present:", true_alleles, "\n")

  accum <- matrix(0, nrow = n_perm, ncol = n_ind)

  for (p in 1:n_perm) {

    order <- sample(n_ind)
    cum_mat <- mat[order, , drop = FALSE]

    for (i in 1:n_ind) {
      # Count alleles present in first i individuals
      if (i == 1) {
        cum_subset <- cum_mat[1, , drop = FALSE]
      } else {
        cum_subset <- cum_mat[1:i, , drop = FALSE]
      }
      
      # Count alleles with at least one occurrence
      accum[p, i] <- sum(colSums(cum_subset) > 0)
    }
  }

  mean_accum <- colMeans(accum)
  sd_accum <- apply(accum, 2, sd)

  # NFDS statistics
  initial_range <- min(5, n_ind)
  tail_range <- min(5, n_ind)
  
  initial_slope <- if(n_ind > 1) {
    mean(diff(mean_accum[1:initial_range]))
  } else {
    NA
  }

  tail_slope <- if(n_ind > tail_range) {
    mean(diff(mean_accum[(n_ind-tail_range+1):n_ind]))
  } else {
    NA
  }

  half <- floor(n_ind / 2)
  
  late_fraction <- if(n_ind > 2 && max(mean_accum) > 0) {
    (max(mean_accum) - mean_accum[half]) / max(mean_accum)
  } else {
    NA
  }

  return(list(
    mean_accum = mean_accum,
    sd_accum = sd_accum,
    initial_slope = initial_slope,
    tail_slope = tail_slope,
    late_fraction = late_fraction,
    true_alleles = true_alleles
  ))
}

############################################################
# 5. Open PDF output
############################################################

pdf(
  "SRK_allele_accumulation_curves.pdf",
  width = 8,
  height = 6
)

############################################################
# 6. Species-level analysis
############################################################

cat("\nRunning species-level analysis\n")

species_res <- allele_accumulation(allele_matrix, "Species")

n_ind <- nrow(allele_matrix)

plot(
  1:n_ind,
  species_res$mean_accum,
  type = "l",
  lwd = 3,
  col = "blue",
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of functional alleles",
  main = paste("Species-level SRK accumulation\n(", species_res$true_alleles, "total alleles)")
)

polygon(
  c(1:n_ind, rev(1:n_ind)),
  c(species_res$mean_accum - species_res$sd_accum,
    rev(species_res$mean_accum + species_res$sd_accum)),
  col = rgb(0,0,1,0.2),
  border = NA
)

lines(1:n_ind, species_res$mean_accum, lwd = 3)

############################################################
# 7. Population-level analysis - CORRECT method
############################################################

pop_counts <- table(geno$Population)
valid_pops <- names(pop_counts[pop_counts >= 5])

cat("\nPopulations retained:", valid_pops, "\n")

pop_results <- list()

for (pop in valid_pops) {

  cat("Processing population:", pop, "\n")

  # CORRECT population subsetting - same as your working code
  target_pop_metadata <- sample_info[sample_info$Pop == pop, ]
  pop_geno <- geno[geno$Individual %in% target_pop_metadata$SampleID, ]
  
  # Extract allele matrix for this population - CORRECT method
  pop_allele_matrix <- as.matrix(pop_geno[, allele_start_col:allele_end_col])
  pop_allele_matrix <- apply(pop_allele_matrix, 2, as.numeric)
  pop_allele_matrix[is.na(pop_allele_matrix)] <- 0
  rownames(pop_allele_matrix) <- pop_geno$Individual

  res <- allele_accumulation(pop_allele_matrix, pop)

  pop_results[[pop]] <- res

  n_ind <- nrow(pop_allele_matrix)

  plot(
    1:n_ind,
    res$mean_accum,
    type = "l",
    lwd = 3,
    col = "red",
    xlab = "Number of individuals sampled",
    ylab = "Cumulative number of functional alleles",
    main = paste("Population", pop, "\n(", res$true_alleles, "alleles,", n_ind, "individuals)")
  )

  polygon(
    c(1:n_ind, rev(1:n_ind)),
    c(res$mean_accum - res$sd_accum,
      rev(res$mean_accum + res$sd_accum)),
    col = rgb(1,0,0,0.2),
    border = NA
  )

  lines(1:n_ind, res$mean_accum, lwd = 3)

}

############################################################
# 8. Close PDF
############################################################

dev.off()

cat("\nPDF written: SRK_allele_accumulation_curves.pdf\n")

############################################################
# 9. Export statistics with CORRECT counts
############################################################

stats <- data.frame(
  Level="Species",
  Population="All",
  N_individuals=nrow(allele_matrix),
  N_alleles=species_res$true_alleles,
  Initial_slope=species_res$initial_slope,
  Tail_slope=species_res$tail_slope,
  Late_fraction=species_res$late_fraction
)

for (pop in names(pop_results)) {

  res <- pop_results[[pop]]
  pop_data <- geno[geno$Population == pop, ]

  stats <- rbind(stats,
    data.frame(
      Level="Population",
      Population=pop,
      N_individuals=nrow(pop_data),
      N_alleles=res$true_alleles,
      Initial_slope=res$initial_slope,
      Tail_slope=res$tail_slope,
      Late_fraction=res$late_fraction
    )
  )
}

write.table(
  stats,
  "SRK_allele_accumulation_stats.tsv",
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("\nStats written: SRK_allele_accumulation_stats.tsv\n")

############################################################
# 10. Print summary with validation
############################################################

print(stats)

# Validation - should match your working code
cat("\n=== VALIDATION ===\n")
cat("Population 27 should now show 17 alleles\n")
cat("Method matches your working R code:\n")
cat("- Alleles from column 2 onwards\n")
cat("- Correct individual matching\n")
cat("- colSums(matrix) > 0 counting method\n")

cat("\nAnalysis complete\n")

