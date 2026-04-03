############################################################
# SRK population genetic summary - FIXED FOR TRUE COUNTS
# Works consistently for ALL populations
############################################################

cat("Starting SRK population genetic summary analysis\n")

############################
# 1. Load files
############################

geno <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

zyg <- read.table(
  "SRK_individual_zygosity.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

meta <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors = FALSE
)

############################
# 2. Keep ingroup only
############################

meta <- meta[meta$Ingroup == 1, ]

geno <- geno[geno$Individual %in% meta$SampleID, ]
zyg  <- zyg[zyg$Individual %in% meta$SampleID, ]

cat("Individuals retained:", nrow(geno), "\n")

############################
# 3. Add population info
############################

geno$Population <- meta$Pop[match(geno$Individual, meta$SampleID)]
zyg$Population  <- meta$Pop[match(zyg$Individual, meta$SampleID)]

############################
# 4. STANDARDIZED genotype processing
############################

# Use the SAME function as accumulation curves for consistency
standardized_genotype_processing <- function(data, exclude_columns = c("Individual", "Population")) {
  # Identify protein columns
  protein_cols <- setdiff(colnames(data), exclude_columns)
  
  cat("Total protein columns detected:", length(protein_cols), "\n")
  
  # Extract and convert to numeric matrix
  protein_data <- data[, protein_cols, drop = FALSE]
  
  # Robust conversion to numeric
  for (i in seq_along(protein_cols)) {
    col_data <- protein_data[, i]
    if (is.character(col_data) || is.factor(col_data)) {
      col_data <- as.character(col_data)
      col_data[is.na(col_data) | col_data == "" | col_data == "NA"] <- "0"
      
      numeric_col <- suppressWarnings(as.numeric(col_data))
      
      if (any(is.na(numeric_col) & col_data != "0")) {
        numeric_col <- ifelse(col_data == "0" | is.na(col_data), 0, 1)
      }
      
      protein_data[, i] <- numeric_col
    } else {
      protein_data[, i] <- ifelse(is.na(protein_data[, i]), 0, protein_data[, i])
    }
  }
  
  # Convert to matrix
  mat <- as.matrix(protein_data)
  mat <- apply(mat, 2, as.numeric)
  mat[is.na(mat)] <- 0
  
  return(list(matrix = mat, protein_cols = protein_cols))
}

# Process the genotype data using standardized method
processed <- standardized_genotype_processing(geno)
geno_matrix <- processed$matrix
allele_cols <- processed$protein_cols

# Add the processed matrix back to geno dataframe
geno[, allele_cols] <- geno_matrix

cat("Genotype matrix processing complete\n")

############################
# 5. Universal population summary function
############################

analyze_population <- function(pop) {

  cat("Analyzing population:", pop, "\n")
  
  g <- geno[geno$Population == pop, ]
  z <- zyg[zyg$Population == pop, ]

  # Extract matrix for this population
  mat <- as.matrix(g[, allele_cols, drop = FALSE])
  
  # Ensure numeric
  mat <- apply(mat, 2, as.numeric)
  mat[is.na(mat)] <- 0

  ########################
  # TRUE allele counts - works for any population
  ########################

  # Handle single individual vs multiple individuals
  if (nrow(mat) == 1) {
    allele_counts <- as.numeric(mat[1, ])
    names(allele_counts) <- colnames(mat)
  } else {
    allele_counts <- colSums(mat, na.rm = TRUE)
  }

  # Only count alleles actually present
  present <- allele_counts > 0
  N_alleles <- sum(present)
  
  cat("  Population", pop, ":", N_alleles, "alleles present in", nrow(mat), "individuals\n")

  ########################
  # Allele frequencies - only for present alleles
  ########################

  if (N_alleles > 0) {
    freqs <- allele_counts[present] / sum(allele_counts[present])
  } else {
    freqs <- numeric(0)
  }

  ########################
  # Effective allele number
  ########################

  if(length(freqs) > 0 && sum(freqs) > 0) {
    Ne <- 1 / sum(freqs^2)
  } else {
    Ne <- NA
  }

  ########################
  # Heterozygosity
  ########################

  if (nrow(z) > 0) {
    prop_het <- mean(z$Zygosity == "Heterozygous", na.rm = TRUE)
    mean_prot <- mean(z$N_distinct_alleles, na.rm = TRUE)
  } else {
    prop_het <- NA
    mean_prot <- NA
  }

  ########################
  # Additional metrics
  ########################
  
  # Count individuals with at least one allele
  if (nrow(mat) == 1) {
    individuals_with_alleles <- ifelse(sum(allele_counts) > 0, 1, 0)
  } else {
    individuals_with_alleles <- sum(rowSums(mat) > 0)
  }
  
  # Most common alleles
  if (N_alleles > 0) {
    top_alleles <- head(sort(allele_counts[present], decreasing = TRUE), 3)
    top_allele_info <- paste(names(top_alleles), "(", top_alleles, ")", collapse = ", ")
  } else {
    top_allele_info <- "None"
  }

  return(data.frame(
    Population = pop,
    N_individuals = nrow(g),
    N_individuals_with_alleles = individuals_with_alleles,
    N_alleles = N_alleles,
    Effective_alleles_Ne = round(Ne, 3),
    Prop_heterozygous = round(prop_het, 3),
    Prop_homozygous = round(1-prop_het, 3),
    Mean_alleles = round(mean_prot, 3),
    Top_alleles = top_allele_info,
    stringsAsFactors = FALSE
  ))
}

############################
# 6. Run analysis for all populations
############################

pops <- sort(unique(geno$Population[!is.na(geno$Population)]))

cat("\nAnalyzing", length(pops), "populations:", paste(pops, collapse = ", "), "\n\n")

# Process all populations with error handling
results <- data.frame()

for (pop in pops) {
  tryCatch({
    pop_result <- analyze_population(pop)
    results <- rbind(results, pop_result)
  }, error = function(e) {
    cat("ERROR analyzing population", pop, ":", e$message, "\n")
    # Create a placeholder row with NAs
    placeholder <- data.frame(
      Population = pop,
      N_individuals = sum(geno$Population == pop, na.rm = TRUE),
      N_individuals_with_alleles = NA,
      N_alleles = NA,
      Effective_alleles_Ne = NA,
      Prop_heterozygous = NA,
      Prop_homozygous = NA,
      Mean_alleles = NA,
      Top_alleles = "ERROR",
      stringsAsFactors = FALSE
    )
    results <<- rbind(results, placeholder)
  })
}

############################
# 7. Save results
############################

write.table(
  results,
  "SRK_population_genetic_summary.tsv",
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("\nSummary table saved to: SRK_population_genetic_summary.tsv\n")

print(results)

############################
# 8. Summary statistics across all populations
############################

valid_results <- results[!is.na(results$N_alleles), ]

if (nrow(valid_results) > 0) {
  cat("\n=== SUMMARY ACROSS ALL POPULATIONS ===\n")
  cat("Populations analyzed:", nrow(valid_results), "\n")
  cat("Total individuals:", sum(valid_results$N_individuals), "\n")
  cat("Allele richness range:", min(valid_results$N_alleles), "-", max(valid_results$N_alleles), "\n")
  cat("Mean alleles per population:", round(mean(valid_results$N_alleles), 1), "\n")
  cat("Populations with >20 alleles:", sum(valid_results$N_alleles > 20), "\n")
  cat("Populations with <5 alleles:", sum(valid_results$N_alleles < 5), "\n")
}

############################
# 9. Create plots for all populations
############################

pdf("SRK_population_genetic_summary.pdf", width=12, height=8)

par(mfrow=c(2,2))

if (nrow(valid_results) > 0) {

  ############ Effective alleles plot
  barplot(
    valid_results$Effective_alleles_Ne,
    names.arg=valid_results$Population,
    col="steelblue",
    ylab="Effective number of alleles (Ne)",
    main="SRK Effective Allele Diversity",
    las=2
  )

  ############ Heterozygosity plot
  barplot(
    valid_results$Prop_heterozygous,
    names.arg=valid_results$Population,
    col="tomato",
    ylab="Proportion heterozygous",
    main="SRK Heterozygosity",
    las=2
  )

  ############ Allele richness plot
  barplot(
    valid_results$N_alleles,
    names.arg=valid_results$Population,
    col="darkgreen",
    ylab="Number of alleles",
    main="SRK Allele Richness",
    las=2
  )

  # Add allele count labels on bars
  text(x = 1:nrow(valid_results), 
       y = valid_results$N_alleles + max(valid_results$N_alleles, na.rm=TRUE) * 0.02, 
       labels = valid_results$N_alleles, 
       cex = 0.8)

  ############ Sample sizes plot
  barplot(
    valid_results$N_individuals,
    names.arg=valid_results$Population,
    col="orange",
    ylab="Number of individuals",
    main="Sample Sizes",
    las=2
  )

} else {
  plot(1, 1, type="n", main="No valid data for plotting")
}

dev.off()

cat("\nPlots saved to: SRK_population_genetic_summary.pdf\n")
cat("\nAnalysis complete - showing true allele counts for all populations\n")

