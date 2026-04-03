############################################################
# SRK zygosity analysis from allele count matrix
############################################################

cat("\nStarting SRK zygosity analysis from allele count matrix\n")

############################################################
# Load allele count matrix
############################################################

geno <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
# Identify allele columns and ensure numeric
############################################################

allele_cols <- setdiff(colnames(geno), "Individual")

geno[, allele_cols] <- lapply(geno[, allele_cols], as.numeric)

############################################################
# Compute per-individual statistics
############################################################

# N_distinct_alleles: number of allele bins with count > 0
N_distinct_alleles <- rowSums(geno[, allele_cols] > 0)

# N_total_proteins: sum of all allele copy counts per individual
N_total_proteins <- rowSums(geno[, allele_cols])

############################################################
# Assign zygosity
############################################################

Zygosity <- ifelse(N_distinct_alleles == 1, "Homozygous", "Heterozygous")

############################################################
# Assign genotype model
# AAAB: 2 distinct alleles, dominant allele has >= 2x copies of the minor allele
# AABB: 2 distinct alleles, roughly equal copy numbers
############################################################

Genotype <- apply(geno[, allele_cols, drop = FALSE], 1, function(row) {
  counts_present <- sort(row[row > 0], decreasing = TRUE)
  n <- length(counts_present)
  if (n == 0) return("Unknown")
  if (n == 1) return("AAAA")
  if (n == 2) {
    if (counts_present[1] >= 2 * counts_present[2]) return("AAAB")
    return("AABB")
  }
  if (n == 3) return("AABC")
  if (n >= 4) return("ABCD")
  return("Unknown")
})

############################################################
# Build allele composition string
# e.g. "Allele_044(2)+Allele_047(2)"
############################################################

Allele_composition <- apply(geno[, allele_cols, drop = FALSE], 1, function(row) {
  present <- sort(row[row > 0])   # sort by allele name (alphabetical)
  if (length(present) == 0) return("None")
  paste(paste0(names(present), "(", present, ")"), collapse = "+")
})

############################################################
# Assemble output
############################################################

zyg <- data.frame(
  Individual         = geno$Individual,
  N_distinct_alleles = N_distinct_alleles,
  N_total_proteins   = N_total_proteins,
  Zygosity           = Zygosity,
  Genotype           = Genotype,
  Allele_composition = Allele_composition,
  stringsAsFactors = FALSE
)

############################################################
# Save results
############################################################

write.table(
  zyg,
  "SRK_individual_zygosity.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Summary
############################################################

cat("\nNumber of individuals:", nrow(zyg), "\n")

cat("\nDistinct allele count distribution:\n")
print(table(zyg$N_distinct_alleles))

cat("\nZygosity distribution:\n")
print(table(zyg$Zygosity))

cat("\nGenotype distribution:\n")
print(table(zyg$Genotype))

############################################################
# Plot
############################################################

pdf("SRK_zygosity_distribution.pdf")

barplot(
  table(zyg$N_distinct_alleles),
  main = "SRK genotype diversity (tetraploid)",
  xlab = "Number of distinct SRK alleles",
  ylab = "Number of individuals",
  col = "gray70"
)

dev.off()

cat("\nOutput files created:\n")
cat("SRK_individual_zygosity.tsv\n")
cat("SRK_zygosity_distribution.pdf\n\n")

cat("Analysis complete\n")
