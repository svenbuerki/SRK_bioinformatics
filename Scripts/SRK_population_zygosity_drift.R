############################################################
# SRK population-level zygosity and drift analysis
# FIXED v2 — numeric conversion corrected
############################################################

cat("\nStarting SRK population zygosity analysis (FIXED v2)\n")

############################################################
# Load files
############################################################

zyg <- read.table(
  "SRK_individual_zygosity.tsv",
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

meta <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors=FALSE
)

geno <- read.table(
  "SRK_individual_genotypes.tsv",
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

############################################################
# Merge metadata
############################################################

zyg <- merge(zyg, meta, by.x="Individual", by.y="SampleID")

geno <- merge(geno, meta, by.x="Individual", by.y="SampleID")

zyg <- zyg[zyg$Ingroup==1,]
geno <- geno[geno$Ingroup==1,]

cat("Individuals retained:", nrow(zyg), "\n")

############################################################
# Detect protein columns
############################################################

exclude_cols <- c(
  "Individual",
  "Pop",
  "Population",
  "Ingroup"
)

protein_cols <- setdiff(colnames(geno), exclude_cols)

cat("Protein columns detected:", length(protein_cols), "\n")

if(length(protein_cols)==0)
stop("No protein columns detected")

############################################################
# FORCE numeric conversion
############################################################

geno[, protein_cols] <- lapply(

  geno[, protein_cols],

  function(x) as.numeric(as.character(x))

)

############################################################
# Population analysis
############################################################

populations <- unique(zyg$Pop)

results <- data.frame()

for(pop in populations){

  sub_z <- zyg[zyg$Pop==pop,]

  sub_g <- geno[geno$Pop==pop,]

  if(nrow(sub_z)<5) next

  ##################################################

  N_ind <- nrow(sub_z)

  prop_het <- sum(sub_z$Zygosity=="Heterozygous") / N_ind

  prop_hom <- sum(sub_z$Zygosity=="Homozygous") / N_ind

  mean_prot <- mean(sub_z$N_proteins)

  ##################################################
  # Numeric matrix
  ##################################################

  mat <- as.matrix(sub_g[, protein_cols])

  ##################################################
  # Allele counts
  ##################################################

  allele_counts <- colSums(mat, na.rm=TRUE)

  N_alleles <- sum(allele_counts>0)

  ##################################################
  # Effective number of alleles
  ##################################################

  total_counts <- sum(allele_counts)

  if(total_counts>0){

    freq <- allele_counts / total_counts

    Ne <- 1 / sum(freq^2)

  } else {

    Ne <- NA

  }

  ##################################################

  results <- rbind(results, data.frame(

    Population=pop,

    N_individuals=N_ind,

    Prop_heterozygous=round(prop_het,3),

    Prop_homozygous=round(prop_hom,3),

    Mean_proteins=round(mean_prot,2),

    N_alleles=N_alleles,

    Effective_alleles_Ne=round(Ne,2)

  ))
}

############################################################
# Save
############################################################

write.table(

results,

"SRK_population_zygosity_drift.tsv",

sep="\t",

quote=FALSE,

row.names=FALSE

)

############################################################
# Plot
############################################################

pdf(

"SRK_population_heterozygosity.pdf",

width=8,

height=6

)

barplot(

results$Prop_heterozygous,

names.arg=results$Population,

col="grey",

ylab="Heterozygosity",

main="SRK heterozygosity by population"

)

dev.off()

############################################################

cat("\nSUCCESS\n")
cat("Output: SRK_population_zygosity_drift.tsv\n")
cat("Plot: SRK_population_heterozygosity.pdf\n\n")
