############################################
# SRK cross compatibility analysis
# Tetraploid-compatible classification
############################################

cat("Starting SRK cross compatibility analysis\n")

## ----------------------------
## Load data
## ----------------------------

meta <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors = FALSE
)

crosses <- read.csv(
  "Crosses_w_seeds_LEPA.csv",
  stringsAsFactors = FALSE
)

geno <- read.table(
  "SRK_individual_genotypes.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

zyg <- read.table(
  "SRK_individual_zygosity.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

## ----------------------------
## Map OccurrenceID → Individual
## ----------------------------

occ_map <- meta[, c("SampleID", "OccurrenceID")]

crosses$Mother <- occ_map$SampleID[
  match(
    crosses$pollenRecieverOccurrenceID,
    occ_map$OccurrenceID
  )
]

crosses$Father <- occ_map$SampleID[
  match(
    crosses$pollenDonorOccurrenceID,
    occ_map$OccurrenceID
  )
]

cat("Crosses total:", nrow(crosses), "\n")

## ----------------------------
## Prepare genotype matrix
## ----------------------------

rownames(geno) <- geno$Individual
geno$Individual <- NULL

geno[] <- lapply(geno, as.numeric)

geno_mat <- as.matrix(geno)

## ----------------------------
## Prepare zygosity lookup
## ----------------------------

rownames(zyg) <- zyg$Individual

## ----------------------------
## Function to count shared
## ----------------------------

count_shared <- function(id1, id2) {

  if (is.na(id1) | is.na(id2)) return(NA)

  if (!(id1 %in% rownames(geno_mat))) return(NA)
  if (!(id2 %in% rownames(geno_mat))) return(NA)

  g1 <- geno_mat[id1, ]
  g2 <- geno_mat[id2, ]

  sum(g1 == 1 & g2 == 1)

}

## ----------------------------
## Storage
## ----------------------------

n <- nrow(crosses)

Shared_alleles <- rep(NA, n)

Mother_zygosity <- rep(NA, n)
Father_zygosity <- rep(NA, n)

Mother_genotype <- rep(NA, n)
Father_genotype <- rep(NA, n)

## ----------------------------
## Main loop
## ----------------------------

for (i in seq_len(n)) {

  mom <- crosses$Mother[i]
  dad <- crosses$Father[i]

  Shared_alleles[i] <- count_shared(mom, dad)

  if (!is.na(mom) && mom %in% rownames(zyg)) {

    Mother_zygosity[i] <- zyg[mom, "Zygosity"]
    Mother_genotype[i] <- zyg[mom, "Genotype"]

  }

  if (!is.na(dad) && dad %in% rownames(zyg)) {

    Father_zygosity[i] <- zyg[dad, "Zygosity"]
    Father_genotype[i] <- zyg[dad, "Genotype"]

  }

}

## ----------------------------
## Attach info
## ----------------------------

crosses$Shared_alleles <- Shared_alleles

crosses$Mother_zygosity <- Mother_zygosity
crosses$Father_zygosity <- Father_zygosity

crosses$Mother_genotype <- Mother_genotype
crosses$Father_genotype <- Father_genotype

crosses$Cross_type <- paste(
  Mother_zygosity,
  "x",
  Father_zygosity
)

## ----------------------------
## Compatibility classification
## ----------------------------

crosses$Compatibility_class <- NA

crosses$Compatibility_class[
  crosses$Shared_alleles == 0
] <- "Fully_compatible"

crosses$Compatibility_class[
  crosses$Shared_alleles == 1
] <- "Partially_compatible"

crosses$Compatibility_class[
  crosses$Shared_alleles >= 2
] <- "Mostly_incompatible"

crosses$Compatibility_class[
  is.na(crosses$Shared_alleles)
] <- "Unknown"

## ----------------------------
## Seeds
## ----------------------------

crosses$Seeds <- suppressWarnings(
  as.numeric(crosses$germplasmQuantityEstimate)
)

crosses$Seeds[is.na(crosses$Seeds)] <- 0

crosses$Success <- ifelse(
  crosses$Seeds > 0,
  1,
  0
)

## ----------------------------
## Final table
## ----------------------------

final <- data.frame(

  CrossID = seq_len(n),

  Mother = crosses$Mother,
  Father = crosses$Father,

  Mother_zygosity = crosses$Mother_zygosity,
  Father_zygosity = crosses$Father_zygosity,

  Mother_genotype = crosses$Mother_genotype,
  Father_genotype = crosses$Father_genotype,

  Cross_type = crosses$Cross_type,

  Shared_alleles = crosses$Shared_alleles,

  Compatibility_class = crosses$Compatibility_class,

  Seeds = crosses$Seeds,

  Success = crosses$Success,

  stringsAsFactors = FALSE

)

## ----------------------------
## Save
## ----------------------------

write.table(

  final,

  "SRK_cross_compatibility.tsv",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

cat("\nAnalysis complete\n")

table(final$Compatibility_class)
