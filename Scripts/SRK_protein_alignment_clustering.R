############################################################
# SRK Dominance Class Prediction - INGROUP FILTERED VERSION
# Filters proteins to only include ingroup individuals first
############################################################

cat("Starting protein alignment-based SRK class analysis (ingroup only)\n")

############################################################
# 0. Set CRAN mirror and handle package installation
############################################################

cat("=== SETTING UP R ENVIRONMENT ===\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))
cat("✓ CRAN mirror set\n")

# Enhanced package installation function
install_if_missing_safe <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      
      tryCatch({
        install.packages(pkg, dependencies = TRUE,
                        repos = "https://cran.rstudio.com/",
                        type = "both")
        
        if (require(pkg, character.only = TRUE, quietly = TRUE)) {
          cat("✓ Successfully installed and loaded:", pkg, "\n")
        } else {
          cat("⚠ Installed but failed to load:", pkg, "\n")
        }
      }, error = function(e) {
        cat("✗ Failed to install", pkg, "\n")
      })
    } else {
      cat("✓", pkg, "already available\n")
    }
  }
}

# Required packages (minimal set)
required_packages <- c("cluster", "seqinr")

# Install/load packages
install_if_missing_safe(required_packages)

# Check if essential packages loaded
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop("Essential package ", pkg, " not available. Please install manually.")
  }
}

cat("✓ All required packages loaded successfully\n")

############################################################
# 1. Load metadata and filter to ingroup
############################################################

cat("\n=== FILTERING TO INGROUP INDIVIDUALS ===\n")

# Load sampling metadata
if (!file.exists("sampling_metadata.csv")) {
  stop("Metadata file not found: sampling_metadata.csv")
}

metadata <- read.csv("sampling_metadata.csv", stringsAsFactors = FALSE)
cat("✓ Metadata loaded with", nrow(metadata), "samples\n")

# Check for required columns
if (!"Ingroup" %in% colnames(metadata)) {
  stop("Column 'Ingroup' not found in sampling_metadata.csv")
}

if (!"SampleID" %in% colnames(metadata)) {
  stop("Column 'SampleID' not found in sampling_metadata.csv")
}

# Filter to ingroup only
ingroup_metadata <- metadata[metadata$Ingroup == 1, ]
cat("✓ Filtered to", nrow(ingroup_metadata), "ingroup individuals\n")

if (nrow(ingroup_metadata) == 0) {
  stop("No ingroup individuals found (Ingroup == 1)")
}

# Get list of ingroup sample IDs
ingroup_samples <- ingroup_metadata$SampleID
cat("Ingroup samples:", length(ingroup_samples), "\n")

############################################################
# 2. Load protein key and filter to ingroup
############################################################

cat("\n=== FILTERING PROTEIN KEY TO INGROUP ===\n")

if (!file.exists("SRK_functional_protein_key.tsv")) {
  stop("Protein key file not found: SRK_functional_protein_key.tsv")
}

protein_key <- read.table("SRK_functional_protein_key.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("✓ Protein key loaded with", nrow(protein_key), "entries\n")

# Extract individual IDs from sequence names
protein_key$Individual <- sapply(strsplit(protein_key$Original_sequence_ID, "_"), 
                                function(x) paste(x[1:2], collapse="_"))

# Filter protein key to ingroup individuals only
ingroup_protein_key <- protein_key[protein_key$Individual %in% ingroup_samples, ]
cat("✓ Filtered protein key to", nrow(ingroup_protein_key), "entries from ingroup\n")

if (nrow(ingroup_protein_key) == 0) {
  stop("No protein sequences found for ingroup individuals")
}

# Get list of ingroup proteins
ingroup_proteins <- unique(ingroup_protein_key$Protein)
cat("✓ Ingroup proteins identified:", length(ingroup_proteins), "\n")

############################################################
# 3. Load and filter protein alignment to ingroup
############################################################

cat("\n=== LOADING AND FILTERING PROTEIN ALIGNMENT ===\n")

# Check if alignment file exists
if (!file.exists("SRK_proteins_aligned_accurate.fasta")) {
  stop("Alignment file not found: SRK_proteins_aligned_accurate.fasta")
}

# Read full alignment
alignment_full <- read.alignment("SRK_proteins_aligned_accurate.fasta", format = "fasta")
cat("✓ Full alignment loaded with", alignment_full$nb, "sequences\n")

# Filter alignment to ingroup proteins only
ingroup_indices <- which(alignment_full$nam %in% ingroup_proteins)

if (length(ingroup_indices) == 0) {
  stop("No ingroup proteins found in alignment file")
}

# Create filtered alignment
alignment <- list(
  nb = length(ingroup_indices),
  nam = alignment_full$nam[ingroup_indices],
  seq = alignment_full$seq[ingroup_indices]
)

cat("✓ Filtered alignment to", alignment$nb, "ingroup sequences\n")
cat("✓ Alignment length:", nchar(alignment$seq[[1]]), "positions\n")

# Convert to matrix format for analysis
alignment_matrix <- matrix(NA, nrow = alignment$nb, ncol = nchar(alignment$seq[[1]]))
rownames(alignment_matrix) <- alignment$nam
colnames(alignment_matrix) <- 1:ncol(alignment_matrix)

# Fill matrix with amino acids
for (i in 1:alignment$nb) {
  alignment_matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))
}

cat("✓ Alignment converted to matrix format\n")
cat("Matrix dimensions:", nrow(alignment_matrix), "x", ncol(alignment_matrix), "\n")

# Remove gap-only columns
gap_cols <- apply(alignment_matrix, 2, function(x) all(x == "-"))
alignment_matrix_clean <- alignment_matrix[, !gap_cols]
cat("✓ Removed", sum(gap_cols), "gap-only columns\n")
cat("Clean matrix dimensions:", nrow(alignment_matrix_clean), "x", ncol(alignment_matrix_clean), "\n")

############################################################
# 4. Calculate protein distance matrices
############################################################

cat("\n=== CALCULATING PROTEIN DISTANCES ===\n")

# Function to calculate amino acid distances
calculate_protein_distances <- function(alignment_matrix) {
  n_seqs <- nrow(alignment_matrix)
  seq_names <- rownames(alignment_matrix)
  
  # Initialize distance matrix
  dist_matrix <- matrix(0, nrow = n_seqs, ncol = n_seqs)
  rownames(dist_matrix) <- colnames(dist_matrix) <- seq_names
  
  cat("Calculating pairwise distances for", n_seqs, "sequences...\n")
  
  # Calculate pairwise distances
  for (i in 1:(n_seqs-1)) {
    for (j in (i+1):n_seqs) {
      seq1 <- alignment_matrix[i, ]
      seq2 <- alignment_matrix[j, ]
      
      # Count differences (excluding gaps)
      valid_positions <- (seq1 != "-") & (seq2 != "-")
      
      if (sum(valid_positions) > 0) {
        differences <- sum(seq1[valid_positions] != seq2[valid_positions])
        distance <- differences / sum(valid_positions)  # Proportion of differences
      } else {
        distance <- 1  # Maximum distance if no valid positions
      }
      
      dist_matrix[i, j] <- dist_matrix[j, i] <- distance
    }
    
    if (i %% 10 == 0) cat("Processed", i, "sequences...\n")
  }
  
  return(dist_matrix)
}

# Calculate distance matrix
dist_matrix <- calculate_protein_distances(alignment_matrix_clean)

cat("✓ Distance matrix calculated\n")
cat("Distance range:", round(min(dist_matrix[dist_matrix > 0]), 4), "to", round(max(dist_matrix), 4), "\n")

############################################################
# 5. Analyze alignment properties
############################################################

cat("\n=== ALIGNMENT ANALYSIS ===\n")

# Calculate alignment statistics using base R
protein_names <- rownames(alignment_matrix_clean)
sequence_lengths <- apply(alignment_matrix_clean, 1, function(x) sum(x != "-"))
gap_counts <- apply(alignment_matrix_clean, 1, function(x) sum(x == "-"))
gap_proportions <- gap_counts / ncol(alignment_matrix_clean)

alignment_stats <- data.frame(
  Protein = protein_names,
  Length = sequence_lengths,
  Gaps = gap_counts,
  Gap_Proportion = gap_proportions,
  stringsAsFactors = FALSE
)

cat("Alignment statistics:\n")
cat("Mean sequence length:", round(mean(alignment_stats$Length)), "\n")
cat("Mean gap proportion:", round(mean(alignment_stats$Gap_Proportion), 3), "\n")

# Identify conserved and variable positions
position_conservation <- apply(alignment_matrix_clean, 2, function(x) {
  non_gap <- x[x != "-"]
  if (length(non_gap) == 0) return(0)
  max_freq <- max(table(non_gap))
  return(max_freq / length(non_gap))
})

conserved_positions <- sum(position_conservation > 0.9)
variable_positions <- sum(position_conservation < 0.5)

cat("Highly conserved positions (>90% identical):", conserved_positions, "\n")
cat("Highly variable positions (<50% identical):", variable_positions, "\n")

############################################################
# 6. K-means clustering on protein distances
############################################################

cat("\n=== K-MEANS CLUSTERING ON PROTEIN DISTANCES ===\n")

# Determine optimal number of clusters
max_k <- min(8, nrow(dist_matrix) - 1)
wss <- numeric(max_k)
sil_scores <- numeric(max_k - 1)

# Calculate WSS and silhouette scores
for (k in 1:max_k) {
  if (k == 1) {
    wss[k] <- sum(dist_matrix^2) / 2
  } else {
    tryCatch({
      kmeans_result <- kmeans(dist_matrix, centers = k, nstart = 25)
      wss[k] <- kmeans_result$tot.withinss
      
      # Calculate silhouette score for k > 1
      if (k > 1) {
        sil <- silhouette(kmeans_result$cluster, as.dist(dist_matrix))
        sil_scores[k-1] <- mean(sil[, 3])
      }
    }, error = function(e) {
      wss[k] <- NA
      if (k > 1) sil_scores[k-1] <- NA
    })
  }
}

# Find optimal k
valid_wss <- wss[!is.na(wss)]
if (length(valid_wss) > 2) {
  wss_diffs <- diff(valid_wss)
  optimal_k_elbow <- which.min(wss_diffs) + 1
} else {
  optimal_k_elbow <- 2
}

valid_sil <- sil_scores[!is.na(sil_scores)]
if (length(valid_sil) > 0) {
  optimal_k_sil <- which.max(valid_sil) + 1
} else {
  optimal_k_sil <- 2
}

cat("Optimal clusters (Elbow method):", optimal_k_elbow, "\n")
cat("Optimal clusters (Silhouette method):", optimal_k_sil, "\n")

# Test multiple k values
k_values <- unique(c(2, 3, 4, optimal_k_elbow, optimal_k_sil))
k_values <- k_values[k_values <= nrow(dist_matrix) & k_values >= 2]

kmeans_results <- list()

for (k in k_values) {
  cat("Running K-means with k =", k, "\n")
  
  tryCatch({
    kmeans_result <- kmeans(dist_matrix, centers = k, nstart = 50)
    sil <- silhouette(kmeans_result$cluster, as.dist(dist_matrix))
    avg_sil <- mean(sil[, 3])
    
    kmeans_results[[paste0("k", k)]] <- list(
      clusters = kmeans_result$cluster,
      centers = kmeans_result$centers,
      silhouette = avg_sil,
      within_ss = kmeans_result$tot.withinss
    )
    
    cat("  Silhouette score:", round(avg_sil, 3), "\n")
  }, error = function(e) {
    cat("  Failed for k =", k, "\n")
  })
}

# Select best k
if (length(kmeans_results) > 0) {
  silhouette_scores <- sapply(kmeans_results, function(x) x$silhouette)
  best_k <- names(kmeans_results)[which.max(silhouette_scores)]
  best_kmeans <- kmeans_results[[best_k]]
  best_k_num <- as.numeric(gsub("k", "", best_k))
  
  cat("✓ Best clustering: k =", best_k_num, "with silhouette score =", round(best_kmeans$silhouette, 3), "\n")
} else {
  stop("No successful clustering results")
}

############################################################
# 7. Hierarchical clustering
############################################################

cat("\n=== HIERARCHICAL CLUSTERING ===\n")

hclust_result <- hclust(as.dist(dist_matrix), method = "ward.D2")

hclust_clusters <- list()
for (k in k_values) {
  tryCatch({
    clusters <- cutree(hclust_result, k = k)
    sil <- silhouette(clusters, as.dist(dist_matrix))
    avg_sil <- mean(sil[, 3])
    
    hclust_clusters[[paste0("k", k)]] <- list(
      clusters = clusters,
      silhouette = avg_sil
    )
    
    cat("Hierarchical clustering k =", k, ", silhouette =", round(avg_sil, 3), "\n")
  }, error = function(e) {
    cat("Hierarchical clustering failed for k =", k, "\n")
  })
}

############################################################
# 8. Create final class assignments
############################################################

cat("\n=== FINAL CLASS ASSIGNMENTS ===\n")

# Create results
final_assignments <- data.frame(
  Protein = rownames(dist_matrix),
  Kmeans_Cluster = paste0("K", best_kmeans$clusters),
  stringsAsFactors = FALSE
)

# Add sequence statistics
match_indices <- match(final_assignments$Protein, alignment_stats$Protein)
final_assignments$Sequence_Length <- alignment_stats$Length[match_indices]
final_assignments$Gap_Proportion <- alignment_stats$Gap_Proportion[match_indices]

# Add hierarchical clustering if available
if (length(hclust_clusters) > 0) {
  hclust_silhouette_scores <- sapply(hclust_clusters, function(x) x$silhouette)
  best_hclust_k <- names(hclust_clusters)[which.max(hclust_silhouette_scores)]
  best_hclust <- hclust_clusters[[best_hclust_k]]
  final_assignments$Hclust_Cluster <- paste0("H", best_hclust$clusters)
}

# Create dominance class assignments
if (best_k_num == 2) {
  final_assignments$Dominance_Class <- ifelse(
    grepl("K1", final_assignments$Kmeans_Cluster), "Class_I", "Class_II"
  )
} else {
  cluster_numbers <- as.numeric(gsub("K", "", final_assignments$Kmeans_Cluster))
  final_assignments$Dominance_Class <- paste0("Class_", cluster_numbers)
}

# Summary
class_summary <- table(final_assignments$Dominance_Class)
cat("✓ Final dominance class distribution (ingroup only):\n")
print(class_summary)

############################################################
# 9. Analyze cluster characteristics
############################################################

cat("\n=== CLUSTER CHARACTERISTICS ===\n")

# Calculate cluster characteristics
unique_classes <- unique(final_assignments$Dominance_Class)
cluster_analysis <- data.frame(
  Dominance_Class = unique_classes,
  N_proteins = NA,
  Mean_length = NA,
  Mean_gaps = NA,
  stringsAsFactors = FALSE
)

for (i in 1:length(unique_classes)) {
  class_data <- final_assignments[final_assignments$Dominance_Class == unique_classes[i], ]
  
  cluster_analysis$N_proteins[i] <- nrow(class_data)
  cluster_analysis$Mean_length[i] <- round(mean(class_data$Sequence_Length, na.rm = TRUE))
  cluster_analysis$Mean_gaps[i] <- round(mean(class_data$Gap_Proportion, na.rm = TRUE), 3)
}

cat("Cluster characteristics:\n")
print(cluster_analysis)

# Calculate separation metrics
within_cluster_dists <- c()
between_cluster_dists <- c()

for (i in 1:(nrow(dist_matrix)-1)) {
  for (j in (i+1):nrow(dist_matrix)) {
    if (best_kmeans$clusters[i] == best_kmeans$clusters[j]) {
      within_cluster_dists <- c(within_cluster_dists, dist_matrix[i, j])
    } else {
      between_cluster_dists <- c(between_cluster_dists, dist_matrix[i, j])
    }
  }
}

separation_ratio <- mean(between_cluster_dists) / mean(within_cluster_dists)
cat("Within-cluster mean distance:", round(mean(within_cluster_dists), 4), "\n")
cat("Between-cluster mean distance:", round(mean(between_cluster_dists), 4), "\n")
cat("Separation ratio:", round(separation_ratio, 2), "\n")

############################################################
# 10. Visualizations
############################################################

cat("\n=== CREATING VISUALIZATIONS ===\n")

tryCatch({
  pdf("SRK_protein_alignment_clustering_ingroup.pdf", width = 16, height = 12)
  
  par(mfrow = c(3, 3))
  
  # Plot 1: Elbow plot
  plot(1:length(wss), wss, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters K", ylab = "Total WSS",
       main = "Elbow Method (Ingroup Only)")
  if (!is.na(optimal_k_elbow)) abline(v = optimal_k_elbow, col = "red", lty = 2)
  
  # Plot 2: Silhouette scores
  if (length(valid_sil) > 0) {
    plot(2:(length(valid_sil)+1), valid_sil, type = "b", pch = 19, frame = FALSE,
         xlab = "Number of clusters K", ylab = "Average Silhouette Score",
         main = "Silhouette Analysis (Ingroup Only)")
    if (!is.na(optimal_k_sil)) abline(v = optimal_k_sil, col = "red", lty = 2)
  }
  
  # Plot 3: Hierarchical clustering dendrogram
  plot(hclust_result, main = "Hierarchical Clustering (Ingroup)",
       xlab = "SRK Proteins", sub = "", cex = 0.7)
  if (exists("best_hclust") && length(hclust_clusters) > 0) {
    best_hclust_num <- as.numeric(gsub("k", "", best_hclust_k))
    rect.hclust(hclust_result, k = best_hclust_num, border = "red")
  }
  
  # Plot 4: Distance heatmap
  if (nrow(dist_matrix) <= 30) {
    heatmap(dist_matrix, main = "Protein Distance Heatmap (Ingroup)",
            col = heat.colors(50), symm = TRUE, cex.main = 0.8)
  } else {
    sample_indices <- sample(nrow(dist_matrix), min(30, nrow(dist_matrix)))
    heatmap(dist_matrix[sample_indices, sample_indices], 
            main = "Protein Distance Heatmap (Ingroup Sample)",
            col = heat.colors(50), symm = TRUE, cex.main = 0.8)
  }
  
  # Plot 5: Distance distributions
  hist(within_cluster_dists, col = rgb(0,0,1,0.5), 
       main = "Distance Distributions (Ingroup)", xlab = "Protein Distance", 
       xlim = range(c(within_cluster_dists, between_cluster_dists)),
       breaks = 20)
  hist(between_cluster_dists, col = rgb(1,0,0,0.5), add = TRUE, breaks = 20)
  legend("topright", c("Within cluster", "Between cluster"), 
         fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
  
  # Plot 6: Sequence length by cluster
  boxplot(final_assignments$Sequence_Length ~ final_assignments$Dominance_Class,
          main = "Sequence Length by Class (Ingroup)", ylab = "Sequence Length",
          xlab = "Dominance Class")
  
  # Plot 7: Gap proportion by cluster
  boxplot(final_assignments$Gap_Proportion ~ final_assignments$Dominance_Class,
          main = "Gap Proportion by Class (Ingroup)", ylab = "Gap Proportion",
          xlab = "Dominance Class")
  
  # Plot 8: Cluster sizes
  barplot(table(final_assignments$Dominance_Class), 
          main = "Cluster Sizes (Ingroup)", ylab = "Number of Proteins")
  
  # Plot 9: Conservation profile
  plot(1:length(position_conservation), position_conservation, type = "l",
       main = "Sequence Conservation Profile (Ingroup)", 
       xlab = "Alignment Position", ylab = "Conservation Score")
  abline(h = 0.9, col = "red", lty = 2)
  abline(h = 0.5, col = "blue", lty = 2)
  legend("topright", c("Highly conserved", "Variable"), 
         col = c("red", "blue"), lty = 2, cex = 0.8)
  
  dev.off()
  cat("✓ Visualizations saved\n")
  
}, error = function(e) {
  cat("⚠ Visualization failed:", e$message, "\n")
})

############################################################
# 11. Export results
############################################################

cat("\n=== EXPORTING RESULTS ===\n")

# Add protein abundance information
protein_summary <- data.frame(
  Protein = unique(ingroup_protein_key$Protein),
  Count = NA,
  N_ingroup_individuals = NA,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(protein_summary)) {
  protein_rows <- ingroup_protein_key[ingroup_protein_key$Protein == protein_summary$Protein[i], ]
  protein_summary$Count[i] <- protein_rows$Count[1]
  protein_summary$N_ingroup_individuals[i] <- length(unique(protein_rows$Individual))
}

# Merge with final results
final_results <- merge(final_assignments, protein_summary, by = "Protein", all.x = TRUE)

# Save main results
write.table(final_results, "SRK_protein_alignment_class_assignments_ingroup.tsv",
           sep = "\t", quote = FALSE, row.names = FALSE)

# Save other files
write.table(dist_matrix, "SRK_protein_distance_matrix_ingroup.tsv",
           sep = "\t", quote = TRUE, row.names = TRUE)

clustering_summary <- data.frame(
  Method = "K-means_on_protein_alignment_ingroup",
  Total_Samples_Original = nrow(metadata),
  Ingroup_Samples = nrow(ingroup_metadata),
  Total_Proteins_Analyzed = nrow(dist_matrix),
  Optimal_K = best_k_num,
  Silhouette_Score = best_kmeans$silhouette,
  Within_Cluster_Distance = mean(within_cluster_dists),
  Between_Cluster_Distance = mean(between_cluster_dists),
  Separation_Ratio = separation_ratio,
  Alignment_Length = ncol(alignment_matrix_clean)
)

write.table(clustering_summary, "SRK_protein_clustering_summary_ingroup.tsv",
           sep = "\t", quote = FALSE, row.names = FALSE)

write.table(cluster_analysis, "SRK_cluster_characteristics_ingroup.tsv",
           sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
# 12. Final summary
############################################################

cat("\n=== INGROUP-FILTERED PROTEIN ANALYSIS SUMMARY ===\n")
cat("Original samples in metadata:", nrow(metadata), "\n")
cat("Ingroup samples:", nrow(ingroup_metadata), "\n")
cat("Ingroup proteins analyzed:", nrow(dist_matrix), "\n")
cat("Optimal number of clusters:", best_k_num, "\n")
cat("Best silhouette score:", round(best_kmeans$silhouette, 3), "\n")
cat("Cluster separation ratio:", round(separation_ratio, 2), "\n")

cat("\nDominance class distribution (ingroup only):\n")
print(class_summary)

cat("\nFiles created:\n")
if (file.exists("SRK_protein_alignment_clustering_ingroup.pdf")) {
  cat("✓ SRK_protein_alignment_clustering_ingroup.pdf\n")
}
cat("✓ SRK_protein_alignment_class_assignments_ingroup.tsv\n")
cat("✓ SRK_protein_distance_matrix_ingroup.tsv\n")
cat("✓ SRK_protein_clustering_summary_ingroup.tsv\n")
cat("✓ SRK_cluster_characteristics_ingroup.tsv\n")

cat("\n✓ Ingroup-filtered protein analysis complete!\n")
