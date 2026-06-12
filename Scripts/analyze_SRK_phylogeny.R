############################################################
# SRK Phylogenetic Analysis and Dominance Class Prediction
# Enhanced with automatic package installation for Mac
############################################################

cat("Starting SRK phylogenetic analysis\n")

############################################################
# 0. Package installation and loading (Mac-optimized)
############################################################

cat("=== CHECKING AND INSTALLING REQUIRED PACKAGES ===\n")

# Function to check and install packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      
      # Try CRAN first
      tryCatch({
        install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
        library(pkg, character.only = TRUE)
        cat("✓ Successfully installed", pkg, "from CRAN\n")
      }, error = function(e) {
        cat("✗ Failed to install", pkg, "from CRAN\n")
        cat("Error:", e$message, "\n")
      })
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# Function to install Bioconductor packages
install_bioconductor_if_missing <- function(packages) {
  # Install BiocManager if not available
  if (!require("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cran.rstudio.com/")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing Bioconductor package:", pkg, "\n")
      tryCatch({
        BiocManager::install(pkg, dependencies = TRUE, ask = FALSE)
        library(pkg, character.only = TRUE)
        cat("✓ Successfully installed", pkg, "from Bioconductor\n")
      }, error = function(e) {
        cat("✗ Failed to install", pkg, "from Bioconductor\n")
        cat("Error:", e$message, "\n")
      })
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# Required CRAN packages
cran_packages <- c(
  "ape",        # Phylogenetic analysis
  "phytools",   # Phylogenetic tools
  "dplyr",      # Data manipulation
  "ggplot2"     # Plotting
)

# Required Bioconductor packages
bioc_packages <- c(
  "ggtree"      # Tree visualization
)

# Install packages
cat("Installing CRAN packages...\n")
install_if_missing(cran_packages)

cat("\nInstalling Bioconductor packages...\n")
install_bioconductor_if_missing(bioc_packages)

# Load all packages
cat("\n=== LOADING PACKAGES ===\n")
required_packages <- c(cran_packages, bioc_packages)

for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓ Loaded:", pkg, "\n")
  } else {
    cat("✗ Failed to load:", pkg, "\n")
    cat("Please try installing manually:\n")
    if (pkg %in% bioc_packages) {
      cat("BiocManager::install('", pkg, "')\n", sep = "")
    } else {
      cat("install.packages('", pkg, "')\n", sep = "")
    }
    stop("Package loading failed")
  }
}

cat("\n✓ All packages loaded successfully!\n\n")

############################################################
# 1. Load tree and data
############################################################

cat("=== LOADING DATA ===\n")

# Check if required files exist
required_files <- c("SRK_phylogeny.contree", "SRK_functional_protein_key.tsv")

for (file in required_files) {
  if (!file.exists(file)) {
    stop("Required file not found: ", file)
  }
}

# Load IQ-TREE results
tryCatch({
  tree <- read.tree("SRK_phylogeny.contree")
  cat("✓ Tree loaded with", Ntip(tree), "tips\n")
}, error = function(e) {
  stop("Failed to load tree file: ", e$message)
})

# Load protein key for metadata
tryCatch({
  protein_key <- read.table("SRK_functional_protein_key.tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("✓ Protein key loaded with", nrow(protein_key), "entries\n")
}, error = function(e) {
  stop("Failed to load protein key file: ", e$message)
})

############################################################
# 2. Root the tree (if possible)
############################################################

cat("\n=== ROOTING TREE ===\n")

tryCatch({
  # Option 1: Midpoint rooting
  tree_rooted <- midpoint.root(tree)
  cat("✓ Tree rooted using midpoint method\n")
  
  # Option 2: If you have outgroup sequences, uncomment and modify:
  # outgroup_tips <- c("outgroup1", "outgroup2")  # specify outgroup names
  # tree_rooted <- root(tree, outgroup_tips)
  # cat("✓ Tree rooted using outgroup method\n")
  
}, error = function(e) {
  cat("⚠ Rooting failed, using original tree:", e$message, "\n")
  tree_rooted <- tree
})

############################################################
# 3. Extract individual/population information
############################################################

cat("\n=== PROCESSING METADATA ===\n")

# Extract metadata from protein key
get_individual_info <- function(protein_key) {
  tryCatch({
    # Extract individual and library info from sequence IDs
    protein_key$Individual <- sapply(strsplit(protein_key$Original_sequence_ID, "_"), 
                                    function(x) paste(x[1:2], collapse="_"))
    
    protein_key$Library <- sapply(strsplit(protein_key$Original_sequence_ID, "_"), 
                                 function(x) x[1])
    
    return(protein_key)
  }, error = function(e) {
    cat("⚠ Warning: Could not extract individual info:", e$message, "\n")
    return(protein_key)
  })
}

protein_metadata <- get_individual_info(protein_key)

# Summarize protein distribution
tryCatch({
  protein_summary <- protein_metadata %>%
    group_by(Protein) %>%
    summarise(
      Count = first(Count),
      N_individuals = n_distinct(Individual),
      N_libraries = n_distinct(Library),
      .groups = 'drop'
    )
  
  cat("✓ Protein summary created\n")
  print(head(protein_summary))
  
}, error = function(e) {
  cat("⚠ Warning: Could not create protein summary:", e$message, "\n")
  protein_summary <- data.frame(
    Protein = unique(protein_metadata$Protein),
    Count = 1,
    N_individuals = 1,
    N_libraries = 1
  )
})

############################################################
# 4. Identify major clades
############################################################

cat("\n=== IDENTIFYING MAJOR CLADES ===\n")

identify_major_clades <- function(tree, min_bootstrap = 70) {
  tryCatch({
    # Get bootstrap values
    bootstrap_values <- as.numeric(tree$node.label)
    bootstrap_values <- bootstrap_values[!is.na(bootstrap_values)]
    
    if (length(bootstrap_values) == 0) {
      cat("⚠ No bootstrap values found in tree\n")
      return(list())
    }
    
    # Find well-supported nodes
    well_supported <- which(as.numeric(tree$node.label) >= min_bootstrap)
    
    # Get descendants for each well-supported node
    major_clades <- list()
    clade_count <- 0
    
    for (i in seq_along(well_supported)) {
      node <- well_supported[i] + Ntip(tree)  # Adjust for tip numbering
      
      if (node <= (Ntip(tree) + Nnode(tree))) {  # Valid node check
        descendants <- extract.clade(tree, node)$tip.label
        
        if (length(descendants) >= 3) {  # Only consider clades with 3+ members
          clade_count <- clade_count + 1
          major_clades[[paste0("Clade_", clade_count)]] <- descendants
        }
      }
    }
    
    return(major_clades)
    
  }, error = function(e) {
    cat("⚠ Error identifying clades:", e$message, "\n")
    cat("Creating simple clades based on tree structure...\n")
    
    # Fallback: create clades based on tree structure
    n_tips <- Ntip(tree)
    if (n_tips >= 6) {
      # Simple division into 2 clades
      tips <- tree$tip.label
      mid <- ceiling(n_tips / 2)
      return(list(
        Clade_1 = tips[1:mid],
        Clade_2 = tips[(mid+1):n_tips]
      ))
    } else {
      return(list())
    }
  })
}

major_clades <- identify_major_clades(tree_rooted, min_bootstrap = 70)

cat("✓ Major clades identified:", length(major_clades), "\n")
for (i in seq_along(major_clades)) {
  cat("  ", names(major_clades)[i], ":", length(major_clades[[i]]), "proteins\n")
}

############################################################
# 5. Assign dominance classes based on phylogeny
############################################################

cat("\n=== ASSIGNING DOMINANCE CLASSES ===\n")

assign_dominance_classes <- function(tree, major_clades) {
  # Create dominance class assignments
  dominance_classes <- data.frame(
    Protein = tree$tip.label,
    Dominance_Class = "Unassigned",
    Clade = "None",
    stringsAsFactors = FALSE
  )
  
  # Assign clades
  for (i in seq_along(major_clades)) {
    clade_name <- names(major_clades)[i]
    clade_members <- major_clades[[i]]
    
    dominance_classes$Clade[dominance_classes$Protein %in% clade_members] <- clade_name
    
    # Simple assignment: Clade_1 = Class I (dominant), Clade_2 = Class II (recessive)
    if (i == 1) {
      dominance_classes$Dominance_Class[dominance_classes$Protein %in% clade_members] <- "Class_I"
    } else if (i == 2) {
      dominance_classes$Dominance_Class[dominance_classes$Protein %in% clade_members] <- "Class_II"
    } else {
      dominance_classes$Dominance_Class[dominance_classes$Protein %in% clade_members] <- paste0("Class_", i)
    }
  }
  
  return(dominance_classes)
}

dominance_assignments <- assign_dominance_classes(tree_rooted, major_clades)

# Summary of dominance classes
class_summary <- table(dominance_assignments$Dominance_Class)
cat("✓ Dominance class distribution:\n")
print(class_summary)

############################################################
# 6. Visualize phylogeny with dominance classes
############################################################

cat("\n=== CREATING VISUALIZATIONS ===\n")

tryCatch({
  # Create color scheme for classes
  n_classes <- length(unique(dominance_assignments$Dominance_Class))
  class_colors <- rainbow(n_classes)
  names(class_colors) <- unique(dominance_assignments$Dominance_Class)
  
  # Basic tree plot
  pdf("SRK_phylogeny_dominance_classes.pdf", width = 12, height = 10)
  
  # Plot 1: Circular tree with dominance classes
  p1 <- ggtree(tree_rooted, layout = "circular") %<+% dominance_assignments +
    geom_tippoint(aes(color = Dominance_Class), size = 2) +
    geom_tiplab(size = 2, offset = 0.01) +
    scale_color_manual(values = class_colors) +
    theme(legend.position = "bottom") +
    ggtitle("SRK Phylogeny with Predicted Dominance Classes")
  
  print(p1)
  
  # Plot 2: Rectangular tree with bootstrap support
  p2 <- ggtree(tree_rooted) %<+% dominance_assignments +
    geom_tippoint(aes(color = Dominance_Class), size = 3) +
    geom_nodelab(aes(label = label), size = 2, hjust = 1.2) +  # Bootstrap values
    geom_tiplab(size = 2.5) +
    scale_color_manual(values = class_colors) +
    theme_tree2() +
    ggtitle("SRK Phylogeny with Bootstrap Support")
  
  print(p2)
  
  # Plot 3: Tree with protein abundance information
  if ("Count" %in% colnames(protein_summary)) {
    protein_abundance <- protein_summary$Count
    names(protein_abundance) <- protein_summary$Protein
    
    p3 <- ggtree(tree_rooted) %<+% dominance_assignments +
      geom_tippoint(aes(color = Dominance_Class, 
                       size = protein_abundance[label]), alpha = 0.7) +
      geom_tiplab(size = 2) +
      scale_color_manual(values = class_colors) +
      scale_size_continuous(name = "Protein\nAbundance", range = c(1, 6)) +
      ggtitle("SRK Phylogeny with Protein Abundance")
    
    print(p3)
  }
  
  dev.off()
  
  cat("✓ Visualizations saved to: SRK_phylogeny_dominance_classes.pdf\n")
  
}, error = function(e) {
  cat("⚠ Error creating visualizations:", e$message, "\n")
  cat("Tree analysis will continue without plots\n")
})

############################################################
# 7. Export results
############################################################

cat("\n=== EXPORTING RESULTS ===\n")

# Save dominance class assignments
write.table(dominance_assignments, "SRK_dominance_class_assignments.tsv",
           sep = "\t", quote = FALSE, row.names = FALSE)
cat("✓ Saved: SRK_dominance_class_assignments.tsv\n")

# Save detailed results with protein metadata
detailed_results <- merge(dominance_assignments, protein_summary, by = "Protein", all.x = TRUE)
write.table(detailed_results, "SRK_phylogeny_detailed_results.tsv", 
           sep = "\t", quote = FALSE, row.names = FALSE)
cat("✓ Saved: SRK_phylogeny_detailed_results.tsv\n")

# Save major clades information
if (length(major_clades) > 0) {
  clade_info <- data.frame(
    Clade = rep(names(major_clades), sapply(major_clades, length)),
    Protein = unlist(major_clades),
    stringsAsFactors = FALSE
  )
  write.table(clade_info, "SRK_major_clades.tsv",
             sep = "\t", quote = FALSE, row.names = FALSE)
  cat("✓ Saved: SRK_major_clades.tsv\n")
}

############################################################
# 8. Summary statistics
############################################################

cat("\n=== PHYLOGENETIC ANALYSIS SUMMARY ===\n")
cat("Total SRK proteins analyzed:", Ntip(tree), "\n")
cat("Major clades identified:", length(major_clades), "\n")
cat("Dominance classes assigned:\n")
print(class_summary)

cat("\nFiles created:\n")
if (file.exists("SRK_phylogeny_dominance_classes.pdf")) {
  cat("✓ SRK_phylogeny_dominance_classes.pdf\n")
}
cat("✓ SRK_dominance_class_assignments.tsv\n")
cat("✓ SRK_phylogeny_detailed_results.tsv\n")
if (length(major_clades) > 0) {
  cat("✓ SRK_major_clades.tsv\n")
}

cat("\nNext steps:\n")
cat("1. Examine phylogeny for clear class divisions\n")
cat("2. Test dominance predictions with crossing data\n")
cat("3. Refine class assignments based on biological knowledge\n")

cat("\n✓ Analysis complete!\n")
