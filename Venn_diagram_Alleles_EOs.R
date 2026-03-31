###
# Allele overlap between EOs
###

# Load data
# Genotype data
geno <- read.csv("modeling/data/filespopulation27/SRK_individual_protein_table.tsv", sep = "\t")

# Sampling
samp <- read.csv("modeling/data/filespopulation27/sampling_metadata.csv")

pop <- unique(samp$Pop)
list_alleles <- list()
for(i in 1:length(pop)){
  # Subset to pop
  ind <- subset(samp$SampleID, samp$Pop == pop[i])
  # Filter geno to ind
  pop_tmp <- subset(geno, geno$Individual %in% ind)
  #Get alleles
  alleles <- unique(pop_tmp$Protein)
  list_alleles[[i]] <- alleles
}
names(list_alleles) <- pop

# Venn diagram
png("Venn_diagram_Alleles_EOs.png")
ggVennDiagram(list(list_alleles$`27`, list_alleles$`67`, list_alleles$`25`, list_alleles$`76`), category.names = c("EO27", "EO67", "EO25", "EO76"))
dev.off()
