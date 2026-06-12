###############################
# SRK SI barrier plot
###############################

# Read cross results
df <- read.table(
  "SRK_cross_shared_alleles.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Ensure numeric
df$Shared_alleles <- as.numeric(df$Shared_alleles)
df$Success <- as.numeric(df$Success)

## ----------------------------
## Fit logistic model (visualization only)
## ----------------------------

m1 <- glm(
  Success ~ Shared_alleles,
  family = binomial,
  data = df
)

## ----------------------------
## Create prediction grid
## ----------------------------

x_pred <- seq(
  min(df$Shared_alleles, na.rm = TRUE),
  max(df$Shared_alleles, na.rm = TRUE),
  length.out = 200
)

pred_df <- data.frame(Shared_alleles = x_pred)

pred_df$predicted_prob <- predict(
  m1,
  newdata = pred_df,
  type = "response"
)

## ----------------------------
## Empirical success rates
## ----------------------------

empirical <- aggregate(
  Success ~ Shared_alleles,
  data = df,
  mean
)

counts <- aggregate(
  Success ~ Shared_alleles,
  data = df,
  length
)

empirical$N <- counts$Success

## ----------------------------
## PDF output
## ----------------------------

pdf("SRK_SI_barrier_plot.pdf", width = 7, height = 6)

plot(
  df$Shared_alleles,
  df$Success,
  xlab = "Number of shared SRK alleles between parents",
  ylab = "Probability of successful cross",
  main = "Functional SI barrier mediated by SRK genotype",
  pch = 16,
  ylim = c(0, 1)
)

# Jitter raw points so overlaps are visible
points(
  jitter(df$Shared_alleles, amount = 0.05),
  jitter(df$Success, amount = 0.02),
  pch = 16
)

# Add fitted logistic curve
lines(
  pred_df$Shared_alleles,
  pred_df$predicted_prob,
  lwd = 3
)

# Add empirical means
points(
  empirical$Shared_alleles,
  empirical$Success,
  pch = 21,
  cex = 2
)

# Add sample size labels
text(
  empirical$Shared_alleles,
  empirical$Success + 0.08,
  labels = paste0("n=", empirical$N),
  cex = 0.9
)

dev.off()

cat("✅ SI barrier plot written to: SRK_SI_barrier_plot.pdf\n")

