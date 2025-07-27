# Example analysis workflow to visualize association results
# Replace placeholder paths, file names, and sample identifiers with your own values

# Set working directory containing VCF files
setwd("/path/to/working/directory")

# Load required libraries
library(vcfR)
library(SNPRelate)
library(gdsfmt)
library(ellipse)
library(qqman)

# Path to input VCF file
vcf_path <- "path/to/input.vcf"

# Convert VCF to GDS format
snpgdsVCF2GDS(vcf_path, "temp.gds")

# Open the generated GDS file
genofile <- snpgdsOpen("temp.gds")

# Retrieve sample identifiers
sample_info <- read.gdsn(index.gdsn(genofile, "sample.id"))

# List of case sample IDs
case_samples <- c("CASE_SAMPLE_01", "CASE_SAMPLE_02")

# Mark each sample as case or control
case_control_vector <- ifelse(sample_info %in% case_samples, "case", "control")

# Perform PCA
pca <- snpgdsPCA(genofile, snp.id = genofile$snp.id, num.thread = 4)

# Variance explained by the first two principal components
pc1_var <- round(pca$varprop[1] * 100, 2)
pc2_var <- round(pca$varprop[2] * 100, 2)

# Build data frame for plotting
pca_df <- data.frame(
  PC1 = pca$eigenvect[, 1],
  PC2 = pca$eigenvect[, 2],
  Group = factor(case_control_vector, levels = c("case", "control"))
)

# Create PCA plot with 95% confidence ellipse
pdf("PCA_plot_with_variance_and_ellipse.pdf")
plot(pca_df$PC1, pca_df$PC2,
     col = ifelse(pca_df$Group == "case", "red", "blue"),
     pch = ifelse(pca_df$Group == "case", 19, 1),
     xlab = paste0("PC1 (", pc1_var, "% variance)"),
     ylab = paste0("PC2 (", pc2_var, "% variance)"))

for (group in levels(pca_df$Group)) {
  points_in_group <- pca_df[pca_df$Group == group, c("PC1", "PC2")]
  if (nrow(points_in_group) > 2) {
    lines(ellipse(cov(points_in_group),
                  centre = colMeans(points_in_group),
                  level = 0.95),
          col = ifelse(group == "case", "red", "blue"))
  }
}
legend("topright", legend = c("Case", "Control"),
       col = c("red", "blue"), pch = c(19, 1))
dev.off()

# Close the GDS file
snpgdsClose(genofile)

# Example association results data frame
# 'result' should include columns: chromosome (chr), position (bp), ID (snp) and p-value (P_Value)
# result <- read.table("path/to/association_results.txt", header = TRUE)

# Create Manhattan plot
png("manhattan_plot.png", width = 1200, height = 600)
manhattan(result, chr = "X.CHROM", bp = "POS", snp = "ID", p = "P_Value",
          col = c("orange3", "blue4"), main = "Manhattan Plot",
          genomewideline = FALSE)
# Optional significance line at -log10(p) = 5
abline(h = 5, col = "red", lwd = 2)
dev.off()

# Create Q-Q plot
png("qq_plot.png", width = 1000, height = 1000)
qq(result$P_Value, main = "Q-Q plot")
dev.off()
