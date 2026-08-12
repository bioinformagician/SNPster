library(biomaRt)


# Install Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install biomaRt for querying Ensembl
BiocManager::install("biomaRt")

# Load the package
library(biomaRt)
library(data.table)

gwas_data = "C:/Users/frezz/Desktop/genome_snp_project/vacation genomics project/gwas_data_filtered.csv"

gwas_data = fread(gwas_data)


snps=gwas_data$rsid

# Connect to the Ensembl database
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")


# Query Ensembl for strand and alleles information
results <- getBM(
  attributes = c("refsnp_id", "allele", "chrom_start", "chrom_strand"),
  filters = "snp_filter",
  values = snps,
  mart = ensembl
)

# Display the results
print(results)
View(results)
head(results)

# 1 means forward strand, -1 means reverse strand


