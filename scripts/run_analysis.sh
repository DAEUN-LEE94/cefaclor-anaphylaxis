#!/bin/bash
set -euo pipefail

##############################################
# STEP 1: VCF Preprocessing
##############################################

# 1-1. Decompress VCF
gzip -dk merged.vcf.gz

# 1-2. Remove chrX, chrY, chrM (use predefined autosomes list)
bcftools view -T chrs_autosome.txt -Oz -o autosome.vcf.gz merged.vcf.gz
gzip -dk autosome.vcf.gz

# 1-3. Apply variant-level quality filters: DP ≥ 10, QUAL ≥ 30, GQ ≥ 30
bcftools filter -i 'INFO/DP>=10 & QUAL>=30 & FORMAT/GQ>=30' autosome.vcf -o merged_qc.vcf

# 1-4. Remove multiallelic variants
bgzip -c merged_qc.vcf > merged_qc.vcf.gz
bcftools index merged_qc.vcf.gz
bcftools filter -i 'N_ALT=1' merged_qc.vcf.gz -o merged_biallelic.vcf

# 1-5. HWE filtering (p > 0.05)
plink2 --vcf merged_biallelic.vcf --make-bed --pheno phenotype.txt --out plink_merged
plink2 --bfile plink_merged --hwe 0.05 --make-bed --out plink_merged_hwe

##############################################
# STEP 2: MAF Filtering (MAF ≥ 0.01)
##############################################

plink2 --bfile plink_merged --maf 0.01 --make-bed --out plink_merged_maf

##############################################
# STEP 3: Allele Frequencies
##############################################

# 3-1. Total
plink2 --bfile plink_merged_maf --freq --out plink_all_AF

# 3-2. Case-only
plink2 --bfile plink_merged_maf --freq --keep case_samples.txt --out plink_case_AF

# 3-3. Control-only
plink2 --bfile plink_merged_maf --freq --keep control_samples.txt --out plink_control_AF

##############################################
# STEP 4: Genotype Counts
##############################################

# 4-1. Total
plink2 --bfile plink_merged_maf --geno-counts --out plink_all_GC

# 4-2. Case-only
plink2 --bfile plink_merged_maf --geno-counts --keep case_samples.txt --out plink_case_GC

# 4-3. Control-only
plink2 --bfile plink_merged_maf --geno-counts --keep control_samples.txt --out plink_control_GC

##############################################
# STEP 5: Final VCF Export for R Integration
##############################################

plink2 --bfile plink_merged_maf --export vcf --out plink_merged_maf_final

##############################################
# STEP 6: Association Analysis
##############################################

# covar_file.txt: must include FID, IID, and sex
# phenotype_file.txt: must include FID, IID, and phenotype (0/1)
plink2 --bfile plink_merged_maf --glm covar --covar covar_file.txt --pheno phenotype_file.txt --out glm_results

