# fSuSiE Annotation and Enrichment Results

## 1. Annotation Files
### ENCODE Annotations & Processing
- **105 BED files** originally curated by Xin He's lab in **hg19**, lifted over to **hg38** (removing redundant annotations).
- **Types of annotations used**:
  - **CS-SNPs with PIP > 0.5**
  - **CS-peak regions (200bp windows)**
  - **TensorQTL qval 0.05 cutoff**

### Processed Annotation Files
| Molecular Trait | Data Type | Annotation File |
|----------------|-----------|----------------------------------|
| fSuSiE haQTL CS Peak | Peak-Based | `ROSMAP_haQTL_cs_effect_ha_peak_annotation.tsv.gz` |
| fSuSiE haQTL CS SNP | SNP-Based | `ROSMAP_haQTL_cs_snp_annotation.tsv.gz` |
| SuSiE haQTL CS topPC1 SNP | SNP-Based | `ROSMAP_haQTL_cs_snp_toppc1_annotation.tsv.gz` |
| TensorQTL haQTL SNP | SNP-Based | `ROSMAP_haQTL_qtl_snp_qval0.05_annotation.tsv.gz` |
| fSuSiE mQTL CS CpG | CpG-Based | `ROSMAP_mQTL_cs_effect_cpg_annotation.tsv.gz` |
| fSuSiE mQTL CS SNP | SNP-Based | `ROSMAP_mQTL_cs_snp_annotation.tsv.gz` |
| SuSiE mQTL CS topPC1 SNP | SNP-Based | `ROSMAP_mQTL_cs_snp_toppc1_annotation.tsv.gz` |
| TensorQTL mQTL SNP（ongoing） | SNP-Based | `ROSMAP_mQTL_qtl_snp_qval0.05_annotation.tsv.gz` |
---

## 2. SLDSC Heritability Enrichment Results
### Background Annotations
- **97 pre-computed annotation tracks** from **Alkes Price lab** were used for SLDSC.
- We performed SLDSC **with all background annotations** and **removing data-related background annotations**.
- Results **remained consistent** between the two versions, so we primarily report outputs **using all annotation backgrounds**.

### SLDSC Enrichment Data Summary
#### Datasets Included in Analysis
- **haQTL:**
  - fSuSiE **haQTL CS Peak**
  - fSuSiE **haQTL CS SNP** (maxPIP, inCS)
  - SuSiE **haQTL CS topPC1 SNP** (maxPIP, inCS)
  - **TensorQTL haQTL SNP**
- **mQTL:**
  - fSuSiE **mQTL CS CpG**
  - fSuSiE **mQTL CS SNP** (maxPIP, inCS)
  - SuSiE **mQTL CS topPC1 SNP** (maxPIP, inCS)
  - **TensorQTL mQTL SNP**

#### CS-SNP Level Format

1. maxPIP: Using PIP as a continuous annotation score, assigning max PIP for the same SNP across different CSs.
2. inCS: Using binary annotation where all CS-SNPs are assigned a value of 1.
3. topPC1: Using the first principal component (TOP PC1) of SNP-level annotations to capture major variation.


### Meta-Analysis Results
- **File:** `sldsc_all_methods_summary_allbg.tsv.gz`
- **Description:** This file contains **meta-results** summarizing **heritability enrichment** across:
  - **Overall studies (57)**
  - **Blood-specific studies (22)**
  - **Brain-specific studies (18)**

### Per-Trait Results
- Each trait's SLDSC results are stored as intermediate files before meta-analysis.
- Key outputs include:
  - `data$<trait_name>$enrichment$enrichment_summary`: Heritability enrichment results for each trait.
  - `data$<trait_name>$enrichment$meta_enrich_mean_sd_input_stats`: Used for meta-analysis to compute **meta enrichment mean and standard deviation**.
  - `data$<trait_name>$enrichment$meta_enrich_pvalue_input_stats`: Used for meta-analysis to compute **meta enrichment p-value**.

### Task Files
- **TASKFILES:** Provide the trait names categorized into **ALL**, **Brain**, and **Blood** study groups.


## 3. Over-Representation Analysis (ORA) Results
### Background & Method
- 105 BED annotation files from Xin He's lab were used.
- ORA was performed on **CS-SNP and peak-level annotations** to assess over-representation using **odds ratio analysis**.
- Statistical Testing: Jackknife resampling was used for **mean and standard error (SE) estimation**.

### **ORA Data Summary**
| Molecular Trait | Data Type | Enrichment File |
|----------------|-----------|----------------------------------|
| fSuSiE haQTL CS Peak | Peak-Based | `haQTL_cs_peak_200bp.range_enrichment_stats_summary.tsv.gz` |
| fSuSiE haQTL CS SNP | SNP-Based | `haQTL_cs_snp_pip0.5.enrichment_results_summary.tsv.gz` |
| SuSiE haQTL CS topPC1 SNP | SNP-Based | `haQTL_cs_snp_toppc1_pip0.5.enrichment_results_summary.tsv.gz` |
| TensorQTL haQTL SNP | SNP-Based | `haQTL_qtl_snp_qvalue0.05.enrichment_results_summary.tsv.gz` |
| fSuSiE mQTL CS SNP | SNP-Based | `mQTL_cs_snp_pip0.5.enrichment_results_summary.tsv.gz` |
| SuSiE mQTL CS topPC1 SNP | SNP-Based | `mQTL_cs_snp_toppc1_pip0.5.enrichment_results_summary.tsv.gz` |
| TensorQTL mQTL SNP | SNP-Based | `mQTL_qtl_snp_qvalue0.05.enrichment_results_summary.tsv.gz` |


## Summary
1. **Annotation Framework:**
   - We used ENCODE annotations from Xin He's lab, lifted to hg38** and filtered redundant features.
   - Analyzed CS-SNPs (PIP > 0.5), CS-peaks (200bp), CpG site, and TensorQTL q-value < 0.05 SNPs.

2. **SLDSC Heritability Enrichment:**
   - 97 annotation backgrounds were used, and trends were consistent whether or not data-related annotations were removed.
   - Meta-analysis of 57 total GWAS, 22 blood-specific GWAS, and 18 brain-specific GWAS.

3. **ORA Functional Enrichment:**
   - CS-SNP and CS-Peak regions showed enrichment trends.
   - Depletion observed in repressor-related annotations, supporting valid enrichment estimates.

## Documentation notebook:
`fSuSiE_functional_enrichment_Feb12.html`