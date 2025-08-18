## Data Resources

This folder hosts (i) raw model outputs from our pipeline for **Topological associating domain and boundary–defined (TADB) regions** in our analysis that contain an AD GWAS signal and are subsequently tested for associations, and (ii) variant-level tables used to generate manuscript visualizations.

## What’s here

### A) Model files (RDS)

Raw outputs for AD-colocalized TADB regions.  
**Filename convention:** `{event_ID}.{region_ID}.fsusie_mixture_normal_top_pc_weights.rds` where (**`event_ID`** — molecular phenotype + cohort context exactly as used in analysis, e.g., `ROSMAP_DLPFC_haQTL`; **`region_ID`** — TADB cis region identifier, e.g., `chr5_85967320_89904257`). Each RDS contains:
- **`fsusie` models** — default output from `fsusie`, containing the essential components of a `susie` object **and** the fSuSiE-specific regional estimated effect.
- **`susie_top_pc` models** — default output from `susie` from the fine-mapping of the top **1–10 phenotype PCs** derived from the phenotype matrix.
- **Metadata** — the TADB region, the list of variants that go into the analysis, and the genomic coordinates for the input epigenetic marks.
- **Result summaries** — (1) variant-level information for variants in a **95% credible set (CS)**; (2) a trimmed-down fSuSiE object supplemented by **correlations** of the CS identified in fSuSiE.

**Notes**  
There are **30** pair-wise colocalizations being done. This release includes **26** model RDS files; **4** of the pair-wise colocalizations are done within the same TADBs. The **18 loci** colocalized with AD reported in the manuscript are among them.

### B) Variant-level tables (bed.gz)

**FunGen_xQTL_epi.bulk.exported.raw.bed.gz**  
- A compilation of all fSuSiE 95% CS variants with non-zero estimated effects.

**FunGen_xQTL_epi.bulk.exported.rosmap.sentinel_snp_global_maf_filtered_cpg_snp_removed.duplicated_removed.zero_effect_removed.bed.gz**  
- Filtered 95% CS variants.  
  **Criteria:**
  1. Keep only credible sets whose sentinel SNP (max PIP) has **MAF > 0.05** (global).
  2. For mQTL, remove credible sets where any SNP shares the **same coordinates** as a CpG site.
  3. If the **same QTL** is identified in multiple **overlapping TADB** regions, retain **only one** instance (deduplicated).

## Column dictionary (variant tables)

_Types shown as **R / Python** for quick reference; storage for list-like fields is **colon-separated** strings._

**Core**
| Column Name | Type (R / Py) | Description |
|---|---|---|
| `#chr` | integer / int | Chromosome number |
| `start` | integer / int | Genomic start coordinate (**0-based**) |
| `end` | integer / int | Genomic end coordinate (**1-based**) |
| `variant_ID` | character / str | Variant ID `chr:pos:ref:alt` |
| `region_ID` | character / str | id for the TADB cis region |
| `event_ID` | character / str | The molecular phenotype and cohort ID |
| `PIP` | double / float | Posterior inclusion probability |

**Frequency / region / set identifiers**
| Column Name | Type (R / Py) | Description |
|---|---|---|
| `maf` | numeric / float | Minor allele frequency in the testing samples |
| `global_maf` | numeric / float | Minor allele frequency in all the ROSMAP samples |
| `TADB_start` | integer / int | Cis-window start position, typically defined by TAD boundary |
| `TADB_end` | integer / int | Cis-window end position, typically defined by TAD boundary |
| `grid_resolution` | integer / int | Number of grid rows used (typically 512 or 1024) |
| `cs_coverage_0.95` | integer / int | CS index if in **95%** CS with **purity > 0.5**, else 0 |
| `cs_id` | character / str | ID of the 95% credible set within this context and TAD |
| `cs_root` | character / str | Unified root name for overlapping CS in the same context |

**fSuSiE effect-profile fields**
| Column Name | Type (R / Py) | Description |
|---|---|---|
| `grid_positions` | character / str | **Colon-separated numeric** grid positions for the effect estimated by fSuSiE |
| `grid_effects` | character / str | **Colon-separated numeric** effect sizes (`beta`) for those positions |
| `epi_mark_positions` | character / str | **Colon-separated numeric** genomic positions for epigenetic marks within the TADB regions |
| `epi_mark_names` | character / str | **Colon-separated** names for epigenetic marks within the TADB regions |
| `epi_mark_effects` | character / str | **Colon-separated numeric (interpolated)** effects for marks; note there is **no 1:1** mapping to `grid_effects` (interpolation from neighboring grid points) |

## Conventions & provenance
- **Genome build:** `<GRCh38>`;
- **Region keys:** `chr_start_end` in **1-based** closed intervals for human-readable labels; model grids use the `grid_resolution` stated above.
