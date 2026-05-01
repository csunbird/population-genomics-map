# PopGen Map — User Guide

**Version 0.6.0 · Phase 5 — Export & Documentation**
_Open-source conservation genomics visualisation platform_

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Quick Start](#2-quick-start)
3. [Input File Formats](#3-input-file-formats)
4. [Phase 1 — Genetic Diversity Mapping](#4-phase-1--genetic-diversity-mapping)
5. [Phase 2 — Population Structure](#5-phase-2--population-structure)
6. [Phase 3 — Inbreeding & Effective Population Size](#6-phase-3--inbreeding--effective-population-size)
7. [Phase 4 — Adaptive Potential](#7-phase-4--adaptive-potential)
8. [Phase 5 — Export & Reporting](#8-phase-5--export--reporting)
9. [Interpreting Results](#9-interpreting-results)
10. [Troubleshooting](#10-troubleshooting)
11. [Glossary](#11-glossary)
12. [Citation & Licence](#12-citation--licence)

---

## 1. Introduction

**PopGen Map** turns population genomic data into interactive geographic maps — designed for conservation practitioners who need to understand the genetic health of wildlife populations, not bioinformaticians who already know how.

Upload a standard VCF file and a spreadsheet of GPS coordinates, and the app produces:

- An interactive map coloured by genetic diversity, structure, inbreeding, effective population size, and climate vulnerability
- Plain-language conservation risk assessments for each population
- Downloadable reports suitable for management plans

**No bioinformatics experience is needed to use the app.** No infrastructure is needed to run it — just R and a web browser.

### What PopGen Map is _not_

PopGen Map is a **visualisation and summary tool**. It does not replace dedicated population genomics pipelines (GATK, STACKS, ADMIXTURE, etc.) for primary analysis. Its inputs are the outputs of those tools.

---

## 2. Quick Start

### Prerequisites

- R 4.2 or later ([r-project.org](https://www.r-project.org))
- A web browser (Chrome, Firefox, or Edge recommended)

### Installation

```r
# 1. Open R or RStudio and navigate to the popgen-map folder
setwd("path/to/popgen-map")

# 2. Install all dependencies (run once)
source("install_packages.R")

# 3. Launch the app
shiny::runApp()
```

The app opens in your browser. If it does not open automatically, visit `http://127.0.0.1:PORT` as shown in the R console.

### First steps with demo data

Click **"Load demo dataset"** in the sidebar. This loads a simulated five-population large felid dataset and demonstrates all Phase 1–5 features immediately — no files required.

> ⚠️ Demo data is simulated. Values are illustrative only.

---

## 3. Input File Formats

### 3.1 Genotype data — VCF file

| Property | Requirement |
|---|---|
| Format | `.vcf` or `.vcf.gz` (gzip-compressed) |
| Variant type | Biallelic SNPs only |
| Variant callers | Any: GATK, STACKS, FreeBayes, DeepVariant, bcftools |
| Sample names | Must match the `sample_id` column in your metadata CSV |
| Size | Up to a few hundred MB; app processes in memory |

**Preparing your VCF:**

If your VCF contains multiallelic sites or indels, filter them first:

```bash
# Keep biallelic SNPs, remove monomorphic sites
bcftools view -m2 -M2 -v snps your_data.vcf.gz -Oz -o filtered.vcf.gz
```

### 3.2 Sample metadata — CSV file

Required columns (column names are case-insensitive; common variants accepted):

| Column | Accepted names | Description |
|---|---|---|
| `sample_id` | `sampleid`, `id`, `individual`, `ind` | Unique sample identifier — must match VCF sample names |
| `population` | `pop`, `pop_id`, `site`, `location` | Population or site label |
| `latitude` | `lat`, `y` | Decimal degrees (WGS84) |
| `longitude` | `lon`, `long`, `x` | Decimal degrees (WGS84) |

**Example:**

```csv
sample_id,population,latitude,longitude
TIG_001,Central Highland,4.56,116.42
TIG_002,Central Highland,4.56,116.42
TIG_003,Northern Forest,6.12,117.83
```

Click **"Metadata CSV template"** in the sidebar to download a pre-filled template.

> Each population must have at least 3 samples for meaningful diversity estimates. Fewer than 6 samples per population will produce wide uncertainty intervals.

### 3.3 Environmental data — CSV file (Phase 4 only)

Required for GEA and genomic offset analyses. One row per population:

| Column | Description |
|---|---|
| `population` | Must match population names in your metadata CSV |
| `BIO1`, `BIO12`, … | WorldClim / CMIP6 bioclimatic variable values |

You can also use the **"Auto-fetch from WorldClim"** option in the sidebar to have the app pull climate variables automatically from population GPS coordinates (requires internet access).

---

## 4. Phase 1 — Genetic Diversity Mapping

### What it computes

| Metric | Formula | What it means |
|---|---|---|
| **He** | 2p(1−p), averaged across loci | Expected heterozygosity — how genetically diverse a population _should_ be under random mating |
| **Ho** | Proportion of heterozygous genotypes | Observed heterozygosity — how diverse individuals actually _are_ |
| **π** | Nei (1978) diploid-corrected He | Nucleotide diversity — average differences per site between two randomly chosen sequences |
| **F** | (He − Ho) / He | Inbreeding coefficient — deviation from Hardy-Weinberg equilibrium |

### Map tab

- **Colour by**: choose He, Ho, π, F, Ne (if estimated), or sample size (n) from the dropdown
- **Marker size**: fixed, proportional to n (sample size), or proportional to He
- **Fit to data**: re-centres and zooms the map to your populations
- **FST network**: shows pairwise genetic differentiation as colour-coded edges between populations
- **ADMIXTURE pies**: overlays ancestry proportion pie charts (requires Q-matrix upload in the ADMIXTURE tab)

### Reading the map

Click any population marker to open its **popup card**, which shows all statistics and a plain-language risk assessment. Hover for a quick summary.

### Conservation thresholds (He)

These thresholds are calibrated on microsatellite and allozyme data (Frankham et al. 2014). SNP-based He filtered at MAF ≥ 0.05 typically runs 20–40% lower than microsatellite He for the same population. **Use as indicative guidance, not hard cut-offs.**

| He range | Status | Guidance |
|---|---|---|
| < 0.10 | 🔴 Critical | High extinction risk; urgent intervention required |
| 0.10 – 0.15 | 🟠 High risk | Low diversity; translocation or breeding programme likely needed |
| 0.15 – 0.20 | 🟡 Moderate risk | Monitoring and genetic rescue assessment recommended |
| 0.20 – 0.30 | 🟢 Lower risk | Moderate diversity; routine monitoring |
| > 0.30 | 🔵 Healthy | High diversity; candidate reference population |

### Statistics tab

A sortable table of all per-population metrics. Use **"Export CSV"** to download the full table.

---

## 5. Phase 2 — Population Structure

### FST Heatmap (Structure tab)

Displays the pairwise Hudson F_ST matrix as a colour-coded heatmap. Each cell shows the genetic differentiation between two populations.

- **Hover** over a cell to see the exact F_ST value and Wright (1978) category
- **Wright categories**: < 0.05 = little | 0.05–0.15 = moderate | 0.15–0.25 = great | > 0.25 = very great differentiation

High F_ST between populations that are geographically close may indicate landscape barriers (roads, rivers, habitat loss). Very high F_ST (> 0.25) combined with low He is a warning sign of long-term isolation.

### PCA (PCA tab)

Principal Component Analysis of individual genotypes. Each dot is an individual, coloured by population. Tight clusters = low within-population diversity; large separation between clusters = high F_ST.

- Use PC1 vs PC2 as default; toggle to PC3 for additional structure
- Outlier individuals may indicate misidentified samples or admixed individuals

### ADMIXTURE (ADMIXTURE tab)

Upload a Q-matrix CSV from an external ADMIXTURE run (K ancestry components per individual). The app renders:

- A **stacked bar chart** of ancestry proportions per population
- **Pie charts** overlaid on the map (toggle with "Show ADMIXTURE pies" in the Map tab)

**Running ADMIXTURE externally:**

```bash
# Convert VCF to PLINK format first
plink --vcf your_data.vcf --make-bed --out your_data --allow-extra-chr

# Run ADMIXTURE for K = 2 to 8, choose K by cross-validation error
for K in 2 3 4 5 6 7 8; do
  admixture --cv your_data.bed $K | tee admixture_K${K}.log
done

# Use the K with lowest CV error
grep "CV error" admixture_K*.log
```

The Q-matrix file should have columns: `sample_id, K1, K2, …, Kn`.

---

## 6. Phase 3 — Inbreeding & Effective Population Size

### F-ROH (F-ROH tab)

Runs of Homozygosity (ROH) are long stretches of identical alleles inherited from a common ancestor. F-ROH = total ROH length / autosome length. It captures **recent** inbreeding (last ~10 generations), unlike the F coefficient which reflects deviation from HWE.

| F-ROH threshold | Approximate relationship |
|---|---|
| ≥ 0.25 | Full-sibling or parent-offspring mating |
| ≥ 0.125 | Half-sibling mating |
| ≥ 0.0625 | First-cousin mating |
| < 0.0625 | Moderate background inbreeding |

> ROH detection requires dense SNP coverage (> 5,000 SNPs evenly distributed across chromosomes). Results are unreliable with RADseq or very sparse panels.

### Effective Population Size — Ne (Ne tab)

Ne is the size of an idealised random-mating population that experiences the same rate of genetic drift as the real population. It is always smaller than census size (N) — typically 10–20% of N in wildlife.

The app estimates Ne using the **single-sample linkage disequilibrium method** (Waples & Do 2008): pairs of loci in the same population that are correlated (r²) in a predictable way under drift.

**IUCN/SSC Ne thresholds:**

| Ne | Status |
|---|---|
| < 50 | Critical — rapid allele fixation; high short-term extinction risk |
| 50 – 99 | High risk — genetic drift dominates |
| 100 – 499 | Concern — long-term viability at risk |
| ≥ 500 | Viable |

> Ne estimates from LD are sensitive to sample size (require n ≥ 10 per population) and linkage pruning. Always report confidence intervals.

---

## 7. Phase 4 — Adaptive Potential

Phase 4 requires environmental data — upload a per-population climate CSV or use auto-fetch in the sidebar.

### Partial RDA (GEA tab)

Redundancy Analysis (RDA) regresses population-level allele frequencies on environmental variables, with geographic position (latitude/longitude) partialled out. This separates **adaptive genetic variation** (driven by environment) from **neutral differentiation** (driven by drift and isolation).

Output:
- **Constrained variance (%)**: how much genetic variation is explained by the environment
- **Top candidate loci**: SNPs with the highest loadings on environmental gradients
- **Biplot**: populations and environmental gradients in ordination space

### LFMM2 (GEA tab)

LFMM2 tests each SNP individually for genotype-environment association, controlling for population structure via K latent factors. Results:

- **Manhattan plot**: −log₁₀(p-value) per SNP, per environmental variable
- **Significant loci**: after Benjamini-Hochberg FDR correction at 5%
- The Genomic Inflation Factor (GIF) calibrates p-values — values far from 1.0 indicate model misspecification (try different K)

### Genomic Offset (Offset tab)

Genomic offset (Fitzpatrick & Keller 2015) measures a population's **climate vulnerability**: the Euclidean distance in constrained ordination space between its current and projected future environmental position.

Larger offset → larger adaptive gap → more vulnerable to climate change.

> Genomic offset is an _exposure_ metric — it does not account for phenotypic plasticity, migration potential, or the rate of adaptation. Combine with demographic data for management decisions.

**Upload requirements:**
- Current environment CSV: one row per population, BIO1, BIO12, etc.
- Future environment CSV: same populations and variables under a climate scenario (e.g. SSP2-4.5 for 2070)

---

## 8. Phase 5 — Export & Reporting

### 8.1 HTML Conservation Report

The **Export tab** generates a self-contained HTML report summarising all analyses. It includes:

- Executive summary cards (populations, individuals, mean He, at-risk count)
- Risk overview with colour-coded badges
- Per-population statistics table
- Pairwise F_ST matrix with differentiation flags
- ADMIXTURE ancestry table (if Q-matrix loaded)
- Inbreeding summary (F-ROH or F coefficient)
- Ne estimates with IUCN risk flags
- GEA candidate loci summary (if run)
- Genomic offset table (if run)
- Automated conservation recommendations
- Full methods and citation list

Click **"Download HTML"** — save the file and open it in any browser. The report is self-contained (no internet required to view it) and can be printed to PDF directly from the browser.

> **Tip:** Before downloading, make sure all the analyses you want included (Ne, GEA, Offset) have been run in their respective tabs.

### 8.2 Statistics CSV

Click **"Download CSV"** in the Export tab (or the Statistics tab) for a flat file of all per-population metrics. Opens directly in Excel, R, or Python.

### 8.3 Map PNG

1. Navigate to the **Map tab**
2. Click the **camera icon** in the top-left corner of the map
3. The current viewport (including markers, FST edges, and ADMIXTURE pies) is saved as a PNG

The camera icon is added by the `leaflet.extras2` package. If it is missing, run:
```r
install.packages("leaflet.extras2")
```
then restart the app.

### 8.4 Chart PNG / SVG

Every Plotly chart (FST heatmap, PCA scatter, ADMIXTURE bar, Ne plot, GEA Manhattan, Offset bar) has a built-in download button:

1. Hover over the chart — the Plotly toolbar appears in the top-right corner
2. Click the **camera icon**
3. A PNG is saved to your Downloads folder

To change the format to SVG or to set a custom filename/size, use the Plotly modebar settings. In R, charts are configured with:
```r
plotly::config(
  p,
  toImageButtonOptions = list(
    format   = "svg",          # "png", "svg", "jpeg", "webp"
    filename = "my_chart",
    width    = 1200,
    height   = 800
  )
)
```

---

## 9. Interpreting Results

### Putting it all together

A complete genetic health assessment for a population requires looking at multiple metrics simultaneously:

| Scenario | Interpretation | Recommended action |
|---|---|---|
| Low He + high F + low Ne | Isolated, inbred population losing diversity fast | Genetic rescue — introduce individuals from a genetically distinct population |
| Low He + low F + moderate Ne | Historical bottleneck, random mating now | Monitor diversity; protect from further fragmentation |
| High F_ST between nearby populations | Landscape barrier to gene flow | Connectivity modelling; corridor restoration |
| High genomic offset + low Ne | Climate vulnerable AND genetically fragile | High priority for assisted migration or ex-situ insurance colony |
| High He + low F + high Ne | Healthy reference population | Maintain; potential donor for genetic rescue elsewhere |

### Common misconceptions

**"Low SNP He means the population is doomed."**
He from SNP data filtered at MAF ≥ 0.05 is systematically lower than microsatellite He. A SNP He of 0.12 may correspond to a microsatellite He of 0.50 or more. Always compare with a reference population from the same dataset, not across studies.

**"F > 0 always means inbreeding is a problem."**
Positive F can also arise from the Wahlund effect (mixing samples from multiple true sub-populations), null alleles in the sequencing, or population subdivision. Confirm with ROH-based F-ROH before drawing management conclusions.

**"High genomic offset means the population will go extinct under climate change."**
Genomic offset measures genetic maladaptation risk only. Populations with high plasticity, short generation times, or access to suitable microhabitats may survive even high offsets. It is one input among many.

---

## 10. Troubleshooting

### App won't start

```
Error in library(vcfR): there is no package called 'vcfR'
```
Run `source("install_packages.R")` to install all dependencies.

### VCF upload fails

- Check the VCF has biallelic SNPs only (filter with bcftools as shown in Section 3.1)
- Confirm the file is valid: `bcftools stats your_data.vcf.gz | head -30`
- Very large VCFs (> 500 MB) may time out — try subsetting to ≤ 50,000 SNPs

### Metadata CSV not matching

```
Error: Samples in VCF not found in metadata CSV: TIG_001, TIG_002
```
The `sample_id` column must exactly match the sample names in the VCF header. Check with:
```bash
bcftools query -l your_data.vcf.gz   # list VCF sample names
```
Then confirm those names appear in your CSV's `sample_id` column.

### Map shows no markers

- Confirm your metadata CSV has valid decimal-degree coordinates (not DMS)
- Latitudes must be between −90 and 90; longitudes between −180 and 180
- South latitudes and west longitudes must be negative (e.g. −3.45, −60.12)

### FST heatmap is all grey / zero

- Hudson F_ST requires at least 2 populations with overlapping loci
- Populations with very small n (< 3) are excluded from FST calculations
- Populations fixed for opposite alleles at every locus will show F_ST = 1.0 (not an error)

### Conservation Report download is empty or shows an error

Ensure `rmarkdown` is installed:
```r
install.packages(c("rmarkdown", "knitr"))
```
If the app is running on a server without internet access, `rmarkdown` may need pre-installed Pandoc. Check with `rmarkdown::pandoc_available()`.

### Map screenshot button missing

Install `leaflet.extras2` and restart the app:
```r
install.packages("leaflet.extras2")
shiny::runApp()
```

---

## 11. Glossary

**Allele frequency** — the proportion of a particular allele variant among all alleles at a locus in a population.

**ADMIXTURE** — a maximum-likelihood method that estimates the proportion of ancestry from K genetic clusters for each individual (Alexander et al. 2009).

**Effective population size (Ne)** — the size of an idealised Wright-Fisher population experiencing the same rate of genetic drift as the observed population. Always ≤ census size.

**Expected heterozygosity (He)** — the probability that two alleles randomly drawn from a population are different. Equivalent to gene diversity. He = 2p(1−p) averaged across loci.

**F_ST** — a measure of genetic differentiation between populations. Ranges 0 (identical allele frequencies) to 1 (completely fixed for different alleles).

**F-ROH** — inbreeding coefficient derived from the proportion of the genome in runs of homozygosity (ROH). Captures recent inbreeding from the last ~10 generations.

**Genomic offset** — RDA-based distance between a population's current and projected future environmental position in constrained ordination space; measures climate vulnerability.

**GEA (Genotype-Environment Association)** — statistical test linking genotype variation to environmental gradients, used to identify candidate loci under local adaptation.

**Inbreeding coefficient (F)** — (He − Ho) / He. Positive F indicates excess homozygosity relative to Hardy-Weinberg expectations.

**LFMM2** — Latent Factor Mixed Model, a mixed-effects regression that controls for population structure (K latent factors) when testing individual SNPs for GEA.

**Nucleotide diversity (π)** — average number of pairwise nucleotide differences per site between randomly chosen sequences. Closely related to He for diploid SNPs with the Nei (1978) sample-size correction.

**Observed heterozygosity (Ho)** — fraction of individuals in a population that carry two different alleles at a given locus, averaged across loci.

**Partial RDA** — Redundancy Analysis with geographic position partialled out; separates environment-driven adaptive variation from neutral structure.

**ROH** — Run of Homozygosity; a long stretch of consecutive homozygous SNP genotypes inherited identical-by-descent from a common ancestor.

**SNP** — Single Nucleotide Polymorphism; a single base-pair position where individuals carry different nucleotides.

**VCF** — Variant Call Format; the standard text file format for storing genotype data from sequencing pipelines.

**Wright-Fisher population** — a theoretical model of a constant-size, randomly-mating, non-overlapping-generation population used as the null baseline for drift and inbreeding calculations.

---

## 12. Citation & Licence

If you use PopGen Map in published research or management documents, please cite:

> CS (2026). *PopGen Map: an open-source browser-accessible platform for conservation genomics visualisation*. Version 0.6.0. [github.com/csunbird/population-genomics-map](https://github.com/csunbird/population-genomics-map)

**Licence:** MIT © 2026 CS — free to use, modify, and distribute with attribution.

**Contact:** seahcheyanne@gmail.com

**Repository:** [github.com/csunbird/population-genomics-map](https://github.com/csunbird/population-genomics-map)

---

_PopGen Map is provided as-is for research and conservation planning purposes. It is not a substitute for expert conservation genetics advice. Management decisions should integrate field surveys, demographic data, and specialist consultation._
