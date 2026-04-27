# popgen-map

**An open-source, browser-accessible platform for visualising population genomic data on interactive real-world maps — designed for conservation practitioners.**

> Phase 1 — Genetic Diversity Mapping  
> Part of the Population Genomics Map project

---

## What it does

Upload a VCF file and a sample coordinates CSV, and the app renders your populations on an interactive map coloured by genetic diversity metrics (He, Ho, π, F). Click any population marker to see its genetic health profile in plain language.

No bioinformatics expertise required to *use* it. No infrastructure required to *run* it — just R and a browser.

---

## Quick start

### 1. Install dependencies

```r
source("install_packages.R")
```

### 2. Launch the app

```r
shiny::runApp()
```

The app opens in your browser. Click **"Load demo dataset"** in the sidebar to explore immediately — no files needed.

---

## Input format

| File | Format | Required columns |
|------|--------|-----------------|
| Genotype data | VCF (`.vcf` or `.vcf.gz`) | Biallelic SNPs, any variant caller |
| Sample metadata | CSV | `sample_id`, `population`, `latitude`, `longitude` |

A template CSV is available via the **"Metadata CSV template"** download button in the app.

---

## Phase 1 metrics

| Metric | Description |
|--------|-------------|
| **He** | Expected heterozygosity — 2p(1−p) averaged across loci |
| **Ho** | Observed heterozygosity — proportion of heterozygous genotypes |
| **π** | Nucleotide diversity — sample-size corrected He |
| **F** | Inbreeding coefficient — (He − Ho) / He |

### Conservation thresholds (He)

| Range | Interpretation |
|-------|---------------|
| < 0.10 | 🔴 Critical — high extinction risk |
| 0.10 – 0.15 | 🟠 High risk — low diversity |
| 0.15 – 0.20 | 🟡 Moderate risk — monitoring recommended |
| 0.20 – 0.30 | 🟢 Lower risk — moderate diversity |
| > 0.30 | 🔵 Healthy — high diversity |

---

## Project structure

```
popgen-map/
├── app.R                      # Entry point
├── global.R                   # Packages, constants, module sourcing
├── ui.R                       # User interface
├── server.R                   # Server logic
├── install_packages.R         # One-time dependency installer
├── R/
│   ├── genomics.R             # VCF parsing, diversity statistics, demo data
│   ├── utils.R                # Metadata parsing, formatting helpers
│   ├── mod_upload.R           # File upload Shiny module
│   ├── mod_map.R              # Interactive Leaflet map module
│   └── mod_stats_panel.R      # Population statistics table module
├── data/
│   └── demo/
│       ├── demo_samples.csv   # Demo sample coordinates (committed)
│       └── generate_demo.R    # Script to regenerate demo VCF
├── www/
│   └── styles.css             # Custom CSS
└── tests/
    └── testthat/
        └── test_genomics.R    # Unit tests
```

---

## Demo dataset

The built-in demo simulates a conservation genetics survey for a fictional large felid species across five populations in a fragmented tropical landscape:

| Population | n | He (target) | Status |
|-----------|---|------------|--------|
| Central Highland | 20 | 0.32 | Source population |
| Northern Forest | 20 | 0.27 | Connected |
| River Delta | 18 | 0.21 | Partially isolated |
| Eastern Ridge | 15 | 0.14 | Fragmented |
| Southern Peatland | 12 | 0.09 | Isolated — bottleneck |

> ⚠️ Demo data is **simulated**. It is not real genotype data for any species.

To regenerate the demo VCF (e.g. after modifying `generate_demo.R`):

```r
source("data/demo/generate_demo.R")
```

---

## Running tests

```r
testthat::test_file("tests/testthat/test_genomics.R")
```

---

## Roadmap

| Phase | Features | Status |
|-------|----------|--------|
| **1** | Diversity mapping (He, Ho, π, F) on Leaflet map | ✅ Current |
| **2** | Population structure (ADMIXTURE pies, FST network) | 🔜 Planned |
| **3** | Inbreeding (ROH / F_ROH) and Ne per population | 🔜 Planned |
| **4** | Adaptive potential (RDA, LFMM2, climate vulnerability) | 🔜 Planned |
| **5** | Polish, multi-user testing, shinyapps.io deployment | 🔜 Planned |

---

## Deployment to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "your-account", token = "...", secret = "...")
rsconnect::deployApp()
```

---

## Contributing

Pull requests welcome. Please open an issue first to discuss significant changes.

Repository: [github.com/csunbird/population-genomics-map](https://github.com/csunbird/population-genomics-map)

---

## Licence

MIT © 2026 CS
