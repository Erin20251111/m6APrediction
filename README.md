---
editor_options: 
  markdown: 
    wrap: 72
---

**m6APrediction: Predict m6A RNA Modifications Using Random Forest
Models**

## Overview

`m6APrediction` is an R package for predicting m6A RNA modification
sites using a trained Random Forest classifier.\
It provides tools to:

-   encode DNA 5-mer sequences into per-position categorical features\
-   perform m6A probability prediction for both single observations and
    large batches\
-   return probability and binary m6A status (`Positive` / `Negative`)

This package is designed for the BIO215 practical sessions and serves as
a template for developing small machine-learning R packages.

The package exposes two main user functions:

-   **`prediction_multiple()`** — run predictions on a `data.frame`\
-   **`prediction_single()`** — convenient wrapper for predicting one
    observation

------------------------------------------------------------------------

## Installation

Install the development version directly from GitHub:

``` r
# If devtools is not installed:
# install.packages("devtools")

devtools::install_github("Erin20251111/m6APrediction")

# Alternatively, using remotes:
# install.packages("remotes")
# remotes::install_github("Erin20251111/m6APrediction")
```

------------------------------------------------------------------------

## Quick start

The package includes an example Random Forest model and input features
in `inst/extdata/`. You can run a full prediction workflow as follows:

``` r
library(m6APrediction)

# Load bundled model and example dataset
rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))

# ---- 1) Predict for multiple observations ----
out <- prediction_multiple(rf_fit, df, positive_threshold = 0.6)
head(out)

# ---- 2) Predict for a single observation ----
r1 <- df[1, ]

prediction_single(
  ml_fit = rf_fit,
  gc_content = r1$gc_content,
  RNA_type = as.character(r1$RNA_type),
  RNA_region = as.character(r1$RNA_region),
  exon_length = r1$exon_length,
  distance_to_junction = r1$distance_to_junction,
  evolutionary_conservation = r1$evolutionary_conservation,
  DNA_5mer = r1$DNA_5mer,
  positive_threshold = 0.5
)
```

------------------------------------------------------------------------

## Model Design

`prediction_multiple()` automatically:

-   checks input feature columns
-   encodes 5-mer sequences into per-nucleotide factor features
    (`nt_pos1`…`nt_pos5`)
-   standardizes factor levels for `RNA_type` and `RNA_region`
-   retrieves the `"Positive"` probability using
    `predict(..., type = "prob")`
-   assigns binary status based on a user-defined threshold

`prediction_single()` is a thin wrapper that converts one set of
features into a one-row data.frame and calls `prediction_multiple()`
internally.

------------------------------------------------------------------------

## Model performance (ROC / PRC)

To visualize classifier performance, include ROC and Precision–Recall
curves generated in Practical 4.

A common location for figures is:

![ROC Curve](man/figures/roc.png)

![PRC Curve](man/figures/prc.png)

and referenced as:

``` markdown
![ROC Curve](man/figures/roc.png)
![PRC Curve](man/figures/prc.png)
```

------------------------------------------------------------------------

## Development Notes

Use the following workflow when modifying the package:

``` r
devtools::document()  # regenerate Rd documentation
devtools::load_all()  # load the development version of the package
devtools::check()     # run R CMD check
```

------------------------------------------------------------------------

## Citation

If you use this package in academic or coursework settings, please cite:

m6APrediction (2025). GitHub: Erin20251111/m6APrediction.
