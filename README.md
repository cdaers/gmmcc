# gmmcc

**Gaussian Mixture Model for Clustering and Classification**

An R package that implements a hybrid classification model combining Logistic Regression (LR) and Gaussian Mixture Model (GMM) chains. Designed for robust, interpretable classification in both binary and multiclass settings.

---

## Installation

To install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("cdaers/gmmcc")
```

---

## Overview

`gmmcc` integrates statistical rigor and iterative refinement by combining LR-based feature selection with chained GMM classifiers. It is particularly suited for structured data where interpretability and performance are both desired.

**Key Features:**

- Logistic Regression + GMM chaining
- Feature selection via significance test or LASSO
- Flexible chain control: `num_chains`, `chain_length`
- Binary and multiclass support
- Automatic graphics output and misclassification analysis
- Built-in model stability test

---

## Quick Start

```r
library(gmmcc)

# Load demo dataset
data(brcancer_diag)

# Run hybrid LR-GMM model
model_result <- model_run(
  brcancer_diag,
  dat_nm = "brcancer_diag",
  lr_params = list(penalty = TRUE, min_features = 5),
  model_saved = TRUE,
  opt_metric = "AUC"
)

# Compare training and test accuracy
accuracy <- model_result$metric_rec$Accuracy
phase <- model_result$metric_rec$phase
boxplot(accuracy ~ phase, main = "Accuracy by Phase")
```

---

## Example: Multiclass Classification

```r
data(sale)

model_result <- model_run(
  sale,
  dat_nm = "sale",
  gmm_args = list(num_chains = 3),
  model_saved = TRUE,
  opt_metric = "AUC"
)

# Visualize accuracy by phase
accuracy <- model_result$metric_rec$Accuracy
phase <- model_result$metric_rec$phase
boxplot(accuracy ~ phase, main = "Sales Forecasting Accuracy")
```

---

## Vignette

A detailed walkthrough is available in the vignette:

```r
browseVignettes("gmmcc")
# or
vignette("gmmcc-intro")
```

---

## Function Reference

| Function            | Purpose                                                |
|---------------------|--------------------------------------------------------|
| `model_run()`       | Full training/testing pipeline                         |
| `model_train()`     | Train LR and/or GMM classifiers                        |
| `model_test()`      | Test on held-out data                                  |
| `gmm_classifier()`  | Core function for GMM-based error correction chains    |
| `stab_test_run()`   | Model stability testing over different train/test splits |

---

## Output

Each `model_run()` returns:

- `metric_rec`: Accuracy, Precision, Recall, F1, AUC by phase and chain
- `template_rec`: Feature templates selected in each chain
- Optional plots saved in the `graphics/` directory

---

## Citation

If you use this package, please cite:

> Hao Wang, Jue Hao. *gmmcc: Generalized Method of Moments for Clustering and Classification*. 2025.

```r
citation("gmmcc")
```

---

## License

MIT License © Hao Wang, Jue Hao

---

## Contact

Feel free to open an issue or contact the authors for questions and collaboration ideas.

