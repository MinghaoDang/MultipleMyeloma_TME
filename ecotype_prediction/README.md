# Reference-based ecotype prediction

This repository provides a reproducible pipeline to assign new samples to the reference immune microenvironment ecotypes.

## Recommended use

The primary classifier uses a sample-by-cell-state frequency matrix. Rows are samples and columns are reference cell states. The classifier aligns the input matrix to the reference feature space and returns:

- predicted reference ecotype
- class probabilities for E1-E5
- missing model features
- optional feature quality-control summaries

An optional DEG-based expression classifier is also included for datasets where matched cell-state frequencies are unavailable. The cell-frequency model is the recommended primary approach for single-cell datasets.

## Folder contents

- `scripts/predict_reference_ecotype.R`: assign ecotypes to a new sample-by-cell-state frequency matrix.
- `scripts/train_reference_ecotype_classifiers.R`: reproduce training/validation of the cell-frequency and DEG-based classifiers.
- `models/reference_ecotype_classifiers.rds`: trained classifier object.
- `data/reference/train.mtx.tsv`: reference sample-by-cell-state frequency matrix.
- `data/reference/train.group_label.tsv`: reference ecotype labels.
- `data/reference/W.tsv`: reference NMF W matrix used for cell-state feature selection.
- `data/reference/train_data.rda`: reference expression object used to reproduce the optional DEG-based classifier.
- `data/deconvolution/CIBERSORTx_reference_signature_matrix.txt`: CIBERSORTx-compatible reference signature matrix for estimating cell-state frequencies from bulk expression data.
- `features/cell_frequency_selected_features.csv`: selected cell states used by the primary classifier.
- `features/logTPM_DEG_selected_genes.csv`: selected genes used by the optional expression classifier.
- `examples/input/example_cell_frequency_matrix.tsv`: example input matrix.
- `examples/output/reference_ecotype_predictions.tsv`: example output.

## R dependencies

The prediction script requires:

```r
install.packages("glmnet")
```

The training/validation script additionally requires:

```r
install.packages(c("caret", "ggplot2", "patchwork"))
BiocManager::install(c("SummarizedExperiment", "limma"))
```

## Predict ecotypes for new samples

Input should be a tab-delimited file with `SampleID` as the first column and cell-state frequencies in the remaining columns:

```bash
Rscript scripts/predict_reference_ecotype.R \
  --input examples/input/example_cell_frequency_matrix.tsv \
  --classifier cell_frequency \
  --model models/reference_ecotype_classifiers.rds \
  --outdir results/prediction
```

For deconvolution-derived bulk frequency matrices, optional filters can be applied to remove sparse or highly concentrated cell states:

```bash
Rscript scripts/predict_reference_ecotype.R \
  --input path/to/new_cell_frequency_matrix.tsv \
  --classifier cell_frequency \
  --outdir results/my_dataset \
  --min_detected_fraction 0.10 \
  --max_sample_share 0.15 \
  --top5_sample_share 0.50
```

To use the optional DEG-based expression classifier, provide a gene-by-sample expression matrix with gene symbols in the first column:

```bash
Rscript scripts/predict_reference_ecotype.R \
  --input path/to/new_logTPM_expression_matrix.tsv \
  --classifier deg_expression \
  --outdir results/my_expression_dataset
```

Output files:

- `reference_ecotype_predictions.tsv`: predicted ecotype and probabilities.
- `reference_ecotype_prediction_counts.csv`: number of samples assigned to each ecotype.
- `input_feature_qc.csv`: feature-level detection and concentration metrics.
- `missing_model_features.txt`: selected model features absent from the input.
- `filtered_model_features.txt`: selected model features removed by optional filters.

## New scRNA-seq datasets

For new scRNA-seq datasets, users should first annotate cells by Seurat label transfer using the reference single-cell dataset. The reference data are available from Zenodo:

https://doi.org/10.5281/zenodo.20299059

This Zenodo record contains the complete reference dataset. Users should split the reference by major cell type and perform label transfer within each cell type to annotate the corresponding cell-state subtypes. After subtype annotation, users should generate a sample-by-cell-state frequency matrix and use `scripts/predict_reference_ecotype.R` to assign reference ecotypes.

## Bulk RNA-seq datasets

For bulk RNA-seq datasets, cell-state frequencies are not directly observed. We provide two possible options:

1. Estimate cell-state frequencies using deconvolution with `data/deconvolution/CIBERSORTx_reference_signature_matrix.txt`, then apply the primary cell-frequency classifier.
2. Apply the optional DEG-based expression classifier directly to expression data.

Because bulk classifier performance depends on the accuracy of deconvolution and cannot be directly evaluated without reference ecotype labels, we recommend trying both approaches when possible and reviewing class probabilities to identify ambiguous assignments.

## Reproduce training and validation

```bash
Rscript scripts/train_reference_ecotype_classifiers.R
```

Outputs are written to `results/train_validate/`.
