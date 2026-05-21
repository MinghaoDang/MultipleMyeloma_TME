suppressPackageStartupMessages({
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] == length(args)) return(default)
  args[idx[1] + 1]
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "scripts/predict_reference_ecotype.R"
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_dir, "models"))) repo_dir <- normalizePath(getwd(), mustWork = FALSE)

input_file <- get_arg("--input", file.path(repo_dir, "examples", "input", "example_cell_frequency_matrix.tsv"))
model_file <- get_arg("--model", file.path(repo_dir, "models", "reference_ecotype_classifiers.rds"))
out_dir <- get_arg("--outdir", file.path(repo_dir, "results", "prediction"))
classifier_type <- get_arg("--classifier", "cell_frequency")
min_detected_fraction <- as.numeric(get_arg("--min_detected_fraction", "0"))
max_sample_share <- as.numeric(get_arg("--max_sample_share", "1"))
top5_sample_share <- as.numeric(get_arg("--top5_sample_share", "1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!classifier_type %in% c("cell_frequency", "deg_expression")) {
  stop("--classifier must be either 'cell_frequency' or 'deg_expression'")
}

clean_feature_names <- function(x) gsub("/", ".", x, fixed = TRUE)

read_frequency_matrix <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("SampleID" %in% colnames(x))
  rownames(x) <- x$SampleID
  x$SampleID <- NULL
  colnames(x) <- clean_feature_names(colnames(x))
  as.data.frame(lapply(x, as.numeric), check.names = FALSE, row.names = rownames(x))
}

read_expression_matrix <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- intersect(colnames(x), c("gene_id", "gene", "Gene", "GeneSymbol", "genesymbols"))
  if (length(gene_col) > 0) {
    rownames(x) <- x[[gene_col[1]]]
    x[[gene_col[1]]] <- NULL
  } else {
    rownames(x) <- x[[1]]
    x[[1]] <- NULL
  }
  as.matrix(as.data.frame(lapply(x, as.numeric), check.names = FALSE, row.names = rownames(x)))
}

normalize_frequency_rows <- function(x) {
  rs <- rowSums(x, na.rm = TRUE)
  out <- x
  idx <- is.finite(rs) & rs > 0
  out[idx, ] <- sweep(x[idx, , drop = FALSE], 1, rs[idx], "/")
  out[!is.finite(as.matrix(out))] <- 0
  out
}

classifier <- readRDS(model_file)

if (classifier_type == "cell_frequency") {
  new_frequency <- read_frequency_matrix(input_file)
  selected_features <- classifier$primary_cell_frequency_model$selected_features

  feature_qc <- data.frame(
    cell_type = colnames(new_frequency),
    detected_fraction = colMeans(new_frequency > 0),
    max_sample_share = apply(new_frequency, 2, function(z) if (sum(z) > 0) max(z) / sum(z) else 1),
    top5_sample_share = apply(new_frequency, 2, function(z) if (sum(z) > 0) sum(sort(z, decreasing = TRUE)[seq_len(min(5, length(z)))]) / sum(z) else 1),
    stringsAsFactors = FALSE
  )
  feature_qc$kept <- feature_qc$detected_fraction >= min_detected_fraction &
    feature_qc$max_sample_share <= max_sample_share &
    feature_qc$top5_sample_share <= top5_sample_share

  x <- normalize_frequency_rows(new_frequency)
  missing_model_features <- setdiff(selected_features, colnames(x))
  filtered_model_features <- intersect(selected_features, feature_qc$cell_type[!feature_qc$kept])

  aligned <- matrix(0, nrow = nrow(x), ncol = length(selected_features), dimnames = list(rownames(x), selected_features))
  usable <- intersect(selected_features, colnames(x))
  usable <- setdiff(usable, filtered_model_features)
  aligned[, usable] <- as.matrix(x[, usable, drop = FALSE])

  fit <- classifier$primary_cell_frequency_model
  write.csv(feature_qc, file.path(out_dir, "input_feature_qc.csv"), row.names = FALSE)
  writeLines(missing_model_features, file.path(out_dir, "missing_model_features.txt"))
  writeLines(filtered_model_features, file.path(out_dir, "filtered_model_features.txt"))
} else {
  expression <- read_expression_matrix(input_file)
  selected_genes <- classifier$optional_logTPM_DEG_model$selected_genes
  missing_model_features <- setdiff(selected_genes, rownames(expression))

  aligned <- matrix(0, nrow = ncol(expression), ncol = length(selected_genes), dimnames = list(colnames(expression), selected_genes))
  usable <- intersect(selected_genes, rownames(expression))
  aligned[, usable] <- t(expression[usable, , drop = FALSE])

  fit <- classifier$optional_logTPM_DEG_model
  writeLines(missing_model_features, file.path(out_dir, "missing_model_genes.txt"))
}

pred <- as.character(predict(fit$model, aligned, type = "class", s = fit$lambda))
prob <- predict(fit$model, aligned, type = "response", s = fit$lambda)[, , 1]

predictions <- data.frame(
  SampleID = rownames(aligned),
  Classifier = classifier_type,
  Predicted_Ecotype = pred,
  Prediction_Probability = apply(prob, 1, max),
  as.data.frame(prob, check.names = FALSE),
  check.names = FALSE
)

write.table(predictions, file.path(out_dir, "reference_ecotype_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(as.data.frame(table(Predicted_Ecotype = factor(predictions$Predicted_Ecotype, levels = classifier$class_levels))), file.path(out_dir, "reference_ecotype_prediction_counts.csv"), row.names = FALSE)

cat("Wrote predictions to:", file.path(out_dir, "reference_ecotype_predictions.tsv"), "\n")
print(table(Predicted_Ecotype = factor(predictions$Predicted_Ecotype, levels = classifier$class_levels)))
