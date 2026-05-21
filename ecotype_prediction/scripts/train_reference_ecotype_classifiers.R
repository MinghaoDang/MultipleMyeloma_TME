suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(caret)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260519)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1] %||% "scripts/train_reference_ecotype_classifiers.R")
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_dir, "data", "reference"))) {
  repo_dir <- normalizePath(getwd(), mustWork = FALSE)
}
out_dir <- file.path(repo_dir, "results", "train_validate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_feature_names <- function(x) gsub("/", ".", x, fixed = TRUE)

read_frequency_matrix <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(x) <- x$SampleID
  x$SampleID <- NULL
  colnames(x) <- clean_feature_names(colnames(x))
  as.data.frame(lapply(x, as.numeric), check.names = FALSE, row.names = rownames(x))
}

read_w <- function(path) {
  w <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(w) <- clean_feature_names(w[[1]])
  w[[1]] <- NULL
  as.matrix(w)
}

macro_f1 <- function(obs, pred) {
  labs <- levels(obs)
  vals <- vapply(labs, function(lab) {
    tp <- sum(obs == lab & pred == lab)
    fp <- sum(obs != lab & pred == lab)
    fn <- sum(obs == lab & pred != lab)
    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    if (is.na(precision) || is.na(recall) || precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)
  }, numeric(1))
  mean(vals)
}

balanced_accuracy <- function(obs, pred) {
  labs <- levels(obs)
  vals <- vapply(labs, function(lab) {
    tp <- sum(obs == lab & pred == lab)
    fn <- sum(obs == lab & pred != lab)
    tn <- sum(obs != lab & pred != lab)
    fp <- sum(obs != lab & pred == lab)
    sens <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    spec <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
    mean(c(sens, spec), na.rm = TRUE)
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

metric_row <- function(model, split, obs, pred) {
  data.frame(model = model, split = split, n = length(obs), accuracy = mean(obs == pred), balanced_accuracy = balanced_accuracy(obs, pred), macro_f1 = macro_f1(obs, pred))
}

make_foldid <- function(y, k = 5) {
  folds <- createFolds(y, k = min(k, min(table(y))), returnTrain = FALSE)
  foldid <- integer(length(y))
  for (i in seq_along(folds)) foldid[folds[[i]]] <- i
  foldid
}

fit_tuned_glmnet <- function(x_train, y_train, alpha_grid = c(0, 0.5, 1)) {
  foldid <- make_foldid(y_train, k = 5)
  best <- NULL
  for (alpha in alpha_grid) {
    cv <- cv.glmnet(x_train, y_train, family = "multinomial", alpha = alpha, type.measure = "class", foldid = foldid, standardize = TRUE, maxit = 1e6)
    err <- min(cv$cvm)
    if (!is.finite(err)) next
    if (is.null(best) || err < best$cv_error) best <- list(alpha = alpha, lambda = cv$lambda.min, cv_error = err)
  }
  fit <- glmnet(x_train, y_train, family = "multinomial", alpha = best$alpha, lambda = best$lambda, standardize = TRUE, maxit = 1e6)
  list(fit = fit, alpha = best$alpha, lambda = best$lambda, cv_error = best$cv_error)
}

select_frequency_features <- function(w, dominant_fraction_cutoff = 0.5, relative_weight_cutoff = 0.3) {
  sorted <- t(apply(w, 1, sort, decreasing = TRUE))
  total <- rowSums(w)
  max_w <- sorted[, 1]
  second_w <- sorted[, 2]
  dominant <- colnames(w)[max.col(w, ties.method = "first")]
  spec <- data.frame(cell_type = rownames(w), dominant_ecotype = dominant, max_weight = max_w, second_weight = second_w, total_weight = total, dominant_fraction = ifelse(total > 0, max_w / total, 0), within_ecotype_relative_weight = NA_real_, stringsAsFactors = FALSE)
  for (ec in unique(spec$dominant_ecotype)) {
    idx <- spec$dominant_ecotype == ec
    spec$within_ecotype_relative_weight[idx] <- spec$max_weight[idx] / max(spec$max_weight[idx])
  }
  spec$selected <- spec$dominant_fraction >= dominant_fraction_cutoff & spec$within_ecotype_relative_weight >= relative_weight_cutoff
  spec
}

select_deg_union <- function(expr_train, y_train, top_n_per_group = 25, min_detected_frac = 0.1) {
  detected <- rowMeans(expr_train > 0) >= min_detected_frac
  x <- expr_train[detected, , drop = FALSE]
  y_train <- droplevels(y_train)
  design <- model.matrix(~ 0 + y_train)
  colnames(design) <- levels(y_train)
  fit <- eBayes(lmFit(x, design))
  selected <- character(0)
  tables <- list()
  for (grp in levels(y_train)) {
    rest <- setdiff(levels(y_train), grp)
    cont <- makeContrasts(contrasts = paste0(grp, "-(", paste(rest, collapse = "+"), ")/", length(rest)), levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cont))
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$gene_id <- rownames(tt)
    tt$ecotype_vs_rest <- grp
    selected <- union(selected, head(tt$gene_id, top_n_per_group))
    tables[[grp]] <- tt
  }
  list(selected = selected, tables = do.call(rbind, tables))
}

ref_dir <- file.path(repo_dir, "data", "reference")
freq_all <- read_frequency_matrix(file.path(ref_dir, "train.mtx.tsv"))
labels_df <- read.delim(file.path(ref_dir, "train.group_label.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
w <- read_w(file.path(ref_dir, "W.tsv"))

common_ids <- intersect(rownames(freq_all), labels_df$SampleID)
freq_all <- freq_all[common_ids, , drop = FALSE]
labels_df <- labels_df[match(common_ids, labels_df$SampleID), , drop = FALSE]
y <- factor(labels_df$Ecotype, levels = paste0("E", 1:5))
names(y) <- rownames(freq_all)
train_idx <- createDataPartition(y, p = 2/3, list = FALSE)[, 1]
validation_idx <- setdiff(seq_along(y), train_idx)

shared_features <- intersect(colnames(freq_all), rownames(w))
freq_all <- freq_all[, shared_features, drop = FALSE]
w <- w[shared_features, , drop = FALSE]
freq_spec <- select_frequency_features(w)
freq_features <- freq_spec$cell_type[freq_spec$selected]
x_freq <- as.matrix(freq_all[, freq_features, drop = FALSE])

freq_fit <- fit_tuned_glmnet(x_freq[train_idx, , drop = FALSE], droplevels(y[train_idx]))
freq_train_pred <- factor(as.character(predict(freq_fit$fit, x_freq[train_idx, , drop = FALSE], type = "class", s = freq_fit$lambda)), levels = levels(y))
freq_valid_pred <- factor(as.character(predict(freq_fit$fit, x_freq[validation_idx, , drop = FALSE], type = "class", s = freq_fit$lambda)), levels = levels(y))

load(file.path(ref_dir, "train_data.rda"))
expr <- assay(train_data, "logTPM")
gene_info <- as.data.frame(rowData(train_data))
expr <- expr[, names(y), drop = FALSE]
if ("protein_coding_nonMtRp" %in% colnames(gene_info)) {
  keep <- as.logical(gene_info$protein_coding_nonMtRp)
  keep[is.na(keep)] <- FALSE
  expr <- expr[keep, , drop = FALSE]
}
deg <- select_deg_union(expr[, train_idx, drop = FALSE], droplevels(y[train_idx]))
deg_genes <- intersect(deg$selected, rownames(expr))
x_expr_train <- t(expr[deg_genes, train_idx, drop = FALSE])
x_expr_valid <- t(expr[deg_genes, validation_idx, drop = FALSE])
expr_fit <- fit_tuned_glmnet(x_expr_train, droplevels(y[train_idx]))
expr_train_pred <- factor(as.character(predict(expr_fit$fit, x_expr_train, type = "class", s = expr_fit$lambda)), levels = levels(y))
expr_valid_pred <- factor(as.character(predict(expr_fit$fit, x_expr_valid, type = "class", s = expr_fit$lambda)), levels = levels(y))

metrics <- rbind(
  metric_row("cell_frequency_specific_contribution_glmnet", "train_2thirds_apparent", y[train_idx], freq_train_pred),
  metric_row("cell_frequency_specific_contribution_glmnet", "validation_1third", y[validation_idx], freq_valid_pred),
  metric_row("logTPM_DEG_glmnet_optional", "train_2thirds_apparent", y[train_idx], expr_train_pred),
  metric_row("logTPM_DEG_glmnet_optional", "validation_1third", y[validation_idx], expr_valid_pred)
)

write.csv(metrics, file.path(out_dir, "train_validation_performance.csv"), row.names = FALSE)
write.csv(freq_spec[freq_spec$selected, ], file.path(out_dir, "cell_frequency_selected_features.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = deg_genes), file.path(out_dir, "logTPM_DEG_selected_genes.csv"), row.names = FALSE)
saveRDS(list(primary_cell_frequency_model = list(model = freq_fit$fit, alpha = freq_fit$alpha, lambda = freq_fit$lambda, cv_error = freq_fit$cv_error, selected_features = freq_features), optional_logTPM_DEG_model = list(model = expr_fit$fit, alpha = expr_fit$alpha, lambda = expr_fit$lambda, cv_error = expr_fit$cv_error, selected_genes = deg_genes), class_levels = levels(y)), file.path(out_dir, "reference_ecotype_classifiers.rds"))

print(metrics)
