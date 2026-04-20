source('global_config.R')

library('broom')
library('emmeans')
library('logistf')
library('patchwork')
library('car')

###################### Figure 3A-C correlation of tumor burden, ITH, proliferating with cell type abundance ######################
data1<-readRDS('sample.level.info.rds')
data2<-readRDS('cell_type_level2.frequency.mtx.rds')

data<-left_join(data1,data2)

data<-data[,c('NewID','Diagnosis',
              'Hyperdiploidy', 't(4;14)', 't(6;14)', 't(11;14)', 't(14;16)', 't(14;20)',
              '1q gain/CKS1B Pos', '17p loss/p53 Del', 'monosomy 13','Ecotype','Tumor.pct','S_G2M.pct','Diversity.Score.PC',
              "CD14_Mono","CD14CD16_Mono",
              "CD16_Mono","CD4T",
              "CD8T","cDC",
              "gdT","Macrophage",
              "MAIT","Mature_B",
              "Neutrophil","NK",
              "pDC","Platelet",
              "Progenitors","Stroma")]


colnames(data)[4:11]<-c('t_4_14','t_6_14','t_11_14','t_14_16','t_14_20','x1q_gain','x17p_loss','mono_13')

data$t_MAF<-data$t_14_16
data$t_MAF[data$t_14_20=='Pos']<-'Pos'

data$t_CCND<-data$t_11_14
data$t_CCND[data$t_6_14=='Pos']<-'Pos'

data<-data[,c('NewID','Diagnosis','Hyperdiploidy','t_CCND','t_4_14','t_MAF','x1q_gain','x17p_loss','mono_13',
              'Ecotype','Tumor.pct','S_G2M.pct','Diversity.Score.PC',
              "CD14_Mono","CD14CD16_Mono",
              "CD16_Mono","CD4T",
              "CD8T","cDC",
              "gdT","Macrophage",
              "MAIT","Mature_B",
              "Neutrophil","NK",
              "pDC","Platelet",
              "Progenitors","Stroma")]

df<-subset(data,Diagnosis!='nBM')

df$Diagnosis <- factor(df$Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM"))

cell_types <- c(
  "CD14_Mono","CD14CD16_Mono","CD16_Mono","CD4T","CD8T",
  "cDC","gdT","Macrophage","MAIT","Mature_B","Neutrophil",
  "NK","pDC","Platelet","Progenitors","Stroma"
)

#### heatmap ####

# exclude RRMM
df_rmRRMM <- df %>%
  filter(Diagnosis %in% c("MGUS", "SMM", "NDMM")) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("MGUS", "SMM", "NDMM")))

# 2. helper functions

safe_spearman <- function(x, y) {
  sub <- data.frame(x = x, y = y) %>% drop_na()
  if (nrow(sub) < 5) {
    return(data.frame(estimate = NA_real_, p = NA_real_, n = nrow(sub)))
  }
  test <- cor.test(sub$x, sub$y, method = "spearman")
  data.frame(
    estimate = unname(test$estimate),   # rho
    p = test$p.value,
    n = nrow(sub)
  )
}

safe_lm_term <- function(df_sub, outcome, predictor, covar) {
  sub <- df_sub %>%
    dplyr::select(all_of(c(outcome, predictor, covar))) %>%
    drop_na()
  
  if (nrow(sub) < 5) {
    return(data.frame(estimate = NA_real_, p = NA_real_, n = nrow(sub)))
  }
  
  fit <- lm(as.formula(paste0("`", outcome, "` ~ `", predictor, "` + `", covar, "`")), data = sub)
  tt <- broom::tidy(fit) %>% filter(term == predictor)
  
  if (nrow(tt) == 0) {
    return(data.frame(estimate = NA_real_, p = NA_real_, n = nrow(sub)))
  }
  
  data.frame(
    estimate = tt$estimate,   # beta
    p = tt$p.value,
    n = nrow(sub)
  )
}


# 3. run all 5 analyses

results <- list()

for (ct in cell_types) {
  
  # 1) Tumor burden
  tmp1 <- safe_spearman(df_rmRRMM[[ct]], df_rmRRMM$Tumor.pct)
  results[[length(results) + 1]] <- data.frame(
    Analysis = "Tumor burden",
    Panel = "Raw",
    CellType = ct,
    estimate = tmp1$estimate,
    p = tmp1$p,
    n = tmp1$n
  )
  
  # 2) S_G2M raw
  tmp2 <- safe_spearman(df_rmRRMM[[ct]], df_rmRRMM$S_G2M.pct)
  results[[length(results) + 1]] <- data.frame(
    Analysis = "S_G2M (raw)",
    Panel = "Raw",
    CellType = ct,
    estimate = tmp2$estimate,
    p = tmp2$p,
    n = tmp2$n
  )
  
  # 3) S_G2M adjusted
  tmp3 <- safe_lm_term(
    df_sub = df_rmRRMM,
    outcome = ct,
    predictor = "S_G2M.pct",
    covar = "Tumor.pct"
  )
  results[[length(results) + 1]] <- data.frame(
    Analysis = "S_G2M (adjusted)",
    Panel = "Adjusted",
    CellType = ct,
    estimate = tmp3$estimate,
    p = tmp3$p,
    n = tmp3$n
  )
  
  # 4) Diversity raw
  tmp4 <- safe_spearman(df_rmRRMM[[ct]], df_rmRRMM$Diversity.Score.PC)
  results[[length(results) + 1]] <- data.frame(
    Analysis = "Diversity (raw)",
    Panel = "Raw",
    CellType = ct,
    estimate = tmp4$estimate,
    p = tmp4$p,
    n = tmp4$n
  )
  
  # 5) Diversity adjusted
  tmp5 <- safe_lm_term(
    df_sub = df_rmRRMM,
    outcome = ct,
    predictor = "Diversity.Score.PC",
    covar = "Tumor.pct"
  )
  results[[length(results) + 1]] <- data.frame(
    Analysis = "Diversity (adjusted)",
    Panel = "Adjusted",
    CellType = ct,
    estimate = tmp5$estimate,
    p = tmp5$p,
    n = tmp5$n
  )
}

res <- bind_rows(results)


# 4. FDR correction within each analysis

res <- res %>%
  group_by(Analysis) %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  ungroup()

# significance label
res <- res %>%
  mutate(
    sig = case_when(
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


# 5. order cell types by Tumor burden effect size

tumor_order <- res %>%
  filter(Analysis == "Tumor burden") %>%
  arrange(estimate) %>%
  pull(CellType)

res$CellType <- factor(res$CellType, levels = tumor_order)

# row order in each panel
res$Analysis <- factor(
  res$Analysis,
  levels = rev(c(
    "Tumor burden",
    "S_G2M (raw)",
    "Diversity (raw)",
    "S_G2M (adjusted)",
    "Diversity (adjusted)"
  ))
)

# split
res_raw <- res %>%
  filter(Panel == "Raw") %>%
  mutate(
    Analysis = factor(Analysis, levels = rev(c("Tumor burden", "S_G2M (raw)", "Diversity (raw)")))
  )

res_adj <- res %>%
  filter(Panel == "Adjusted") %>%
  mutate(
    Analysis = factor(Analysis, levels = rev(c("S_G2M (adjusted)", "Diversity (adjusted)")))
  )


# 6. choose color scale limits
#    use symmetric raw limits separately for rho and beta panels

raw_lim <- max(abs(res_raw$estimate), na.rm = TRUE)
adj_lim <- max(abs(res_adj$estimate), na.rm = TRUE)


# 7. plot raw panel

p_raw <- ggplot(res_raw, aes(x = CellType, y = Analysis, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sig), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-raw_lim, raw_lim),
    name = "rho"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  ggtitle("Raw")


# 8. plot adjusted panel

p_adj <- ggplot(res_adj, aes(x = CellType, y = Analysis, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sig), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-adj_lim, adj_lim),
    name = "beta"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  ggtitle("Adjusted for tumor burden")

# 9. combine
p_combined <- p_raw / p_adj + plot_layout(heights = c(3, 2))

ggsave("heatmap_raw_adjusted_split_panels_noRRMM.pdf", p_combined, width = 12, height = 6.5)


#### dot plot ####

# 2. run three analyses within each disease stage

results <- list()

for (ct in cell_types) {
  for (stage in levels(df$Diagnosis)) {
    
    sub <- df %>% filter(Diagnosis == stage)
    
    # -------- Tumor burden correlation --------
    sub1 <- sub %>% drop_na(Tumor.pct, all_of(ct))
    if (nrow(sub1) >= 5) {
      fit1 <- cor.test(sub1[[ct]], sub1$Tumor.pct, method = "spearman")
      results[[length(results) + 1]] <- data.frame(
        CellType = ct,
        Stage = stage,
        Feature = "Tumor",
        estimate = unname(fit1$estimate),   # rho
        p = fit1$p.value,
        n = nrow(sub1)
      )
    }
    
    # -------- Proliferation adjusted for tumor --------
    sub2 <- sub %>% drop_na(S_G2M.pct, Tumor.pct, all_of(ct))
    if (nrow(sub2) >= 5) {
      fit2 <- lm(sub2[[ct]] ~ S_G2M.pct + Tumor.pct, data = sub2)
      tmp2 <- tidy(fit2) %>% filter(term == "S_G2M.pct")
      results[[length(results) + 1]] <- data.frame(
        CellType = ct,
        Stage = stage,
        Feature = "Proliferation",
        estimate = tmp2$estimate,           # beta
        p = tmp2$p.value,
        n = nrow(sub2)
      )
    }
    
    # -------- Diversity adjusted for tumor --------
    sub3 <- sub %>% drop_na(Diversity.Score.PC, Tumor.pct, all_of(ct))
    if (nrow(sub3) >= 5) {
      fit3 <- lm(sub3[[ct]] ~ Diversity.Score.PC + Tumor.pct, data = sub3)
      tmp3 <- tidy(fit3) %>% filter(term == "Diversity.Score.PC")
      results[[length(results) + 1]] <- data.frame(
        CellType = ct,
        Stage = stage,
        Feature = "Diversity",
        estimate = tmp3$estimate,           # beta
        p = tmp3$p.value,
        n = nrow(sub3)
      )
    }
  }
}

res <- bind_rows(results)


# 3. multiple testing correction

res <- res %>%
  group_by(Feature, Stage) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()


# 4. make points visible
#   - standardize effect size within each Feature
#   - set minimum dot size

res <- res %>%
  group_by(Feature) %>%
  mutate(
    effect_z = as.numeric(scale(estimate)),
    effect_z = ifelse(is.na(effect_z), 0, effect_z)
  ) %>%
  ungroup() %>%
  mutate(
    dot_size = pmax(-log10(p_adj), 1.5),   # minimum visible size
    sig = case_when(
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Feature = factor(Feature, levels = c("Diversity", "Proliferation", "Tumor")),
    Stage = factor(Stage, levels = c("MGUS", "SMM", "NDMM", "RRMM"))
  )


plot_df <- res %>%
  mutate(
    Stage = factor(Stage, levels = rev(c("MGUS", "SMM", "NDMM", "RRMM"))),
    Feature = factor(Feature, levels = rev(c("Diversity", "Proliferation", "Tumor"))),
    logFDR = -log10(p_adj),
    logFDR_plot = pmax(logFDR, 1.2),
    sig = p_adj < 0.05
  ) %>%
  group_by(Feature) %>%
  mutate(
    effect_z = as.numeric(scale(estimate)),
    effect_z = ifelse(is.na(effect_z), 0, effect_z)
  ) %>%
  ungroup()

cell_order <- c("CD4T","gdT","NK","Platelet","Mature_B","cDC",
                "CD16_Mono","MAIT","CD8T","CD14_Mono","CD14CD16_Mono",
                "Neutrophil","pDC","Macrophage","Stroma","Progenitors")

plot_df$CellType <- factor(plot_df$CellType, levels = cell_order)


plot_df$effect_z[plot_df$effect_z>3]<-3

p_stage <- ggplot(plot_df, aes(x = CellType, y = Stage)) +
  
  geom_point(
    aes(size = logFDR, fill = effect_z),
    shape = 21,
    color = "grey60",
    stroke = 0.3
  ) +
  
  geom_point(
    data = plot_df %>% filter(sig == "TRUE"),
    aes(size = logFDR, fill = effect_z),
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  facet_wrap(~ Feature, ncol = 1)+
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  scale_size_continuous(
    name = "-log10(FDR)",
    range = c(2, 15),
    breaks = c(2.5, 7.5, 12.5)
  ) +
  #scale_size(range = c(2, 15)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("dotplot_threepanel_sig_outline.pdf", p_stage, width = 11.5, height = 7.2)


#### scatter plot ####

meta_cols <- c(
  "NewID","Diagnosis","Hyperdiploidy","t_CCND","t_4_14","t_MAF",
  "x1q_gain","x17p_loss","mono_13","Ecotype","Tumor.pct"
)

diag_order <- c('MGUS','SMM','NDMM','RRMM')

cell_types <- setdiff(colnames(df), meta_cols)
cell_types <- setdiff(cell_types, "Neutrophil")

cell_order <- c(
  "Progenitors",
  "Stroma",
  "Macrophage",
  "pDC",
  "NK",
  "gdT",
  "CD4T",
  "CD16_Mono",
  "CD14CD16_Mono",
  "CD8T",
  "Platelet",
  "MAIT",
  "CD14_Mono",
  "cDC",
  "Mature_B"
)

df_long <- df %>%
  pivot_longer(
    cols = all_of(cell_types),
    names_to = "CellType",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(Tumor.pct), !is.na(Abundance))

df_long$CellType <- factor(df_long$CellType, levels = cell_order)

df_long <- df_long %>%
  mutate(panel = paste(Diagnosis, CellType, sep = " | "))

panel_levels <- unlist(lapply(diag_order, function(dx) {
  paste(dx, cell_order, sep = " | ")
}))

df_long$panel <- factor(df_long$panel, levels = panel_levels)

cor_df <- df_long %>%
  group_by(Diagnosis, CellType, panel) %>%
  group_modify(~{
    dd <- .x
    if (nrow(dd) < 5) {
      return(data.frame(rho = NA_real_, p.value = NA_real_, n = nrow(dd)))
    }
    
    test <- cor.test(dd$Tumor.pct, dd$Abundance, method = "spearman")
    
    data.frame(
      rho = unname(test$estimate),
      p.value = test$p.value,
      n = nrow(dd)
    )
  }) %>%
  ungroup() %>%
  group_by(Diagnosis) %>%
  mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(sig_panel = !is.na(FDR) & FDR < 0.05)

df_long2 <- df_long %>%
  left_join(
    cor_df %>% dplyr::select(Diagnosis, CellType, panel, rho, p.value, FDR, sig_panel),
    by = c("Diagnosis", "CellType", "panel")
  )

label_df <- df_long2 %>%
  group_by(Diagnosis, CellType, panel, rho, FDR, sig_panel) %>%
  summarise(
    x = quantile(Tumor.pct, 0.95, na.rm = TRUE),
    y = quantile(Abundance, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(sig_panel) %>%
  mutate(
    label = paste0(
      "rho = ", sprintf("%.2f", rho),
      "\nFDR = ", format(FDR, digits = 2, scientific = TRUE)
    )
  )

label_df$panel <- factor(label_df$panel, levels = panel_levels)

df_sig <- df_long2 %>% filter(sig_panel)
df_ns  <- df_long2 %>% filter(!sig_panel | is.na(sig_panel))

stage_colors <- c(
  "MGUS" = "#1c75bc",
  "SMM"  = "#fbb040",
  "NDMM" = "#d01c8b",
  "RRMM" = "#8856a7"
)


p <- ggplot() +
  geom_point(
    data = df_ns,
    aes(x = Tumor.pct, y = Abundance),
    color = "grey75", size = 0.8, alpha = 0.8
  ) +
  geom_smooth(
    data = df_ns,
    aes(x = Tumor.pct, y = Abundance, group = 1),
    method = "lm", se = FALSE,
    color = "grey65", linewidth = 0.6
  ) +
  geom_point(
    data = df_sig,
    aes(x = Tumor.pct, y = Abundance, color = Diagnosis),
    size = 0.9, alpha = 0.9
  ) +
  geom_smooth(
    data = df_sig,
    aes(x = Tumor.pct, y = Abundance, color = Diagnosis, group = 1),
    method = "lm", se = FALSE,
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ panel,
    scales = "free",
    ncol = length(cell_order)
  ) +
  
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.3,
    hjust = 1, vjust = 1,
    fontface = "bold"
  ) +
  
  scale_color_manual(values = stage_colors) +
  
  labs(
    x = "Tumor burden (%)",
    y = "Cell type abundance"
  ) +
  
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    #strip.text = element_blank(),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.25, "lines")
  )


ggsave(
  "Supplementary_scatterplots_stage_by_celltype_sigPanelHighlight.pdf",
  p, width = 18, height = 5
)

###################### Figure 3D correlation of tumor burden with cell status abundance ######################
data1<-readRDS('sample.level.info.rds')
data2<-readRDS('cell_subset.frequency.mtx.rds')

data<-left_join(data1,data2)

data<-data[,c(2,6,15:20,22:24,32,34,40,42:120)]


colnames(data)[4:11]<-c('t_4_14','t_6_14','t_11_14','t_14_16','t_14_20','x1q_gain','x17p_loss','mono_13')

data$t_MAF<-data$t_14_16
data$t_MAF[data$t_14_20=='Pos']<-'Pos'

data$t_CCND<-data$t_11_14
data$t_CCND[data$t_6_14=='Pos']<-'Pos'

data<-data[,c('NewID','Diagnosis','Hyperdiploidy','t_CCND','t_4_14','t_MAF','x1q_gain','x17p_loss','mono_13',
              'Ecotype','Tumor.pct',colnames(data)[14:93])]

df<-subset(data,Diagnosis!='nBM')

df$Diagnosis <- factor(df$Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM"))

cell_types <- colnames(df)[14:91]

# 2. run three analyses within each disease stage

results <- list()

for (ct in cell_types) {
  for (stage in levels(df$Diagnosis)) {
    
    sub <- df %>% filter(Diagnosis == stage)
    
    # -------- Tumor burden correlation --------
    sub1 <- sub %>% drop_na(Tumor.pct, all_of(ct))
    if (nrow(sub1) >= 5) {
      fit1 <- cor.test(sub1[[ct]], sub1$Tumor.pct, method = "spearman")
      results[[length(results) + 1]] <- data.frame(
        CellType = ct,
        Stage = stage,
        Feature = "Tumor",
        estimate = unname(fit1$estimate),   # rho
        p = fit1$p.value,
        n = nrow(sub1)
      )
    }
    
    # -------- Proliferation adjusted for tumor --------
    sub2 <- sub %>% drop_na(S_G2M.pct, Tumor.pct, all_of(ct))
    if (nrow(sub2) >= 5) {
      fit2 <- lm(sub2[[ct]] ~ S_G2M.pct + Tumor.pct, data = sub2)
      tmp2 <- tidy(fit2) %>% filter(term == "S_G2M.pct")
      results[[length(results) + 1]] <- data.frame(
        CellType = ct,
        Stage = stage,
        Feature = "Proliferation",
        estimate = tmp2$estimate,           # beta
        p = tmp2$p.value,
        n = nrow(sub2)
      )
    }
    
    # -------- Diversity adjusted for tumor --------
    sub3 <- sub %>% drop_na(Diversity.Score.PC, Tumor.pct, all_of(ct))
    if (nrow(sub3) >= 5) {
      fit3 <- lm(sub3[[ct]] ~ Diversity.Score.PC + Tumor.pct, data = sub3)
      tmp3 <- tidy(fit3) %>% filter(term == "Diversity.Score.PC")
      results[[length(results) + 1]] <- data.frame(
        CellType = ct,
        Stage = stage,
        Feature = "Diversity",
        estimate = tmp3$estimate,           # beta
        p = tmp3$p.value,
        n = nrow(sub3)
      )
    }
  }
}

res <- bind_rows(results)

# 3. multiple testing correction

res <- res %>%
  group_by(Feature, Stage) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()


# 4. make points visible
#   - standardize effect size within each Feature
#   - set minimum dot size

res <- res %>%
  group_by(Feature) %>%
  mutate(
    effect_z = as.numeric(scale(estimate)),
    effect_z = ifelse(is.na(effect_z), 0, effect_z)
  ) %>%
  ungroup() %>%
  mutate(
    dot_size = pmax(-log10(p_adj), 1.5),   # minimum visible size
    sig = case_when(
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Feature = factor(Feature, levels = c("Diversity", "Proliferation", "Tumor")),
    Stage = factor(Stage, levels = c("MGUS", "SMM", "NDMM", "RRMM"))
  )

# optional: order cell types by strongest overall signal
cell_order <- res %>%
  group_by(CellType) %>%
  summarise(max_signal = max(abs(effect_z), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_signal) %>%
  pull(CellType)

res$CellType <- factor(res$CellType, levels = cell_order)

df_plot <- res %>%
  mutate(
    Stage = factor(Stage, levels = rev(c("MGUS", "SMM", "NDMM", "RRMM"))),
    Feature = factor(Feature, levels = rev(c("Diversity", "Proliferation", "Tumor"))),
    logFDR = -log10(p_adj),
    logFDR_plot = pmax(logFDR, 1.2),
    sig = p_adj < 0.05
  ) %>%
  group_by(Feature) %>%
  mutate(
    effect_z = as.numeric(scale(estimate)),
    effect_z = ifelse(is.na(effect_z), 0, effect_z)
  ) %>%
  ungroup()


meta<-readRDS('MM.TME.235samples.TME.meta.data.rds')
aa<-unique(meta[,c('cell.type.level2','cell.status')])

df_plot<-left_join(df_plot,aa,by=c('CellType'='cell.status'))
df_plot<-subset(df_plot,cell.type.level2!='Neutrophil')

df_plot$cell.type.level2<-as.character(df_plot$cell.type.level2)

df_plot$cell.type.level2[df_plot$cell.type.level2%in%c('CD14_Mono','CD14CD16_Mono','CD16_Mono')]<-'Monocytes'
df_plot$cell.type.level2[df_plot$cell.type.level2%in%c('Platelet','Stroma')]<-'Other'
df_plot$cell.type.level2[df_plot$cell.type.level2%in%c('gdT','MAIT')]<-'OtherT'
df_plot$cell.type.level2[df_plot$cell.type.level2%in%c('cDC','pDC')]<-'DC'

df_plot$cell.type.level2<-factor(df_plot$cell.type.level2,levels=c("CD4T","OtherT","NK","CD8T",
                                                                   "Mature_B","Monocytes","Macrophage",
                                                                   "DC","Other","Progenitors"))
df_plot$effect_z[df_plot$effect_z>3]<-3
df_plot$effect_z[df_plot$effect_z< -3]<- -3

p_stage <- ggplot(df_plot, aes(x = CellType, y = Stage)) +
  
  geom_point(
    aes(size = logFDR, fill = effect_z),
    shape = 21,
    color = "grey60",
    stroke = 0.3
  ) +
  
  geom_point(
    data = df_plot %>% filter(sig == "TRUE"),
    aes(size = logFDR, fill = effect_z),
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  facet_grid(rows=vars(Feature),cols=vars(cell.type.level2), scales = "free_x", space = "free_x")+
  #facet_wrap(~ Feature, ncol = 1)+
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  scale_size_continuous(
    name = "-log10(FDR)",
    range = c(2, 10),
    breaks = c(2.5, 7.5, 12.5)
  ) +
  #scale_size(range = c(2, 15)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90)
  )

ggsave("cell_subset_dotplot_threepanel_sig_outline.pdf", p_stage, width = 25, height = 6)


###################### Figure 3E RRMM prior CART TME cell type level2 compare ######################
meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  select(NewID, cell.type.level2, Diagnosis) %>%
  add_count(NewID, name = "n_cells") %>%
  filter(
    n_cells >= 200,
    Diagnosis %in% c("NDMM", "RRMM"),
    NewID != "RRMM01"
  ) %>%
  select(-n_cells)

group_info <- readRDS("RRMM.treatment_history.rds") %>%
  select(NewID, `Prior CAR-T`, `Prior lines`) %>%
  rename(
    group = `Prior CAR-T`,
    Prior_lines = `Prior lines`
  )

data <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup() %>%
  left_join(group_info, by = "NewID") %>%
  mutate(
    Prior_lines = coalesce(Prior_lines, 0),
    Prior_lines = if_else(Prior_lines > 3, ">3", "<=3"),
    group = coalesce(group, "NDMM"),
    group = factor(group, levels = c("NDMM", "No", "Yes"))
  ) %>%
  as.data.frame()

library('ggsignif')

plotx<-signif_plot(data,
                   x='group',
                   y='percentage',
                   split='cell.type.level2',
                   signif.cutoff=0.05,
                   ncol=8,
                   fill_by='group',
                   fill_color=c("NDMM" = "#d01c8b",'Yes'='#a6d96a','No'='#fee08b'),
                   color_by='Prior_lines',
                   color=c('<=3'='grey60','>3'='grey30'),
                   ylab='% (of TME cells)')
ggsave('TME.cell.type.level2.compare.recent.CART.pdf',plotx,width=24,height=12)


###################### Figure 3F RRMM anti-CD38 TME cell type level2 compare ######################
meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  select(NewID, cell.type.level2, Diagnosis) %>%
  add_count(NewID, name = "n_cells") %>%
  filter(
    n_cells >= 200,
    Diagnosis %in% c("NDMM", "RRMM"),
    NewID != "RRMM01"
  ) %>%
  select(-n_cells)

group_info <- readRDS("RRMM.treatment_history.rds") %>%
  select(NewID, `time interval between sample collection and most recent anti-CD38`, `Prior lines`, `Recent exposure group`) %>%
  rename(
    group = `time interval between sample collection and most recent anti-CD38`,
    Prior_lines = `Prior lines`,
    recent = `Recent exposure group`
  )

data <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  tidyr::complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup() %>%
  left_join(group_info, by = "NewID") %>%
  mutate(
    Prior_lines = coalesce(Prior_lines, 0),
    Prior_lines = if_else(Prior_lines > 3, ">3", "<=3"),
    group = coalesce(group, "NDMM"),
    group = case_when(
      recent == "anti-CD38" ~ "Recent_exposure",
      Diagnosis == "RRMM" & group == "NDMM" ~ "No-anti-CD38",
      TRUE ~ group
    ),
    group = factor(group, levels = c("NDMM", "No-anti-CD38", "Recent_exposure", "<6mo", ">6mo"))
  ) %>%
  as.data.frame()

plotx<-signif_plot(data,
                   x='group',
                   y='percentage',
                   split='cell.type.level2',
                   signif.cutoff=0.05,
                   ncol=8,
                   fill_by='group',
                   color_by='Prior_lines',
                   color=c('<=3'='grey60','>3'='grey30'),
                   fill_color=c("NDMM" = "#d01c8b",'Recent_exposure'='#1a9850','<6mo'='#a6d96a','>6mo'='#fee08b','No-anti-CD38'='#f46d43'),
                   ylab='% (of TME cells)')
ggsave('TME.cell.type.level2.compare.anti_CD38_exposure.pdf',plotx,width=24,height=12)

###################### Figure 3G RRMM anti-CD38 TME lymphoid cell subset compare ######################
meta <- readRDS("MM.TME.235samples.TME.meta.data.rds")
meta = meta[,c('NewID','cell.type.level2','cell.status','Diagnosis')]

group_info <- readRDS("RRMM.treatment_history.rds") %>%
  select(NewID, `time interval between sample collection and most recent anti-CD38`, `Prior lines`, `Recent exposure group`) %>%
  rename(
    group = `time interval between sample collection and most recent anti-CD38`,
    Prior_lines = `Prior lines`,
    recent = `Recent exposure group`
  )

for(celltype in c('CD8T','CD4T','NK','Mature_B')){ # 'Macrophage','Neutrophil','DC','CD16_Mono','CD14_Mono',
  
  meta.use <- meta %>% 
    filter(cell.type.level2==celltype) %>%
    select(NewID, cell.status, Diagnosis) %>%
    add_count(NewID, name = "n_cells") %>%
    filter(
      n_cells > 50,
      Diagnosis %in% c("NDMM", "RRMM"),
      NewID != "RRMM01"
    ) %>%
    select(-n_cells)
  
  data <- meta.use %>%
    count(NewID, cell.status, name = "Freq") %>%
    tidyr::complete(NewID, cell.status, fill = list(Freq = 0)) %>%
    left_join(meta.use %>% distinct(NewID, Diagnosis), by = "NewID") %>%
    group_by(NewID) %>%
    mutate(
      count = sum(Freq),
      percentage = 100 * Freq / count
    ) %>%
    ungroup() %>%
    left_join(group_info, by = "NewID") %>%
    mutate(
      Prior_lines = coalesce(Prior_lines, 0),
      Prior_lines = if_else(Prior_lines > 3, ">3", "<=3"),
      group = coalesce(group, "NDMM"),
      group = case_when(
        recent == "anti-CD38" ~ "Recent_exposure",
        Diagnosis == "RRMM" & group == "NDMM" ~ "No-anti-CD38",
        TRUE ~ group
      ),
      group = factor(group, levels = c("NDMM", "No-anti-CD38", "Recent_exposure", "<6mo", ">6mo"))
    ) %>%
    as.data.frame()
  
  plotx<-signif_plot(data,
                     x='group',
                     y='percentage',
                     split='cell.status',
                     signif.cutoff=0.05,
                     ncol=11,
                     fill_by='group',
                     color_by='Prior_lines',
                     color=c('<=3'='grey60','>3'='grey30'),
                     fill_color=c("NDMM" = "#d01c8b",'Recent_exposure'='#1a9850','<6mo'='#a6d96a','>6mo'='#fee08b','No-anti-CD38'='#f46d43'),
                     ylab=paste0('% (of ',celltype,' cells)'))
  
  ggsave(paste0(celltype,'.cell.status.compare.anti_CD38_exposure.pdf'),plotx,width=30,height=5)
  
}

###################### Figure 3G RRMM anti-CD38 TME myeloid cell subset compare ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  filter(cell.type.level2%in%c('Macrophage','Neutrophil','cDC','pDC','CD16_Mono','CD14_Mono','CD14CD16_Mono'))  %>%
  select(NewID, cell.status, Diagnosis) %>%
  add_count(NewID, name = "n_cells") %>%
  filter(
    n_cells > 100,
    Diagnosis %in% c("NDMM", "RRMM"),
    NewID != "RRMM01"
  ) %>%
  select(-n_cells)

group_info <- readRDS("RRMM.treatment_history.rds") %>%
  select(NewID, `time interval between sample collection and most recent anti-CD38`, `Prior lines`, `Recent exposure group`) %>%
  rename(
    group = `time interval between sample collection and most recent anti-CD38`,
    Prior_lines = `Prior lines`,
    recent = `Recent exposure group`
  )

data <- meta %>%
  count(NewID, cell.status, name = "Freq") %>%
  tidyr::complete(NewID, cell.status, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup() %>%
  left_join(group_info, by = "NewID") %>%
  mutate(
    Prior_lines = coalesce(Prior_lines, 0),
    Prior_lines = if_else(Prior_lines > 3, ">3", "<=3"),
    group = coalesce(group, "NDMM"),
    group = case_when(
      recent == "anti-CD38" ~ "Recent_exposure",
      Diagnosis == "RRMM" & group == "NDMM" ~ "No-anti-CD38",
      TRUE ~ group
    ),
    group = factor(group, levels = c("NDMM", "No-anti-CD38", "Recent_exposure", "<6mo", ">6mo"))
  ) %>%
  as.data.frame()

plotx<-signif_plot(data,
                   x='group',
                   y='percentage',
                   split='cell.status',
                   signif.cutoff=0.05,
                   ncol=14,
                   fill_by='group',
                   color_by='Prior_lines',
                   color=c('<=3'='grey60','>3'='grey30'),
                   fill_color=c("NDMM" = "#d01c8b",'Recent_exposure'='#1a9850','<6mo'='#a6d96a','>6mo'='#fee08b','No-anti-CD38'='#f46d43'),
                   ylab=paste0('% (of ','Myeloid',' cells)'))

ggsave(paste0('Myeloid','.cell.status.compare.anti_CD38_exposure.pdf'),plotx,width=35,height=10)


###################### Figure S5A TME cell type level2 compare between primary cytogenetics ######################

cyto_level <- c("HY", "t_CCND", "t_MAF", "t(4;14)")

clinical_info <- readRDS("sample.level.info.rds") %>%
  as.data.frame() %>%
  dplyr::select(NewID, Primary_Cytogenetics) %>%
  mutate(
    Primary_Cytogenetics = dplyr::recode(
      Primary_Cytogenetics,
      "HY;t(11;14)" = "t(11;14)",
      "HY;t(4;14)"  = "t(4;14)",
      "HY;t(14;16)" = "t(14;16)"
    )
  )

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  dplyr::select(NewID, cell.type.level2, Diagnosis) %>%
  as.data.frame() %>%
  left_join(clinical_info, by = "NewID") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(
    n_cells >= 200,
    Primary_Cytogenetics %in% c("HY", "t(11;14)", "t(14;16)", "t(14;20)", "t(4;14)", "t(6;14)")
  ) %>%
  mutate(
    Primary_Cytogenetics = case_when(
      Primary_Cytogenetics %in% c("t(11;14)", "t(6;14)")   ~ "t_CCND",
      Primary_Cytogenetics %in% c("t(14;16)", "t(14;20)") ~ "t_MAF",
      TRUE ~ Primary_Cytogenetics
    ),
    Primary_Cytogenetics = factor(Primary_Cytogenetics, levels = cyto_level)
  ) %>%
  select(-n_cells)

data <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup() %>%
  left_join(
    meta %>% distinct(NewID, Diagnosis, Primary_Cytogenetics),
    by = "NewID"
  )

P<-ggplot(data, aes(x=Primary_Cytogenetics, y=percentage) ) + 
  geom_boxplot(alpha = 1, show.legend = T,width = 0.75, aes(fill=Primary_Cytogenetics), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=2, aes(color=Diagnosis)) +
  scale_color_manual(values=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  scale_fill_manual(values=c('t(4;14)'='#8dd3c7','t_CCND'='#bebada','t_MAF'='#80b1d3','HY'='#ff7f00'))+
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  facet_wrap(.~cell.type.level2,ncol=8,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)") 
ggsave('TME.compare.by.cytogenetics.pdf',P,width=25,height=12)


kruskal_results <- data %>%
  group_by(cell.type.level2) %>%
  summarise(p_value = kruskal.test(percentage ~ Primary_Cytogenetics)$p.value)


###################### Figure S5B correlation of cytogenetics with cell type abundance heatmap ######################
meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  filter(Diagnosis!='nBM') %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c('MGUS','SMM','NDMM','RRMM'))) %>%
  select(NewID, cell.type.level2, Diagnosis) %>%
  as.data.frame()

clinical.info <- readRDS("sample.level.info.rds") %>%
  as.data.frame() %>%
  dplyr::select(
    NewID,
    Hyperdiploidy,
    `t(4;14)`,
    `t(6;14)`,
    `t(11;14)`,
    `t(14;16)`,
    `t(14;20)`,
    `1q gain/CKS1B Pos`,
    `17p loss/p53 Del`,
    `monosomy 13`
  ) %>%
  dplyr::rename(
    t_4_14    = `t(4;14)`,
    t_6_14    = `t(6;14)`,
    t_11_14   = `t(11;14)`,
    t_14_16   = `t(14;16)`,
    t_14_20   = `t(14;20)`,
    x1q_gain  = `1q gain/CKS1B Pos`,
    x17p_loss = `17p loss/p53 Del`,
    mono_13   = `monosomy 13`
  ) %>%
  mutate(
    t_MAF  = if_else(t_14_20 == "Pos", "Pos", t_14_16),
    t_CCND = if_else(t_6_14 == "Pos", "Pos", t_11_14)
  ) %>%
  dplyr::select(
    NewID,
    Hyperdiploidy,
    t_CCND,
    t_4_14,
    t_MAF,
    x1q_gain,
    x17p_loss,
    mono_13
  )

meta<-left_join(meta,clinical.info,by=c('NewID'='NewID'))

data <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  tidyr::complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = Freq / count
  ) %>%
  ungroup() %>%
  left_join(
    meta %>%
      distinct(
        NewID, Diagnosis, Hyperdiploidy,
        t_CCND, t_4_14, t_MAF,
        x1q_gain, x17p_loss, mono_13
      ),
    by = "NewID"
  )

library(limma)
library(car)

features <- c(
  "Hyperdiploidy", "t_CCND", "t_4_14", "t_MAF",
  "x1q_gain", "x17p_loss", "mono_13"
)

dat <- data %>%
  filter(
    count >= 200,
    !is.na(percentage),
    !is.na(Diagnosis)
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis),
    percentage2 = pmin(pmax(percentage, 0.025), 0.975),
    logit_prop = car::logit(percentage2)
  )

run_one_feature <- function(feature) {
  dat %>%
    filter(.data[[feature]] %in% c("Pos", "Neg")) %>%
    mutate(group = factor(.data[[feature]], levels = c("Neg", "Pos"))) %>%
    group_by(cell.type.level2) %>%
    group_modify(~{
      fit <- lm(logit_prop ~ group + Diagnosis, data = .x)
      broom::tidy(fit) %>%
        filter(term == "groupPos") %>%
        mutate(
          feature = feature,
          n = nrow(.x),
          n_pos = sum(.x$group == "Pos"),
          n_neg = sum(.x$group == "Neg")
        )
    }) %>%
    ungroup()
}

results <- map_dfr(features, run_one_feature) %>%
  mutate(
    FDR_all = p.adjust(p.value, method = "BH")
  ) %>%
  group_by(feature) %>%
  mutate(
    FDR_by_feature = p.adjust(p.value, method = "BH"),
    odds_ratio = exp(estimate),
    direction = ifelse(estimate > 0, "higher in Pos", "lower in Pos")
  ) %>%
  ungroup()


heatmap_df <- results %>%
  mutate(
    sig = case_when(
      FDR_by_feature < 0.001 ~ "***",
      FDR_by_feature < 0.01  ~ "**",
      FDR_by_feature < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

heatmap_df$feature<-factor(heatmap_df$feature,levels=rev(c("Hyperdiploidy","t_CCND","t_4_14","t_MAF","x1q_gain","x17p_loss","mono_13")))

plotx<-ggplot(heatmap_df, aes(y = feature, x = cell.type.level2, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sig), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "Log-odds\nratio"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave('TME.compare.by.cytogenetics.heatmap.pdf',plotx,width=10,height=5)


###################### Figure S5C correlation of cytogenetics with cell status abundance ######################

df <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  filter(
    Diagnosis != "nBM",
    !cell.type.level2 %in% c("gdT", "MAIT", "Platelet", "Stroma")
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM")),
    cell.type.level2 = case_when(
      cell.type.level2 %in% c(
        "CD14_Mono", "CD14CD16_Mono", "CD16_Mono",
        "Macrophage", "cDC", "Neutrophil", "pDC"
      ) ~ "Myeloid",
      TRUE ~ cell.type.level2
    ),
    cell.type.level2 = factor(
      cell.type.level2,
      levels = c("CD4T", "CD8T", "NK", "Mature_B", "Myeloid", "Progenitors")
    )
  ) %>%
  select(NewID, cell.type.level2, cell.status, Diagnosis) %>%
  count(NewID, cell.type.level2, cell.status, Diagnosis, name = "count") %>%
  group_by(NewID, cell.type.level2, Diagnosis) %>%
  complete(cell.status, fill = list(count = 0)) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup()


clinical.info <- readRDS("sample.level.info.rds") %>%
  as.data.frame() %>%
  dplyr::select(
    NewID,
    Hyperdiploidy,
    `t(4;14)`,
    `t(6;14)`,
    `t(11;14)`,
    `t(14;16)`,
    `t(14;20)`,
    `1q gain/CKS1B Pos`,
    `17p loss/p53 Del`,
    `monosomy 13`
  ) %>%
  dplyr::rename(
    t_4_14    = `t(4;14)`,
    t_6_14    = `t(6;14)`,
    t_11_14   = `t(11;14)`,
    t_14_16   = `t(14;16)`,
    t_14_20   = `t(14;20)`,
    x1q_gain  = `1q gain/CKS1B Pos`,
    x17p_loss = `17p loss/p53 Del`,
    mono_13   = `monosomy 13`
  ) %>%
  mutate(
    t_MAF  = if_else(t_14_20 == "Pos", "Pos", t_14_16),
    t_CCND = if_else(t_6_14 == "Pos", "Pos", t_11_14)
  ) %>%
  dplyr::select(
    NewID,
    Hyperdiploidy,
    t_CCND,
    t_4_14,
    t_MAF,
    x1q_gain,
    x17p_loss,
    mono_13
  )


df<-left_join(df,clinical.info,by=c('NewID'='NewID'))


features <- c(
  "Hyperdiploidy", "t_CCND", "t_4_14",
  "t_MAF", "x1q_gain", "x17p_loss", "mono_13"
)


df2 <- df %>%
  filter(
    total >= 50,
    !is.na(proportion),
    !is.na(Diagnosis)
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis),
    proportion2 = pmin(pmax(proportion, 0.025), 0.975),
    logit_prop = car::logit(proportion2)
  )

run_model <- function(feature) {
  
  df2 %>%
    filter(.data[[feature]] %in% c("Pos", "Neg")) %>%
    mutate(group = factor(.data[[feature]], levels = c("Neg", "Pos"))) %>%
    group_by(cell.type.level2, cell.status) %>%
    group_modify(~{
      dd <- .x
      
      n_sample <- n_distinct(dd$NewID)
      n_pos <- n_distinct(dd$NewID[dd$group == "Pos"])
      n_neg <- n_distinct(dd$NewID[dd$group == "Neg"])
      
      ## QC 过滤
      if (n_sample < 10) return(data.frame())
      if (n_pos < 3 || n_neg < 3) return(data.frame())
      
      ## Diagnosis 有 >=2 个水平时纳入协变量，否则退化为 group-only model
      fit <- tryCatch(
        {
          if (n_distinct(dd$Diagnosis) >= 2) {
            lm(logit_prop ~ group + Diagnosis, data = dd)
          } else {
            lm(logit_prop ~ group, data = dd)
          }
        },
        error = function(e) NULL
      )
      
      if (is.null(fit)) return(data.frame())
      
      out <- broom::tidy(fit) %>%
        filter(term == "groupPos")
      
      if (nrow(out) == 0) return(data.frame())
      
      out %>%
        mutate(
          feature = feature,
          n_sample = n_sample,
          n_pos = n_pos,
          n_neg = n_neg
        )
    }) %>%
    ungroup()
}

results <- map_dfr(features, run_model)

results <- results %>%
  group_by(feature) %>%
  mutate(
    FDR = p.adjust(p.value, method = "BH"),
    odds_ratio = exp(estimate),
    direction = case_when(
      estimate > 0 ~ "higher in Pos",
      estimate < 0 ~ "lower in Pos",
      TRUE ~ "no change"
    )
  ) %>%
  ungroup()

feature_order <- c(
  "mono_13", "x17p_loss", "x1q_gain",
  "t_MAF", "t_4_14", "t_CCND", "Hyperdiploidy"
)

plot_df <- results %>%
  mutate(
    feature = factor(feature, levels = feature_order),
    sig = case_when(
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    estimate_plot = pmin(pmax(estimate, -0.4), 0.4)
  )

status_order <- plot_df %>%
  distinct(cell.type.level2, cell.status) %>%
  arrange(cell.type.level2, cell.status) %>%
  pull(cell.status)

plot_df <- plot_df %>%
  mutate(cell.status = factor(cell.status, levels = status_order))


p <- ggplot(plot_df, aes(x = cell.status, y = feature, fill = estimate_plot)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    data = plot_df %>% filter(sig != ""),
    aes(label = sig),
    size = 4,
    fontface = "bold"
  ) +
  facet_grid(. ~ cell.type.level2, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B73D9",
    mid = "#F2F2F2",
    high = "#F4A090",
    midpoint = 0,
    limits = c(-0.4, 0.4),
    name = "Log-odds\nratio"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing.x = unit(0.4, "lines")
  )

ggsave("cytogenetic_cellstate_heatmap_adjusted.pdf", p, width = 20, height = 5)

###################### Figure S5D correlation of cytogenetics with cell status abundance (separate by diagnosis) ######################
df <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  filter(
    Diagnosis != "nBM",
    !cell.type.level2 %in% c("gdT", "MAIT", "Platelet", "Stroma")
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM")),
    cell.type.level2 = case_when(
      cell.type.level2 %in% c(
        "CD14_Mono", "CD14CD16_Mono", "CD16_Mono",
        "Macrophage", "cDC", "Neutrophil", "pDC"
      ) ~ "Myeloid",
      TRUE ~ cell.type.level2
    ),
    cell.type.level2 = factor(
      cell.type.level2,
      levels = c("CD4T", "CD8T", "NK", "Mature_B", "Myeloid", "Progenitors")
    )
  ) %>%
  select(NewID, cell.type.level2, cell.status, Diagnosis) %>%
  count(NewID, cell.type.level2, cell.status, Diagnosis, name = "count") %>%
  group_by(NewID, cell.type.level2, Diagnosis) %>%
  complete(cell.status, fill = list(count = 0)) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup()


clinical.info <- readRDS("sample.level.info.rds") %>%
  as.data.frame() %>%
  dplyr::select(
    NewID,
    Hyperdiploidy,
    `t(4;14)`,
    `t(6;14)`,
    `t(11;14)`,
    `t(14;16)`,
    `t(14;20)`,
    `1q gain/CKS1B Pos`,
    `17p loss/p53 Del`,
    `monosomy 13`
  ) %>%
  dplyr::rename(
    t_4_14    = `t(4;14)`,
    t_6_14    = `t(6;14)`,
    t_11_14   = `t(11;14)`,
    t_14_16   = `t(14;16)`,
    t_14_20   = `t(14;20)`,
    x1q_gain  = `1q gain/CKS1B Pos`,
    x17p_loss = `17p loss/p53 Del`,
    mono_13   = `monosomy 13`
  ) %>%
  mutate(
    t_MAF  = if_else(t_14_20 == "Pos", "Pos", t_14_16),
    t_CCND = if_else(t_6_14 == "Pos", "Pos", t_11_14)
  ) %>%
  dplyr::select(
    NewID,
    Hyperdiploidy,
    t_CCND,
    t_4_14,
    t_MAF,
    x1q_gain,
    x17p_loss,
    mono_13
  )


df<-left_join(df,clinical.info,by=c('NewID'='NewID'))


features <- c(
  "Hyperdiploidy", "t_CCND", "t_4_14",
  "t_MAF", "x1q_gain", "x17p_loss", "mono_13"
)

feature_order <- c(
  "Hyperdiploidy", "t_CCND", "t_4_14",
  "t_MAF", "x1q_gain", "x17p_loss", "mono_13"
)

diag_order <- c("MGUS", "SMM", "NDMM", "RRMM")

df2 <- df %>%
  filter(
    total >= 50,
    !is.na(proportion),
    !is.na(Diagnosis)
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = diag_order),
    proportion2 = pmin(pmax(proportion, 0.025), 0.975),
    logit_prop = car::logit(proportion2)
  )

run_model_byDx <- function(feature) {
  df2 %>%
    filter(.data[[feature]] %in% c("Pos", "Neg")) %>%
    mutate(group = factor(.data[[feature]], levels = c("Neg", "Pos"))) %>%
    group_by(Diagnosis, cell.type.level2, cell.status) %>%
    group_modify(~{
      dd <- .x
      
      n_sample <- n_distinct(dd$NewID)
      n_pos <- n_distinct(dd$NewID[dd$group == "Pos"])
      n_neg <- n_distinct(dd$NewID[dd$group == "Neg"])
      
      if (n_sample < 10) return(data.frame())
      if (n_pos < 3 || n_neg < 3) return(data.frame())
      
      fit <- tryCatch(
        lm(logit_prop ~ group, data = dd),
        error = function(e) NULL
      )
      
      if (is.null(fit)) return(data.frame())
      
      out <- broom::tidy(fit) %>%
        filter(term == "groupPos")
      
      if (nrow(out) == 0) return(data.frame())
      
      out %>%
        mutate(
          feature = feature,
          n_sample = n_sample,
          n_pos = n_pos,
          n_neg = n_neg
        )
    }) %>%
    ungroup()
}

results_dx <- map_dfr(features, run_model_byDx) %>%
  group_by(Diagnosis, feature) %>%
  mutate(
    FDR = p.adjust(p.value, method = "BH"),
    odds_ratio = exp(estimate)
  ) %>%
  ungroup()

plot_df <- results_dx %>%
  mutate(
    feature = factor(feature, levels = feature_order),
    Diagnosis = factor(Diagnosis, levels = diag_order),
    sig = case_when(
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    estimate_plot = pmin(pmax(estimate, -0.4), 0.4),
    feature_diag = paste(feature, Diagnosis, sep = " | ")
  )

feature_diag_levels <- unlist(lapply(feature_order, function(f) {
  paste(f, diag_order, sep = " | ")
}))

plot_df$feature_diag <- factor(plot_df$feature_diag, levels = rev(feature_diag_levels))

status_order <- plot_df %>%
  distinct(cell.type.level2, cell.status) %>%
  arrange(cell.type.level2, cell.status) %>%
  pull(cell.status)

plot_df <- plot_df %>%
  mutate(cell.status = factor(cell.status, levels = status_order))

gap_rows <- data.frame(
  Diagnosis = NA,
  cell.type.level2 = unique(plot_df$cell.type.level2)[1],
  cell.status = levels(plot_df$cell.status)[1],
  estimate = NA,
  p.value = NA,
  FDR = NA,
  odds_ratio = NA,
  feature = NA,
  n_sample = NA,
  n_pos = NA,
  n_neg = NA,
  sig = "",
  estimate_plot = NA,
  feature_diag = paste0(feature_order, " | gap")
)

gap_levels <- c()
for (f in feature_order) {
  gap_levels <- c(gap_levels, paste(f, diag_order, sep = " | "), paste0(f, " | gap"))
}
gap_levels <- rev(gap_levels)

plot_df2 <- bind_rows(plot_df, gap_rows) %>%
  mutate(feature_diag = factor(feature_diag, levels = gap_levels))

p <- ggplot(plot_df2, aes(x = cell.status, y = feature_diag, fill = estimate_plot)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(
    data = plot_df2 %>% filter(sig != ""),
    aes(label = sig),
    size = 2.8,
    fontface = "bold"
  ) +
  facet_grid(. ~ cell.type.level2, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B73D9",
    mid = "#F2F2F2",
    high = "#F4A090",
    midpoint = 0,
    limits = c(-0.4, 0.4),
    na.value = "white",
    name = "Log-odds\nratio"
  ) +
  scale_y_discrete(
    labels = function(x) gsub(" \\| gap", "", x)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing.x = unit(0.25, "lines"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

ggsave("heatmap_stratified_grouped_by_feature_optimized.pdf", p, width = 20, height = 10)
