source('global_config.R')

library('broom')
library('emmeans')
library('logistf')
library('patchwork')
library('car')
library('scales')

###################### Figure 6A correlation of ecotype signatures with tumor burden, ITH, proliferating ######################
data<-readRDS('sample.level.info.rds')

data<-data[,c('NewID','Diagnosis',
              'Tumor.pct','S_G2M.pct','Diversity.Score.PC',
              'Sig_E1','Sig_E2','Sig_E3','Sig_E4','Sig_E5')]

df<-subset(data,Diagnosis!='nBM')

df$Diagnosis <- factor(df$Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM"))

cell_types <- c(
  'Sig_E1','Sig_E2','Sig_E3','Sig_E4','Sig_E5'
)

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
    Stage = factor(Stage, levels = c("MGUS", "SMM", "NDMM", "RRMM")),
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


plot_df$effect_z[plot_df$effect_z>3]<-3

p_stage <- ggplot(plot_df, aes(y = CellType, x = Stage)) +
  
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
  facet_wrap(~ Feature, ncol = 3)+
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  scale_size_continuous(
    name = "-log10(FDR)",
    range = c(2, 12.5),
    breaks = c(2.5, 7.5, 12.5)
  ) +
  #scale_size(range = c(2, 15)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("Ecotype.sig.tumor.corr.dotplot.pdf", p_stage, width = 5, height = 6.5)


###################### Figure 6B correlation of ecotype signatures to cytogenetics tumor burden ######################

data<-readRDS('sample.level.info.rds')

data.use<-data[data$Non_Eryth.cell.count>200 & !is.na(data$Sig_E1) & data$Diagnosis!='nBM',c('Diagnosis','Tumor.pct','Sig_E1','Sig_E2','Sig_E3','Sig_E4','Sig_E5')]

data.use<-melt(data.use,id=c('Diagnosis','Tumor.pct'))

plotx <- ggscatter(data.use, y = "value", x = "Tumor.pct",color='Diagnosis',
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE) + 
  stat_cor(method = "spearman") +
  #facet_grid(cols=vars(variable),rows=vars(Diagnosis),scales = 'free') +
  facet_wrap(vars(Diagnosis,variable),nrow=4,ncol=5,scales = 'free') +
  scale_color_manual(values=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20))

ggsave('correlation.of.NMF.signature.with.tumor_burden.pdf',plotx,width=17,height=15)

###################### Figure 6C tumor burden across ecotypes by diagnosis ######################

data<-readRDS('sample.level.info.rds')

data.use<-data[data$Non_Eryth.cell.count>200 & !is.na(data$Sig_E1) & data$Diagnosis!='nBM',c('Diagnosis','Tumor.pct','Ecotype')]
data.use$Ecotype<-factor(data.use$Ecotype,levels=paste0('E',1:5))

P<-ggplot(data.use, aes(x=Ecotype, y=Tumor.pct) ) + 
  geom_boxplot(show.legend = T,width = 0.75, aes(fill=Ecotype), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=1.5) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = combn(levels(factor(data.use$Ecotype)),2, simplify = F),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)
                        vjust=-0.15,
                        textsize=6,
                        size=0.75,
                        step_increase = 0.15)+
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  facet_wrap(vars(Diagnosis),ncol=4,scales='free')+
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="tumor cell proportion") 

###################### Figure 6D, 6E correlation of ecotype signatures with cytogenetics ######################

data<-readRDS('sample.level.info.rds')
data<-data[,c('NewID','Diagnosis','Hyperdiploidy', 't(4;14)', 't(6;14)', 't(11;14)', 't(14;16)', 't(14;20)',
              '1q gain/CKS1B Pos', '17p loss/p53 Del', 'monosomy 13','Sig_E1','Sig_E2','Sig_E3','Sig_E4','Sig_E5')]

colnames(data)<-c('NewID','Diagnosis','Hyperdiploidy','t_4_14','t_6_14','t_11_14','t_14_16','t_14_20','x1q_gain','x17p_loss','mono_13','Sig_E1','Sig_E2','Sig_E3','Sig_E4','Sig_E5')

data$t_MAF<-data$t_14_16
data$t_MAF[data$t_14_20=='Pos']<-'Pos'

data$t_CCND<-data$t_11_14
data$t_CCND[data$t_6_14=='Pos']<-'Pos'

data<-data[,c('NewID','Diagnosis','Hyperdiploidy','t_CCND','t_4_14','t_MAF','x1q_gain','x17p_loss','mono_13','Sig_E1','Sig_E2','Sig_E3','Sig_E4','Sig_E5')]
data<-subset(data,Diagnosis!='nBM')

data<-melt(data,id.vars = c('NewID', 'Diagnosis', 'Hyperdiploidy', 't_CCND', 't_4_14', 't_MAF', 'x1q_gain', 'x17p_loss', 'mono_13'))

library(limma)
library(car)


features <- c(
  "Hyperdiploidy", "t_CCND", "t_4_14", "t_MAF",
  "x1q_gain", "x17p_loss", "mono_13"
)


run_one_feature <- function(feature) {
  data %>%
    filter(.data[[feature]] %in% c("Pos", "Neg")) %>%
    mutate(group = factor(.data[[feature]], levels = c("Neg", "Pos"))) %>%
    group_by(variable) %>%
    group_modify(~{
      fit <- lm(value ~ group + Diagnosis, data = .x)
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

heatmap_df$feature<-factor(heatmap_df$feature,levels=c("Hyperdiploidy","t_CCND","t_4_14","t_MAF","x1q_gain","x17p_loss","mono_13"))

plotx<-ggplot(heatmap_df, aes(y = feature, x = variable, fill = estimate)) +
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

ggsave('Ecotype.sig.cytogenetics.corr.heatmap.pdf',plotx,width=5,height=4)


P<-ggplot(subset(data,x17p_loss!='NA'), aes(x=x17p_loss, y=value) ) + 
  geom_boxplot(show.legend = T,width = 0.75, aes(fill=x17p_loss), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=1.5, aes(color=Diagnosis)) +
  stat_compare_means(
    comparisons = list(c('Pos','Neg')),
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE
  )+
  scale_color_manual(values=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  scale_fill_manual(values=c('Pos'='#ff7f00','Neg'='#8dd3c7'))+
  facet_grid(.~variable)+
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="Signatures") 
ggsave('correlation.of.NMF.signature.with.17q_loss.pdf',P,width=15,height=5)





###################### Figure S12A adjust distribution of cytogenetics in ecotypes ######################

dat <- readRDS("sample.level.info.rds") %>%
  as.data.frame() %>%
  dplyr::select(
    NewID, Diagnosis, Ecotype, Hyperdiploidy,
    `t(4;14)`, `t(6;14)`, `t(11;14)`, `t(14;16)`, `t(14;20)`,
    `1q gain/CKS1B Pos`, `17p loss/p53 Del`, `monosomy 13`
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
    Diagnosis, Ecotype,
    Hyperdiploidy, t_CCND, t_4_14, t_MAF, x1q_gain, x17p_loss, mono_13
  ) %>%
  filter(Diagnosis != "nBM", !is.na(Ecotype)) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM")),
    Ecotype = factor(Ecotype, levels = c("E1", "E2", "E3", "E4", "E5"))
  )

event_order <- c(
  "Hyperdiploidy", "t_CCND", "t_4_14", "t_MAF",
  "x1q_gain", "x17p_loss", "mono_13"
)

# rare events to use Firth
rare_events <- c("t_4_14", "t_MAF")

fit_event <- function(event, n_boot = 300, seed = 1) {
  d <- dat %>%
    dplyr::select(Diagnosis, Ecotype, all_of(event)) %>%
    filter(!is.na(.data[[event]])) %>%
    mutate(y = if_else(.data[[event]] == "Pos", 1, 0)) %>%
    droplevels()
  
  if (nrow(d) == 0 || length(unique(d$y)) < 2) return(NULL)
  
  # global p-value: full vs reduced model
  full_glm <- glm(y ~ Ecotype + Diagnosis, data = d, family = binomial())
  red_glm  <- glm(y ~ Diagnosis, data = d, family = binomial())
  p <- anova(red_glm, full_glm, test = "LRT")$`Pr(>Chi)`[2]
  
  if (!(event %in% rare_events)) {
    out <- emmeans(full_glm, ~ Ecotype, type = "response") %>%
      as.data.frame() %>%
      transmute(
        Event = event,
        Ecotype,
        Adj_Prob = prob,
        CI_low = asymp.LCL,
        CI_high = asymp.UCL,
        p_value = p
      )
    return(out)
  }
  
  # Firth for rare events
  set.seed(seed)
  fit_firth <- logistf(y ~ Ecotype + Diagnosis, data = d)
  
  newdat <- expand.grid(
    Ecotype = levels(d$Ecotype),
    Diagnosis = levels(d$Diagnosis)
  )
  newdat$pred <- predict(fit_firth, newdata = newdat, type = "response")
  
  point_est <- newdat %>%
    group_by(Ecotype) %>%
    summarise(
      Adj_Prob = mean(pred),
      .groups = "drop"
    )
  
  boot_fun <- function() {
    idx <- sample(seq_len(nrow(d)), replace = TRUE)
    db <- droplevels(d[idx, , drop = FALSE])
    
    if (length(unique(db$y)) < 2) {
      return(rep(NA_real_, length(levels(d$Ecotype))))
    }
    
    fit_b <- tryCatch(
      logistf(y ~ Ecotype + Diagnosis, data = db),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      return(rep(NA_real_, length(levels(d$Ecotype))))
    }
    
    nd <- expand.grid(
      Ecotype = levels(d$Ecotype),
      Diagnosis = levels(d$Diagnosis)
    )
    nd$pred <- predict(fit_b, newdata = nd, type = "response")
    
    nd %>%
      group_by(Ecotype) %>%
      summarise(Adj_Prob = mean(pred), .groups = "drop") %>%
      arrange(match(Ecotype, levels(d$Ecotype))) %>%
      pull(Adj_Prob)
  }
  
  boot_mat <- replicate(n_boot, boot_fun())
  
  point_est %>%
    mutate(
      CI_low = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
      CI_high = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
      Event = event,
      p_value = p
    ) %>%
    dplyr::select(Event, Ecotype, Adj_Prob, CI_low, CI_high, p_value)
}

# run all events
res <- map_dfr(event_order, fit_event, n_boot = 300, seed = 1)

# BH correction at event level
p_tab <- res %>%
  distinct(Event, p_value) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

plot_all <- res %>%
  left_join(p_tab, by = c("Event", "p_value")) %>%
  mutate(
    Event = factor(Event, levels = event_order),
    Ecotype = factor(Ecotype, levels = c("E1", "E2", "E3", "E4", "E5")),
    facet_lab = paste0(
      Event,
      "\nGlobal LRT p=", sprintf("%.3g", p_value),
      ", q=", sprintf("%.3g", p_adj_BH)
    )
  )

p_all <- ggplot(plot_all, aes(x = Ecotype, y = Adj_Prob)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.12) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Ecotype",
    y = "Adjusted probability (95% CI)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold")
  )

p_all
