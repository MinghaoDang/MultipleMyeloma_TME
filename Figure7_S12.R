source('global_config.R')

library('broom')
library('emmeans')
library('logistf')
library('patchwork')
library('car')

###################### Figure 7A correlation of Diagnosis with Ecotype ######################
df <- readRDS("sample.level.info.rds") %>%
  filter(
    Diagnosis != "nBM",
    !is.na(Ecotype)) %>%
  transform(
    Diagnosis = factor(Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM"))
  ) %>%
  dplyr::select(Diagnosis, Ecotype)

res<-chisq.test(table(df$Diagnosis, df$Ecotype))
ratio<-res$observed/res$expected;
ratio[ratio>2]<-2
pheatmap::pheatmap(t(ratio),fontsize_col = 10,
                   color = c(colorRampPalette(c('#313695','white'))(40),
                             #colorRampPalette(c('#e0f3f8','#ffffbf'))(2),
                             colorRampPalette(c('white','#ff7f00'))(50)),
                   #scale ='row',
                   border_color = NA,
                   treeheight_row =15,
                   treeheight_col =30,
                   cluster_rows = F,
                   show_rownames = T,
                   show_colnames = T,
                   cluster_cols = F
);
res$p.value

data_plot <- df %>%
  dplyr::count(Diagnosis, Ecotype, name = "Freq") %>%
  dplyr::group_by(Diagnosis) %>%
  dplyr::mutate(percentage = 100 * Freq / sum(Freq)) %>%
  dplyr::ungroup()

plotx<-ggplot(data_plot, aes(x=Diagnosis,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  #facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

###################### Figure 7B correlation of prior lines of treatment with Ecotype ######################

df <- readRDS("sample.level.info.rds") %>%
  transmute(
    Prior_Lines_of_Therapies = if_else(`Prior Lines of Therapies` <= 3, "Less", "More"),
    Ecotype = if_else(Ecotype %in% c("E1", "E2"), "E1+E2", Ecotype),
    Diagnosis = Diagnosis
  ) %>%
  filter(Diagnosis == "RRMM", !is.na(Ecotype), !is.na(Prior_Lines_of_Therapies))

res<-chisq.test(table(df$Prior_Lines_of_Therapies, df$Ecotype))
ratio<-res$observed/res$expected;
ratio[ratio>2]<-2
pheatmap::pheatmap(t(ratio),fontsize_col = 10,
                   color = c(colorRampPalette(c('#313695','white'))(50),
                             #colorRampPalette(c('#e0f3f8','#ffffbf'))(2),
                             colorRampPalette(c('white','#ff7f00'))(50)),
                   #scale ='row',
                   border_color = NA,
                   treeheight_row =15,
                   treeheight_col =30,
                   cluster_rows = F,
                   show_rownames = T,
                   show_colnames = T,
                   cluster_cols = F
);
res$p.value

data_plot <- df %>%
  dplyr::count(Prior_Lines_of_Therapies, Ecotype, name = "Freq") %>%
  dplyr::group_by(Prior_Lines_of_Therapies) %>%
  dplyr::mutate(percentage = 100 * Freq / sum(Freq)) %>%
  dplyr::ungroup()

plotx<-ggplot(data_plot, aes(x=Prior_Lines_of_Therapies,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1+E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  #facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx


###################### Figure 7C correlation of prior anti-CD38 with Ecotype ######################

df <- readRDS("RRMM.treatment_history.rds") %>%
  transmute(
    NewID,
    group = case_when(
      `Recent exposure group` == "anti-CD38" ~ "Recent",
      is.na(`time interval between sample collection and most recent anti-CD38`) ~ "No",
      TRUE ~ `time interval between sample collection and most recent anti-CD38`
    )
  ) %>%
  left_join(
    readRDS("sample.level.info.rds") %>%
      transmute(
        NewID,
        Ecotype = if_else(Ecotype %in% c("E1", "E2"), "E1+E2", Ecotype)
      ),
    by = "NewID"
  )

df$group<-factor(df$group,levels=c('No','Recent','<6mo','>6mo'))

res<-chisq.test(table(df$group, df$Ecotype))
ratio<-res$observed/res$expected;
ratio[ratio>2]<-2
pheatmap::pheatmap(t(ratio),fontsize_col = 10,
                   color = c(colorRampPalette(c('#313695','white'))(50),
                             #colorRampPalette(c('#e0f3f8','#ffffbf'))(2),
                             colorRampPalette(c('white','#ff7f00'))(50)),
                   #scale ='row',
                   border_color = NA,
                   treeheight_row =15,
                   treeheight_col =30,
                   cluster_rows = F,
                   show_rownames = T,
                   show_colnames = T,
                   cluster_cols = F
);
res$p.value

data_plot <- df %>%
  dplyr::count(group, Ecotype, name = "Freq") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(percentage = 100 * Freq / sum(Freq)) %>%
  dplyr::ungroup()

plotx<-ggplot(data_plot, aes(x=group,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1+E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  #facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx


###################### Figure S12B ecotype diagnosis piechart ######################

data_plot <- readRDS("sample.level.info.rds") %>%
  filter(!is.na(Ecotype)) %>%
  count(Diagnosis, Ecotype, name = "Freq") %>%
  group_by(Ecotype) %>%
  mutate(
    percentage = round(100 * Freq / sum(Freq), 2),
    label_text = paste0(round(percentage, 1), "%")
  ) %>%
  ungroup()


plotx<-ggplot(data_plot, aes(x = 2, y = percentage, fill = Diagnosis)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # makes the hole
  facet_grid(. ~ Ecotype) +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(
    nBM  = "#4dac26",
    MGUS = "#1c75bc",
    SMM  = "#fbb040",
    NDMM = "#d01c8b",
    RRMM = "#8856a7"
  )) +
  theme_void() +
  theme(legend.position = "bottom")

plotx

###################### Figure S12C immunotherapy response ######################

df <- readRDS("sample.level.info.rds")

aa <- df %>%
  filter(Diagnosis == "RRMM", !is.na(`response to immunotherapy`)) %>%
  transmute(
    Ecotype,
    response = if_else(
      `response to immunotherapy` %in% c("Short-Term Responder", "Non-responder"),
      "Non/Short-responder",
      `response to immunotherapy`
    )
  )

chisq.test(table(aa$Ecotype, aa$response))

data1 <- aa %>%
  count(Ecotype, response, name = "Freq")

data2 <- data1 %>%
  group_by(response) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()


plotx<-ggplot(data2, aes(x=response,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  #facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

library('ggalluvial')
plotx<-ggplot(data1,
              aes(axis1 = Ecotype, axis2 = response, y = Freq)) +
  geom_alluvium(aes(fill = response), width = 1/12, alpha = 0.55, color = NA) +
  geom_stratum(width = 1/3, aes(fill = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6','Responder'='#00bfc4','Non/Short-responder'='#f8766d'))+
  scale_x_discrete(limits = c("Group", "Ecotype"), expand = c(0.1, 0.1)) +
  theme_minimal()

plotx

###################### Figure D12D correlation of other clinical info with Ecotype ######################

df <- readRDS("sample.level.info.rds") %>%
  mutate(
    Age_group = case_when(
      Age <= 50 ~ "<=50",
      Age <= 60 ~ "<=60",
      Age <= 70 ~ "<=70",
      Age <= 80 ~ "<=80",
      TRUE ~ ">80"
    )
  ) %>%
  dplyr::rename(
    progression = `Progression from precursor`
  )

# gender

chisq.test(t(table(df[,c('Ecotype','Gender')])))

data_plot <- df %>%
  filter(!is.na(Ecotype)) %>%
  dplyr::select(Ecotype,Gender) %>%
  count(Ecotype, Gender, name = "Freq") %>%
  group_by(Gender) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x=Gender,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))
plotx

# age
chisq.test(t(table(df[,c('Ecotype','Age_group')])))

data_plot <- df %>%
  filter(!is.na(Ecotype)) %>%
  dplyr::select(Ecotype,Age_group) %>%
  count(Ecotype, Age_group, name = "Freq") %>%
  group_by(Age_group) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x=Age_group,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

# SMM progression

df.sub <- df %>%
  filter(Diagnosis == "SMM") %>%
  mutate(progression = if_else(is.na(progression), "No", progression))

chisq.test(t(table(df.sub[,c('Ecotype','progression')])))

data_plot <- df.sub %>%
  filter(!is.na(Ecotype)) %>%
  dplyr::select(Ecotype,progression) %>%
  count(Ecotype, progression, name = "Freq") %>%
  group_by(progression) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x=progression,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

# R-ISS stage

df.sub <- df %>%
  filter(Diagnosis %in% c('NDMM','RRMM'), !is.na(Ecotype), !is.na(`R-ISS`)) %>%
  dplyr::rename(R_ISS=`R-ISS`)

chisq.test(t(table(df.sub[,c('Ecotype','R_ISS')])))

data_plot <- df.sub %>%
  filter(!is.na(Ecotype)) %>%
  dplyr::select(Ecotype,R_ISS) %>%
  count(Ecotype, R_ISS, name = "Freq") %>%
  group_by(R_ISS) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x=R_ISS,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

# RRMM prior ASCT

df.sub <- df %>%
  filter(Diagnosis %in% c('RRMM'), !is.na(`Prior ASCT`)) %>%
  dplyr::rename(Prior_ASCT=`Prior ASCT`)

chisq.test(t(table(df.sub[,c('Ecotype','Prior_ASCT')])))

data_plot <- df.sub %>%
  filter(!is.na(Ecotype)) %>%
  dplyr::select(Ecotype,Prior_ASCT) %>%
  count(Ecotype, Prior_ASCT, name = "Freq") %>%
  group_by(Prior_ASCT) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x=Prior_ASCT,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

# RRMM prior Prior_CAR-T/BiTE

df.sub <- df %>%
  filter(Diagnosis %in% c('RRMM'), !is.na(`Prior CAR-T/TCE`)) %>%
  dplyr::rename(Prior_CART_BiTE=`Prior CAR-T/TCE`)

chisq.test(t(table(df.sub[,c('Ecotype','Prior_CART_BiTE')])))

data_plot <- df.sub %>%
  filter(!is.na(Ecotype)) %>%
  dplyr::select(Ecotype,Prior_CART_BiTE) %>%
  count(Ecotype, Prior_CART_BiTE, name = "Freq") %>%
  group_by(Prior_CART_BiTE) %>%
  mutate(
    count = sum(Freq),
    percentage = 100 * Freq / count
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x=Prior_CART_BiTE,fill=Ecotype,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  labs(x='',y = 'Proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 0,size = 10,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx

###################### Figure 7D, S12E, S12F MMRF NMF heatmap, alluvial ######################

meta <- readRDS("per_aliquot_metadata_MMRF_IA.rds") %>%
  filter(VJ_INTERVAL == "Baseline") %>%
  mutate(
    d_tx_induction_cat = case_when(
      d_tx_induction_cat == "chemo_imid_pi_steroid" ~ "imid_pi_steroid",
      d_tx_induction_cat == "chemo_pi_steroid" ~ "pi_steroid",
      TRUE ~ d_tx_induction_cat
    ),
    d_pt_race_1 = case_when(
      d_pt_race_1 %in% c("asian_nos", "unknown") ~ "Others/unknown",
      d_pt_race_1 == "black_african_american" ~ "black",
      TRUE ~ d_pt_race_1
    ),
    d_pt_race_1 = factor(d_pt_race_1, levels = c("white", "black", "Others/unknown")),
    d_dx_amm_iss_stage = as.character(d_dx_amm_iss_stage),
    d_dx_amm_ecog = case_when(
      as.character(d_dx_amm_ecog) %in% c("3", "4") ~ ">=3",
      TRUE ~ as.character(d_dx_amm_ecog)
    ),
    d_dx_amm_ecog = factor(d_dx_amm_ecog, levels = c("0", "1", "2", ">=3"))
  )


data <- readRDS("MMRF.263baseline.meta.data.rds") %>%
  add_count(public_id, name = "n_sample") %>%
  add_count(predicted.id, name = "n_celltype") %>%
  filter(n_sample >= 200, n_celltype >= 100) %>%
  count(public_id, predicted.id, name = "Freq") %>%
  complete(public_id, predicted.id, fill = list(Freq = 0)) %>%
  group_by(public_id) %>%
  mutate(
    count = sum(Freq),
    percentage = Freq / count
  )

ratio <- data %>%
  dplyr::select(public_id,predicted.id,percentage) %>%
  tidyr::pivot_wider(
    names_from = predicted.id,
    values_from = percentage,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("public_id") %>%
  as.matrix()

scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
scale_ratio <- as.data.frame(scale_ratio) %>% t()

ranks <- 2:10
estim.coad <- nmf(scale_ratio, ranks, nrun=500,method = "lee")
plot(estim.coad)


z_ratio <- scale(ratio) %>% as.data.frame() %>% t()

nmf_list<-readRDS('MMRF.transferred_label.res_k2-10.rds')


nmf.rank<-nmf_list[[3]]
library(dplyr)
library(tibble)

B <- basis(nmf.rank)

group_map <- c(`1` = "5", `2` = "3", `3` = "2", `4` = "1")
group_levels <- c("1", "2", "3", "5")

sample.group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "samples")),
  group = predict(nmf.rank, what = "samples")
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = group_levels)
  ) %>%
  left_join(
    meta %>% dplyr::select(public_id, progression_group, davies_based_risk),
    by = c("SampleID" = "public_id")
  ) %>%
  mutate(
    progression_group = factor(progression_group, levels = c("NP", "P", "RP", "Inc")),
    davies_based_risk = factor(
      davies_based_risk,
      levels = c("standard_risk", "high_risk", "not_calculable", "no_risk_data")
    )
  )

sample.order <- with(sample.group, order(group, progression_group, davies_based_risk))

cell.group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "features")),
  group = predict(nmf.rank, what = "features"),
  value = apply(B, 1, max)
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = group_levels)
  )

cell.order <- with(cell.group, order(group, -value))

z_ratio_2 <- z_ratio[cell.group$SampleID[cell.order], sample.group$SampleID[sample.order]]

cell.group.split <- cell.group$group[cell.order]
sample.group.split <- sample.group.2$group[sample.order]

cell.use <- cell.group$SampleID[cell.group$value > 0.02]
keep <- rownames(z_ratio_2) %in% cell.use

z_ratio_3 <- z_ratio_2[keep, ]
cell.group.split <- cell.group.split[keep]

anno.use <- sample.group %>%
  column_to_rownames("SampleID") %>%
  select(-group)

top_ann <- HeatmapAnnotation(
  df = anno.use[colnames(z_ratio_3), ],
  which = "column",
  annotation_name_side = "right",
  col = list(
    progression_group = c(
      NP = "#ffffcc",
      P = "#feb24c",
      RP = "#bd0026",
      Inc = "#878787"
    ),
    davies_based_risk = c(
      standard_risk = "#9e9ac8",
      high_risk = "#54278f",
      not_calculable = "#737373",
      no_risk_data = "#737373"
    )
  ),
  show_legend = TRUE
)

plot_matrix<-z_ratio_3/4
plot_matrix[plot_matrix>0.5]<- 0.5
plot_matrix[plot_matrix< -0.5]<- -0.5

pdf('MMRF.transferred_label.sample_NMF_k4.heatmap.pdf',height=15,width=25)
ht_list<-Heatmap(plot_matrix, name = "ratio",
                 top_annotation = top_ann,
                 row_split = cell.group.split,
                 column_split = sample.group.split,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = F,
                 cluster_columns = F,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 5))
draw(ht_list, annotation_legend_side = "right")


dev.off()

# progression group alluvial
data_progression <- sample.group %>%
  mutate(Ecotype = paste0("E", group)) %>%
  filter(progression_group %in% c("NP", "RP")) %>%
  transmute(
    Ecotype,
    progression_group = as.character(progression_group)
  )

chisq.test(table(data_progression$Ecotype, data_progression$progression_group))

data_plot <- data_progression %>%
  count(Ecotype, progression_group, name = "Freq")

plotx<-ggplot(data_plot,
              aes(axis1 = Ecotype, axis2 = progression_group, y = Freq)) +
  geom_alluvium(aes(fill = Ecotype), width = 1/12, alpha = 0.55, color = NA) +
  geom_stratum(width = 1/3, aes(fill = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E5'='#cab2d6','NP'='#ffffcc','RP'='#bd0026'))+
  scale_x_discrete(limits = c("Ecotype", "Progression"), expand = c(0.1, 0.1)) +
  theme_minimal()

plotx

# risk group alluvial
data_risk <- sample.group %>%
  mutate(Ecotype = paste0("E", group)) %>%
  filter(davies_based_risk %in% c('high_risk','standard_risk')) %>%
  transmute(
    Ecotype,
    davies_based_risk = as.character(davies_based_risk)
  )

chisq.test(table(data_risk$Ecotype, data_risk$davies_based_risk))

data_plot <- data_risk %>%
  count(Ecotype, davies_based_risk, name = "Freq")

data_plot$davies_based_risk<-factor(data_plot$davies_based_risk,levels=c('standard_risk','high_risk'))

plotx<-ggplot(data_plot,
              aes(axis1 = Ecotype, axis2 = davies_based_risk, y = Freq)) +
  geom_alluvium(aes(fill = Ecotype), width = 1/12, alpha = 0.55, color = NA) +
  geom_stratum(width = 1/3, aes(fill = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E5'='#cab2d6','standard_risk'='#9e9ac8','high_risk'='#54278f'))+
  scale_x_discrete(limits = c("Ecotype", "Risk"), expand = c(0.1, 0.1)) +
  theme_minimal()

plotx

saveRDS(sample.group,'MMRF.transferred_label.k4.NMF.sample_group.rds')

###################### Figure 7E, 7F, S12G survival ######################

library(survminer);
library(survival);

meta <- readRDS("per_aliquot_metadata_MMRF_IA.rds") %>%
  filter(VJ_INTERVAL == "Baseline") %>%
  mutate(
    d_tx_induction_cat = case_when(
      d_tx_induction_cat == "chemo_imid_pi_steroid" ~ "imid_pi_steroid",
      d_tx_induction_cat == "chemo_pi_steroid" ~ "pi_steroid",
      TRUE ~ d_tx_induction_cat
    ),
    d_pt_race_1 = case_when(
      d_pt_race_1 %in% c("asian_nos", "unknown") ~ "Others/unknown",
      d_pt_race_1 == "black_african_american" ~ "black",
      TRUE ~ d_pt_race_1
    ),
    d_pt_race_1 = factor(d_pt_race_1, levels = c("white", "black", "Others/unknown")),
    d_dx_amm_iss_stage = as.character(d_dx_amm_iss_stage),
    d_dx_amm_ecog = case_when(
      as.character(d_dx_amm_ecog) %in% c("3", "4") ~ ">=3",
      TRUE ~ as.character(d_dx_amm_ecog)
    ),
    d_dx_amm_ecog = factor(d_dx_amm_ecog, levels = c("0", "1", "2", ">=3"))
  )

sample.group<-readRDS('MMRF.transferred_label.k4.NMF.sample_group.rds')

aa<-left_join(sample.group[,c('SampleID','Ecotype')],meta,by=c('SampleID'='public_id'))

#### PFS ####
sur=Surv(aa$ttcpfs, aa$censpfs);
fit<-survfit(sur~Ecotype,data = aa);

medi<-surv_median(fit)$median

ggsurvplot(
  fit,
  data = aa,
  size = 1,                 # change line size
  conf.int = F,          # Add confidence interval
  pval = TRUE,              # Add p-value
  palette = c("Ecotype=E1"="#a6cee3","Ecotype=E2"="#b2df8a","Ecotype=E3"="#fb9a99","Ecotype=E5"="#cab2d6"),
  risk.table = T,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  xlab = 'Time in days',
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
) -> g

g$plot = g$plot + #labs(title=k) +
  #annotate('text',x=100*0.9,y=0.95,label=i,size=5) +
  annotate('text',x=100*0.9,y=0.85,label=paste('med_1','=',round(medi[1],1)),size=3) +
  annotate('text',x=100*0.9,y=0.8,label=paste('med_2','=',round(medi[2],1)),size=3) +
  annotate('text',x=100*0.9,y=0.75,label=paste('med_3','=',round(medi[3],1)),size=3) +
  annotate('text',x=100*0.9,y=0.7,label=paste('med_4','=',round(medi[4],1)),size=3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
g$table = g$table + labs(title='Number at risk') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

pdf(paste0('MMRF.discovery.baseline.PFS.',Sys.Date(),'.pdf'),height = 6,width = 4.5)
print(g)
dev.off()

sample.group<-rbind(sample.group,c('MMRF_1720',NA,NA,NA,NA))
sample.group<-rbind(sample.group,c('MMRF_1977',NA,NA,NA,NA))

aa<-left_join(sample.group[,c('SampleID','Ecotype')],meta[meta$cohort=='discovery',],by=c('SampleID'='public_id'))
aa$davies_based_risk[aa$davies_based_risk%notin%c('high_risk','standard_risk')]<-NA
aa$davies_based_risk<-factor(aa$davies_based_risk,levels=c('standard_risk','high_risk'))

bigmodel <-coxph( Surv(ttcpfs,censpfs) ~  d_dx_amm_age + d_pt_sex + d_pt_race_1 + d_dx_amm_bmi + d_dx_amm_iss_stage + d_amm_tx_asct_1st + d_tx_induction_cat + d_dx_amm_ecog + Ecotype, data =  aa);
ggforest(bigmodel, data=aa) -> g
ggsave(paste0('MMRF.discovery.baseline.multi.variable.PFS.',Sys.Date(),'.pdf'),g,height=10,width=8)


#### OS ####
sur=Surv(aa$ttcos, aa$censos);
fit<-survfit(sur~Ecotype,data = aa);

medi<-surv_median(fit)$median

ggsurvplot(
  fit,
  data = aa,
  size = 1,                 # change line size
  conf.int = F,          # Add confidence interval
  pval = TRUE,              # Add p-value
  palette = c("Ecotype=E1"="#a6cee3","Ecotype=E2"="#b2df8a","Ecotype=E3"="#fb9a99","Ecotype=E5"="#cab2d6"),
  risk.table = T,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  xlab = 'Time in days',
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
) -> g

g$plot = g$plot + #labs(title=k) +
  annotate('text',x=100*0.9,y=0.85,label=paste('med_1','=',round(medi[1],1)),size=3) +
  annotate('text',x=100*0.9,y=0.8,label=paste('med_2','=',round(medi[2],1)),size=3) +
  annotate('text',x=100*0.9,y=0.75,label=paste('med_3','=',round(medi[3],1)),size=3) +
  annotate('text',x=100*0.9,y=0.7,label=paste('med_4','=',round(medi[4],1)),size=3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
g$table = g$table + labs(title='Number at risk') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

pdf(paste0('MMRF.discovery.baseline.OS.',Sys.Date(),'.pdf'),height = 6,width = 4.5)
print(g)
dev.off()


aa<-left_join(sample.group[,c('SampleID','Ecotype')],meta[meta$cohort=='discovery',],by=c('SampleID'='public_id'))
aa$davies_based_risk[aa$davies_based_risk%notin%c('high_risk','standard_risk')]<-NA
aa$davies_based_risk<-factor(aa$davies_based_risk,levels=c('standard_risk','high_risk'))

bigmodel <-coxph( Surv(ttcos,censos) ~  d_dx_amm_age + d_pt_sex + d_pt_race_1 + d_dx_amm_bmi + d_dx_amm_iss_stage + d_amm_tx_asct_1st + d_tx_induction_cat + d_dx_amm_ecog + Ecotype, data =  aa);
ggforest(bigmodel, data=aa) -> g
ggsave(paste0('MMRF.discovery.baseline.multi.variable.OS.',Sys.Date(),'.pdf'),g,height=10,width=8)

