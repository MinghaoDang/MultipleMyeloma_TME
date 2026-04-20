library('ggplot2')
library('ggpubr')
library('ggrepel')
library('dplyr')
library('Seurat')
library('RColorBrewer')
library('cowplot')
library('pheatmap')
library('ggplotify')
library('readxl')
library('viridis')
library('future')
library('grid')
library('ComplexHeatmap')
library('data.table')
library('gridExtra')
library('purrr')
library('ggsignif')
library('tidyr')

`%notin%`<-Negate(`%in%`)


cell.type.color<-c('Mature_B'='#1f78b4','Pre/Immature_B'='#9ecae1','Pro_B'='#78c679','Prolif_Pro/Pre_B'='#FFBF00',
                   'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67','Prolif_T_NK'='#8c96c6',
                   'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c',
                   'cDC'='#D947E4','pDC'='#dc9ae9','Pro_cDC2'='#e08214','Prolif_pDC'='#3288bd',
                   'Platelet'='#6a3d9a','Pro_Mk'='#2c7fb8',
                   'BaEoMa'='#C20088',
                   'Ortho_E'='#ce1256','Poly_E'='#e7298a','Baso_E'='#df65b0','Pro_E'='#c994c7',
                   'HSC'='#bd0026','GMP'='#b15928',
                   'Endothelial'='#41ab5d','Fibroblast'='#c51b7d')


cell.type.level2.color<-c('Progenitors'='#f4a900',
                          'CD4T'='#936BD8','CD8T'='#33a02c','gdT'='#FFA8BB','MAIT'='#84DF67','NK'='#a6cee3',
                          'Mature_B'='#1f78b4',
                          'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','cDC'='#D947E4',
                          'Neutrophil'='#e31a1c',
                          'pDC'='#dc9ae9',
                          'Platelet'='#6a3d9a','Erythrocytes'='#fc4e2a',
                          'Stroma'='#b15928')

CD8T_color<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#FF7F00","#E31A1C",
              "#FDBF6F","#CAB2D6","#6A3D9A","#FFFF99")
names(CD8T_color)<-c('CD8T_Tn','CD8T_Tcm','CD8T_Tem','CD8T_Trm','CD8T_Teff_CCL4','CD8T_Teff_MHCII','CD8T_Teff_Fos/Jun','CD8T_Teff_TGFB1','CD8T_Teff_FGFBP2','CD8T_Teff_NKlike','CD8T_Tisg')


cell.status.color<-c('HSC'="#E41A1C", 'GMP'="#377EB8", 'Prolif_pDC'="#4DAF4A", 'Prolif_Pro_B'="#984EA3", 'Pro_B'="#FF7F00", 'Prolif_Pre_B'="#FFFF33", 'Pre/Immature_B'="#A65628", 'Prolif_CD4T'="#F781BF", 'Prolif_CD8T'="#999999", 'Prolif_NK'="#66C2A5", 'Pro_Mk'="#FC8D62", 'Pro_E'="#8DA0CB", 'BaEoMa'="#E78AC3",
                     'CD4T_Tn'="#A6CEE3", 'CD4T_Fos/Jun'="#1F78B4", 'CD4T_LGALS1'="#B2DF8A", 'CD4T_Th17'="#33A02C", 'CD4T_Tfh'="#FB9A99", 'CD4T_Treg'="#E31A1C", 'CD4T_CTL_GZMK'="#FDBF6F", 'CD4T_CTL_NKG7'="#FF7F00", 'CD4T_Tisg'="#CAB2D6",
                     'CD8T_Tn'="#6A3D9A", 'CD8T_Tcm'="#228B22", 'CD8T_Tem'="#B15928", 'CD8T_Trm'="#8DD3C7", 'CD8T_Teff_CCL4'="#FFFFB3", 'CD8T_Teff_MHCII'="#BEBADA", 'CD8T_Teff_Fos/Jun'="#FB8072", 'CD8T_Teff_TGFB1'="#80B1D3", 'CD8T_Teff_FGFBP2'="#FDB462", 'CD8T_Teff_NKlike'="#B3DE69", 'CD8T_Tisg'="#FCCDE5",
                     'NK_CD16high'="#66CDAA", 'NK_CD16high_MHCII'="#BC80BD", 'NK_GZMK'="#CCEBC5", 'NK_GZMK_SELL'="#FFED6F", 'NK_GZMK_TOX2'="#FFB300", 'NK_JUND'="#803E75", 'NK_Tisg'="#FF6800",
                     'Bcell_IgMhi_TransB'="#A6BDD7", 'Bcell_Naive'="#FF7F7F", 'Bcell_Isg'="#99C2A2", 'Bcell_Fos/Jun'="#DAB6C4", 'Bcell_Mem'="#9932CC",
                     'CD14_Mono_CCL3_IL1B'="#1B9E77", 'CD14_Mono_S100A8'="#D95F02", 'CD14_Mono_MAP3K8'="#7570B3", 'CD14_Mono_ISG15'="#E7298A", 'CD14_Mono_MHCII'="#66A61E",
                     'CD16_Mono_MT'="#E6AB02", 'CD16_Mono_ISG15'="#A6761D", 'CD16_Mono_CDKN1C'="#8B0000", 'CD16_Mono_CCL3'="#A1C9F4", 'CD16_Mono_JUND'="#FFB482",
                     'Mac_FTL'="#FF9F9B", 'Mac_APOE'="#8DE5A1", 'Mac_S100A4'="#5F9EA0", 'Mac_FRMD4B'="#D2691E", 'Mac_FOS_CCL3'="#DA70D6",
                     'cDC1'="#6B8E23", 'cDC2_MNDA'="#9370DB", 'cDC2_VEGFA'="#BA55D3", 'LAMP3_DC'="#A52A2A", 'pDC'="#1E90FF",
                     'Neutro_CAMP'="#B0E0E6", 'Neutro_MMP9'="#DDA0DD", 'Neutro_CXCL8'="#FA8072", 'Neutro_NEAT1'="#9ACD32", 'Neutro_Isg'="#20B2AA", 'Neutro_ABCA13'="#87CEFA", 'Neutro_TENM1'="#D2691E")


theme_umap<-function(){
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
}

theme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(), 
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9), 
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9), 
      axis.ticks = element_line(color = "white", size  =  0.2), 
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)), 
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)), 
      axis.ticks.length = unit(0.3, "lines"),  
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"), 
      legend.key = element_rect(color = "white",  fill = "black"), 
      legend.key.size = unit(1.2, "lines"), 
      legend.key.height = NULL, 
      legend.key.width = NULL,     
      legend.text = element_text(size = base_size*0.8, color = "white"), 
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"), 
      legend.position = "none", 
      legend.text.align = NULL, 
      legend.title.align = NULL, 
      legend.direction = "vertical", 
      legend.box = NULL,
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA), 
      panel.border = element_rect(fill = NA, color = "white"), 
      ##panel.grid.major = element_line(color = "grey35"), 
      panel.grid.major = element_blank(), 
      ##panel.grid.minor = element_line(color = "grey20"), 
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(0.5, "lines"),  
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"), 
      strip.text.x = element_text(size = base_size*0.8, color = "white"), 
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90), 
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"), 
      plot.title = element_text(size = base_size*1.2, color = "white"), 
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

theme_vlnplot<-function(){
  theme(strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.text.y = element_text(size=12, color="black",angle = 0),
        axis.text.x = element_text(angle = 90,size = 12,vjust=0.5,hjust=0.95),
        text = element_text(size=20),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA))
}

theme_clean<-function(){
  theme(axis.text.x = element_text(size = 12,color='black'),
        axis.text.y = element_text(size = 12,angle=0,color='black'),
        #axis.line = element_blank(),axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
}

signif_plot = function(data,x,y,split_by=NULL,signif.cutoff=0.05, height.factor=0.23,
                       ncol, fill_by=NULL, fill_color='black',
                       color_by=NULL,color='black', ylab) {
  
  # Get full range of y-values
  y_rng = range(data[,y])
  
  # Generate a list of three plots, one for each Species (these are the facets)
  plot_list = lapply(split(data, data[,split_by]), function(d) {
    
    # Get pairs of x-values for current facet
    pairs = combn(levels(data[,x]),2, simplify = F)
    
    # Run wilcox test on every pair
    w.tst =  pairs %>% 
      map_df(function(lv) { 
        p.value = wilcox.test(d$percentage[d[,x]==lv[1]], d$percentage[d[,x]==lv[2]])$p.value
        data.frame(levs=paste(lv, collapse=" "), p.value)
      })
    
    # Record number of significant p.values. We'll use this later to adjust the top of the
    # y-range of the plots
    num_signif = sum(w.tst$p.value <= signif.cutoff)
    
    # Plot significance levels only for combinations with p <= signif.cutoff
    p = ggplot(d, aes_string(x=x, y=y)) +
      geom_boxplot(alpha = 1, show.legend = T,width = 0.75, aes_string(fill=fill_by), outlier.shape=NA) + 
      geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes_string(color=color_by)) +
      facet_grid(as.formula(paste("~", split_by)), scales="free", space="free_x") +
      scale_y_continuous(expand = expansion(mult = c(0.15)))+
      geom_signif(test="wilcox.test", comparisons = pairs[which(w.tst$p.value <= signif.cutoff)],
                  map_signif_level = F,            
                  vjust=0.1,
                  textsize=5,
                  size=0.75,
                  step_increase = 0.15) +
      scale_fill_manual(values=fill_color)+
      scale_color_manual(values=color)+
      theme(strip.text.x = element_text(size=14, color="black", face="bold"),
            axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_blank(),
            axis.title=element_blank(),legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            text = element_text(size=20))
    
    # Return the plot and the number of significant p-values
    return(list(num_signif, p))
  })
  
  # Get the highest number of significant p-values across all three "facets"
  max_signif = max(sapply(plot_list, function(x) x[[1]]))
  
  # Lay out the three plots as facets (one for each Species), but adjust so that y-range is same
  # for each facet. Top of y-range is adjusted using max_signif.
  grid.arrange(grobs=lapply(plot_list, function(x) x[[2]]), 
               ncol=ncol, left=ylab)
}