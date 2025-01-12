library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)


source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')



#######################################
#############  FIGURE 6A  #############
#######################################
M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/ec_subclust.rds")

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Names',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()

#######################################
#############  FIGURE 6B  #############
#######################################
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_dot.pdf'), width=4.5, height=4)
DotPlot(M1, c("VenousScore1","ArterialScore1","LymphScore1",
  "CapillaryScore1","EndoScore1"),col.min=0)+xlab('Marker Score')+
  scale_x_discrete(labels=c('1','2','3','4','5'))
dev.off()


#######################################
#############  FIGURE 6C  #############
#######################################


M1 <- SetIdent(M1, value = "Names")



all_percent_cell <- as.data.frame(table(M1@active.ident)/length(M1@active.ident)*100)

NF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="NF")@active.ident)/length(subset(M1,group=="NF")@active.ident)*100),type = "NF")
NF_percent_cell$sum <- (rev(cumsum(rev(NF_percent_cell$Freq))) - NF_percent_cell$Freq/2)/100
NF_percent_cell$Freq <- NF_percent_cell$Freq/100


pRV_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="pRV")@active.ident)/length(subset(M1,group=="pRV")@active.ident)*100),type = "pRV")
pRV_percent_cell$sum <- (rev(cumsum(rev(pRV_percent_cell$Freq))) - pRV_percent_cell$Freq/2)/100
pRV_percent_cell$Freq <- pRV_percent_cell$Freq/100


RVF_percent_cell <- cbind(as.data.frame(table(subset(M1,group=="RVF")@active.ident)/length(subset(M1,group=="RVF")@active.ident)*100),type = "RVF")
RVF_percent_cell$sum <- (rev(cumsum(rev(RVF_percent_cell$Freq))) - RVF_percent_cell$Freq/2)/100
RVF_percent_cell$Freq <- RVF_percent_cell$Freq/100


percent_cell_df <- rbind(NF_percent_cell,pRV_percent_cell,RVF_percent_cell)

pdf('~/Downloads/hdWGCNA_TOM/EC_subclust_prev_stacked.pdf',width=5,height=5)
ggplot(percent_cell_df, aes(fill=Var1, y=Freq, x=type,label=round(sum,1))) +  geom_bar(position="stack", stat="identity",width=0.6) + theme_classic() + xlab("Disease State") + ylab("Frequency") + labs(fill="Cell type",color='black') + theme(text = element_text(size=20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),legend.text=element_text(color="black")) + scale_y_continuous(expand=c(0,0)) + geom_label_repel(aes(type,sum,label=scales::percent(round(Freq,2))),fill=NA,nudge_x=0.5,direction="y")
dev.off()


#Prevalence comparisons
cells <- table(M1@active.ident,M1@meta.data$patient)
cells <- sweep(cells,2,colSums(cells),'/')
cells <- data.frame(cells)


#ggboxplot(cells[77:1,],x="Var2",y="Freq",fill="Var2",group="Var2")+theme_classic() + theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),legend.title=element_text(size=16),legend.text=element_text(size=16),text=element_text(color='black'),axis.text=element_text(color='black')) + labs(color='Group',x="Disease",y='Frequency') + facet_wrap(~Var1,ncol=7)


cells$group = rep(c("RVF","pRV","RVF","NF","pRV","pRV","RVF","NF","NF","pRV","NF"),each=5)


#my_comparisons <- list( c("NF", "pRV"),c("pRV", "RVF"),c("NF", "RVF"))

library(ggpubr)
pdf('~/Downloads/hdWGCNA_TOM/EC_clust_freq.pdf',width=12.5,height=5)
p <- ggboxplot(cells[55:1,],x="group",y="Freq",fill="group",group="group")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=16),
  axis.text.y=element_text(size=16),
  axis.title.x=element_text(size=16),
  axis.title.y=element_text(size=16),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  text=element_text(color='black'),
  axis.text=element_text(color='black')) + 
  labs(color='Group',x="Disease",y='Frequency') + 
  facet_wrap(~Var1,ncol=7) + 
  #stat_compare_means(aes(group=group),comparisons=my_comparisons,method="t.test",ref.group="NF")+
  stat_compare_means(aes(group=group),method="anova")
p
dev.off()


#######################################
#############  FIGURE 6D  #############
#######################################
DefaultAssay(M1) <- "RNA"


M1 <- SetupForWGCNA(
  M1,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "ec" # the name of the hdWGCNA experiment
)

M1 <- MetacellsByGroups(seurat_obj = M1,
  group.by = c("Names"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Names' # set the Idents of the metacell seurat object
)

M1 <- NormalizeMetacells(M1)


M1 <- SetDatExpr(
  M1,
  group_name = c("Capillary","Arterial","Venous"), # the name of the group of interest in the group.by column
  group.by='Names', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

M1 <- TestSoftPowers(
  M1,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

library(patchwork)
plot_list <- PlotSoftPowers(M1)
wrap_plots(plot_list, ncol=2)

M1 <- ConstructNetwork(
  M1,
  tom_name = 'ec_net' # name of the topoligical overlap matrix written to disk
)


M1 <- ScaleData(M1)
M1 <- ModuleEigengenes(
 M1,
 group.by.vars="patient"
)


MEs <- GetMEs(M1,harmonized=TRUE)

M1 <- ModuleConnectivity(
  M1,
  group.by = 'Names', group_name =  c("Capillary","Arterial","Venous")
)

modules <- GetModules(M1)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

#Blue is CM - M2
#Yellow module EC specific - M4
#Green is endocardial - M5
#Brown is SM - M3
#Turqouise is EC specific - M1
#Red is mural cells and fibroblast - M6
#Black is ECs - M7

# add hMEs to Seurat meta-data:
M1@meta.data <- cbind(M1@meta.data, MEs)

mapping <- labels2colors(1:100)
mods_num <- paste0('M',match(mods,mapping))

colNames <- colnames(M1@meta.data)
colNames[match(mods,colNames)] <- mods_num
colnames(M1@meta.data) <- colNames

# compute the module UMAPs
M1 <- RunModuleUMAP(
  M1,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.3,
  spread=2,
  target_weight=0.1,
  supervised=TRUE
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(M1)

# plot with ggplot
plot_df <- umap_df

# compute coordinates for cluster labels
centroid_df <- data.frame()
for(cur_cluster in unique(plot_df[['module']])){
    cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
    df <- data.frame(
      cluster = cur_cluster,
      UMAP1 = mean(cur_meta$UMAP1),
      UMAP2 = mean(cur_meta$UMAP2)
    )
  centroid_df <- rbind(centroid_df, df)
}

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    ggrastr::rasterise(geom_point(
     color=umap_df$color,
     size=umap_df$kME*2
   ), dpi=500, scale=0.5) +
    umap_theme() +
    theme(
      plot.margin = margin(0,0,0,0),
      plot.title = element_text(hjust=0.5)
    ) + ggtitle('Bulk RNASeq') +
    ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=3)


hub_genes <- GetHubGenes(M1, 3)

  # add annotation
anno_genes <- hub_genes$gene_name
plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')

plot_df_anno <- subset(plot_df, anno != '')
p <-  plot_df %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
    ggrastr::rasterise(
      geom_point(
        inherit.aes=FALSE,
        data=plot_df,
        aes(x=UMAP1, y=UMAP2, color=module),
        color=plot_df$color,
        size=plot_df$kME*2,
      ), dpi=500, dpi_scale=0.5) +
    geom_point(
      inherit.aes = FALSE,
      data = plot_df_anno,
      shape=21, color='black',
      fill=plot_df_anno$color,
      size=plot_df_anno$kME*2,
      aes(x=UMAP1, y=UMAP2, fill=module)
    ) +
    # add labels
    ggrepel::geom_text_repel(data = centroid_df, label=paste0('M',match(centroid_df$cluster,mapping)), color='black', max.overlaps=Inf, size=3, fontface='bold') +
    geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=3) +
    umap_theme() + NoLegend() +
    coord_equal() +
    theme(
      plot.margin = margin(0,0,0,0)
    )

pdf('~/Downloads/hdWGCNA_TOM/EC_hdWGCNA.pdf',width=6,height=6)
print(p)
dev.off()


library(enrichR)
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
M1 <- RunEnrichr(
  M1,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 10000 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(M1)

EnrichrBarPlot(
  M1,
  outdir = "~/Downloads/hdWGCNA_TOM/scEC_subclust_enrichr_plot", 
  n_terms = 5,
  plot_size = c(5,4), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

#######################################
#############  FIGURE 6E  #############
#######################################


modules <- GetModules(M1)
color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  rename(c(group=module, colour=color))
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(M1)
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")

# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]

order_mod <- c(1,4,5,7,2,3,6)

selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mapping[order_mod]
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, " \\s*\\([^\\)]+\\)", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 45)

selected_terms$Term <- factor(
  as.character(selected_terms$Term),
  levels = rev(unique(as.character(selected_terms$Term)))
)

selected_terms$wrap <- factor(
  as.character(selected_terms$wrap),
  levels = rev(unique(as.character(selected_terms$wrap)))
)

library(viridis)

# Reactome Term dot plot
p <- selected_terms %>%
  ggplot(aes(x = group, y = wrap, color =logp, size=log(Combined.Score))) +
  geom_point() +
  scale_color_stepsn(colors=rev(magma(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )


mapping <- labels2colors(1:100)


color_df <- modules %>% subset(module!='grey') %>%
  select(c(module, color)) %>% distinct %>%
  rename(c(group=module, colour=color))
color_df <- subset(color_df,colour %in% mods)


color_df$group<-paste0('M',match(color_df$group,mapping))

color_df$group <- factor(
  as.character(color_df$group),
  levels = paste0('M',order_mod)
)


c_vect <- color_df$colour
names(c_vect) <- color_df$group



# make the colorbar as its own heatmap
color_df$var <- 1
colorbar <- color_df %>%
  ggplot(aes(x=group, y=var, fill=group)) +
  geom_tile() +
  scale_fill_manual(values=c_vect) +
  coord_equal() +
  NoLegend() + RotatedAxis() +
  theme(
    plot.title=element_blank(),
    axis.line=element_blank(),
    axis.ticks.y =element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.margin=margin(0,0,0,0),
  )


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_GO_terms.pdf'), width=8, height=7)
p / colorbar 
dev.off()




#######################################
#############  FIGURE 6F  #############
#######################################


pdf('~/Downloads/hdWGCNA_TOM/EC_celltype_mods.pdf',width=6,height=3)
p <- DotPlot(M1, features=sort(mods_num), group.by = 'Names',col.min=0)
p
dev.off()


pdf('~/Downloads/hdWGCNA_TOM/EC_group_mods.pdf',width=6,height=2)
p <- DotPlot(M1, features=sort(mods_num), group.by = 'group',col.min=0)
p
dev.off()

up_RVF = c('M1','M4','M5','M7')
down_RVF = c('M2','M3','M6')


#######################################
#############  FIGURE 6G  #############
#######################################

pdf('~/Downloads/hdWGCNA_TOM/EC_feature.pdf',width=4,height=6)
FeaturePlot(M1,c('M1','M4','M7'),label=T,min.cutoff=0,ncol=1)
dev.off()

#######################################
#############  FIGURE 6H  #############
#######################################
M1 <- SetIdent(M1, value = "group")

modules <- GetModules(M1)
modules$module <- match(modules$module,mapping)
combined_set <- data.frame()
mods_idx <- c(1,4,5,7,2,3,6)
for (i in mods_idx){
  key_genes <- subset(modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(M1)]
  gene_set <- FindMarkers(M1, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
  gene_set<-subset(gene_set,p_val_adj<0.05)
  gene_set$module <- paste0('M',i)
  gene_set$color <- mapping[i]
  if (length(combined_set) == 0){
    combined_set <- gene_set
  }
  else {
    combined_set <- rbind(combined_set,gene_set)
  }
}


library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
keyvals <- combined_set$color
names(keyvals)  <- combined_set$module



M1_subset <- subset(combined_set,color=='turquoise')
M1_subset_up <- subset(M1_subset,avg_log2FC>0)
M1_subset_down <- subset(M1_subset,avg_log2FC<0)
keyvals <- M1_subset$color
names(keyvals)  <- M1_subset$module


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'RV_EC_module_volcano_aM1.pdf'), width=8, height=6)

EnhancedVolcano(M1_subset,lab=rownames(M1_subset),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip()
dev.off()




cat(rownames(M1_subset_up),sep='\n')
cat(rownames(M1_subset_down),sep='\n')


M4_subset <- subset(combined_set,module=='M4')
M4_subset_up <- subset(M4_subset,avg_log2FC>0)
M4_subset_down <- subset(M4_subset,avg_log2FC<0)


cat(rownames(M4_subset_up),sep='\n')
cat(rownames(M4_subset_down),sep='\n')

#saveRDS(M1,'~/Downloads/hdWGCNA_TOM/EC_hdWGCNA_by_celltype.rds')


DefaultAssay(M1) <- "SCT"

#VlnPlot(subset(M1,Names=='Arterial'),'MECOM',group.by='group')
#VlnPlot(M1,'VEGFC',group.by='group')
#VlnPlot(M1,'VASH1',group.by='group')
#VlnPlot(M1,'FLT1',group.by='group')
#VlnPlot(M1,'SP100',group.by='group')
#VlnPlot(M1,'SEMA6A',group.by='group')
#VlnPlot(M1,'PPP1R16B',group.by='group')
#VlnPlot(M1,'PLXND1',group.by='group')


#Get all terms for + regulation of angiogenesis from GO
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0045765', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="AngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0045766', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="PositiveAngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0016525', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="NegativeAngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0060978', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="CoronaryAngioGenes")

gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0001525', mart = ensembl)
angio_genes <- unique(gene.data$hgnc_symbol)
angio_genes <- intersect(angio_genes,rownames(M1))
M1 <- AddModuleScore(M1,list(angio_genes),name="AllAngioGenes")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_all_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'AllAngioGenes1')
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_coronary_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'CoronaryAngioGenes1')
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_negative_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'NegativeAngioGenes1')
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_positive_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'PositiveAngioGenes1')
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_regulatory_angio_genes.pdf'), width=8, height=6)
VlnPlot(M1,'AngioGenes1')
dev.off()


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_combined_angio_genes.pdf'), width=10, height=4)
VlnPlot(M1,c('AllAngioGenes1','CoronaryAngioGenes1','PositiveAngioGenes1','NegativeAngioGenes1'),group.by='group',ncol=4,pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_MECOM.pdf'), width=3, height=3)
VlnPlot(subset(M1,Names=='Arterial'),'SMAD1',group.by='group',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_MECOM_feat.pdf'), width=4, height=3)
DotPlot(M1,'MECOM',group.by='Names',col.min=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_SMAD1.pdf'), width=3, height=3)
VlnPlot(subset(M1,Names=='Capillary'),'SMAD1',group.by='group',pt.size=0)
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'EC_SMAD1_feat.pdf'), width=4, height=3)
DotPlot(M1,'SMAD1',group.by='Names',col.min=0)
dev.off()



#######################################
#############  FIGURE 6I  #############
#######################################






















DefaultAssay(M1) <- "RNA"


M1 <- SetupForWGCNA(
  M1,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "ec_group" # the name of the hdWGCNA experiment
)

M1 <- MetacellsByGroups(seurat_obj = M1,
  group.by = c("group"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'group' # set the Idents of the metacell seurat object
)

M1 <- NormalizeMetacells(M1)


M1 <- SetDatExpr(
  M1,
  group_name = c("NF","pRV","RVF"), # the name of the group of interest in the group.by column
  group.by='group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

M1 <- TestSoftPowers(
  M1,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

library(patchwork)
plot_list <- PlotSoftPowers(M1)
wrap_plots(plot_list, ncol=2)

M1 <- ConstructNetwork(
  M1,
  tom_name = 'ec_net_group' # name of the topoligical overlap matrix written to disk
)


M1 <- ScaleData(M1)
M1 <- ModuleEigengenes(
 M1,
 group.by.vars="patient"
)


MEs <- GetMEs(M1,harmonized=TRUE)

M1 <- ModuleConnectivity(
  M1,
  group.by = 'group', group_name =  c("NF","pRV","RVF")
)

modules <- GetModules(M1)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
saveRDS(M1,'~/Downloads/hdWGCNA_TOM/EC_hdWGCNA_by_group.rds')



