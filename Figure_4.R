library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')


#######################################
#############  FIGURE 4A  #############
#######################################

seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_bulk2sn_projection.rds')


seurat_ref <- SetActiveWGCNA(seurat_ref , 'bulk2sn')


mapping <- labels2colors(1:100)
MEs <- GetMEs(seurat_ref, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods_num <- paste0('M',match(mods,mapping))
modules_up <- c('M2','M12','M28')
modules_down <- c('M10','M25','M26')



colnames(MEs)<-paste0('M',match(colnames(MEs),mapping))
seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs)
seurat_ref <- SetIdent(seurat_ref, value = "Names")


consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]


library(dplyr)
score_calc <- consensus_modules %>% group_by(module) %>% group_split()
module_colors <- unique(unlist(lapply(score_calc,'[[','module')))
module_colors <- paste0('M',match(module_colors,mapping))


rm(seurat_ref)
gc()
seurat_ref<-readRDS('/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/cm_subclust.rds')
#seurat_ref<-readRDS('~/Downloads/hdWGCNA_TOM/RV_data.rds')

new.cluster.ids <- c("Cm1","Cm2","Cm3","Cm4","Cm5","Cm6","Cm7","Cm8","Cm9","Cm10")
names(new.cluster.ids) <- levels(seurat_ref)
seurat_ref <- RenameIdents(seurat_ref, new.cluster.ids)


idx <- rownames(MEs) %in% colnames(seurat_ref)
seurat_ref$Subnames <- seurat_ref@active.ident
seurat_ref$SubNames_Groups <- paste(seurat_ref$Subnames,seurat_ref$group,sep='_')

seurat_ref <- SetIdent(seurat_ref, value = "Subnames")



seurat_ref@meta.data <- cbind(seurat_ref@meta.data, MEs[idx,])



seurat_ref <- AddModuleScore(seurat_ref,lapply(score_calc,'[[','gene_name'),name="module_score")


cols_current <- colnames(seurat_ref@meta.data)
cols_current[startsWith(colnames(seurat_ref@meta.data),'module_score')] <- paste0('module_',module_colors)
colnames(seurat_ref@meta.data) <- cols_current

#Dot Plot of enrichment cell type of CM enriched modules

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_dot_subclust_up.pdf'), width=5, height=5)

p <- DotPlot(seurat_ref,paste0('module_',modules_up),group.by='Subnames',dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_dot_subclust_down.pdf'), width=5, height=5)

p <- DotPlot(seurat_ref,paste0('module_',modules_down),group.by='Subnames',dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


seurat_ref <- SetIdent(seurat_ref, value = "group")
my_levels <- c("NF","pRV","RVF")
Idents(seurat_ref) <- factor(Idents(seurat_ref), levels= my_levels)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_RV_trend_CM.pdf'), width=4.5, height=2)

p <- DotPlot(seurat_ref,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()

seurat_ref <- SetIdent(seurat_ref, value = "Subnames")


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'sc_seurat_RV_trend_CM_subnames.pdf'), width=4.5, height=5)

p <- DotPlot(seurat_ref,paste0('module_',
  c('M2','M12','M25','M26','M10','M28')),dot.min=0,col.min=0,col.max=2) +
  RotatedAxis() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
p
dev.off()


#######################################
#############  FIGURE 4B  #############
#######################################
Idents(seurat_ref) <- "group"

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_Consensus_Sigs")


#Run enrichment by cell type
combined_set <- data.frame()
combined_output <- data.frame()
bulk_modules <- consensus_modules
bulk_modules$module <- match(consensus_modules$module,mapping)
mods_idx <- list(2,12,10,25,26,28)
cell_types <- unique(seurat_ref$Names)
comparison <- list(c("RVF","NF"),c("RVF","pRV"))
for (i in mods_idx){
  for (j in cell_types){
    for (k in comparison){
      key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
      key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]

      gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(k[1]), ident.2 = paste0(k[2]),features=key_genes)
      
      gene_set<-subset(gene_set,p_val_adj<0.05)
      if (length(rownames(gene_set))==0){next}
      gene_set$module <- paste0('M',i)
      gene_set$color <- mapping[i]
      gene_set$comparison <- paste0(k[1],'_',k[2])
      gene_set$celltype <- j

      if (length(combined_set) == 0){
        combined_set <- gene_set
      }
      else {
        combined_set <- rbind(combined_set,gene_set)
      }

      gene_enrich <- subset(gene_set,avg_log2FC<0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'down'
            combined_output <- rbind(combined_output, cur_df)
          }
      }

      gene_enrich <- subset(gene_set,avg_log2FC>0)
      enriched <- enrichR::enrichr(rownames(gene_enrich), dbs)
      Sys.sleep(5)
      for(db in names(enriched)){
          cur_df <- enriched[[db]]
          if (nrow(cur_df) > 1){
            cur_df$db <- db
            cur_df$module <- paste0('M',i)
            cur_df$celltype <- j
            cur_df$comparison <- paste0(k[1],'_',k[2])
            cur_df$color <- mapping[i]
            cur_df$direction <- 'up'
            combined_output <- rbind(combined_output, cur_df)
          }
      }
    }
  }
}


###RVF vs NF

#Up
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ (GO).*", "")

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


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_by_cluster_terms_cell_type_up_RVF_vs_NF.pdf'), width=6, height=4)
p / colorbar 
dev.off()



###RVF vs NF

#Down
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ (GO).*", "")

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


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_by_cluster_terms_cell_type_down_RVF_vs_NF.pdf'), width=6, height=4)
p / colorbar 
dev.off()



#Both
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
#selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(2,12,10,25,26,28)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1 <- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
#idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_1,]


selected_terms$group <- factor(
  as.character(selected_terms$module_celltype),
  levels = unique(selected_terms$module_celltype)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ (GO).*", "")

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


color_df <- data.frame(group = selected_terms$group,colour = selected_terms$color)



color_df$group <- factor(
  as.character(color_df$group),
  levels = unique(color_df$group)
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_by_cluster_terms_cell_type_both_RVF_vs_NF.pdf'), width=6, height=4)
p / colorbar 
dev.off()


#######################################
#############  FIGURE 4C  #############
#######################################
#Deep dive M2

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs NF
Idents(seurat_ref) <- "group"

combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
  gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M12_genes_up <- rownames(subset(combined_set,module=="M12" & avg_log2FC>0))
M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M12_genes_down <- rownames(subset(combined_set,module=="M12" & avg_log2FC<0))

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

enriched <- enrichr(M2_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M2_enrichr_up.pdf',width=5,height=2.5)
p1<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p1
dev.off()


enriched <- enrichr(M12_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M12_enrichr_up.pdf',width=5,height=2.5)
p2<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p2
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M12_enrichr_up.pdf',width=5,height=5)
p1/p2
dev.off()

enriched <- enrichr(M2_genes_down, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M2_enrichr_down.pdf',width=5,height=2.5)
p3<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p3
dev.off()


enriched <- enrichr(M12_genes_down, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M12_enrichr_down.pdf',width=5,height=2.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p4
dev.off()





bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs NF
Idents(seurat_ref) <- "group"

combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
  gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "pRV",features=key_genes)
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

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M12_genes_up <- rownames(subset(combined_set,module=="M12" & avg_log2FC>0))
M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M12_genes_down <- rownames(subset(combined_set,module=="M12" & avg_log2FC<0))

dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

enriched <- enrichr(M2_genes_up, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_RVF_vs_pRV_M2_enrichr_up.pdf',width=5,height=2.5)
p5<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p5
dev.off()


enriched <- enrichr(M2_genes_down, dbs)
enriched[[4]] <- subset(enriched[[4]],Adjusted.P.value<0.05)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_RVF_vs_pRV_M2_enrichr_down.pdf',width=5,height=2.5)
p6<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p6
dev.off()

pdf('~/Downloads/hdWGCNA_TOM/CM_RV_RVF_vs_NF_and_RVF_vs_pRV_M2_enrichr_up.pdf',width=5,height=4.5)
p1/p5
dev.off()



#######################################
#############  FIGURE 4D  #############
#######################################


library(EnhancedVolcano)

temp_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
temp_set <- subset(temp_set,module == "M2")
keyvals <- temp_set$color
names(keyvals)  <- temp_set$module
temp_set$p_val_adj[temp_set$p_val_adj < 1e-50] = 1e-50
temp_set$avg_log2FC[temp_set$avg_log2FC < -5] = -4.99
temp_set$avg_log2FC[temp_set$avg_log2FC > 5] = 5



pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_all_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(temp_set,lab=rownames(temp_set),
  x='avg_log2FC',y='p_val_adj',xlim=c(-5,5),
  FCcutoff = 0.1) + coord_flip()
dev.off()

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_ACTA1_vln.pdf'), width=3, height=2.8)
VlnPlot(seurat_ref,'ACTA1',group.by='group')
dev.off()



#######################################
#############  FIGURE 4E  #############
#######################################

#Deep dive M12

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs NF
Idents(seurat_ref) <- "group"

combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
  gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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


temp_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
temp_set <- subset(temp_set,module == "M12")
keyvals <- temp_set$color
names(keyvals)  <- temp_set$module
temp_set$p_val_adj[temp_set$p_val_adj < 1e-50] = 1e-50
temp_set$avg_log2FC[temp_set$avg_log2FC < -5] = -4.99
temp_set$avg_log2FC[temp_set$avg_log2FC > 5] = 5

#LRRC39,TNNI3,BMPR2,MYH7,TNNI3K

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_M12_key_vln.pdf'), width=6, height=3)

VlnPlot(seurat_ref,c('TNNI3','MYH7','BMPR2','TNNI3K'),group.by='group',ncol=4)

dev.off()


#######################################
#############  FIGURE 4F  #############
#######################################

library(GeneOverlap)

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules <- subset(bulk_modules,color %in% mods)

sc_modules <- read.csv("~/Downloads/hdWGCNA_TOM/sc_heart_modules.csv")
sc_modules <- subset(sc_modules,str_replace(module,"\\-.*", "")
 == "CM")

bulk_modules <- subset(bulk_modules,gene_name %in% sc_modules$gene_name)
sc_modules <- subset(sc_modules,gene_name %in% bulk_modules$gene_name)

idx<-match(bulk_modules$gene_name,sc_modules$gene_name)

bulk2sc_mods<-data.frame(gene_name=bulk_modules$gene_name,bulk_module=bulk_modules$module,sc_module=sc_modules[idx,2])
bulk2sc_mods$bulk_module <- factor(bulk2sc_mods$bulk_module)
bulk2sc_mods$sc_module <- factor(bulk2sc_mods$sc_module)


#table(bulk2sc_mods$sc_module,match(bulk2sc_mods$bulk_module,mapping))
cont_tbl <- table(bulk2sc_mods$sc_module,bulk2sc_mods$bulk_module)
cont_tbl <- t(cont_tbl) / colSums(cont_tbl)

#cat(subset(bulk2sc_mods,bulk_module == 'tan' & sc_module == 'CM-M1')$gene_name,sep='\n')
#cat(subset(bulk2sc_mods,bulk_module == 'tan' & sc_module == 'CM-M5')$gene_name,sep='\n')
#cat(subset(bulk2sc_mods,bulk_module == 'tan' & sc_module == 'CM-M6')$gene_name,sep='\n')

cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M10')$gene_name,sep='\n')
cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M4')$gene_name,sep='\n')
cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M5')$gene_name,sep='\n')
cat(subset(bulk2sc_mods,bulk_module == 'blue' & sc_module == 'CM-M6')$gene_name,sep='\n')

cat(subset(bulk2sc_mods,sc_module == 'CM-M6' & (bulk_module == 'blue' | bulk_module == 'tan' | bulk_module == 'skyblue'))$gene_name,sep='\n')
cat(subset(bulk2sc_mods,sc_module == 'CM-M5' & (bulk_module == 'blue' | bulk_module == 'tan' | bulk_module == 'skyblue'))$gene_name,sep='\n')

cat(subset(bulk2sc_mods,bulk_module == 'blue')$gene_name,sep='\n')

VlnPlot(seurat_ref,'module_M2',pt.size=0,group.by='group')

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#key_genes <- subset(bulk_modules,module %in% c(10,25,26))$gene_name
#key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
#Idents(seurat_ref) <- "group"


#gene_set <- FindMarkers(seurat_ref, ident.1 = "NF", ident.2 = "RVF",features=key_genes,only.pos=T)
#gene_set<-subset(gene_set,p_val_adj<0.05)


#key_genes <- subset(bulk_modules,module %in% c(2))$gene_name
#key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
#Idents(seurat_ref) <- "group"


#gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes,only.pos=T)
#ene_set<-subset(gene_set,p_val_adj<0.05)

#key_genes <- subset(bulk_modules,module %in% c(12))$gene_name
#key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
#Idents(seurat_ref) <- "group"


#gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes,only.pos=T)
#gene_set<-subset(gene_set,p_val_adj<0.05)


#RVF vs NF
Idents(seurat_ref) <- "group"

combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
  gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_all.pdf'), width=8, height=6)

EnhancedVolcano(combined_set,lab=rownames(combined_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip()
dev.off()


sub_set <- sub_set[!grepl('MT-',rownames(sub_set)),]
sub_set <- subset(combined_set,combined_set$module != 'M2')
keyvals <- sub_set$color
names(keyvals)  <- sub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_no_M2_all.pdf'), width=8, height=6)

EnhancedVolcano(sub_set,lab=rownames(sub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip()
dev.off()

subsub_set <- subset(sub_set,sub_set$module != 'M12')
keyvals <- subsub_set$color
names(keyvals)  <- subsub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_mito_mods_all.pdf'), width=8, height=6)

EnhancedVolcano(subsub_set,lab=rownames(subsub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip(xlim=c(-3, 1)) 
dev.off()

#RVF vs pRV
combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
  gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "pRV",features=key_genes)
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

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_all_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(combined_set,lab=rownames(combined_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip()
dev.off()


sub_set <- sub_set[!grepl('MT-',rownames(sub_set)),]
sub_set <- subset(combined_set,combined_set$module != 'M2')
keyvals <- sub_set$color
names(keyvals)  <- sub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_no_M2_all_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(sub_set,lab=rownames(sub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip()
dev.off()

subsub_set <- subset(sub_set,sub_set$module != 'M12')
keyvals <- subsub_set$color
names(keyvals)  <- subsub_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_mito_mods_all_RVF_vs_pRV.pdf'), width=8, height=6)

EnhancedVolcano(subsub_set,lab=rownames(subsub_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip(xlim=c(-3, 1)) 
dev.off()



##############################################
##############################################
#### Figure 4G
##############################################
##############################################

pdf('~/Downloads/hdWGCNA_TOM/Xenium/CM_Xenium_Trends.pdf',width=20,height=5)

p1 <- VlnPlot(Xenium.cm,'sct_NPPA',group.by='group', pt.size = 0) + NoLegend()
p2 <- VlnPlot(Xenium.cm,'NPPB',group.by='group', pt.size = 0) + NoLegend()
p3 <- VlnPlot(Xenium.cm,'MYH6',group.by='group', pt.size = 0) + NoLegend()
p4 <- VlnPlot(Xenium.cm,'MYH7',group.by='group', pt.size = 0) + NoLegend()
p5 <- VlnPlot(Xenium.cm,'TNNI3',group.by='group', pt.size = 0) + NoLegend()
p6 <- VlnPlot(Xenium.cm,'ANKRD1',group.by='group', pt.size = 0) + NoLegend()
p7 <- VlnPlot(Xenium.cm,'BMPR2',group.by='group', pt.size = 0) + NoLegend()
p8 <- VlnPlot(Xenium.cm,'PANK1',group.by='group', pt.size = 0) + NoLegend()
p9 <- VlnPlot(Xenium.cm,'TMEM65',group.by='group', pt.size = 0) + NoLegend()


p1 | p2 | p4 | p5 | p6 | p7 | p8 | p9 | p3 
dev.off()

##############################################
##############################################
#### Figure 4H
##############################################
##############################################



#M1 <- readRDS(file = "~/Downloads/hdWGCNA_TOM/cm_new_subclust.rds")
#Xenium.cm <- readRDS('~/Downloads/hdWGCNA_TOM/Xenium/cm_minimalist.rds')

















####### OLD















#######################################
#############  FIGURE 3B  #############
#######################################

seurat_obj <- readRDS('~/Downloads/RV_bulkRNASeq_seurat.rds')

#modules <- GetModules(seurat_obj)
#color_df <- modules %>% subset(module!='grey') %>%
#  select(c(module, color)) %>% distinct %>%
#  rename(c(group=module, colour=color))
#mods <- levels(modules$module)
#mods <- mods[mods!='grey']

mods <- mapping[c(2,12,28,10,25,26)]



#Enrichments
library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022')

# compute GO terms:
enrich_list <- list()
seurat_obj <- RunEnrichr(seurat_obj, dbs=dbs,max_genes = 10000)



# helper function to wrap text
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

combined_output <- GetEnrichrTable(seurat_obj)
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,module %in% mapping[c(2,12,28,10,25,26)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ R-HSA.*", "")

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
  levels = paste0('M',match(mods,mapping))
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_reactome_terms.pdf'), width=8, height=7)
p / colorbar 
dev.off()






combined_output <- GetEnrichrTable(seurat_obj)
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,module %in% mapping[c(2,12,28,10,25,26)])


# subset selected terms
selected_terms <- subset(selected_terms, P.value < 0.05)
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_5 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2,idx_top_1+3,idx_top_1+4))

selected_terms<-selected_terms[idx_top_5,]


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = mods
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove Reactome Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\ .*", "")

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
  levels = paste0('M',match(mods,mapping))
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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_chea_terms.pdf'), width=8, height=7)
p / colorbar 
dev.off()


#######################################
#############  FIGURE 3E  #############
#######################################

bulk_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
bulk_modules$module <- match(bulk_modules$module,mapping)

#RVF vs NF
Idents(seurat_ref) <- "group"

combined_set <- data.frame()
mods_idx <- c(2,12,28,10,25,26)
for (i in mods_idx){
  key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
  key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]
  gene_set <- FindMarkers(seurat_ref, ident.1 = "RVF", ident.2 = "NF",features=key_genes)
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

M2_genes_up <- rownames(subset(combined_set,module=="M2" & avg_log2FC>0))
M12_genes_up <- rownames(subset(combined_set,module=="M12" & avg_log2FC>0))
M2_genes_down <- rownames(subset(combined_set,module=="M2" & avg_log2FC<0))
M12_genes_down <- rownames(subset(combined_set,module=="M12" & avg_log2FC<0))


library(EnhancedVolcano)

combined_set <- combined_set[!grepl('MT-',rownames(combined_set)),]
M2_set <- subset(combined_set,module=="M2")

keyvals <- M2_set$color
names(keyvals)  <- M2_set$module

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'CM_module_volcano_M2.pdf'), width=8, height=6)

EnhancedVolcano(M2_set,lab=rownames(combined_set),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,
  colCustom=keyvals) + coord_flip()
dev.off()









#######################################
#############  FIGURE 3F  #############
#######################################


dbs <- c("ChEA_2022","WikiPathway_2023_Human","Reactome_2016","GO_Biological_Process_2023")

library(enrichR)
library(forcats)

parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

enriched <- enrichr(M2_genes_up, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M2_enrichr_up.pdf',width=5,height=2.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p4
dev.off()


enriched <- enrichr(M12_gene_up, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M12_enrichr_up.pdf',width=5,height=2.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p4
dev.off()

enriched <- enrichr(M2_genes_down, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M2_enrichr_down.pdf',width=5,height=2.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p4
dev.off()


enriched <- enrichr(M12_genes_down, dbs)
pdf('~/Downloads/hdWGCNA_TOM/CM_RV_M12_enrichr_down.pdf',width=5,height=2.5)
p4<- ggplot(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + 
  ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + 
  ggtitle('GO Biological Process Up') + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched[[4]][order(enriched[[4]]$Combined.Score,decreasing=T),][rev(1:5),]$Term," \\(GO"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256)))
p4
dev.off()































#######################################
#############  FIGURE 3E  #############
#######################################


seurat_ref <- readRDS('~/Downloads/hdWGCNA_TOM/scWGCNA_all_celltypes.rds')
seurat_ref<-SetActiveWGCNA(seurat_ref, "CM")

# get MEs from seurat object
MEs <- GetMEs(seurat_ref)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
meta <- seurat_ref@meta.data
seurat_ref@meta.data <- cbind(meta, MEs)


# plot with Seurat's DotPlot function
p <- DotPlot(seurat_ref, features=mods, group.by = 'group')

seurat_ref <- RunModuleUMAP(
  seurat_ref,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.3,
  spread=2,
  target_weight=0.1,
  supervised=TRUE
)



