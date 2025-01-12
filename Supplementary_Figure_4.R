library(Seurat)
library(hdWGCNA)
library(ggeasy)
library(harmony)



source('~/Downloads/hdWGCNA_TOM/spatial_functions.R')

#######################################
#############  FIGURE S5A  ############
#######################################

M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/fb_subclust.rds")

new.cluster.ids <- c("Fb1","Fb2","Fb3","Fb4","Fb5","Fb6","Fb7")
names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Subnames',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()



#######################################
#############  FIGURE S5B  ############
#######################################

M1 <- AddModuleScore(M1, features=list(c("SCN7A","ACSM3")),assay="SCT",name="Fb1Score")
M1 <- AddModuleScore(M1, features=list(c("KAZN","CNTNAP2","C7")),assay="SCT",name="Fb2Score")
M1 <- AddModuleScore(M1, features=list(c("LTBP2","PLXDC2","ELN","NOX4","FGF14")),assay="SCT",name="Fb3Score")
M1 <- AddModuleScore(M1, features=list(c("PCOLCE2","FBN1","MFAP5","CREB5")),assay="SCT",name="Fb4Score")
M1 <- AddModuleScore(M1, features=list(c("GRID2","NAMPT","NR4A3","NR4A1")),assay="SCT",name="Fb5Score")
M1 <- AddModuleScore(M1, features=list(c("SERPINE1","DEC1","TNC","FN1")),assay="SCT",name="Fb6Score")
M1 <- AddModuleScore(M1, features=list(c("SYN3","TIMP3")),assay="SCT",name="Fb7Score")

#Fb1 (SCN7A, ACSM3), Fb2 (KAZN, CNTNAP2, C7), Fb3 (LTBP2, PLXDC2, ELN, NOX4, FGF14)
#Fb4 (PCOLCE2, FBN1, MFAP5, CREB5), Fb5 (GRID2, NAMPT, NR4A3, NR4A1)
#Fb6 (SERPINE1, DEC1 TNC, FN1), Fb7 (SYN3, TMIP3)

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_sn_Dot.pdf'), width=6, height=5)

DotPlot(M1,features=c("Fb1Score1",
	"Fb2Score1","Fb3Score1","Fb4Score1",
	"Fb5Score1","Fb6Score1","Fb7Score1"),
col.min=0,col.max=2)+ylab('Condition')+xlab('Marker Score')+
scale_x_discrete(labels=c('Clust 1','Clust 2','Clust 3',
	'Clust 4','Clust 5','Clust 6','Clust 7'))
dev.off()


a <- FindAllMarkers(M1)
b <- split( a , f = a$cluster )
c <-  lapply(b,subset,avg_log2FC>0 & p_val_adj < 0.05)
d <- lapply(c,'[[','gene')


M1 <- AddModuleScore(M1,d,name='ClustMarkers')

DotPlot(M1,features=c("ClustMarkers1",
	"ClustMarkers2","ClustMarkers3","ClustMarkers4",
	"ClustMarkers5","ClustMarkers6","ClustMarkers7"),
col.min=0,col.max=2)+ylab('Condition')+xlab('Marker Score')+
scale_x_discrete(labels=c('Clust 1','Clust 2','Clust 3',
	'Clust 4','Clust 5','Clust 6','Clust 7'))

#######################################
#############  FIGURE S5C  ############
#######################################

library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2021','GO_Molecular_Function_2021','LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down', 'WikiPathway_2021_Human', 'KEGG_2021_Human')

# compute GO terms:
enrich_list <- list()

combined_output <- data.frame()
for (i in 1:length(unique(M1$Subnames))) {
    cur_mod <- unique(M1$Subnames)[i]
    cur_info <- subset(a, cluster == cur_mod & avg_log2FC>0 & p_val_adj < 0.05)
    cur_info <- cur_info[, c("gene")]
    
    enriched <- enrichR::enrichr(cur_info, dbs)
    Sys.sleep(5)
    for (db in names(enriched)) {
        cur_df <- enriched[[db]]
        if (nrow(cur_df) > 1) {
            cur_df$db <- db
            cur_df$module <- cur_mod
            combined_output <- rbind(combined_output, cur_df)
        }
    }
}

M1 <- SetupForWGCNA(M1,wgcna_name='temp')
M1 <- SetEnrichrTable(M1, combined_output)

outdir = '~/Downloads/hdWGCNA_TOM/scFB_subclust_enrichr_plot'


wrapText <- function(x, len) {
        sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), 
            USE.NAMES = FALSE)
}

 enrichr_df <- GetEnrichrTable(M1)

for (i in 1:length(unique(M1$Subnames))) {
    cur_mod <- unique(M1$Subnames)[i]
    cur_terms <- subset(enrichr_df, module == cur_mod)
    print(cur_mod)
    #cur_color <- modules %>% subset(module == cur_mod) %>% 
    #    .$color %>% unique %>% as.character
    #if (!is.null(plot_bar_color)) {
    #    cur_color <- plot_bar_color
    #}
    if (nrow(cur_terms) == 0) {
        next
    }
    cur_terms$wrap <- wrapText(cur_terms$Term, 45)
    plot_list <- list()
    for (cur_db in dbs) {
        plot_df <- subset(cur_terms, db == cur_db) %>% top_n(5, 
            wt = Combined.Score)
        text_color = "black"
        cur_color = 'green'

        plot_df$Combined.Score <- log(plot_df$Combined.Score)
        lab <- "Enrichment log(combined score)"
        x <- 0.2

        plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Combined.Score, 
            y = reorder(wrap, Combined.Score))) + geom_bar(stat = "identity", 
            position = "identity", color = "white", fill = cur_color) + 
            geom_text(aes(label = wrap), x = x, color = text_color, 
              size = 3.5, hjust = "left") + ylab("Term") + 
            xlab(lab) + ggtitle(cur_db) + theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.title = element_blank(), 
            axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            plot.title = element_text(hjust = 0.5))
    }
    pdf(paste0(outdir, "/", cur_mod, ".pdf"), width = 5, 
        height = 4)
    for (plot in plot_list) {
        print(plot)
    }
    dev.off()
}





#Plot enrichments by cluster
selected_terms <- subset(combined_output,db=="GO_Biological_Process_2023")
selected_terms <- subset(selected_terms, P.value < 0.05)



# subset selected terms
idx_top_1 <- match(unique(selected_terms$module),selected_terms$module)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


selected_terms$group <- factor(
  as.character(selected_terms$module),
  levels = levels(selected_terms$module)
)


# set max pval
quantile(-log(selected_terms$P.value), 0.95)
max_p <- 10

selected_terms$logp <- -log(selected_terms$P.value)
selected_terms$logp <- ifelse(selected_terms$logp > max_p, max_p, selected_terms$logp)

# remove GO Term ID
library(stringr)
selected_terms$Term <- str_replace(selected_terms$Term, "\\(GO.*", "")

selected_terms <- selected_terms %>%
  arrange(group)


selected_terms$wrap <- wrapText(selected_terms$Term, 35)

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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_by_cluster_GO.pdf'), width=6, height=8)
p
dev.off()


#######################################
#############  FIGURE S5D  ############
#######################################
pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'MyoFB_dot.pdf'), width=4.5, height=4)
DotPlot(M1,c('ACTA2','CDH11','TAGLN','SLIT3',
  'MINDY2','MYO1B','LIMS2','GARS'),group.by='group',dot.min=0,col.min=0,col.max=2) +
  coord_flip() + ylab('')+ xlab('')+
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  theme(
    panel.border = element_rect(size=1,fill=NA, color='black'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
) 
dev.off()


M1 <- AddModuleScore(M1, features=list(c('ACTA2','CDH11','TAGLN','SLIT3',
  'MINDY2','MYO1B','LIMS2','GARS')),assay="SCT",name="MyoFBScore")

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'MyoFB_Vln.pdf'), width=4.5, height=4)
VlnPlot(M1,'MyoFBScore1',split.by='group',pt.size=0)
dev.off()



#######################################
#############  FIGURE S5E  ############
#######################################

consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]

mapping <- labels2colors(1:100)

bulk_modules <- consensus_modules
bulk_modules$module <- match(consensus_modules$module,mapping)

seurat_ref <- M1

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


#Run enrichment by cell type
Idents(seurat_ref) <- "SubNames_Groups"
combined_set <- data.frame()
combined_output <- data.frame()

cell_types <- unique(seurat_ref$Subnames)
comparison <- list(c("RVF","NF"),c("RVF","pRV"))
for (j in cell_types){
  for (k in comparison){
    key_genes <- rownames(seurat_ref)

    gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes)
    
    gene_set<-subset(gene_set,p_val_adj<0.05)
    if (length(rownames(gene_set))==0){next}
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
          cur_df$celltype <- j
          cur_df$comparison <- paste0(k[1],'_',k[2])
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
          cur_df$celltype <- j
          cur_df$comparison <- paste0(k[1],'_',k[2])
          cur_df$direction <- 'up'
          combined_output <- rbind(combined_output, cur_df)
        }
    }
  }
}


#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_chea_terms_cell_type_up_RVF_vs_NF.pdf'), width=4, height=4)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_chea_terms_cell_type_down_RVF_vs_NF.pdf'), width=4, height=4)
p 
dev.off()



#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_chea_terms_cell_type_up_RVF_vs_pRV.pdf'), width=4, height=4)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_chea_terms_cell_type_down_RVF_vs_pRV.pdf'), width=4, height=4)
p 
dev.off()








#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_reactome_terms_cell_type_up_RVF_vs_NF.pdf'), width=8, height=7)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_reactome_terms_cell_type_down_RVF_vs_NF.pdf'), width=8, height=7)
p 
dev.off()



#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_reactome_terms_cell_type_up_RVF_vs_pRV.pdf'), width=8, height=7)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'FB_reactome_terms_cell_type_down_RVF_vs_pRV.pdf'), width=8, height=7)
p 
dev.off()

#######################################
#############  FIGURE S5F  ############
#######################################


M1 <- readRDS(file = "/Volumes/Extreme\ SSD/Final_Analysis/CellTypes/pc_sm_subclust.rds")


new.cluster.ids <- c("Sm","Pc","Pc","Sm")

names(new.cluster.ids) <- levels(M1)
M1 <- RenameIdents(M1, new.cluster.ids)

M1$Subnames <- M1@active.ident
M1$SubNames_Groups <- paste(M1$Subnames,M1$group,sep='_')

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PC_SM_snUMAP.pdf'), width=5, height=5)
PlotEmbedding(M1,group.by='Subnames',point_size=1,plot_under=TRUE,plot_theme=umap_theme()+NoLegend(),raster_dpi=400,raster_scale=0.5)
dev.off()


#######################################
#############  FIGURE S5G  ############
#######################################
#M5, M20, M11


consensus_modules <- read.csv("~/Downloads/hdWGCNA_TOM/bulk_heart_modules.csv")
consensus_modules <- consensus_modules[,1:3]
consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_ref))
# remove duplicate gene names
consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]

mapping <- labels2colors(1:100)

bulk_modules <- consensus_modules
bulk_modules$module <- match(consensus_modules$module,mapping)

seurat_ref <- M1

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


#Run enrichment by cell type
Idents(seurat_ref) <- "SubNames_Groups"
combined_set <- data.frame()
combined_output <- data.frame()

mods_idx <- c(5,11,20)
cell_types <- unique(seurat_ref$Subnames)
comparison <- list(c("RVF","NF"),c("RVF","pRV"))
for (i in mods_idx){
  for (j in cell_types){
    for (k in comparison){
      key_genes <- subset(bulk_modules,module %in% c(i))$gene_name
      key_genes <- key_genes[key_genes %in% rownames(seurat_ref)]

      gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes)
      
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


#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(5,11,20)])


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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_reactome_terms_cell_type_up_RVF_vs_NF.pdf'), width=8, height=7)
p / colorbar 
dev.off()


#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")
selected_terms <- subset(selected_terms,color %in% mapping[c(5,11,20)])


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


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_reactome_terms_cell_type_down_RVF_vs_NF.pdf'), width=8, height=7)
p / colorbar 
dev.off()








seurat_ref <- M1

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','Reactome_2022', 'ChEA_2022',"LINCS_L1000_Chem_Pert_up")


#Run enrichment by cell type
Idents(seurat_ref) <- "SubNames_Groups"
combined_set <- data.frame()
combined_output <- data.frame()

cell_types <- unique(seurat_ref$Subnames)
comparison <- list(c("RVF","NF"),c("RVF","pRV"))
for (j in cell_types){
  for (k in comparison){
    key_genes <- rownames(seurat_ref)

    gene_set <- FindMarkers(seurat_ref, ident.1 = paste0(j,"_",k[1]), ident.2 = paste0(j,"_",k[2]),features=key_genes)
    
    gene_set<-subset(gene_set,p_val_adj<0.05)
    if (length(rownames(gene_set))==0){next}
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
          cur_df$celltype <- j
          cur_df$comparison <- paste0(k[1],'_',k[2])
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
          cur_df$celltype <- j
          cur_df$comparison <- paste0(k[1],'_',k[2])
          cur_df$direction <- 'up'
          combined_output <- rbind(combined_output, cur_df)
        }
    }
  }
}


#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_chea_terms_cell_type_up_RVF_vs_NF.pdf'), width=4, height=4)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_chea_terms_cell_type_down_RVF_vs_NF.pdf'), width=4, height=4)
p 
dev.off()



#Up
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_chea_terms_cell_type_up_RVF_vs_pRV.pdf'), width=4, height=4)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="ChEA_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_chea_terms_cell_type_down_RVF_vs_pRV.pdf'), width=4, height=4)
p 
dev.off()








#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_reactome_terms_cell_type_up_RVF_vs_NF.pdf'), width=8, height=7)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_NF")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_reactome_terms_cell_type_down_RVF_vs_NF.pdf'), width=8, height=7)
p 
dev.off()



#Up
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="up")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_reactome_terms_cell_type_up_RVF_vs_pRV.pdf'), width=8, height=7)
p 
dev.off()

#Down
selected_terms <- subset(combined_output,db=="Reactome_2022")
selected_terms <- subset(selected_terms,direction=="down")
selected_terms <- subset(selected_terms,comparison=="RVF_pRV")

# subset selected terms
selected_terms <- subset(selected_terms, Adjusted.P.value < 0.05)
selected_terms$module_celltype <- paste0(selected_terms$module,'_',selected_terms$celltype)
idx_top_1<- match(unique(selected_terms$module_celltype),selected_terms$module_celltype)
idx_top_3 <- sort(c(idx_top_1,idx_top_1+1,idx_top_1+2))

selected_terms<-selected_terms[idx_top_3,]


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
    plot.margin = margin(0,0,0,0),
    panel.grid = element_line(size=0.25, color='lightgrey')
  )

pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_reactome_terms_cell_type_down_RVF_vs_pRV.pdf'), width=8, height=7)
p 
dev.off()

#######################################
#############  FIGURE S5H  ############
#######################################
a<-FindMarkers(seurat_ref,ident.1=c('Pc_RVF','Sm_RVF'),ident.2=c('Pc_pRV','Sm_pRV'))


pdf(paste0('~/Downloads/hdWGCNA_TOM/', 'PCSM_volcano_RVF_vs_pRV.pdf'), width=4, height=8)

EnhancedVolcano(a,lab=rownames(a),
  x='avg_log2FC',y='p_val_adj',
  FCcutoff = 0.1,pCutoff=0.05,xlim=c(-6,6),ylim=c(0,27))
dev.off()


