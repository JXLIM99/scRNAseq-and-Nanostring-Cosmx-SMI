#R script for Fig5, FigS5 and Supplementary Table 2
#Brief Experiment Background: 
#Rag1-/- mice were subcutaneously injected with B16F10 cells. 
#Tumor were harvested when tumor size reached around 800mm3 and subjected to Nanostring Cosmx Analysis
#######################################################################################################

#Libraries
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(reshape2)
library(BiocManager)
library(openxlsx)
library(dplyr)
library(HGNChelper)
library(biomaRt)
library(scRNAtoolVis)

#Section 1: Data Import
Rag.seurat <- LoadNanostring(data.dir = "C:/Users/User/Desktop/spatial/30 May 2024/Data Import", #address to raw data files
                           fov = "tumorfov", assay = "Nanostring")

#Section 2: Normalisation and Clustering
Rag.seurat <- SCTransform(Rag.seurat, assay = "Nanostring",
                                         clip.range = c(-10,10), verbose = FALSE)
Rag.seurat <- RunPCA(Rag.seurat, npcs=5)
Rag.seurat <- RunUMAP(Rag.seurat, dims = 1:5)
Rag.seurat <- FindNeighbors(Rag.seurat, reduction = "pca", dims = 1:5)
Rag.seurat <- FindClusters(Rag.seurat, 
                          resolution = c(0.2),
                          verbose = F)

library(clustree)
clustree(Rag.seurat, prefix = "SCT_snn_res.")
Idents(object = Rag.seurat) <- "SCT_snn_res.0.2"

#Section 3: Marker Identification
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list <-  gene_sets_prepare("CosMx_Ref_Matrix.xlsx") #available in Supplementary Table 2a
es.max <-  sctype_score(scRNAseqData = Rag.seurat@assays[["SCT"]]@scale.data, 
                        scaled = TRUE, 
                        gs = gs_list$gs_positive, 
                        gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(Rag.seurat@meta.data$SCT_snn_res.0.2), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Rag.seurat@meta.data[Rag.seurat@meta.data$SCT_snn_res.0.2==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Rag.seurat@meta.data$Rag.seurat==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

Rag.seurat@meta.data$cellidentity = ""

for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  Rag.seurat@meta.data$cellidentity[Rag.seurat@meta.data$SCT_snn_res.0.2== j] = as.character(cl_type$type[1])
}

#Merged Umap (Fig5a)
Idents(Rag.seurat)<- "cellidentity"
colours_frequency <- c("#FFF085","#FCA917","#F65D63","#FF4D96","#C75EAC","#98217C","#601652")  
DimPlot(Rag.seurat, reduction = "umap", 
        label = F, repel = TRUE, pt.size = 0.5,
        cols = colours_frequency)

#Umap of Rag1ko reconstituted with either WT Teff:WT Treg or WT Teff:PD1KO Treg (Fig5a)
DimPlot(Rag.seurat, reduction = "umap", 
        label = F, repel = TRUE, pt.size = 0.5,
        cols = colours_frequency,
        split.by = "Group")

#Since Naive CD4+T cells containing heterogenous population of immune cells hence was renamed to Immune Cells
Rag.seurat@meta.data$cellidentity[Rag.seurat@meta.data$cellidentity == 'Naive CD4+ T cells'] <- 'Immune cells'

#Frequency of 1st clustering (Fig5b)
table <- table(Rag.seurat@meta.data$Group, Rag.seurat@meta.data$cellidentity)
Frequency <- as.data.frame(table)
colnames(Frequency) <- c("Group", "CellIdentity", "Percentage")
colours_frequency <- c("#FFF085","#FCA917","#F65D63","#98217C","#601652")  
Frequency$CellIdentity <- factor(Frequency$CellIdentity, 
                                  levels = c("Basal cell", "Immune cells", "Schwann cell", "Melanoma tumour cell", "Melanocyte"))

ggplot(Frequency, aes(fill=CellIdentity, y=Percentage, x=Group)) + 
geom_bar(position="fill", stat="identity") +
scale_fill_manual(values=colours_frequency) +
theme(axis.text.x = element_text(angle = 0, hjust = 1),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(size = 0.5))

#Sub-clustering of Immune cells 
Idents(Rag.seurat) <- "cellidentity"
IC.seurat <- subset(Rag.seurat, idents = "Immune cells")
#repeat previous steps with IC.seurat starting from SCTransform in Section 2 till cell identity allocation in Section 3; npcs 20 was used in function "RunPCA" and resolution 1.6 were used in function "FindCluster".

#Redefine Foxp3+ in immune cell cluster as Treg
grep('Foxp3', rownames(IC.seurat@assays[["SCT"]]$counts))
length(which(IC.seurat@assays[["SCT"]]$counts[rownames(IC.seurat@assays[["SCT"]])[378], ] != 0)) 
IC.seurat@meta.data$cellidentity[which(IC.seurat@assays[["SCT"]]$counts[rownames(IC.seurat@assays[["SCT"]])[378], ] != 0)] <- 'Treg'

#Frequency of subclustering (Fig5c)
table <- table(IC.seurat@meta.data$Group, IC.seurat@meta.data$cellidentity)
Frequency <- as.data.frame(table)
colnames(Frequency) <- c("Group", "CellIdentity", "Percentage")
FrequencyCell.Identity <- factor(Frequency$CellIdentity, 
                                  levels = c("Monocyte", "DC", "ILC", "Mast cell", "Treg", 
                                             "NK", "Macrophage", "CD4+ T cell"))
colours_frequencyIC <- c("#c34bff","#dc4dff","#C5AAF5","#A3CBF1","#6e73f7","#3F59A7","#182D90","#79BFA1")

ggplot(Frequency, aes(fill=CellIdentity, y=Percentage, x=Group)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=colours_frequencyIC) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5))

#Remap subcluster's identity to original seurat object
Idents(IC.seurat) <- "cellidentity"
Idents(Rag.seurat) <- "cellidentity"
Rag.seurat$cellidentity <- as.character(Idents(Rag.seurat))
Rag.seurat$cellidentity[Cells(IC.seurat)] <- paste(Idents(IC.seurat))

#Dotplot of WT Treg and Pd1ko Treg (Fig5d-e, S5c-d)
Rag.seurat.Treg <- subset(Rag.seurat, idents = "Treg")
Rag.seurat.Treg$group.celltype <- paste(Rag.seurat.Treg$Group,Rag.seurat.Treg$cellidentity, sep = "_")
Idents(Rag.seurat.Treg) <- "group.celltype"

jjDotPlot(integrate.filtered.rna,
          gene = c("Pdcd1","Cd274", "Tigit", "Ctla4", "Tnfrsf18", "Havcr2"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'group.celltype')

jjDotPlot(integrate.filtered.rna,
          gene = c("Tpi1","Eno1/b", "Apoe", "Adipoq", "Pparg"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'group.celltype')

jjDotPlot(integrate.filtered.rna,
          gene = c("Ccl11", "Ccl17", "Ccl19",	"Ccl2",	"Ccl20",	"Ccl21a/b/d",	"Ccl22",	"Ccl26",	"Ccl28",	"Ccl3",
                   "Ccl4",	"Ccl5",	"Ccl8",	"Ccl9","Cx3cl1","Cxcl1",	"Cxcl10",	"Cxcl12",	"Cxcl13",	"Cxcl14",	"Cxcl16",	"Cxcl17",	"Cxcl2",	
                   "Cxcl3",	"Cxcl5",	"Cxcl9"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'group.celltype')

jjDotPlot(integrate.filtered.rna,
          gene = c("Ccr1",	"Ccr10",	"Ccr2",	"Ccr5",	"Ccr7", "Cx3cr1","Cxcr1",	"Cxcr3",	"Cxcr4",	"Cxcr5",	"Cxcr6"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'group.celltype')

#Violin plot of each cell cluster (Fig5i-l)
Idents(Rag.seurat) <- "group.celltype"
VlnPlot(object = Rag.seurat, 
        features = c("Ifng"), pt.size = 0.5, #also includes Fig.5i:Areg, Il2; Fig.5j: Il10, Osm, Reg1; Fig5k Adgrf5; Fig.5l: Gzmk, Gzmb, Csk
        log= T, adjust = 1, alpha = 1, 
        idents = c("RagWT_CD4+ T cell", "RagKO_CD4+ T cell"), #also includes Fig.5j: DC; Fig.5k: Macrophage,Fig.5l: NK
        fill.by = "feature") + theme(axis.text.x = element_text(angle = 90)) + 
  stat_compare_means(comparisons = list(c("RagWT_CD4+ T cell", "RagKO_CD4+ T cell")), method = "wilcox.test")

#Violin plot of Tnfsf8 expression in all clusters (FigS5f)
Idents(Rag.seurat) <- "Group"
RagWT.seurat <- subset (Rag.seurat, idents = "RagWT")
RagKO.seurat <- subset (Rag.seurat, idents = "RagKO")

Idents(RagWT.seurat) <- "cellidentity"
Idents(RagKO.seurat) <- "cellidentity"

colours_frequencyALL <- c("#FFF085","#F65D63","#98217C","#601652","#c34bff",
                          "#dc4dff","#C5AAF5","#A3CBF1","#6e73f7","#3F59A7","#182D90","#3F8377","#79BFA1","#ABDECC")

VlnPlot(object = RagWT.seurat, 
        features = c("Tnfsf8"), pt.size = 0.5,
        log= T, adjust = 30, alpha = 1,
        cols = colours_frequencyALL,
        y.max = 5, same.y.lims = T,
        fill.by = "feature") + theme(axis.text.x = element_text(angle = 90))

VlnPlot(object = RagKO.seurat, 
        features = c("Tnfsf8"), pt.size = 0.5,
        log= T, adjust = 30, alpha = 1,
        cols = colours_frequencyALL,
        y.max = 5, same.y.lims = T,
        fill.by = "feature") + theme(axis.text.x = element_text(angle = 90))

#Section 4: Enrichment Analysis (Fig5f)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(biomaRt)

Idents(Rag.seurat) <- Rag.seurat@meta.data$group.celltype
DE_Treg <- FindMarkers(Rag.seurat, 
                        ident.1 = "RagKO_Treg", 
                        ident.2 = "RagWT_Treg", 
                        verbose = FALSE)

mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
results <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                 filters = "external_gene_name", 
                 values = DE_Treg$gene,
                 mart = mart)
results %>% head

DE_Treg <- DE_Treg %>% 
  left_join(., results, by = c("gene" = "external_gene_name")) %>% 
  filter(!is.na(entrezgene_id)) %>% 
  filter(!is.na(p_val_adj))

DE_Treg <- dplyr::filter(DE_Treg, p_val < 0.05, avg_log2FC > 0.585)

ego<- enrichGO(gene = DE_Treg$entrezgene_id,
                      keyType = "ENTREZID",
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

library(enrichplot)
options(enrichplot.colours = c("red","blue"))
dotplot(ego, showCategory=5) #available in Supplementary Table 2b

#Section 5: FOV Analysis (FigS5b)
Idents(Rag.seurat) <- "FOV"
ImageDimPlot(Rag.seurat, fov = "tumorfov", axes = T) #to locate coordinates

crop.FOV20 <- Crop(Rag.seurat[["tumorfov"]], 
                     x=c(62500,67000), y=c(30500,34500)) 
Rag.seurat[["zoomFOV20"]] <- crop.FOV20
DefaultBoundary(Rag.seurat[["zoomFOV20"]])<- "segmentation"

Idents(Rag.seurat) <- Rag.seurat@meta.data$cellidentity
Rag.seurat <- SetIdent(Rag.seurat, 
                       value = factor(Idents(Rag.seurat), 
                               levels = c("Basal cell", "Schwann cell", "Melanocyte", "Melanoma Tumour Cell", "Monocyte", "DC", "ILC", "Mast cell", "Treg", 
                               "NK", "Macrophage", "Lymphatic endothelial cell", "CD4+ T cell", "Endothelial")))

colours_mol <- c("red", "salmon", "darksalmon", "orange", "brown", "lightgreen", "green", "darkgreen", "lightblue", "blue", "darkblue", "purple", "violet","pink","grey", "darkgrey")

ImageDimPlot(Rag.seurat, axes = TRUE, fov = "zoomFOV20", 
             mols.size = 1, 
             border.color = "black",
             dark.background = F,
             cols = colours_frequencyALL,
             molecules = c("Cd68", "Foxp3", "Cd274", "Entpd1", "Mcpt8", "Itgae", "Cd81",
                           "Mrc1", "Msr1", "Eomes", "Itgax", "Klrb1", "Notch2", "Cd3d", "Cd4"),
             mols.cols = colours_mol,
             coord.fixed = T)

#Section 6: Analysis of Treg through physical interaction on FOVs
#Add CellID to metadata
Rag.seurat@meta.data$CellID <- paste(Rag.seurat@meta.data$FOV, Rag.seurat@meta.data$CellID, sep = "_")

#Re-distribute Treg into different groups with cell ID (Supplementary Table 2c) through visualisation of physical interaction on FOV
TregInteraction <- read.xlsx("All Treg Spatial Interaction.xlsx") #available in Supplementary Table 2d
Rag.seurat.Treg@meta.data$Interaction <- TregInteraction$Interaction

#1.Treg Metabolic and Oxygenation State (Fig5m)
Idents(Rag.seurat.Treg) <- Rag.seurat.Treg@meta.data$Interaction
Rag.seurat.Treg$group.interaction <- paste(Rag.seurat.Treg$Group,Rag.seurat.Treg$Interaction, sep = "_")
Idents(Rag.seurat.Treg) <- Rag.seurat.Treg@meta.data$group.interaction

Tregoxg <- DotPlot(Rag.seurat.Treg, group.by = "group.interaction", 
                   features = c("Hif1a", "Cox4i2", "Tpi1", "Eno1/b"))
Tregoxg <- Tregoxg$data
write.xlsx(Tregoxg_df, "Tregoxg_df.xlsx", rowNames = T, colNames = T) 

Tregoxg_annotation <- read.xlsx("Tregoxggroup.xlsx")
Tregoxggroup <- Tregoxg_annotation$Interaction

Tregoxg_heatmap <- read.xlsx("Tregoxg_input.xlsx")
rownames(Tregoxg_heatmap)<- Tregoxg_heatmap$Gene
Tregoxg_heatmap$Gene <- NULL

group_colors <- list(Group = c("WT Teff: WT Treg" = "#00BFC4", "WT Teff: Pd1-/- Treg" = "#F8766D"))
Tregoxyg_annotation <- data.frame(Group = Tregoxg_annotation$Group)
rownames(Tregoxyg_annotation) <- Tregoxg_annotation$Interaction
Tregoxg_heatmap[Tregoxg_heatmap < 0] <- 0

pheatmap(Tregoxg_heatmap, cluster_rows = F, cluster_cols = F, scale = "none",
         color = colorRampPalette(c("white","Salmon"))(100), 
         main = "Treg Oxygenation State", border_color = "black", 
         cellwidth = 10, cellheight = 10,
         annotation_col = Tregoxyg_annotation,
         annotation_colors = group_colors)

library(pheatmap)
#2.Short range and Long range heatmap (WT) (Fig.S5e)
Idents(Rag.seurat.Treg) <- "Group"
Rag.seurat.Treg.WT <- subset(Rag.seurat.Treg, idents = "RagWT")
Idents(Rag.seurat.Treg) <- "Interaction"
WTTregInteraction.marker <- FindAllMarkers(Rag.seurat.Treg.WT, only.pos = T) 
write.xlsx(WTTregInteraction.marker, "WTSRLRgene.xlsx", rowNames = T, colNames = T) #available in Supplementary Table 2e

WTannotation_SRLR <- read.xlsx("WTSRLRgene.xlsx") 
WTgene_SRLR <- WTannotation_SRLR$Gene

Idents(Rag.seurat.Treg.WT) <- Rag.seurat.Treg.WT@meta.data$Interaction
WTSRLR <- DotPlot(Rag.seurat.Treg.WT, group.by = "Interaction",
                    features = WTgene_SRLR)
WTSRLR_df <- WTSRLR$data
write.xlsx(WTSRLR_df, "WTSRLR_df.xlsx", rowNames = T, colNames = T)

WTSRLR_heatmap <- read.xlsx("WTSRLR_input.xlsx") 
rownames(WTSRLR_heatmap)<- WTSRLR_heatmap$Gene
WTSRLR_heatmap$Gene <- NULL

WTSRLR_colors <- list(Signaling = c("SR" = "#C7B1DF", "LR" = "#876F98"))
WTSRLR_annotation <- data.frame(Signaling = WTannotation_SRLR$Signaling)
rownames(WTSRLR_annotation) <- WTannotation_SRLR$Gene
WTSRLR_heatmap[WTSRLR_heatmap < 0] <- 0

pheatmap(WTSRLR_heatmap, cluster_rows = F, cluster_cols = F, scale = "none",
         color = colorRampPalette(c("#141743", "lightblue", "white", "salmon", "darkred"))(1000),
         main = "WT Treg Communication", border_color = "black", 
         cellwidth = 15, cellheight = 10,
         annotation_row = WTSRLR_annotation,
         annotation_colors = WTSRLR_colors)

#Short range and Long range heatmap (PD1KO) (Fig.S5f)
Idents(Rag.seurat.Treg) <- "Group"
Rag.seurat.Treg.PD1KO <- subset(Rag.seurat.Treg, idents = "RagKO")
Idents(Rag.seurat.Treg) <- "Interaction"
PD1KOTregInteraction.marker <- FindAllMarkers(Rag.seurat.Treg.PD1KO, only.pos = T)
write.xlsx(PD1KOTregInteraction.marker, "KOSRLRgene.xlsx", rowNames = T, colNames = T) #available in Supplementary Table 2f

PD1KOannotation_SRLR <- read.xlsx("KOSRLRgene.xlsx")
PD1KOgene_SRLR <- PD1KOannotation_SRLR$Gene

Idents(Rag.seurat.Treg.PD1KO) <- Rag.seurat.Treg.PD1KO@meta.data$Interaction
PD1KOSRLR <- DotPlot(Rag.seurat.Treg.PD1KO, group.by = "Interaction",
                  features = PD1KOgene_SRLR)
PD1KOSRLR_df <- PD1KOSRLR$data
write.xlsx(PD1KOSRLR_df, "PD1KOSRLR_df.xlsx", rowNames = T, colNames = T)

PD1KOSRLR_heatmap <- read.xlsx("PD1KOSRLR_input.xlsx") #available in Supplementary Table 2f
rownames(PD1KOSRLR_heatmap )<- PD1KOSRLR_heatmap$Gene
PD1KOSRLR_heatmap$Gene <- NULL

SRLR_colors <- list(Signaling = c("SR" = "#C7B1DF", "LR" = "#876F98"))
PD1KOSRLR_annotation <- data.frame(Signaling = PD1KOannotation_SRLR$Signaling)
rownames(PD1KOSRLR_annotation) <- PD1KOannotation_SRLR$Gene
PD1KOSRLR_heatmap[PD1KOSRLR_heatmap < 0] <- 0

pheatmap(PD1KOSRLR_heatmap, cluster_rows = F, cluster_cols = F, scale = "none",
         color = colorRampPalette(c("#141743", "lightblue", "white", "salmon", "darkred"))(1000),
         main = "Pd1-/- Treg Communication", border_color = "black", 
         cellwidth = 15, cellheight = 10,
         annotation_row = PD1KOSRLR_annotation,
         annotation_colors = SRLR_colors)
