#R script for Fig5, FigS7 and Supplementary Table 4
#Brief Experiment Background: 
#Rag1-/- mice were subcutaneously injected with B16F10 cells.
#Rag1-/- mice were reconstituted with WT Teff along with either WT Treg or Pd1-/- Treg.
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
Rag.seurat <- LoadNanostring(data.dir = "C:/Users/User/Desktop/spatial/Data Import", #address to raw data files
                           fov = "tumorfov", assay = "Nanostring")

#Section 2: Normalization and Clustering
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
gs_list <-  gene_sets_prepare("CosMx_Ref_Matrix.xlsx") #available in Supplementary Table 4a
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

#Merged Umap (Fig5a, Left panel)
Idents(Rag.seurat)<- "cellidentity"
colours_frequency <- c("#FFF085","#FCA917","#F65D63","#FF4D96","#C75EAC","#98217C","#601652")  
DimPlot(Rag.seurat, reduction = "umap", 
        label = F, repel = TRUE, pt.size = 0.5,
        cols = colours_frequency)

#Umap of Rag1ko reconstituted with either WT Teff:WT Treg or WT Teff:PD1KO Treg (Fig5a, Middle and right panel)
DimPlot(Rag.seurat, reduction = "umap", 
        label = F, repel = TRUE, pt.size = 0.5,
        cols = colours_frequency,
        split.by = "Group")

#Since Naive CD4+T cells containing heterogenous population of immune cells hence was renamed to Immune Cells
Rag.seurat@meta.data$cellidentity[Rag.seurat@meta.data$cellidentity == 'Naive CD4+ T cells'] <- 'Immune cells'

#Frequency of 1st clustering (Fig5b)
table <- table(Rag.seurat@meta.data$Group, Rag.seurat@meta.data$cellidentity)
Frequency <- as.data.frame(table)
colnames(Frequency) <- c("Group", "CellIdentity", "Frequency")

colours_frequency <- c("#FFF085","#FCA917","#F65D63","#98217C","#601652")  
Frequency$CellIdentity <- factor(Frequency$CellIdentity, 
                                  levels = c("Basal cell", "Immune cells", "Schwann cell", "Tumor Cell", "Melanocyte"))

Frequency$Group <- factor(Frequency$Group, 
                                 levels = c("WT", "KO"))

ggplot(Frequency, aes(fill=CellIdentity, y=Frequency, x=Group)) + 
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
Immunecells <- subset(IC.seurat, idents = c("CD4+ T cell", "DC" , "ILC",
                                            "Macrophage", "Mast cell",  "Monocyte", "NK", "Treg"))
table <- table(Immunecells@meta.data$Group, Immunecells@meta.data$cellidentity)
Frequency <- as.data.frame(table)
colnames(Frequency) <- c("Group", "CellIdentity", "Frequency")
Frequency$CellIdentity <- factor(Frequency$CellIdentity, 
                                  levels = c("Monocyte", "DC", "ILC", "Mast cell", "Treg", 
                                             "NK", "Macrophage", "CD4+ T cell"))
colours_frequencyIC <- c("#c34bff","#dc4dff","#C5AAF5","#A3CBF1","#6e73f7","#3F59A7","#182D90","#79BFA1")
Frequency$Group <- factor(Frequency$Group, 
                          levels = c("WT", "KO"))

ggplot(Frequency, aes(fill=CellIdentity, y=Frequency, x=Group)) + 
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

#key markers for cell identity (Fig.S7c)
jjDotPlot(Rag.seurat,
          gene = c("Cd4",
                   "Eomes", "Foxp3", "Tbx21", "Gata3", "Rora", "Irf4", "Maf",
                   "Itgax", "Itgam", "Klrb1", "Cd163", "Cd68","Cd14", "Col1a1", "Col1a2",
                   "Cd80", "Cd86", "Gpnmb", "Apoe", "Becn1", "Rxra",
                   "Ctla4", "Entpd1", "Havcr2", "Lag3", "Pdcd1", "Tigit",
                   "Cd28", "Cd40", "Icos", "Tnfrsf1b", "Tnfrsf4", "Tnfrsf9", "Tnfrsf13b", "Tnfrsf17", "Tnfrsf18",
                   "Ccr7", "Sell", "Tcf7",
                   "Mki67", "Pcna", 
                   "Il1a", "Il1b", "Il2", "Il6", "Il12a", "Il12b", "Il15", "Il17a", "Il18", "Ifng", "Tgfb1", "Tgfb3",
                   "Gzma", "Gzmb", "Gzmk", "Nkg7", "Prf1", "S100b", "Cd34", "Vcam1"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'cellidentity')

#Violin plot of Il10 in DC and Gzmb in NK with clustered Wilcoxon rank-sum test applied (Fig5f-g)
Rag.seurat$group.celltype <- paste(Rag.seurat$Group,Rag.seurat$cellidentity, sep = "_")
Idents(Rag.seurat) <- "group.celltype"
VlnPlot(object = Rag.seurat, 
        features = c("Il10"), pt.size = 0.5,
        log= F, adjust = 20, alpha = 1, 
        idents = c("RagWT_DC", "RagKO_DC"), #also visualized Gzmb for RagWT_NK vs RagKO_NK
        fill.by = "feature") + theme(axis.text.x = element_text(angle = 90)) 

library(clusrank)
library(dplyr)
DC <- subset(Rag.seurat, idents = "DC")
gene <- "Il10" 
expression_data <- FetchData(object = DC, vars = c("Group", "Sample", gene))
clusWilcox.test(Il10 ~ Group + cluster(Sample), data = expression_data, method = "ds")

NK <- subset(Rag.seurat, idents = "NK")
gene <- "Gzmb" 
expression_data <- FetchData(object = NK, vars = c("Group", "Sample", gene))
clusWilcox.test(Gzmb ~ Group + cluster(Sample), data = expression_data, method = "ds")

#Violin plot of Tnfsf8 (CD30L) expression in all clusters (FigS7e, Supplemental Table 4f)
#with clustered Wilcoxon rank-sum test corrected by FDR applied
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

gene <- "Tnfsf8"
expression_data <- FetchData(object = RagWT.seurat, vars = c("cellidentity", "Sample", gene))
with(expression_data,pairwise.clusWilcox.test(
  expression_data$Tnfsf8,
  group = cellidentity,
  cluster = Sample,
  method = "ds",
  p.adjust.method = "fdr"))

VlnPlot(object = RagKO.seurat, 
        features = c("Tnfsf8"), pt.size = 0.5,
        log= T, adjust = 30, alpha = 1,
        cols = colours_frequencyALL,
        y.max = 5, same.y.lims = T,
        fill.by = "feature") + theme(axis.text.x = element_text(angle = 90)) 

gene <- "Tnfsf8"
expression_data <- FetchData(object = RagKO.seurat, vars = c("cellidentity", "Sample", gene))
with(expression_data,pairwise.clusWilcox.test(
  expression_data$Tnfsf8,
  group = cellidentity,
  cluster = Sample,
  method = "ds",
  p.adjust.method = "fdr"))

#Section 4: Enrichment Analysis (FigS7d)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(biomaRt)

Rag.seurat$group.celltype <- paste(Rag.seurat$Group, Rag.seurat$cellidentity, sep = "_")
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
dotplot(ego, showCategory=5)

#Section 5: FOV Analysis (FigS7b, 5j)
Idents(Rag.seurat) <- "FOV"
ImageDimPlot(Rag.seurat, fov = "tumorfov", axes = T) #to locate coordinates

crop.FOV20 <- Crop(Rag.seurat[["tumorfov"]], 
                     x=c(62500,67000), y=c(30500,34500)) 
Rag.seurat[["zoomFOV20"]] <- crop.FOV20
DefaultBoundary(Rag.seurat[["zoomFOV20"]])<- "segmentation"

Idents(Rag.seurat) <- Rag.seurat@meta.data$cellidentity
Rag.seurat <- SetIdent(Rag.seurat, 
                       value = factor(Idents(Rag.seurat), 
                               levels = c("Basal cell", "Schwann cell", "Melanocyte", "Tumor Cell", "Monocyte", "DC", "ILC", "Mast cell", "Treg", 
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

#Re-distribute Treg into different groups with cell ID (Supplementary Table 4c) through visualization of physical interaction on FOV
TregInteraction <- read.xlsx("All Treg Spatial Interaction.xlsx") #available in Supplementary Table 4d
Rag.seurat.Treg <- 
Rag.seurat.Treg@meta.data$Interaction <- TregInteraction$Interaction

#2.Short range and Long range (Fig.5e, Supplementary Table 4d,e)
Idents(Rag.seurat.Treg) <- "Group"
Rag.seurat.Treg.WT <- subset(Rag.seurat.Treg, idents = "RagWT")
Rag.seurat.Treg.PD1KO <- subset(Rag.seurat.Treg, idents = "RagKO")                             
DEWTgene_SRLR <- FindMarkers(Rag.seurat.Treg.WT, only.pos = T) #Supplementary Table 4d
DEKOgene_SRLR <- FindMarkers(Rag.seurat.Treg.PD1KO, only.pos = T) #Supplementary Table 4e
