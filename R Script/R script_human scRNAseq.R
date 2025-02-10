#R script for Fig6, FigS8 and Supplementary Table 6
#Brief Experiment Background: 
#n=3 healthy controls and n=3 stage IV melanoma patients 
#PBMCs were subjected to BD Rhapsody single-cell RNA seq.
###############################################################################
#packages 
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("Matrix")
install.packages("RCurl")
install.packages("scales")
install.packages("data.table")
install.packages("readxl")
install.packages("BiocManager")
install.packages("ggpubr")
install.packages("Seurat")

BiocManager::install("ensembldb")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("biomaRt")
BiocManager::install("enrichplot")
BiocManager::install("AnnotationHub")

library(BiocManager)
library(Matrix)
library(RCurl)
library(readxl)
library(scales)
library(ggpubr)
library(AnnotationHub)
library(ensembldb)
library(HGNChelper)
library(clusterProfiler)
library(AnnotationDbi)
library(biomaRt)
library(enrichplot)
library(ggplot2)

#Section 1 data import
library(tidyverse)
library(Seurat)
library(data.table)

-expMat_C1 <- ReadMtx(mtx = "C1 matrix.mtx", cells = "C1 barcodes.tsv", features = "C1 features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)

expMat_C1[1:5,1:5]

rna_C1 <- Seurat::CreateSeuratObject(counts = expMat_C1,
                                     min.cells = 0, min.features = 0, assay = "RNA")
rna_C1@meta.data %>% head()

smk_C1 <- fread(file = "C1 sample tag.csv", sep = ",", header = TRUE) %>%
  data.frame(row.names = 1)

smk_C1[1:5,]

rna_C1 <- AddMetaData(object = rna_C1, metadata = smk_C1)
rna_C1@meta.data %>% head()

expMat_C2 <- ReadMtx(mtx = "C2 matrix.mtx", cells = "C2 barcodes.tsv", features = "C2 features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)
expMat_C2[1:5,1:5]
rna_C2 <- Seurat::CreateSeuratObject(counts = expMat_C2,
                                     min.cells = 0, min.features = 0, assay = "RNA")
rna_C2@meta.data %>% head()
smk_C2 <- fread(file = "C2 sample tag.csv", sep = ",", header = TRUE) %>%
  data.frame(row.names = 1)
smk_C2[1:5,]

rna_C2 <- AddMetaData(object = rna_C2, metadata = smk_C2)
rna_C2@meta.data %>% head()

expMat_C3 <- ReadMtx(mtx = "C3 matrix.mtx", cells = "C3 barcodes.tsv", features = "C3 features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)

expMat_C3[1:5,1:5]
rna_C3 <- Seurat::CreateSeuratObject(counts = expMat_C3,
                                     min.cells = 0, min.features = 0, assay = "RNA")
rna_C3@meta.data %>% head()

smk_C3 <- fread(file = "C3 sample tag.csv", sep = ",", header = TRUE) %>%
  data.frame(row.names = 1)

smk_C3[1:5,]

rna_C3 <- AddMetaData(object = rna_C3, metadata = smk_C3)
rna_C3@meta.data %>% head()

Idents(rna_C2) <- rna_C2@meta.data$Sample_Name
Healthy_C2 <- subset(rna_C2, idents = "HC2")
Melanoma_C2 <- subset(rna_C2, idents = "MC2")

Healthy_C3$Sample_Name <- "HC"
Healthy_C2$Sample_Name <- "HC"
Healthy_C1$Sample_Name <- "HC"
Melanoma_C1$Sample_Name <- "MC"
Melanoma_C2$Sample_Name <- "MC"
Melanoma_C3$Sample_Name <- "MC"

MCHC <- merge(Healthy_C1, y = c(Healthy_C2, Healthy_C3, Melanoma_C1, Melanoma_C2, Melanoma_C3),
              add.cell.ids = c("HC1", "HC2", "HC3", "MC1", "MC2", "MC3"))

#Section 2: Quality Control
MCHC@meta.data %>%
  group_by(Sample_Name) %>%
  tally() %>%
  mutate(pct = n/sum(n)) %>%
  ggplot(aes(x=Sample_Name, y=n, fill= Sample_Name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = n, 
                label= paste0(n, "\n", scales::percent(pct)),
                vjust= -0.5, size= 2.5), show.legend = FALSE) +
  theme_classic()

#Remove multiplet and undetermines
MCHC <- subset(MCHC, subset = Sample_Name %in% c("Multiplet", "Undetermined"), invert = T)

# Quality control is to filter out cell labels that are of low quality, e.g. dead cells
MCHC@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nCount_RNA, fill= Sample_Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  ggtitle("nCount_RNA")

MCHC@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nFeature_RNA, fill= Sample_Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  ggtitle("nFeature_RNA")

MCHC@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500, colour = "red") +
  geom_hline(yintercept = 300, colour = "red") +
  facet_wrap(~Sample_Name) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Seurat object 
filtered.MCHC <- subset(x = MCHC, subset = (nCount_RNA >= 1200) & (nCount_RNA <= 60000) & (nFeature_RNA >= 50) & (nFeature_RNA <= 200))
counts <- GetAssayData(JoinLayers(object = filtered.MCHC, slot = "counts"))
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 1
filtered_counts <- counts[keep_genes, ]

filtered.MCHC <- CreateSeuratObject(filtered_counts,meta.data = filtered.MCHC@meta.data)
filtered.MCHC@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nCount_RNA, fill= Sample_Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  ggtitle("nCount_RNA")

filtered.MCHC@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nFeature_RNA, fill= Sample_Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  ggtitle("nFeature_RNA")

filtered.MCHC@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500, colour = "red") +
  geom_hline(yintercept = 300, colour = "red") +
  facet_wrap(~Sample_Name) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Cell Number 
MCHC@meta.data %>% 
  group_by(Sample_Name) %>% 
  dplyr::summarise(cell_number = length(Sample_Name))

filtered.MCHC@meta.data %>% 
  group_by(Sample_Name) %>% 
  dplyr::summarise(cell_number = length(Sample_Name))

#Normalize
filtered.MCHC <- NormalizeData(filtered.MCHC, 
                               normalization.method = "LogNormalize")

# split the dataset into a list of two seurat objects (hc and mc)
split.filtered.MCHC <- SplitObject(filtered.MCHC, split.by = "Sample_Name")

# normalize and identify variable features for each dataset independently
split.filtered.MCHC <- lapply(X = split.filtered.MCHC, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Select features that are repeatedly variable across datasets for integration
integ.features <- SelectIntegrationFeatures(object.list = split.filtered.MCHC, 
                                            nfeatures = 2000) 

# Select the most variable features to use for integration
anchors <- FindIntegrationAnchors(object.list = split.filtered.MCHC,
                                  anchor.features = integ.features)


integrate.filtered.MCHC <- IntegrateData(anchorset = anchors)
integrate.filtered.MCHC <- ScaleData(integrate.filtered.MCHC, verbose = FALSE)
integrate.filtered.MCHC <- RunPCA(integrate.filtered.MCHC, npcs = 50, verbose = FALSE)

ElbowPlot(integrate.filtered.MCHC, ndims = 50)
integrate.filtered.MCHC <- RunUMAP(integrate.filtered.MCHC, 
                                   reduction = "pca", 
                                   dims = 1:20) 

integrate.filtered.MCHC <- FindNeighbors(integrate.filtered.MCHC, 
                                         reduction = "pca", 
                                         dims = 1:20)

DimPlot(integrate.filtered.MCHC, group.by = "Sample_Name")

#Section 3: Clustering
# Determine the K-nearest neighbor graph
integrate.filtered.MCHC <- FindNeighbors(object = integrate.filtered.MCHC, dims = 1:20)
# Determine the clusters for various resolutions  
integrate.filtered.MCHC<- FindClusters(object = integrate.filtered.MCHC,
                                       resolution = c(0.2, 0.4, 0.6,  0.8, 1, 1.2, 1.4),
                                       verbose = F)
integrate.filtered.MCHC@meta.data %>% 
  dplyr::select(contains("integrat")) %>% 
  map_int(~ unique(.x) %>% length)

install.packages("clustree")
library(clustree)

clustree(integrate.filtered.MCHC, prefix = "integrated_snn_res.")
Idents(object = integrate.filtered.MCHC) <- "integrated_snn_res.0.6"

#UMAP visualization
DimPlot(integrate.filtered.MCHC,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

#Section 4: Marker Identification
library(tidyverse)
library(Seurat)
library(data.table)
library(ggpubr)
library(HGNChelper)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list <-  gene_sets_prepare("/Users/teganmctaggart/Desktop/unzipped rhapsody files/ScTypeDB_short.xlsx", "Immune system")
es.max <-  sctype_score(scRNAseqData = integrate.filtered.MCHC@assays[["integrated"]]@scale.data, 
                        scaled = TRUE, 
                        gs = gs_list$gs_positive, 
                        gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(integrate.filtered.MCHC@meta.data$integrated_snn_res.0.6), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(integrate.filtered.MCHC@meta.data[integrate.filtered.MCHC@meta.data$integrated_snn_res.0.6==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(integrate.filtered.MCHC@meta.data$integrate.filtered.MCHC==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
integrate.filtered.MCHC@meta.data$customclassif = ""

for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]}

integrate.filtered.MCHC@meta.data$celltypes = ""

for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  integrate.filtered.MCHC@meta.data$customclassif[integrate.filtered.MCHC@meta.data$integrated_snn_res.0.6 == j] = as.character(cl_type$type[1]) }

integrate.filtered.MCHC$celltypes <- integrate.filtered.MCHC$customclassif

#BubbleHeatmap
install.packages('devtools')
devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

#(Fig.S8a)
jjDotPlot(integrate.filtered.MCHC,
          gene = c("Tim3-HAVCR2-AHS0016-pAbO","IgM-IGHM-AHS0198-pAbO","IgD-IGHD-AHS0058-pAbO","HLA-DR-CD74-AHS0035-pAbO","GITR-TNFRSF18-AHS0104-pAbO","CXCR6-CXCR6-AHS0148-pAbO","CXCR5-CXCR5-AHS0039-pAbO","CD8:SK1-CD8A-AHS0228-pAbO","CD4:SK3-CD4-AHS0032-pAbO","CD3:UCHT1-CD3E-AHS0231-pAbO","CD28:L293-CD28-AHS0138-pAbO","CD27:M-T271-CD27-AHS0025-pAbO","CD279:EH12-1-PDCD1-AHS0014-pAbO","CD278-ICOS-AHS0012-pAbO","CD272-BTLA-AHS0052-pAbO","CD25:2A3-IL2RA-AHS0026-pAbO","CD19:SJ25C1-CD19-AHS0030-pAbO","CD196-CCR6-AHS0034-pAbO","CD183-CXCR3-AHS0031-pAbO","CD16:3G8-FCGR3A-AHS0053-pAbO","CD161:HP-3G10-KLRB1-AHS0205-pAbO","CD14:MPHIP9-CD14-AHS0037-pAbO","CD137-TNFRSF9-AHS0003-pAbO","CD134:ACT35-TNFRSF4-AHS0013-pAbO", "CD127-IL7R-AHS0028-pAbO","CD11c:B-LY6-ITGAX-AHS0056-pAbO", "CCR7-CCR7-AHS0273-pAbO", "CD62L:DREG-56-SELL-AHS0049-pAbO", "CD45RA:HI100-PTPRC-AHS0009-pAbO", "CD56:NCAM16.2-NCAM1-AHS0019-pAbO"),
          dot.col = c('darkblue', 'white', 'darkred'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'cellidentity')

#(Fig.S8c)
jjDotPlot(integrate.filtered.MCHC,
          gene = c("CTLA4","TIGIT","ENTPD1","HAVCR2","LAG3","CD28","CD40","CD80","CD86","ICOS","TNFRSF13C","TNFRSF17","TNFRSF25","TNFRSF4","TNFSF10","TNFSF13","TNFSF13B","TNFSF14","CCL20","CCL3","CCL4","CCL5","CXCL16","CXCL2","CXCL3","CXCL5","CCR10","CCR1","CCR2","CCR3","CCR4","CCR5","CCR8","CCR9","CX3CR1","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","GZMA","GZMK","GZMH","GZMB","GNLY","PRF1","NKG7","BAX","BCL2","BIRC3","FAS","FASLG","MKI67","PCNA","IFNG","IL12A","IL13","IL15","IL18","IL2","IL22","IL5","IL6","TGFB1","TGFB3","TNF"),
          dot.col = c('darkblue', 'white', 'darkred'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'cellidentity')

#Cell identity changed based on AbSeq expression and top5 gene heatmap 
sctype_scores[9, "type"] <- "Tregs"
sctype_scores[16, "type"] <- "Other"
sctype_scores[14, "type"] <- "CD11c+ CD14+"
sctype_scores[12, "type"] <- "CD11c+ CD16+"
sctype_scores[7, "type"] <- "CD161+ NK"
sctype_scores[10, "type"] <- "Tc1 CD8"
sctype_scores[13, "type"] <- "Tc17 CD8"
sctype_scores[2, "type"] <- "B cells"
sctype_scores[8, "type"] <- "gammadelta"

#(Fig.6a-b)
DimPlot(integrate.filtered.MCHC, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') + NoLegend()
DimPlot(integrate.filtered.MCHC, reduction = "umap", label = TRUE, repel = TRUE, split.by = 'Sample_Name') + NoLegend()

#Removing AbSeq
counts_trial <- GetAssayData(JoinLayers(integrate.filtered.MCHC, assay = "RNA"))
counts_trial <- counts_trial[-(which(rownames(counts_trial) %in% c("IgD-IGHD-AHS0058-pAbO", "IgM-IGHM-AHS0198-pAbO", "Tim3-HAVCR2-AHS0016-pAbO", "CD45RA:HI100-PTPRC-AHS0009-pAbO", "CD4:SK3-CD4-AHS0032-pAbO", "CD56:NCAM16.2-NCAM1-AHS0019-pAbO", "CD62L:DREG-56-SELL-AHS0049-pAbO", "CD8:SK1-CD8A-AHS0228-pAbO", "CXCR5-CXCR5-AHS0039-pAbO", "CXCR6-CXCR6-AHS0148-pAbO", "GITR-TNFRSF18-AHS0104-pAbO", "HLA-DR-CD74-AHS0035-pAbO", "CD196-CCR6-AHS0034-pAbO", "CD19:SJ25C1-CD19-AHS0030-pAbO", "CD25:2A3-IL2RA-AHS0026-pAbO", "CD272-BTLA-AHS0052-pAbO", "CD278-ICOS-AHS0012-pAbO", "CD279:EH12-1-PDCD1-AHS0014-pAbO", "CD27:M-T271-CD27-AHS0025-pAbO", "CD28:L293-CD28-AHS0138-pAbO", "CD3:UCHT1-CD3E-AHS0231-pAbO", "CCR7-CCR7-AHS0273-pAbO", "CD11c:B-LY6-ITGAX-AHS0056-pAbO", "CD127-IL7R-AHS0028-pAbO", "CD134:ACT35-TNFRSF4-AHS0013-pAbO", "CD137-TNFRSF9-AHS0003-pAbO", "CD14:MPHIP9-CD14-AHS0037-pAbO", "CD161:HP-3G10-KLRB1-AHS0205-pAbO", "CD16:3G8-FCGR3A-AHS0053-pAbO", "CD183-CXCR3-AHS0031-pAbO"))),]
noabseq <- subset(integrate.filtered.MCHC, features = rownames(counts_trial))

if ('IgD-IGHD-AHS0058-pAbO' %in% rownames(noabseq)) {
  print("'IgD-IGHD-AHS0058-pAbO' gene is still present in the Seurat object.")
} else {
  print("'IgD-IGHD-AHS0058-pAbO' gene has been successfully removed from the Seurat object.")
}

#Heatmap of top 5 genes (Fig.S8b)
Idents(noabseq) <- "celltypes"
pbmc.markers <- FindAllMarkers(JoinLayers(noabseq, only.pos = TRUE))
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
candidate_genes <- top5$gene

DoHeatmap(JoinLayers(noabseq), group.by= "celltypes", features = candidate_genes)

#Frequency of cell types per sample (Fig.6c)
table <- table(integrate.filtered.MCHC@meta.data$integrated_snn_res.0.6, integrate.filtered.MCHC@meta.data$Sample_Tag, integrate.filtered.MCHC@meta.data$celltypes)
integrate.filtered.MCHC.df <- as.data.frame(table)
colnames(integrate.filtered.MCHC.df) <- c("cluster", "sample_tag", "id", "freq")

integrate.filtered.MCHC.df$Freq_Percentage <- with(integrate.filtered.MCHC.df, (freq / tapply(freq, sample_tag, sum)[sample_tag]))

ggplot(integrate.filtered.MCHC.df, aes(x = sample_tag, y = Freq_Percentage, fill = factor(id))) +
  geom_bar(stat = "identity") +  # Use stat = "identity" to ensure bars represent the actual values
  labs(x = "Donors", y = "Frequency (%)", title = "Cluster Frequencies per Sample") +
  scale_y_continuous(labels = scales::percent_format()) +  # Format y-axis labels as percentages
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#CHAPTER 6 differential expression
library(tidyverse)
library(Seurat)
library(data.table)
library(ggpubr)

integrate.filtered.MCHC$cohort.celltypes <- paste(integrate.filtered.MCHC$Sample_Name, integrate.filtered.MCHC$celltypes, sep = "_")

#(Supplementary Table 6b)
#Differential expression analysis
DefaultAssay(noabseq) <- "RNA"
Idents(noabseq) <- "cohort.celltypes"
tregs.de <- FindMarkers(JoinLayers(noabseq),
                          ident.1 = "MC_Tregs", 
                          ident.2 = "HC_Tregs", 
                          verbose = FALSE)
