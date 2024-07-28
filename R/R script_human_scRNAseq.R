#R script for Fig6, FigS6 and Supplementary Table 3
#Brief Experiment Background: 
#n=3 healthy controls and n=3 stage IV melanoma patients 
#PBMCs were subjected to BD Rhapsody single-cell RNA seq.
###############################################################################

#Section 1 data import
library(tidyverse)
library(Seurat)
library(data.table)

expMat_C1 <- ReadMtx(mtx = "C1 matrix.mtx", cells = "C1 barcodes.tsv", features = "C1 features.tsv",
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

#packages 
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

DimPlot(integrate.filtered.MCHC, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') + NoLegend()
DimPlot(integrate.filtered.MCHC, reduction = "umap", label = TRUE, repel = TRUE, split.by = 'Sample_Name') + NoLegend()

#BubbleHeatmap
install.packages('devtools')
devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

#(Fig.6d)
jjDotPlot(integrate.filtered.MCHC,
          gene = c("CTLA4","TIGIT","ENTPD1","HAVCR2","LAG3","CD28","CD40","CD80","CD86","ICOS","TNFRSF13C","TNFRSF17","TNFRSF25","TNFRSF4","TNFSF10","TNFSF13","TNFSF13B","TNFSF14","CCL20","CCL3","CCL4","CCL5","CXCL16","CXCL2","CXCL3","CXCL5","CCR10","CCR1","CCR2","CCR3","CCR4","CCR5","CCR8","CCR9","CX3CR1","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","GZMA","GZMK","GZMH","GZMB","GNLY","PRF1","NKG7","BAX","BCL2","BIRC3","FAS","FASLG","MKI67","PCNA","IFNG","IL12A","IL13","IL15","IL18","IL2","IL22","IL5","IL6","TGFB1","TGFB3","TNF"),
          dot.col = c('darkblue', 'white', 'darkred'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'cellidentity')

#(Fig.S6a)
jjDotPlot(integrate.filtered.MCHC,
          gene = c("Tim3-HAVCR2-AHS0016-pAbO","IgM-IGHM-AHS0198-pAbO","IgD-IGHD-AHS0058-pAbO","HLA-DR-CD74-AHS0035-pAbO","GITR-TNFRSF18-AHS0104-pAbO","CXCR6-CXCR6-AHS0148-pAbO","CXCR5-CXCR5-AHS0039-pAbO","CD8:SK1-CD8A-AHS0228-pAbO","CD4:SK3-CD4-AHS0032-pAbO","CD3:UCHT1-CD3E-AHS0231-pAbO","CD28:L293-CD28-AHS0138-pAbO","CD27:M-T271-CD27-AHS0025-pAbO","CD279:EH12-1-PDCD1-AHS0014-pAbO","CD278-ICOS-AHS0012-pAbO","CD272-BTLA-AHS0052-pAbO","CD25:2A3-IL2RA-AHS0026-pAbO","CD19:SJ25C1-CD19-AHS0030-pAbO","CD196-CCR6-AHS0034-pAbO","CD183-CXCR3-AHS0031-pAbO","CD16:3G8-FCGR3A-AHS0053-pAbO","CD161:HP-3G10-KLRB1-AHS0205-pAbO","CD14:MPHIP9-CD14-AHS0037-pAbO","CD137-TNFRSF9-AHS0003-pAbO","CD134:ACT35-TNFRSF4-AHS0013-pAbO", "CD127-IL7R-AHS0028-pAbO","CD11c:B-LY6-ITGAX-AHS0056-pAbO", "CCR7-CCR7-AHS0273-pAbO", "CD62L:DREG-56-SELL-AHS0049-pAbO", "CD45RA:HI100-PTPRC-AHS0009-pAbO", "CD56:NCAM16.2-NCAM1-AHS0019-pAbO"),
          dot.col = c('darkblue', 'white', 'darkred'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'cellidentity')

#Cell identity changed based on specific gene and AbSeq expression from heatmap 
#example 
sctype_scores[9, "type"] <- "Tregs"
sctype_scores[16, "type"] <- "Other"
sctype_scores[14, "type"] <- "CD11c+ CD14+"
sctype_scores[12, "type"] <- "CD11c+ CD16+"
sctype_scores[7, "type"] <- "CD161+ NK"
sctype_scores[10, "type"] <- "Tc1 CD8"

#Removing AbSeq
counts_trial <- GetAssayData(JoinLayers(integrate.filtered.MCHC, assay = "RNA"))
counts_trial <- counts_trial[-(which(rownames(counts_trial) %in% c("IgD-IGHD-AHS0058-pAbO", "IgM-IGHM-AHS0198-pAbO", "Tim3-HAVCR2-AHS0016-pAbO", "CD45RA:HI100-PTPRC-AHS0009-pAbO", "CD4:SK3-CD4-AHS0032-pAbO", "CD56:NCAM16.2-NCAM1-AHS0019-pAbO", "CD62L:DREG-56-SELL-AHS0049-pAbO", "CD8:SK1-CD8A-AHS0228-pAbO", "CXCR5-CXCR5-AHS0039-pAbO", "CXCR6-CXCR6-AHS0148-pAbO", "GITR-TNFRSF18-AHS0104-pAbO", "HLA-DR-CD74-AHS0035-pAbO", "CD196-CCR6-AHS0034-pAbO", "CD19:SJ25C1-CD19-AHS0030-pAbO", "CD25:2A3-IL2RA-AHS0026-pAbO", "CD272-BTLA-AHS0052-pAbO", "CD278-ICOS-AHS0012-pAbO", "CD279:EH12-1-PDCD1-AHS0014-pAbO", "CD27:M-T271-CD27-AHS0025-pAbO", "CD28:L293-CD28-AHS0138-pAbO", "CD3:UCHT1-CD3E-AHS0231-pAbO", "CCR7-CCR7-AHS0273-pAbO", "CD11c:B-LY6-ITGAX-AHS0056-pAbO", "CD127-IL7R-AHS0028-pAbO", "CD134:ACT35-TNFRSF4-AHS0013-pAbO", "CD137-TNFRSF9-AHS0003-pAbO", "CD14:MPHIP9-CD14-AHS0037-pAbO", "CD161:HP-3G10-KLRB1-AHS0205-pAbO", "CD16:3G8-FCGR3A-AHS0053-pAbO", "CD183-CXCR3-AHS0031-pAbO"))),]
noabseq <- subset(integrate.filtered.MCHC, features = rownames(counts_trial))

if ('IgD-IGHD-AHS0058-pAbO' %in% rownames(noabseq)) {
  print("'IgD-IGHD-AHS0058-pAbO' gene is still present in the Seurat object.")
} else {
  print("'IgD-IGHD-AHS0058-pAbO' gene has been successfully removed from the Seurat object.")
}

#Heatmap of top 5 genes (Fig.S6b)
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

#Differential expression analysis
DefaultAssay(noabseq) <- "RNA"
Idents(noabseq) <- "cohort.celltypes"
tregs.de <- FindMarkers(JoinLayers(noabseq),
                          ident.1 = "MC_Tregs", 
                          ident.2 = "HC_Tregs", 
                          verbose = FALSE)

#Volcano pf HC Tregs compared to Melanoma Tregs(Fig.6e, S6c, S6d)
tregs.de$diffexpressed <- "No Changed"
tregs.de$diffexpressed[tregs.de$avg_log2FC > 1 & tregs.de$p_val < 0.05] <- "Upregulated"
tregs.de$diffexpressed[tregs.de$avg_log2FC < -1 & tregs.de$p_val < 0.05] <- "Downregulated"

tregs.de$DElabel <- NA
tregs.de$gene <- row.names(tregs.de)
tregs.de$DElabel[tregs.de$diffexpressed != "NA"] <- tregs.de$gene[tregs.de$diffexpressed != "NA"]

library(ggrepel)
ggplot(data=tregs.de, aes(x=avg_log2FC, y=-log10(p_val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color ="#999999") +
  geom_point(aes(size = -log10(p_val), color = -log10(p_val))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.99),
        legend.justification = c(0,1)) +
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("darkblue", "blue", "red", "red", "darkred")) +
  scale_size_continuous(range = c(1,7)) +
  guides(col = guide_colorbar(title = "-log10(p_val)"),
         size = "none") +
  geom_text(aes(label = DElabel, color = -log10(p_val),
                size = 3, hjust = 1, vjust = 1.7)) +
  xlab("avg_log2FC") +
  ylab("-log10(p_val)")

#Enrichment Analysis (Fig.6f)
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("gseaplot2")
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(biomaRt)
library(gseaplot2)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                 filters = "external_gene_name", 
                 values = rownames(tregs.de),
                 mart = mart)
results %>% head

tregs.de$gene <- rownames(tregs.de) 

tregs.de<- left_join(tregs.de, results, by = c("gene" = "external_gene_name")) %>% 
  filter(!is.na(entrezgene_id)) %>% 
  filter(!is.na(p_val))

tregs.de.sig <- dplyr::filter(tregs.de, p_val < 0.05, abs(avg_log2FC) > 1)

ego <- enrichGO(gene = tregs.de.sig$entrezgene_id,
                    keyType = "ENTREZID",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
dotplot(ego, showCategory=10) #fig. 6f

#GSEA (Fig.6g)
tregs.de <- tregs.de[order(-tregs.de$avg_log2FC), ]
gene_list <- tregs.de$avg_log2FC
names(gene_list) <- tregs.de$entrezgene_id

gse <- gseGO(geneList=gene_list,
             ont = "BP",
             keyType = "ENTREZID",
             minGSSize = 3,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = "org.Hs.eg.db",
             pAdjustMethod = "none")

gse.table <- as.data.frame(gse) %>% 
  mutate(geneSetID = 1:dim(.)[1])

gseaplot2(gse, geneSetID = c(1,5)) #Result of GSEA of Treg available in Supplementary Table 1d, Fig.2k

# CellChat (Fig. 7)
library(devtools)
devtools::install_github("sqjin/CellChat")
BiocManager::install("ComplexHeatmap")
BiocManager::install("BiocNeighbors")

library(CellChat)
library(Seurat)
library(reticulate)
library(tidyverse)
library(Seurat)
library(BiocManager)
library(Matrix)
library(RCurl)
library(readxl)
library(scales)
library(ggpubr)
library(data.table)
library(AnnotationHub)
library(ensembldb)
library(HGNChelper)
library(clusterProfiler)
library(AnnotationDbi)
library(biomaRt)
library(enrichplot)
library(ggplot2)

#Cell Chat object HC
Idents(integrate.filtered.MCHC) <- "Sample_Name"
seurat_objectHC <- subset(integrate.filtered.MCHC, idents = c("HC"))
data.inputHC <-GetAssayData(JoinLayers(seurat_objectHC, assay ="RNA", slot = "data"))
metaHC <- data.frame(group = seurat_objectHC$customclassif, row.names = names(labels))
cellchatHC <- createCellChat(object = data.inputHC, meta =  metaHC, group.by = "group")
#Cell Chat object MC
seurat_objectMC <- subset(integrate.filtered.MCHC, idents = c("MC"))
data.inputMC <-GetAssayData(JoinLayers(seurat_objectMC, assay ="RNA", slot = "data"))
metaMC <- data.frame(group = seurat_objectMC$customclassif, row.names = names(labels))
cellchatMC <- createCellChat(object = data.inputMC, meta =  metaMC, group.by = "group")

#Cellchat databases
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB.human)
#show structure of databases
dplyr::glimpse(CellChatDB$interaction)
#set all cellchat db for cell-cell communications 
CellChatDB.use <- CellChatDB.human

#Inference of cell-cell communication network analysis
library(CellChat)
library(patchwork)
library(circlize)
options(stringsAsFactors = FALSE)

cellchatHC@DB <- CellChatDB.use
cellchatMC@DB <- CellChatDB.use

cellchatHC<-subsetData(cellchatHC)
cellchatMC<-subsetData(cellchatMC)

cellchatHC <- identifyOverExpressedGenes(cellchatHC)
cellchatHC<-identifyOverExpressedInteractions(cellchatHC)
cellchatMC <- identifyOverExpressedGenes(cellchatMC)
cellchatMC<-identifyOverExpressedInteractions(cellchatMC)

#Project gene expression data onto protein-protein interaction (PPI)
cellchatHC<- projectData(cellchatHC, PPI.human)
cellchatMC<- projectData(cellchatMC, PPI.human)

#Compute the communication probability and infer cellular communication network
cellchatHC <- computeCommunProb(cellchatHC, raw.use = FALSE)
cellchatMC <- computeCommunProb(cellchatMC, raw.use = FALSE)

#Filter out the cell-cell communication
cellchatHC<-filterCommunication(cellchatHC, min.cells = 3)
cellchatMC<-filterCommunication(cellchatMC, min.cells = 3)

#Infer the cell-cell communication at a signaling pathway level
cellchatHC<-computeCommunProbPathway(cellchatHC)
cellchatMC<-computeCommunProbPathway(cellchatMC)

#Calculate the aggregated cell-cell communication network
cellchatHC<-aggregateNet(cellchatHC)
cellchatHC@net$count
cellchatHC@net$weight

cellchatMC<-aggregateNet(cellchatMC)
cellchatMC@net$count
cellchatMC@net$weight

#Visualize the aggregated cell-cell communication network
groupSize <-as.numeric(table(cellchatHC@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatHC@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchatHC@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")

#Examine the signalling sent from each cell group, circle plot (Fig. 7c)
mat <- cellchatHC@net$weight
par(mfrow= c(2,4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames (mat))
  mat2 [i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}

#Visualisation of cell-cell communication network
cellchatHC@netP[["pathways"]]
cellchatMC@netP[["pathways"]]

extractEnrichedLR(cellchatHC, signaling = c(cellchatHC@netP[["pathways"]]),
                  geneLR.return = TRUE)
extractEnrichedLR(cellchatMC, signaling = c(cellchatMC@netP[["pathways"]]),
                  geneLR.return = TRUE)

#Visualize the contribution of each LR pairs to the communication network
netAnalysis_contribution(cellchatHC,
                         signaling = c(cellchatHC@netP[["pathways"]]),
                         title = "Contribution of each LR pairs")
netAnalysis_contribution(cellchatMC,
                         signaling = c(cellchatMC@netP[["pathways"]]),
                         title = "Contribution of each LR pairs")

#Chord diagram (Fig 7d-e)
netVisual_chord_gene(cellchatMC, sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), targets.use = c(16),
                     lab.cex = 0.5, legend.pos.x = 15)
netVisual_chord_gene(cellchatMC, sources.use = c(16), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                     lab.cex = 0.5, legend.pos.x = 15)

library(RColorBrewer)
pathwaysHC <- cellchatHC@netP[["pathways"]]
pathwaysMC <- cellchatMC@netP[["pathways"]]

#Bubble plot of LR pairs (Fig. 7b)
netVisual_bubble(cellchatHC, sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), targets.use = c(16),
                 signaling = c(pathwaysHC),
                 remove.isolate = F)
netVisual_bubble(cellchatHC, sources.use = c(16), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
                 signaling = c(pathwaysHC),
                 remove.isolate = F)
