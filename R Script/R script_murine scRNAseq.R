#R script for Fig2, FigS4, FigS5 and Supplementary Table 3
#Brief Experiment Background: 
#C57BL/6 WT and Pd1-/- mice were subcutaneously injected with B16F10 cells. 
#TILs were harvested on d14 and subjected to BD Rhapsody single-cell RNA seq.
###############################################################################

#Section 1: Data Import
library(tidyverse)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library (HGNChelper)

expMat <- ReadMtx (mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv", 
                   cell.column = 1, feature.column = 2, cell.sep = "\t", feature.sep= "\t", 
                   skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE)
expMat[1:5,1:5]

rna <- Seurat::CreateSeuratObject(counts = expMat, min.cells = 0, min.features = 0, assay = "RNA")
rna@meta.data %>% head

smk <- fread(file = "Sampletag.csv", sep = ",", header=TRUE) %>% 
data.frame(row.names = 1)

smk[1:5,]

rna <- AddMetaData(object = rna, metadata = smk)
rna@meta.data %>% head()

rna@meta.data %>% 
group_by(Sample_Name) %>% 
tally() %>% 
mutate(pct = n/sum(n)) %>% 
ggplot(aes(x=Sample_Name, y=n, fill= Sample_Name)) + 
geom_bar(stat = "identity") +
geom_text(aes(y = n,
              label = paste0(n, "\n", scales::percent(pct)),
              vjust= -0.5, size = 2.5), show.legend = FALSE) + 
theme_classic()

rna <- subset(rna, subset = Sample_Name %in% c("Multiplet", "Undetermined"), invert = T)

#Section 2: Quality Control
rna@meta.data %>% 
ggplot(aes(color=Sample_Name, x=nCount_RNA, fill= Sample_Name)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
ggtitle("nCount_RNA")

rna@meta.data %>% 
ggplot(aes(color=Sample_Name, x=nFeature_RNA, fill= Sample_Name)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
ggtitle("nFeature_RNA")

rna@meta.data %>% 
ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
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

filtered.rna <- subset(x = rna, subset = (nCount_RNA >= 10) & (nFeature_RNA >= 20) & (nFeature_RNA <= 300))

# Extract counts
counts <- GetAssayData(object = rna, slot = "counts")

# Output a logical matrix specifying for each gene whether or not 
# there are more than zero counts per cell
nonzero <- counts > 0
 
# Sum all TRUE values and return TRUE 
# if equal or more than 5 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 1
 
# Only keep those genes
filtered_counts <- counts[keep_genes, ]

# Create Seurat object
filtered.rna <- CreateSeuratObject(filtered_counts, meta.data = rna@meta.data)
filtered.rna@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nCount_RNA, fill= Sample_Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  ggtitle("nCount_RNA")

filtered.rna@meta.data %>% 
  ggplot(aes(color=Sample_Name, x=nFeature_RNA, fill= Sample_Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  ggtitle("nFeature_RNA")

# Cell Count before filter
rna@meta.data %>% 
group_by(Sample_Name) %>% 
dplyr::summarise(cell_number = length(Sample_Name))

# Cell Count after filter
filtered.rna@meta.data %>% 
  group_by(Sample_Name) %>% 
  dplyr::summarise(cell_number = length(Sample_Name))

#Normalization
filtered.rna <- NormalizeData(filtered.rna, 
                              normalization.method = "LogNormalize")

# split the dataset into a list of two seurat objects (ko and wt)
split.filtered.rna <- SplitObject(filtered.rna, split.by = "Sample_Name")

# normalize and identify variable features for each dataset independently
split.filtered.rna <- lapply(X = split.filtered.rna, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select the most variable features to use for integration
integ.features <- SelectIntegrationFeatures(object.list = split.filtered.rna, 
                                            nfeatures = 2000) 
anchors <- FindIntegrationAnchors(object.list = split.filtered.rna,
                                  anchor.features = integ.features)

integrate.filtered.rna <- IntegrateData(anchorset = anchors)
integrate.filtered.rna <- ScaleData(integrate.filtered.rna, verbose = FALSE)
integrate.filtered.rna <- RunPCA(integrate.filtered.rna, npcs = 50, verbose = FALSE)

ElbowPlot(integrate.filtered.rna, ndims = 50)

integrate.filtered.rna <- RunUMAP(integrate.filtered.rna, 
                                  reduction = "pca", 
                                  dims = 1:10)

integrate.filtered.rna <- FindNeighbors(integrate.filtered.rna, 
                                        reduction = "pca", 
                                        dims = 1:10)

#Section 3: Clustering
# Determine the K-nearest neighbor graph
integrate.filtered.rna <- FindNeighbors(object = integrate.filtered.rna, dims = 1:10)

# Determine the clusters for various resolutions                                
integrate.filtered.rna <- FindClusters(object = integrate.filtered.rna,
                                       resolution = c(0.2, 0.4, 0.6,  0.8, 1, 1.2, 1.4),
                                       verbose = F)

integrate.filtered.rna@meta.data %>% 
dplyr::select(contains("integrat")) %>% 
map_int(~ unique(.x) %>% length)

library(clustree)
clustree(integrate.filtered.rna, prefix = "integrated_snn_res.")

Idents(object = integrate.filtered.rna) <- "integrated_snn_res.0.6"

#UMAP visualization
DimPlot(integrate.filtered.rna,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

#UMAP of clusters by sample_name
DimPlot(integrate.filtered.rna, 
        label = FALSE, 
        split.by = "Sample_Name") 

#Section 4: Marker Identification
library(tidyverse)
library(Seurat)
library(data.table)
library(ggpubr)
library(HGNChelper)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list <-  gene_sets_prepare("ScTypeDB_short.xlsx") #Available in Supplementary Table 3a
es.max <-  sctype_score(scRNAseqData = integrate.filtered.rna@assays[["integrated"]]@scale.data, 
scaled = TRUE, 
gs = gs_list$gs_positive, 
gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(integrate.filtered.rna@meta.data$integrated_snn_res.0.6), function(cl){
es.max.cl = sort(rowSums(es.max[ ,rownames(integrate.filtered.rna@meta.data[integrate.filtered.rna@meta.data$integrated_snn_res.0.6==cl, ])]), decreasing = !0)
head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(integrate.filtered.rna@meta.data$integrate.filtered.rna==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
integrate.filtered.rna@meta.data$cellidentity = ""

for(j in unique(sctype_scores$cluster)){
cl_type = sctype_scores[sctype_scores$cluster==j,]; 
integrate.filtered.rna@meta.data$cellidentityf[integrate.filtered.rna@meta.data$integrated_snn_res.0.6 == j] = as.character(cl_type$type[1]) }

Idents(integrate.filtered.rna)<- "cellidentity"                                                                                        
DimPlot(integrate.filtered.rna, reduction = "umap",label = F, repel = TRUE) is #Fig2a
DimPlot(integrate.filtered.rna, reduction = "umap",label = F, repel = TRUE, split.by = "Sample_Name") #Fig.2b

#Visualize cell numbers in each clusters in each sample
FetchData(integrate.filtered.rna, 
          vars = c("ident", "Sample_Name")) %>%
  dplyr::count(ident, Sample_Name) %>%
  tidyr::spread(ident, n)

#Umap plot of Cd4 and Foxp3 expression (Fig.2c)
FeaturePlot(integrate.filtered.rna, 
            reduction = "umap", 
            features = c("Cd4", "Foxp3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE)

#Dotplot (Fig.2d, S4c)
library(scRNAtoolVis)
jjDotPlot(integrate.filtered.rna,
         gene = c("Cd4","Cd8a","Cd8b1", 
                  "Eomes", "Foxp3", "Tbx21", "Gata3", "Rora","Rorc", "Bcl6", "Irf4", "Maf",
                  "Ctla4", "Entpd1", "Havcr2", "Lag3", "Pdcd1", "Tigit",
                  "Cd28", "Cd40", "Icos", "Tnfrsf1b", "Tnfrsf4", "Tnfrsf8", "Tnfrsf9", "Tnfrsf13b", "Tnfrsf17", "Tnfrsf18", "Tnfrsf25",
                  "Ccr7", "Lef1", "Sell", "Tcf7",
                  "Mki67", "Pcna", 
                  "Il1a", "Il1b", "Il2", "Il4", "Il6", "Il10", "Il12a", "Il12b", "Il13", "Il15", "Il17a", "Il17f", "Il18", "Il21", "Il33", "Ifng", "Tgfb1", "Tgfb3",
                  "Gzma", "Gzmb", "Gzmk", "Nkg7", "Prf1"),
         dot.col = c('blue', 'white', 'red'),
         tile.geom = F,
         ytree = F,
         rescale.min = -2,
         rescale.max = 2,
         midpoint = 0,
         id = 'cellidentity')

#Since cluster with cell identity "Effector CD8+ T cells" co-expressing both CD4 and CD8 hence rename it as "T Effectors" and proceed to sub-clustering
integrate.filtered.rna@meta.data$cellidentity[integrate.filtered.rna@meta.data$cellidentity == 'Effector CD8+ T cells'] <- 'T Effectors'
TEffectors <- subset(integrate.filtered.rna, idents = "T Effectors") 
#repeat previous steps with TEffectors starting from ScaleData in Section 2 till cell identity allocation in Section 4; resolution 0.8 was used in function "FindCluster".

#Violin plot of WT Tregs vs PD1ko Tregs, with clustered Wilcoxon rank-sum test applied (Fig.S4b)
Idents(integrate.filtered.rna) <- "cellidentity"
Treg <- subset(integrate.filtered.rna, idents = "Treg")
Idents(Treg) <- "cellidentity.group"
Treg <- SetIdent(Treg, 
                 value = factor(Idents(Treg), 
                 levels = c("wt_Treg", "ko_Treg")))
VlnPlot(object = Treg, 
        features = c("Tnfrsf8"), #also visualized Tnfrsf18, Lrrc32, Ctla4, Tigit
        pt.size = 0.5,
        log= F, adjust =0.4, alpha = 1, 
        fill.by = "feature",
        cols = c("#FCDEDD", "#BC6F7E"),
        idents = c("wt_Treg", "ko_Treg"))

library(clusrank)
library(dplyr)

gene <- "Tnfrsf8" #also tested Tnfrsf18, Lrrc32, Ctla4, Tigit
expression_data <- FetchData(object = Treg, vars = c("Sample_Name", "Sample_Tag", gene))
clusWilcox.test(Tnfrsf8 ~ Sample_Name + cluster(Sample_Tag), data = expression_data, method = "ds")

# Average expression of a specific gene in different individual (Fig.2h)
Idents(Treg) <- "Sample_Tag"
avg_expr <- FetchData(Treg, vars = c("Tnfrsf8", "Ctla4", "Lrrc32", "Tigit", "Tnfrsf18"), layer = "data")
avg_expr$Sample_Tag <- Treg$Sample_Tag

avg_expr_summary <- avg_expr %>%
  group_by(Sample_Tag) %>%
  summarise(across(c("Tnfrsf8", "Ctla4", "Lrrc32", "Tigit", "Tnfrsf18"), mean), .groups = "drop")

#Section 5: Differential Expression Analysis
# libraries for this section
library(tidyverse)
library(Seurat)
library(data.table)
library(ggpubr)

#Removing AbSeq
all_genes <- rownames(integrate.filtered.rna)
genes_to_remove <- grep("pAb", all_genes, value = TRUE)
counts_pAb <- GetAssayData(JoinLayers(integrate.filtered.rna, assay = "RNA"))
counts_pAb <- counts_pAb[-(which(rownames(counts_pAb) %in% genes_to_remove)),]
Idents(integrate.filtered.rna) <- "features"
noabseq <- subset(JoinLayers(integrate.filtered.rna, 
                             features = setdiff(rownames(integrate.filtered.rna),genes_to_remove)))
noabseq$cellidentity.group <- paste(noabseq$Sample_Name, noabseq$cellidentity, sep = "_")

#Differential gene expression analysis
Idents(noabseq) <- "cellidentity.group"
Treg.DE <- FindMarkers(noabseq, 
                       ident.1 = "ko_Treg",
                       ident.2 = "wt_Treg",
                       verbose = FALSE) #also performed on CD8Teff and CD4Teff

library(openxlsx)
write.xlsx(Treg.DE, "TregDE.xlsx", rowNames = T) #available in Supplementary Table 3b-d for Treg, CD8Teff &CD4Teff

#Volcano Plot (Fig.2f, S4c-d)
DEresult <- read.xlsx('TregDE.xlsx') #also visualized CD8 Teff or CD4 Teff; available in Supplementary Table 3b-d
DEresult$DElabel <- NA
DEresult <- DEresult %>%
  mutate(gene_type = case_when(avg_log2FC >= 0.585 & p_val_adj <= 0.05 ~ "up",
                               avg_log2FC <= -0.585 & p_val_adj <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
DEresult %>%
  count(gene_type)

cols <- c("up" = "red", "down" = "blue", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 1)

sig_il_genes <- DEresult %>%
  filter(Gene %in% c("Tnfrsf8","Tnfrsf18", "Ctla4", "Tigit", "Lrrc32","Pdcd1"))

up_il_genes <- DEresult %>%
  filter(Gene %in% c("Tnfrsf8"))

ns_il_genes <- DEresult %>%
  filter(Gene %in% c("Tnfrsf18", "Ctla4", "Tigit", "Lrrc32"))

down_il_genes <- DEresult %>%
  filter(Gene == "Pdcd1")

DEresult <- DEresult%>%
  mutate(gene_type = fct_relevel(gene_type, "up", "down")) 

ggplot(data = DEresult,
       aes(x = avg_log2FC,
           y = -log10(p_val_adj))) + 
  geom_point(aes(colour = gene_type), 
             alpha = 0.3,
             shape = 19,
             size = 2) + 
  geom_point(data = up_il_genes,
             shape = 21,
             size = 3, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = down_il_genes,
             shape = 21,
             size = 3, 
             fill = "darkblue", 
             colour = "black") +
  geom_point(data = ns_il_genes,
             shape = 21,
             size = 3, 
             fill = "darkgrey", 
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.585, 0.585),
             linetype = "dashed") +
  geom_label_repel(data = sig_il_genes, 
                   aes(label = Gene),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  scale_x_continuous(breaks = c(seq(-6, 4, 2)),     
                     limits = c(-6, 4))  +
  labs(title = "Gene expression changes in Pd1-/- Treg versus WT Treg",
       x = "log2(FC)",
       y = "-log10(adjusted P-value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

#Section 6: Enrichment Analysis of Tregs
# libraries for this section
library(tidyverse)
library(Seurat)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(biomaRt)
DEresults <- read.xlsx('TregDE.xlsx')
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
results <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
filters = "external_gene_name", 
values = DEresults$Gene, #available in Supplementary Table 3b
mart = mart)

results %>% head

#keep genes with matched Entrez ID
rownames(DEresult) <- DEresult$Gene
DEresult$Gene <- NULL
DEresult <- DEresult %>% 
  rownames_to_column("Gene") %>% 
  left_join(., results, by = c("Gene" = "external_gene_name")) %>% 
  filter(!is.na("entrezgene_id")) %>% 
  filter(!is.na(p_val_adj))

#GO analysis available in Supplementary Table 3e
#Isolate enriched genes in Pd1-/- Treg
DEresult_KO <- dplyr::filter(DEresult, p_val_adj < 0.05 & avg_log2FC > 0.585)
ego <- enrichGO(gene = DEresult_KO$entrezgene_id,
keyType = "ENTREZID",
OrgDb = org.Mm.eg.db,
ont = "BP",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05,
readable = TRUE)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
options(enrichplot.colours = c("red","blue"))
dotplot(ego,showCategory=8)

#GSEA analysis
DEresult <- DEresult[order(-DEresult$avg_log2FC), ]
gene_list <- DEresult$avg_log2FC
names(gene_list) <- DEresult$entrezgene_id

# perform GSEA
gse <- gseGO(geneList=gene_list,
ont = "BP",
keyType = "ENTREZID",
OrgDb = "org.Mm.eg.db",
pAdjustMethod = "BH",
pvalueCutoff = 1)

# convert the results into a dataframe
gse.table <- as.data.frame(gse) %>% 
mutate(geneSetID = 1:dim(.)[1])

#GSEA plot of "Regulation of Cytokine Production"
gseaplot2(gse, geneSetID = c(12)) #Result of GSEA of Treg available in Supplementary Table 3f, Fig.2g

#Section 7: CellChat Analysis
library(CellChat)
library(reticulate)

seurat_objectWT <- subset(integrate.filtered.rna, idents = c("wt"))
seurat_objectKO <- subset(integrate.filtered.rna, idents = c("ko"))

data.inputWT <-GetAssayData(JoinLayers(seurat_objectWT, assay ="RNA", slot = "data"))
data.inputKO <-GetAssayData(JoinLayers(seurat_objectKO, assay ="RNA", slot = "data"))

metaWT <- data.frame(group = seurat_objectWT$cellidentity, row.names = names(labels))
metaKO <- data.frame(group = seurat_objectKO$cellidentity, row.names = names(labels))

cellchatWT <- createCellChat(object = data.inputWT, meta =  metaWT, group.by = "group")
cellchatKO <- createCellChat(object = data.inputKO, meta =  metaKO, group.by = "group")

CellChatDB.mouse <- CellChatDB.mouse
cellchatWT@DB <- CellChatDB.mouse
cellchatKO@DB <- CellChatDB.mouse

#Pre-processing the expression data
#subset the expression data to use less RAM
CellchatKO.use<-subsetData(cellchatKO)
CellchatWT.use<-subsetData(cellchatWT)

#Pre-processing the expression data
CellchatWT.use<-identifyOverExpressedGenes(CellchatWT.use)
CellchatWT.use<-identifyOverExpressedInteractions(CellchatWT.use)

CellchatKO.use<-identifyOverExpressedGenes(CellchatKO.use)
CellchatKO.use<-identifyOverExpressedInteractions(CellchatKO.use)

#Project gene expression data onto protein-protein interaction (PPI)
CellchatKO.use<- projectData(CellchatKO.use, PPI.mouse)
CellchatWT.use<- projectData(CellchatWT.use, PPI.mouse)

#Compute the communication probability and infer cellular communication network
CellchatKO.use<- computeCommunProb(CellchatKO.use, raw.use = FALSE)
CellchatWT.use<- computeCommunProb(CellchatWT.use, raw.use = FALSE)

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
CellchatWT.use<-filterCommunication(CellchatWT.use, min.cells = 3)
CellchatKO.use<-filterCommunication(CellchatKO.use, min.cells = 3)

#Infer the cell-cell communication at a signaling pathway level
CellchatKO.use<-computeCommunProbPathway(CellchatKO.use)
CellchatWT.use<-computeCommunProbPathway(CellchatWT.use)

#Calculate the aggregated cell-cell communication network
CellchatKO.use<-aggregateNet(CellchatKO.use)
cellchat@net$count
cellchat@net$weight
CellchatWT.use<-aggregateNet(CellchatWT.use)
cellchat@net$count
cellchat@net$weight

#examine the signalling sent from each cell group (Fig.S5c,d)
matKO.use <- CellchatKO.use@net$weight
par(mfrow= c(2,5), xpd = TRUE)
for (i in 1:nrow(matKO.use)) {
  mat2 <-matrix(0, nrow = nrow(matKO.use), ncol = ncol(matKO.use), dimnames = dimnames (matKO.use))
  mat2 [i, ] <- matKO.use[i, ]
  netVisual_circle(mat2, vertex.weight = groupSizeKO.use, weight.scale = T,
                   edge.weight.max = max(matKO.use), title.name = rownames(matKO.use)[i])}

matWT.use <- CellchatWT.use@net$weight
par(mfrow= c(2,5), xpd = TRUE)
for (i in 1:nrow(matWT.use)) {
  mat2 <-matrix(0, nrow = nrow(matWT.use), ncol = ncol(matWT.use), dimnames = dimnames (matWT.use))
  mat2 [i, ] <- matWT.use[i, ]
  netVisual_circle(mat2, vertex.weight = groupSizeWT.use, weight.scale = T,
                   edge.weight.max = max(matWT.use), title.name = rownames(matWT.use)[i])}

CellchatKO.use@netP[["pathways"]]
CellchatWT.use@netP[["pathways"]]
extractEnrichedLR(CellchatKO.use, signaling = c(CellchatKO.use@netP[["pathways"]]),
                  geneLR.return = TRUE)

#Chord diagram: define source and target cell types (Fig.2m)
netVisual_chord_gene(CellchatKO.use, sources.use = c(9), targets.use = c(1,2,3,4,5,6,7,8,9),
                     lab.cex = 0.5, legend.pos.x = 15)
netVisual_chord_gene(CellchatWT.use, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(9),
                     lab.cex = 0.5, legend.pos.x = 15)

#Systematic analysis of cell-cell communication networks
library(NMF)
library(ggalluvial)
library(CellChat)

#1.compute the network centrality scores
CellchatWT.use <- netAnalysis_computeCentrality(CellchatWT.use, slot.name = "netP")
CellchatKO.use <- netAnalysis_computeCentrality(CellchatKO.use, slot.name = "netP")

#Heatmap to visualise dominant cell types for each signaling pathway (Fig.S5a-b)
netAnalysis_signalingRole_heatmap(CellchatWT.use, pattern = "outgoing", height = 11) 
netAnalysis_signalingRole_heatmap(CellchatWT.use, pattern = "incoming", height = 11) 
netAnalysis_signalingRole_heatmap(CellchatKO.use, pattern = "outgoing", height = 11) 
netAnalysis_signalingRole_heatmap(CellchatKO.use, pattern = "incoming", height = 11) 

#2.Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate
#Identify and visualize incoming communication pattern of secreting cells (Fig.2k-l)
selectK(CellchatWT.use, pattern = "incoming")
ellchatWT.use <- identifyCommunicationPatterns(ellchatWT.use, pattern = "outgoing",
                                          k = 4, width = 5, height = 9)

selectK(CellchatKO.use, pattern = "incoming")
CellchatKO.use <- identifyCommunicationPatterns(CellchatKO.use, pattern = "outgoing",
                                          k = 4, width = 5, height = 9)
