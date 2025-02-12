#R script for Fig6, S9 and Supplemental table 6
#This script contains the analysis flow for two data sets as follows:
#1. Treatment naive melanoma patients from GSE72056
#2. Anti-PD1 baseline/treated melanoma patients from GSE120575

###############################################################################
#1. Treatment naive melanoma patients from GSE72056
#Libraries
library(Seurat)
library(Matrix)
library(readr)
library(dplyr)
library(ggpubr)
library(openxlsx)
library(scRNAtoolVis)

#1.Data import
data <- read.delim("processed files.txt", header = TRUE, sep = "\t")

#2. Create MetaData from raw data
df_transposed <- t(data)
colnames(df_transposed) <- as.character(df_transposed[1, ])
df_transposed <- df_transposed[-1, ]
df_transposed <- as.data.frame(df_transposed)

meta_data <- data.frame(
  CellID = rownames(df_transposed),
  CellType = df_transposed$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`,
  Individual = df_transposed$tumor
)

#3. Create Expression Data from raw data
expressiondata <- data
expressiondata <- expressiondata[-1,] #repeat another 2 times
rownames(expressiondata) <- expressiondata$Cell
expressiondata$Cell <- NULL
expressiondata <- 2^expressiondata - 1

#4. Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = expressiondata, meta.data = meta_data)

#5. Change Individual and CellType to integer
seurat_obj$CellType <- as.integer(seurat_obj$CellType)
seurat_obj$Individual <- as.integer(seurat_obj$Individual)
Idents(seurat_obj) <- seurat_obj$Individual

#6. Interest individuals with Tregs
TN_so <- subset(seurat_obj, idents = c("79", "80", "82", "53"))

#6.Clustering
TN_so <- FindVariableFeatures(TN_so)
TN_so <- RunPCA(TN_so, npcs = 50, verbose = FALSE)
ElbowPlot(TN_so, ndims = 50)

TN_so <- RunUMAP(TN_so, 
                reduction = "pca", 
                dims = 1:10)

# Determine the K-nearest neighbor graph
TN_so <- FindNeighbors(TN, 
                      reduction = "pca", 
                      dims = 1:10)

# Determine the clusters for various resolutions                                
TN_so <- FindClusters(object = TN_so,
                      resolution = c(0.2, 0.4, 0.6,  0.8, 1, 1.2, 1.4),
                      verbose = F)

TN_so@meta.data %>% 
  dplyr::select(contains("snn")) %>% 
  map_int(~ unique(.x) %>% length)

# Determine the DGE between each clusters with resolution 0.2 (Supplemental Table 6c)
Idents(TN_so) <- "RNA_snn_res.0.2"
DE_Cluster <- FindAllMarkers(TN_so, only.pos = T) 

#Cell annotation based on the result of DGE analysis (p<0.0001)
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '0'] <- 'Melanoma'
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '1'] <- 'T'
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '2'] <- 'B'
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '3'] <- 'Macrophage'
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '4'] <- 'Endo'
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '5'] <- 'CAFs'
TN_so@meta.data$celltypes[TN_so@meta.data$RNA_snn_res.0.2 == '6'] <- 'NK'

#Subdivision of 'T' into CD8+, CD3+CD8+CD4+, CD4+ T and Tregs (CD4+FoxP3+)
Idents(TN_so) <- "celltypes"
T_seurat <- subset(TN_so, idents = "T")
TN_so@meta.data$celltypes[TN_so@meta.data$celltypes == 'T'] <- 'CD4-CD8-T'

grep("^CD8A$", rownames(T_seurat@assays[["RNA"]]$counts))
length(which(T_seurat@assays[["RNA"]]$counts[rownames(T_seurat@assays[["RNA"]])[18230], ] != 0)) 
T_seurat@meta.data$celltypes[which(T_seurat@assays[["RNA"]]$counts[rownames(T_seurat@assays[["RNA"]])[18230], ] != 0)] <- 'CD8 T'

grep("^CD8B$", rownames(T_seurat@assays[["RNA"]]$counts))
length(which(T_seurat@assays[["RNA"]]$counts[rownames(T_seurat@assays[["RNA"]])[13022], ] != 0)) 
T_seurat@meta.data$celltypes[which(T_seurat@assays[["RNA"]]$counts[rownames(T_seurat@assays[["RNA"]])[13022], ] != 0)] <- 'CD8 T'

Idents(Tirosh_so) <- "celltypes"
T_seurat <- subset(Tirosh_so, idents = "T")
grep("^CD4$", rownames(T_seurat@assays[["RNA"]]$counts))
length(which(T_seurat@assays[["RNA"]]$counts[rownames(T_seurat@assays[["RNA"]])[7414], ] != 0)) 
T_seurat@meta.data$celltypes[which(T_seurat@assays[["RNA"]]$counts[rownames(T_seurat@assays[["RNA"]])[7414], ] != 0)] <- 'CD4 T'

Idents(T_seurat) <- "celltypes"
TN_so$celltypes[Cells(T_seurat)] <- paste(Idents(T_seurat))

Idents(TN_so) <- 'celltypes'
CD4_seurat <- subset(TN_so, idents = "CD4 T")
grep("^FOXP3$", rownames(CD4_seurat@assays[["RNA"]]$counts))
length(which(CD4_seurat@assays[["RNA"]]$counts[rownames(CD4_seurat@assays[["RNA"]])[11928], ] != 0)) 
CD4_seurat@meta.data$celltypes[which(CD4_seurat@assays[["RNA"]]$counts[rownames(CD4_seurat@assays[["RNA"]])[11928], ] != 0)] <- 'Treg'
Idents(CD4_seurat) <- "celltypes"
TN_so$celltypes[Cells(CD4_seurat)] <- paste(Idents(CD4_seurat))

DimPlot(TN_so,
        reduction = "umap",
        label = F,
        pt.size = 0.5,
        label.size = 6,
        group.by = "celltypes") #Fig6d

#Frequency (Fig.6e)
table_celltypes <- table(TN_so@meta.data$Individual, TN_so@meta.data$celltypes)
Frequency_celltypes <- as.data.frame(table_celltypes)
Frequency_celltypes$Percentage <- round(Frequency_celltypes$Freq / sum(Frequency_celltypes$Freq) * 100, 1)
colnames(Frequency_celltypes) <- c("Individual", "CellType", "Frequency")

ggplot(Frequency_celltypes, aes(fill=CellType, y=Frequency, x=Individual)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.line = element_line(size = 0.5),
  ) + 
  theme_minimal()

#Marker Identification (Fig.S9a)
jjDotPlot(TN_so,
          gene = c("CD3D", "CD4", "CD8A", "CD19",
                   "EOMES", "FOXP3", "TBX21", "GATA3", "RORA", "RORC", "BCL6", "IRF4", "MAF",
                   "CTLA4", "ENTPD1", "HAVCR2", "LAG3", "PDCD1", "TIGIT",
                   "CD28", "CD40", "ICOS", "TNFRSF1B", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFRSF13B", "TNFRSF17", "TNFRSF18", "TNFRSF25",
                   "CCR7", "LEF1", "SELL", "TCF7",
                   "MKI67", "PCNA", 
                   "IL1A", "IL1B", "IL2", "IL4", "IL6", "IL10", "IL12A", "IL12B", "IL13", "IL15", "IL17A", "IL17F", "IL18", "IL21", "IL33", "IFNG", "TGFB1", "TGFB3",
                   "GZMA", "GZMB", "GZMK", "NKG7", "PRF1",
                   "CD163", "CD14",
                   "PECAM1", "VWF",
                   "FAP",
                   "MIA", "TYR"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'celltypes')

#Analysis of scRNAseq data derived from: 
#our in-house BD Rhapsody scRNAseq of PBMCs from healthy donor (HC) and treatment naive melanoma patients (Mel), 
#and tumor from melanoma patients (publicly mined data), related to Fig6f-i

TN_so$Sample_Name<- "TIL"
TN_so$Sample_Tag <- TN_so$Individual
common_genes <- intersect(rownames(TN_so), rownames(integrate.filtered.MCHC)) #integrate.filtered.MCHC contains our in-house human scRNAseq data
obj1 <- subset(TN_so, features = common_genes)
obj2 <- subset(integrate.filtered.MCHC, features = common_genes)

merged_obj <- merge(obj1, y = obj2, add.cell.ids = c("TIL", "PBMC"))
merged_obj$cellidentity.group <- paste(merged_obj$Sample_Name, merged_obj$celltypes, sep = "_")
Idents(merged_obj) <- "cellidentity.group"
merged_Treg <- subset(merged_obj, idents= c("TIL_Treg", "HC_Tregs", "Mel_Tregs"))

Idents(merged_Treg) <- "cellidentity.group"
VlnPlot(object = merged_Treg, 
        features = c("TNFRSF8"), pt.size = 0.8, #also visualized CTLA4, LRRC32, TIGIT
        log= F, adjust = 30, alpha = 1, fill.by = "feature",y.max = 10,
        group.by = "cellidentity.group")  

gene <- "Tnfsf8" #also tested CTLA4, LRRC32, TIGIT
expression_data <- FetchData(object = merged_Treg, vars = c("Sample_Name", "Sample_Tag", gene))
with(expression_data,pairwise.clusWilcox.test(
  expression_data$Tnfsf8,
  group = Sample_Name,
  cluster = Sample_Tag,
  method = "ds",
  p.adjust.method = "fdr"))

# Average expression of a specific gene in different individual within HC vs Mel vs TIL (Fig.6f-i)
avg_expr <- FetchData(merged_Treg, vars = c("TNFRSF8", "CTLA4", "LRRC32", "TIGIT"), layer = "data")
avg_expr$Sample_Tag <- merged_Treg$Sample_Tag

avg_expr_summary <- avg_expr %>%
  group_by(Sample_Tag) %>%
  summarise(across(c("TNFRSF8", "CTLA4", "LRRC32", "TIGIT"), mean), .groups = "drop")

################################################################################
#2. Anti-PD1 baseline/treated melanoma patients from GSE120575
#Libraries
library(Seurat)
library(Matrix)
library(readr)
library(dplyr)
library(ggpubr)
library(openxlsx)
library(scRNAtoolVis)

#1.Data import and create metadata
expression_data <- read.delim("processed files.txt", header = TRUE, sep = "\t")
colnames(expression_data) <- c(colnames(expression_data)[-1], NA)
expression_data <- expression_data[, -ncol(expression_data)]

metadata <- expression_data[1,]
metadata <- t(metadata)
metadata <- as.data.frame(metadata)

#Include clinical information of patients into metadata based on the "Patient Characteristics.xlsx" in the GSE120575.
metadata$response <- ifelse(metadata$PatientID %in% 
                              c("Pre_P1", "Pre_P7", "Pre_P8", "Pre_P24", "Pre_P26", "Pre_P28", "Pre_P29",
                                "Pre_P33", "Pre_P35", "Post_P1", "Post_P4", "Post_P5_2", "Post_P7",
                                "Post_P8", "Post_P17", "Post_P19", "Post_P21"), "responder", "non-responder")

metadata$treatment <- ifelse(metadata$PatientID %in% 
                               c("Pre_P4","Pre_P7", "Pre_P8", "Pre_P26", "Pre_P28","Post_P4", "Post_P7", 
                                 "Post_P8", "Post_P13", "Post_P28", "Post_P28_2"), "aCTLA4 n aPD1", 
                             ifelse(metadata$PatientID %in% 
                                      c("Pre_P1", "Pre_P6"), "aCTLA4", "aPD1"))

#2. Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = expressiondata, meta.data = metadata)

#3. Filter for patient samples related to the study, excluding mixed-response samples
Idents(seuratobj) <- "PatientID"
aPD1_so <- subset(seurat_obj, idents = c("Pre_P2", "Pre_P6",  "Pre_P15", "Pre_P20", "Pre_P25", "Pre_P27", "Pre_P31", 
                                         "Post_P2", "Post_P6", "Post_P10", "Post_P11", "Post_P12", "Post_P14",
                                         "Post_P15", "Post_P16", "Post_P18","Post_P20","Post_P23","Post_P30",
                                          "Post_P17", "Post_P19", "Post_P21"))

#4.Clustering
aPD1_so <- FindVariableFeatures(aPD1_so)
aPD1_so<- RunPCA(aPD1_so, npcs = 50, verbose = FALSE)
ElbowPlot(aPD1_so, ndims = 50)

aPD1_so <- RunUMAP(aPD1_so, 
                  reduction = "pca", 
                  dims = 1:20)

# Determine the K-nearest neighbor graph
aPD1_so <- FindNeighbors(aPD1_so, 
                        reduction = "pca", 
                        dims = 1:20)

# Determine the clusters for various resolutions                                
aPD1_so <- FindClusters(object = aPD1_so,
                       resolution = c(0.1, 0.2), 
                       verbose = F)

aPD1_so@meta.data %>% 
  dplyr::select(contains("snn")) %>% 
  map_int(~ unique(.x) %>% length)

# Determine the DGE between each clusters with resolution 0.1 (Supplemental Table 6d)
Idents(aPD1_so) <- "RNA_snn_res.0.1"
DE_Cluster <- FindAllMarkers(aPD1_so, only.pos = T) 

FeaturePlot(aPD1_so,
            cols = c("lightgrey", "blue"),
            pt.size = 0.1,
            features = c("CD4", "FOXP3")) #Fig.S9c

#Treg annotation based on CD4 and FOXP3
Idents(aPD1_so) <- 'RNA_snn_res.0.1'
CD4_seurat <- subset(aPD1_so, idents = "2")
CD4_seurat$celltypes <- "CD4 T"
grep("^FOXP3$", rownames(CD4_seurat@assays[["RNA"]]$counts))
length(which(CD4_seurat@assays[["RNA"]]$counts[rownames(CD4_seurat@assays[["RNA"]])[678], ] != 0)) 
CD4_seurat@meta.data$celltypes[which(CD4_seurat@assays[["RNA"]]$counts[rownames(CD4_seurat@assays[["RNA"]])[678], ] != 0)] <- 'Treg'
Idents(CD4_seurat) <- "celltypes"
aPD1_so$celltypes <-""
aPD1_so$celltypes[Cells(CD4_seurat)] <- paste(Idents(CD4_seurat))

#Cell annotation
aPD1_so@meta.data$celltypes[aPD1_so@meta.data$RNA_snn_res.0.1 == '0'] <- 'CD8 T'
aPD1_so@meta.data$celltypes[aPD1_so@meta.data$RNA_snn_res.0.1 == '1'] <- 'NK' 
aPD1_so@meta.data$celltypes[aPD1_so@meta.data$RNA_snn_res.0.1 == '3'] <- 'DC'
aPD1_so@meta.data$celltypes[aPD1_so@meta.data$RNA_snn_res.0.1 == '4'] <- 'B'
aPD1_so@meta.data$celltypes[aPD1_so@meta.data$RNA_snn_res.0.1 == '5'] <- 'Mono and Mac'

DimPlot(aPD1_so,
        reduction = "umap",
        label = F,
        pt.size = 0.3,
        label.size = 6,
        group.by = "celltypes") #Fig.S9b

#Frequency (Fig.S9c)
table_celltypes <- table(aPD1_so@meta.data$PatientID, aPD1_so@meta.data$celltypes)
Frequency_celltypes <- as.data.frame(table_celltypes)
Frequency_celltypes$Percentage <- round(Frequency_celltypes$Freq / sum(Frequency_celltypes$Freq)*100, 1)
colnames(Frequency_celltypes) <- c("PatientID", "CellType", "Frequency", "Percentage (%)")

ggplot(Frequency_celltypes, aes(fill=CellType, y=Frequency, x=PatientID)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.line = element_line(size = 0.5),
  ) + 
  theme_minimal()

#Marker Identification (Fig.S9d)
jjDotPlot(aPD1_so,
          gene = c("CD3D", "CD4", "CD8A", "CD19", "CD1D", "CD80", "CD86",
                   "CD163", "CD14", "CD83","NCAM1", "ITGAX", "ITGAM", "IL7R", "CD68",
                   "EOMES", "FOXP3", "TBX21", "GATA3", "RORA", "RORC", "BCL6", "IRF4", "MAF",
                   "CTLA4", "ENTPD1", "HAVCR2", "LAG3", "PDCD1", "TIGIT",
                   "CD28", "CD40", "ICOS", "TNFRSF1B", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFRSF13B", "TNFRSF17", "TNFRSF18", "TNFRSF25",
                   "CCR7", "LEF1", "SELL", "TCF7",
                   "MKI67", "PCNA", 
                   "IL1A", "IL1B", "IL2", "IL4", "IL6", "IL10", "IL12A", "IL12B", "IL13", "IL15", "IL17A", "IL17F", "IL18", "IL21", "IL33", "IFNG", "TGFB1", "TGFB3",
                   "GZMA", "GZMB", "GZMK", "NKG7", "PRF1"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = 'celltypes')

#Violin plot visualizing gene expression in Tregs in baseline vs post aPD1 therapy (Fig.S9f-j)
#Clustered Wilcoxon rank-sum test was applied
Idents(aPD1_so) <- "celltypes"
Treg_so <- subset(aPD1_so, idents = "Treg")
Idents(Treg_so) <- "PatientID"
Nir_prepost <- subset(Treg_so, idents = c("Pre_P15", "Pre_P2", "Pre_P6","Pre_P20","Pre_P25", "Pre_P27", "Pre_P31",
                                          "Post_P2", "Post_P6", "Post_P10","Post_P11", "Post_P12", "Post_P14",
                                          "Post_P16", "Post_P23", "Post_P30", "Post_P15", "Post_P18", "Post_P20"))

Nir_prepost$preorpost <- ifelse(Nir_prepost$PatientID %in% 
                                  c("Pre_P15", "Pre_P2", "Pre_P6","Pre_P20","Pre_P25", "Pre_P27", "Pre_P31"), 
                                "Pre", "Post")

Nir_prepost$preorpost <- factor(Nir_prepost$preorpost, 
                                levels = c("Pre", "Post"))

Idents(Nir_prepost) <- "preorpost"
VlnPlot(object = Nir_prepost, 
        features = c("TNFRSF8"), pt.size = 0.8, #also visualized CTLA4, TNFRSF18, LRRC32, TIGIT
        log= F, adjust =30, alpha = 1, fill.by = "feature",y.max = 10,
        group.by = "preorpost",
        cols = c("#11B2B6","#F15F5B"),
        idents = c("Pre", "Post")) 

library(clusrank)
library(dplyr)

gene <- "TNFRSF8" #also tested CTLA4, TNFRSF18, LRRC32, TIGIT
expression_data <- FetchData(object = Nir_prepost, vars = c("preorpost", "PatientID", gene))
clusWilcox.test(TNFRSF8 ~ preorpost + cluster(PatientID), data = expression_data, method = "ds")

# Average expression of a specific gene in different individual within baseline vs post aPD1 therapy cohort (Fig.6k-o)
avg_expr <- FetchData(Nir_prepost, vars = c("TNFRSF8", "CTLA4", "LRRC32", "TIGIT", "TNFRSF18"), layer = "data")
avg_expr$PatientID <- Nir_prepost$PatientID

avg_expr_summary <- avg_expr %>%
  group_by(PatientID) %>%
  summarise(across(c("TNFRSF8", "CTLA4", "LRRC32", "TIGIT", "TNFRSF18"), mean), .groups = "drop")

# Average expression of a specific gene in different individual within responder vs non-responder (Fig.6p,S9k-n)
Idents(Treg_so) <- "PatientID"
Nir_response <- subset(Treg_so, idents = c("Post_P17", "Post_P19", "Post_P21",
                                          "Post_P2", "Post_P6", "Post_P10","Post_P11", "Post_P12", "Post_P14",
                                          "Post_P16", "Post_P23", "Post_P30", "Post_P15", "Post_P18", "Post_P20"))

avg_expr <- FetchData(Nir_response, vars = c("TNFRSF8", "CTLA4", "LRRC32", "TIGIT", "TNFRSF18"), layer = "data")
avg_expr$PatientID <- Nir_response$PatientID

avg_expr_summary <- avg_expr %>%
  group_by(PatientID) %>%
  summarise(across(c("TNFRSF8", "CTLA4", "LRRC32", "TIGIT", "TNFRSF18"), mean), .groups = "drop")
