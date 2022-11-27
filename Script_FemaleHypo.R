library(devtools)
library(Seurat)
library(Matrix)
library(xlsx)
library(umap)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)

## Setup the Seurat Object- intact female
ifemale.data <- Read10X(data.dir = "../P29_intact_F/outs/filtered_feature_bc_matrix/")
ifemale.atleastone <- apply(ifemale.data, 2, function(x) sum(x>0))
hist(ifemale.atleastone, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
ifemale.tmp <- apply(ifemale.data, 1, function(x) sum(x>0))
table(ifemale.tmp>=3)
ifemale.keep <- ifemale.tmp>=3
ifemale.tmp <- ifemale.data[ifemale.keep,]
ifemale.atleastone <- apply(ifemale.tmp, 2, function(x) sum(x>0))
summary(ifemale.atleastone)
rownames(x = ifemale.data) <- gsub(pattern = '_', replacement = '-', x = rownames(x = ifemale.data))
ifemale <- CreateSeuratObject(counts = ifemale.data, min.cells = 3, min.features = 200, project = "Intact_female")

## Setup the Seurat Object-eb300 female
efemale.data <- Read10X(data.dir = "../P29_EB300_F/outs/filtered_feature_bc_matrix/")
efemale.atleastone <- apply(efemale.data, 2, function(x) sum(x>0))
hist(efemale.atleastone, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
efemale.tmp <- apply(efemale.data, 1, function(x) sum(x>0))
table(efemale.tmp>=3)
efemale.keep <- efemale.tmp>=3
efemale.tmp <- efemale.data[efemale.keep,]
efemale.atleastone <- apply(efemale.tmp, 2, function(x) sum(x>0))
summary(efemale.atleastone)
rownames(x = efemale.data) <- gsub(pattern = '_', replacement = '-', x = rownames(x = efemale.data))
efemale <- CreateSeuratObject(counts = efemale.data, min.cells = 3, min.features = 200, project = "EB300_female")

## Set up WT object
ctrl <- CreateSeuratObject(counts = ifemale.data, project = "Intact_female", min.cells = 5)
ctrl$stim <- "Intact"
ctrl <- subset(x = ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(object = ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000)

## Set up KO object
stim <- CreateSeuratObject(counts = efemale.data, project = "EB300_female", min.cells = 5)
stim$stim <- "EB300"
stim <- subset(x = stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(object = stim, verbose = FALSE)
stim <- FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 2000)

# Find Anchors
anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = combined) <- "integrated"

# Standard workflow for visualization and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.7)
DefaultAssay(object = combined) <- "RNA"

## Plotting
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(object = combined, reduction = "umap", split.by = "stim")

FeaturePlot(object = combined, features = c("Msx2"),  
            cols = c("grey", "blue"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = combined, features = c("Amh"),  
            cols = c("grey", "blue"), blend.threshold = 0.5, pt.size = 0.5)

## cluster markers
C0.markers <- FindConservedMarkers(object = combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
head(x = C0.markers)
C1.markers <- FindConservedMarkers(object = combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
head(x = C1.markers)
C2.markers <- FindConservedMarkers(object = combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
head(x = C2.markers)
C3.markers <- FindConservedMarkers(object = combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
head(x = C3.markers)
C4.markers <- FindConservedMarkers(object = combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
head(x = C4.markers)
C5.markers <- FindConservedMarkers(object = combined, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
head(x = C5.markers)
C6.markers <- FindConservedMarkers(object = combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(x = C6.markers)
C7.markers <- FindConservedMarkers(object = combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(x = C7.markers)
C8.markers <- FindConservedMarkers(object = combined, ident.1 = 8, grouping.var = "stim", verbose = FALSE)
head(x = C8.markers)
C9.markers <- FindConservedMarkers(object = combined, ident.1 = 9, grouping.var = "stim", verbose = FALSE)
head(x = C9.markers)
C10.markers <- FindConservedMarkers(object = combined, ident.1 = 10, grouping.var = "stim", verbose = FALSE)
head(x = C10.markers)
C11.markers <- FindConservedMarkers(object = combined, ident.1 = 11, grouping.var = "stim", verbose = FALSE)
head(x = C11.markers)
C12.markers <- FindConservedMarkers(object = combined, ident.1 = 12, grouping.var = "stim", verbose = FALSE)
head(x = C12.markers)
C13.markers <- FindConservedMarkers(object = combined, ident.1 = 13, grouping.var = "stim", verbose = FALSE)
head(x = C13.markers)
C14.markers <- FindConservedMarkers(object = combined, ident.1 = 14, grouping.var = "stim", verbose = FALSE)
head(x = C14.markers)
C15.markers <- FindConservedMarkers(object = combined, ident.1 = 15, grouping.var = "stim", verbose = FALSE)
head(x = C15.markers)
C16.markers <- FindConservedMarkers(object = combined, ident.1 = 16, grouping.var = "stim", verbose = FALSE)
head(x = C16.markers)
C17.markers <- FindConservedMarkers(object = combined, ident.1 = 17, grouping.var = "stim", verbose = FALSE)
head(x = C17.markers)
C18.markers <- FindConservedMarkers(object = combined, ident.1 = 18, grouping.var = "stim", verbose = FALSE)
head(x = C18.markers)
C19.markers <- FindConservedMarkers(object = combined, ident.1 = 19, grouping.var = "stim", verbose = FALSE)
head(x = C19.markers)
C20.markers <- FindConservedMarkers(object = combined, ident.1 = 20, grouping.var = "stim", verbose = FALSE)
head(x = C20.markers)
C21.markers <- FindConservedMarkers(object = combined, ident.1 = 21, grouping.var = "stim", verbose = FALSE)
head(x = C21.markers)

## Identification of cell types
FeaturePlot(object = combined, features = c("Snap25", "Agt", "Slc32a1", "Slc17a6"),  
            cols = c("grey", "blue"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = combined, features = c("Pdgfra", "Fyn", "Mobp", "Olig1"),  
            cols = c("grey", "blue"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = combined, features = c("Tmem119", "Tmem212", "Mgp", "Pf4"),
            cols = c("grey", "blue"), blend.threshold = 0.5, pt.size = 0.5)

## Identification of KISS cells
FeaturePlot(object = combined, features = c("Kiss1"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = combined, features = c("Kiss1", "Pomc", "Esr1", "Esr2"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1)
FeaturePlot(object = combined, features = c("Vim", "Tac3", "Cfap126", "Crym"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1)
combined$stim <- factor(x = combined$stim, levels = c('Intact', 'EB300'))

# UMAP ploting with cluster numbers
DimPlot(object = combined, reduction = "umap", label = TRUE, label.size = 5)
DimPlot(object = combined, reduction = "umap", split.by = "stim")
DimPlot(object = combined, reduction = "umap", split.by = "stim", label = TRUE, label.size = 5)
FeaturePlot(object = combined, features = c("Kiss1"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 1.5, label = TRUE, label.size = 3)
DotPlot(object = combined, features = c("Tac3", "Pdyn", "Kiss1"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6)
FeaturePlot(object = combined, features = c("Hsd17b4"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
FeaturePlot(object = combined, features = c("Hsd17b8"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
FeaturePlot(object = combined, features = c("Hsd17b1"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
FeaturePlot(object = combined, features = c("Hsd17b7"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
VlnPlot(combined, features = c("Hsd17b4"), split.by = "stim", 
        pt.size = 0.1, combine = FALSE)
VlnPlot(combined, features = c("Hsd17b8"), split.by = "stim", 
        pt.size = 0.1, combine = FALSE)
VlnPlot(combined, features = c("Hsd17b1"), split.by = "stim", 
        pt.size = 0.1, combine = FALSE)
VlnPlot(combined, features = c("Hsd17b7"), split.by = "stim", 
        pt.size = 0.1, combine = FALSE)

## Sub-clustering (Kiss1+ cluster; 6)
C_Kiss1 <- subset(x = combined, idents = c("6"))
C_Kiss1
C_Kiss1$stim <- factor(x = C_Kiss1$stim, levels = c('Intact', 'EB300'))

## Plotting Sub-cluster (Kiss1+ cluster; 6)
p1 <- DimPlot(object = C_Kiss1, reduction = "umap", group.by = "stim", label.size = 5)
p2 <- DimPlot(object = C_Kiss1, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1)
plot_grid(p1, p2)
DimPlot(object = C_Kiss1, reduction = "umap", split.by = "stim", label.size = 5)
FeaturePlot(object = C_Kiss1, features = c("Kiss1", "Vim", "Pomc", "Esr1"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
table(C_Kiss1@active.ident)
VlnPlot(C_Kiss1, features = c("Kiss1"), 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Tac3"), 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Pdyn"), 
        pt.size = 0.5, combine = FALSE)
RidgePlot(object = C_Kiss1, feature = c("Kiss1"), group.by = "stim")
DotPlot(object = combined, features = c("Kiss1", "Vim", "Pomc", "Esr1"), split.by = "stim")
DotPlot(object = combined, features = c("Kiss1", "Esr1"), split.by = "stim")
DotPlot(object = combined, features = c("Esr2", "Esr1", "Kiss1"))
DotPlot(object = combined, features = c("Gper1", "Esrrg", "Esrra", "Esr2", "Esr1"))
DotPlot(object = combined, features = c("Agrp", "Npy", "Pomc", "Kiss1"))

FeatureScatter(object = combined, feature1 = "Kiss1", feature2 = "Vim")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Esr1")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Npy")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Pomc")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Ar")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Slc7a11", group.by = "stim")
FeatureScatter(object = combined, feature1 = "Kiss1", feature2 = "Npy")
FeatureScatter(object = combined, feature1 = "Kiss1", feature2 = "Agrp")
FeatureScatter(object = combined, feature1 = "Kiss1", feature2 = "Pomc")

C_Kiss1.Kiss1 <- subset(x = C_Kiss1, subset = Kiss1 > 1)
C_Kiss1.Kiss1
C_Kiss1.Kiss1$stim <- factor(x = C_Kiss1.Kiss1$stim, levels = c('Intact', 'EB300'))
table(C_Kiss1.Kiss1@active.ident)
table(C_Kiss1.Kiss1@meta.data$stim)

table(combined@meta.data$stim)
table(C_Kiss1@meta.data$stim)

FeaturePlot(object = C_Kiss1, features = c("Cox6a2"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = C_Kiss1, features = c("Gstm1", "Gsta4"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = C_Kiss1, features = c("Gstm2", "Gsta1"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
VlnPlot(C_Kiss1, features = c("Cox6a2"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Gstm1"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Gsta4"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Gstm2"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Gsta1"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Esrra"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)
VlnPlot(C_Kiss1, features = c("Esrrg"), group.by = "stim", 
        pt.size = 0.5, combine = FALSE)

DotPlot(object = combined, features = c("Gstm1", "Gsta4", "Gstm2", "Gsta1"))

FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Gstm1", group.by = "stim")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Gsta4", group.by = "stim")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Gstm2", group.by = "stim")
FeatureScatter(object = C_Kiss1, feature1 = "Kiss1", feature2 = "Gsta1", group.by = "stim")
FeatureScatter(object = C_Kiss1, feature1 = "Pomc", feature2 = "Esrra", group.by = "stim")

# Differential gene expression 
Idents(object = C_Kiss1) <- "stim"
avg.C_Kiss1 <- log1p(x = AverageExpression(object = C_Kiss1, verbose = FALSE)$RNA)
avg.C_Kiss1$gene <- rownames(x = avg.C_Kiss1)
Idents(object = C_Kiss1) <- "celltype.stim"
C_Kiss1$celltype <- Idents(object = C_Kiss1)
C_Kiss1$celltype.stim <- paste(Idents(object = C_Kiss1), C_Kiss1$stim, sep = "_")
head(C_Kiss1$celltype.stim)
Idents(object = C_Kiss1) <- "celltype.stim"
response.C_Kiss1 <- FindMarkers(object = C_Kiss1, ident.1 = c("celltype.stim_Intact"), ident.2 = c("celltype.stim_EB300"), verbose = FALSE)
head(x = response.C_Kiss1, n = 15)

genes.to.label =response.C_Kiss1 %>%
  tibble::rownames_to_column("gene_name") %>% #gene names was in rownames, need to put in a seprate column
  dplyr::filter(p_val_adj<0.05) %>%  #select only significant genes #keep only significant genes
  arrange(desc(abs(avg_logFC))) %>%  #sort based on the obsolute FC value, from highest to lowesr
  head(n=20) %>%  #look at top 10
  pull(gene_name) # extract gene names
genes.to.label

table(response.C_Kiss1$p_val_adj<0.05)
# only mark genes with FDR p-value < 0.05
sig_genes <- response.C_Kiss1[response.C_Kiss1$p_val_adj<0.05,]  
#make common coloum "gene" for joining
sig_genes$gene <- rownames(sig_genes) 
#need to merge test results and epression table together so everything is in the same order
avg.C_Kiss1 <- avg.C_Kiss1 %>% left_join(sig_genes) 

avg.C_Kiss1$status <- "NoSig"
avg.C_Kiss1$status[avg.C_Kiss1$avg_logFC>(0) & avg.C_Kiss1$p_val_adj<0.05] <- "Up"   
avg.C_Kiss1$status[avg.C_Kiss1$avg_logFC<(0) & avg.C_Kiss1$p_val_adj<0.05] <- "Down"
avg.C_Kiss1$status <- factor(avg.C_Kiss1$status, levels = c("NoSig", "Up", "Down"))
#re-arrange the df so "NoSig" points got drew first
avg.C_Kiss1 <- avg.C_Kiss1 %>% arrange(status) 
#re-order the levels to make the legend show in ggplot
avg.C_Kiss1$status <- factor(avg.C_Kiss1$status, levels = c( "Up","NoSig", "Down")) 
#for label points to work
rownames(avg.C_Kiss1) <- avg.C_Kiss1$gene 

pp1 <- ggplot(avg.C_Kiss1, aes(EB300, Intact, color = status)) + geom_point() + ggtitle("DEG")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Plot - DEG genes
FeaturePlot(object = combined, features = c("Vim"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
FeaturePlot(object = combined, features = c("Rorb"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
FeaturePlot(object = combined, features = c("Sfrp2"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")
FeaturePlot(object = combined, features = c("Serpinf1"), 
            cols = c("lightgrey", "red"), blend.threshold = 0.5, pt.size = 0.1, split.by = "stim")