library(Seurat)
library(dplyr)
library(Matrix)

#### USER INPUT #####
gene_threshold <- c(200, 2500) # [min, max] count
mito_threshold <- c(-Inf, 0.05) # [min, max] frequency

data_location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/merged_counts_M1_M2.csv" #@TODO add ability to do 10X sparse data
cell.label.location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/sample_labels.csv"
# cell.label.location <- NULL
##### LOAD DATA #####
# km.data <- Read10X(data.dir = data_dir)
km.data <- read.csv(file = data_location, header = FALSE)
cell.label <- read.csv(file = cell.label.location, header = FALSE)
cell.label <- cell.label[-1]
cell.label <- t(cell.label)
rownames(cell.label) <- colnames(km.data)
new_rownames <- make.names(km.data[ , 1], unique = TRUE)
rownames(km.data) <- new_rownames
km.data <- km.data[ ,-1]

#
km <- CreateSeuratObject(raw.data = km.data, min.cells = 3, min.genes = 200, project = "GATA_WT_MERGED")

##### QC #####

# UMIs mapping to mitochondrial genes
# note: check to make sure that mt genes are denoted as mt.atp and not as mt-atp
mito.genes <- grep(pattern = "^mt\\.", x = rownames(x = km@data), value = TRUE)
percent.mito <- Matrix::colSums(km@raw.data[mito.genes, ]) / Matrix::colSums(km@raw.data)

km <- AddMetaData(object = km, metadata = percent.mito, col.name = "percent.mito")

km <- AddMetaData(object = km, metadata = cell.label, col.name = "cell.label")
VlnPlot(object = km, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1,2))
GenePlot(object = km, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = km, gene1 = "nUMI", gene2 = "nGene")

km <- FilterCells(object = km, subset.names = c("nGene", "percent.mito"),
                  low.thresholds = c(gene_threshold[1], mito_threshold[1]),
                  high.thresholds = c(gene_threshold[2], mito_threshold[2]))

##### NORMALIZE #####
km <- NormalizeData(object = km, normalization.method = "LogNormalize", scale.factor = 1e4)

##### VARIABLE GENES #####
km <- FindVariableGenes(object = km, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
 #@TODO check these parameters to see if they make sense for our data

print("Number of variable genes:", quote = FALSE)
length(x = km@var.genes)

##### SCALE DATA AND REGRESSION #####

km <- ScaleData(object = km, vars.to.regress = c("nUMI", "percent.mito")) #@SPEED

##### PCA DIMENSIONALITY REDUCTION #####
km <- RunPCA(object = km, pc.genes = km@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)

# PrintPCA(object = km, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
# VizPCA(object = km, pcs.use = 1:2)
PCAPlot(object = km, dim.1 = 1, dim.2 = 2)

# Jackstraw method
km <- JackStraw(object = km, num.replicate = 100, do.print = FALSE) #@SPEED
 
JackStrawPlot(object = km, PCs = 1:20)

# PC Elbow plot
PCElbowPlot(object = km)

##### CLUSTERING #####
num_of_pca_components_to_use = 10
km <- FindClusters(object = km, reduction.type = "pca", dims.use = 1:num_of_pca_components_to_use, resolution = 0.6, print.output = 0, save.SNN = TRUE)
km <- FindClusters(object = km, reduction.type = "pca", dims.use = 1:num_of_pca_components_to_use, resolution = 1.2, print.output = 0, save.SNN = TRUE)


##### SAVE OBJECT #####

save(km, file = "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/km_merged_v1.Robj")

##### tSNE #####
km <- RunTSNE(object = km, dims.use = 1:num_of_pca_components_to_use, do.fast = TRUE)
TSNEPlot(object = km, do.label = TRUE)

TSNEPlot(object = km, do.label = TRUE, group.by = "cell.label")

##### FIND ALL MARKERS OF CLUSTERS #####
HSCcluster.markers <- FindMarkers(object = km, ident.1 = "thrombocytes", min.pct = 0.25)
write.table(HSCcluster.markers,
            file = "/Users/jon/data_analysis/10X_scRNAseq/seurat/gfp_markers.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            na = "\t")


km.markers <- FindAllMarkers(object = km, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top20 <- km.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)


write.table(top20,
            file = "/Users/jon/data_analysis/10X_scRNAseq/seurat/top20_markers.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            na = "\t")


top10 <- km.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top5 <- km.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

##### VIOLIN PLOT OF USER GENES ACROSS CLUSTERS #####
violin_gene_list <- c("gene1", "gene2") 
VlnPlot(object = km, features.plot = violin_gene_list)

##### FEATURE TSNE PLOTS #####
tsne_feature_list <- c("gata1b")
colors_to_use <- c("grey", "blue") #scaled from first color in list to second color in list
FeaturePlot(object = km, features.plot = tsne_feature_list, cols.use = colors_to_use, reduction.use = "tsne")

##### HEATMAP OF MARKER GENES #####
DoHeatmap(object = km, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE)

##### IDENTIFY CLUSTERS #####
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
new.cluster.ids <- c("progenitor_1",
                     "progenitor_2",
                     "progenitor_3",
                     "erythrocyte_1",
                     "erythrocyte_2",
                     "B_cell",
                     "progenitor_4",
                     "progenitor_5",
                     "neutrophil_1",
                     "macrophage",
                     "neutrophil_2",
                     "kidney_progenitors_1",
                     "kidney_mucin",
                     "thrombocytes",
                     "kidney_progenitors_2")
km@ident <- plyr::mapvalues(x = km@ident, from = current.cluster.ids, to = new.cluster.ids)

tsne_feature_list <- c("npm1a",
                       "alas2",
                       "lyz",
                       "marco",
                       "meis1b",
                       "dap1b",
                       "EGFP")
colors_to_use <- c("lightgrey", "blue") #scaled from first color in list to second color in list
FeaturePlot(object = km, features.plot = tsne_feature_list, cols.use = colors_to_use, reduction.use = "tsne", pt.size = 0.5)


