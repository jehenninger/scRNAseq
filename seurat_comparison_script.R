library(Seurat)
library(dplyr)
library(Matrix)

#### USER INPUT #####
gene_threshold <- c(200, 2500) # [min, max] count
mito_threshold <- c(-Inf, 0.05) # [min, max] frequency

sample_A_location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/JH_M2/JH_M2.counts.tr.csv" #@TODO add ability to do 10X sparse data
sample_B_location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/JH_M1/JH_M1.counts.tr.csv"

#cell.label.location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/sample_labels.csv"
# cell.label.location <- NULL
##### LOAD DATA #####
# km.data <- Read10X(data.dir = data_dir)

# Set up sample A
sample_A.data <- read.csv(file = sample_A_location, header = FALSE)
new_rownames <- make.names(sample_A.data[ , 1], unique = TRUE)
rownames(sample_A.data) <- new_rownames
sample_A.data <- sample_A.data[ ,-1]
sample_A <- CreateSeuratObject(raw.data = sample_A.data, project = "GATA", min.cells = 3)
sample_A@meta.data$genotype <- "WT"
sample_A <- FilterCells(sample_A, subset.names = "nGene", low.thresholds = gene_threshold[1], high.thresholds = gene_threshold[2])
sample_A <- NormalizeData(sample_A)
sample_A <- ScaleData(sample_A, display.progress = T)

# Set up sample B
sample_B.data <- read.csv(file = sample_B_location, header = FALSE)
new_rownames <- make.names(sample_B.data[ , 1], unique = TRUE)
rownames(sample_B.data) <- new_rownames
sample_B.data <- sample_B.data[ ,-1]
sample_B <- CreateSeuratObject(raw.data = sample_B.data, project = "GATA", min.cells = 3)
sample_B@meta.data$genotype <- "GATA"
sample_B <- FilterCells(sample_B, subset.names = "nGene", low.thresholds = gene_threshold[1], high.thresholds = gene_threshold[2])
sample_B <- NormalizeData(sample_B)
sample_B <- ScaleData(sample_B, display.progress = T)

#Gene selection for input to CCA
sample_A <- FindVariableGenes(sample_A, do.plot = F)
sample_B <- FindVariableGenes(sample_B, do.plot = F)
g.1 <- head(rownames(sample_A@hvg.info), 1000)
g.2 <- head(rownames(sample_B@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(sample_A@scale.data))
genes.use <- intersect(genes.use, rownames(sample_B@scale.data))

# Canonical correlation analysis
combined <- RunCCA(sample_A, sample_B, genes.use = genes.use, num.cc = 30, add.cell.id1 = "_WT", add.cell.id2 = "_GATA")

p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "genotype",
              pt.size = 0.5, do.return = TRUE)

p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "genotype",
              do.return = TRUE)

plot_grid(p1,p2)

PrintDim(object = combined, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

p3 <- MetageneBicorPlot(combined, grouping.var = "genotype", dims.eval = 1:30, display.progress = FALSE)

combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "genotype",
                          dims.align = 1:20)

combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
combined <- FindClusters(combined, reduction.type = "cca.aligned",
                         resolution = 0.6, dims.use = 1:20)

p1 <- TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "genotype")
p2 <- TSNEPlot(combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1,p2)

conserved_markers <- list()
for(i in 1:19){
  conserved_markers[i] <- FindConservedMarkers(combined, ident.1 = i, grouping.var = "genotype", 
                                              print.bar = FALSE)
}
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

##### SAVE OBJECT #####

save(km, file = "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/km_merged_v1.Robj")

##### tSNE #####
km <- RunTSNE(object = km, dims.use = 1:num_of_pca_components_to_use, do.fast = TRUE)
TSNEPlot(object = km, do.label = TRUE)

TSNEPlot(object = km, do.label = TRUE, group.by = "cell.label")

##### FIND ALL MARKERS OF CLUSTERS #####
HSCcluster.markers <- FindMarkers(object = km, ident.1 = "HSCs thrombocytes2", min.pct = 0.25)

km.markers <- FindAllMarkers(object = km, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- km.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

##### VIOLIN PLOT OF USER GENES ACROSS CLUSTERS #####
violin_gene_list <- c("gene1", "gene2") 
VlnPlot(object = km, features.plot = violin_gene_list)

##### FEATURE TSNE PLOTS #####
tsne_feature_list <- c("gata1b")
colors_to_use <- c("grey", "blue") #scaled from first color in list to second color in list
FeaturePlot(object = km, features.plot = tsne_feature_list, cols.use = colors_to_use, reduction.use = "tsne")

##### HEATMAP OF MARKER GENES #####
DoHeatmap(object = km, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

##### IDENTIFY CLUSTERS #####
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids <- c("myeloid2",
                     "progenitors",
                     "myeloid1",
                     "erythrocytes2",
                     "erythrocytes1",
                     "B-cells",
                     "macrophage1",
                     "HSCs thrombocytes1",
                     "macrophage2",
                     "neutrophil",
                     "kidney",
                     "unknown",
                     "HSCs thrombocytes2")
km@ident <- plyr::mapvalues(x = km@ident, from = current.cluster.ids, to = new.cluster.ids)
