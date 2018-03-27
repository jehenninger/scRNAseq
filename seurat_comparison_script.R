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


save(combined, file = "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/km_comparison.Robj")

## this doesn't work
# conserved_markers <- list()
# count <- 1
# for(i in 0:18){
#   conserved_markers[count] <- FindConservedMarkers(combined, ident.1 = i, grouping.var = "genotype", 
#                                               print.bar = FALSE)
#   count <- count + 1
# }


##### Quantify clusters by genotype #####
combined <- SetAllIdent(object = combined, id = "celltype")
count_table <- table(combined@ident, combined@meta.data$genotype)
freq_table <- prop.table(x = count_table, margin = 2)
count_table
freq_table

tsne_feature_list <- c("myca",
                       "npm1a",
                       "alas2",
                       "lyz",
                       "marco",
                       "kdrl",
                       "meis1b",
                       "dap1b",
                       "lck",
                       "dusp2",
                       "zap70")
colors_to_use <- c("lightgrey", "blue") #scaled from first color in list to second color in list
FeaturePlot(object = combined, features.plot = tsne_feature_list, cols.use = colors_to_use, reduction.use = "tsne", pt.size = 0.5)


##### Name clusters #####
new.ident <- c("neutrophil_1",
               "progenitor_1",
               "erythrocytes_1",
               "erythrocytes_2",
               "T",
               "B_1",
               "progenitor_2",
               "neutrophil_2",
               "B_2",
               "macrophage",
               "vascular",
               "unknown_1",
               "unknown_2",
               "progenitor_3",
               "unknown_3",
               "HSC_1",
               "unknown_4",
               "HSC_2",
               "unknown_5")

for (i in 0:18) {
  combined <- RenameIdent(object = combined, old.ident.name = i, 
                                 new.ident.name = new.ident[i + 1])
}

##### Compare clusters between groups #####
combined@meta.data$celltype.genotype <- paste0(combined@ident, "_", 
                                                  combined@meta.data$genotype)
combined <- StashIdent(combined, save.name = "celltype")
combined <- SetAllIdent(combined, id = "celltype.genotype")

combined <- SetAllIdent(combined, id = "genotype")
overall_gata_vs_wt <- FindMarkers(combined, ident.1 = "GATA", ident.2 = "WT", print.bar = TRUE)


neutrophil_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "neutrophil_1_GATA", ident.2 = "neutrophil_1_WT", print.bar = FALSE)

neutrophil_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "neutrophil_2_GATA", ident.2 = "neutrophil_2_WT", print.bar = FALSE)

progenitor_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "progenitor_1_GATA", ident.2 = "progenitor_1_WT", print.bar = FALSE)
progenitor_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "progenitor_2_GATA", ident.2 = "progenitor_2_WT", print.bar = FALSE)
progenitor_3_gata_vs_wt <- FindMarkers(combined, ident.1 = "progenitor_3_GATA", ident.2 = "progenitor_3_WT", print.bar = FALSE)

hsc_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "HSC_1_GATA", ident.2 = "HSC_1_WT", print.bar = FALSE)
hsc_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "HSC_2_GATA", ident.2 = "HSC_2_WT", print.bar = FALSE)

erythrocytes_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "erythrocytes_1_GATA", ident.2 = "erythrocytes_1_WT", print.bar = FALSE)
erythrocytes_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "erythrocytes_2_GATA", ident.2 = "erythrocytes_2_WT", print.bar = FALSE)

compared_feature_list <- c("lyz", "mpx")
FeatureHeatmap(combined, features.plot = compared_feature_list, group.by = "genotype", pt.size = 2, key.position = "top", max.exp = 3)




##### Scatter plots of group comparison #####
neutrophils <- SubsetData(combined, ident.use = c(0,7), subset.raw = T)
neutrophils <- SetAllIdent(neutrophils, id = "genotype")
avg.t.cells <- log1p(AverageExpression(t.cells, show.progress = FALSE))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- SubsetData(immune.combined, ident.use = "CD14 Mono", subset.raw = T)
cd14.mono <- SetAllIdent(cd14.mono, id = "stim")
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, show.progress = FALSE))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label1 = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1")
genes.to.label2 = c("IFIT2", "IFIT1")
genes.to.label3 = c("CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelUR(p1, genes = c(genes.to.label1, genes.to.label2), avg.t.cells, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p1 <- LabelUL(p1, genes = genes.to.label3, avg.t.cells, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelUR(p2, genes = c(genes.to.label1, genes.to.label3), avg.cd14.mono, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = genes.to.label2, avg.cd14.mono, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
plot_grid(p1, p2)