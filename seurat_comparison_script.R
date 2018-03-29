library(Seurat)
library(dplyr)
library(Matrix)


#### USER INPUT #####
gene_threshold <- c(200, 2500) # [min, max] count
mito_threshold <- c(-Inf, 0.05) # [min, max] frequency

sample_A_location <- "/Users/jon/data_analysis/indrop_single_cell/JH001_transposed.csv" #@TODO add ability to do 10X sparse data
sample_B_location <- "/Users/jon/data_analysis/indrop_single_cell/JH-003_transposed.csv"

output_data_path <- "/Users/jon/data_analysis/indrop_single_cell/seurat/km_comparison/"

#cell.label.location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/sample_labels.csv"
# cell.label.location <- NULL
##### LOAD DATA #####
# km.data <- Read10X(data.dir = data_dir)

# sample_locations <- c("")
sample_labels <- c("Runx_WKM", "Zbow_WKM")
num_of_samples <- length(sample_locations)

# samples <- list()
# for(j in 1:num_of_samples){
#   samples[[j]] <- read.csv(file = sample_locations[j], header = FALSE)
#   new_rownames <- make.names(samples[[j]][ , 1], unique = TRUE)
#   rownames(samples[[j]]) <- new_rownames
#   samples[[j]] <- samples[[j]][ ,-1]
#   samples[[j]] <- CreateSeuratObject(raw.data = samples[[j]], project = "PILOT", min.cells = 3)
#   mito.genes <- grep(pattern = "^mt\\.", x = rownames(x = samples[[j]]@data), value = TRUE)
#   percent.mito <- Matrix::colSums(samples[[j]]@raw.data[mito.genes, ]) / Matrix::colSums(samples[[j]]@raw.data)
#   samples[[j]] <- AddMetaData(object = samples[[j]], metadata = percent.mito, col.name = "percent.mito")
#   samples[[j]]@meta.data$genotype <- sample_labels[j]
#   samples[[j]] <- FilterCells(samples[[j]], subset.names = c("nGene", "percent.mito"),
#                           low.thresholds = c(gene_threshold[1], mito_threshold[1]),
#                           high.thresholds = c(gene_threshold[2], mito_threshold[2]))
#   samples[[j]] <- NormalizeData(samples[[j]])
#   samples[[j]] <- ScaleData(samples[[j]], display.progress = T)
# }

# # Set up sample A
sample_A.data <- read.csv(file = sample_A_location, header = FALSE)
new_rownames <- make.names(sample_A.data[ , 1], unique = TRUE)
rownames(sample_A.data) <- new_rownames
sample_A.data <- sample_A.data[ ,-1]
sample_A <- CreateSeuratObject(raw.data = sample_A.data, project = "PILOT", min.cells = 3)
sample_A.mito.genes <- grep(pattern = "^mt\\.", x = rownames(x = sample_A@data), value = TRUE)
sample_A.percent.mito <- Matrix::colSums(sample_A@raw.data[sample_A.mito.genes, ]) / Matrix::colSums(sample_A@raw.data)
sample_A <- AddMetaData(object = sample_A, metadata = sample_A.percent.mito, col.name = "percent.mito")
sample_A@meta.data$genotype <- sample_labels[1]
sample_A <- FilterCells(sample_A, subset.names = c("nGene", "percent.mito"),
                        low.thresholds = c(gene_threshold[1], mito_threshold[1]),
                        high.thresholds = c(gene_threshold[2], mito_threshold[2]))
sample_A <- NormalizeData(sample_A)
sample_A <- ScaleData(sample_A, display.progress = T)

# Set up sample B
sample_B.data <- read.csv(file = sample_B_location, header = FALSE)
new_rownames <- make.names(sample_B.data[ , 1], unique = TRUE)
rownames(sample_B.data) <- new_rownames
sample_B.data <- sample_B.data[ ,-1]
sample_B <- CreateSeuratObject(raw.data = sample_B.data, project = "PILOT", min.cells = 3)
sample_B.mito.genes <- grep(pattern = "^mt\\.", x = rownames(x = sample_B@data), value = TRUE)
sample_B.percent.mito <- Matrix::colSums(sample_B@raw.data[sample_B.mito.genes, ]) / Matrix::colSums(sample_B@raw.data)
sample_B <- AddMetaData(object = sample_B, metadata = sample_B.percent.mito, col.name = "percent.mito")
sample_B@meta.data$genotype <- sample_labels[2]
sample_B <- FilterCells(sample_B, subset.names = c("nGene", "percent.mito"),
                        low.thresholds = c(gene_threshold[1], mito_threshold[1]),
                        high.thresholds = c(gene_threshold[2], mito_threshold[2]))
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

##### STOP AND CHOOSE NUMBER OF DIMENSIONS FOR ALIGNING DATA #####

num_dims_to_align <- 15
combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "genotype",
                          dims.align = 1:num_dims_to_align)

combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:num_dims_to_align, do.fast = T)
combined <- FindClusters(combined, reduction.type = "cca.aligned",
                         resolution = 0.6, dims.use = 1:num_dims_to_align)

num_of_clusters <- length(unique(combined@meta.data$res.0.6))

p1 <- TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "genotype")
p2 <- TSNEPlot(combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1,p2)


save(combined, file = paste0(output_data_path, "km_samples_comparison.Robj"))

conserved_markers <- list()
count <- 1
for(i in 0:(num_of_clusters - 1)){
  conserved_markers[[count]] <- FindConservedMarkers(combined, ident.1 = i, grouping.var = "genotype", print.bar = TRUE)
  write.table(conserved_markers[[count]],
              file = paste0(output_data_path,
                            "cluster_",i,"_conserved_markers.tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = NA,
              na = "\t")
  count <- count + 1
  
}

##### STOP HERE AND DECIDE ON CLUSTER NOMENCLATURE #####
##### Name clusters #####
new.ident <- c("neutrophil",
               "erythrocyte_2",
               "erythrocyte_3",
               "T_cells",
               "progenitor_2",
               "progenitor_1",
               "B_cell_1",
               "B_cell_2",
               "vascular",
               "kidney_mucin",
               "erythrocyte_1",
               "progenitor_4",
               "macrophage_1",
               "thrombocytes",
               "macrophage_2",
               "kidney_prox",
               "progenitor_5",
               "kidney_progenitor")

for (i in 0:(length(new.ident)-1)) {
  combined <- RenameIdent(object = combined, old.ident.name = i, 
                          new.ident.name = new.ident[i + 1])
}

# write markers to cluster-specific labels
count <- 1
for(i in 0:(num_of_clusters - 1)){
  write.table(conserved_markers[[count]],
              file = paste0("/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/markers/",
                            new.ident[count],"_conserved_markers.tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = NA,
              na = "\t")
  count <- count + 1
  
}

# create 'celltype' ID
combined@meta.data$celltype.genotype <- paste0(combined@ident, "_", 
                                               combined@meta.data$genotype)
combined <- StashIdent(combined, save.name = "celltype")


##### Quantify clusters by genotype #####
combined <- SetAllIdent(object = combined, id = "celltype")
count_table <- table(combined@ident, combined@meta.data$genotype)
freq_table <- prop.table(x = count_table, margin = 2)
count_table
freq_table
write.table(count_table,
            file = paste0("/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/",
                          "genotype_cluster_counts.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            na = "\t")

write.table(freq_table,
            file = paste0("/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/",
                          "genotype_cluster_percentages.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            na = "\t")


tsne_feature_list <- c("myca",
                       "npm1a",
                       "alas2",
                       "lyz",
                       "marco",
                       "meis1b",
                       "dap1b",
                       "lck",
                       "dusp2")
colors_to_use <- c("lightgrey", "blue") #scaled from first color in list to second color in list
FeaturePlot(object = combined, features.plot = tsne_feature_list, cols.use = colors_to_use, reduction.use = "tsne", pt.size = 0.5)



##### Compare clusters between groups #####
combined <- SetAllIdent(combined, id = "celltype.genotype")
params_to_compare <- sort(unique(combined@meta.data$celltype.genotype))
num_of_params <- length(params_to_compare)

output <- list()
count <- 1
for(k in seq(from = 1, to = num_of_params, by = 2)){
  output[[count]] <- FindMarkers(combined, ident.1 = params_to_compare[k], ident.2 = params_to_compare[k+1], print.bar = TRUE)
  write.table(output[[count]],
              file = paste0("/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/cluster_comparisons/",
                            params_to_compare[k],"_vs_", params_to_compare[k+1],".tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = NA,
              na = "\t")
  count <- count + 1
}

##### OLD #######
combined <- SetAllIdent(combined, id = "genotype")
overall_gata_vs_wt <- FindMarkers(combined, ident.1 = "GATA", ident.2 = "WT", print.bar = TRUE)
write.table(overall_gata_vs_wt,
            file = "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-1/project/seurat/overall_gata_vs_wt.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            na = "\t")

# heatmap the top genes changed in GATA vs WT per cluster
top_genes <- sort(overall_gata_vs_wt$avg_logFC, decreasing = TRUE, index.return = TRUE)
top_genes <- rownames(overall_gata_vs_wt[top_genes$ix,])
top_genes <- top_genes[1:20]

bottom_genes <- sort(overall_gata_vs_wt$avg_logFC, decreasing = FALSE, index.return = TRUE)
bottom_genes <- rownames(overall_gata_vs_wt[bottom_genes$ix,])
bottom_genes <- bottom_genes[1:20]

genes_to_heatmap <- c(top_genes, bottom_genes)
DoHeatmap(object = combined, genes.use = genes_to_heatmap, group.by = "genotype", slim.col.label = TRUE, remove.key = TRUE)


neutrophil_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "neutrophil_1_GATA", ident.2 = "neutrophil_1_WT", print.bar = FALSE)

neutrophil_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "neutrophil_2_GATA", ident.2 = "neutrophil_2_WT", print.bar = FALSE)

progenitor_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "progenitor_1_GATA", ident.2 = "progenitor_1_WT", print.bar = FALSE)
progenitor_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "progenitor_2_GATA", ident.2 = "progenitor_2_WT", print.bar = FALSE)
progenitor_3_gata_vs_wt <- FindMarkers(combined, ident.1 = "progenitor_3_GATA", ident.2 = "progenitor_3_WT", print.bar = FALSE)

hsc_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "HSC_1_GATA", ident.2 = "HSC_1_WT", print.bar = FALSE)
hsc_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "HSC_2_GATA", ident.2 = "HSC_2_WT", print.bar = FALSE)

erythrocytes_1_gata_vs_wt <- FindMarkers(combined, ident.1 = "erythrocytes_1_GATA", ident.2 = "erythrocytes_1_WT", print.bar = FALSE)
erythrocytes_2_gata_vs_wt <- FindMarkers(combined, ident.1 = "erythrocytes_2_GATA", ident.2 = "erythrocytes_2_WT", print.bar = FALSE)

compared_feature_list <- c("dap1b", "meis1b")
FeatureHeatmap(combined, features.plot = compared_feature_list, group.by = "genotype", pt.size = 2, key.position = "top", max.exp = 3)

# Heatmap of gene differences
combined <- SetAllIdent(combined, id = "celltype.genotype")
DoHeatmap(object = combined, genes.use = rownames(overall_gata_vs_wt))


subset <- SubsetData(combined, ident.use = c("progenitor_1_GATA", "progenitor_1_WT"), subset.raw = T)
subset <- SetAllIdent(subset, id = "genotype")
DoHeatmap(object = subset, genes.use = rownames(overall_gata_vs_wt))


#Scatter plots of group comparison
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