#import data, make rows = cells and columns = genes

#filter out cells with fewer than X (1000) UMIs

#row normalize the gene expression

#filter genes with mean expression < X (0.1) and fano factor < Y (3)
# mean_filter = mean(E,1) >= Ecutoff;
# var_filter = var(E,1) ./ (mean(E,1)+0.0001) >= Vcutoff;
# gene_filter = all([mean_filter; var_Filter],1);
# Efiltered = E(:, gene_filter);

#Z-score and do PCA with 20-25 PCs
# use scale() for z-score and prcomp() for PCA

# Do T-SNE
# tsne_out <- Rtsne(data_sample, dims = 2, theta = 0.5, perplexity = 10, pca = FALSE)
# 
# plot(tsne_out$Y, pch = 21, bg = color, col = NULL)