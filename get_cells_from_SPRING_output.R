data_location <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-2/project/DMSO_Flk+/DMSO_Flk+.counts.tr.csv"
cells_to_get <- read.table("/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-2/project/DMSO_Flk+/selected_cells (3).txt")
save_as <- "/Users/jon/data_analysis/scRNAseq inDrop second try with Elliott/20180201-Fast-2/project/DMSO_Flk+/selected_cell_output.csv"

num_of_columns <- read.table(data_location, sep = ",", nrows = 1)
num_of_columns <- ncol(num_of_columns)

cells_to_get <- unlist(cells_to_get)
cells_to_get <- cells_to_get + 2 #add 2 because first column is gene name and cells_to_get is 0-indexed, but R is 1-indexed

cc <- rep('NULL', num_of_columns)
cc_rownames <- rep('NULL', num_of_columns)
cc_rownames[1] <- 'character'

cc[cells_to_grab] <- 'integer'
output <- read.table(data_location, sep = ",", colClasses = cc)
output_rownames <- read.table(data_location, sep = ",", colClasses = cc_rownames)

rownames(output) <- unlist(output_rownames)

write.table(output, file = save_as,
            quote = FALSE, sep = ",", row.names = TRUE, col.names = FALSE)
