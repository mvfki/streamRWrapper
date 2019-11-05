This set of test data came from real experiment, which means, it has a large amount of single cells and realistic expression information. In fact, this set of cells mainly came from one cell type, but two types were put here for testing the functions. 
### Result Preview
- The 2D UMAP visualization
![image](http://github.com/mvfki/streamRWrapper/raw/master/testData_real/RStream_UMAP_2D.png)  
  
- SingleCellExperiment object converted from AnnData object
```
# After running the pipeline
sce <- adata2sce(adata)
sce
## class: SingleCellExperiment
## dim: 12288 2386
## metadata(0):
## assays(2): counts stream_matrix
## rownames(12288): Mrpl15 Lypla1 ... DHRSX CAAA01147332.1
## rowData names(4): geneSymbol gene_ids n_counts n_cells
## colnames(2386): AAACCCACAGTTCTAG-1 AAACGAACAAAGCGTG-1 ...
##   TTTGTTGCACATCATG-1 TTTGTTGGTATTGAGA-1
## colData names(5): barcodes label label_color n_counts n_genes
## reducedDimNames(4): STREAM_var_genes STREAM_X_mlle STREAM_X_dr
##   STREAM_X_vis_umap
## spikeNames(0):
```

