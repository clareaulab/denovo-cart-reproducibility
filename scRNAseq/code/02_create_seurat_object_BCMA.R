library(Seurat)
library(dplyr)
library(data.table)
setwd("/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq")

## Load BCMA experiment matrix
bcma_mat = Read10X_h5("./data/gex/20250923_CAR5p_BCMA_GEX_rna.h5")
bcma_so = CreateSeuratObject(counts = bcma_mat, assay = "RNA")
## Load BCMA barcode assignments
bcma_barcodes = read.delim("./output/20250923_CAR5p_BCMA_hashing_assignment.tsv", header = FALSE, row.names = 1)
bcma_so@meta.data$HTO = bcma_barcodes[rownames(bcma_so@meta.data), 1]
bcma_so@meta.data$binder_name = bcma_barcodes[rownames(bcma_so@meta.data), 2]

# QC-metrics for RNA
DefaultAssay(bcma_so) = "RNA"
bcma_so[["percent.mt"]] = PercentageFeatureSet(bcma_so, pattern = "^MT-") # Mito %
bcma_so[["percent.ribo"]] = PercentageFeatureSet(bcma_so, pattern = "^RPL|^RPS") # Ribosomal %
bcma_so_pre_qc_plot = VlnPlot(
  bcma_so,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo"),
  ncol = 4, group.by = "binder_name", pt.size = 0.1
)

# Apply Filtering
bcma_so_filtered = subset(
  bcma_so,
  subset = nFeature_RNA > 500 & 
    nCount_RNA < 30000 &
    nFeature_RNA < 5000 & 
    percent.mt < 20 &
    (!is.na(binder_name))
)

# seurat machinery goes brr
bcma_so_filtered = bcma_so_filtered %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(assay="RNA",reduction_name="pca") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_rna")

## Assign Cell cycle scores
bcma_so_filtered = CellCycleScoring(bcma_so_filtered, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

## Find clusters
bcma_so_filtered =  bcma_so_filtered %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.2)

## Inspect clusters
DimPlot(bcma_so_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
        )

## Identify cluster markers
bcma_all_markers = FindAllMarkers(bcma_so_filtered)

Idents(bcma_so_filtered) = "binder_name"
FindMarkers(
  bcma_so_filtered,
  ident.1 = "BCMA_B4_I59_nonint_1399",
  ident.2 = "BCMA_A2_nonint_0366",
  #ident.2 = "BCMA_Abecma",
  logfc.threshold = 0.1,
  min.pct = 0.01,
  only.pos = TRUE
)

## Visualize markers
FeaturePlot(
  bcma_so_filtered,
  features = c("CD3E","GZMB","GZMA","LTB","CCL5","TNFRSF18","IL2","TNF","MS4A1","PRAME"),
  reduction = "umap_rna",
  ncol = 5
)

## Inspection reveals cluster 6 and 7 are likely contaminated tumor cells
## We will remove them and re-perform clustering
bcma_so_tumor_filtered = subset(
  bcma_so_filtered,
  subset = ((seurat_clusters != 6) & (seurat_clusters != 7))
)

# Re-perform seurat machinery with now less cells
bcma_so_tumor_filtered = bcma_so_tumor_filtered %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(assay="RNA",reduction_name="pca") %>%
  RunHarmony("Phase") %>%
  RunUMAP(dims = 1:30, reduction ="pca", reduction.name = "umap_rna") %>%
  RunUMAP(dims = 1:30, reduction ="harmony", reduction.name = "umap_harmony")

bcma_so_tumor_filtered = CellCycleScoring(
  bcma_so_tumor_filtered, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes
)

## Find clusters
bcma_so_tumor_filtered =  bcma_so_tumor_filtered %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.2)

## Inspect clusters
DimPlot(bcma_so_tumor_filtered,
        #reduction = "umap_rna",
        reduction = "umap_harmony",
        group.by = c("seurat_clusters","Phase","binder_name")
)

#Idents(bcma_so_tumor_filtered) = "binder_name"
Idents(bcma_so_tumor_filtered) = "seurat_clusters"

FindMarkers(
  bcma_so_tumor_filtered,
  #ident.1 = "BCMA_B11_I59_int_1525",
  #ident.1 = "BCMA_561726_WT",
  ident.1= 0,
  ident.2 = 6,
  #ident.2 = "BCMA_561726_WT",
  #ident.2 = "BCMA_Abecma",
  logfc.threshold = 0.1,
  min.pct = 0.01,
  only.pos = TRUE
)

FeaturePlot(bcma_so_tumor_filtered,features=c("GZMA","GZMB","IL2"),split.by = "binder_name")

## Save the filtered seurat object for further analysis
saveRDS(bcma_so_tumor_filtered,"./data/seurat_objects/bcma_so_tumor_filtered.rds")

## Save the variable genes
bcma_so_variable_genes = VariableFeatures(bcma_so_tumor_filtered)
writeLines(bcma_so_variable_genes, "./data/seurat_objects/bcma_so_tumor_filtered_variable_genes.txt")

## Write the outputs for starCAT calculation
bcma_filtered_counts = bcma_so_tumor_filtered@assays$RNA$counts
data_dir = "/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq/data/starCAT/BCMA_input/"
writeMM(bcma_filtered_counts, paste0(data_dir, 'matrix.mtx'))
library(R.utils)
gzip(paste0(data_dir, 'matrix.mtx'))

# Output cell barcodes
bcma_barcodes = colnames(bcma_filtered_counts)
write_delim(as.data.frame(bcma_barcodes), paste0(data_dir, 'barcodes.tsv'),col_names = FALSE)
gzip(paste0(data_dir, 'barcodes.tsv'))

# Output feature names
bcma_gene_names = rownames(bcma_filtered_counts)
bcma_features = data.frame("gene_id" = bcma_gene_names,"gene_name" = bcma_gene_names,type = "Gene Expression")
write_delim(as.data.frame(bcma_features),delim = "\t", paste0(data_dir, 'features.tsv'),
            col_names = FALSE)
gzip(paste0(data_dir, 'features.tsv'))

## Command ran
## starcat -r "TCAT.V1" -c /home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq/data/starCAT/BCMA_input/matrix.mtx.gz --output-dir /home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq/data/starCAT/BCMA_output --name "bcma"
