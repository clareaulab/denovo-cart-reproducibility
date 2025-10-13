library(Seurat)
library(dplyr)
library(data.table)
setwd("/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq")

## Load CD22 experiment matrix
cd22_oe_mat = Read10X_h5("./data/gex/20250923_CAR5p_CD22_K562OE_GEX_rna.h5")
cd22_oe_so = CreateSeuratObject(counts = cd22_oe_mat, assay = "RNA")
cd22_oe_barcodes = read.delim("./output/20250923_CAR5p_CD22_K562OE_hashing_assignment.tsv", header = FALSE, row.names = 1)
cd22_oe_so@meta.data$HTO = cd22_oe_barcodes[rownames(cd22_oe_so@meta.data), 1]
cd22_oe_so@meta.data$binder_name = cd22_oe_barcodes[rownames(cd22_oe_so@meta.data), 2]

# QC-metrics for RNA
DefaultAssay(cd22_oe_so) = "RNA"
cd22_oe_so[["percent.mt"]] = PercentageFeatureSet(cd22_oe_so, pattern = "^MT-") # Mito %
cd22_oe_so[["percent.ribo"]] = PercentageFeatureSet(cd22_oe_so, pattern = "^RPL|^RPS") # Ribosomal %
cd22_oe_so_pre_qc_plot = VlnPlot(
  cd22_oe_so,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo"),
  ncol = 4, group.by = "binder_name", pt.size = 0.1
)
cd22_oe_so_pre_qc_plot

# Apply Filtering
cd22_oe_so_filtered = subset(
  cd22_oe_so,
  subset = nFeature_RNA > 500 & 
    nCount_RNA < 30000 &
    nFeature_RNA > 500 & 
    nFeature_RNA < 5000 & 
    percent.mt < 20 &
    (!is.na(binder_name))
)


# seurat machinery goes brr
cd22_oe_so_filtered = cd22_oe_so_filtered %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(assay="RNA",reduction_name="pca") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_rna")

cd22_oe_so_filtered = CellCycleScoring(cd22_oe_so_filtered, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

## Find clusters
cd22_oe_so_filtered =  cd22_oe_so_filtered %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.2)

## Inspect clusters
DimPlot(cd22_oe_so_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)

## Identify cluster markers
#bcma_all_markers = FindAllMarkers(bcma_so_filtered)

Idents(cd22_oe_so_filtered) = "binder_name"
Idents(cd22_oe_so_filtered) = "seurat_clusters"

# FindMarkers(
#   cd22_oe_so_filtered,
#   ident.1 = "CD22_i61_2",
#   ident.2 = "CD22_m971",
#   logfc.threshold = 0.1,
#   min.pct = 0.01,
#   only.pos = TRUE
# )

## Visualize markers
FeaturePlot(
  cd22_oe_so_filtered,
  features = c("CD3E","GZMB","GZMA","LTB","CCL5","TNFRSF18","IL2","TNF","MS4A1","PRAME"),
  reduction = "umap_rna",
  ncol = 5
)

## Inspection reveals cluster 5 is likely contaminated tumor cells
## We will remove them and re-perform clustering
cd22_oe_so_tumor_filtered = subset(
  cd22_oe_so_filtered,
  subset = (seurat_clusters != 5)
)

# Reinspect the clusters, we need to regress out the cell-cycle effects
DimPlot(cd22_oe_so_tumor_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)

bcma_variable_genes = readLines("./data/seurat_objects/bcma_so_tumor_filtered_variable_genes.txt")

# Re-perform seurat machinery with now less cells
cd22_oe_so_tumor_filtered = cd22_oe_so_tumor_filtered %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(assay="RNA",reduction_name="pca", features = bcma_variable_genes) %>%
  RunHarmony("Phase") %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap_rna") %>%
  RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony")

cd22_oe_so_tumor_filtered = CellCycleScoring(
  cd22_oe_so_tumor_filtered, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes
)

## Find clusters
cd22_oe_so_tumor_filtered =  cd22_oe_so_tumor_filtered %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.2)

## Inspect clusters
DimPlot(cd22_oe_so_tumor_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)

Idents(cd22_oe_so_tumor_filtered) = "binder_name"
#Idents(bcma_so_tumor_filtered) = "seurat_clusters"

FindMarkers(
  cd22_oe_so_tumor_filtered,
  ident.1 = "CD22_A11_int_1411",
  ident.2 = "CD22_m971",
  logfc.threshold = 0.1,
  min.pct = 0.01,
  only.pos = TRUE
)

## Visualize markers
FeaturePlot(
  cd22_oe_so_tumor_filtered,
  features = c("CD3E","GZMB","GZMA","LTB","CCL5","TNFRSF18","IL2","TNF","MS4A1","PRAME"),
  reduction = "umap_rna",
  ncol = 5
)

saveRDS(cd22_oe_so_tumor_filtered,"./data/seurat_objects/cd22_oe_so_tumor_filtered.rds")

## Write the outputs for starCAT calculation
library(R.utils)
cd22_oe_filtered_counts = cd22_oe_so_tumor_filtered@assays$RNA$counts
data_dir = "/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq/data/starCAT/CD22_OE_input/"
writeMM(cd22_oe_filtered_counts, paste0(data_dir, 'matrix.mtx'))
gzip(paste0(data_dir, 'matrix.mtx'))

# Output cell barcodes
cd22_oe_barcodes = colnames(cd22_oe_filtered_counts)
write_delim(as.data.frame(cd22_oe_barcodes), paste0(data_dir, 'barcodes.tsv'),col_names = FALSE)
gzip(paste0(data_dir, 'barcodes.tsv'))

# Output feature names
cd22_oe_gene_names = rownames(cd22_oe_filtered_counts)
cd22_oe_features = data.frame("gene_id" = cd22_oe_gene_names,"gene_name" = cd22_oe_gene_names,type = "Gene Expression")
write_delim(as.data.frame(cd22_oe_features),delim = "\t", paste0(data_dir, 'features.tsv'),
            col_names = FALSE)
gzip(paste0(data_dir, 'features.tsv'))

