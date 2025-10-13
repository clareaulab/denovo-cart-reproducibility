library(Seurat)
library(dplyr)
library(data.table)
setwd("/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq")

## Load CD22 experiment matrix
cd22_rpmi_mat = Read10X_h5("./data/gex/20250923_CAR5p_CD22_RPMI8226_GEX_rna.h5")
cd22_rpmi_so = CreateSeuratObject(counts = cd22_rpmi_mat, assay = "RNA")
cd22_rpmi_barcodes = read.delim("./output/20250923_CAR5p_CD22_RPMI8226_hashing_assignment.tsv", header = FALSE, row.names = 1)
cd22_rpmi_so@meta.data$HTO = cd22_rpmi_barcodes[rownames(cd22_rpmi_so@meta.data), 1]
cd22_rpmi_so@meta.data$binder_name = cd22_rpmi_barcodes[rownames(cd22_rpmi_so@meta.data), 2]


# QC-metrics for RNA
DefaultAssay(cd22_rpmi_so) = "RNA"
cd22_rpmi_so[["percent.mt"]] = PercentageFeatureSet(cd22_rpmi_so, pattern = "^MT-") # Mito %
cd22_rpmi_so[["percent.ribo"]] = PercentageFeatureSet(cd22_rpmi_so, pattern = "^RPL|^RPS") # Ribosomal %
cd22_rpmi_so_pre_qc_plot = VlnPlot(
  cd22_rpmi_so,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo"),
  ncol = 4, group.by = "binder_name", pt.size = 0.1
)
cd22_rpmi_so_pre_qc_plot

# Apply Filtering
cd22_rpmi_so_filtered = subset(
  cd22_rpmi_so,
  subset = nFeature_RNA > 500 & 
    nCount_RNA > 500 &
    nFeature_RNA > 500 &
    percent.mt < 20 &
    (!is.na(binder_name))
)

## RNA counts are reallllly low 
VlnPlot(
  cd22_rpmi_so_filtered,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 4, group.by = "binder_name", pt.size = 0.1, log = TRUE
)

# seurat machinery goes brr anyways
cd22_rpmi_so_filtered = cd22_rpmi_so_filtered %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(assay="RNA",reduction_name="pca") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_rna")

cd22_rpmi_so_filtered = CellCycleScoring(cd22_rpmi_so_filtered, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

## Find clusters
cd22_rpmi_so_filtered =  cd22_rpmi_so_filtered %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.2)

## Inspect clusters
DimPlot(cd22_rpmi_so_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)

Idents(cd22_rpmi_so_filtered) = "binder_name"
Idents(cd22_rpmi_so_filtered) = "seurat_clusters"

FindMarkers(
  cd22_rpmi_so_filtered,
  ident.1 = "CD22_i61_2",
  ident.2 = "CD22_m971",
  #ident.1 = 2,
  #ident.2 = 0,
  logfc.threshold = 0.1,
  min.pct = 0.01,
  only.pos = TRUE
)

## Visualize markers
FeaturePlot(
  cd22_rpmi_so_filtered,
  features = c("CD3E","GZMB","GZMA","LTB","CCL5","TNFRSF18","IL2","TNF","MS4A1","MKI67"),
  reduction = "umap_rna",
  ncol = 5
)

VlnPlot(cd22_rpmi_so_filtered,features = c("GZMB","IL2", "IFNG", "CD69", "CD25"),group.by = "binder_name",log=TRUE)

saveRDS(cd22_rpmi_so_filtered,"./data/seurat_objects/cd22_rpmi_so_filtered.rds")

## Write the outputs for starCAT calculation
library(R.utils)
cd22_rpmi_filtered_counts = cd22_rpmi_so_filtered@assays$RNA$counts
data_dir = "/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq/data/starCAT/CD22_RPMI_input/"
writeMM(cd22_rpmi_filtered_counts, paste0(data_dir, 'matrix.mtx'))
gzip(paste0(data_dir, 'matrix.mtx'))

# Output cell barcodes
cd22_rpmi_barcodes = colnames(cd22_rpmi_filtered_counts)
write_delim(as.data.frame(cd22_rpmi_barcodes), paste0(data_dir, 'barcodes.tsv'),col_names = FALSE)
gzip(paste0(data_dir, 'barcodes.tsv'))

# Output feature names
cd22_rpmi_gene_names = rownames(cd22_rpmi_filtered_counts)
cd22_rpmi_features = data.frame("gene_id" = cd22_rpmi_gene_names,"gene_name" = cd22_rpmi_gene_names,type = "Gene Expression")
write_delim(as.data.frame(cd22_rpmi_features),delim = "\t", paste0(data_dir, 'features.tsv'),
            col_names = FALSE)
gzip(paste0(data_dir, 'features.tsv'))

