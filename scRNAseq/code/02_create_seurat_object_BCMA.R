library(Seurat)
library(dplyr)
library(data.table)
library(harmony)

## Load BCMA experiment matrix
bcma_mat = Read10X_h5("../data/gex/20250923_CAR5p_BCMA_GEX_rna.h5")
bcma_so = CreateSeuratObject(counts = bcma_mat, assay = "RNA")
## Load BCMA barcode assignments
bcma_barcodes = read.delim("../output/20250923_CAR5p_BCMA_hashing_assignment.tsv", header = FALSE, row.names = 1)
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

# ## Find clusters. Over-cluster so we can assign CD4/CD8
# bcma_so_tumor_filtered =  bcma_so_tumor_filtered %>%
#   FindNeighbors() %>%
#   FindClusters(resolution = 3)

## Find clusters. Over-cluster so we can assign CD4/CD8
bcma_so_tumor_filtered =  bcma_so_tumor_filtered %>%
  FindNeighbors() %>%
  FindClusters(resolution = 2)


## Inspect clusters
DimPlot(bcma_so_tumor_filtered,
        reduction = "umap_rna",
        #reduction = "umap_harmony",
        group.by = c("seurat_clusters","Phase","binder_name")
)

## Compute Average CD4/ CD8A expression per cluster
bcma_cd_agg_expression = AggregateExpression(bcma_so_tumor_filtered,features=c("CD4","CD8A","CD8B"))$RNA
bcma_cd_agg_expression = as.data.frame(t(bcma_cd_agg_expression))
bcma_cd_agg_expression$CD8 = bcma_cd_agg_expression$CD8A + bcma_cd_agg_expression$CD8B
bcma_cd_agg_expression$CD8CD4_Ratio = bcma_cd_agg_expression$CD8 / bcma_cd_agg_expression$CD4
bcma_cd_agg_expression[["seurat_cluster_name"]] = rownames(bcma_cd_agg_expression)
bcma_cd_agg_expression[["seurat_clusters"]] = as.numeric(substr(rownames(bcma_cd_agg_expression),start=2,stop=3))

cutoff_val = 1
bcma_cd8_cd4_ratio_cutoff_plot = ggplot(bcma_cd_agg_expression, aes(x = CD4, y = CD8, color = seurat_cluster_name)) + 
  geom_point(size = 4) +
  geom_abline(slope = cutoff_val, linetype = "dashed") +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_minimal()
bcma_cd8_cd4_ratio_cutoff_plot

# ## Test majotiy voting
# bcma_so_tumor_filtered_counts  <- GetAssayData(bcma_so_tumor_filtered, layer = "counts")
# bcma_so_tumor_filtered_scaled  <- GetAssayData(bcma_so_tumor_filtered, layer = "scale.data")
# 
# cd8_module <- c("CD8A","CD8B","GZMK","GZMB","NKG7","CCL5","CTSW")
# cd4_module <- c("CD4","IL7R","TCF7","MAL","LTB")
# bcma_so_tumor_filtered <- AddModuleScore(bcma_so_tumor_filtered, features = list(cd8_module), name = "CD8_mod")
# bcma_so_tumor_filtered <- AddModuleScore(bcma_so_tumor_filtered, features = list(cd4_module), name = "CD4_mod")
# 
# bcma_so_tumor_filtered_meta <- bcma_so_tumor_filtered@meta.data %>%
#   mutate(
#     cd8_det = (bcma_so_tumor_filtered_counts["CD8A", ] > 0) | (bcma_so_tumor_filtered_counts["CD8B", ] > 0),
#     cd4_det =  bcma_so_tumor_filtered_counts["CD4",  ] > 0
#   )
# 
# cluster_votes <- bcma_so_tumor_filtered_meta %>%
#   group_by(seurat_clusters) %>%
#   summarise(
#     vote_raw    = as.integer((mean(bcma_so_tumor_filtered_counts["CD8A", cur_group_rows()]) + mean(bcma_so_tumor_filtered_counts["CD8B", cur_group_rows()])) >= mean(bcma_so_tumor_filtered_counts["CD4", cur_group_rows()])),
#     vote_pct    = as.integer(mean(cd8_det) >= mean(cd4_det)),
#     vote_module = as.integer(mean(CD8_mod1) >= mean(CD4_mod1)),
#     cd8_votes   = vote_raw + vote_pct + vote_module,
#     lineage     = if_else(cd8_votes >= 2, "CD8", "CD4"),
#     .groups = "drop"
#   )
# 
# cluster_map <- setNames(cluster_votes$lineage, cluster_votes$seurat_clusters)
# bcma_so_tumor_filtered$lineage  <- cluster_map[as.character(bcma_so_tumor_filtered$seurat_clusters)]


## Assign cluster based on cd8/cd4 cutoffs and propogate label to all cells in cluster
meta_data_temp = bcma_so_tumor_filtered@meta.data %>% dplyr::select(-one_of(c("CD8CD4_Ratio","cd_type")))
meta_data_temp$cell = rownames(meta_data_temp)
meta_data_temp$seurat_cluster_name = paste0("g",meta_data_temp$seurat_clusters)
meta_data_temp$seurat_cluster_name = as.character(meta_data_temp$seurat_cluster_name)
merged_meta = merge(meta_data_temp,bcma_cd_agg_expression,by="seurat_cluster_name",all.x=TRUE)
rownames(merged_meta) = merged_meta$cell
merged_meta_ordered = merged_meta[rownames(bcma_so_tumor_filtered@meta.data),]
merged_meta_ordered$cd_type = case_when(
  merged_meta_ordered$CD8CD4_Ratio < cutoff_val ~ "CD4+",
  merged_meta_ordered$CD8CD4_Ratio >= cutoff_val ~ "CD8+"
)
# cd8_cutoff = 250
# merged_meta_ordered$cd_type = case_when(
#   merged_meta_ordered$CD8 >= cd8_cutoff ~ "CD8+",
#   merged_meta_ordered$CD8 < cd8_cutoff ~ "CD4+"
# )
merged_meta_ordered
## Add metadata
bcma_so_tumor_filtered = AddMetaData(bcma_so_tumor_filtered,metadata = merged_meta_ordered$CD8CD4_Ratio,col.name = 'CD8CD4_Ratio')
bcma_so_tumor_filtered = AddMetaData(bcma_so_tumor_filtered,metadata = merged_meta_ordered$cd_type,col.name = 'cd_type')

DimPlot(bcma_so_tumor_filtered,reduction="umap_rna",group.by="seurat_clusters")

cutoff_plot = DimPlot(bcma_so_tumor_filtered,reduction="umap_rna",group.by="cd_type")
cutoff_plot

FeaturePlot(bcma_so_tumor_filtered,features = c("CD8A","CD8B","CD4","GZMB")) | 
  cutoff_plot |
  DimPlot(bcma_so_tumor_filtered,group.by = "seurat_clusters") 
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

table(bcma_so_tumor_filtered@meta.data$cd_type,bcma_so_tumor_filtered@meta.data$binder_name)

FeaturePlot(bcma_so_tumor_filtered,features=c("CD4","CD8A")) | DimPlot(bcma_so_tumor_filtered,group.by="cd_type")

bcma_so_tumor_filtered@meta.data %>%
  dplyr::count(binder_name, cd_type) %>%
  group_by(binder_name) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  arrange(binder_name, cd_type)

FeaturePlot(bcma_so_tumor_filtered,features=c("GZMA","GZMB","IL2"),split.by = "binder_name")

## Save the filtered seurat object for further analysis
saveRDS(bcma_so_tumor_filtered,"../data/seurat_objects/bcma_so_tumor_filtered.rds")

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
