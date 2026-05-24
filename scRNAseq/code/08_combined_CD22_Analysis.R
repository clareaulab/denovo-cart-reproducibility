library(Seurat)
library(dplyr)
library(data.table)

## Load RDS
cd22_so_oe_tumor_filtered = readRDS("../data/seurat_objects/cd22_oe_so_tumor_filtered.rds")
cd22_rpmi_tumor_filtered = readRDS("../data/seurat_objects/cd22_rpmi_so_filtered.rds")

cd22_so_oe_tumor_filtered$experiment = "OE"
cd22_rpmi_tumor_filtered$experiment = "RPMI"

## Merge files
cd22_runs_combined = merge(
  cd22_so_oe_tumor_filtered,
  y = cd22_rpmi_tumor_filtered,
  add.cell.ids = c("OE", "RPMI"),
  project = "CD22Combined"
)
cd22_runs_combined = JoinLayers(cd22_runs_combined)
cd22_runs_combined

# seurat machinery goes brr
cd22_runs_combined = cd22_runs_combined %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(assay="RNA",reduction_name="pca") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_rna")

## Cellcycle scoring
cd22_runs_combined = CellCycleScoring(cd22_runs_combined, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# Run Harmony
cd22_runs_combined = cd22_runs_combined %>%
  RunHarmony("Phase") %>%
  RunUMAP(dims = 1:30, reduction ="pca", reduction.name = "umap_rna") %>%
  RunUMAP(dims = 1:30, reduction ="harmony", reduction.name = "umap_harmony")


DimPlot(cd22_runs_combined,group.by = c("experiment","binder_name"),reduction="umap_harmony")# | 
  #FeaturePlot(cd22_runs_combined,c("GZMB"),reduction="umap_harmony")


## Find clusters
cd22_runs_combined =  cd22_runs_combined %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.25)

## Inspect clusters
DimPlot(cd22_runs_combined,
        reduction = "umap_harmony",
        group.by = c("seurat_clusters","experiment","binder_name")
)

test_markers = FindAllMarkers(
  cd22_runs_combined,
  min.pct = 0.1,
  test.use = "wilcox"
)

test_markers %>% filter(cluster == 6)

