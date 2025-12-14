library(Seurat)
library(dplyr)
library(data.table)
cd22_runs_combined = merge(cd22_so_oe_tumor_filtered, y = cd22_rpmi_tumor_filtered, add.cell.ids = c("OE", "RPMI"), project = "CD22Combined")
cd22_runs_combined
