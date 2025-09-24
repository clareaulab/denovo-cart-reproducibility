library(data.table)
library(Signac)
library(Seurat)
library(Matrix)
library(dplyr)
library(BuenColors)

#sample_i = "20250923_CAR5p_BCMA"
#sample_i = "20250923_CAR5p_CD22_K562OE"
sample_i = "20250923_CAR5p_CD22_RPMI8226"

# import data
genes <- fread(paste0("../data/hashing/",sample_i,"_HTO_kboutput.genes.txt"), header = FALSE)[[1]]
barcodes <- fread(paste0("../data/hashing/",sample_i,"_HTO_kboutput.barcodes.txt"), header = FALSE)[[1]]
mtx <- fread(paste0("../data/hashing/",sample_i,"_HTO_kboutput.mtx"), skip = 3)
hto_mat <- sparseMatrix(
  i = c(mtx[[1]],length(barcodes)),
  j = c(mtx[[2]],length(genes)),
  x = c(mtx[[3]],0)
)
rownames(hto_mat) <- paste0(barcodes, "-1")
colnames(hto_mat) <- genes

data.frame(colSums(hto_mat))

# Import 
counts <- Read10X_h5(paste0("../data/gex/",sample_i,"_GEX_rna.h5"))
dim(counts)
common <- (intersect(colnames(counts), rownames(hto_mat)))
mat_ss <- hto_mat[common,]
hto_mat_ss <- t(mat_ss[,c(1:5)])

#Create Seurat object
so <- CreateSeuratObject(counts = counts[,common])

# Add hto
so[["HTO"]] <- CreateAssayObject(counts = hto_mat_ss[,colnames(so)] +1)

# analysis
so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR", margin =1)
so <- HTODemux(so, assay = "HTO", positive.quantile = 0.99)
HTOHeatmap(so, assay = "HTO", ncells = 5000)
table(so@meta.data$HTO_classification.global)
df <- so@meta.data[so@meta.data$hash.ID %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"),]

if(FALSE){
  #20250923_CAR5p_BCMA
  df$CAR <- case_when(
    df$hash.ID == "HTO1" ~ "BCMA_B4_I59_nonint_1399",
    df$hash.ID == "HTO2" ~ "BCMA_Abecma",
    df$hash.ID == "HTO3" ~ "BCMA_561726_WT",
    df$hash.ID == "HTO4" ~ "BCMA_B11_I59_int_1525",
    df$hash.ID == "HTO5" ~ "BCMA_A2_nonint_0366"
  )
  #"20250923_CAR5p_CD22_K562OE"
  df$CAR <- case_when(
    df$hash.ID == "HTO1" ~ "CD22_m971",
    df$hash.ID == "HTO2" ~ "CD22_i61_2",
    df$hash.ID == "HTO3" ~ "CD22_A3_nonint_4434",
    df$hash.ID == "HTO4" ~ "CD22_A11_int_1411",
    df$hash.ID == "HTO5" ~ "CD22_B1_nonint_2261"
  )
  #20250923_CAR5p_CD22_RPMI8226
  df$CAR <- case_when(
    df$hash.ID == "HTO1" ~ "CD22_m971",
    df$hash.ID == "HTO2" ~ "CD22_i61_2",
    df$hash.ID == "HTO3" ~ "CD22_B1_nonint_2261",
    df$hash.ID == "HTO4" ~ "CD22_A3_nonint_4434",
    df$hash.ID == "HTO5" ~ "CD22_A11_int_1411"
  )
}

write.table(df[,c("HTO_classification", "CAR")], file = paste0("../output/", sample_i, "_hashing_assignment.tsv"),
            row.names = TRUE, sep = "\t", quote = FALSE, col.names = FALSE)
