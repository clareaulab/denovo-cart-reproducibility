library(Seurat)
library(dplyr)
library(data.table)
library(BuenColors)
setwd("/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq")

## Load data
bcma_so_tumor_filtered = readRDS("./data/seurat_objects/bcma_so_tumor_filtered.rds")

selected_markers_bcma = VlnPlot(bcma_so_tumor_filtered, features = c("IL2","IFNG","GZMA","GZMB"),ncol=4)
cowplot::ggsave2("./plots/BCMA/selected_markers.png",selected_markers_bcma,dpi=300,width=12,height=4)

## Inspect clusters
bcma_dim_plot = DimPlot(bcma_so_tumor_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)
bcma_dim_plot
cowplot::ggsave2("./plots/BCMA/Dimplot.png",bcma_dim_plot,dpi=300,width=12,height=4)

bcma_so_tumor_filtered@meta.data$binder_name = factor(bcma_so_tumor_filtered@meta.data$binder_name,levels = c(
  "BCMA_Abecma","BCMA_561726_WT","BCMA_B11_I59_int_1525","BCMA_A2_nonint_0366","BCMA_B4_I59_nonint_1399"
))
## Draw the dots properly
bcma_dots_simple = DimPlot(bcma_so_tumor_filtered, label = FALSE, group.by = c( "binder_name"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  scale_color_manual(values = c("#FFB81C","#D91E18","#8B0000","#9966CC","#00BFFF")) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")
bcma_dots_simple
cowplot::ggsave2(bcma_dots_simple, file = "./plots/BCMA/umap_base.png", width = 6, height = 6, dpi = 300)

bcma_dots_simple

## Draw the features
mk_plot <- function(so,gene){
  pu <- FeaturePlot(so, features = c(gene),  
                    pt.size = 0.1, max.cutoff = "q90") + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}
bcma_feature_plots = cowplot::plot_grid(
  mk_plot(bcma_so_tumor_filtered,"CD3E"),
  mk_plot(bcma_so_tumor_filtered,"IFNG"),
  mk_plot(bcma_so_tumor_filtered,"TNF"),
  ncol=2,scale=1
)
cowplot::ggsave2(bcma_feature_plots,file="./plots/BCMA/umap_features.png",width=6*2,height=6*2, dpi=600)

## Split activated from non-activated
Idents(bcma_so_tumor_filtered)="binder_name"
bcma_so_tumor_filtered@meta.data$binder_is_activated = bcma_so_tumor_filtered@meta.data$binder_name %in% c("BCMA_561726_WT","BCMA_B11_I59_int_1525","BCMA_Abecma")

split_activate_plot = DimPlot(bcma_so_tumor_filtered,reduction="umap_harmony",split.by = "binder_is_activated")

split_activate_plot

bcma_marker_data = FetchData(bcma_so_tumor_filtered,
                             vars = c("IL2","IFNG","TNF","GZMA","GZMB","binder_name"),
                             slot = "data"
                             ) %>%
  pivot_longer(
    cols = c("IL2", "IFNG", "GZMA", "GZMB"),
    names_to = "gene",
    values_to = "val"
  )

bcma_markers_boxplot = ggplot(bcma_marker_data, aes(x = binder_name, y = val)) +
  geom_violin(aes(fill=binder_name)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) +
  pretty_plot(fontsize = 8) +
  L_border() + theme(legend.position = "none") +
  facet_wrap(~gene)
bcma_markers_boxplot

FeaturePlot(bcma_so_tumor_filtered,features=c("IL2","IFNG","TNF"),ncol=3, reduction = "umap_harmony")

split_activate_plot | FeaturePlot(bcma_so_tumor_filtered,features=c("IL2","IFNG"), reduction = "umap_harmony")
  
FeaturePlot(
  bcma_so_tumor_filtered,
  features = colnames(bcma_so_starCAT)
)
Idents(bcma_so_tumor_filtered) = "binder_name"

## Find Markers
bcma_so_tumor_filtered_clust_markers = FindAllMarkers(bcma_so_tumor_filtered,assay="RNA",group.by="seurat_clusters")
bcma_so_tumor_filtered_binder_markers = FindAllMarkers(bcma_so_tumor_filtered,assay="RNA",group.by="binder_name")

## Show markers by cluster
top_bcma_markers_by_clust = bcma_so_tumor_filtered_clust_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

Idents(bcma_so_tumor_filtered) = "seurat_clusters"
top_bcma_markers_by_clust_heatmap = DoHeatmap(
  bcma_so_tumor_filtered,
  features = top_bcma_markers_by_clust$gene
) + NoLegend()
top_bcma_markers_by_clust_heatmap
ggsave("./plots/top_bcma_markers_by_rna_clusters.png",top_bcma_markers_by_clust_heatmap,dpi=300,height=8,width=12)

## Show markers by binder names
top_bcma_markers_by_binders = bcma_so_tumor_filtered_binder_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

Idents(bcma_so_tumor_filtered) = "binder_name"
top_bcma_markers_by_binder_heatmap = DoHeatmap(
  bcma_so_tumor_filtered,
  features = top_bcma_markers_by_binders$gene
) + NoLegend()
top_bcma_markers_by_binder_heatmap
ggsave("./plots/top_bcma_markers_by_binders.png",top_bcma_markers_by_binder_heatmap,dpi=300,height=8,width=12)

## Runs fgsea enrichment analysis
run_fgsea_hallmark = function(markers_df, species = "Homo sapiens", top_filter = NULL) {
  hallmark_sets = msigdbr(species = species, category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  pathways = split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
  ranks = markers_df$avg_log2FC
  names(ranks) = markers_df$gene
  ranks = sort(ranks, decreasing = TRUE)
  fgsea_res = fgsea(pathways = pathways, stats = ranks, nperm = 1000) %>%
    arrange(padj)
  if (!is.null(top_filter)) {
    fgsea_res <- fgsea_res %>% head(top_filter)
  }
  return(fgsea_res)
}

run_fgsea_GO = function(markers_df, species = "Homo sapiens", top_filter = NULL) {
  go_sets = msigdbr(species = species, category = "C5", subcategory = "GO:BP") %>%
    dplyr::select(gs_name, gene_symbol)
  pathways = split(go_sets$gene_symbol, go_sets$gs_name)
  ranks = markers_df$avg_log2FC
  names(ranks) = markers_df$gene
  ranks = sort(ranks, decreasing = TRUE)
  fgsea_res = fgsea(pathways = pathways, stats = ranks, nperm = 1000) %>%
    dplyr::arrange(padj)
  if (!is.null(top_filter)) {
    fgsea_res <- fgsea_res %>% head(top_filter)
  }
  return(fgsea_res)
}

## Helper function to plot fgsea enrichment result
plot_fgsea_bar = function(fgsea_res, top_n = 20) {
  df = fgsea_res %>%
    mutate(
      sig_cat = case_when(
        padj < 0.05 ~ "padj < 0.05",
        pval < 0.05 ~ "pval < 0.05",
        TRUE ~ "N.S."
      )
    ) %>%
    arrange(desc(NES)) %>%
    head(top_n) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  p <- ggplot(df, aes(x = NES, y = pathway, fill = sig_cat)) +
    geom_col() +
    scale_fill_manual(
      values = c("padj < 0.05" = "red", "pval < 0.05" = "orange", "N.S." = "gray70")
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway",
      fill = "Significance"
      #title = "FGSEA Hallmark Enrichment"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

## Function that performs differential expression analysis
run_de_and_volcano = function(seurat_obj,
                              sample_col = "sample_name",
                              sample1,
                              sample2,
                              subset_col = NULL,
                              subset_value = NULL,
                              assay = "RNA",
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                              title_prefix = "Markers") {
  
  # ---- Conditional subsetting ----
  if (!is.null(subset_col) && !is.null(subset_value)) {
    subset_obj <- subset(seurat_obj, subset = !!as.name(subset_col) == subset_value)
    plot_title_suffix <- paste0(": ", subset_value)
  } else {
    subset_obj <- seurat_obj
    plot_title_suffix <- ""
  }
  print("setting idents")
  # set identities
  Idents(subset_obj) <- subset_obj[[sample_col, drop = TRUE]]
  # ---- Run differential expression ----
  print("Running differential")
  markers <- FindMarkers(
    subset_obj,
    ident.1 = sample1,
    ident.2 = sample2,
    assay = assay,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = "wilcox"
  )
  print("markers finfihsed")
  print(head(markers))
  # ---- Map features to gene names for ATAC/chromVAR assays ----
  if (assay %in% c("ATAC", "chromvar")) {
    print("handling motifs")
    # Get the named list of motif IDs and gene symbols from the ATAC assay's motifs slot
    motif_list <- seurat_obj[["ATAC"]]@motifs@motif.names
    # Create a new column with human-readable labels using a simple vector lookup
    markers$gene <- unlist(motif_list[rownames(markers)])
    #print(markers)
    # Use motif name if the gene symbol is missing
    markers$gene[is.na(markers$gene)] <- rownames(markers)[is.na(markers$gene)]
    print("motig name used")
  } else {
    # For RNA and other assays, the gene column is already the label
    markers$gene <- rownames(markers)
  }
  
  ## Guard against empty DEG table
  if (nrow(markers) == 0) {
    message(
      paste0(
        "No markers passed the criteria (logfc.threshold = ", logfc.threshold,
        ") for ", subset_value, " in ", assay, " assay. Returning an empty plot."
      )
    )
    # Create an empty ggplot object to serve as a placeholder plot.
    empty_plot <- ggplot() + 
      ggtitle(paste0(title_prefix, plot_title_suffix, " (", sample1, " vs ", sample2, ")")) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12)
      )
    
    return(list(markers = markers, volcano_plot = empty_plot))
  }
  
  # ---- Volcano plot ----
  print("PLotting volcanoes")
  print("marker before volcano")
  #print(markers)
  p <- EnhancedVolcano(
    markers,
    lab = markers$gene,
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 0.25,
    title = paste0(title_prefix, plot_title_suffix, " (", sample1, " vs ", sample2, ")"),
    subtitle = NULL,
    legendPosition = "right",
    labSize = 3.0,
    pointSize = 2.0,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'adjusted p-value')
  )
  
  return(list(markers = markers, volcano_plot = p))
}

## For a single data modality, perform DEG analysis and FGSEA on markers 
compare_samples_with_enrichment = function(seurat_obj,
                                           sample_col = "sample_name",
                                           sample1,
                                           sample2,
                                           subset_col = "predicted.celltype.l2",
                                           subset_value,
                                           assay = "RNA",
                                           min.pct = 0.1,
                                           logfc.threshold = 0,
                                           title_prefix = "Markers",
                                           top_fgsea = 50) {
  # run DE and volcano
  print("running run_de_and_volcano")
  de_out <- run_de_and_volcano(
    seurat_obj = seurat_obj,
    sample_col = sample_col,
    sample1 = sample1,
    sample2 = sample2,
    subset_col = subset_col,
    subset_value = subset_value,
    assay = assay,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    title_prefix = title_prefix
  )
  
  # run fgsea
  print("Running fgsea")
  fgsea_res <- run_fgsea_hallmark(de_out$markers)
  fgsea_res_go <- run_fgsea_GO(de_out$markers)
  print("plottig fgsea")
  # plot fgsea
  fgsea_plot <- plot_fgsea_bar(fgsea_res, top_n = top_fgsea)
  fgsea_plot_go <- plot_fgsea_bar(fgsea_res_go, top_n = top_fgsea)
  print("Done plotting fgsea")
  return(list(
    markers = de_out$markers,
    volcano_plot = de_out$volcano_plot,
    fgsea = fgsea_res,
    fgsea_go = fgsea_res_go,
    fgsea_plot = fgsea_plot,
    fgsea_plot_GO = fgsea_plot_go
  ))
}

## Perform DEG analysis between all de novo binders to ABECMA
de_novo_binders = c("BCMA_561726_WT","BCMA_B4_I59_nonint_1399","BCMA_B11_I59_int_1525","BCMA_A2_nonint_0366")
for (binder in de_novo_binders) {
  denovo_abecma_comparison = compare_samples_with_enrichment(
    bcma_so_tumor_filtered,
    sample_col="binder_name",
    sample1=binder,
    sample2="BCMA_Abecma",
    subset_col=NULL,
    subset_value=NULL
  )
  cowplot::ggsave2(paste0("./plots/BCMA/volcano_",binder,"_vs_BCMA_Abecma.pdf"),denovo_abecma_comparison$volcano_plot,dpi=300,width=8,height=4)
  cowplot::ggsave2(paste0("./plots/BCMA/volcano_",binder,"_vs_BCMA_Abecma.png"),denovo_abecma_comparison$volcano_plot,dpi=300,width=8,height=4)
  cowplot::ggsave2(paste0("./plots/BCMA/fgsea_",binder,"_vs_BCMA_Abecma.pdf"),denovo_abecma_comparison$fgsea_plot,dpi=300,width=8)
  cowplot::ggsave2(paste0("./plots/BCMA/fgsea_go_",binder,"_vs_BCMA_Abecma.pdf"),denovo_abecma_comparison$fgsea_plot_GO,dpi=300,width=24,height=24)
  write.csv(denovo_abecma_comparison$markers,paste0("./data/degs/degs_",binder,"_vs_BCMA_Abecma.csv"))
}

## Perform DEG analysis between all evolved binders to parental
evolved_binders = c("BCMA_B4_I59_nonint_1399","BCMA_B11_I59_int_1525","BCMA_A2_nonint_0366")
#evolved_binders = c("BCMA_B11_I59_int_1525")
for (binder in evolved_binders) {
  evolved_parental_comparison = compare_samples_with_enrichment(
    bcma_so_tumor_filtered,
    sample_col="binder_name",
    sample1=binder,
    sample2="BCMA_561726_WT",
    subset_col=NULL,
    subset_value=NULL
  )
  cowplot::ggsave2(paste0("./plots/BCMA/volcano_",binder,"_vs_BCMA_561726_WT.pdf"),evolved_parental_comparison$volcano_plot,dpi=300,width=8,height=4)
  cowplot::ggsave2(paste0("./plots/BCMA/volcano_",binder,"_vs_BCMA_561726_WT.png"),evolved_parental_comparison$volcano_plot,dpi=300,width=8,height=4)
  cowplot::ggsave2(paste0("./plots/BCMA/fgsea_",binder,"_vs_BCMA_561726_WT.pdf"),evolved_parental_comparison$fgsea_plot,dpi=300,width=8,height=8)
  write.csv(evolved_parental_comparison$markers,paste0("./data/degs/degs_",binder,"_vs_BCMA_561726_WT.csv"))
}

## Compared good evolved to no activation unevolved
evolved_good_vs_evolved_bad_comparison = compare_samples_with_enrichment(
  bcma_so_tumor_filtered,
  sample_col="binder_name",
  sample1="BCMA_B11_I59_int_1525",
  sample2="BCMA_B4_I59_nonint_1399",
  subset_col=NULL,
  subset_value=NULL
)
cowplot::ggsave2(paste0("./plots/BCMA/volcano_","BCMA_B11_I59_int_1525","_vs_BCMA_B4_I59_nonint_1399.pdf"),evolved_good_vs_evolved_bad_comparison$volcano_plot,dpi=300,width=8,height=4)
cowplot::ggsave2(paste0("./plots/BCMA/volcano_","BCMA_B11_I59_int_1525","_vs_BCMA_B4_I59_nonint_1399.png"),evolved_good_vs_evolved_bad_comparison$volcano_plot,dpi=300,width=8,height=4)
cowplot::ggsave2(paste0("./plots/BCMA/fgsea_","BCMA_B11_I59_int_1525","_vs_BCMA_B4_I59_nonint_1399.pdf"),evolved_good_vs_evolved_bad_comparison$fgsea_plot,dpi=300,width=8,height=8)
write.csv(evolved_good_vs_evolved_bad_comparison$markers,paste0("./data/degs/degs_","BCMA_B11_I59_int_1525","_vs_BCMA_B4_I59_nonint_1399.csv"))

## Compare the other one
evolved_good_vs_evolved_bad_comparison_2 = compare_samples_with_enrichment(
  bcma_so_tumor_filtered,
  sample_col="binder_name",
  sample1="BCMA_B11_I59_int_1525",
  sample2="BCMA_A2_nonint_0366",
  subset_col=NULL,
  subset_value=NULL
)
cowplot::ggsave2(paste0("./plots/BCMA/volcano_","BCMA_B11_I59_int_1525","_vs_BCMA_A2_nonint_0366.pdf"),evolved_good_vs_evolved_bad_comparison_2$volcano_plot,dpi=300,width=8,height=4)
cowplot::ggsave2(paste0("./plots/BCMA/volcano_","BCMA_B11_I59_int_1525","_vs_BCMA_A2_nonint_0366.png"),evolved_good_vs_evolved_bad_comparison_2$volcano_plot,dpi=300,width=8,height=4)
cowplot::ggsave2(paste0("./plots/BCMA/fgsea_","BCMA_B11_I59_int_1525","_vs_BCMA_A2_nonint_0366.pdf"),evolved_good_vs_evolved_bad_comparison_2$fgsea_plot,dpi=300,width=8,height=8)
write.csv(evolved_good_vs_evolved_bad_comparison_2$markers,paste0("./data/degs/degs_","BCMA_B11_I59_int_1525","_vs_BCMA_A2_nonint_0366.csv"))

## Also compare difference between tonic vs. no activation
tonic_vs_nonactive_comparison = compare_samples_with_enrichment(
  bcma_so_tumor_filtered,
  sample_col="binder_name",
  sample1="BCMA_B4_I59_nonint_1399",
  sample2="BCMA_A2_nonint_0366",
  subset_col=NULL,
  subset_value=NULL
)
cowplot::ggsave2(paste0("./plots/BCMA/volcano_","BCMA_B4_I59_nonint_1399","_vs_BCMA_A2_nonint_0366.pdf"),tonic_vs_nonactive_comparison$volcano_plot,dpi=300,width=8,height=4)
cowplot::ggsave2(paste0("./plots/BCMA/volcano_","BCMA_B4_I59_nonint_1399","_vs_BCMA_A2_nonint_0366.png"),tonic_vs_nonactive_comparison$volcano_plot,dpi=300,width=8,height=4)
cowplot::ggsave2(paste0("./plots/BCMA/fgsea_","BCMA_B4_I59_nonint_1399","_vs_BCMA_A2_nonint_0366.pdf"),tonic_vs_nonactive_comparison$fgsea_plot,dpi=300,width=8,height=8)
write.csv(tonic_vs_nonactive_comparison$markers,paste0("./data/degs/degs_","BCMA_B4_I59_nonint_1399","_vs_BCMA_A2_nonint_0366.csv"))


## Load the overexpressed genes from evolved to wildtype as binder activation module score
evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_561726_WT.csv")
tonic_vs_nonactive = read.csv("./data/degs/degs_BCMA_B4_I59_nonint_1399_vs_BCMA_A2_nonint_0366.csv")
#evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_B4_I59_nonint_1399.csv")
#evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_A2_nonint_0366.csv")

# tcat_ref = t(read.table("./data/starCAT/TCAT.V1.reference.tsv")) %>% as.data.frame()
# binder_activation_genes = tcat_ref %>% filter(Cytotoxic > 0) %>% arrange(-Cytotoxic) %>% rownames()
# binder_activation_genes = binder_activation_genes#[1:500]

binder_activation_genes = evolved_vs_parental %>%
  filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>%
  pull(gene)
binder_activation_genes
binder_tonic_genes = tonic_vs_nonactive %>%
  filter(avg_log2FC > 0.75, p_val_adj < 0.05) %>%
  pull(gene)
binder_tonic_genes

## Calculate Module Score
effector_genes = c(
  "TNFRSF9","CTLA4","TIGIT","PRDM1",
  "CD40LG","KLRG1","TBX21","ZBTB32"
)

exhaustion_genes = c(
  "BATF3","NR4A2","PDCD1","BHLHE40",
  "HAVCR2","TOX","NR4A3","ENTPD1","LAG3"
)
function_genes = c(
  "CSF1","CCL20","IL2","TNF","PRF1",
  "IFNG","GZMB","CCL1","CCL3","CCL4"
)
memory_genes = c(
  "TCF7","LEF1","IL7R","CXCR3","KLF2","CD27"
)

general_activation_genes = c(
  "CD69","ICOS","IL2","IFNG","TNFRSF9",
  "TNF","PRF1","GZMB","FASLG"
)

# VlnPlot(
#   bcma_so_tumor_filtered,
#   #features=c("TNFRSF9","CD28")
#   features=c("CD69","ICOS","IL2","IFNG","TNFRSF9",
#              "TNF","PRF1","GZMB","FASLG")
# )


module_gene_list = list(
  "effector" = effector_genes,
  "exhaustion" = exhaustion_genes,
  "function" = function_genes,
  "memory" = memory_genes,
  "activation" = general_activation_genes,
  "binder_activation" = binder_activation_genes,
  "binder_tonic" = binder_tonic_genes
)

bcma_so_tumor_filtered = AddModuleScore(
  object = bcma_so_tumor_filtered,
  features = module_gene_list,
  name = c("effector","exhaustion","function","memory", "activation","binder_activation","binder_tonic")
)

VlnPlot(bcma_so_tumor_filtered, pt.size = 0, 
        features = c("effector1","exhaustion2","function3","binder_activation6","binder_tonic7"),
        group.by = "binder_name") & geom_boxplot(outlier.shape = NA,)

## Visualize the DEG results
highlight_volcano_genes = c("IL2","IFNG","ENTPD1","SELL")
rownames(evolved_vs_parental) = evolved_vs_parental$gene
evolved_vs_parental[highlight_volcano_genes,]

evolved_vs_parental_capped = evolved_vs_parental %>% mutate(p_val_adj_capped = pmax(p_val_adj,10**-150))
bcma_volcano_highlight = evolved_vs_parental_capped %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj_capped), color = gene %in% highlight_volcano_genes)) + 
  geom_point(size = 0.5) +
  pretty_plot(fontsize = 7) + L_border() + 
  labs(x = "log2FC (B5.I0/B5)", y = "-log10(padj)") +
  scale_color_manual(values = c("black", "firebrick")) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 150)) +
  theme(legend.position = "none")

bcma_volcano_highlight
cowplot::ggsave2(bcma_volcano_highlight, file = "./plots/BCMA/volcano_bcma_evolved_vs_wildtype.pdf", width = 1.5, height = 1.5)



library(BuenColors)
# replacement_dict = c(
#   "BCMA_Abecma" = "Abecma", "BCMA_561726_WT" = "B5", "BCMA_B11_I59_int_1525" = "B5.I0",
#   "BCMA_A2_nonint_0366" = "", "BCMA_B4_I59_nonint_1399" = ""
#   )
bcma_so_tumor_filtered@meta.data$binder_name = factor(bcma_so_tumor_filtered@meta.data$binder_name,levels = c(
  "BCMA_Abecma","BCMA_561726_WT","BCMA_B11_I59_int_1525","BCMA_A2_nonint_0366","BCMA_B4_I59_nonint_1399"
))

bcma_module_score_boxplot = ggplot(bcma_so_tumor_filtered@meta.data,
                                   aes(x = binder_name, y = binder_activation6, fill=binder_name)) + 
  geom_violin(aes(fill=binder_name)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6) + 
  pretty_plot(fontsize = 8) + 
  L_border() + theme(legend.position = "none") +
  stat_compare_means(
      comparisons = list(c("BCMA_B11_I59_int_1525","BCMA_Abecma"))
  ) #+
  #scale_fill_manual(values = c("orange", "black", "dodgerblue3", "firebrick", "firebrick"))
bcma_module_score_boxplot

library(ggpubr)
#c("function3","activation5","Cytotox1")
p = VlnPlot(bcma_so_tumor_filtered,features = "binder_activation6",group.by = "binder_name", pt.size = 0) +
  geom_boxplot(outlier.shape = NA)
# p = p + stat_compare_means(
#   comparisons = list(c("BCMA_B11_l59_int_1525","BCMA_Abecma"))
# )
p

VlnPlot(bcma_so_tumor_filtered,features = c("IL2","IFNG","GZMB","TNF"),group.by = "binder_name")

DoHeatmap(bcma_so_tumor_filtered,features = unlist(module_gene_list))


library(msigdbr)
m_df = msigdbr(species = "Homo sapiens", category = "C5") # C5 is the Gene Ontology category
cytotoxicity_geneset = m_df %>% filter(gs_name == "GOBP_T_CELL_MEDIATED_CYTOTOXICITY")
cytotoxicity_gene_list = unique(cytotoxicity_geneset$gene_symbol)

cytotoxicity_module_list = list(
  Cytotoxicity_Module = cytotoxicity_gene_list
)

bcma_so_tumor_filtered = AddModuleScore(
  object = bcma_so_tumor_filtered,
  features = cytotoxicity_module_list,
  name = 'Cytotox'
)

VlnPlot(bcma_so_tumor_filtered,features = c("GZMB","IL2"),group.by = "binder_name")

## Make the Average Expression Plot

# t_cell_features = c("TCF1","FOXO1","BCL2","CCR7",
#                     "GZMB","PRF1","IFNG","IL2",
#                     "PDCD1","LAG3","TOX","PRDM1",
#                     "TIM3","TCF7","TNF"
#                     )

# t_cell_features = c("TCF1","FOXO1","BCL2","CCR7",
#                     "GZMB","PRF1","IFNG","IL2",
#                     "PDCD1","LAG3","TOX","PRDM1",
#                     "TIM3","TCF7","TNF"
# )

t_cell_features = c("CSF1", "CCL20", "IL2", "TNF", "PRF1", "IFNG", "GZMB",  "CCL1", "CCL3", "CCL4")


t_cell_features_3 = c("TCF1","BCL2","GZMB","IFNG","IL2","TOX","PDCD1","PRDM1","CD69")
t_cell_features_2 = c("GZMA","GZMB","IL2","IFNG","TNF",
                      "TCF7","CCR7","BCL2","FOXO1","PRF1",
                      "TBX21","PRDM1","CD45RA","PDCD1","LAG3",
                      "TOX","EGR2","CD69","IL2RA","CD28","CD137","NFKB","PGC1A","HIF1A",
                      "EOMES","IL21")

avg_expression <- AverageExpression(
  object = bcma_so_tumor_filtered,
  features = t_cell_features,      # The genes we want to plot (all 10 in this example)
  group.by = "binder_name",       # Grouping variable from metadata
  assays = "RNA",                 # Use the 'RNA' assay (normalized data)
  slot = "data"                   # Use the 'data' slot (normalized log-transformed data)
)

# Extract the expression matrix (it's a list, we need the RNA element)
avg_expression_matrix <- avg_expression$RNA

p1 = FeaturePlot(bcma_so_tumor_filtered,features=c("PLAUR"))
p2 = DimPlot(bcma_so_tumor_filtered,group.by = "binder_name")
cowplot::plot_grid(p1,p2)

# --- 3. Prepare Data for ComplexHeatmap ---

# Optional: Scale the data before plotting the heatmap for better visualization
# Scaling transforms values to mean=0 and SD=1, making gene-to-gene comparisons easier.
scaled_matrix <- t(scale(t(avg_expression_matrix)))
# Define a color palette
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Create the heatmap
library(ComplexHeatmap)
heatmap_plot <- Heatmap(
  scaled_matrix,
  name = "Scaled Average Expression", # Name for the legend
  
  # Clustering/Aesthetics for Rows (Genes)
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_title = "Genes",
  
  # Clustering/Aesthetics for Columns (Binder Groups)
  cluster_columns = FALSE, # We usually don't cluster groups defined by metadata
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Average Expression by Binder Name",
  
  # Color
  col = col_fun,
  
  # Optional: Draw borders
  border = TRUE,
  
  # Optional: Adjust cell size (default is good for small matrix)
  width = unit(5, "cm"),
  height = unit(10, "cm")
)

# Print the heatmap
print(heatmap_plot)

VlnPlot(bcma_so_tumor_filtered,features=c("GZMB","IL2","PDCD1","IFNG","PRDM1"))

ggplot(bcma_so_tumor_filtered@meta.data,aes(y=Cytotox1,x=binder_name)) +
  geom_jitter()

DotPlot(bcma_so_tumor_filtered,features = t_cell_features)

## Unused ASA scores

## Load starCAR outputs
bcma_so_starCAT_programs = read.table(paste0("./data/starCAT/BCMA_output/bcma.rf_usage_normalized.txt"))
bcma_so_starCAT_scores = read.table(paste0("./data/starCAT/BCMA_output/bcma.scores.txt"))
bcma_so_tumor_filtered = bcma_so_tumor_filtered %>%
  AddMetaData(metadata = bcma_so_starCAT_programs) %>%
  AddMetaData(metadata = bcma_so_starCAT_scores)

## Plot the proportion of cells predicted to have underwent antigen-specific activation
bcma_asa_binary_data = bcma_so_tumor_filtered@meta.data %>% # <-- Start the pipe with the extracted data.frame
  dplyr::count(binder_name, ASA_binary) %>%
  dplyr::group_by(binder_name) %>%
  dplyr::mutate(Proportion = n / sum(n))# %>%

bcma_asa_plot_df = bcma_asa_binary_data %>%
  filter(ASA_binary == "True") %>%
  group_by(binder_name) %>%
  mutate(
    Total_N = sum(bcma_asa_binary_data$n[bcma_asa_binary_data$binder_name == binder_name]), # Recalculate Total_N using original data frame
    Success_N = n
  ) %>%
  ungroup() %>%
  mutate(
    SE = sqrt(Proportion * (1 - Proportion) / Total_N),
    CI_lower = Proportion - 1.96 * SE,
    CI_upper = Proportion + 1.96 * SE
  ) %>%
  mutate(
    CI_lower = pmax(0, CI_lower),
    CI_upper = pmin(1, CI_upper)
  )

bcma_asa_plot = ggplot(bcma_asa_plot_df, aes(x = binder_name, y = Proportion*100)) +
  # Add bar geometry
  geom_col(fill = "white", color = "black", width = 0.7) +
  # Add error bars using the calculated confidence intervals
  geom_errorbar(aes(ymin = CI_lower*100, ymax = CI_upper*100),
                width = 0.2, # Width of the horizontal lines at the end of the error bar
                color = "black",
                linewidth = 0.8) +
  pretty_plot() + L_border() + ylim(0,100) +
  labs(x="Binder Name",y="% Antigen-specific Activated Cells")
bcma_asa_plot
cowplot::ggsave2("./plots/BCMA/percent_antigen_specfic_activation.pdf",bcma_asa_plot,dpi=300,width=8,height=3)
cowplot::ggsave2("./plots/BCMA/percent_antigen_specfic_activation.png",bcma_asa_plot,dpi=300,width=8,height=3)

## Make a heatmap of the means
functional_modules = c(
  "CellCycle.G2M","CellCycle.S","CellCycle.Late.S","Cytoskeleton","Cytotoxic",
  "ISG","Exhaustion","ASA","Proliferation","Translation","HLA"
)
functional_modules = c(colnames(bcma_so_starCAT_programs),"ASA","Proliferation")
test_plot_data <- bcma_so_tumor_filtered@meta.data %>%
  dplyr::select(all_of("binder_name"), all_of(functional_modules)) %>%
  group_by(binder_name) %>% 
  summarise(across(all_of(functional_modules), mean, .names = "Mean_{.col}"),
            .groups = 'drop')
heatmap_matrix_raw <- test_plot_data %>%
  as.data.frame() %>%
  column_to_rownames(var = "binder_name") %>%
  t() %>%
  as.matrix()
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_raw)))
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(
  heatmap_matrix_scaled,
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  # Titles and Legend
  name = "Z-Score",
  column_title = "Mean Functional Modules by Binder Name",
  
  # Colors
  col = col_fun,
  
  # Text appearance
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10)
)

# --- 4. Draw the Heatmap ---
draw(ht)

asa_marker_bcma = VlnPlot(bcma_so_tumor_filtered, features = c("ASA","Proliferation"),ncol=2)
cowplot::ggsave2("./plots/BCMA/asa_markers.png",asa_marker_bcma,dpi=300,width=8,height=4)

