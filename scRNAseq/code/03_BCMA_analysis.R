library(Seurat)
library(dplyr)
library(data.table)
library(BuenColors)

## Load data
bcma_so_tumor_filtered = readRDS("../data/seurat_objects/bcma_so_tumor_filtered.rds")

bcma_so_tumor_filtered@meta.data$binder_name_simple = dplyr::recode(
  bcma_so_tumor_filtered@meta.data$binder_name,
  !!!c("BCMA_Abecma"="Abecma","BCMA_561726_WT"="B5","BCMA_B11_I59_int_1525"="B5.I0",
       "BCMA_A2_nonint_0366"="B5.N6","BCMA_B4_I59_nonint_1399"="B5.N9")
)
bcma_so_tumor_filtered@meta.data$binder_name_simple = factor(
  bcma_so_tumor_filtered@meta.data$binder_name_simple,
  levels = c(
    "Abecma","B5","B5.I0","B5.N6","B5.N9"
  ))

Idents(bcma_so_tumor_filtered) = "binder_name_simple"
selected_markers_bcma = VlnPlot(bcma_so_tumor_filtered, features = c("IL2","IFNG","GZMA","GZMB"),ncol=4,split.by = "cd_type")
selected_markers_bcma
cowplot::ggsave2("./plots/BCMA/selected_markers.png",selected_markers_bcma,dpi=300,width=12,height=4)

## Check CD type proportions
bcma_cd_prop_df = table(bcma_so_tumor_filtered@meta.data$binder_name_simple,bcma_so_tumor_filtered@meta.data$cd_type) %>% as.data.frame()
bcma_cd_prop_df = bcma_cd_prop_df %>% group_by(Var1) %>% mutate(proportion = Freq / sum(Freq)) %>%ungroup()
bcma_cd_prop_plot = ggplot(bcma_cd_prop_df, aes(x = Var1, y = proportion, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values=c("CD4+"="dodgerblue3","CD8+"="firebrick")) +
  scale_y_continuous(labels = scales::percent,expand = c(0,0)) +
  pretty_plot() + L_border() +
  theme(axis.title = element_blank(),legend.position = "none",text = element_text(size=7))
  #labs(y = "Proportion of T Cells", x = "Binder Name", fill = "CD Type") +

cowplot::ggsave2("../plots/BCMA/bcma_cd_prop_plot.pdf",bcma_cd_prop_plot,dpi=300,width=1.7,height=1.3)

bcma_cd_prop_df

# Get CD4 proportion per binder, normalized to B5.N6
ref <- "B5.N6"

bcma_cd4_se_df <- bcma_so_tumor_filtered@meta.data %>%
  mutate(is_cd4 = as.integer(cd_type == "CD4+")) %>%
  group_by(binder_name_simple) %>%
  summarise(
    n       = n(),
    prop    = mean(is_cd4),
    se      = sqrt(prop * (1 - prop) / n),
    .groups = "drop"
  ) %>%
  mutate(
    ref_prop = prop[binder_name_simple == ref],
    norm     = prop / ref_prop,
    norm_ci  = 1.96 * se / ref_prop
  )

bcma_cd4_se_df$binder_name_simple = factor(bcma_cd4_se_df$binder_name_simple,levels=c(
  "Abecma","B5","B5.I0","B5.N6","B5.N9"
))
#"Abecma","B5","B5.I0","B5.N6","B5.N9"

bcma_cd4_norm_plot <- ggplot(bcma_cd4_se_df %>% filter(binder_name_simple != ref),
                             aes(x = binder_name_simple, y = norm)) +
  #geom_errorbar(aes(ymin = norm - norm_ci, ymax = norm + norm_ci), width = 0.2) +
  geom_errorbar(aes(ymin = norm - norm_ci, ymax = norm + norm_ci), width = 0.5) +
  geom_point(size = 0.8) +
  ylim(0.9, 1.28) +
  geom_hline(yintercept = 1, linetype = 2) +
  pretty_plot() + L_border() + #labs(y = "CD4 Prop. / (B5.N6 CD4 Prop.)") +
  theme(axis.title = element_blank(), text = element_text(size = 7))
bcma_cd4_norm_plot
cowplot::ggsave2("../plots/BCMA/bcma_cd4_prop_norm_plot.pdf",bcma_cd4_norm_plot,dpi=300,width=1.3,height=1.3)

library(purrr)

binders <- levels(bcma_cd4_se_df$binder_name_simple)
cell_df  <- bcma_so_tumor_filtered@meta.data %>%
  mutate(is_cd4 = as.integer(cd_type == "CD4+")) %>%
  filter(binder_name_simple %in% binders)

# All pairwise prop.test
pairwise_df <- combn(binders, 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    d <- cell_df %>% filter(binder_name_simple %in% pair) %>%
      group_by(binder_name_simple) %>%
      summarise(x = sum(is_cd4), n = n(), .groups = "drop")
    p <- prop.test(x = d$x, n = d$n)$p.value
    tibble(group1 = pair[1], group2 = pair[2], pval = p)
  }) %>%
  mutate(
    padj  = p.adjust(pval, method = "BH"),
    label = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns")
  )

pairwise_df

## Inspect clusters
bcma_dim_plot = DimPlot(bcma_so_tumor_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)
bcma_dim_plot
cowplot::ggsave2("./plots/BCMA/Dimplot.png",bcma_dim_plot,dpi=300,width=12,height=4)

## Check CD3
bcma_so_tumor_filtered@meta.data$binder_name = factor(bcma_so_tumor_filtered@meta.data$binder_name,levels = c(
  "BCMA_Abecma","BCMA_561726_WT","BCMA_B11_I59_int_1525","BCMA_A2_nonint_0366","BCMA_B4_I59_nonint_1399"
))
## Draw the dots properly
bcma_dots_simple = DimPlot(bcma_so_tumor_filtered, label = FALSE, group.by = c( "binder_name"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  #scale_color_manual(values = c("#FFB81C","#D91E18","#8B0000","#9966CC","#00BFFF")) + 
  scale_color_manual(values = c("#FFB81C","#D91E18","#8B0000","#9966CC","#2CA02C")) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")
bcma_dots_simple
cowplot::ggsave2(bcma_dots_simple, file = "../plots/BCMA/umap_base_color_updated.png", width = 6, height = 6, dpi = 300)


bcma_dots_simple

# ## Check on other targets
# featuers_to_check = c(
#   "PLAUR","CD33","IL3RA","CD47","CD70","CLEC12A","HAVCR2","FLT3","CD38","BST1","CD200","LILRB4","CD70",
# )
# set2 = c(
#   "BST"
# )
# FeaturePlot(bcma_so_tumor_filtered,features = featuers_to_check)

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
bcma_feature_plots
cowplot::ggsave2(bcma_feature_plots,file="./plots/BCMA/umap_features.png",width=6*2,height=6*2, dpi=600)

## Make more feature plots for supplements
bcma_feature_plots_extended = cowplot::plot_grid(
  mk_plot(bcma_so_tumor_filtered,"CD4"),
  mk_plot(bcma_so_tumor_filtered,"CD8A"),
  mk_plot(bcma_so_tumor_filtered,"GZMB"),
  ncol=2,scale=1
)
bcma_feature_plots_extended
cowplot::ggsave2(bcma_feature_plots_extended,file="../plots/BCMA/umap_features_extended.png",width=6*2,height=6*2, dpi=600)


## Split activated from non-activated
Idents(bcma_so_tumor_filtered)="binder_name"
bcma_so_tumor_filtered@meta.data$binder_is_activated = bcma_so_tumor_filtered@meta.data$binder_name %in% c("BCMA_561726_WT","BCMA_B11_I59_int_1525","BCMA_Abecma")

split_activate_plot = DimPlot(bcma_so_tumor_filtered,reduction="umap_harmony",split.by = "binder_is_activated")

split_activate_plot

bcma_marker_data = FetchData(bcma_so_tumor_filtered,
                             vars = c("IL2","IFNG","TNF","GZMA","GZMB","binder_name"),
                             #slot = "data"
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

## Perform DEG analysis between all evolved binders to parental (Per cell types)
evolved_binders = c("BCMA_B4_I59_nonint_1399","BCMA_B11_I59_int_1525","BCMA_A2_nonint_0366")
for (binder in evolved_binders) {
  for (cd_type in c("CD4+","CD8+")) {
    #sub_so = subset(bcma_so_tumor_filtered,cd_type==cd_type)
    clean_name = substr(cd_type,1,3)
    evolved_parental_comparison = compare_samples_with_enrichment(
      bcma_so_tumor_filtered,
      sample_col="binder_name",
      sample1=binder,
      sample2="BCMA_561726_WT",
      subset_col="cd_type",
      subset_value=cd_type
    )
    cowplot::ggsave2(paste0("../plots/BCMA/volcano_",binder,"_vs_BCMA_561726_WT_",clean_name,".pdf"),evolved_parental_comparison$volcano_plot,dpi=300,width=8,height=4)
    cowplot::ggsave2(paste0("../plots/BCMA/volcano_",binder,"_vs_BCMA_561726_WT_",clean_name,".png"),evolved_parental_comparison$volcano_plot,dpi=300,width=8,height=4)
    cowplot::ggsave2(paste0("../plots/BCMA/fgsea_",binder,"_vs_BCMA_561726_WT_",clean_name,".pdf"),evolved_parental_comparison$fgsea_plot,dpi=300,width=8,height=8)
    write.csv(evolved_parental_comparison$markers,paste0("../data/degs/degs_",binder,"_vs_BCMA_561726_WT_",clean_name,".csv"))
  }
}

## Check if CD4/CD8 effect sizes are consistent
bcma_evolved_cd4_degs = read.csv("../data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_561726_WT_CD4.csv")
colnames(bcma_evolved_cd4_degs) <- paste(colnames(bcma_evolved_cd4_degs), "CD4", sep = "_")
bcma_evolved_cd8_degs = read.csv("../data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_561726_WT_CD8.csv")
colnames(bcma_evolved_cd8_degs) <- paste(colnames(bcma_evolved_cd8_degs), "CD8", sep = "_")
bcma_evolved_degs_combined = merge(bcma_evolved_cd4_degs,bcma_evolved_cd8_degs,by.x="gene_CD4",by.y="gene_CD8",all.x=TRUE)
bcma_evolved_degs_combined

bcma_evolved_degs_combined$sig_in_cd4 = bcma_evolved_degs_combined$p_val_adj_CD4 < 0.05
bcma_evolved_degs_combined_annotate_subset = bcma_evolved_degs_combined %>% filter(
  (avg_log2FC_CD4 < -1.25) | (avg_log2FC_CD4 > 1) | (avg_log2FC_CD8 > 1.3) 
)
effect_size_comparison_plot = ggplot(bcma_evolved_degs_combined,aes(x=avg_log2FC_CD4,y=avg_log2FC_CD8)) +
  geom_point(aes(color=sig_in_cd4)) + pretty_plot() + L_border() +
  geom_abline(slope=1,linetype="dashed",color="gray") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text_repel(data=bcma_evolved_degs_combined_annotate_subset,aes(label=gene_CD4)) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  labs(x="log2FC(Evolved CD4/ Parental CD4)",y="log2FC(Evolved CD8/ Parental CD8)") +
  theme(legend.position = "none")

effect_size_comparison_plot

ggsave("../plots/BCMA/CD4_CD8_effect_size_comparison_plot.pdf",effect_size_comparison_plot,dpi=300)

## Make a new one that only contains sig hits
bcma_evolved_degs_combined_sig_only = bcma_evolved_degs_combined %>% filter(
  (p_val_adj_CD4 < 0.05) | (p_val_adj_CD8 < 0.05)
)
bcma_evolved_degs_combined_sig_only$to_highlight = bcma_evolved_degs_combined_sig_only$gene_CD4 %in% c(
  "IL2","ENTPD1","IFNG","PDCD1"
)
effect_size_comparison_plot_sig_only = ggplot(bcma_evolved_degs_combined_sig_only,aes(x=avg_log2FC_CD4,y=avg_log2FC_CD8)) +
  geom_point(aes(color=to_highlight)) + pretty_plot() + L_border() +
  geom_abline(slope=1,linetype="dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_smooth(method = "lm", se = TRUE,) +
  #geom_text_repel(data=bcma_evolved_degs_combined_sig_only %>% filter(to_highlight),aes(label=gene_CD4)) +
  scale_color_manual(values=c("TRUE"="firebrick3","FALSE"="black")) +
  labs(x="log2FC(Evolved CD4/ Parental CD4)",y="log2FC(Evolved CD8/ Parental CD8)") +
  theme(legend.position = "none",axis.title = element_blank(),text = element_text(size=7),
        axis.line = element_line(linewidth = 0.5))
effect_size_comparison_plot_sig_only
ggsave("../plots/BCMA/effect_size_comparison_plot_sig_only.pdf",effect_size_comparison_plot_sig_only,dpi=300,width=1.5,height=1.4)


VlnPlot(
  subset(bcma_so_tumor_filtered, binder_name_simple == "B5" | binder_name_simple == "B5.I0"),
  features = c("IL24"),
  split.by = c("cd_type")
  )

VlnPlot(
  subset(bcma_so_tumor_filtered, cd_type == "CD4+" | cd_type == "CD8+"),
  features = c("IL24"),
  split.by = c("binder_name_simple")
)

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
evolved_vs_parental = read.csv("../data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_561726_WT.csv")
tonic_vs_nonactive = read.csv("../data/degs/degs_BCMA_B4_I59_nonint_1399_vs_BCMA_A2_nonint_0366.csv")
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
  "memory" = memory_genes
  #"activation" = general_activation_genes,
  #"binder_activation" = binder_activation_genes,
  #"binder_tonic" = binder_tonic_genes
)

bcma_so_tumor_filtered = AddModuleScore(
  object = bcma_so_tumor_filtered,
  features = module_gene_list,
  name = c("effector","exhaustion","function","memory")#, "activation","binder_activation","binder_tonic")
)

## Visualize differnece in module scores by cell types

bcma_cdtype_module_scores = bcma_so_tumor_filtered@meta.data[,c("binder_name_simple","effector1","exhaustion2","cd_type")]
bcma_cdtype_module_scores_test_results = bcma_cdtype_module_scores %>%
  group_by(binder_name_simple) %>%
  summarise(
    # Means per CD type
    mean_effector1_CD4  = mean(effector1[cd_type == "CD4+"], na.rm = TRUE),
    mean_effector1_CD8  = mean(effector1[cd_type == "CD8+"], na.rm = TRUE),
    mean_exhaustion2_CD4 = mean(exhaustion2[cd_type == "CD4+"], na.rm = TRUE),
    mean_exhaustion2_CD8 = mean(exhaustion2[cd_type == "CD8+"], na.rm = TRUE),
    # Wilcoxon p-values (CD4+ vs CD8+)
    p_effector1   = wilcox.test(
      effector1[cd_type == "CD4+"],
      effector1[cd_type == "CD8+"]
    )$p.value,
    p_exhaustion2 = wilcox.test(
      exhaustion2[cd_type == "CD4+"],
      exhaustion2[cd_type == "CD8+"]
    )$p.value,
    .groups = "drop"
  )
bcma_cdtype_module_scores_test_results[,c("binder_name_simple","p_effector1","p_exhaustion2")]


bcma_cd4_cd8_effector_plot <- ggplot(bcma_cdtype_module_scores,
                                     aes(x = binder_name_simple, y = effector1, fill = cd_type)) +
  geom_violin(position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA) +
  #scale_fill_manual(values=c("CD4+"="dodgerblue3","CD8+"="firebrick")) +
  #scale_fill_manual(values=c("CD4+"="dodgerblue3","CD8+"="firebrick")) +
  scale_fill_manual(values = c("CD4+" = "#2E86AB", "CD8+" = "#E07B39")) +
  pretty_plot() + L_border() +
  theme(axis.title = element_blank(), legend.position = "none")

bcma_cd4_cd8_effector_plot
ggsave("../plots/BCMA/bcma_cd4_cd8_effector_plot.pdf",bcma_cd4_cd8_effector_plot,dpi=300,width=2.6,height=1.3)

bcma_cd4_cd8_exhaustion_plot = ggplot(bcma_cdtype_module_scores, aes(x = binder_name_simple, y = exhaustion2, fill = cd_type)) +
  geom_violin(position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA) +
  #scale_fill_manual(values=c("CD4+"="dodgerblue3","CD8+"="firebrick")) +
  scale_fill_manual(values = c("CD4+" = "#2E86AB", "CD8+" = "#E07B39")) +
  pretty_plot() + L_border() +
  theme(axis.title = element_blank(), legend.position = "none")
  
bcma_cd4_cd8_exhaustion_plot
ggsave("../plots/BCMA/bcma_cd4_cd8_exhaustion_plot.pdf",bcma_cd4_cd8_exhaustion_plot,dpi=300,width=2.6,height=1.3)


VlnPlot(bcma_so_tumor_filtered, pt.size = 0, 
        features = c("effector1","exhaustion2","function3","memory4","activation5"),
        group.by = "binder_name") & geom_boxplot(outlier.shape = NA,)

VlnPlot(bcma_so_tumor_filtered, pt.size = 0, 
        features = c("effector1","exhaustion2"),
        group.by = "binder_name_simple",split.by = "cd_type") & geom_boxplot(outlier.shape = NA,)

binder_tonic_genes

VlnPlot(bcma_so_tumor_filtered,features=c("Treg"))

## Plot exhaustion and function
bcma_module_score_subset = bcma_so_tumor_filtered@meta.data[,c("binder_name_simple","effector1","exhaustion2")] %>%
  pivot_longer(c(exhaustion2, effector1), names_to = "Metric", values_to = "Score") %>%
  mutate(binder_name_simple = factor(binder_name_simple))

library(gghalves)
effector_data <- subset(bcma_module_score_subset, Metric == "effector1")
exhaustion_data <- subset(bcma_module_score_subset, Metric == "exhaustion2")

##
bcma_metric_group_means = bcma_module_score_subset %>% group_by(binder_name_simple,Metric) %>% summarize(Mean_Score = mean(Score, na.rm = TRUE)) %>% arrange(Metric)
bcma_metric_group_means_pvals = bcma_module_score_subset %>%
  group_by(Metric) %>%
  summarise(
    t_test_results = list(
      pairwise.t.test(
        x = Score, 
        g = binder_name_simple, 
        p.adjust.method = "bonferroni"
      ) %>% tidy()
    ),
    .groups = "drop"
  ) %>%
  unnest(t_test_results)
bcma_metric_group_means_pvals

bcma_exhaustion_vs_effector_scores = ggplot(bcma_module_score_subset, aes(x = binder_name_simple, y = Score, fill = Metric)) +
  geom_half_violin(
    data = subset(bcma_module_score_subset, Metric == "effector1"), side = "l",
    position = position_nudge(x = -0.0), trim = TRUE#, alpha = 0.7
  ) +
  geom_half_violin(
    data = subset(bcma_module_score_subset, Metric == "exhaustion2"), side = "r",
    position = position_nudge(x = 0.0), trim = TRUE#, alpha = 0.7
  ) +
  # geom_boxplot(
  #   aes(group = interaction(binder_name_simple, Metric)), width = 0.1,
  #   outlier.shape = NA, position = position_dodge(width = 0.1)
  # ) +
  scale_fill_manual(
    values = c(effector1 = "#1f78b4", exhaustion2 = "#e31a1c"),
    labels = c("Activation (Left)","Exhaustion (Right)")
  ) +
  labs(
    x = "Binder Name",
    y = "Module Score",
    fill = "Metric"
  ) +
  pretty_plot() + L_border() + theme(legend.position="none",axis.title.x = element_blank())
bcma_exhaustion_vs_effector_scores
cowplot::ggsave2("../plots/BCMA/bcma_module_score_comparison.pdf",bcma_exhaustion_vs_effector_scores,dpi=300,height=1.6,width=2.4)


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

## Specifically visualize the genes highlughted int he volcano

library(BuenColors)
bcma_so_tumor_filtered@meta.data$binder_name_simple = dplyr::recode(
  bcma_so_tumor_filtered@meta.data$binder_name,
  !!!c("BCMA_Abecma"="Abecma","BCMA_561726_WT"="B5","BCMA_B11_I59_int_1525"="B5.I0",
       "BCMA_A2_nonint_0366"="B5.N6","BCMA_B4_I59_nonint_1399"="B5.N9")
)
bcma_so_tumor_filtered@meta.data$binder_name_simple = factor(
  bcma_so_tumor_filtered@meta.data$binder_name_simple,
  levels = c(
    "Abecma","B5","B5.I0","B5.N6","B5.N9"
  ))
#bcma_so_tumor_filtered@meta.data$IL2_counts = GetAssayData(object = bcma_so_tumor_filtered, assay = "RNA", slot = "data")[c("IL2"), ]
#bcma_so_tumor_filtered@meta.data$IFNG_counts = GetAssayData(object = bcma_so_tumor_filtered, assay = "RNA", slot = "data")[c("IFNG"), ]
bcma_so_tumor_filtered@meta.data$IL2_counts = GetAssayData(object = bcma_so_tumor_filtered, assay = "RNA")[c("IL2"), ]
bcma_so_tumor_filtered@meta.data$IFNG_counts = GetAssayData(object = bcma_so_tumor_filtered, assay = "RNA")[c("IFNG"), ]

## Alternative color scheme
bcma_color_mapping <- c(
  "m971" = "#FFB81C",#"#8B0000",
  "B5" = "#D91E18",
  "B5.I0" = "#8B0000",
  "B5.N6" = "#9966CC",
  "B5.N9" = "#00BFFF"
)

bcma_il2_boxplot = ggplot(bcma_so_tumor_filtered@meta.data,
                                        aes(x = binder_name_simple, y = IL2_counts, fill=binder_name_simple)) + 
  geom_violin(aes(fill=binder_name_simple),linewidth = 0.3) + geom_jitter(size=0.01) +
  pretty_plot(fontsize = 8) + 
  L_border() + 
  labs(x="Binder Name",y="IL2") +
  scale_fill_manual(values = bcma_color_mapping) +
  theme(legend.position = "none", axis.title.x=element_blank())

bcma_ifng_boxplot = ggplot(bcma_so_tumor_filtered@meta.data,
                          aes(x = binder_name_simple, y = IFNG_counts, fill=binder_name_simple)) + 
  geom_violin(aes(fill=binder_name_simple),linewidth = 0.3) + geom_jitter(size=0.01) +
  pretty_plot(fontsize = 8) + 
  L_border() + 
  labs(x="Binder Name",y="IFNG") +
  scale_fill_manual(values = bcma_color_mapping) +
  theme(legend.position = "none", axis.title.x=element_blank())

bcma_il2_ifng_violin_plot = cowplot::plot_grid(bcma_il2_boxplot,bcma_ifng_boxplot,ncol=1)
bcma_il2_ifng_violin_plot
cowplot::ggsave2("../plots/BCMA/bcma_il2_ifng_violin.pdf",bcma_il2_ifng_violin_plot,dpi=300,width=1.8,height=1.8)
cowplot::ggsave2("../plots/BCMA/bcma_il2_ifng_violin.png",bcma_il2_ifng_violin_plot,dpi=300,width=1.8,height=1.8)


VlnPlot(bcma_so_tumor_filtered,features=c("IL2","IFNG"))

# combined_module_score_plots = cowplot::plot_grid(
#   cd22_oe_module_score_boxplot,cd22_rpmi_module_score_boxplot,ncol=1
# )
# combined_module_score_plots
# cowplot::ggsave2("../plots/CD22_OE_RPMI_Activation_Scores.pdf",combined_module_score_plots,dpi=300,width=1.8,height=1.8)
# 

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
bcma_so_starCAT_programs = read.table(paste0("../data/starCAT/BCMA_output/bcma.rf_usage_normalized.txt"))
bcma_so_starCAT_scores = read.table(paste0("../data/starCAT/BCMA_output/bcma.scores.txt"))
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

