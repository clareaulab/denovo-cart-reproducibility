library(Seurat)
library(dplyr)
library(data.table)
setwd("/home/chuh/protein_design/denovo-cart-reproducibility/scRNAseq")

## Load data
cd22_rpmi_tumor_filtered = readRDS("./data/seurat_objects/cd22_rpmi_so_filtered.rds")

## Load starCAT outputs
cd22_rpmi_starCAT_programs = read.table(paste0("./data/starCAT/CD22_RPMI_output/CD22_RPMI.rf_usage_normalized.txt"))
cd22_rpmi_starCAT_scores = read.table(paste0("./data/starCAT/CD22_RPMI_output/CD22_RPMI.scores.txt"))
cd22_rpmi_tumor_filtered = cd22_rpmi_tumor_filtered %>%
  AddMetaData(metadata = cd22_rpmi_starCAT_programs) %>%
  AddMetaData(metadata = cd22_rpmi_starCAT_scores)

## Plot the proportion of cells predicted to have underwent antigen-specific activation
cd22_rpmi_asa_binary_data = cd22_rpmi_tumor_filtered@meta.data %>% # <-- Start the pipe with the extracted data.frame
  dplyr::count(binder_name, ASA_binary) %>%
  dplyr::group_by(binder_name) %>%
  dplyr::mutate(Proportion = n / sum(n))# %>%

cd22_rpmi_asa_plot_df = cd22_rpmi_asa_binary_data %>%
  filter(ASA_binary == "True") %>%
  group_by(binder_name) %>%
  mutate(
    Total_N = sum(cd22_rpmi_asa_binary_data$n[cd22_rpmi_asa_binary_data$binder_name == binder_name]), # Recalculate Total_N using original data frame
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

cd22_rpmi_asa_plot = ggplot(cd22_rpmi_asa_plot_df, aes(x = binder_name, y = Proportion*100)) +
  # Add bar geometry
  geom_col(fill = "white", color = "black", width = 0.7) +
  # Add error bars using the calculated confidence intervals
  geom_errorbar(aes(ymin = CI_lower*100, ymax = CI_upper*100),
                width = 0.2, # Width of the horizontal lines at the end of the error bar
                color = "black",
                linewidth = 0.8) +
  pretty_plot() + L_border() + ylim(0,100) +
  labs(x="Binder Name",y="% Antigen-specific Activated Cells")
cd22_rpmi_asa_plot
cowplot::ggsave2("./plots/CD22_RPMI/percent_antigen_specfic_activation.pdf",cd22_rpmi_asa_plot,dpi=300,width=8,height=3)
cowplot::ggsave2("./plots/CD22_RPMI/percent_antigen_specfic_activation.png",cd22_rpmi_asa_plot,dpi=300,width=8,height=3)

## Plot module scores
functional_modules = c(colnames(cd22_rpmi_starCAT_programs),"ASA","Proliferation")
test_plot_data <- cd22_rpmi_tumor_filtered@meta.data %>%
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

VlnPlot(cd22_rpmi_tumor_filtered,features = c("ASA","Proliferation"))

cd22_so_rpmi_asa_only = subset(cd22_rpmi_tumor_filtered, subset = ASA_binary == "True")
VlnPlot(cd22_so_rpmi_asa_only,features = c("ASA","Proliferation"))

## Inspect clusters
cd22_rpmi_dim_plot = DimPlot(cd22_rpmi_tumor_filtered,
        reduction = "umap_rna",
        group.by = c("seurat_clusters","Phase","binder_name")
)
cd22_rpmi_dim_plot
cowplot::ggsave2("./plots/CD22_RPMI/Dimplot.png",cd22_rpmi_dim_plot,dpi=300,width=12,height=4)

## Show selected markers
selected_markers_cd22_rpmi = VlnPlot(cd22_rpmi_tumor_filtered, features = c("IL2","IFNG","GZMA","GZMB"),ncol=4)
selected_markers_cd22_rpmi
cowplot::ggsave2("./plots/CD22_RPMI/selected_markers.png",selected_markers_cd22_rpmi,dpi=300,width=12,height=4)

asa_marker_cd22_rpmi = VlnPlot(cd22_rpmi_tumor_filtered, features = c("ASA","Proliferation"),ncol=2)
asa_marker_cd22_rpmi
cowplot::ggsave2("./plots/CD22_RPMI/asa_markers.png",asa_marker_cd22_rpmi,dpi=300,width=8,height=4)

## Find Markers
cd22_rpmi_tumor_filtered_clust_markers = FindAllMarkers(cd22_rpmi_tumor_filtered,assay="RNA",group.by="seurat_clusters")
cd22_rpmi_tumor_filtered_binder_markers = FindAllMarkers(cd22_rpmi_tumor_filtered,assay="RNA",group.by="binder_name")

## Show markers by cluster
top_cd22_rpmi_markers_by_clust = cd22_rpmi_tumor_filtered_clust_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

Idents(cd22_rpmi_tumor_filtered) = "seurat_clusters"
top_cd22_rpmi_markers_by_clust_heatmap = DoHeatmap(
  cd22_rpmi_tumor_filtered,
  features = top_cd22_rpmi_markers_by_clust$gene
) + NoLegend()
top_cd22_rpmi_markers_by_clust_heatmap
ggsave("./plots/top_cd22_rpmi_markers_by_rna_clusters.png",top_cd22_rpmi_markers_by_clust_heatmap,dpi=300,height=8,width=12)

## Show markers by binder names
top_cd22_rpmi_markers_by_binders = cd22_rpmi_tumor_filtered_binder_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

Idents(cd22_rpmi_tumor_filtered) = "binder_name"
top_cd22_rpmi_markers_by_binder_heatmap = DoHeatmap(
  cd22_rpmi_tumor_filtered,
  features = top_cd22_rpmi_markers_by_binders$gene
) + NoLegend()
top_cd22_rpmi_markers_by_binder_heatmap
ggsave("./plots/top_cd22_rpmi_markers_by_binders.png",top_cd22_rpmi_markers_by_binder_heatmap,dpi=300,height=8,width=12)

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
      fill = "Significance",
      title = "FGSEA Hallmark Enrichment"
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
  print("plottig fgsea")
  # plot fgsea
  fgsea_plot <- plot_fgsea_bar(fgsea_res, top_n = top_fgsea)
  print("Done plotting fgsea")
  return(list(
    markers = de_out$markers,
    volcano_plot = de_out$volcano_plot,
    fgsea = fgsea_res,
    fgsea_plot = fgsea_plot
  ))
}

## Perform DEG analysis between all de novo binders to ABECMA
de_novo_binders = c("CD22_i61_2","CD22_B1_nonint_2261","CD22_A3_nonint_4434","CD22_A11_int_1411")
for (binder in de_novo_binders) {
  denovo_m971_comparison = compare_samples_with_enrichment(
    cd22_rpmi_tumor_filtered,
    sample_col="binder_name",
    sample1=binder,
    sample2="CD22_m971",
    subset_col=NULL,
    subset_value=NULL
  )
  cowplot::ggsave2(paste0("./plots/CD22_RPMI/volcano_",binder,"_vs_CD22_m971.pdf"),denovo_m971_comparison$volcano_plot,dpi=300,width=8)
  cowplot::ggsave2(paste0("./plots/CD22_RPMI/volcano_",binder,"_vs_CD22_m971.png"),denovo_m971_comparison$volcano_plot,dpi=300,width=8)
  cowplot::ggsave2(paste0("./plots/CD22_RPMI/fgsea_",binder,"_vs_CD22_m971.pdf"),denovo_m971_comparison$fgsea_plot,dpi=300,width=8)
  write.csv(denovo_m971_comparison$markers,paste0("./data/degs/degs_",binder,"_vs_CD22_m971.csv"))
}

## Perform DEG analysis between all evolved binders to parental
cd22_evolved_binders = c("CD22_B1_nonint_2261","CD22_A3_nonint_4434","CD22_A11_int_1411")
for (binder in cd22_evolved_binders) {
  cd22_evolved_parental_comparison = compare_samples_with_enrichment(
    cd22_rpmi_tumor_filtered,
    sample_col="binder_name",
    sample1=binder,
    sample2="CD22_i61_2",
    subset_col=NULL,
    subset_value=NULL
  )
  cowplot::ggsave2(paste0("./plots/CD22_RPMI/volcano_",binder,"_vs_CD22_i61_2.pdf"),cd22_evolved_parental_comparison$volcano_plot,dpi=300)
  cowplot::ggsave2(paste0("./plots/CD22_RPMI/volcano_",binder,"_vs_CD22_i61_2.png"),cd22_evolved_parental_comparison$volcano_plot,dpi=300)
  cowplot::ggsave2(paste0("./plots/CD22_RPMI/fgsea_",binder,"_vs_CD22_i61_2.pdf"),cd22_evolved_parental_comparison$fgsea_plot,dpi=300)
  write.csv(cd22_evolved_parental_comparison$markers,paste0("./data/degs/degs_",binder,"_vs_CD22_i61_2.csv"))
}

## Load the overexpressed genes from evolved to wildtype as binder activation module score
evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_561726_WT.csv")
#evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_B4_I59_nonint_1399.csv")
#evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_A2_nonint_0366.csv")
#evolved_vs_parental = read.csv("./data/degs/degs_BCMA_B11_I59_int_1525_vs_BCMA_Abecma.csv")

# tcat_ref = t(read.table("./data/starCAT/TCAT.V1.reference.tsv")) %>% as.data.frame()
# binder_activation_genes = tcat_ref %>% filter(Cytotoxic > 0) %>% arrange(-Cytotoxic) %>% rownames()
# binder_activation_genes = binder_activation_genes#[1:500]


binder_activation_genes = evolved_vs_parental %>%
  filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>%
  pull(gene)
binder_activation_genes

tonic_vs_nonactive = read.csv("./data/degs/degs_BCMA_B4_I59_nonint_1399_vs_BCMA_A2_nonint_0366.csv")
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

module_gene_list = list(
  "effector" = effector_genes,
  "exhaustion" = exhaustion_genes,
  "function" = function_genes,
  "memory" = memory_genes,
  "activation" = general_activation_genes,
  "binder_activation" = binder_activation_genes,
  "binder_tonic" = binder_tonic_genes
)

cd22_rpmi_tumor_filtered = AddModuleScore(
  object = cd22_rpmi_tumor_filtered,
  features = module_gene_list,
  name = c("effector","exhaustion","function","memory", "activation","binder_activation","binder_tonic")
)

VlnPlot(cd22_rpmi_tumor_filtered, pt.size = 0,
        features = c("effector1","exhaustion2","function3","memory4","binder_activation6","binder_tonic7"),
        group.by = "binder_name")  & geom_boxplot(outlier.shape = NA,)


library(BuenColors)
cd22_rpmi_tumor_filtered@meta.data$binder_name_simple = dplyr::recode(
  cd22_rpmi_tumor_filtered@meta.data$binder_name,
  !!!c("CD22_m971"="m971","CD22_i61_2"="D1","CD22_B1_nonint_2261"="D1.N1",
       "CD22_A3_nonint_4434"="D1.N0","CD22_A11_int_1411"="D1.I3")
)
cd22_rpmi_tumor_filtered@meta.data$binder_name_simple = factor(
  cd22_rpmi_tumor_filtered@meta.data$binder_name_simple,
  levels = c(
    "m971","D1","D1.N0","D1.N1","D1.I3"
  ))

## Alternative color scheme
cd22_color_mapping <- c(
  "m971" = "#d62728",#"#8B0000",
  "D1" = "#2CA02C",
  "D1.N0" = "#176017",
  "D1.N1" = "#FFD700",
  "D1.I3" = "orange"
)


cd22_rpmi_module_score_boxplot = ggplot(cd22_rpmi_tumor_filtered@meta.data,
                                      aes(x = binder_name_simple, y = binder_activation6, fill=binder_name_simple)) + 
  geom_violin(aes(fill=binder_name_simple),linewidth = 0.3) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA, width = 0.6,lwd=0.1) + 
  pretty_plot(fontsize = 8) + 
  L_border() + theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("m971","D1.N0"),c("m971","D1"),c("D1","D1.N0")),size = 1.5
  ) +
  labs(x="Binder Name",y="Activation Score") +
  scale_fill_manual(values = cd22_color_mapping)
#scale_fill_manual(values = c("orange", "black", "dodgerblue3", "firebrick", "firebrick"))
cd22_rpmi_module_score_boxplot
cowplot::ggsave2("../plots/CD22_RPMI/minibinder_module_score.pdf",cd22_rpmi_module_score_boxplot,dpi=300,width=1.8,height=1.6)


combined_module_score_plots = cowplot::plot_grid(
  cd22_oe_module_score_boxplot,cd22_rpmi_module_score_boxplot,ncol=1
)
combined_module_score_plots
cowplot::ggsave2("../plots/CD22_OE_RPMI_Activation_Scores.pdf",combined_module_score_plots,dpi=300,width=1.8,height=1.8)



VlnPlot(cd22_rpmi_tumor_filtered,features = c("activation5"),group.by = "binder_name")
VlnPlot(cd22_rpmi_tumor_filtered,features = c("FASLG"),group.by = "binder_name")


library(msigdbr)
# m_df = msigdbr(species = "Homo sapiens", category = "C5") # C5 is the Gene Ontology category
# cytotoxicity_geneset = m_df %>% filter(gs_name == "GOBP_T_CELL_MEDIATED_CYTOTOXICITY")
# cytotoxicity_gene_list = unique(cytotoxicity_geneset$gene_symbol)

m_df = msigdbr(species = "Homo sapiens", category = "H") # C5 is the Gene Ontology category
il2_stat5_gene_set = m_df %>% filter(gs_name == "HALLMARK_IL2_STAT5_SIGNALING")
il2_stat5_gene_list = unique(il2_stat5_gene_set$gene_symbol)
#cytotoxicity_geneset = m_df %>% filter(gs_name == "GOBP_T_CELL_MEDIATED_CYTOTOXICITY")
#cytotoxicity_gene_list = unique(cytotoxicity_geneset$gene_symbol)

# 
cytotoxicity_module_list = list(
  il2_stat5_module = il2_stat5_gene_list,
  custom_set = c("MT2A","MT1E","MTBP","MYB")
  #Cytotoxicity_Module = cytotoxicity_gene_list
)

cd22_rpmi_tumor_filtered = AddModuleScore(
  object = cd22_rpmi_tumor_filtered,
  features = cytotoxicity_module_list,
  name = 'Cytotox'
)

VlnPlot(cd22_rpmi_tumor_filtered,features = c("Cytotox1","Cytotox2"),group.by = "binder_name")

# 
# VlnPlot(cd22_rpmi_tumor_filtered,features = c("Cytotox1","GZMB","IL2"),group.by = "binder_name")

