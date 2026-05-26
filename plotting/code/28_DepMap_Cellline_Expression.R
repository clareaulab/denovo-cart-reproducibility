library(dplyr)
library(BuenColors)
# 
# depmap_df = read.csv("../data/DepMap_celline_subset.csv")
# depmap_df_melt = depmap_df %>% melt()
# 
# depmap_df_melt$Cellline = factor(
#   depmap_df_melt$Cellline,
#   levels=c(
#     "RAMOS","RAJI",
#     "K562","U937",
#     "MDAMB468",
#     "RPMI8226",
#     "SET2","A549","HELA"
#   )
# )
# 
# depmap_expr_plot <- ggplot(depmap_df_melt,aes(x = Cellline, y = variable, fill = value)) +
#   geom_tile(color = "black", linewidth = 0.25) +
#   scale_fill_gradientn(
#     colors = jdb_palette("solar_rojos"),
#     limits = c(0, max(depmap_df_melt$value)),
#     na.value = "grey90"
#   ) +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   pretty_plot(fontsize = 8) +
#   L_border() +
#   theme(
#     axis.title.x    = element_blank(),
#     axis.title.y    = element_blank(),
#     axis.text = element_text(size=5)
#     #axis.ticks.x   = element_blank(),
#     #axis.ticks.y   = element_blank()
#     #legend.position = "none"
#   )
# 
# depmap_expr_plot
# 
# depmap_expr_plot_clean = depmap_expr_plot + theme(legend.position = "none")
# 
# ggsave("../plots/DepMap_celline_expr_plot.pdf",depmap_expr_plot,dpi=300,width=3.5,height=1)
# ggsave("../plots/DepMap_celline_expr_plot_clean.pdf",depmap_expr_plot_clean,dpi=300,width=3.5,height=1)

llab_celllines = read.csv("/home/chuh/protein_design/denovo-cart-reproducibility/plotting/data/llab_cellline_tpm.csv")
llab_celllines_melt = llab_celllines %>% melt()
llab_celllines_melt$variable = factor(
  llab_celllines_melt$variable,
  levels=c(
    "RAMOS","Raji",
    "DG75","NALM6",
    "RPMI8226",
    "K562","U937",
    "MDAMB468",
    "SET2","A549","HELA",
    "MOLT4"
  )
)
llab_celllines_melt$X = factor(
  llab_celllines_melt$X,
  levels=c("CD22","CD19","TNFRSF17")
)

## Cap the TPM to emphasize which lines ar epositive
cap_val = 150
llab_celllines_melt$value_capped = pmin(llab_celllines_melt$value, cap_val)
llab_cellline_expr_plot = ggplot(llab_celllines_melt,aes(y = X, x = variable, fill = value_capped)) +
  geom_tile(color = "black", linewidth = 0.25) +
  scale_fill_gradientn(
    colors = jdb_palette("solar_rojos"),
    limits = c(0, cap_val),
    #limits = c(0, max(log2(llab_celllines_melt$value+1))),
    na.value = "grey90"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  pretty_plot(fontsize = 8) +
  L_border() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text = element_text(size=5)
  )

llab_cellline_expr_plot
llab_cellline_expr_plot_clean = llab_cellline_expr_plot + theme(legend.position = "none")

ggsave("../plots/llab_cellline_expr_plot.pdf",llab_cellline_expr_plot,dpi=300,width=4.5,height=1.1)
ggsave("../plots/llab_cellline_expr_plot_clean.pdf",llab_cellline_expr_plot_clean,dpi=300,width=4.5,height=1.1)
