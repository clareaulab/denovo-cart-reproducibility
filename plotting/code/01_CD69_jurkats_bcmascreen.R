library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

fread("../data/jurkat_data_BCMA_screen.txt") %>%
  mutate(class = case_when(
    (bK562_BCMA - cCAR_alone > 25) &  (bK562_BCMA / cCAR_alone > 2) ~ "Good",
    cCAR_alone > 25 ~ "Tonic", 
    TRUE ~ "bad"
  )) %>% 
  arrange(desc(class), bK562_BCMA - cCAR_alone ) %>%
  mutate(BINDER = factor(as.character(BINDER), levels = as.character(BINDER))) %>%
  reshape2::melt(id.vars = c("BINDER", "class")) %>%
  ggplot(aes(x = BINDER, y = variable, fill = value)) + 
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() + theme(axis.text.x=element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") -> p1
p1
cowplot::ggsave2(p1, file = "../plots/grid_cd69_jurkats_bcma.pdf", width = 3.7, height = 1.3)

process_line <- function(what1){
  oldf <- fread("../data/primary_bcma_screen_cytokine.tsv") %>%
    filter(line == what1)
  
  mean_df <- oldf %>%
    reshape2::melt(id.vars = c("what", "line", "replicate")) %>%
    group_by(variable, line, what) %>%
    dplyr::summarize(meanV = mean(value), count = n())
  mean_df
  
}

fread("../data/primary_bcma_screen_cytokine.tsv") %>%
  reshape2::melt(id.vars = c("what", "line", "replicate")) %>%
  group_by(variable, line, what) %>%
  dplyr::summarize(meanV = mean(value), count = n()) %>%
  mutate(line = factor(as.character(line), levels = rev(c("CARalone", "K562","RPMI8226", "Raji",  "K562BCMA")))) %>%
  mutate(variable = factor(as.character(variable), levels = (c("UTD", "no_binder", "298275_mpnn", "561726_mpnn", "abecma")))) %>%
  ggplot(aes(x = variable, y = line, fill = meanV)) + 
  geom_tile(color = "black") + facet_wrap(~what)+
  scale_fill_gradientn(colors = jdb_palette("solar_blues"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + theme(axis.text.x=element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))  + theme(legend.position = "none") -> p1
p1
cowplot::ggsave2(p1, file = "../plots/grid_cd69_primary_cytokines.pdf", width = 3.7, height = 2)

fread("../data/primary-bcma-car-data.tsv") %>%
  mutate(CAR = factor(as.character(CAR), levels = c("UTD", "B3", "B1", "B4", "B2", "B5", "ABECMA"))) %>%
  mutate(Line = factor(as.character(Line), levels = rev(c("CAR_alone", "K562", "RPMI_8226", "Raji", "K562_BCMA")))) %>%
  ggplot(aes(x = CAR, y = Line, fill = CD69)) + 
  geom_tile(color = "black") + facet_grid(~Donor, scales = "free_x", space = "free_x") + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + #theme(axis.text.x=element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + theme(legend.position = "none") -> p3

cowplot::ggsave2(p3, file = "../plots/grid_cd69_donors_bcma.pdf", width = 2.6, height = 1.6)
