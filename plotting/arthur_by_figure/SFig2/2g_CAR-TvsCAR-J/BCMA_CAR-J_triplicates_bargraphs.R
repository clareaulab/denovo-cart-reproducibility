library(BuenColors)
library(dplyr)
library(data.table)
library(tidyr)
library(ggbeeswarm)

order_cells <- c("K562.Parental", "K562.BCMA",  "RPMI.8226", "Raji")
order_binders <- c("No Binder", "Abecma", "561726_original")

# Find all columns EXCEPT "CAR" and "Replicate" and pivot them
fread("../data/BCMA_CAR-J_CD69_triplicates.tsv") %>% data.frame() %>%
  tidyr::pivot_longer(
    cols = !c("CAR", "Replicate"),  # <--- CORRECTED THIS LINE
    names_to = "variable",          # Sets the name for the column that holds the old column names
    values_to = "value"             # Sets the name for the column that holds the data values
  ) %>%
  filter(variable %in% order_cells) %>%
  filter(CAR %in% order_binders) %>% 
  mutate(CAR = factor(as.character(CAR), levels = order_binders)) %>% 
  mutate(variable = factor(as.character(variable), levels = order_cells)) -> ready_df

mean_df<- ready_df %>% group_by(CAR, variable) %>% summarize(value = mean(value))

p1 <- ready_df %>%
  ggplot(aes(x = CAR, y = value, fill = CAR)) + 
  geom_bar(data = mean_df, stat = "identity", color = 'black') + 
  geom_quasirandom() + facet_wrap(~variable, nrow = 1) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + labs(x = "", y = "") +
  scale_fill_manual(values = c("black", "#ff0000", "gray", "#84C8FF")) + 
  theme(legend.position = "none")
p1
cowplot::ggsave2(p1, file = "../plots/BCMA_parental_CD69_CAR-J_triplicates.pdf", width = 4.5, height = 1.6)


# do pairwise test
t_test_go <- function(cell_line1, CAR1, CAR2){
  cldf <- ready_df %>% filter(variable == cell_line1) 
  vec1 <- cldf %>% filter(CAR == CAR1) %>% pull(value)
  vec2 <- cldf %>% filter(CAR == CAR2) %>% pull(value)
  (t.test(vec1, vec2))$p.value
}

t_test_go("K562.Parental", "No Binder", "Abecma")
t_test_go("K562.Parental", "No Binder", "WT")
t_test_go("K562.Parental", "WT", "4434")

t_test_go("K562.CD22", "m971", "WT")
t_test_go("K562.CD22", "m971", "4434")
t_test_go("K562.CD22", "WT", "4434")

t_test_go("RPMI.8226", "m971", "WT")
t_test_go("RPMI.8226", "m971", "4434")
t_test_go("RPMI.8226", "WT", "4434")

t_test_go("Raji", "m971", "WT")
t_test_go("Raji", "WT", "4434")
t_test_go("Raji", "m971", "4434")

p1

