library(BuenColors)
library(dplyr)
library(data.table)
library(ggbeeswarm)

order_cells <- c("K562.Parental", "K562.CD22",  "RPMI.8226", "Raji")
order_binders <- c("No Binder", "m971", "WT", "2261")

fread("../data/CD22-cd69-11.tsv") %>% data.frame() %>%
  reshape2::melt(id.vars = c("CAR", "Replicate")) %>%
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
  scale_fill_manual(values = c("lightgrey", "red", "darkgrey", "firebrick")) + 
  theme(legend.position = "none")

cowplot::ggsave2(p1, file = "../plots/CD22_triplicate_CD69.pdf", width = 4.5, height = 1.6)


# do pairwise test
t_test_go <- function(cell_line1, CAR1, CAR2){
  cldf <- ready_df %>% filter(variable == cell_line1) 
  vec1 <- cldf %>% filter(CAR == CAR1) %>% pull(value)
  vec2 <- cldf %>% filter(CAR == CAR2) %>% pull(value)
  (t.test(vec1, vec2))$p.value
}

t_test_go("K562.Parental", "m971", "WT")
t_test_go("K562.Parental", "WT", "2261")
t_test_go("K562.Parental", "m971", "2261")

t_test_go("K562.CD22", "m971", "WT")
t_test_go("K562.CD22", "WT", "2261")
t_test_go("K562.CD22", "m971", "2261")

t_test_go("RPMI.8226", "m971", "WT")
t_test_go("RPMI.8226", "WT", "2261")
t_test_go("RPMI.8226", "m971", "2261")

t_test_go("Raji", "m971", "WT")
t_test_go("Raji", "WT", "2261")
t_test_go("Raji", "m971", "2261")

