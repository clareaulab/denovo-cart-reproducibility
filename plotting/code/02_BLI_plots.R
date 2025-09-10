library(BuenColors)
library(dplyr)
library(readxl)

raw_in <- read_xlsx("../data/BLI_BCMA_minibinders_raw.xlsx") %>% data.frame()

melt_df <- data.frame(
  time = raw_in[[1]],
  binder = c(raw_in[[6]], raw_in[[11]], raw_in[[15]], raw_in[[20]]),
  a1nM = c(raw_in[[2]], raw_in[[7]], raw_in[[12]], raw_in[[16]]),
  b10nM = c(raw_in[[3]], raw_in[[8]], raw_in[[13]], raw_in[[17]]),
  c100nM = c(raw_in[[4]], raw_in[[9]], raw_in[[14]], raw_in[[18]]),
  d1000nM = c(raw_in[[5]], raw_in[[10]], rep(NA, length(raw_in[[14]])), raw_in[[19]])
) %>%
  reshape2::melt(id.vars = c("binder", "time")) 

px <- melt_df %>% mutate(variable = factor(as.character(variable), rev(c("d1000nM", "c100nM", "b10nM", "a1nM")))) %>%
  arrange(desc(variable)) %>%
  filter(value > -0.02) %>%
  filter(binder != "B4_903444") %>%
  filter(binder != "B3_961593") %>%
  ggplot(., aes(x = time, y = value, color= variable)) +
  geom_point(size = 0.3) + pretty_plot(fontsize = 8) +
  facet_wrap(~binder)  + scale_color_viridis_d()
cowplot::ggsave2(px, file = "../plots/fig1e_BLI_two.pdf", width = 3.5, height = 1.4)
