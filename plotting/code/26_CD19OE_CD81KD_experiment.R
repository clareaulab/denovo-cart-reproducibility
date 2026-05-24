library(BuenColors)
library(dplyr)
library(data.table)
library(cowplot)

pair_table = read.csv("../data/CD19_CD81_KD_K562_CD69_coculturevalues.csv")

pair_table$cd19_oe_gain = pair_table$K562_CD19OE-pair_table$CAR.only
pair_table$cd19oe_cd81ko_gain = pair_table$K562_CD19OE_CD81_KD-pair_table$CAR.only

pair_table = pair_table %>% mutate(Binder_name = case_when(
  Binder_name == "l112" ~ "C1",
  Binder_name == "l115" ~ "C2",
  Binder_name == "l115_D2" ~ "C2.N2",
  TRUE ~ Binder_name
))

pair_table_long = pair_table %>% 
  pivot_longer(
    cols = c(CAR.only,K562_CD19OE, K562_CD19OE_CD81_KD), 
    names_to = "Condition", 
    values_to = "Value"
  )
pair_table_long

pair_table_long_renamed = pair_table_long %>% mutate(
  Condition = case_when(
    Condition == "CAR.only" ~ "CAR alone",
    Condition == "K562_CD19OE" ~ "K562\n(CD19 OE)",
    Condition == "K562_CD19OE_CD81_KD" ~ "K562\n(CD19 OE, CD81 KD)",
    TRUE ~ Condition
  )
)

pair_table_long_renamed$Binder_name = factor(
  pair_table_long_renamed$Binder_name,
  levels = c(
    "UTD","NB","FMC63","C1","C2","C2.N2"
  )
)
pair_table_long_renamed$Condition = factor(
  pair_table_long_renamed$Condition,
  levels = c(
    "CAR alone","K562\n(CD19 OE)","K562\n(CD19 OE, CD81 KD)"
  )
)

pair_table_palette = c(
  "UTD"="darkgray",
  "NB"="gray",
  "FMC63"="red",
  "C1"="dodgerblue",
  "C2"="dodgerblue2",
  "C2.N2"="dodgerblue4"
)

pair_table_long_renamed_plot = ggplot(pair_table_long_renamed, aes(x=Condition,y=Value,group=Binder_name)) +
  geom_line(size = 0.8,aes(color=Binder_name)) +
  geom_point(aes(color = Binder_name), size = 3) + 
  pretty_plot() + L_border() +
  scale_color_manual(values=pair_table_palette) +
  labs(y="%CD69+ Cells",color="Binder Name")

pair_table_long_renamed_plot

## Plot gains only
pair_table_gain_long = pair_table %>% 
  pivot_longer(
    cols = c(cd19_oe_gain, cd19oe_cd81ko_gain), 
    names_to = "Condition", 
    values_to = "Value"
  )

pair_table_gain_long

pair_table_gain_long_renamed = pair_table_gain_long %>% mutate(
  Condition = case_when(
    Condition == "cd19_oe_gain" ~ "K562 CD19 OE",
    Condition == "cd19oe_cd81ko_gain" ~ "K562 CD19 OE\n(CD81 KD)",
    TRUE ~ Condition
  )
)

pair_table_gain_long_renamed$Binder_name = factor(
  pair_table_gain_long_renamed$Binder_name,
  levels = c(
    "UTD","NB","FMC63","C1","C2","C2.N2"
  ))

pair_table_gain_long_renamed_plot = ggplot(pair_table_gain_long_renamed, aes(x=Condition,y=Value,group=Binder_name)) +
  geom_line(size = 0.8,aes(color=Binder_name)) +
  geom_point(aes(color = Binder_name), size = 3) + 
  pretty_plot() + L_border() +
  scale_color_manual(values=pair_table_palette) +
  labs(y="Activity Gain",color="Binder Name")

pair_table_gain_long_renamed_plot
cowplot::ggsave2("../plots/CD19OE_CD81KD_pairplot.png",pair_table_gain_long_renamed_plot,dpi=300)

## Another way of showing the same data



