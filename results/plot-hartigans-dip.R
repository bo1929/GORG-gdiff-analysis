library(dplyr); library(ggplot2); library(vroom); library(latex2exp); library(tibble)
library(scales); library(ggpubr); library(cowplot); library(rstatix); library(purrr); library(tidyr)
library(tidyquant)

dipstat <- vroom("ani-comparison/gdiff-estimates-all/dipstat-gdiff-k23-w23-h9-l333-n1000-frac50.tsv.gz") %>% select(genome_a, genome_b, dip)
dist <- vroom("ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h9-l333-n1000-frac50-chisq06635-p66.tsv.gz") %>% select(genome_a, genome_b, distance)
df <- merge(dipstat, dist)

df_bi <- df %>% 
  mutate(distance_bin = cut(distance, c(0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.3, 0.35))) %>% 
  group_by(distance_bin) %>% 
  slice_max(order_by = dip, n = 1, with_ties = FALSE)

df %>%
  ggplot() +
  aes(y=dip, x=distance) +
  geom_point(alpha=0.15, size=0.5) +
  geom_point(data=df_bi, color="red", size=2) +
  theme_bw() +
  labs(x="Distance", y="Hartigans' dip statistic")

if (F) {
roll <- vroom("bimodal-pair-roll.tsv", comment = "#")
bi <- rbind(
  df_bi %>% rename(query=genome_a, reference=genome_b),
  df_bi %>% rename(query=genome_b, reference=genome_a)
)
# i <- as.integer(length(unique(dfr$seq))/6)
i <- 150
dfr <- merge(roll, bi)
s <- as.list((names(sort(table((dfr$seq)), decreasing = T))))[i:(i+5)]
dfr %>% filter(seq %in% s) %>%
  ggplot() +
  facet_wrap(~seq, scale="free") +
  aes(x=start/1e3, y=d, color=distance_bin) +
  # geom_point() +
  geom_ma(linetype=1, n = 1) + 
  theme_bw() +
  scale_color_brewer(palette="Paired", direction = -1) +
  scale_x_continuous(name="Coordinate (Kb)") +
  labs(y="Distance")
}

gt <- vroom("../anib-groundtruth.csv")
gt <- rbind(
  gt %>% rename(genome_a = `Genome A`, genome_b = `Genome B`),
  gt %>% rename(genome_b = `Genome A`, genome_a = `Genome B`)
)
nrow(df)
nrow(gt)
df <- merge(df, gt)
nrow(df)

mt <- vroom("../gorg-tropics_sags.tsv")
df <- merge(df, mt %>% rename(genome_a = SAG))
df %>%
  ggplot() +
  aes(y=dip, x=(ani_alignment_coverage_ab+ani_alignment_coverage_ba)/2*100) +
  # aes(y=dip, x=`Orthologous gene fraction, %`, color=1-`ANIgenome, %`/100) +
  geom_point(aes(color=Sample), alpha=0.15, size=0.5) +
  # geom_point(data=df %>% filter(genome_b == "AG-894-M23", genome_a=="AG-899-N21"), color="red") +
  # geom_point(data=df %>% filter(genome_b == "AG-893-E15", genome_a=="AG-899-N21"), color="red") +
  # geom_point(data=df %>% filter(genome_b == "AG-892-D10", genome_a=="AG-911-I02"), color="red") +
  # scale_color_viridis_c() +
  theme_bw() +
  labs(x="Alignment fraction (%)", y="Hartigans' dip statistic", color="D")

df %>%
  ggplot() +
  aes(y=dip, x=1-`ANIgenome, %`/100, color=(ani_alignment_coverage_ab+ani_alignment_coverage_ba)/2*100) +
  geom_point(alpha=0.15, size=0.5) +
  theme_bw() +
  scale_color_viridis_c() +
  labs(x="D", y="Hartigans' dip statistic", color="Alignment fraction (%)")
