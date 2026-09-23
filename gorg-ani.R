library(dplyr); library(ggplot2); library(vroom); library(latex2exp)
library(scales)
library(ggpubr)

mc = c("#809D6F", "#D0D55C", "#9D4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")
colm = c("genome_a", "genome_b", "ani_true")
df_meta <- read.delim("all_pairs.tsv", comment.char = "#", header = FALSE, col.names=colm, colClasses = c("character", "character", "numeric")) |>
  mutate(
    pair = paste(pmin(genome_a, genome_b), pmax(genome_a, genome_b), sep = "|"),
    # ani_bin = cut(ani_true, breaks = c(55, 60, 65, 70, 75, 80, 85, 90, 95, 97.5, 99.9, 100), include.lowest = TRUE)
    ani_bin = cut(ani_true, breaks = c(55, 60, 80, 90, 95, 99, 100), include.lowest = TRUE)
  )
n_pairs <- nrow(df_meta)
df_meta <- df_meta %>% 
  mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
  mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
  select(!c(genome_a, genome_b)) %>%
  rename(genome_a = genome_x, genome_b = genome_y)
df_r <- rbind(
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v4.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100),
  vroom("results/ani-comparison/mash-estimates-all/distances/mash-sensitive.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct, distance),
  vroom("results/ani-comparison/skani-estimates-all/distances/skani-slow.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100),
  vroom("results/ani-comparison/fastani-estimates-all/distances/fastani-frag1000.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100),
  vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h9-l500-n1000-frac50-chisq10828-p66.tsv") %>%
  # vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h9-l333-n1000-frac50-chisq06635-p66.tsv") %>%
    select(genome_a, genome_b, distance) %>% mutate(method = "gdiff", ani_est = (1-distance)*100)
) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
  mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
  select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
df_r <- df_r %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
df <- merge(df_r, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)

df %>% filter(ani_true < 99) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.15, 0.25, 0.4), include.lowest = TRUE)) %>%
  ggplot() +
  # aes(x=ani_bin, fill=method, y=d_est/d_true, color=method) +
  aes(x=ani_bin, fill=method, y=(100-ani_est)/(100-ani_true)-1, color=method) +
  # aes(x=ani_bin, fill=method, y=((100-ani_est)-(100-ani_true))/(100-ani_true), color=method) +
  geom_boxplot(outliers = F, outlier.size = 0.6, position = position_dodge(0.8), color = "grey10") +
  # stat_summary(aes(group=method), geom="line") +
  # stat_summary(aes(group=method)) +
  # facet_wrap(~ani_bin, scale="free") +
  geom_hline(yintercept = 0, linetype=2) +
  # coord_cartesian(ylim = c(0.95, 1.15)) +
  # scale_y_log10() +
  theme_bw() + labs(x="ANIb", y=TeX(r'(Error (%))'), title="251,534 pairs from GORG") +
  scale_fill_manual(values = mc) +
  scale_color_manual(values = mc) + scale_y_continuous(labels=percent)

df %>%
  filter(ani_true > 99) %>%
  mutate(ani_bin = cut((100-ani_true)/100, c(0, 0.001, 0.005, 0.01), include.lowest = T)) %>%
  group_by(ani_bin, method) %>%
  summarize(z=mean(abs((100-ani_true)-(100-ani_est))/100)) %>%
  ggplot() +
  aes(fill=z, y=ani_bin, x=method) +
  geom_tile() +
  geom_label(aes(label=round(z, 4), color=z<0.005), show.legend = F) +
  labs(fill=TeX(r'(${\hat{D}-D_{ANIb}}$)'), x="Method", y=TeX(r'(${D_{ANIb}$})')) +
  # geom_abline(linetype="dashed") +
  theme_cowplot(font_size = 10) +
  scale_fill_viridis_c(option = "B", limits = c(0, 0.01)) +
  scale_color_manual(values = c("black", "white")) +
  theme(axis.text.x.bottom = element_text(angle=0.45), legend.title = element_text(vjust=2))

df %>% filter(ani_true < 99) %>%
  # mutate(ani_bin = cut(ani_true, c(60, 75, 90, 95, 99))) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.15, 0.25, 0.4), include.lowest = TRUE)) %>%
  ggplot() +
  # aes(x=ani_bin, fill=method, y=d_est/d_true, color=method) +
  # aes(x=ani_bin, fill=method, y=(100-ani_est)/(100-ani_true), color=method) +
  aes(x=method, fill=method, y=((100-ani_est)-(100-ani_true))/(100-ani_true), color=method) +
  geom_violin(outliers = F, trim = T, scale = "width", outlier.size = 0.2, color = "grey10") +
  # stat_summary(aes(group=method), geom="line") +
  stat_summary(color="black") +
  facet_wrap(~ani_bin, scale="free", nrow=1) +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim = c(0.95, 1.15)) +
  # scale_y_log10() +
  theme_bw() + labs(x="ANIb", y=TeX(r'(Error)'), title="251,534 pairs from GORG") +
  scale_color_manual(values = mc) + scale_y_continuous(labels=percent) +
  scale_fill_manual(values = mc) + theme(axis.text.x = element_blank()) +
  labs(x="", fill="color") + coord_cartesian(ylim = c(-0.60, 0.60)) 

ggscatter(
  df %>% filter(ani_true > 90),
  x = "ani_true", y = "ani_est", color = "method",
  add = "reg.line", alpha=0.05, conf.int = TRUE
  ) +
  stat_cor(aes(color = method), method = "pearson", show.legend = F) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  geom_abline(linewidth=2, alpha=1, linetype="dashed") +
  coord_cartesian(x=c(90, 100), y=c(90, 100)) +
  labs(x="ANIb", y=TeX(r'(${\hat{ANI}}$)'))

df %>% filter(ani_true > 90) %>% ggplot() +
  aes(x=ani_true, y=ani_est, color=method) +
  stat_cor(method = "pearson", show.legend = F) +
  geom_point(alpha=0.05) +
  stat_smooth(method = "lm", linewidth=1.5) +
  geom_abline(linewidth=1, alpha=1, linetype="dashed") +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  theme_bw() +
  labs(x="ANIb", y=TeX(r'(${\hat{ANI}}$)'))

df %>% filter(ani_true < 99) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4), include.lowest = TRUE)) %>%
  ggplot() +
  aes(x=ani_bin, y=(100-ani_est)/(100-ani_true)-1, color=method) +
  stat_summary(aes(group=method), geom="line") +
  stat_summary(aes(group=method)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() +
  labs(x=TeX(r'($D_{ANIb})'), y=TeX(r'(Error (%))'), title="251,534 pairs from GORG", color="Method") +
  scale_fill_manual(values = mc) +
  scale_color_manual(values = mc)

merge(
  df %>% mutate(ani_bin=cut(ani_true, c(65, 70, 75, 80, 100))) %>% group_by(method, ani_bin) %>% summarise(c=n()),
  df_meta %>% mutate(ani_bin=cut(ani_true, c(65, 70, 75, 80, 100))) %>% group_by(ani_bin) %>% summarise(t=n()), all.y=T
) %>% complete(method, ani_bin, fill = list(c = 0, t = 1)) %>%
  ggplot() +
  aes(x=ani_bin, fill=method, y=c/t) +
  geom_col(color="gray10", position = position_dodge2(width = 0.8, padding = 0.15), na.rm = F) +
  # geom_line(aes(group=method, color=method)) +
  theme_bw() + labs(x="ANIb", y="Pairs with an estimate (%)", title="251,534 pairs from GORG") +
  scale_fill_manual(values = mc) +
  # scale_color_manual(values = mc) +
  scale_y_continuous(labels = percent)

df %>%
  group_by(ani_true, method) %>% summarise(nx=n()) %>%
  group_by(method) %>%
  arrange(ani_true) %>%
  mutate(nt = cumsum(nx)/n_pairs) %>%
  mutate(ani_true=round(ani_true, 1)) %>%
  group_by(ani_true, method) %>%
  summarize(nt=mean(nt)) %>%
  ggplot() +
  aes(x=ani_true, y=nt, color=method) +
  geom_line(aes(group=method, linetype = method), linewidth=2) +
  # stat_smooth(aes(group=method), m) +
  theme_bw() + labs(x="ANIb", y="Pairs with an estimate", title="251,534 pairs from GORG") +
  scale_color_manual(values = mc) +
  scale_linetype_manual(values = c(1, 1, 2, 3, 1)) +
  scale_y_continuous(labels = percent)

# dashing2 configuration comparison
if (T) {
df_dashing <- rbind(
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v1.tsv") %>% mutate(method="v1"),
  vroom("../results/ani-comparison/dashing2-estimates-all/distances/dashing2-v2.tsv") %>% mutate(method="v2"),
  vroom("../results/ani-comparison/dashing2-estimates-all/distances/dashing2-v3.tsv") %>% mutate(method="v3"),
  vroom("../results/ani-comparison/dashing2-estimates-all/distances/dashing2-v4.tsv") %>% mutate(method="v4"),
  vroom("../results/ani-comparison/dashing2-estimates-all/distances/dashing2-v5.tsv") %>% mutate(method="v5")
) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
  mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
  select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
df_dashing <- df_dashing %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
df_dashing <- merge(df_dashing, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
df_dashing %>%
  ggplot() +
  aes(x=ani_bin, fill=method, y=ani_est/ani_true, color=method) +
  stat_summary(aes(group=method), geom="line") +
  stat_summary() +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(0.9, 1.2)) +
  theme_bw() + labs(x="ANIb", y=TeX(r'(${\hat{ANI}}/{{ANIb}}$)'), title="251,534 pairs from GORG") +
  scale_color_brewer(palette="Set1")
}

# mash configuration comparison
if (T) {
  df_mash <- rbind(
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-sensitive.tsv.gz") %>% mutate(method="sensitive"),
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-long-k.tsv.gz") %>% mutate(method="long-k"),
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-large-sketch.tsv.gz") %>% mutate(method="large-sketch"),
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-default.tsv.gz") %>% mutate(method="default")
  ) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_mash <- df_mash %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_mash <- merge(df_mash, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_mash %>%
    ggplot() +
    aes(x=ani_bin, fill=method, y=ani_est/ani_true, color=method) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 1) +
    coord_cartesian(ylim = c(0.9, 1.2)) +
    theme_bw() + labs(x="ANIb", y=TeX(r'(${\hat{ANI}}/{{ANIb}}$)'), title="251,534 pairs from GORG") +
    scale_color_brewer(palette="Set1")
}

# skani configuration comparison
if (T) {
  df_skani <- rbind(
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-slow.tsv.gz") %>% mutate(method="slow"),
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-sensitive.tsv.gz") %>% mutate(method="sensitive"),
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-fast.tsv.gz") %>% mutate(method="fast"),
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-default.tsv.gz") %>% mutate(method="default")
  ) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_skani <- df_skani %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_skani <- merge(df_skani, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_skani %>%
    ggplot() +
    aes(x=ani_bin, fill=method, y=ani_est/ani_true, color=method) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 1) +
    coord_cartesian(ylim = c(0.9, 1.2)) +
    theme_bw() + labs(x="ANIb", y=TeX(r'(${\hat{ANI}}/{{ANIb}}$)'), title="251,534 pairs from GORG") +
    scale_color_brewer(palette="Set1")
}

# fastani configuration comparison
if (T) {
  df_fastani <- rbind(
    vroom("results/ani-comparison/fastani-estimates-all/distances/fastani-frag1000.tsv.gz") %>% mutate(method="--fragLen 1000 --minFraction 0.1"),
    vroom("results/ani-comparison/fastani-estimates-all/distances/fastani-frag3000.tsv.gz") %>% mutate(method="--fragLen 3000 --minFraction 0.1")
  ) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100)  %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_fastani <- df_fastani %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_fastani <- merge(df_fastani, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_fastani %>%
    ggplot() +
    aes(x=ani_bin, fill=method, y=ani_est/ani_true, color=method) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 1) +
    coord_cartesian(ylim = c(0.9, 1.2)) +
    theme_bw() + labs(x="ANIb", y=TeX(r'(${\hat{ANI}}/{{ANIb}}$)'), title="251,534 pairs from GORG") +
    scale_color_brewer(palette="Set1")
}

# gdiff configuration comparison
if (T) {
  df_gdiff <- rbind(
    vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h10-l500-n1000-frac50-chisq06635-p66.tsv") %>% mutate(method="k23-w23-h10-l500-n1000-frac50"),
    vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h11-l500-n1000-frac50-chisq03841-p66.tsv") %>% mutate(method="k23-w23-h11-l500-n1000-frac50"),
    vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h9-l333-n1000-frac50-chisq06635-p66.tsv") %>% mutate(method="k23-w23-h9-l333-n1000-frac50"),
    vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h9-l500-n1000-frac33-chisq10828-p66.tsv") %>% mutate(method="k23-w23-h9-l500-n1000-frac33"),
    vroom("results/ani-comparison/gdiff-estimates/distances/gdiff-k23-w23-h9-l500-n1000-frac50-chisq10828-p66.tsv") %>% mutate(method="k23-w23-h9-l500-n1000-frac50")
    ) %>% select(method, genome_a, genome_b, distance) %>% mutate(ani_est = 100-distance*100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_gdiff <- df_gdiff %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_gdiff <- merge(df_gdiff, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_gdiff %>% filter(ani_true < 99) %>%
    ggplot() +
    aes(x=ani_bin, fill=method, y=(100-ani_est)/(100-ani_true)-1) +
    # facet_wrap(~ani_bin, scale="free") +
    stat_summary(aes(group=method, color=method), geom="line") +
    # geom_boxplot(color="gray20", outliers = F) +
    stat_summary(aes(color=method), fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 0) +
    # coord_cartesian(ylim = c(0.9, 1.2)) +
    theme_bw() +
    labs(x="ANIb", y=TeX(r'(${\hat{ANI}}/{{ANIb}}$)'), title="251,534 pairs from GORG") +
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1")
}
