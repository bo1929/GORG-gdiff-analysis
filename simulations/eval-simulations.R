library(vroom)
library(rstatix)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(latex2exp)
library(cowplot)
library(fitdistrplus)

detach("package:fitdistrplus", unload = TRUE)
detach("package:MASS", unload = TRUE)
make_symmetric_max <- function(df) {
  df <- rbind(df %>% filter(genome_a > genome_b), df %>% filter(genome_b > genome_a) %>%rename(genome_a=genome_b, genome_b=genome_a)) %>%
    group_by(genome_a, genome_b) %>% summarize(ani_pct=max(ani_pct, na.rm=T))
  df <- df %>% filter(grepl("gnd", genome_a)) 
  return(df)
}
make_symmetric_mean <- function(df) {
  df <- rbind(df %>% filter(genome_a > genome_b), df %>% filter(genome_b > genome_a) %>%rename(genome_a=genome_b, genome_b=genome_a)) %>%
    group_by(genome_a, genome_b) %>% summarize(ani_pct=mean(ani_pct, na.rm=T))
  df <- df %>% filter(grepl("gnd", genome_a)) 
  return(df)
}

df_fastani <- vroom("results/distances/fastani-frag1000.tsv") %>% 
  select(genome_a, genome_b, ani_pct)
df_fastani <- make_symmetric_max(df_fastani) %>% mutate(method="fastANI")
df_mash <- vroom("results/distances/mash-k25-s10000.tsv") %>% 
  select(genome_a, genome_b, ani_pct) %>% mutate(method="mash") %>%
  filter(grepl("gnd", genome_a))
df_skani <- vroom("results/distances/skani-slow-min-af0.tsv") %>% 
  select(genome_a, genome_b, ani_pct) %>% mutate(method="skani") %>%
  filter(grepl("gnd", genome_a))
df_dashing2 <- vroom("results/distances/dashing2-v4.tsv") %>% 
  select(genome_a, genome_b, ani_pct) %>% mutate(method="dashing2") %>%
  filter(grepl("gnd", genome_a))
# df_gdiff <- vroom("results/distances/gdiff-frac50-k23-w23-h9-l500-n1000-delta3-chisq10828-p66.tsv") %>% 
df_gdiff <- vroom("results/distances/gdiff-frac50-k23-w23-h9-l333-n1000-delta3-chisq06635-p66.tsv") %>% 
  mutate(ani_pct=100-100*distance) %>% 
  select(genome_a, genome_b, ani_pct) %>%
  mutate(method="gdiff")

pairs <- df_mash %>% select(genome_a, genome_b)
dfm <- vroom("metadata.tsv") %>% rename(genome_b=seed, genome_a=key)
dfm <- merge(pairs, dfm)
dfm <- dfm  %>% 
  mutate(missing_percent = if_else(is.na(missing_percent), 0, missing_percent)) %>%
  mutate(s = if_else(is.na(s), 0, s)) %>%
  mutate(ratevar=if_else(grepl("a22", level), "low", "high")) %>%
  mutate(true_ani_bin=cut(true_ani, c(80, 90, 95, 99.9, 100))) %>%
  mutate(missing_percent_bin=cut(missing_percent, c(0, 1, 5, 10, 15, 20), include.lowest = T)) %>%
  mutate(true_ani_bin=cut(true_ani, c(80, 90, 95, 99.9, 100))) %>%
  mutate(missing_percent_bin=cut(missing_percent, c(0, 1, 5, 10, 15, 20), include.lowest = T))

df <- merge(rbind(df_mash, df_skani, df_gdiff, df_fastani, df_dashing2), dfm)
df 

# df %>% filter(missing_percent == 0) %>%
#   filter(ratevar == "high") %>%
#   ggplot() +
#   aes(x=true_ani, y=ani_pct, color=method) +
#   geom_point() +
#   theme_bw() + labs(x="True ANI", y="Estimated")
mc = c("#809D6F", "#D0D55C", "#9D4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")

df %>%
  filter(ratevar == "high") %>%
  filter(true_ani <= 99.9) %>%
  filter (missing_percent <=15) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, scales = "free_y") +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  geom_hline(yintercept = 0, linetype = 14) +
  aes(x=missing_percent_bin, y=((100-ani_pct)-(100-true_ani))/(100-true_ani)) +
  # geom_boxplot(outliers = F) +
  stat_summary(aes(color=method)) +
  stat_summary(geom="line", aes(group=method, color=method), show.legend = FALSE, position = position_dodge2()) +
  # labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Missing portion (%)", color="Method") +
  labs(y=TeX(r'(Error (%))'), x="Missing portion (%)", color="Method") +
  theme_bw() +
  scale_y_log10() +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)

df %>%
  filter(ratevar == "high") %>%
  filter(true_ani > 99.9) %>%
  filter (missing_percent == 0) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, scales = "free_y") +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  # geom_hline(yintercept = 0, linetype = 14) +
  aes(x=method, y=(ani_pct-true_ani)) +
  geom_boxplot(aes(fill=method), outliers = T) +
  # stat_summary(aes(color=method)) +
  # stat_summary(geom="line", aes(group=method, color=method), show.legend = FALSE, position = position_dodge2()) +
  # labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Missing portion (%)", color="Method") +
  labs(y=TeX(r'(${\hat{D}}-{D}$)'), fill="Method", x="Method") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)

df %>%
  filter(missing_percent == 0) %>% filter(true_ani < 99.9) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, scales = "free_y") +
  aes(x=ratevar, y=(100-ani_pct)/(100-true_ani), fill=method) +
  geom_hline(yintercept = 1, linetype = 14) +
  geom_boxplot(outliers = F) +
  labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Rate variation") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)

df %>%
  filter(missing_percent == 0) %>%
  ggplot() +
  facet_wrap(~ratevar) +
  aes(x=(100-true_ani)/100, y=(100-ani_pct)/100, color=method) +
  geom_point(alpha=0.15) +
  stat_smooth() +
  labs(y=TeX(r'(${\hat{D}}$)'), x=TeX(r'(${{D}}$)')) +
  geom_abline(linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  coord_cartesian(xlim=c(0, 0.18), ylim=c(0, 0.18))

df %>%
  filter(missing_percent == 0) %>% filter(true_ani > 99) %>%
  ggplot() +
  facet_wrap(~ratevar) +
  aes(x=(100-true_ani)/100, y=(100-ani_pct)/100, color=method, shape=method) +
  geom_point(alpha=0.1) +
  stat_summary_bin() +
  labs(y=TeX(r'(${\hat{D}}$)'), x=TeX(r'(${{D}}$)')) +
  geom_abline(linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)

df %>%
  filter(missing_percent == 0) %>% filter(true_ani > 99) %>%
  mutate(true_ani_bin = cut((100-true_ani)/100, c(0, 0.001, 0.005, 0.01), include.lowest = T)) %>%
  group_by(true_ani_bin, method, ratevar) %>%
  summarize(z=mean(abs((100-true_ani)-(100-ani_pct))/100)) %>%
  ggplot() +
  facet_wrap(~ratevar) +
  aes(fill=z, y=true_ani_bin, x=method) +
  geom_tile() +
  geom_label(aes(label=round(z, 4), color=z<0.0004), show.legend = F) +
  labs(fill=TeX(r'(${\hat{D}-D}$)'), x="Method", y=TeX(r'(D)')) +
  # geom_abline(linetype="dashed") +
  theme_cowplot(font_size = 10) +
  scale_fill_viridis_c() +
  scale_color_manual(values = c("black", "white")) +
  theme(axis.text.x.bottom = element_text(angle=0.45))
  

# Comparing different configurations of gdiff
if (F) {
df_gdiff <- merge(
  rbind(
    vroom("results/distances/gdiff-frac50-k23-w23-h9-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff-1"),
    vroom("results/distances/gdiff-frac50-k23-w23-h10-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff-2"),
    vroom("results/distances/gdiff-frac50-k23-w23-h11-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff-3"),
    vroom("results/distances/gdiff-frac33-k23-w23-h9-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff-4"),
    vroom("results/distances/gdiff-frac50-k23-w23-h9-l333-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff-5")
  )%>% 
    mutate(ani_pct=100-100*distance) %>% 
    select(genome_a, genome_b, ani_pct, method),
  dfm
)
df_gdiff %>%
  filter(ratevar == "high") %>% filter(true_ani != 100) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = 14) +
  aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
  geom_boxplot(outliers = F) +
  # stat_summary(aes(color=method)) +
  # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
  labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired")
}

# Comparing different configurations of skani
if (F) {
  df_skani <- merge(
    rbind(
      vroom("results/distances/skani-c30-m300-slow.tsv") %>%  mutate(method="skani-a"), 
      vroom("results/distances/skani-robust-min-af0-c70.tsv") %>%  mutate(method="skani-b"),
      vroom("results/distances/skani-slow-min-af0.tsv") %>%  mutate(method="skani-c")
    )%>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_skani %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    ggplot() +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = 14) +
    aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
    theme_bw() +
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired")
}

# Comparing different configurations of mash
if (F) {
  df_mash <- merge(
    rbind(
      vroom("results/distances/mash-k19-s10000.tsv") %>%  mutate(method="mash-k19-s10000"), 
      vroom("results/distances/mash-k21-s1000.tsv") %>%  mutate(method="mash-k21-s1000"), 
      vroom("results/distances/mash-k21-s10000.tsv") %>%  mutate(method="mash-k21-s10000"),
      vroom("results/distances/mash-k23-s10000.tsv") %>%  mutate(method="mash-k23-s10000")
    )%>% 
      mutate(ani_pct=100-100*distance) %>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_mash %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    ggplot() +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = 14) +
    aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
    theme_bw() +
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired")
}

# Comparing different configurations of fastani
if (F) {
  df_fastani <- merge(
    rbind(
      vroom("results/distances/fastani-frag1000.tsv") %>%  mutate(method="fastani-frag1000"), 
      vroom("results/distances/fastani-frag3000.tsv") %>%  mutate(method="fastani-frag3000")
    )%>% 
      # mutate(ani_pct=100-100*distance) %>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_fastani %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    ggplot() +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = 14) +
    aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
    theme_bw() +
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired")
}

# Comparing different configurations of dashing2
if (F) {
  df_dashing2 <- merge(
    rbind(
      vroom("results/distances/dashing2-v4.tsv") %>%  mutate(method="dashing2-v4"), 
      vroom("results/distances/dashing2-default.tsv") %>%  mutate(method="dashing2-default"), 
      vroom("results/distances/dashing2-containment.tsv") %>%  mutate(method="dashing2-containment"), 
      vroom("results/distances/dashing2-symmetric.tsv") %>%  mutate(method="dashing2-symmetric")
    )%>% 
      mutate(ani_pct=100-100*distance) %>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_dashing2 %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    ggplot() +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = 14) +
    aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
    theme_bw() +
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired")
}

# See cv_analysis for a better coefficient of variance analysis:
if (F) {
  dfc <- vroom("results-ORFvANI/results-combined.csv") %>% filter(grepl("_p0_s0", genome_b))
  dfc <- dfc %>% mutate(genome_b=gsub(genome_b, pattern="_gnd00000_.*", replacement = ""))
  dfc <- dfc %>% mutate(genome_a=gsub(genome_a, pattern="_c.*", replacement=""))
  dfcm <- merge(
    dfm %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")),
    rbind(
      # df_gdiff_mean %>% select(genome_a, genome_b, ani_pct, method) %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")),
      df_gdiff %>% select(genome_a, genome_b, ani_pct, method) %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")),
      dfc %>% mutate(method="BLAST") %>% rename(ani_pct=mean_pident) %>% select(genome_b, genome_a, ani_pct, method),
      df_skani %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")),
      df_fastani %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")),
      # df_gdiff_median %>% select(genome_a, genome_b, ani_pct, method) %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")),
      df_mash %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement=""))
    ), by=c("genome_a", "genome_b")
  ) %>%  mutate(missing_percent = if_else(is.na(missing_percent), 0, missing_percent)) %>%
    mutate(s = if_else(is.na(s), 0, s)) 
  
  dfcm %>% filter(grepl("_p", genome_a)) %>%
    ggplot() +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
    # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
    geom_hline(yintercept = 1, linetype = 14) +
    aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
    theme_bw()
  
  merge(
    rbind(
      dfc %>% select(genome_a, genome_b, stdev_pident, mean_pident) %>% mutate(method="BLAST") %>% rename(ani_pct=mean_pident),
      dfs_gdiff %>% filter(grepl("a5", genome_a)) %>% mutate(genome_a=gsub(genome_a, pattern="_a5", replacement="")) %>% group_by(genome_a, genome_b) %>% summarise(ani_pct=100-100*mean(d, na.rm=T), stdev_pident=sd(100-100*d, na.rm=T)) %>% mutate(method="gdiff")),
    dfm %>% filter(grepl("a5", genome_a)) %>%mutate(genome_a=gsub(genome_a, pattern="_a5", replacement=""))
  ) %>%
    filter(grepl("_p", genome_a)) %>% filter(true_ani !=100) %>%
    # filter(missing_percent == 0) %>% # select(true_ani)
    ggplot() +
    aes(x=true_ani, y=stdev_pident/(100-ani_pct)*100, color=method) +
    # stat_summary()+
    geom_point(alpha=0.4)+
    stat_smooth()+
    facet_wrap(~missing_percent_bin) +
    # stat_ecdf() + 
    # geom_abline()+
    geom_hline(yintercept = 1/sqrt(5)*100)+
    theme_bw() + labs(x="ANI", y="Coefficient of variance", title="Without outlier removal")+
    coord_cartesian(ylim=c(0, 100))
}

# Artifacts from sample inspection:
if (F) {
  reconcile_directional <- function(df) {
    df <- df %>%
      mutate(
        X = pmax(genome_a, genome_b),
        Y = pmin(genome_a, genome_b),
        direction = if_else(genome_a == X, "XY", "YX")
      )
    
    keep_cols <- setdiff(names(df), c("genome_a", "genome_b", "X", "Y", "direction"))
    
    result <- df %>%
      group_by(config, X, Y) %>%
      group_modify(~ {
        grp <- .x
        
        xy <- grp %>% filter(direction == "XY") %>% arrange(d) %>% select(all_of(keep_cols[keep_cols %in% names(grp)]))
        yx <- grp %>% filter(direction == "YX") %>% arrange(d) %>% select(all_of(keep_cols[keep_cols %in% names(grp)]))
        
        n_x <- nrow(xy)
        n_y <- nrow(yx)
        n_max <- max(n_x, n_y, 1)
        
        pad <- function(block, n) {
          if (nrow(block) < n) {
            extra <- n - nrow(block)
            na_rows <- block[rep(NA_integer_, extra), , drop = FALSE]
            block <- bind_rows(block, na_rows)
          }
          block
        }
        
        xy_p <- pad(xy, n_max)
        yx_p <- pad(yx, n_max)
        
        chosen <- map_dfr(seq_len(n_max), function(i) {
          row_xy <- xy_p[i, , drop = FALSE]
          row_yx <- yx_p[i, , drop = FALSE]
          d_xy <- row_xy$d
          d_yx <- row_yx$d
          
          if (is.na(d_xy) && is.na(d_yx)) {
            out <- row_xy
            out[] <- NA
          } else if (is.na(d_xy)) {
            out <- row_yx
          } else if (is.na(d_yx)) {
            out <- row_xy
          } else if (d_xy <= d_yx) {
            out <- row_xy
          } else {
            out <- row_yx
          }
          out
        })
        
        chosen$N_X <- n_x
        chosen$N_Y <- n_y
        chosen
      }) %>%
      ungroup() %>%
      rename(genome_a = X, genome_b = Y) %>%
      relocate(config, genome_a, genome_b, N_X, N_Y)
    
    result
  }
  dfs_gdiff <- vroom("results/samples/gdiff-frac50-k23-w23-l500-n1000-delta3.tsv")
  # dfs_gdiff <- reconcile_directional(dfs_gdiff)
  df_gdiff_mean <- dfs_gdiff %>%
    group_by(genome_a, genome_b) %>%
    summarise(n=n(), r=sum(!is.na(n))/n(), ani_pct=100-100*mean(d, na.rm = T))
  df_gdiff_median <- dfs_gdiff %>%
    group_by(genome_a, genome_b) %>%
    summarise(n=n(), r=sum(!is.na(n))/n(), ani_pct=100-100*median(d, na.rm = T))
  
  df <- merge(
    rbind(
      df_gdiff,
      df_mash,
      df_skani,
      df_fastani,
      df_gdiff,
      make_symmetric_mean(df_gdiff_mean) %>% select(genome_a, genome_b, ani_pct)  %>% mutate(method="gdiff-mean-pair-avg"),
      make_symmetric_mean(df_gdiff_median) %>% select(genome_a, genome_b, ani_pct) %>% mutate(method="gdiff-median-pair-avg"),
      make_symmetric_max(df_gdiff_mean) %>% select(genome_a, genome_b, ani_pct)  %>% mutate(method="gdiff-mean-pair-max"),
      make_symmetric_max(df_gdiff_median) %>% select(genome_a, genome_b, ani_pct) %>% mutate(method="gdiff-median-pair-max")
    ),
    dfm
  )
  df  
  
  df %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    filter(grepl("gdiff", method)) %>%
    ggplot() +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
    # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
    geom_hline(yintercept = 1, linetype = 14) +
    aes(x=missing_percent_bin, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    labs(y=TeX(r'(${\hat{ANI}}/{ANI}$)'), x="Missing portion (%)") +
    theme_bw()
  
  df %>% filter(grepl("gdiff", method)) %>%
    filter(missing_percent == 0) %>% filter(true_ani != 100) %>%
    ggplot() +
    geom_hline(yintercept = 1, linetype = 14) +
    facet_wrap(~true_ani_bin, scales = "free_y") +
    aes(x=ratevar, y=(ani_pct)/(true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    labs(y=TeX(r'(${\hat{GND}}/{GND}$)'), x="Rate variation") +
    theme_bw()
  
  # Couple of examples distributions for debugging:
  df_wi <- df %>%
    filter(missing_percent == 0) %>% filter(ratevar == "low") %>% filter(true_ani != 100) %>%
    group_by(true_ani_bin, method) %>%
    slice_max(abs(1-(ani_pct)/(true_ani)), n = 2) %>%
    filter(method %in% c("gdiff"))#  %>% rename(genome_a=genome_b, genome_b=genome_a)
  genome_wia <- df_wi$genome_a
  genome_wib <- df_wi$genome_b
  df_inspect <- merge(dfs_gdiff %>% filter((genome_a %in% genome_wia) & (genome_b %in% genome_wib)), df_wi)
  df_inspect %>% mutate(d=if_else(!is.finite(d), 0.5, d)) %>%
    ggplot() +
    aes(x=d) +
    facet_wrap(~true_ani, scale="free") +
    stat_bin() +
    geom_vline(aes(xintercept=(100-ani_pct)/100), color="blue", label="Estimated GND") +
    geom_vline(aes(xintercept=(100-true_ani)/100), color="red", label="True GND")

}