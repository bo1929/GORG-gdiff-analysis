library(vroom)
library(rstatix)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(latex2exp)
library(cowplot)
library(ggpubr)
library(vroom)

df <- vroom("blastn-comparison/genes_mapped.tsv")
names(df)[1] <- "Genome A"
names(df)[2] <- "Genome B"

df_ani <- df_ani <- vroom("../anib-groundtruth.csv")
df_ani <- rbind(df_ani %>% filter(`Genome A` < `Genome B`), df_ani %>% filter(`Genome A` >= `Genome B`) %>% rename(`Genome A`=`Genome B`, `Genome B`=`Genome A`))

merge(df, df_ani) %>% 
  filter(blastn_aln_len > 200) %>%
  ggplot() +
  aes(x=blastn_d, y=gdiff_d, color=`ANIgenome, %`) +
  geom_point(alpha=0.25) +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.y = 0.29) +
  geom_abline(color="red") +
  scale_color_viridis_c(end = 1) +
  # geom_text(
  #  data = function(d) d %>% 
  #    filter(is.finite(blastn_d), is.finite(gdiff_d)) %>%
  #    summarize(mae = mean(abs(((blastn_d) - (gdiff_d)))),  .groups = "drop") %>% 
  #    mutate(x = 0.0, y = 0.25),
  #  aes(x=x, y=y, label=paste0("MAE = ", round(mae, 2), "")), hjust=0, vjust=0, show.legend = F, inherit.aes = F
  # ) +
  theme_minimal_grid(font_size = 9) +
  panel_border(color = "black", size = 0.5) +
  coord_cartesian(xlim=c(0, 0.3), ylim=c(0, 0.3)) +
  # labs(x="BLASTn", y=TeX(r'(gdiff $D$)'), color="ANI", title="Distances between annotated genes in marine microbes")
  labs(x="BLASTn", y="gdiff", color="ANI", title="Distances between annotated genes")
ggsave("./gene-dist-compare.pdf", width=3.75, height = 3.45)
