library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(rafalib)
library(Rtsne)
library(umap)
library(hash)
library(R.matlab)


theme_set(theme_gray() +
            theme(
              axis.line = element_line(size=0.5),
              panel.background = element_rect(fill=NA,size=rel(20)),
              panel.grid.minor = element_line(colour = NA),
              axis.text = element_text(size=16),
              axis.title = element_text(size=18)
            )
)

library(dplyr)
library(viridis)

data <- read.csv("./input_of_R/forebrain_half_dropout_cell_sum.csv")
pdf(file=sprintf("../results/plot/0503/supp/forebrain_do_cell_sum.pdf"), width=12, height=6)
# sample size
sample_size = data %>% group_by(do_rate) %>% summarize(num=n())

# Plot
data %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(do_rate, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=cell_sum, fill=do_rate)) +
    geom_violin(width=1) +
    geom_boxplot(width=0.15, color="grey", alpha=0.7) +
    scale_fill_viridis(discrete = TRUE) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("forebrain_half_dropout_cell_sum") +
    xlab("")

dev.off()
