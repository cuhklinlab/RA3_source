# ~/anaconda3/bin/R
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(rafalib)
library(Rtsne)
library(umap)
library(hash)
library(gridExtra)

theme_set(theme_gray() +
            theme(
              axis.line = element_line(size=0.5),
              panel.background = element_rect(fill=NA,size=rel(20)),
              panel.grid.minor = element_line(colour = NA),
              axis.text = element_text(size=16),
              axis.title = element_text(size=18)
            )
)

res_table <- read.csv("./input_of_R/time.csv")
# head(res_table)
# res_table$peak.selection.k <- as.character(res_table$peak.selection.k)
# res_table$H1drop_pct <- as.character(res_table$H1drop_pct)
# res_table$K3 <- as.character(res_table$K3)

all_methods = c("scABC","SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC","RA3")


data <- which(res_table$DR.methods%in%all_methods)
data <- res_table[data,]
data$DR.methods <- factor(data$DR.methods, levels=all_methods)
data$Dataset <- factor(data$Dataset, levels=c("MPP_LMPP_CLP","donor_BM0828","forebrain_half","ALL_blood","tissue4_forebrain", "sc_simulation_tissue4_forebrain_p81190_n10000","sc_simulation_tissue4_forebrain_p81190_n20000","sc_simulation_tissue4_forebrain_p81190_n50000"))
data$Dim.Reduction.Time <- data$Dim.Reduction.Time.ms / 1000
data$Dim.Reduction.Time.log <- log10(data$Dim.Reduction.Time)

data$Dim.Reduction.Time.log[which(data$Dim.Reduction.Time.ms==0)] <- 0

pdf(file=sprintf("../results/plot/1021/Q2_time_bar.pdf"), width=20, height=8)
data1 <- data[which(data$Dataset%in%c("MPP_LMPP_CLP","donor_BM0828","forebrain_half","ALL_blood")),]
res1 <- (ggplot(data1) + geom_bar(aes(x=Dataset, y=Dim.Reduction.Time.log, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(-0.5,3.5))  )

data2 <- data[which(data$Dataset%in%c("tissue4_forebrain", "sc_simulation_tissue4_forebrain_p81190_n10000","sc_simulation_tissue4_forebrain_p81190_n20000","sc_simulation_tissue4_forebrain_p81190_n50000")),]
res2 <- (ggplot(data2) + geom_bar(aes(x=Dataset, y=Dim.Reduction.Time.log, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(0,5))  )
grid.arrange(res1, res2, nrow = 2)


dev.off()


