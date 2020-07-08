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

res_table <- read.csv("./input_of_R/res_multidonor_0503_truncate.csv")
# head(res_table)
res_table$Peak.selection.k <- as.character(res_table$Peak.selection.k)
res_table$DR.methods <- factor(res_table$DR.methods, levels=c("scABC","SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC","RA3"))

# unique(res_table$Dataset)
    
compare_method <- c("scABC","SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC")
plot_line <- function(dataName_ind,clusterMethod){
  # dataName <- 'ALL_blood'
  # coul <- brewer.pal(8, "Dark2")
  # coul <- coul[c(1,2,3,4,5,6,8)]
  # scale_color_manual(values = coul)
  

  data <- which(res_table$Clustering.methods==clusterMethod & res_table$Dataset.number%in%dataName_ind & res_table$Peak.selection.k==3 & ((res_table$DR.methods %in% compare_method) & (res_table$Dimension.of.latent.space == 10) | (res_table$DR.methods == "RA3") & (res_table$Dimension.of.latent.space == 5) | (res_table$DR.methods == 'scABC') & (res_table$Dimension.of.latent.space == 999))
)
  data <- res_table[data,]
  ari <- (ggplot(data, aes(x=Dataset.number, y=ARI, colour=DR.methods)) + geom_line(size=1.5) + scale_color_brewer(palette="Paired") + ylim(min(data$ARI,0),1) )
  hom <- (ggplot(data, aes(x=Dataset.number, y=Homogeneity, colour=DR.methods)) + geom_line(size=1.5) + scale_color_brewer(palette="Paired") + ylim(min(data$Homogeneity,0),1) )
  nmi <- (ggplot(data, aes(x=Dataset.number, y=NMI, colour=DR.methods)) + geom_line(size=1.5) + scale_color_brewer(palette="Paired") + ylim(min(data$NMI,0),1) )
  ami <- (ggplot(data, aes(x=Dataset.number, y=AMI, colour=DR.methods)) + geom_line(size=1.5) + scale_color_brewer(palette="Paired") + ylim(min(data$AMI,0),1) )
  # mean_se for calculating mean and standard error of different. Refer to  https://ggplot2.tidyverse.org/reference/ and https://stackoverflow.com/questions/44872951/how-do-i-add-se-error-bars-to-my-barplot-in-ggplot2


  # ari <- (ggplot(res_table[which(res_table$Dataset==dataName),], mapping=aes(x=Peak.selection.k, y=ARI, fill=DR.methods)) + stat_summary(geom="bar", fun=mean, position="dodge") +  stat_summary(geom="errorbar", fun.data=mean_se, position="dodge") + scale_fill_brewer(palette="Paired") + ylim(0,1) + labs(title=dataName) + theme(legend.position="none") )
  # hom <- (ggplot(res_table[which(res_table$Dataset==dataName),], mapping=aes(x=Peak.selection.k, y=Homogeneity, fill=DR.methods)) + stat_summary(geom="bar", fun=mean, position="dodge") +  stat_summary(geom="errorbar", fun.data=mean_se, position="dodge") + scale_fill_brewer(palette="Paired") + ylim(0,1) + labs(title=dataName) + theme(legend.position="none") )
  # nmi <- (ggplot(res_table[which(res_table$Dataset==dataName),], mapping=aes(x=Peak.selection.k, y=NMI, fill=DR.methods)) + stat_summary(geom="bar", fun=mean, position="dodge") +  stat_summary(geom="errorbar", fun.data=mean_se, position="dodge") + scale_fill_brewer(palette="Paired") + ylim(0,1) + labs(title=dataName) + theme(legend.position="none") )
  # ami <- (ggplot(res_table[which(res_table$Dataset==dataName),], mapping=aes(x=Peak.selection.k, y=AMI, fill=DR.methods)) + stat_summary(geom="bar", fun=mean, position="dodge") +  stat_summary(geom="errorbar", fun.data=mean_se, position="dodge") + scale_fill_brewer(palette="Paired") + ylim(0,1) + labs(title=dataName) )
  grid.arrange(ari, hom, nmi, ami, nrow = 1)
}

pdf(file="../results/plot/0503/line_multidonor_fig4.pdf", width=42, height=6)
# mypar(mfrow=c(1, 1))
dataName_ind = c(7:17)
clusterMethod <- "Louvain clustering"
plot_line(dataName_ind,clusterMethod)
dev.off()