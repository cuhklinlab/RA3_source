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
plot_bar <- function(clusterMethod,mainDRdim,ourDRdim){
  dataName = 'donor_BM0828'
  data_ind <- which(res_table$Clustering.methods==clusterMethod & 
    res_table$Dataset==dataName & res_table$Peak.selection.k==3 & 
    ( ((res_table$DR.methods %in% compare_method) & (res_table$Dimension.of.latent.space==mainDRdim)) | 
      ((res_table$DR.methods=="RA3") & (res_table$Dimension.of.latent.space==ourDRdim)) | 
      ((res_table$DR.methods=='scABC') & (res_table$Dimension.of.latent.space==999)) )
  )
  data <- res_table[data_ind,]
  for (DRmethod in compare_method){
    if (!(DRmethod %in% data$DR.methods)){
      data <- data[c(1, 1:dim(data)[1]),]
      data[1,]$DR.methods <- DRmethod
      data[1,]$ARI <- 0
      data[1,]$Homogeneity <- 0
      data[1,]$NMI <- 0
      data[1,]$AMI <- 0
    }
  }
  for (dataName in c('MPP_LMPP_CLP','ALL_blood','GMvsHL','forebrain_half','tissue4_forebrain')){
    data_ind <- which(res_table$Clustering.methods==clusterMethod & 
      res_table$Dataset==dataName & res_table$Peak.selection.k==3 & 
      ( ((res_table$DR.methods %in% compare_method) & (res_table$Dimension.of.latent.space==mainDRdim)) | 
        ((res_table$DR.methods=="RA3") & (res_table$Dimension.of.latent.space==ourDRdim)) | 
        ((res_table$DR.methods=='scABC') & (res_table$Dimension.of.latent.space==999)) )
    )
    data_tmp <- res_table[data_ind,]
    for (DRmethod in compare_method){
      if (!(DRmethod %in% data_tmp$DR.methods)){
        data_tmp <- data_tmp[c(1, 1:dim(data_tmp)[1]),]
        data_tmp[1,]$DR.methods <- DRmethod
        data_tmp[1,]$ARI <- 0
        data_tmp[1,]$Homogeneity <- 0
        data_tmp[1,]$NMI <- 0
        data_tmp[1,]$AMI <- 0
      }
    }
    data <- rbind(data, data_tmp)
  }

  data$Dataset <- factor(data$Dataset, levels=c('MPP_LMPP_CLP','donor_BM0828','ALL_blood','GMvsHL','forebrain_half','tissue4_forebrain'))

  cat(dim(data))
  ari <- (ggplot(data) + geom_bar(aes(x=Dataset, y=ARI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$ARI,0),1)) )
  hom <- (ggplot(data) + geom_bar(aes(x=Dataset, y=Homogeneity, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$Homogeneity,0),1)) )
  nmi <- (ggplot(data) + geom_bar(aes(x=Dataset, y=NMI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$NMI,0),1)) )
  ami <- (ggplot(data) + geom_bar(aes(x=Dataset, y=AMI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$AMI,0),1)) )
  grid.arrange(ari, hom, nmi, ami, ncol = 1)
  cat(unique(data$Peak.selection.k))
}


mainDRdim <- 10
ourDRdim <- 5
pdf(file=sprintf("../results/plot/0503/bar_multidonor_fig3-2_dim%d-%d.pdf",mainDRdim,ourDRdim), width=24, height=15)
# mypar(mfrow=c(1, 1))
clusterMethod <- "Louvain clustering"
plot_bar(clusterMethod,mainDRdim,ourDRdim)
dev.off()
