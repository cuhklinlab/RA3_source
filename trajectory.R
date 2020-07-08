suppressPackageStartupMessages({
library(Rtsne)
library(umap)
library(R.matlab)
library(ggplot2)
library(RColorBrewer)
library(hash)
library(slingshot)
})
# display.brewer.all(colorblindFriendly=TRUE)
theme_set(theme_gray() +
            theme(
              axis.line = element_line(size=0.5),
              panel.background = element_rect(fill=NA,size=rel(20)),
              panel.grid.minor = element_line(colour = NA),
              axis.text = element_text(size=16),
              axis.title = element_text(size=18)
            )
)

plot_res <- function(data_type, data_name, bulk_name, peak_selection_k, K2, multi_donor){
  
  if (multi_donor == T) {
    task_name <- sprintf('Our_%s_peak%03d_dim%d_multidonor%s',data_name,100*peak_selection_k,K2,bulk_name);
  }else{
    task_name <- sprintf('Our_%s_peak%03d_dim%d%s',data_name,100*peak_selection_k,K2,bulk_name);
  }
  data <- readMat(sprintf('~/data/resultForComparison/%s.mat',task_name))
  
  label_file_name <- sprintf('/home/zhangwenyu/dataForComparison_clustering/clusters_our/%s_clusters.tsv',task_name)
  cluster_assign <- read.csv(label_file_name, sep='\t')
  data$cell.labels <- cluster_assign$label
  data$cell.labels <- unlist(data$cell.labels)

  palette_name <- "Dark2"
  
  if (data_type=='blood') {
    legend_cellType_color <- read.csv("./input_of_R/blood/cellTypeAll1_U_color.csv", header=FALSE)
    legend_cellType_color <- as.matrix(legend_cellType_color)
    legend_cellType_color_rgb <- c()
    for (i in 1:nrow(legend_cellType_color)){
      legend_cellType_color_rgb <- c( legend_cellType_color_rgb, rgb(legend_cellType_color[i,1],
                                                                     legend_cellType_color[i,2],
                                                                     legend_cellType_color[i,3]) ) 
    }
    cellTypeAll_U <- read.csv("./input_of_R/blood/cellTypeAll1_U.csv", header=FALSE)
    cellTypeAll_U <- c(as.matrix(cellTypeAll_U))
    color_map <- hash(cellTypeAll_U, legend_cellType_color_rgb)
    cells_colors_rgb <- values(color_map, keys=data$cell.labels )
  }else{
    
    cellTypeAll_U <- unique(data$cell.labels)
    # legend_cellType_color_rgb <- colorRampPalette(brewer.pal(8, name=palette_name))(length(unique(data$cell.labels)))
    legend_cellType_color_rgb <- brewer.pal(length(unique(data$cell.labels)), name=palette_name)
    color_map <- hash(cellTypeAll_U, legend_cellType_color_rgb)
    cells_colors_rgb <- values(color_map, keys=data$cell.labels )
  }

  donor_shape <- rep.int(20,length(data$cell.labels))
  if (multi_donor == T) {
    donor_file_name <- sprintf('~/data/resultForComparison/%s_donor_labels.mat',data_name)
    data$donor.labels <- readMat(donor_file_name)$donor.label
    data$donor.labels <- unlist(data$donor.labels)
    donor_index <- as.numeric( factor(data$donor.labels) )
    donor_shape <- donor_index
    donor_shape[which(donor_index==1)] <- 4
    donor_shape[which(donor_index==2)] <- 16
    donor_shape[which(donor_index==3)] <- 15
    donor_shape[which(donor_index==4)] <- 17
    donor_shape[which(donor_index==5)] <- 3
    donor_shape[which(donor_index==6)] <- 8
    donor_shape[which(donor_index==7)] <- 11
    
  }
  

  if (data$K2.left[1]==0){
    DRres_mat <- t(data$H.hat)[,1:(data$K1[1])]
  }else{
    DRres_mat1 <- t(data$H.hat)[,1:(data$K1[1])]
    # DRres_mat2 <- t(data$H.hat)[, (data$sparse.index.left + data$K1[1])]
    DRres_mat2 <- t(data$H2.trun)
    DRres_mat  <- cbind(DRres_mat1, DRres_mat2)
  }
      
  set.seed(23)
  tsne_H_hat_12  <- Rtsne( DRres_mat )$Y


  
  pdf(file=sprintf("../results/plot/0512/case/traj_%s_pred.pdf",task_name), width=16, height=4)
  par(mfrow=c(1,4),  pty="s", mar=c(1, 3, 1, 1) + 0.1, cex.axis=1.5)

  if (data_name %in% c('ALL_blood','tissue4_forebrain','tissue2_forebrain','tissue3_forebrain')){
    cex = 0.6
  }
  else{
    cex = 1
  }
  lin1 <- getLineages(DRres_mat, cluster_assign$louvain)
  # lin1 <- getLineages(DRres_mat, cluster_assign$louvain, start.clus= '3')
  crv1 <- getCurves(lin1)


  sds.new <- embedCurves(crv1, tsne_H_hat_12, shrink.method='density')
  plot(tsne_H_hat_12, col = cells_colors_rgb, pch = donor_shape, cex=cex,xlab="",ylab="")
  lines(sds.new, lwd = 3)


  
  dev.off()
}

################################

data_type <- 'blood'
data_name <- 'donor_BM0828'
multi_donor <- F
bulk_name <- ''
peak_selection_k <- 3
K2 <- 5
plot_res(data_type, data_name, bulk_name, peak_selection_k, K2, multi_donor)
