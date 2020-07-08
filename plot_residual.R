library(rafalib)
library(Rtsne)
library(umap)
library(R.matlab)
library(ggplot2)
library(RColorBrewer)
library(hash)

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


plot_jitter <- function(data, label, xaxis, main, pch=3, cex=1, xlabs="Components", ylabs="H_hat value", jitter_prop=0.8, col1=brewer.pal(8, name="Set2")[8], col2=brewer.pal(8, name="Set2")[2], norm=F){
  if (norm){
    data <- scale(data) 
  } 
  x_jitter <- matrix(rep(1:ncol(data), nrow(data)) + jitter_prop*(runif(length(data))-0.5), nrow=nrow(data), byrow=T) 
  ## 1 is common, 2 is rare
  cols <- rep(NA, length(label))
  cols[which(label==1)] <- col1
  cols[which(label==2)] <- col2
  ##
  plot(as.numeric(x_jitter), as.numeric(data), type="n", xaxt="n", main=main, xlab=xlabs, ylab=ylabs)
  for (i in 1:ncol(data)){
    points(x_jitter[,i], data[,i], pch=pch, col=cols, cex=cex)
  }
  axis(1, at=1:ncol(data),labels=xaxis)
}


plot_res <- function(data_type, data_name, bulk_name, peak_selection_k, K2, multi_donor){
  
  if (multi_donor == T) {
    task_name <- sprintf('Our_%s_peak%03d_dim%d_multidonor%s_residual',data_name,100*peak_selection_k,K2,bulk_name);
  }else{
    task_name <- sprintf('Our_%s_peak%03d_dim%d%s_residual',data_name,100*peak_selection_k,K2,bulk_name);
  }
  data <- readMat(sprintf('~/data/RA3_test/%s.mat',task_name))
  print(task_name)
  
  label_file_name <- sprintf('/home/zhangwenyu/dataForComparison_clustering/clusters_our/%s_clusters.tsv',sprintf('Our_%s_peak%03d_dim%d%s',data_name,100*peak_selection_k,K2,''))#clusters_32to35, clusters_our
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
  

  DRres_mat <- data$score.res

  
  pdf(file=sprintf("../results/plot/0512/case/%s.pdf",task_name), width=16, height=4)
  par(mfrow=c(1,4), cex.axis=1.5)

  if (data_name %in% c('ALL_blood','tissue4_forebrain','tissue2_forebrain','tissue3_forebrain')){
    cex = 0.6
  }
  else{
    cex = 1
  }
  for (seed in 110:110){
    set.seed(seed)
    tsne_residual  <- Rtsne( DRres_mat )$Y
    plot(tsne_residual,    col=cells_colors_rgb, pch=donor_shape, cex=cex,xlab="",ylab="")

    set.seed(seed)
    tsne_residual_bulk  <- Rtsne( cbind(DRres_mat, data$bulk.proj) )$Y
    plot(tsne_residual_bulk,    col=cells_colors_rgb, pch=donor_shape, cex=cex,xlab="",ylab="")

  }
  

  plot(1:2, 1:2, type="n")
  legend("bottomright", legend=cellTypeAll_U, col=legend_cellType_color_rgb, pch=donor_shape, bty = "n", cex=cex)
  
  dev.off()
}

###############################
data_type <- 'blood'
data_name <- 'donor_BM0828'
multi_donor <- F
bulk_name <- '_Bulkblood_RSHscMppLmppCmp'
peak_selection_k <- 3
K2 <- 5
plot_res(data_type, data_name, bulk_name, peak_selection_k, K2, multi_donor)
