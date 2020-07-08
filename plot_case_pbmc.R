library(rafalib)
library(Rtsne)
library(umap)
library(R.matlab)
library(ggplot2)
library(RColorBrewer)
library(hash)
library(gridExtra)

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
# theme_set(theme_gray())

legend.col <- function(col, lev){
   
  opar <- par
   
  n <- length(col)
   
  bx <- par("usr")
   
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
  bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
   
  xx <- rep(box.cx, each = 2)
   
  par(xpd = TRUE)
  for(i in 1:n){
   
  yy <- c(box.cy[1] + (box.sy * (i - 1)),
  box.cy[1] + (box.sy * (i)),
  box.cy[1] + (box.sy * (i)),
  box.cy[1] + (box.sy * (i - 1)))
  polygon(xx, yy, col = col[i], border = col[i])
   
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
  ylim = c(min(lev), max(lev)),
  yaxt = "n", ylab = "",
  xaxt = "n", xlab = "",
  frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}

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


plot_res <- function(data_type, data_name, bulk_name, peak_selection_k, K2, multi_donor, gene_score, gene_name, data){
  
  if (multi_donor == T) {
    task_name <- sprintf('Our_%s_peak%03d_dim%d_multidonor%s',data_name,100*peak_selection_k,K2,bulk_name);
  }else{
    task_name <- sprintf('Our_%s_peak%03d_dim%d%s',data_name,100*peak_selection_k,K2,bulk_name);
  }
  cat(task_name,'\n')
  # data <- readMat(sprintf('~/data/resultForComparison/%s.mat',task_name))
  
  score <- gene_score[which(gene_score$X==gene_name),2:5336]
  score <- as.numeric(score)

  cell_barcode <- read.csv('~/data/dataForComparison_0415/pbmc5k_cell_barcode.txt', sep='\t', header=FALSE)
  cat(sum(gsub("-", ".", cell_barcode[1:5335,]) == colnames(gene_score)[2:5336]),'\n')

  if (data$K2.left[1]==0){
    DRres_mat <- t(data$H.hat)[,1:(data$K1[1])]
  }else{
    DRres_mat1 <- t(data$H.hat)[,1:(data$K1[1])]
    # DRres_mat2 <- t(data$H.hat)[, (data$sparse.index.left + data$K1[1])]
    DRres_mat2 <- t(data$H2.trun)
    DRres_mat  <- cbind(DRres_mat1, DRres_mat2)
  }
  set.seed(23)
  tsne_sc_pca    <- Rtsne(data$sc.pca.10 )$Y
  set.seed(23)
  tsne_bulk_proj <- Rtsne(data$bulk.proj )$Y
  set.seed(23)
  tsne_H_hat_1   <- Rtsne( t(data$H.hat)[,1:(data$K1)] )$Y
  set.seed(23)
  tsne_H_hat_12  <- Rtsne( DRres_mat )$Y
  # tsne_H_hat_12  <- Rtsne( t(data$H.hat)[,1:(data$K1+data$K2)] )$Y
  set.seed(23)
  umap_sc_pca    <- umap(data$sc.pca.10)$layout
  set.seed(23)
  umap_bulk_proj <- umap(data$bulk.proj)$layout
  set.seed(23)
  umap_H_hat_1  <- umap( t(data$H.hat)[,1:(data$K1)] )$layout
  set.seed(23)
  umap_H_hat_12 <- umap( DRres_mat )$layout
  # umap_H_hat_12 <- umap( t(data$H.hat)[,1:(data$K1+data$K2)] )$layout
  
  
  cols = brewer.pal(9, "RdBu")
  pal = colorRampPalette(c(cols[6],cols[1]))
  cells_colors_rgb = pal(100)[as.numeric(cut(score,breaks=100))]

  pdf(file=sprintf("../results/plot/0503/case/%s_%s.pdf",task_name,gene_name), width=16, height=8)
  par(mfrow=c(2,4))

  plot(tsne_sc_pca,    col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(tsne_bulk_proj, col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(tsne_H_hat_1,   col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(tsne_H_hat_12,  col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(umap_sc_pca,    col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(umap_bulk_proj, col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(umap_H_hat_1,   col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  plot(umap_H_hat_12,  col=cells_colors_rgb, pch=20)
  legend.col(col = pal(100), lev = score)
  
  par(mfrow=c(2,4))
  
  for (i in 1:nrow(data$H.hat)){
    plot(1:ncol(data$H.hat),data$H.hat[i,], col=cells_colors_rgb, pch=20, ylab=i)
    legend.col(col = pal(100), lev = score)
  }
  
  
  dev.off()
}


data_type <- 'PBMC'
data_name <- 'pbmc5k'
multi_donor <- F
bulk_name <- ''
peak_selection_k <- 0
K2 <- 5

if (multi_donor == T) {
  task_name <- sprintf('Our_%s_peak%03d_dim%d_multidonor%s',data_name,100*peak_selection_k,K2,bulk_name);
}else{
  task_name <- sprintf('Our_%s_peak%03d_dim%d%s',data_name,100*peak_selection_k,K2,bulk_name);
}
data <- readMat(sprintf('~/data/resultForComparison/%s.mat',task_name))

for (gene_name in c('S100A12', 'MS4A1', 'GAPDH')){
  cat(gene_name, '\n')
  plot_res(data_type, data_name, bulk_name, peak_selection_k, K2, multi_donor, gene_score, gene_name, data)
}

