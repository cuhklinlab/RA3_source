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


plot_res <- function(data_type, data_name, bulk_name, peak_selection_k, Dim, multi_donor, gene_score, gene_name, DRres_mat, DRmethod){
  
  set.seed(23)
  tsne_H_hat_12  <- Rtsne( DRres_mat )$Y

    score <- gene_score[which(gene_score$X==gene_name),2:5336]
    score <- as.numeric(score)
    cols = brewer.pal(9, "RdBu")
    pal = colorRampPalette(c(cols[6],cols[1]))
    cells_colors_rgb = pal(100)[as.numeric(cut(score,breaks=100))]
    plot(tsne_H_hat_12,  col=cells_colors_rgb, pch=20, main=sprintf('%s_%s',DRmethod,gene_name))
    legend.col(col = pal(100), lev = score)
}

data_type <- 'PBMC'
data_name <- 'pbmc5k_10xPeak'
multi_donor <- F
bulk_name <- '_Bulkpbmc5k'
peak_selection_k <- 3
Dim <- 10


pdf(file="../results/plot/0503/case/pbmc5k.pdf", width=24, height=4)
par(mfrow=c(1,6))
gene_names <- c('S100A12','MS4A1','PROM1')

for (gene_name in gene_names){
  score <- gene_score[which(gene_score$X==gene_name),2:5336]
  if (nrow(score) == 0){
    print(gene_name)
  }

  for (DRmethod in c("SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC","RA3")){

    if (DRmethod!="RA3"){
      task_name <- sprintf('%s_%s_peak%03d_dim%d',DRmethod,data_name,100*peak_selection_k,Dim)
      rds_file_name <- sprintf('/home/zhangwenyu/dataForComparison_RDS/10xpbmc/%s.rds',task_name)
      if (!file.exists(rds_file_name)) {
        print(sprintf('===[DOES NOT EXIST]=== %s', rds_file_name))
        plot(1:2, 1:2, type="n")
        next
      } 
      DRres <- readRDS(rds_file_name)
      DRres_mat <- t(DRres$FM)
    }else{
      task_name <- sprintf('%s_%s_peak%03d_dim%d%s',DRmethod,data_name,100*peak_selection_k,5,bulk_name)
      mat_file_name <- sprintf('~/data/resultForComparison/%s.mat',task_name)
      if (!file.exists(mat_file_name)) {
        print(sprintf('===[DOES NOT EXIST]=== %s', task_name))
        next
      } 
      DRres <- readMat(mat_file_name)
      if (DRres$K2.left[1]==0){
        DRres_mat <- t(DRres$H.hat)[,1:(DRres$K1[1])]
      }else{
        DRres_mat1 <- t(DRres$H.hat)[,1:(DRres$K1[1])]
        # DRres_mat2 <- t(DRres$H.hat)[, (DRres$sparse.index.left + DRres$K1[1])]
        DRres_mat2 <- t(DRres$H2.trun)
        DRres_mat  <- cbind(DRres_mat1, DRres_mat2)
      }
    }

    print(task_name)
    plot_res(data_type, data_name, bulk_name, peak_selection_k, Dim, multi_donor, gene_score, gene_name, DRres_mat, DRmethod)
  }
}

dev.off()