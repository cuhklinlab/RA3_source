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

plot_scatter <- function(dims.reduce, cell.labels, data_type='NULL', cell.colors=NULL, donor.labels=NULL, multi_donor=F, perplexity=30,DRmethod,data_name){
  
  # cell.labels <- cluster_assign$label
  # dims.reduce <- t(DRres$FM)
  palette_name <- "Dark2"
  
  if (data_type=='blood') {
    # cells_colors_rgb <- c()
    # for (i in 1:nrow(cell.colors)){
    #   cells_colors_rgb <- c( cells_colors_rgb, rgb(cell.colors[i,1], cell.colors[i,2], cell.colors[i,3]) ) 
    # }
    ## plot the celltype color
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
    cells_colors_rgb <- values(color_map, keys=cell.labels )
    
    # plot(1:2, 1:2, type="n")
    # legend("bottomright", legend=cellTypeAll_U, col=legend_cellType_color_rgb, pch=20, bty = "n")
    type_ind <- c()
    for (i in 1:length(cellTypeAll_U)){
      if (cellTypeAll_U[i] %in% unique(cell.labels) & !(cellTypeAll_U[i] %in% c('GMP1low','GMP2mid','GMP3high'))){
        type_ind <- c( type_ind, i)
      }
    }
    cellTypeAll_U <- cellTypeAll_U[type_ind]
    legend_cellType_color_rgb <- legend_cellType_color_rgb[type_ind]
    
  }else{
    
    cellTypeAll_U <- unique(cell.labels)
    legend_cellType_color_rgb <- brewer.pal(length(unique(cell.labels)), name=palette_name)
    color_map <- hash(cellTypeAll_U, legend_cellType_color_rgb)
    cells_colors_rgb <- values(color_map, keys=cell.labels )

    # plot(1:2, 1:2, type="n")
    # legend("bottomright", legend=cellTypeAll_U, col=legend_cellType_color_rgb, pch=20, bty = "n")
    levels(cellTypeAll_U) <- c(levels(cellTypeAll_U), c('AC','EX CPN','EX CThPN','EX SCPN','IN','MG','OC'))
    cellTypeAll_U[(cellTypeAll_U) == "Astrocytes"] <- 'AC'
    cellTypeAll_U[(cellTypeAll_U) == "Ex. neurons CPN"] <- 'EX CPN'
    cellTypeAll_U[(cellTypeAll_U) == "Ex. neurons CThPN"] <- 'EX CThPN'
    cellTypeAll_U[(cellTypeAll_U) == "Ex. neurons SCPN"] <- 'EX SCPN'
    cellTypeAll_U[(cellTypeAll_U) == "Inhibitory neurons"] <- 'IN'
    cellTypeAll_U[(cellTypeAll_U) == "Microglia"] <- 'MG'
    cellTypeAll_U[(cellTypeAll_U) == "Oligodendrocytes"] <- 'OC'
    sort_index <- order(cellTypeAll_U)
    cellTypeAll_U <- cellTypeAll_U[sort_index]
    legend_cellType_color_rgb <- legend_cellType_color_rgb[sort_index]
  }
  
  donor_shape <- rep.int(20,length(cell.labels))
  if (multi_donor == T) {
    donor.labels <- unlist(donor.labels)
    donor_index <- as.numeric( factor(donor.labels) )
    donor_shape <- donor_index
    donor_shape[which(donor_index==1)] <- 4
    donor_shape[which(donor_index==2)] <- 16
    donor_shape[which(donor_index==3)] <- 15
    donor_shape[which(donor_index==4)] <- 17
    donor_shape[which(donor_index==5)] <- 3
    donor_shape[which(donor_index==6)] <- 8
    donor_shape[which(donor_index==7)] <- 11
    
    # pdf(file="../results/plot/0512/MPPLMPPCLP_donor_shape.pdf")
    # plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    # legend("center", legend=unique(donor.labels), col="black", pch=unique(donor_shape), bty="n", cex=1.47)
    # dev.off()
  }
  
  set.seed(23)
  tsne_DRres    <- Rtsne(dims.reduce, perplexity=perplexity)$Y
  # set.seed(23)
  # umap_DRres    <- umap(dims.reduce)$layout
  
  
  if (DRmethod=="SCALE"){
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("center", legend=cellTypeAll_U, col=legend_cellType_color_rgb, pch=20, bty="n", cex=2, title=final_data_name(data_name))#, inset=c(-0.5,0)   par(xpd=TRUE)for outside
  }

  if (data_name %in% c('ALL_blood','tissue4_forebrain','tissue2_forebrain','tissue3_forebrain')){
    plot(tsne_DRres,    col=cells_colors_rgb, pch=donor_shape,xlab="",ylab="", cex = 0.6)#, main=file_name
  }
  else{
    plot(tsne_DRres,    col=cells_colors_rgb, pch=donor_shape,xlab="",ylab="")#, main=file_name
  }
  
}


final_data_name <- function(data_name){
  if (data_name=="MPP_LMPP_CLP"){return("CLP/LMPP/MPP")}
  if (data_name=="donor_BM0828"){return('Donor BM0828')}
  if (data_name=="ALL_blood"){return("Human bone marrow")}
  if (data_name=="GMvsHL"){return("GM/HL")}
  if (data_name=="GMvsHek"){return("GM/HEK")}
  if (data_name=="InSilico"){return("InSilico mixture")}
  if (data_name=="forebrain_half"){return("Mouse forebrain")}
  if (data_name=="forebrain_half_do025"){return("Mouse forebrain\n(dropout rate=25%)")}
  if (data_name=="tissue4_forebrain"){return("MCA mouse\nprefrontal cortex")}
  if (data_name=="tissue1_forebrain"){return("MCA mouse\ncerebellum")}
  if (data_name=="tissue2_forebrain"){return("MCA mouse\nwhole brain 1")}
  if (data_name=="tissue3_forebrain"){return("MCA mouse\nwhole brain 2")}
}





intv = '_interval'
interval = seq(0.05,0.25,0.05)
mainDRdim <- '10'
ourDRdim <- '5'
plot_size = 3.8
multi_donor = T;
# multi_donor = F;
if (multi_donor == T){
  pdf(file=sprintf("../results/plot/1021/Q6_scatter_theta_multidonor%s.pdf",intv), width=length(interval)*plot_size, height=3*plot_size)
}else{
  pdf(file=sprintf("../results/plot/1021/Q6_scatter_theta%s.pdf",intv), width=length(interval)*plot_size, height=3*plot_size)
}
tsne_perplexity <- 30
# datasets <- unique(res_table$Dataset)
# dataset_ind <- 0
peakR <- 3

par(mfrow=c(3,length(interval)), pty="s", mar=c(1, 3, 1, 1) + 0.1, cex.axis=1.5) 
for (data_name in c('ALL_blood','forebrain_half','tissue4_forebrain')){
for (theta in interval){
      DRmethod = "RA3"
      if (DRmethod=="RA3"){
        DRdim <- ourDRdim
        if ((multi_donor == T) & (data_name %in% c('ALL_blood','MPP_CMP_CLP','MPP_LMPP_CLP','all_forebrain','all_forebrain_RmEx'))){
          file_name <- sprintf('%s_%s_peak%03d_dim%s_multidonor_theta%03d',DRmethod,data_name,peakR*100,DRdim, round(100*theta))
        }else{
          file_name <- sprintf('%s_%s_peak%03d_dim%s_theta%03d',DRmethod,data_name,peakR*100,DRdim, round(100*theta))
        }
        label_file_name <- sprintf('/home/zhangwenyu/dataForComparison_clustering/202010/clusters/%s_clusters.tsv',file_name)
        mat_file_name <- sprintf('~/data/RA3/R1_results/%s.mat',file_name)
        if (!file.exists(mat_file_name)) {
          print(sprintf('===[DOES NOT EXIST]=== %s', file_name))
          next
        } 
        DRres <- readMat(mat_file_name)
        if (DRres$K2.left[1]==0){
          DRres_mat <- t(DRres$H.hat)[,1:(DRres$K1[1])]
        }else{
          DRres_mat1 <- t(DRres$H.hat)[,1:(DRres$K1[1])]
          DRres_mat2 <- t(DRres$H2.trun)
          DRres_mat  <- cbind(DRres_mat1, DRres_mat2)
        }
      }else{
        DRdim <- mainDRdim
        file_name <- sprintf('%s_%s_peak%03d_dim%s_theta%03d',DRmethod,data_name,peakR*100,DRdim, round(100*theta))
        label_file_name <- sprintf('/home/zhangwenyu/dataForComparison_clustering/202010/clusters/%s_clusters.tsv',file_name)
        rds_file_name <- sprintf('/home/zhangwenyu/dataForComparison_RDS/202010/%s.rds',file_name)
        if (!file.exists(rds_file_name)) {
          print(sprintf('===[DOES NOT EXIST]=== %s', rds_file_name))

          file_name <- sprintf('%s_%s_peak%03d_dim%s','cisTopicCGS',data_name,peakR*100,DRdim)
          label_file_name <- sprintf('/home/zhangwenyu/dataForComparison_clustering/clusters/%s_clusters.tsv',file_name)
          cluster_assign <- read.csv(label_file_name, sep='\t')
          cell.labels = cluster_assign$label

          if (dataset_ind < 4) {
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
            
            type_ind <- c()
            for (i in 1:length(cellTypeAll_U)){
              if (cellTypeAll_U[i] %in% unique(cell.labels) & !(cellTypeAll_U[i] %in% c('GMP1low','GMP2mid','GMP3high'))){
                type_ind <- c( type_ind, i)
              }
            }
            cellTypeAll_U <- cellTypeAll_U[type_ind]
            legend_cellType_color_rgb <- legend_cellType_color_rgb[type_ind]
          
          }else{
            palette_name <- "Dark2"
            cellTypeAll_U <- unique(cell.labels)
            legend_cellType_color_rgb <- brewer.pal(length(unique(cell.labels)), name=palette_name)

            levels(cellTypeAll_U) <- c(levels(cellTypeAll_U), c('AC','EX CPN','EX CThPN','EX SCPN','IN','MG','OC'))
            cellTypeAll_U[(cellTypeAll_U) == "Astrocytes"] <- 'AC'
            cellTypeAll_U[(cellTypeAll_U) == "Ex. neurons CPN"] <- 'EX CPN'
            cellTypeAll_U[(cellTypeAll_U) == "Ex. neurons CThPN"] <- 'EX CThPN'
            cellTypeAll_U[(cellTypeAll_U) == "Ex. neurons SCPN"] <- 'EX SCPN'
            cellTypeAll_U[(cellTypeAll_U) == "Inhibitory neurons"] <- 'IN'
            cellTypeAll_U[(cellTypeAll_U) == "Microglia"] <- 'MG'
            cellTypeAll_U[(cellTypeAll_U) == "Oligodendrocytes"] <- 'OC'
            sort_index <- order(cellTypeAll_U)
            cellTypeAll_U <- cellTypeAll_U[sort_index]
            legend_cellType_color_rgb <- legend_cellType_color_rgb[sort_index]
          }

          plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
          legend("center", legend=cellTypeAll_U, col=legend_cellType_color_rgb, pch=20, bty="n", cex=2, title=final_data_name(data_name))
          plot(-5:5, -5:5, type="n", ylab='',xlab='')
          next
        } 
        DRres <- readRDS(rds_file_name)
        DRres_mat <- t(DRres$FM)
      }

      if (sum(is.na(DRres$FM)) > 0) {
        print(sprintf('[NaN in DRres] %s', file_name))
        next
      } 
      print(file_name)

      if (!file.exists(label_file_name)) {
        print(sprintf('===[DOES NOT EXIST]=== %s', label_file_name))
        next
      } 
      cluster_assign <- read.csv(label_file_name, sep='\t')
      

      if (data_name %in% c('MPP_LMPP_CLP','donor_BM0828','ALL_blood')){
        plot_scatter(DRres_mat, cluster_assign$label, data_type='blood', perplexity=tsne_perplexity, DRmethod=DRmethod,data_name=data_name) # using color of source paper
      }else{
        plot_scatter(DRres_mat, cluster_assign$label, perplexity=tsne_perplexity, DRmethod=DRmethod,data_name=data_name) # using automatic color
      }
    # }
  }
}
dev.off()

