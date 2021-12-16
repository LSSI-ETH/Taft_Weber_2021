
#' this function generates prediction value heatmaps
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param thres_range numeric vector of threshold determine the upper and lower limit
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")
#' @return heatmaps of prediction values in working directory, legend is exported in separate file

pred_value_heatmap<-function(data, distance=c(1,2), threshold=0.5,thres_range=c(0.25,0.75), model=c("RF","RNN"), variant=c("alpha","beta","gamma","kappa","wuhan","new")){
  heatmap<-list()
  data_dis<-list()
  data_model<-data_all<- list()
  ace_list<-ab_list<-all_list<-list()
  data1<-list()
  for (v in 1:length(variant) ){
    range<-1:ncol(data[[v]])
    colnames(data[[v]])<-str_remove_all(colnames(data[[v]]),pattern = "_isotonic")
    colnames(data[[v]])<-str_remove_all(colnames(data[[v]]),pattern = "_sigmoid")
    colnames(data[[v]])<-str_remove_all(colnames(data[[v]]),pattern = "_uncalibrated")
    for( t in 1:length(threshold)){
      ##filter out ace2 non-binding
      ace2<-.binary(threshold[t], data1 = as.data.frame(data[[v]]$RF_ACE2), data2 = as.data.frame(data[[v]]$RNN_ACE2), disagree = 0)
      data1[[v]]<-data[[v]][ace2==1,]

      for(j in 1:length(distance)){

        #row selection + partial col selection----
        data_dis[[j]] <-subset(x = data1[[v]],subset = data1[[v]][,str_detect(colnames(data1[[v]]),pattern="ist")]==distance[j],select = range[str_detect(colnames(data1[[v]]),pattern="RF")|str_detect(colnames(data1[[v]]),pattern="RNN")])

        #split model
        for(k in 1:length(model)){
          modelselect<-stringr::str_detect(string = colnames(data_dis[[j]]),pattern = model[k])
          data_model[[k]]<-subset(x=data_dis[[j]],select = modelselect)
        }

        data_all[[j]]<-data_model

        #把每个data_dist 弄成rf rnn挨在一起的。
        #combine by col(antibody)

        for(col in 2:ncol(data_model[[1]])){
          for(k in 1:length(model)){
            uno<-data.frame(data_model[[k]][,col])
            colnames(uno)<-colnames(data_model[[k]])[col]
            if(k==1) duo<-uno else
              duo<-cbind(duo,uno)
          }
          if(col==2) tri<-duo else
            tri<-cbind(tri,duo)
        }
        data_dis[[j]]<-tri
      }
      #TODO: do binary first and find ranking order
      for(i in 1:length(distance)){
        ##note ab escaping: reverse binding value, value means escape. T means yes
        bi<-.binary(th=threshold[t],data1=data_all[[i]][[1]],data2=data_all[[i]][[2]], disagree = -1)
        #how many 0 == how many escapes
        n<-apply(bi[,-1],FUN = .sum0, MARGIN = 1)
        combi<-c()
        for(row in 1:nrow(bi)) combi[row]<-paste(colnames(bi)[which(bi[row,]==1)],collapse = "_")

        rank<-order(n,combi)

        bi3<-.binary(th=c(0.25,0.75),data1=data_all[[i]][[1]],data2=data_all[[i]][[2]])
        n<-nrow(data_all[[i]][[1]])
        ##TODO 根据距离换图长度
        if(distance[i]==1) {
          Y<-F
          h<-500}

        if(distance[i]==2) {
          h<-1000
          Y<-T}


        # pdf(file=paste("bi3_all" , variant[v] ,"dist",distance[i], "threshold", threshold[t],"n",n,".pdf",sep = "_")
        #     ,
        #     width = 5,
        #     height = h
        # )
        #
        # pheatmap( mat = bi3[rank,2:ncol(data_model[[1]])],
        #           border=F,
        #           cellwidth = 10,
        #           cluster_rows = F, cluster_cols = F,
        #           breaks = c(-2.1,-1.2,-0.3,0.6,1.5),
        #           color = c("#999999","#FFFFBF","#4575B4","#D73027"),
        #           legend_breaks = c(-1.2,-0.3,0.6,1.5)-0.45,
        #           legend_labels=c("Disagree","NA","False","True"),
        #           show_rownames = F,
        #           show_colnames = Y,legend = F
        # )
        #
        # dev.off()
        #
        # pdf(file=paste("bi_all" , variant[v] ,"dist",distance[i], "threshold", threshold[t],"n",n,".pdf",sep = "_")
        #     ,
        #     width = 5,
        #     height = h
        # )
        #
        # pheatmap( mat = bi[rank,2:ncol(data_model[[1]])],
        #           border=F,
        #           cellwidth = 10,
        #           cluster_rows = F, cluster_cols = F,
        #           breaks = c(-1,-0.1,0.8,1.7),
        #           color = c("#FFFFBF","#4575B4","#D73027"),
        #           legend_breaks = c(-0.1,0.8,1.7),
        #           legend_labels=c("NA","False","True"),
        #           show_rownames = F,
        #           show_colnames = Y,
        #           legend = F
        # )
        #
        # dev.off()


        #do heat map with data_dis continuous

        # pdf(file =  paste("heatmap_reorder" , variant[v] ,"dist",distance[i],"threshold",threshold[t],"n",n,".pdf",sep = "_"),
        #     width = 5,
        #     height =h,
        # )
        #
        # pheatmap( mat = data_dis[[i]][rank,],
        #           gaps_col = seq(2,2*ncol(data_model[[1]])-1,2),
        #           border=F,
        #           cellwidth = 10,
        #           cluster_rows = F, cluster_cols = F,breaks=seq(0,1,0.01),
        #           show_rownames = F,
        #           show_colnames = Y,
        #           legend = F
        #
        # )
        #
        # dev.off()

        png(filename =  paste("heatmap_reorder" , variant[v] ,"dist",distance[i],"threshold",threshold[t],"n",n,".png",sep = "_"),
            units = "px",
            width = 500,
            height =h,
        )

        pheatmap::pheatmap( mat = data_dis[[i]][rank,],
                  gaps_col = seq(2,2*ncol(data_model[[1]])-1,2),
                  border=F,
                  cellwidth = 10,
                  cluster_rows = F, cluster_cols = F,breaks=seq(0,1,0.01),
                  show_rownames = F,
                  show_colnames = Y,
                  legend = F

        )

        dev.off()

      }
    }
  }
}
