




#' this function extracts sequences that models don't agree on binding decision.
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")
#' @return dataframe with disagreement sequence and its identity data, variant, distance, threshold, antibody
disagree_sequence_extract<-function(data,variant=c("alpha","beta","gamma","kappa","wuhan","new"),threshold=0.5, distance=c(1,2),model=c("RF","RNN")){
  #----------
  data_dis1<-list()
  extract_list<-list()
  heatmap<-list()
  data_dis<-list()
  data_model<-data_all<- list()
  for (v in 1:length(variant) ){
    range<-1:ncol(data[[v]])
    for( t in 1:length(threshold)){

      for(j in 1:length(distance)){

        #row selection + partial col selection
        data_dis[[j]] <-subset(x=data[[v]],subset = data[[v]][,str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = range[str_detect(colnames(data[[v]]),pattern="RF")|str_detect(colnames(data[[v]]),pattern="RNN")])
        # data_dis1[[j]]<-subset(x=data[[v]],subset = data[[v]][,str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = c(1,2,13))
        data_dis1[[j]]<-subset(x=data[[v]],subset = data[[v]][,str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = c(1,2))
        #split model
        for(k in 1:length(model)){
          modelselect<-stringr::str_detect(string = colnames(data_dis[[j]]),pattern = model[k])
          data_model[[k]]<-subset(x=data_dis[[j]],select = modelselect)
        }

        data_all[[j]]<-data_model

      }
      #TODO: do .binary first and find ranking order
      for(i in 1:length(distance)){
        ##note ab escaping: reverse binding value, value means escape. T means yes
        bi<-.binary(th=threshold[t],data1=data_all[[i]][[1]],data2=data_all[[i]][[2]],disagree = -1)
        bi<-bi[bi[,1]==1,]

        for(e in 2:ncol(bi)){
          extract<-data_dis1[[i]][,1][bi[,1]==1&bi[,e]==-1]#ace>th&NA
          el<-length(extract_list)
          extract_list[[el+1]]<-extract
          if(length(extract)==0){
            extract_list[[el+1]]<-"None"
          }
          names(extract_list)[el+1]<-paste(variant[v],distance[i],threshold[t],colnames(bi)[e],sep = ",")
        }
      }

    }
  }


  extract<-c()
  ##reshape
  for(i in 1:length(extract_list)){
    if(length(extract_list[[i]])==0) next
    if(length(extract)==0) extract<-data.frame(sequence=extract_list[[i]],group=rep(names(extract_list)[i],length(extract_list[[i]])))
    if(length(extract)>0) extract<-rbind(extract,data.frame(sequence=extract_list[[i]],group=rep(names(extract_list)[i],length(extract_list[[i]]))))
  }

  meta<-as.data.frame(str_split(extract$group,pattern = ",",n = 4,simplify = T))
  colnames(meta)<-c("variant","dist","th","ab")
  ab<-str_split (meta$ab,pattern = "_",n = 3,simplify = T)
  ab<-paste(ab[,1],ab[,2],sep="_")
  meta$ab<-ab

  extract<-cbind(extract,meta)
return(extract)

}
