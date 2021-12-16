
#' this function returns all combination of antibody escape, can be used before the stacked bar plot function for color set up.
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param thres_range numeric vector of threshold determine the upper and lower limit
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")

escape_combi<-function(data, distance=c(1,2), threshold=0.5, thres_range=c(0.25,0.75), model=c("RF","RNN"), variant=c("alpha","beta","gamma","kappa","wuhan","new"),n_antibody=4){
  stack_list<-list()
  total_list<-list()
  all_combi<-c()
  nAbEscape<-0:n_antibody
  for (v in 1:length(variant) ){
    range<-1:ncol(data[[v]])
    data_dis<-list()
    data_model<-list()
    data_all<-list()
    for(j in 1:length(distance)){

      #row selection + partial col selection
      data_dis[[j]] <-subset(x=data[[v]],subset = data[[v]][,str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = range[str_detect(colnames(data[[v]]),pattern="RF")|str_detect(colnames(data[[v]]),pattern="RNN")])

      #col selection model


      for(k in 1:length(model)){

        modelselect<-stringr::str_detect(string = colnames(data_dis[[j]]),pattern = model[k])
        data_model[[k]]<-subset(x=data_dis[[j]],select = modelselect)

      }

      data_all[[j]]<-data_model
    }

    for( i in 1:length(distance)){

      for (t in 1:length(threshold)){
        #note ab escaping: reverse binding value, value means escape. T means yes

        bi0<-.binary(th=threshold[t],data1=data_all[[i]][[1]],data2=data_all[[i]][[2]],disagree = 0)
        bi<-.binary(th=c(0.25,0.75),data1=data_all[[i]][[1]],data2=data_all[[i]][[2]])
        bi<-bi[bi0[,1]==1,] ##filter out ace non-binding
        bi<-as.data.frame(bi)
        n<-apply(bi[,-1],FUN = .sum0,MARGIN = 1)
        combi<-c()

        for(row in 1:nrow(bi)) combi[row]<-paste(colnames(bi)[which(bi[row,]==0)],collapse = "_")

        all_combi<-c(all_combi,combi)

        stack<-data.frame(bi,n,combi)
        options(dplyr.summarise.inform = FALSE)
        stack<-dplyr::summarise(dplyr::group_by(stack,n,combi),count = dplyr::n())

        total <- stack %>%
          dplyr::group_by(n) %>%
          dplyr::summarise(total = sum(count))

        stack$combi[stack$combi==""]<-"None"


        if(any(nAbEscape %in% unique(stack$n)==F)) {
          extra<-nrow(stack)+1
          x<-nAbEscape[!(nAbEscape %in% unique(stack$n))]
          stack[extra,1]<-stack[extra,3]<-0
          stack[extra,2]<-"None"
          x<-data.frame(x,rep(0,length(x)))
          colnames(x)<-c("n","total")
          total<- rbind(total,x)
        }


        l<-length(stack_list)
        stack_list[[l+1]]<-stack
        total_list[[l+1]]<-total
        names(stack_list)[l+1]<-paste("variant",variant[v], "dist=", distance[i], "threshold=" ,threshold[t],sep = "_")



      }
    }
  }
  #find all combis
  all_combi<-unique(all_combi)
  a<-all_combi
  a<-str_remove_all(a,pattern = "_uncalibrated")
  a<-str_remove_all(a,pattern = "_isotonic")
  a<-str_remove_all(a,pattern = "_sigmoid")
  a<-str_remove_all(a,pattern = "RF_")
  a<-str_remove_all(a,pattern = "ACE2_")
  a<-unique(a)
  a[which(a=="")]<-"None"

  return(a)
}
