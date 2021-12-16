#' this function generates stacked bar plots for number of antibody escape in working directory.
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param thres_range numeric vector of threshold determine the upper and lower limit
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")
#' @param colr a named character vector indicates the color of the bar plot.
#' @param n_antibody total number of antibody in data, default is 4.
#' @return this function generates stacked bar plots for number of antibody escape in working directory. the legend is exported in separate file.
escape_stackedbar<-function(data, n_antibody=4, colr=  c(None="#999999",LY16="#8cb0e6",LY555="#205cb6",REGN33="#79ccec",REGN87="#1c366a",LY16_LY555="#338104",LY16_REGN33="#c5ff60",LY16_REGN87="#4d5026",LY555_REGN33="#88fd41",REGN33_REGN87="#49836e",LY555_REGN87="#fae80b",LY16_LY555_REGN33="#ab6624",LY16_LY555_REGN87="#f6b50a",LY16_REGN33_REGN87="#673900",LY555_REGN33_REGN87="#ff6900",LY16_LY555_REGN33_REGN87="#da0000")
,distance=c(1,2), threshold=0.5,thres_range=c(0.25,0.75), model=c("RF","RNN"), variant=c("alpha","beta","gamma","kappa","wuhan","new")){

  stack_list<-list()
  all_combi<-c()
  total_list<-list()
  nAbEscape<-0:n_antibody

  for (v in 1:length(variant) ){
    range<-1:ncol(data[[v]])
    for(j in 1:length(distance)){

      #row selection + partial col selection
      data_dis[[j]] <-subset(x=data[[v]],subset = data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = range[stringr::str_detect(colnames(data[[v]]),pattern="RF")|stringr::str_detect(colnames(data[[v]]),pattern="RNN")])

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
  a<-stringr::str_remove_all(a,pattern = "_uncalibrated")
  a<-stringr::str_remove_all(a,pattern = "_isotonic")
  a<-stringr::str_remove_all(a,pattern = "_sigmoid")
  a<-stringr::str_remove_all(a,pattern = "RF_")
  a<-stringr::str_remove_all(a,pattern = "ACE2_")
  a<-unique(a)
  a[which(a=="")]<-"None"
  #a %in% names(colr)
  colr_s<-colr[names(colr) %in% a]
  for(l in 1:length(stack_list)){
    a<-stack_list[[l]]$combi
    a<-stringr::str_remove_all(a,pattern = "_uncalibrated")
    a<-stringr::str_remove_all(a,pattern = "_isotonic")
    a<-stringr::str_remove_all(a,pattern = "_sigmoid")
    a<-stringr::str_remove_all(a,pattern = "RF_")
    a<-stringr::str_remove_all(a,pattern = "ACE2_")
    stack_list[[l]]$combi<-a

    p<-ggplot2::ggplot(stack_list[[l]], aes(fill=combi, x=n,y=count)) +
      ggplot2::geom_bar(position="stack", stat="identity",width = 0.7)+
      ggplot2:: xlab("n of antibody escape") + ggplot2::ylab("sequences count") +
      ggplot2::scale_fill_manual(values=colr)+
      ggplot2::theme_classic()+
      ggplot2::theme(legend.position="none",text=ggplot2::element_text(size=20),axis.text= ggplot2::element_text(size=30), axis.title=ggplot2::element_text(size=30), plot.margin =ggplot2:: margin(1,1,1,1, "cm"),)+
      ggplot2::geom_text(ggplot2::aes(n, total, label = total, fill = NULL, vjust = -0.5), size=10, data = total_list[[l]])


    ggplot2::ggsave(filename =  paste0("stack_",names(stack_list)[l],".pdf"), path = wd, device = "pdf", plot =p,height = 12,width =12)
  }

  p<-ggplot2::ggplot(stack_list[[l]], ggplot2::aes(fill=combi, x=n,y=count)) +
    ggplot2::geom_bar(position="stack", stat="identity",width = 0.7)+
    ggplot2::xlab("n of antibody escape") + ggplot2::ylab("sequences count") +
    ggplot2::scale_fill_manual(values=colr_s)+
    ggplot2::theme_classic()+
    ggplot2::theme(text=ggplot2::element_text(size=20),axis.text= ggplot2::element_text(size=30), axis.title=ggplot2::element_text(size=30), plot.margin = ggplot2::margin(1,1,1,1, "cm"),)+
    ggplot2::geom_text(ggplot2::aes(n, total, label = total, fill = NULL, vjust = -0.5),size=10, data = total_list[[l]])
  ggplot2::ggsave(filename =  paste0("stack_legend_",names(stack_list)[l],".pdf"), path = wd, device = "pdf", plot =p,height = 8,width = 8)
}









