

#' this function returns a list of igraph network depicting mutational relation between variants of different mutational distances.
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param thres_range numeric vector of threshold determine the upper and lower limit
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")
#' @return this function returns a list of igraph network depicting mutational relation between variants of different mutational distances.
escape_network<-function(data, distance=c(1,2), threshold=0.5,thres_range=c(0.25,0.75), model=c("RF","RNN"), variant=c("alpha","beta","gamma","kappa","wuhan","new")){
  igraph_list<-list()
  data_dis_result<-list()

  for (v in 1:length(variant) ){

    data[[v]]$No<-1:nrow(data[[v]])
    range<-1:ncol(data[[v]])
    for(j in 1:length(distance)){

      #row selection + partial col selection
      data_dis[[j]] <-subset(x=data[[v]],subset = data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = range[stringr::str_detect(colnames(data[[v]]),pattern="RF")|stringr::str_detect(colnames(data[[v]]),pattern="RNN")])

      data_dis_result[[j]]<-subset(x=data[[v]],subset = data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern="ist")]==distance[j],select = c(1,13,14))
      #col selection model
      for(k in 1:length(model)){

        modelselect<-stringr::str_detect(string = colnames(data_dis[[j]]),pattern = model[k])
        data_model[[k]]<-subset(x=data_dis[[j]],select = modelselect)

      }
      data_all[[j]]<-data_model
    }

    box<-data_dis_result

    for (t in 1:length(threshold)){

      data_dis_result<-box

      for( i in 1:length(distance)){


        #note ab escaping: reverse binding value, value means escape. T means yes


        bi0<-binary(th=threshold[t],data1=data_all[[i]][[1]],data2=data_all[[i]][[2]],disagree = 0)
        bi<-binary(th=c(0.25,0.75),data1=data_all[[i]][[1]],data2=data_all[[i]][[2]])

        bi<-bi[bi0[,1]==1,] ##filter out ace non-binding

        data_dis_result[[i]]<- data_dis_result[[i]][bi0[,1]==1,]


        bi<-as.data.frame(bi)
        n<-apply(bi[,-1],FUN = sum0,MARGIN = 1)

        #stack------
        # combi<-c()
        #
        # for(row in 1:nrow(bi)) combi[row]<-paste(colnames(bi)[which(bi[row,]==1)],collapse = "_")
        #
        # stack<-data.frame(bi,n,combi)
        # stack<-stack[bi[,1]==0,] ##filter out ace non-binding
        # stack<-dplyr::summarise(dplyr::group_by(stack,n,combi),count = dplyr::n())
        #
        # stack$combi[stack$combi==""]<-"None"
        #
        # l<-length(stack_list)
        # stack_list[[l+1]]<-stack
        # names(stack_list)[l+1]<-paste("variant",variant[v], "dist=", distance[i], "threshold=" ,threshold[t],sep = "_")
        #result----
        data_dis_result[[i]]$n<-n*10000
        data_dis_result[[i]]$ace<-bi0[bi0[,1]==1,]
        #filter here( 出错原因不明，改成在上面filter了)
        #data_dis_result[[i]]<-data_dis_result[[i]][bi[,1]==0,]

      }

      igraph_index<-c()
      #for root- dist1 - dist2 link
      for (r in 1:nrow(data_dis_result[[1]])){

        # a<-paste0(stringr::str_split(data_dis_result[[1]]$variant[r],pattern = ",",simplify = T),collapse = "_")
        #
        # b<-stringr::str_split(data_dis_result[[2]]$variant,pattern = ",",simplify = T)
        #
        # c<-apply(X = b,FUN = fun, MARGIN =1 )
        #
        # x<-data_dis_result[[1]]$No[r]
        #
        # y<-data_dis_result[[2]]$No[stringr::str_detect(string = c, pattern = a)]
        a<-data_dis_result[[1]]$No[r]
        b<-data_dis_result[[2]]$No[stringr::str_detect(data_dis_result[[2]]$mutation,data_dis_result[[1]]$mutation[r])]

        if(length(b)==0) vec<-c() else
          vec<-c(as.vector(t(as.matrix(data.frame(a,b)))))
        #vec<-c(data_dis_result[[1]]$n[r],a, as.vector(t(as.matrix(data.frame(a,b)))))

        igraph_index<-c(igraph_index,vec)
      }
      #for root-dist2
      # for(r in 1:nrow(data_dis_result[[2]])){
      #
      #   vec<-c(data_dis_result[[2]]$n[r], data_dis_result[[2]]$No[r])
      #   igraph_index<-c(igraph_index,vec)
      # }


      label<-unique(igraph_index)

      size<-rep(15,length(label))
      size[label>=10000]<-30

      No<-match(igraph_index,label) #reindexing

      result_all<- rbind(data_dis_result[[1]],data_dis_result[[2]])
      #igraph_n<-result_all$n[match(igraph_index,result_all$No)]


      igraph_index_df<-data.frame(No=igraph_index)
      igraph_index_df$col<-1:nrow(igraph_index_df)
      igraph_index_df<-merge(igraph_index_df,result_all,by="No",all.x = T)
      igraph_index_df<-igraph_index_df[order(igraph_index_df$col),]
      igraph_n<-igraph_index_df$n
      #igraph_n igraph_index --> unique format
      igraph_n<-igraph_n[match(label,igraph_index)]



      color<-rep("white",length(label))
      #color for root
      # color[label==0]<-"grey30"
      # color[label==10000]<-"royalblue4"
      # color[label==20000]<-"darkgreen"
      # color[label==30000]<-"gold4"
      # color[label==40000]<-"darkred"

      #color for sequence node
      color[igraph_n==0]<-"grey80"
      color[igraph_n==10000]<-"royalblue1"
      color[igraph_n==20000]<-"green"
      color[igraph_n==30000]<-"gold"
      color[igraph_n==40000]<-"red"


      # color[which(label<=cut)]<-"grey80"
      # color[which(label==10000)+1]<-"royalblue1"
      # color[which(label==20000)+1]<-"green"
      # color[which(label==30000)+1]<-"gold"
      # color[which(label==40000)+1]<-"red"

      # cut<-max(data_dis_result[[1]]$No)
      tag<-label
      # tag[label>cut&label<10000]<-"2"
      # tag[label<=cut&label>0]<-"1"
      # tag[label==10000]<-"1mAb"
      # tag[label==20000]<-"2mAb"
      # tag[label==30000]<-"3mAb"
      # tag[label==40000]<-"4mAb"

      tag[label %in% data_dis_result[[1]]$No]<-1
      tag[label %in% data_dis_result[[2]]$No]<-2


      g<-igraph::graph(No)

      igraph::V(g)$label<-tag
      #igraph::V(g)$shape<-'circle'
      igraph::V(g)$size<-size
      igraph::V(g)$color<-color

      name<-paste(variant[v],threshold[t],sep = "_")
      l<-length(igraph_list)
      igraph_list[[l+1]]<-g
      names(igraph_list)[l+1]<-name

    }
  }
  return(igraph_list)
}


