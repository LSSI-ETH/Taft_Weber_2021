
#' this function generates 2 series of pie charts: pie3 and pie4.
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param thres_range numeric vector of threshold determine the upper and lower limit
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")
#' @return pdf images in working directory. pie3: piecharts, category: "escape" both model,"unsure","binding","disagree". pie4: specify what are the different decision made by the 2 models in case of disagreement and unsure.

pred_value_pie<-function(data, distance=c(1,2), threshold=0.5,thres_range=c(0.25,0.75), model=c("RF","RNN"), variant=c("alpha","beta","gamma","kappa","wuhan","new")){
  pie3_list<-list()
  pie4_list<-list()
  ace_list<-list()
  total_list<-list()
  stack_list<-list()
  all_combi<-c()

  for (v in 1:length(variant) ){
    range<-1:ncol(data[[v]])
    data_dis<-list()
    data_model<-list()
    data_all<-list()
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


        l<-length(pie3_list)

        ace2 <-.binary_pie(th=threshold[t],data1=data_all[[i]][[1]],data2=data_all[[i]][[2]],specific = F,ab_split = T)[,1]
        ace2<-c(ace2[1],ace2[2]+ace2[3])
        names(ace2)<-c("binding","non-binding")
        bi<-.binary(th=threshold[t],data1=data_all[[i]][[1]],data2=data_all[[i]][[2]], disagree = -1)

        a2b<-bi[,1]==1

        mab <-.binary_pie(th=thres_range, data1=data_all[[i]][[1]][a2b,], data2=data_all[[i]][[2]][a2b,], specific = F,ab_split = T)
        ace_list[[l+1]]<-ace2
        pie3_list[[l+1]]<-mab[,-1]
        pie4_list[[l+1]]<-.binary_pie(th=thres_range,data1=data_all[[i]][[1]],data2=data_all[[i]][[2]],specific = T,ab_split = T)[,-1]

        names(pie3_list)[l+1]<-names(pie4_list)[l+1]<-paste(variant[v] ,"dist",distance[i],"n",nrow(data_all[[i]][[1]][a2b,]),sep="_")
        names(ace_list)[l+1]<-paste(variant[v] ,"dist",distance[i],"n",nrow(data_all[[i]][[1]][a2b,]),sep="_")

      }
    }
  }

  for(l in 1:length(pie3_list)){




    pdf(file=paste("pie3",names(pie3_list)[l],".pdf",sep = "_"),
        width = 10,
        height = 40
    )
    layout(matrix(c(1,2,3,4),nr=4,byrow=T))

    #pie(pie3_list[[l]][,1], labels = round(100*pie3_list[[l]][,1]/sum(pie3_list[[l]][,1]), 1),main = paste("ACE2"), col = c("#4575B4","#D73027"))
    # pie(pie3_list[[l]][,2], labels = round(100*pie3_list[[l]][,2]/sum(pie3_list[[l]][,2]), 1), main = paste("LY16"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    # pie(pie3_list[[l]][,3], labels = round(100*pie3_list[[l]][,3]/sum(pie3_list[[l]][,3]), 1), main = paste("LY555"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    # pie(pie3_list[[l]][,4], labels = round(100*pie3_list[[l]][,4]/sum(pie3_list[[l]][,4]), 1), main = paste("REGN33"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    # pie(pie3_list[[l]][,5], labels = round(100*pie3_list[[l]][,5]/sum(pie3_list[[l]][,5]), 1), main = paste("REGN87"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    #
    # nab<-1
    # pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), main = paste("LY16"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    #
    # nab<-nab+1
    # pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), main = paste("LY555"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    #
    # nab<-nab+1
    # pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), main = paste("REGN33"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))
    #
    # nab<-nab+1
    # pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), main = paste("REGN87"), col = c("#4575B4","#FFFFBF","#D73027","#999999"))


    nab<-1
    pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), col = c("#4575B4","#FFFFBF","#D73027","#999999"),cex=5)

    nab<-nab+1
    pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), col = c("#4575B4","#FFFFBF","#D73027","#999999"),cex=5)

    nab<-nab+1
    pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), col = c("#4575B4","#FFFFBF","#D73027","#999999"),cex=5)

    nab<-nab+1
    pie(pie3_list[[l]][,nab], labels = round(100*pie3_list[[l]][,nab]/sum(pie3_list[[l]][,nab]), 1), col = c("#4575B4","#FFFFBF","#D73027","#999999"),cex=5)

    dev.off()

    pdf(file=paste("ace",names(ace_list)[l],".pdf",sep = "_"),
        width = 10,
        height = 10
    )

    pie(ace_list[[l]], labels = round(100*ace_list[[l]]/sum(ace_list[[l]]), 1), col = c("#D73027","#4575B4"),cex=5)

    dev.off()

    pdf(file=paste("pie4",names(pie4_list)[l],".pdf",sep = "_"),
        width = 10,
        height = 40
    )

    layout(matrix(c(1,2,3,4),nr=4,byrow=T))

    pie(pie4_list[[l]][,1], labels = round(100*pie4_list[[l]][,1]/sum(pie4_list[[l]][,1]), 1), main = paste(names(pie4_list)[l],"LY16"), col =RColorBrewer::brewer.pal(9,"Set3"))

    pie(pie4_list[[l]][,2], labels = round(100*pie4_list[[l]][,2]/sum(pie4_list[[l]][,2]), 1), main = paste(names(pie4_list)[l],"LY555"), col = RColorBrewer::brewer.pal(9,"Set3"))

    pie(pie4_list[[l]][,3], labels = round(100*pie4_list[[l]][,3]/sum(pie4_list[[l]][,3]), 1), main = paste(names(pie4_list)[l],"REGN33"), col = RColorBrewer::brewer.pal(9,"Set3"))

    pie(pie4_list[[l]][,4], labels = round(100*pie4_list[[l]][,4]/sum(pie4_list[[l]][,4]), 1), main = paste(names(pie4_list)[l],"REGN87"), col =RColorBrewer::brewer.pal(9,"Set3"))

    dev.off()

  }

  pdf(file=paste("pie3_legend.pdf",sep = "_"),
      width = 10,
      height = 10
  )
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("bottom", rownames(pie3_list[[l]]),cex = 0.8, fill =  c("#4575B4","#FFFFBF","#D73027","#999999") )
  dev.off()

  pdf(file=paste("ace_legend.pdf",sep = "_"),
      width = 10,
      height = 10
  )
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("bottom", names(ace_list[[l]]),cex = 0.8, fill  = c("#D73027","#4575B4") )
  dev.off()

  pdf(file=paste("pie4_legend.pdf",sep = "_"),
      width = 10,
      height = 10
  )
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("bottom", rownames(pie4_list[[l]]),cex = 0.8, fill = RColorBrewer::brewer.pal(9,"Set3"))

  dev.off()



}
