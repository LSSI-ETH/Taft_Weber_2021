
#input 2 dataframe of 2 model, col should have the same order, return binary decision with 2 model agreement, disaggrement can be set to same as false or other value.
# th threshold
# disagree  the value to return

.binary<-function(th, data1, data2, disagree=-1){
  if(length(th)==1){
    bi1<-data1>th
    bi2<-data2>th
    
    c<-bi1&bi2
    c[bi1!=bi2]<-disagree
    
    bi<-matrix(as.numeric(c),dim(bi1))
    
    
    colnames(bi)<-colnames(bi1)
    rownames(bi)<-rownames(bi1)
  }
  
  if(length(th)==2){
    zuo1<-data1<th[1]
    zuo2<-data2<th[1]
    you1<-data1>th[2]
    you2<-data2>th[2]
    
    bi_matrix3<-list()
    
    #negative
    bi_matrix3[[1]]<-matrix(as.numeric(zuo1&zuo2),dim(you1))
    #positive
    bi_matrix3[[2]]<-matrix(as.numeric(you1&you2),dim(you1))
    #disagree
    bi_matrix3[[3]]<-matrix(as.numeric((zuo1&you2)|(you1&zuo2)),dim(you1))
    #unsure
    bi_matrix3[[4]]<-!(bi_matrix3[[1]]|bi_matrix3[[2]]|bi_matrix3[[3]])
    
    num<-c(0,1,-2,-1)
    
    bi<-matrix(data = 0,nrow = nrow(you1),ncol = ncol(you1))
    
    for(i in 1:length(bi_matrix3)){
      bi_matrix3[[i]][which(bi_matrix3[[i]]==1)]<-num[i]
      bi<-bi+bi_matrix3[[i]]
    }
    colnames(bi)<-colnames(data1)
    rownames(bi)<-rownames(data1)
    
  }

  return(bi)
}
