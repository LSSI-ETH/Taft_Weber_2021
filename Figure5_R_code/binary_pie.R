.binary_pie<-function(th, data1, data2,specific=F,ab_split=F){
  
  if(length(th)==1){
    
    bi1<-data1>th
    bi2<-data2>th
    di<-c()
    
    if(specific==T){
      di[1]<-sum(bi1&bi2)
      di[2]<- sum(bi1==T&bi2==F)
      di[3]<- sum(bi1==F&bi2==T)
      di[4]<-sum(!(bi1&bi2))
      di<-data.frame(group=c("binding","RF=T","RNN=T","escape"),value=di)
      
      ma<-matrix(0,4,ncol(data1))
      colnames(ma)<-colnames(data1)
      rownames(ma)<-c("binding","RF=T","RNN=T","escape")
      
      for(i in 1:ncol(data1)){
        unit<-c()
        bi1<-data1[,i]>th
        bi2<-data2[,i]>th
        unit[1]<-sum(bi1&bi2)
        unit[2]<- sum(bi1==T&bi2==F)
        unit[3]<- sum(bi1==F&bi2==T)
        unit[4]<-sum(!(bi1&bi2))
        ma[,i] <-unit
      }
    }
    if(specific==F){
      di[1]<-sum(bi1&bi2)
      di[2]<- sum(bi1==T&bi2==F)+sum(bi1==F&bi2==T)
      di[3]<- sum(!(bi1&bi2))

      di<-data.frame(group=c("binding","disagree","escape"),value=di)
      
      ma<-matrix(0,3,ncol(data1))
      colnames(ma)<-colnames(data1)
      rownames(ma)<-c("binding","disagree","escape")
      
      for(i in 1:ncol(data1)){
        unit<-c()
        bi1<-data1[,i]>th
        bi2<-data2[,i]>th
        unit[1]<-sum(bi1&bi2)
        unit[2]<- sum(bi1==T&bi2==F)+sum(bi1==F&bi2==T)
        unit[3]<- sum(!(bi1&bi2))
        ma[,i] <-unit
      }
    }

    
    
    
    
    
    
    
  }
  
  if(length(th)==2){
    xi1<-data1<th[1]
    xi2<-data2<th[1]
    da1<-data1>th[2]
    da2<-data2>th[2]
    mi1<-da1==F&xi1==F
    mi2<-da2==F&xi2==F
    
    
    you<-sum(da1&da2)
    zhong<-sum(mi1|mi2)
    zuo<-sum(xi1&xi2)
    cuo<-sum((da1&xi2)|(xi1&da2))
    
    group<-c("escape","unsure","binding","disagree")
    value<-c(zuo, zhong, you,cuo)
    di<-data.frame(group,value)
    
    ma<-matrix(0,4,ncol(data1))
    colnames(ma)<-colnames(data1)
    
    for(i in 1:ncol(data1)){
      xi1<-data1[,i]<th[1]
      xi2<-data2[,i]<th[1]
      da1<-data1[,i]>th[2]
      da2<-data2[,i]>th[2]
      mi1<-da1==F&xi1==F
      mi2<-da2==F&xi2==F
      you<-sum(da1&da2)
      zhong<-sum(mi1|mi2)
      zuo<-sum(xi1&xi2)
      cuo<-sum((da1&xi2)|(xi1&da2))
      ma[,i]<-c(zuo, zhong, you, cuo)
    }
    rownames(ma)<-group
    
    
    #more specific:
    if(specific==T){
      
      xi1<-data1<th[1]
      xi2<-data2<th[1]
      da1<-data1>th[2]
      da2<-data2>th[2]
      mi1<-da1==F&xi1==F
      mi2<-da2==F&xi2==F
      
      
      you<-sum(da1&da2)
      zhong<-sum(mi1|mi2)
      zuo<-sum(xi1&xi2)
      cuo<-sum((da1&xi2)|(xi1&da2))
      rnnt<-sum(xi1&da2)
      rft<-sum(da1&xi2)
      
      
      group<-c("negative","RF unsure RNN=T","RF unsure RNN=F","RNN unsure RF=T","RNN unsure RF=F","both unsure","positive","disagree RNN=T","disagree RF=T")
      value<-c(zuo,sum(mi1&xi2),sum(mi1&da2),sum(mi2&xi1),sum(mi2&da1),sum(mi1&mi2),you,rnnt,rft)
      
      di<-data.frame(group,value)
      ma<-matrix(0,9,ncol(data1))
      colnames(ma)<-colnames(data1)
      for (i in 1:ncol(data1)) {
        xi1<-data1[,i]<th[1]
        xi2<-data2[,i]<th[1]
        da1<-data1[,i]>th[2]
        da2<-data2[,i]>th[2]
        mi1<-da1==F&xi1==F
        mi2<-da2==F&xi2==F
        you<-sum(da1&da2)
        zhong<-sum(mi1|mi2)
        zuo<-sum(xi1&xi2)
        cuo<-sum((da1&xi2)|(xi1&da2))
        rnnt<-sum(xi1&da2)
        rft<-sum(da1&xi2)
        value<-c(zuo,sum(mi1&da2),sum(mi1&xi2),sum(mi2&da1),sum(mi2&xi1),sum(mi1&mi2),you,rnnt,rft)
        ma[,i]<-value
      }
      rownames(ma)<-group
    }
  }
  
  if(ab_split==F) output<-di
  if(ab_split==T) output<-ma
  return(output)
}