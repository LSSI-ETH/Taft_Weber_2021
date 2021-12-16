


#' this function read in prediction output from models and split data with different calibration methods into different list items.
#' @param data_path file path of the csv data.
#' @param dataframe_list name of the list of R dataframe, each dataframe should be of 1 variant
#' @param model a vector of model names, the names should be written in the column names with underscore connecting to the antibody name, such as RF_LY16
#' @param variant a vector of variant names.
#' @param cat vector of calibration method names. this name should also be written in the column name with underscore, such as RF_LY16_isotonic
#' @return returns a list of data spited by variant and different calivbration method.

omni_read<-function(data_path, dataframe_list, model=c("RF","RNN"), variant=c("alpha","beta","gamma","kappa","wuhan","new"), cat="isotonic"){
  data<-list()
  if(!missing(data_path)){
    data_path<-list.files(data_path,full.names = T)


    for(i in 1:length(data_path)){
      data[[i]]<-read.csv(file=data_path[i])
    }
  }

  if(missing(data_path)){
    data<-dataframe_list
  }


  data_spl<-list()

  for(v in 1:length(data)){

    common<-data[[v]][,stringr::str_detect(string = colnames(data[[v]]),pattern = "equence")|stringr::str_detect(string = colnames(data[[v]]),pattern = "dist")]

    RF_ACE2<-data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern = "RF_ACE2")]

    data_cali<-list()

    for(c in 1:length(cat)){
      data_cali[[c]]<-data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern = cat[c])]
    }

    RNN_ACE2<-data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern = "RNN_ACE2")]

    data_RNN<-data[[v]][,stringr::str_detect(colnames(data[[v]]),pattern = "RNN")&!(stringr::str_detect(colnames(data[[v]]),pattern = "RNN_ACE2"))]

    mutation<-data[[v]][,stringr::str_detect(string = colnames(data[[v]]),pattern = "mutation")]

    for (i in 1:length(data_cali)){
      #colnames(data_cali[[i]])<-colnames(data_cali[[i]])
      l<-length(data_spl)
      data_spl[[l+1]]<- cbind(common,RF_ACE2,data_cali[[i]],RNN_ACE2,data_RNN,mutation)
      names(data_spl)[l+1]<-paste(variant[v],cat[i])
    }

  }


  data<-data_spl
  variant<-names(data)

  for(v in 1:length(data)){
    a<-colnames(data[[v]])
    a<-stringr::str_remove_all(a,pattern = "_uncalibrated")
    a<-stringr::str_remove_all(a,pattern = "_isotonic")
    a<-stringr::str_remove_all(a,pattern = "_sigmoid")
    colnames(data[[v]])<-a
  }
  return(data)
}







