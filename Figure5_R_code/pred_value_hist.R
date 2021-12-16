
#' this function draw histograms of predicted value.
#' @param data output of omni_read function, list of data separated by variant and their calibration method in prediction models
#' @param distance numeric vector of mutational distances from root variant sequence, which are wanted to be included in the analysis.
#' @param threshold numeric vector of threshold, used for ACE2 binding screening. predicted value > threshold -> binding. every threshold will result in a full set of plots.
#' @param model character vector of model names. default c("RF","RNN")
#' @param variant character vector of variant names. default c("alpha","beta","gamma","kappa","wuhan","new")
pred_value_hist<-function(data, variant=c("alpha","beta","gamma","kappa","wuhan","new"),model=c("RF","RNN"),distance=c(1,2)){
  for(m in 1:length(model)){
    for(v in 1:length(variant)){
      range<-1:ncol(data[[v]])
      pdf(file=paste0("hist_" , paste0(model[m],"_"), paste(variant[v],sep = "_"), ".pdf"),
          width = 10,
          height = 5
      )
      layout(matrix(c(1,2,3,6,4,5),nr=2,byrow=T))
      for (j in range[str_detect(colnames(data[[v]]),pattern=model[m])]){
        a<-stringr::str_replace_all(colnames(data[[v]])[j],pattern = paste0(model[m],"_"),replacement = "")
        subdist<-which(data[[v]]$dist_aa %in% distance)
        hist(data[[v]][subdist,j],breaks=20, main = paste(a,sep = "_"),xlab = "predicted value")
      }
      dev.off()

    }
  }
}
