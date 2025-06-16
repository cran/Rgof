#' This function checks whether the correct methods have been requested
#' @param  doMethods ="all" Which methods should be included?
#' @param  Continuous is data continuous
#' @param  WithWeights with weights?
#' @keywords internal
#' @return TRUE or FALSE
#' @export 
test_methods=function(doMethods, Continuous, WithWeights) {
    if(doMethods[1]%in%c("all","default")) return(FALSE)
    if(Continuous & !WithWeights) 
       methods=c("KS","Kuiper","CvM","AD","W","ZA","ZK","ZC",
                 "Wassp1","ES-l-P","ES-s-P","EP-l-P","EP-s-P",
                          "ES-l-P","ES-s-P","EP-l-P","EP-s-P")
    if(Continuous & WithWeights) 
       methods=c("KS","Kuiper","CvM","AD")
    if(!Continuous) 
       methods=c("KS","Kuiper","CvM","AD","W", "Wassp1","l-P","s-P")
    Good=TRUE
    for(i in seq_along(doMethods)) {
      if(!(doMethods[i]%in%methods)) {Good=FALSE;break}
    }
    if(Good) return(FALSE)
    message(paste0(doMethods[i]," is not an included Method!"))
    if(Continuous & !WithWeights) {
         message("For continuous data without weights included methods are")
         message("Method               Code")
         message("Kolmogorov-Smirnov   KS")
         message("Kuiper               Kuiper")
         message("Cramer-vonMises      CvM")
         message("Anderson-Darling     AD")
         message("Watson               W")
         message("Zhang's tests        ZA, ZK and ZC")
         message("Wasserstein          Wassp1")
         message("Chi square tests     ES-l-P, ES-s-P, EP-l-P, EP-s-P")
         message("                     ES-l-L, ES-s-L, EP-l-L, EP-s-L")
    }
    if(Continuous & WithWeights) {
      message("For continuous data with weights included methods are")
      message("Method               Code")
      message("Kolmogorov-Smirnov   KS")
      message("Kuiper               Kuiper")
      message("Cramer-vonMises      CvM")
      message("Anderson-Darling     AD")
    }
    if(!Continuous) {
      message("For discrete data included methods are")
      message("Method               Code")
      message("Kolmogorov-Smirnov   KS")
      message("Kuiper               Kuiper")
      message("Cramer-vonMises      CvM")
      message("Anderson-Darling     AD")
      message("Watson               W")
      message("Wasserstein          Wassp1")
      message("Chi square tests     l-P, s-P, l-L, s-L")
    }
    TRUE
}
