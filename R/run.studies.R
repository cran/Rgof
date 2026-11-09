#' Power Comparisons
#'
#' This function runs the case studies included in the package and compares the
#' power of a new test to those included.
#'
#' For details consult vignette("Rgof","Rgof")
#'
#' @param studies =1:60, either the names of the studies, or their numbers. If missing all the studies are run.
#' @param types =1:5, either the types of studies, or their numbers. If missing all types are run.
#' @param data_type ="cont", continuous or discrete data
#' @param TS routine to calculate test statistics. If Hybrid=TRUE it has to have the
#'        format required by R2sample (a two-sample test).
#' @param TSextra list passed to TS.
#' @param With.p.value =FALSE does user supplied routine return p values?
#' @param n desired sample size.
#' @param alpha =0.05  type I error
#' @param B = 1000
#' @param SuppressMessages =TRUE
#' @param RerunIncludedCases =FALSE, set to true to verify powers in Rgof::power_study
#'        (or power_study_hybrid, power_study_hybrid_nMC5)
#' @param seed =NULL, optional seed for reproducibility. Case study i in the list of
#'        studies/types that are run uses the seed seed+i, both for the new test and
#'        (if they are recomputed) for the included methods.
#' @param maxProcessor number of processors used in gof_power or hybrid_power. If missing
#'        it is chosen automatically. Results with a given seed are reproducible only
#'        with the same number of processors.
#' @param Null =FALSE, if TRUE the data is generated under the null hypothesis (by rnull),
#'        so the result is the type I error.
#' @param Hybrid =FALSE, if TRUE the goodness-of-fit/two-sample hybrid tests
#'        (hybrid_power, requires package R2sample) are run instead of the
#'        goodness-of-fit tests.
#' @param nMC =1, size of the Monte Carlo data set for the hybrid tests, as a multiple of
#'        the sample size. Stored results exist for nMC=1 (power_study_hybrid) and
#'        nMC=5 (power_study_hybrid_nMC5).
#' @return A matrix of power values.
#' @examples
#' TSextra=list(nbins=5, statistic=FALSE) # Use 5 bins and calculate p values
#' run.studies(1:2, 1:2, TS=myTS_cont, TSextra=TSextra, With.p.value=TRUE, B=100)
#' @export
run.studies=function(studies=1:60, types=1:5, data_type="cont",
                     TS, TSextra, With.p.value =FALSE, n,
                     alpha=0.05, B=1000, SuppressMessages=TRUE,
                     RerunIncludedCases=FALSE, seed=NULL, maxProcessor,
                     Null=FALSE, Hybrid=FALSE, nMC=1) {

# check if old version is assumed
  if(is.function(studies)) {
    stop("the arguments list of the function run.studies has changed.\n
          Please check the help file or the vignette",
         call.=FALSE)
  }
  if(Hybrid) check_R2sample("run.studies(..., Hybrid=TRUE)")
  noMP=missing(maxProcessor)
  NewTest=ifelse(missing(TS), FALSE, TRUE)
  if(missing(TSextra)) TSextra=list(aaa=NULL)

# stored results: power_study, or for the hybrid tests power_study_hybrid (nMC=1)
# or power_study_hybrid_nMC5 (nMC=5), if they exist
  stored_name=if(!Hybrid) "power_study"
              else if(nMC==1) "power_study_hybrid"
              else if(nMC==5) "power_study_hybrid_nMC5"
              else ""
  stored_all=tryCatch(getExportedValue("Rgof", stored_name), error=function(e) NULL)

# power of the included methods (TS missing) or of a new test (TS given) for
# case study tmp
  find_power=function(tmp, TS, seed) {
    a=list(alpha=alpha, B=B, SuppressMessages=SuppressMessages, seed=seed)
    if(!noMP) a$maxProcessor=maxProcessor
    if(!missing(TS)) a=c(a, list(TS=TS, TSextra=TSextra, With.p.value=With.p.value))
    if(Hybrid) {
      a=c(a, list(rnull=tmp$rnull, ralt=tmp$ralt, param_alt=tmp$apar,
                  vals=tmp$vals, phat=tmp$phat, nMC=nMC))
      pwr=do.call(Rgof::hybrid_power, a)
    }
    else {
      a=c(a, list(list.with.everything=tmp, param_alt=tmp$apar))
      pwr=do.call(Rgof::gof_power, a)
    }
    if(is.matrix(pwr)) pwr=stats::setNames(pwr[1, ], colnames(pwr))
    # a user-supplied statistic may come without names
    if(is.null(names(pwr)))
      names(pwr)=if(length(pwr)==1) "NewTest" else paste0("NewTest", seq_along(pwr))
    pwr
  }

  list.of.studies=Rgof::case.studies(ReturnCaseNames = TRUE)
  if(is.numeric(studies[1])) studies=list.of.studies[studies]
  all_types=names(Rgof::power_study[[1]])
  if(is.numeric(types)) types=all_types[types]
  studytype=cbind(rep(studies, length(types)),
                      rep(types, each=length(studies)))
  k=NULL
  for(i in 1:nrow(studytype)) {
    if(studytype[i,1]%in%list.of.studies[1:20] && studytype[i,2]=="mixed")
       k=c(k, i)
    if(studytype[i,1]%in%list.of.studies[-c(1:20)]) {
       studytype[i,2]="mixed"
    }
  }
  if(!is.null(k)) studytype=studytype[-k, , drop=FALSE]
  studytype=unique(studytype)
  rows=paste0(studytype[,1], "/", ifelse(Null, "Null", studytype[,2]))
  Old_n=FALSE
  if(missing(n)) {
    Old_n=TRUE
    n=rep(0, nrow(studytype))
    for(i in 1:nrow(studytype)) {
       if(studytype[i, 2]=="mixed")
         n[i]=Rgof::case_studies_sample_sizes[[data_type]][[2]][studytype[i,1]]
       else
         n[i]=Rgof::case_studies_sample_sizes[[data_type]][[1]][studytype[i,1], studytype[i,2]]
    }
  }
  else n=rep(n, length.out=nrow(studytype))
  out_new_test=list()
  out_old_test=list()
  for(i in 1:nrow(studytype)) {
      tmp <- Rgof::case.studies(
        which = studytype[i,1],
        which_diff = studytype[i,2],
        n = n[i],
        data_type = data_type
      )
      if(data_type=="cont") tmp$vals=NA
      if(Null) {  # data generated under the null hypothesis
        rn=tmp$rnull
        tmp$ralt=function(a) rn()
      }
      s=if(is.null(seed)) NULL else seed+i
      if(NewTest) out_new_test[[i]]=find_power(tmp, TS, s)
      # stored results; under the null hypothesis in attr(, "Null"), if present
      stored=NULL
      if(!is.null(stored_all))
        stored=if(Null) attr(stored_all, "Null")[[data_type]]
               else stored_all[[data_type]][[studytype[i,2]]]
      if(Old_n&&alpha==0.05&&(!RerunIncludedCases)&&
         !is.null(stored)&&studytype[i, 1]%in%rownames(stored)&&
         !all(is.na(stored[studytype[i, 1], ]))) {
          out_old_test[[i]]=stored[studytype[i, 1], ]
      }
      else out_old_test[[i]]=find_power(tmp, seed=s)
  }
  # combine the results, methods are the columns
  as_matrix=function(l) {
    if(length(l)==0) return(NULL)
    nm=unique(unlist(lapply(l, names)))
    m=matrix(NA_real_, length(l), length(nm), dimnames=list(rows, nm))
    for(i in seq_along(l)) m[i, names(l[[i]])]=l[[i]]
    m
  }
  out = cbind(as_matrix(out_new_test), as_matrix(out_old_test))
  a1=apply(out, 1, rank)
  message("Average rank of a method:")
  print(round(100*sort(apply(a1,1,mean), decreasing = TRUE)/nrow(a1),1))
  out
}
