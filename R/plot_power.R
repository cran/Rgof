#' This function draws the power graph, with curves sorted by the mean power and smoothed for easier reading.
#' @param  pwr  a matrix of power values, usually from gof_power or hybrid_power 
#' @param  xname Name of variable on x axis
#' @param  title (Optional) title of graph
#' @param  Smooth =TRUE lines are smoothed for easier reading. Smoothing needs enough
#'         alternatives; with fewer than 10 values of the parameter the powers are
#'         joined by straight lines instead.
#' @param  span  =0.75 span of the loess smoother (as a fraction of the points)
#' @return plt, an object of class ggplot.
#' @export 

plot_power=function(pwr, xname=" ", title, Smooth=TRUE, span=0.75) {

  if(inherits(pwr, "Rgof_power")) pwr <- pwr$power
  if(!is.matrix(pwr))
    stop("plot_power requires power estimates for at least two alternative parameter values",
         call.=FALSE)

# For CRAN CMD check
  x=NULL
  y=NULL
  Method=NULL
# sort methods by their average power
  mu=apply(pwr, 2, mean) 
  lvls = colnames(pwr)[order(mu, decreasing = TRUE)]
# create data frame  
  df=data.frame(x=rep(as.numeric(rownames(pwr)), ncol(pwr)),
       y=100*c(pwr),
       Method=factor(rep(colnames(pwr), each=nrow(pwr)),  
                     levels=lvls,
                     ordered = TRUE))
# colors: the standard palette has only 8 colors, with more methods use distinct ones
  cols = if(length(mu)<=8) seq_along(mu) else grDevices::hcl.colors(length(mu), "Dark 3")
# create ggplot graphics                     
  plt=ggplot2::ggplot(data=df, ggplot2::aes(x=x, y=y, color=Method))+
     ggplot2::xlab(xname)+
     ggplot2::ylab("Power")+
     ggplot2::scale_color_manual(labels=lvls,
                                 values=cols,
                                 name='Method')
  if(!missing(title)) 
     plt=plt+ggplot2::ggtitle(title)
  # loess needs enough points per curve, with only a few alternatives it fails
  # and the graph is empty
  if(Smooth && nrow(pwr)>=10) plt=plt+
    ggplot2::geom_smooth(formula = y ~ x, method="loess", se=FALSE, span=span)
  else  plt=plt+ggplot2::geom_line()+ggplot2::geom_point()
  plt
}
