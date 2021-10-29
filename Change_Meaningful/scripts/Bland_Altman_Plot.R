Bland_Altman_Raw <- function(S1, S2, S1_is_true_model=FALSE, center=TRUE ){
  
  x<-c()
  y<-c()
  
  if(S1_is_true_model==FALSE){
  for(i in 1:length(S1))
  {
    x<-c(x,(S1[i]+S2[i])/2)
    y<-c(y,(S1[i]-S2[i]))
  }
    ylab="Score Difference"
    xlab="Score Mean"
    }
  else{
   for(i in 1:length(S1))
    {
      x<-c(x,S1[i])
      y<-c(y,(S1[i]-S2[i]))
   } 
    ylab="Score Error" 
    xlab="True Score"
  }
  
  yscale<-scale(y, center=center)
  upper_limit<-min(y[which(yscale> 1.96)])
  lower_limit<-max(y[which(yscale< -1.96)])
  if(is.numeric(center)){mean<-center}
  else{mean<-mean(y, na.rm=TRUE)}
  
  BA<-list("XY"=data.frame("x"=x,"y"=y), "ylab"=ylab, "xlab"=xlab, 
           "mean"=mean,
           "upper limit"=upper_limit,
           "lower limit"=lower_limit)
  
  return(BA)
}

Bland_Altman_Plot <- function(BA, density=FALSE){

  if(density==FALSE){
    plot(BA$XY, ylab=BA$ylab, xlab=BA$xlab, main="Bland-Altman Plot")
    abline(h=0)
    abline(h=BA$mean, col="blue")
    abline(h=BA$`lower limit`, col="red")
    abline(h=BA$`upper limit`, col="red")
  }
  else{
    library(ggplot2)
    BA_plot<-ggplot(BA$XY,aes(x=x,y=y)) + geom_point(alpha = 0.3)
    BA_plot + geom_hline(yintercept=c(0,BA$mean,BA$lower,BA$upper), color=c("black", "blue","red","red")) + 
    labs(title="Bland-Altman Plot", y=BA$ylab, x=BA$xlab)
    
  }
    
}


