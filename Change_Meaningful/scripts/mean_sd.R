mean_sd<-function(Data, na.rm=TRUE)
{return(paste0(round(mean(Data,na.rm=na.rm),2),"(",round(sd(Data,na.rm=na.rm),2),")"))}

mean_sd_n<-function(Data, na.rm=TRUE)
{return(paste0(round(mean(Data,na.rm=na.rm),2),"(",round(sd(Data,na.rm=na.rm),2),"), n=",length(!is.na(Data))))}