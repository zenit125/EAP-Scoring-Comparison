Generate_Random_Data_from_Response_Probabilities<-function(N, most_probable_response, mean=50, sd=10)
{
  Tscore<-round(rnorm(N, mean, sd))
  Tscore[which(Tscore<10)]<-10
  Tscore[which(Tscore>90)]<-90
  random_data<-data.frame(matrix(NA,ncol=length(most_probable_response),nrow=length(Tscore)))
  for(i in 1:length(Tscore))
  {
    random_data[i,]
    for(j in 1:length(most_probable_response))
    {
      random_data[i,j]<-as.numeric(sample(
        names(most_probable_response[[j]]),
        size=1,
        prob=most_probable_response[[j]][which(rownames(most_probable_response[[j]])==Tscore[i]),]))
    }
  }
  names(random_data)<-names(most_probable_response)
  random_data$Tscore<-Tscore
  return(random_data)
}