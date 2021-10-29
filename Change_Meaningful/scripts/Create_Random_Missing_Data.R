Create_Random_Missing_Data<-function(Data, Variables, number_of_missing, mean_replacement=FALSE)
{
  if(mean_replacement==FALSE){
  for(i in 1:nrow(Data))
    {
    variables_to_remove<-c(which((names(Data) %in% Variables[sample(1:length(Variables), number_of_missing, replace=F)])))
    Data[i,variables_to_remove]<-NA
  }}
  
  if(mean_replacement==TRUE){
  for(i in 1:nrow(Data))
  {
    variables_to_remove<-c(which((names(Data) %in% Variables[sample(1:length(Variables), number_of_missing, replace=F)])))
    Data[i,variables_to_remove]<-round(mean(as.numeric(Data[i,Variables],na.rm=TRUE)))
  }}
  
  return(Data)
}