Create_Random_Missing_Data_With_Mean_Replacement<-function(Data, Variables, number_of_missing)
{
  for(i in 1:nrow(Data))
    {
    variables_to_remove<-c(which((names(Data) %in% Variables[sample(1:length(Variables), number_of_missing, replace=F)])))
    Data[i,variables_to_remove]<-mean(as.numeric(Data[i,Variables],na.rm=TRUE))
    }
  return(Data)
}