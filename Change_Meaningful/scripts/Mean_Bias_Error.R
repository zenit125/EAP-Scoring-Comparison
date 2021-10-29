Mean_Bias_Error <- function(Model_data, Observed_data){
  Model_data_trim<-Model_data[which(!is.na(Model_data)&!is.na(Observed_data))]
  Observed_data_trim<-Observed_data[which(!is.na(Model_data)&!is.na(Observed_data))]
  Bias_Residuals<-abs(Observed_data_trim-Model_data_trim)/Model_data_trim
  return(sum(Bias_Residuals)/length(Model_data_trim))
}