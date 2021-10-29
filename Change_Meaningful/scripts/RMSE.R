RMSE <- function(Model_data, Observed_data, round_digits=0){

  Model_Observed_Error<-c()

      Model_Unique<-unique(Model_data)
      RMSE_Error_by_Tscore<-c(rep(NA,times=length(Model_Unique)))
      names(RMSE_Error_by_Tscore)<-Model_Unique

      for(Tscore in Model_Unique)
      {
        Model_Observed_Error_by_Tscore<-(Tscore-Observed_data[which(Model_data==Tscore)])
        RMSE_Error_by_Tscore[which(names(RMSE_Error_by_Tscore)==Tscore)]<-sqrt(mean(Model_Observed_Error_by_Tscore^2, na.rm=TRUE))
        Model_Observed_Error<-c(Model_Observed_Error,Model_Observed_Error_by_Tscore)
      }
      
      Scale_RMSE<-sqrt(mean(Model_Observed_Error^2, na.rm = TRUE))
      RMSE_Error_by_Tscore<-RMSE_Error_by_Tscore[order(names(RMSE_Error_by_Tscore))]
      
      if(round_digits==0){
        Return_Data<-list(
          "Scale level RMSE"=Scale_RMSE,
          "RMSE by Tscore list"=RMSE_Error_by_Tscore,
          "RMSE by Tscore range"=paste(range(RMSE_Error_by_Tscore),collapse="-"),
          "RMSE by Tscore SD"=sd(RMSE_Error_by_Tscore),
          "Summary"=paste0(Scale_RMSE,
                           "(",sd(RMSE_Error_by_Tscore),") ",
                           paste(range(RMSE_Error_by_Tscore),collapse="-")))
        }
      else{
      Return_Data<-list(
        "Scale level RMSE"=round(Scale_RMSE,round_digits),
        "RMSE by Tscore list"=round(RMSE_Error_by_Tscore,round_digits),
        "RMSE by Tscore range"=paste(range(round(RMSE_Error_by_Tscore,round_digits) ),collapse="-"),
        "RMSE by Tscore SD"=round(sd(RMSE_Error_by_Tscore),round_digits),
        "Summary"=paste0(round(Scale_RMSE,round_digits),
                         "(",round(sd(RMSE_Error_by_Tscore), round_digits),") ",
                         paste(round(range(RMSE_Error_by_Tscore),round_digits),collapse="-"))
        
        )
        }
  
  
  return(Return_Data)
}