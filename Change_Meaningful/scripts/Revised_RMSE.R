
sqrt(mean((m - o)^2))

Model_data<-Scored_data$Sim$Depression$`True TScore`
Observed_data<-Scored_data$Sim$Depression$RSSS_scoring_dfs$TScore$SF_4a
Observed_data<-Scored_data$Missing$Depression$`1`$RSSS_scoring_dfs$TScore$SF_4a
Observed_data<-Scored_data$Missing$Depression$`2`$RSSS_scoring_dfs$TScore$SF_4a
Observed_data<-Scored_data$Missing$Depression$`3`$RSSS_scoring_dfs$TScore$SF_4a

Observed_data<-Scored_data$Sim$Depression$ThetaSEeap_scoring_dfs$TScore$SF_4a
Observed_data<-Scored_data$Missing$Depression$`1`$ThetaSEeap_scoring_dfs$TScore$SF_4a
Observed_data<-Scored_data$Missing$Depression$`2`$ThetaSEeap_scoring_dfs$TScore$SF_4a
Observed_data<-Scored_data$Missing$Depression$`3`$ThetaSEeap_scoring_dfs$TScore$SF_4a

Model_Observed_Error<-c()
RMSE_Error_by_Tscore<-c()

for(Tscore in min(Model_data):max(Model_data))
{
  Model_Observed_Error_by_Tscore<-(Tscore-Observed_data[which(Model_data==Tscore)])
  RMSE_Error_by_Tscore[Tscore]<-
  sqrt(mean(Model_Observed_Error_by_Tscore^2, na.rm=TRUE))
  Model_Observed_Error<-c(Model_Observed_Error,Model_Observed_Error_by_Tscore)
}

plot(RMSE_Error_by_Tscore, xlab="T score", ylab="RMSE Error", main="Depression 4a \n ThetaSEeap Scoring")
points(RMSE_Error_by_Tscore, col="red")
points(RMSE_Error_by_Tscore, col="green")
points(RMSE_Error_by_Tscore, col="blue")

Scale_RMSE<-sqrt(mean(Model_Observed_Error^2, na.rm = TRUE))
Scale_RMSE;