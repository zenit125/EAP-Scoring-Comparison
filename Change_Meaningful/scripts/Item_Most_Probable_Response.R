Item_Most_Probable_Response<-function(IPAR_dataframe, Tscore_range=NA, Theta_range=NA){
  if(is.na(Tscore_range) & is.na(Theta_range)){
  Tscore_range=c(10,90)
  Theta_range=c(-4,4)}
  else{
    if(!is.na(Tscore_range)){Theta_range<-c((min(Tscore_range)-50)/10,(max(Tscore_range)-50)/10)}
    if(!is.na(Theta_range & is.na(Tscore_range))){Tscore_range<-c((min(Theta_range)+50)/10,(max(Theta_range)+50)/10)}
  }
  
  Theta<-seq(min(Theta_range), max(Theta_range), .1)
  Tscore<-seq(min(Tscore_range), max(Tscore_range), 1)
  most_probable_response<-list()
  for(m in 1:nrow(IPAR_dataframe))
  {
    IPAR<-IPAR_dataframe[m,]
    ThetaProbMatrix <- matrix(nrow = length(Theta), ncol = (IPAR$NCAT-1)) 
    #i is the rows, which are theta increments
    for (i in 1:length(Theta)) {
      #j is columns or response the categories
      for (j in 1:(IPAR$NCAT-1)) {
        ThetaProbMatrix[i, j] <- exp(IPAR$a * (Theta[i] - IPAR[,paste0("cb",j)])) / (1 + exp(IPAR$a * (Theta[i] - IPAR[,paste0("cb",j)])))} } 
    
        ThetaProbMatrix <- cbind((1 - ThetaProbMatrix[, 1]), ThetaProbMatrix)
    for(k in 2:(IPAR$NCAT-1)){ThetaProbMatrix[,k]<-ThetaProbMatrix[,k]-ThetaProbMatrix[,k+1]}
    
    ThetaProbMatrix_df<-as.data.frame(ThetaProbMatrix)
    ThetaProbMatrix_df[ThetaProbMatrix_df<0]<-0
    
    names(ThetaProbMatrix_df)<-c(1:IPAR$NCAT);
    rownames(ThetaProbMatrix_df)<-Tscore
    most_probable_response[[rownames(IPAR_dataframe)[m]]]<-ThetaProbMatrix_df
  }
  return(most_probable_response)
}

