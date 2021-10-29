

table_building <- function(table_structure, RMSE_analysis, Short_Form_Variables, Data, trim_na_col=FALSE ){
  SF_lengths=sapply(Short_Form_Variables, length)
  max_ni=max(c(max(SF_lengths),max(Data$CAT_NI)))
  
  if(table_structure=="Single_Timepoint_Scoring_Error"){
    output=matrix(NA,nrow=6,dimnames = list(c("NI","SF","LU-SF","PR-SF","CAT","n (CAT)"),c()))
    for(ni in 1:max_ni){  
      SF_number=sum(SF_lengths %in% ni)
      output2=matrix(NA,nrow=6,ncol=max(1,SF_number), dimnames=dimnames(output))
      output2["NI",]=rep(ni,max(1,SF_number))
      if(SF_number!=0){
        SFs=which(SF_lengths==ni)
        SFs_names=names(Short_Form_Variables)[SFs]
        output2["SF",]=SFs_names
        output2["LU-SF",]=sapply(SFs_names, function(SF_name){
          RMSE(Data$Tscore, Data[,paste0(SF_name,"_RSSS_Tscore")],2)[[RMSE_analysis]]})
        output2["PR-SF",]=sapply(SFs_names, function(SF_name){
          RMSE(Data$Tscore, Data[,paste0(SF_name,"_PatRes_Tscore")],2)[[RMSE_analysis]]})
        }
      if(ni %in% Data$CAT_NI){
        output2["CAT",1]=RMSE(Data$Tscore[which(Data$CAT_NI==ni)], Data$CAT_Tscore[which(Data$CAT_NI==ni)],2)[[RMSE_analysis]]
        output2["n (CAT)",1]=length(which(Data$CAT_NI==ni))
        }
      output=cbind(output,output2)
    }
    if(trim_na_col){output=output[,which(!is.na(output["CAT",])|!is.na(output["SF",]))]}
    return(output)
}}