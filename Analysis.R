#load in the RSSS and ThetaSEeap scoring scripts
source("scripts/Create_Random_Missing_Data.R")
source("scripts/Create_Random_Missing_Data_With_Mean_Replacement.R")

#Load in the 'RMSE' script from 'Psychometric Scripts' repo on github
source("https://raw.githubusercontent.com/zenit125/Psychometric_Scripts/main/RMSE.R")

#Load in the 'Bland_Altman_Plot' script from 'Psychometric Scripts' repo on github
source("https://raw.githubusercontent.com/zenit125/Psychometric_Scripts/main/Bland_Altman_Plot.R")

#load in the 'table_building' function from the locati
source("~/EAP_Scoring_Comparison/table_building.R")

#Read in previously generated and scored simulated data
Data=readRDS("Data.rds")

Short_Form_Variables=readRDS("~/EAP_Scoring_Comparison/Short_Form_Variables.RDS")

#mapply analyses
#***************************************
#**Single Timepoint (scoring error)*
#***************************************
####One table for each domain (tables, in a list)
#####SF by LuSF/PrSF & NI by CAT 
#######RMSE Across Tscore lvls - RMSE(True Tscore vs Estimated Tscore) 
#######RMSE Range by tscore lvl - RMSE(True Tscore vs Estimated Tscore) 
#######RMSE SD across tscore lvl - RMSE(True Tscore vs Estimated Tscore) 

Single_Timepoint_Scoring_Error<-
  sapply(names(Data), simplify=FALSE, USE.NAMES=TRUE, function(domain){
    sapply(names(Data[[domain]]), simplify=FALSE, USE.NAMES=TRUE, function(Tcenter){
      Baseline_Data=Data[[domain]][[Tcenter]][["baseline"]]
       list("RMSE Summary"= table_building("Single_Timepoint_Scoring_Error","Summary",Short_Form_Variables[[domain]], Baseline_Data),
            "RMSE Across All Tscores"= table_building("Single_Timepoint_Scoring_Error","Scale level RMSE",Short_Form_Variables[[domain]], Baseline_Data),
            "RMSE Range Across Tscores"= table_building("Single_Timepoint_Scoring_Error","RMSE by Tscore range", Short_Form_Variables[[domain]], Baseline_Data),
            "RMSE SD Across Tscores"= table_building("Single_Timepoint_Scoring_Error","RMSE by Tscore SD",Short_Form_Variables[[domain]] , Baseline_Data))      
    })    
  })


output=sapply(names(Single_Timepoint_Scoring_Error), simplify=FALSE, USE.NAMES=TRUE, function(domain){
    output=matrix(NA,nrow=0, ncol=ncol(Single_Timepoint_Scoring_Error[[domain]][[1]][["RMSE Summary"]]))
    for(Tcenter in names(Single_Timepoint_Scoring_Error[[domain]])){
      output=rbind(output, rep(Tcenter,ncol(output)))
      output=rbind(output,Single_Timepoint_Scoring_Error[[domain]][[Tcenter]][["RMSE Summary"]])}
    output
  })

max_ncol=max(sapply(output, function(x){dim(x)[2]}))
output2=matrix(NA, nrow=0, ncol=max_ncol)
for(domain in names(output)){
  output2=rbind(output2, rep(domain, max_ncol))
  output2=rbind(output2, cbind(output[[domain]],matrix(NA, nrow=nrow(output[[domain]]), ncol=(max_ncol-ncol(output[[domain]])))))} 

write.csv(output2, paste0(getwd(),"/output/Single_Timepoint_Scoring_Error.csv"))


#Group_Level_Baseline_Means
Group_Level_Baseline_Means<-sapply(names(Data), function(domain){
  SFs=names(Short_Form_Variables[[domain]])
  SFs_length=length(SFs)
  
  results_output=matrix(NA, nrow=(SFs_length*2+3), ncol=2, dimnames = list(c(),c("SF/CAT","LU/PR")))
  results_output[,1]=c(NA,NA,sapply(SFs,function(x){rep(x,2)}),"CAT")
  results_output[,2]=c(NA,NA,sapply(SFs,function(x){c("LU","PR")}),NA)
  
  for(Tcenter in names(Data[[domain]])){
  
      results_output2=matrix(NA,nrow=dim(results_output)[1],ncol=1)
      results_output2[1,]=Tcenter
      results_output2[2,]="baseline"
      results_output2[3:dim(results_output)[1],]=
        sapply(3:dim(results_output)[1], simplify = "array", function(Scoring_index){
          Score_var=NA
          if(results_output[Scoring_index,"SF/CAT"]!="CAT")
          {Score_var=paste0(results_output[Scoring_index,"SF/CAT"],ifelse(results_output[Scoring_index,"LU/PR"]=="LU","_RSSS_Tscore","_PatRes_Tscore"))}
          else{Score_var="CAT_Tscore"}
          mean=round(mean(Data[[domain]][[Tcenter]][["baseline"]][[Score_var]]),2)
          SD=round(sd(Data[[domain]][[Tcenter]][["baseline"]][[Score_var]]),2)
          paste0(mean, "(", SD, ")")
          
          
        })
      results_output<-cbind(results_output,results_output2)
    
  }
  results_output
})

output=matrix(NA,nrow=0, ncol=ncol(Group_Level_Baseline_Means[[1]]))
for(domain in names(Group_Level_Baseline_Means)){
  output=rbind(output, rep(domain, ncol(output)))
  output=rbind(output,Group_Level_Baseline_Means[[domain]])
}
write.csv(output,paste0(getwd(),"/output/Group_Level_Baseline_Means.csv"))





###Group-level Change
####One table for each domain
#####SF/CAT (&LU/SF) by center and change amount
#this is all based on full sample size - can resample and add ranges/RMSE's

Group_Level_Change<-sapply(names(Data), function(domain){
  SFs=names(Short_Form_Variables[[domain]])
  SFs_length=length(SFs)
  
  results_output=matrix(NA, nrow=(SFs_length*2+3), ncol=2, dimnames = list(c(),c("SF/CAT","LU/PR")))
  results_output[,1]=c(NA,NA,sapply(SFs,function(x){rep(x,2)}),"CAT")
  results_output[,2]=c(NA,NA,sapply(SFs,function(x){c("LU","PR")}),NA)
  
  for(Tcenter in names(Data[[domain]])){
    Tchange_number=length(which(names(Data[[domain]][[Tcenter]])!="baseline"))
    for(Tchange in names(Data[[domain]][[Tcenter]])[which(names(Data[[domain]][[Tcenter]])!="baseline")]){
      results_output2=matrix(NA,nrow=dim(results_output)[1],ncol=1)
      results_output2[1,]=Tcenter
      results_output2[2,]=Tchange
      results_output2[3:dim(results_output)[1],]=
      sapply(3:dim(results_output)[1], simplify = "array", function(Scoring_index){
        Score_var=NA
          if(results_output[Scoring_index,"SF/CAT"]!="CAT")
          {Score_var=paste0(results_output[Scoring_index,"SF/CAT"],ifelse(results_output[Scoring_index,"LU/PR"]=="LU","_RSSS_Tscore","_PatRes_Tscore"))}
          else{Score_var="CAT_Tscore"}
          mean_diff=round(mean(Data[[domain]][[Tcenter]][[Tchange]][[Score_var]])-mean(Data[[domain]][[Tcenter]][["baseline"]][[Score_var]]),2)
          SD_diff=round(sd(Data[[domain]][[Tcenter]][[Tchange]][[Score_var]]-sd(Data[[domain]][[Tcenter]][["baseline"]][[Score_var]])),2)
          paste0(mean_diff, ", d=", round(mean_diff/SD_diff,2))
          
          
        })
      results_output<-cbind(results_output,results_output2)
    }
  }
  results_output
})

output=matrix(NA,nrow=0, ncol=ncol(Group_Level_Change[[1]]))
for(domain in names(Group_Level_Change)){
  output=rbind(output, rep(domain, ncol(output)))
  output=rbind(output,Group_Level_Change[[domain]])
}
write.csv(output,paste0(getwd(),"/output/Group_Level_Change.csv"))

###Individual-level Change


#move the domain level up, centering down
#make sure that everything is reported on tscore metric
#insert sleep measures & cross-check calibrations for sleep&social


Individual_Level_Change<-sapply(names(Data), function(domain){
  SFs=names(Short_Form_Variables[[domain]])
  SFs_length=length(SFs)
  
  results_output=matrix(NA, nrow=(SFs_length*2+3), ncol=2, dimnames = list(c(),c("SF/CAT","LU/PR")))
  results_output[,1]=c(NA,NA,sapply(SFs,function(x){rep(x,2)}),"CAT")
  results_output[,2]=c(NA,NA,sapply(SFs,function(x){c("LU","PR")}),NA)
  
  for(Tcenter in names(Data[[domain]])){
    Tchange_number=length(which(names(Data[[domain]][[Tcenter]])!="baseline"))
    for(Tchange in names(Data[[domain]][[Tcenter]])[which(names(Data[[domain]][[Tcenter]])!="baseline")]){
      results_output2=matrix(NA,nrow=dim(results_output)[1],ncol=1)
      results_output2[1,]=Tcenter
      results_output2[2,]=Tchange
      results_output2[3:dim(results_output)[1],]=
        sapply(3:dim(results_output)[1], simplify = "array", function(Scoring_index){
          Score_var=NA
          if(results_output[Scoring_index,"SF/CAT"]!="CAT")
          {Score_var=paste0(results_output[Scoring_index,"SF/CAT"],ifelse(results_output[Scoring_index,"LU/PR"]=="LU","_RSSS_Tscore","_PatRes_Tscore"))}
          else{Score_var="CAT_Tscore"}
          #mean(Data[[domain]][[Tcenter]][[Tchange]]$Tscore-Data[[domain]][[Tcenter]][["baseline"]]$Tscore)
          #Error=0
          #if(abs(as.numeric(Tchange))==2){Error=1}
          #if(abs(as.numeric(Tchange))==5){Error=2}
          #if(abs(as.numeric(Tchange))==8){Error=3}
          
          #Individual_Change=abs(c(Data[[domain]][[Tcenter]][[Tchange]]$Tscore-Data[[domain]][[Tcenter]][["baseline"]]$Tscore)
          Individual_Change=3
          Theshold_percent=round(100*sum((Individual_Change) <= abs(c(Data[[domain]][[Tcenter]][[Tchange]][[Score_var]]-Data[[domain]][[Tcenter]][["baseline"]][[Score_var]])))/nrow(Data[[domain]][[Tcenter]][["baseline"]]),2)
        })
      results_output<-cbind(results_output,results_output2)
    }
  }
  results_output
})

output=matrix(NA,nrow=0, ncol=ncol(Individual_Level_Change[[1]]))
for(domain in names(Individual_Level_Change)){
  output=rbind(output, rep(domain, ncol(output)))
  output=rbind(output,Individual_Level_Change[[domain]])
}
write.csv(output,paste0(getwd(),"/output/Individual_Level_Change.csv"))

