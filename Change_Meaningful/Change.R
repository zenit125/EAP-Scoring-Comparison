#load in the RSSS and ThetaSEeap scoring scripts
source('scripts/RSSS.R')
source("scripts/ThetaSEeap.R")
source("scripts/CAT_sim_missing.R")
source("scripts/Item_Most_Probable_Response.R")
source("scripts/Generate_Random_Data_from_Response_Probabilities.R")
source("scripts/Read_CUE_from_HM_Box.R")
source("scripts/Create_Random_Missing_Data.R")
source("scripts/Create_Random_Missing_Data_With_Mean_Replacement.R")
source("scripts/RMSE.R")
source("scripts/Bland_Altman_Plot.R")

#read in object full of item parameters
Item_Parameters<-readRDS("Item_Parameters.rds")

###########################################################
###Short Forms w/items & item parameters###################
###########################################################

#need to add sleep

Short_Form_Variables<-list(
  "Physical_Function"=list(
     "4a"=c("PFA11", "PFA21", "PFA23", "PFA53"),
     "6b"=c("PFA11", "PFA21", "PFA23", "PFA53", "PFC12", "PFB1"),
     "8b"=c("PFA11", "PFA21", "PFA23", "PFA53", "PFC12", "PFB1", "PFA5", "PFA4"),
     "10a"=c("PFA1", "PFC36r1", "PFC37", "PFA5", "PFA3", "PFA11", "PFA16r1", "PFB26", "PFA55", "PFC45r1"),
     "10b"=c("PFA11", "PFA56", "PFA21", "PFA53", "PFA9", "PFB28r1", "PFA1", "PFA6", "PFB3", "PFB44"),
     "20a"=c("PFA11", "PFA12", "PFA16r1", "PFA34", "PFA38", "PFA51", "PFA55", "PFA56", "PFB19r1", "PFB22", "PFB24", "PFB26", "PFC45r1", "PFC46", "PFA1", "PFA3", "PFA5", "PFC12", "PFC36r1", "PFC37")),
  "Pain_Interference"=list(
    "4a"=c("PAININ9", "PAININ22", "PAININ31", "PAININ34"),
    "6a"=c("PAININ9", "PAININ22", "PAININ31", "PAININ34", "PAININ12", "PAININ36"),
    "6b"=c("PAININ3", "PAININ8", "PAININ9", "PAININ10", "PAININ14", "PAININ26"),
    "8a"=c("PAININ9", "PAININ22", "PAININ31", "PAININ34", "PAININ12", "PAININ36", "PAININ3", "PAININ13")),
  "Fatigue"=list(
    "4a"=c("HI7", "AN3", "FATEXP41", "FATEXP40"),
    "6a"=c("HI7", "AN3", "FATEXP41", "FATEXP40", "FATEXP35", "FATIMP49"),
    "7a"=c("FATEXP20", "FATEXP5", "FATEXP18", "FATIMP33", "FATIMP30", "FATIMP21", "FATIMP40"),
    "8a"=c("HI7", "AN3", "FATEXP41", "FATEXP40", "FATEXP35", "FATIMP49", "FATIMP3", "FATIMP16"),
    "13a"=c("HI7", "HI12", "AN1", "AN2", "AN3", "AN4", "AN5", "AN7", "AN8", "AN12", "AN14", "AN15", "AN16")),
  "Ability"=list(
    "4a"=c("SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS"),
    "6a"=c("SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS", "SRPPER15_CaPS", "SRPPER28r1"),
    "8a"=c("SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS", "SRPPER15_CaPS", "SRPPER28r1", "SRPPER14r1", "SRPPER26_CaPS")),
  "Anxiety"=list(
    "4a"=c("EDANX01", "EDANX40", "EDANX41", "EDANX53"),
    "6a"=c("EDANX01", "EDANX40", "EDANX41", "EDANX53", "EDANX46", "EDANX07"),
    "7a"=c("EDANX01", "EDANX05", "EDANX30", "EDANX40", "EDANX46", "EDANX53", "EDANX54"),
    "8a"=c("EDANX01", "EDANX40", "EDANX41", "EDANX53", "EDANX46", "EDANX07", "EDANX05", "EDANX54")),
  "Depression"=list(
    "4a"=c("EDDEP04","EDDEP06","EDDEP29","EDDEP41"),
    "6a"=c("EDDEP04","EDDEP06","EDDEP29","EDDEP41","EDDEP22","EDDEP36"),
    "8a"=c("EDDEP04","EDDEP06","EDDEP29","EDDEP41","EDDEP22","EDDEP36","EDDEP05","EDDEP09"),
    "8b"=c("EDDEP04", "EDDEP05", "EDDEP06", "EDDEP17", "EDDEP22", "EDDEP29", "EDDEP36", "EDDEP41")))

############################################################
###Big Ole Data File #######################################
############################################################

#length of sim data:
sim_data_length<-1e3

#set seed
set.seed(125)

###Tscore centers (i.e., means) of data
Tscore_Center=c(
  "50 mean"=50,
  "60 mean"=60,
  "40 mean"=40)

Tscore_Change=c(
  "+8"=8,
  "+5"=5,
  "+2"=2,
  "baseline"=0,#Required for center (mean=mean+0)
  "-2"=-2,
  "-5"=-5,
  "-8"=-8)

#Create a list of Item Most Probably Responses- run this only once, then reference later in the sapply statements
Item_Most_Probable_Response_List<-sapply(names(Short_Form_Variables), simplify=FALSE, USE.NAMES=TRUE, function(z){
Item_Most_Probable_Response(Item_Parameters[[z]][,2:ncol(Item_Parameters[[z]])])})

#Simulated Data & Scores
start_time<-Sys.time()
Data<-sapply(names(Short_Form_Variables), simplify=FALSE, USE.NAMES=TRUE, function(domain){
  domain<-sapply(Tscore_Center, simplify=FALSE, USE.NAMES=TRUE, function(Tcenter){
    Tcenter<-sapply(Tscore_Change, simplify = FALSE, USE.NAMES = TRUE, function(Tchange){
        
      Random_Data=Generate_Random_Data_from_Response_Probabilities(sim_data_length, Item_Most_Probable_Response_List[[domain]], mean=Tcenter+Tchange, sd=10)
      print(paste("RD:",domain,Tcenter,Tchange,dim(Random_Data)[1]))
      for(SF in names(Short_Form_Variables[[domain]])){
        RSSS_scores=RSSS_scoring_wrapper(Item_Parameters[[domain]], Short_Form_Variables[[domain]][[SF]], Random_Data)
          names(RSSS_scores)=paste0(SF,"_RSSS_",names(RSSS_scores))
        PatRes_scores=ThetaSEeap_wrapper(Item_Parameters[[domain]], Short_Form_Variables[[domain]][[SF]], Random_Data)
          names(PatRes_scores)=paste0(SF,"_PatRes_",names(PatRes_scores))
        Random_Data=cbind(Random_Data, RSSS_scores, PatRes_scores)}    
      print(paste("SF:",domain,Tcenter,Tchange,dim(Random_Data)[1]))
      CAT_scores=CAT_sim_missing(Item_Parameters[[domain]],Random_Data[,rownames(Item_Parameters[[domain]])])
       names(CAT_scores)=paste0("CAT_",names(CAT_scores))

      print(paste("CAT:",domain,Tcenter,Tchange,dim(CAT_scores)[1]))
      Random_Data=cbind(Random_Data, CAT_scores)      
           Random_Data
    })
  })
})
end_time<-Sys.time()
print((end_time-start_time))

#Data<-readRDS("Data.rds")

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

