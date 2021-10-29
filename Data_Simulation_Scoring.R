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
Item_Parameters<-readRDS("~/Scripts/Item_Parameters.rds")

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

saveRDS(Short_Form_Variables, "~/EAP_Scoring_Comparison/Short_Form_Variables.RDS")

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


#Load 'Calculate_Response_Probabilities' script from Github
source("https://raw.githubusercontent.com/zenit125/Psychometric_Scripts/main/Calculate_Response_Probabilities")

#Load 'Generate_Random_Data_from_Response_Probabilities' script from Github
source("https://raw.githubusercontent.com/zenit125/Psychometric_Scripts/main/Generate_Random_Data_from_Response_Probabilities.R")


Response_Probabilities=Calculate_Response_Probabilities(ipar)
Sim_Data=Generate_Random_Data_from_Response_Probabilities(N=1000, Response_Probabilities, mean=50, sd=10)


#Create a list of Item Most Probably Responses- run this only once, then reference later in the sapply statements
Response_Probabilities_List<-sapply(names(Short_Form_Variables), simplify=FALSE, USE.NAMES=TRUE, function(z){
Calculate_Response_Probabilities(Item_Parameters[[z]][,2:ncol(Item_Parameters[[z]])])})

#Load in "RSSS" scoring script from github
source("https://raw.githubusercontent.com/zenit125/Psychometric_Scripts/main/RSSS.R")

#Load in "ThetaSEeap" scoring script from github
source("https://raw.githubusercontent.com/zenit125/Psychometric_Scripts/main/ThetaSEeap.R")

#Load in "CAT_sim_missing" CAT simulation script
source("~/EAP_Scoring_Comparison/CAT_sim_missing.R")

#Simulated Data & Scores
Data<-sapply(names(Short_Form_Variables), simplify=FALSE, USE.NAMES=TRUE, function(domain){
  domain<-sapply(Tscore_Center, simplify=FALSE, USE.NAMES=TRUE, function(Tcenter){
    Tcenter<-sapply(Tscore_Change, simplify = FALSE, USE.NAMES = TRUE, function(Tchange){
        
      Random_Data=Generate_Random_Data_from_Response_Probabilities(sim_data_length, Response_Probabilities_List[[domain]], mean=Tcenter+Tchange, sd=10)
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

saveRDS(Data, "Data.rds")

