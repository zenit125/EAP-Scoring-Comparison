#load in the RSSS and ThetaSEeap scoring scripts
source('scripts/RSSS.R')
source("scripts/ThetaSEeap.R")
source("scripts/Item_Most_Probable_Response.R")
source("scripts/Generate_Random_Data_from_Response_Probabilities.R")
source("scripts/Read_CUE_from_HM_Box.R")
source("scripts/Create_Random_Missing_Data.R")
source("scripts/Create_Random_Missing_Data_With_Mean_Replacement.R")
source("scripts/RMSE.R")
source("scripts/Bland_Altman_Plot.R")

#set working directory to read in basic scripts
setwd('/cloud/project')

#read in object full of item parameters
Item_Parameters<-readRDS("/cloud/project/Item_Parameters.rds")


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


###Level 1, type of data, Sim/Real and Dataset type, centering or sample

Scored_data<-list(
  "Simulated Data"=list(
    "0 Theta Centered"=list(),
    "1 Theta Centered"=list(),
    "2 Theta Centered"=list()),
  "Real Data"=list(
    "MyHealth"=list()))

###Level 2, Profiles and datasets NEEDS TO BE TAILORED TO REAL DATA, WHICH DOESN'T HAVE ALL PROFILES

Scored_data<-lapply(Scored_data, function(x){
  x<-lapply(x, function(y)  {
    y<-sapply(names(Short_Form_Variables), simplify=FALSE, USE.NAMES=TRUE, function(z)
    {Random_Data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Item_Most_Probable_Response(Item_Parameters[[z]][,3:ncol(Item_Parameters[[z]])]), mean=50, sd=10)
      z<-lapply(Short_Form_Variables[[z]], function(a){
      vars=c("missing_item_1","missing_item_2","missing_item_3",
             "RSSS_0","RSSS_1","RSSS_2","RSSS_3",
             "ThetaSEeap_0","ThetaSEeap_1","ThetaSEeap","ThetaSEeap_3")
      data.frame(Random_Data[,a],"True_Tscore"=Random_Data$Tscore, matrix(ncol=length(vars),dimnames=list(NULL,vars)))
      })})})})

#redo random data generation based on 0-2 SD out (may have to look at name & change lapply -> sapply)
#make 1-3 random missing items from "a" (short form variables)
#do pattern response and RSSS scoring four times, across levels of missing data
#do RMSE's and Mean Bias's against True Tscore



