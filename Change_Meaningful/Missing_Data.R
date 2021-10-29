#set working directory to read in basic scripts
setwd('/cloud/project')

#read in object full of item parameters
Item_Parameters<-readRDS("/cloud/project/Item_Parameters.rds")

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

Short_Form_Variables<-list()

Scored_data<-list("Sim"=list(),"Real"=list(), "Missing"=list())

#length of sim data:
sim_data_length<-1e3

#set seed
set.seed(125)


########################################################################################################################
###Depression###########################################################################################################
########################################################################################################################
Depression_IPAR<- Item_Parameters$Depression
Depression_Response_Probabilities<-Item_Most_Probable_Response(Depression_IPAR[,3:ncol(Depression_IPAR)])
Depression_sim_data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Depression_Response_Probabilities)
Depression_sim_Tscore<-Depression_sim_data$Tscore
Depression_sim_data<-Depression_sim_data[,which(names(Depression_sim_data)!="Tscore")]

for(n_missing in 0:3){
  
  #4a-EDDEP04, EDDEP06, EDDEP29, EDDEP41
  Depression_4a_variables<-c("EDDEP04","EDDEP06","EDDEP29","EDDEP41")
  ##Raw_Score_Conversion_RSSS
  Depression_4a_RSSS_scoring<-RSSS_scoring_wrapper(Depression_IPAR, Depression_4a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Depression_sim_data, Variables=Depression_4a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Depression_4a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Depression_IPAR, Depression_4a_variables, Create_Random_Missing_Data(Data=Depression_sim_data, Variables=Depression_4a_variables, number_of_missing = n_missing));
  
  #6a-EDDEP04, EDDEP06, EDDEP29, EDDEP41, EDDEP22, EDDEP36
  Depression_6a_variables<-c("EDDEP04","EDDEP06","EDDEP29","EDDEP41","EDDEP22","EDDEP36")
  ##Raw_Score_Conversion_RSSS
  Depression_6a_RSSS_scoring<-RSSS_scoring_wrapper(Depression_IPAR, Depression_6a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Depression_sim_data, Variables=Depression_6a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Depression_6a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Depression_IPAR, Depression_6a_variables, Create_Random_Missing_Data(Data=Depression_sim_data, Variables=Depression_6a_variables, number_of_missing = n_missing));
  #8a-EDDEP04, EDDEP06, EDDEP29, EDDEP51, EDDEP22, EDDEP36, EDDEP05, EDDEP09
  
  Depression_8a_variables<-c("EDDEP04","EDDEP06","EDDEP29","EDDEP41","EDDEP22","EDDEP36","EDDEP05","EDDEP09")
  ##Raw_Score_Conversion_RSSS
  Depression_8a_RSSS_scoring<-RSSS_scoring_wrapper(Depression_IPAR, Depression_8a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Depression_sim_data, Variables=Depression_8a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Depression_8a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Depression_IPAR, Depression_8a_variables, Create_Random_Missing_Data(Data=Depression_sim_data, Variables=Depression_8a_variables, number_of_missing = n_missing));
  
  #8b-EDDEP04, EDDEP05, EDDEP06, EDDEP17, EDDEP22, EDDEP29, EDDEP36, EDDEP41 
  Depression_8b_variables<-c("EDDEP04", "EDDEP05", "EDDEP06", "EDDEP17", "EDDEP22", "EDDEP29", "EDDEP36", "EDDEP41")
  ##Raw_Score_Conversion_RSSS
  Depression_8b_RSSS_scoring<-RSSS_scoring_wrapper(Depression_IPAR, Depression_8b_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Depression_sim_data, Variables=Depression_8b_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Depression_8b_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Depression_IPAR, Depression_8b_variables, Create_Random_Missing_Data(Data=Depression_sim_data, Variables=Depression_8b_variables, number_of_missing = n_missing));
  
  if(n_missing==0){
    #Merge scored sim data for correlations 
    Scored_data$Sim[["Depression"]]=list(
      "True TScore"=Depression_sim_Tscore,
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Depression_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Depression_6a_ThetaSEeap_scoring$TScore, "SF_8a"=Depression_8a_ThetaSEeap_scoring$TScore, "SF_8b"=Depression_8b_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Depression_4a_ThetaSEeap_scoring$SE, "SF_6a"=Depression_6a_ThetaSEeap_scoring$SE, "SF_8a"=Depression_8a_ThetaSEeap_scoring$SE, "SF_8b"=Depression_8b_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8b_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$TScore, "SF_6a"=Depression_6a_RSSS_scoring$TScore, "SF_8a"=Depression_8a_RSSS_scoring$TScore, "SF_8b"=Depression_8b_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$SE, "SF_6a"=Depression_6a_RSSS_scoring$SE, "SF_8a"=Depression_8a_RSSS_scoring$SE, "SF_8b"=Depression_8b_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8b_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Depression_6a_RSSS_scoring$raw_sum_score, "SF_8a"=Depression_8a_RSSS_scoring$raw_sum_score, "SF_8b"=Depression_8b_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$SD, "SF_6a"=Depression_6a_RSSS_scoring$SD, "SF_8a"=Depression_8a_RSSS_scoring$SD, "SF_8b"=Depression_8b_RSSS_scoring$SD)))
    
    Scored_data$Missing[["Depression"]]<-list()}
  
  else{
    #Merge scored sim data for correlations 
    Scored_data$Missing$Depression[[as.character(n_missing)]]=list(
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Depression_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Depression_6a_ThetaSEeap_scoring$TScore, "SF_8a"=Depression_8a_ThetaSEeap_scoring$TScore, "SF_8b"=Depression_8b_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Depression_4a_ThetaSEeap_scoring$SE, "SF_6a"=Depression_6a_ThetaSEeap_scoring$SE, "SF_8a"=Depression_8a_ThetaSEeap_scoring$SE, "SF_8b"=Depression_8b_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8b_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$TScore, "SF_6a"=Depression_6a_RSSS_scoring$TScore, "SF_8a"=Depression_8a_RSSS_scoring$TScore, "SF_8b"=Depression_8b_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$SE, "SF_6a"=Depression_6a_RSSS_scoring$SE, "SF_8a"=Depression_8a_RSSS_scoring$SE, "SF_8b"=Depression_8b_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Depression_sim_Tscore, Depression_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Depression_sim_Tscore, Depression_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Depression_sim_Tscore, Depression_8a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8b"=RMSE(Depression_sim_Tscore, Depression_8b_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Depression_sim_Tscore, Depression_8b_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Depression_6a_RSSS_scoring$raw_sum_score, "SF_8a"=Depression_8a_RSSS_scoring$raw_sum_score, "SF_8b"=Depression_8b_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Depression_4a_RSSS_scoring$SD, "SF_6a"=Depression_6a_RSSS_scoring$SD, "SF_8a"=Depression_8a_RSSS_scoring$SD, "SF_8b"=Depression_8b_RSSS_scoring$SD))
        )}
}


Short_Form_Variables[["Depression"]]=list("SF_4a"=Depression_4a_variables, "SF_6a"=Depression_6a_variables, "SF_8a"=Depression_8a_variables, "SF_8b"=Depression_8b_variables)
rm(list=ls()[which(startsWith(ls(), "Depression")==TRUE)])


########################################################################################################################
###Anxiety##############################################################################################################
########################################################################################################################

Anxiety_IPAR<-Item_Parameters$Anxiety
Anxiety_Response_Probabilities<-Item_Most_Probable_Response(Anxiety_IPAR[,3:ncol(Anxiety_IPAR)])
Anxiety_sim_data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Anxiety_Response_Probabilities)
Anxiety_sim_Tscore<-Anxiety_sim_data$Tscore
Anxiety_sim_data<-Anxiety_sim_data[,which(names(Anxiety_sim_data)!="Tscore")]

for(n_missing in 0:3){
  #4a- EDANX01, EDANX40, EDANX41, EDANX53
  Anxiety_4a_variables<-c("EDANX01", "EDANX40", "EDANX41", "EDANX53")
  ##Raw_Score_Conversion_RSSS
  Anxiety_4a_RSSS_scoring<-RSSS_scoring_wrapper(Anxiety_IPAR, Anxiety_4a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Anxiety_sim_data, Variables=Anxiety_4a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Anxiety_4a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Anxiety_IPAR, Anxiety_4a_variables, Create_Random_Missing_Data(Data=Anxiety_sim_data, Variables=Anxiety_4a_variables, number_of_missing = n_missing));
  
  #6a- EDANX01, EDANX40, EDANX41, EDANX53, EDANX46, EDANX07
  Anxiety_6a_variables<-c("EDANX01", "EDANX40", "EDANX41", "EDANX53", "EDANX46", "EDANX07")
  ##Raw_Score_Conversion_RSSS
  Anxiety_6a_RSSS_scoring<-RSSS_scoring_wrapper(Anxiety_IPAR, Anxiety_6a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Anxiety_sim_data, Variables=Anxiety_6a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Anxiety_6a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Anxiety_IPAR, Anxiety_6a_variables, Create_Random_Missing_Data(Data=Anxiety_sim_data, Variables=Anxiety_6a_variables, number_of_missing = n_missing));
  
  #7a- EDANX01, EDANX05, EDANX30, EDANX40, EDANX46, EDANX53, EDANX54
  Anxiety_7a_variables<-c("EDANX01", "EDANX05", "EDANX30", "EDANX40", "EDANX46", "EDANX53", "EDANX54")
  ##Raw_Score_Conversion_RSSS
  Anxiety_7a_RSSS_scoring<-RSSS_scoring_wrapper(Anxiety_IPAR, Anxiety_7a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Anxiety_sim_data, Variables=Anxiety_7a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Anxiety_7a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Anxiety_IPAR, Anxiety_7a_variables, Create_Random_Missing_Data(Data=Anxiety_sim_data, Variables=Anxiety_7a_variables, number_of_missing = n_missing));
  
  #8a- EDANX01, EDANX40, EDANX41, EDANX53, EDANX46, EDANX07, EDANX05, EDANX54
  Anxiety_8a_variables<-c("EDANX01", "EDANX40", "EDANX41", "EDANX53", "EDANX46", "EDANX07", "EDANX05", "EDANX54")
  ##Raw_Score_Conversion_RSSS
  Anxiety_8a_RSSS_scoring<-RSSS_scoring_wrapper(Anxiety_IPAR, Anxiety_8a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Anxiety_sim_data, Variables=Anxiety_8a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Anxiety_8a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Anxiety_IPAR, Anxiety_8a_variables, Create_Random_Missing_Data(Data=Anxiety_sim_data, Variables=Anxiety_8a_variables, number_of_missing = n_missing));
  
  
  if(n_missing==0){
    #Merge scored sim data for correlations 
    Scored_data$Sim[["Anxiety"]]=list(
      "True TScore"=Anxiety_sim_Tscore,
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Anxiety_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Anxiety_6a_ThetaSEeap_scoring$TScore, "SF_7a"=Anxiety_7a_ThetaSEeap_scoring$TScore, "SF_8a"=Anxiety_8a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Anxiety_4a_ThetaSEeap_scoring$SE, "SF_6a"=Anxiety_6a_ThetaSEeap_scoring$SE, "SF_7a"=Anxiety_7a_ThetaSEeap_scoring$SE, "SF_8a"=Anxiety_8a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_7a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$TScore, "SF_6a"=Anxiety_6a_RSSS_scoring$TScore, "SF_7a"=Anxiety_7a_RSSS_scoring$TScore, "SF_8a"=Anxiety_8a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$SE, "SF_6a"=Anxiety_6a_RSSS_scoring$SE, "SF_7a"=Anxiety_7a_RSSS_scoring$SE, "SF_8a"=Anxiety_8a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_7a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Anxiety_6a_RSSS_scoring$raw_sum_score, "SF_7a"=Anxiety_7a_RSSS_scoring$raw_sum_score, "SF_8a"=Anxiety_8a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$SD, "SF_6a"=Anxiety_6a_RSSS_scoring$SD, "SF_7a"=Anxiety_7a_RSSS_scoring$SD, "SF_8a"=Anxiety_8a_RSSS_scoring$SD)))
    
    Scored_data$Missing[["Anxiety"]]<-list()}
  
  else{
    #Merge scored sim data for correlations 
    Scored_data$Missing$Anxiety[[as.character(n_missing)]]=list(
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Anxiety_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Anxiety_6a_ThetaSEeap_scoring$TScore, "SF_7a"=Anxiety_7a_ThetaSEeap_scoring$TScore, "SF_8a"=Anxiety_8a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Anxiety_4a_ThetaSEeap_scoring$SE, "SF_6a"=Anxiety_6a_ThetaSEeap_scoring$SE, "SF_7a"=Anxiety_7a_ThetaSEeap_scoring$SE, "SF_8a"=Anxiety_8a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_7a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$TScore, "SF_6a"=Anxiety_6a_RSSS_scoring$TScore, "SF_7a"=Anxiety_7a_RSSS_scoring$TScore, "SF_8a"=Anxiety_8a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$SE, "SF_6a"=Anxiety_6a_RSSS_scoring$SE, "SF_7a"=Anxiety_7a_RSSS_scoring$SE, "SF_8a"=Anxiety_8a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Anxiety_sim_Tscore, Anxiety_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Anxiety_sim_Tscore, Anxiety_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Anxiety_sim_Tscore, Anxiety_7a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Anxiety_sim_Tscore, Anxiety_8a_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8a"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_7a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Anxiety_sim_Tscore, Anxiety_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Anxiety_6a_RSSS_scoring$raw_sum_score, "SF_7a"=Anxiety_7a_RSSS_scoring$raw_sum_score, "SF_8a"=Anxiety_8a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Anxiety_4a_RSSS_scoring$SD, "SF_6a"=Anxiety_6a_RSSS_scoring$SD, "SF_7a"=Anxiety_7a_RSSS_scoring$SD, "SF_8a"=Anxiety_8a_RSSS_scoring$SD)))}
}

Short_Form_Variables[["Anxiety"]]=list("SF_4a"=Anxiety_4a_variables , "SF_6a"=Anxiety_6a_variables, "SF_7a"=Anxiety_7a_variables, "SF_8a"=Anxiety_8a_variables)
rm(list=ls()[which(startsWith(ls(), "Anxiety")==TRUE)])


########################################################################################################################
###Ability to Participant in Social Roles###############################################################################
########################################################################################################################

Ability_IPAR<-Item_Parameters$Ability
Ability_Response_Probabilities<-Item_Most_Probable_Response(Ability_IPAR[,2:ncol(Ability_IPAR)])
Ability_sim_data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Ability_Response_Probabilities)
Ability_sim_Tscore<-Ability_sim_data$Tscore
Ability_sim_data<-Ability_sim_data[,which(names(Ability_sim_data)!="Tscore")]


for(n_missing in 0:3){
  #4a- SRPPER11_CaPS, SRPPER18_CaPS, SRPPER23_CaPS, SRPPER46_CaPS
  Ability_4a_variables<-c("SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS")
  ##Raw_Score_Conversion_RSSS
  Ability_4a_RSSS_scoring<-RSSS_scoring_wrapper(Ability_IPAR, Ability_4a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Ability_sim_data, Variables=Ability_4a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Ability_4a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Ability_IPAR, Ability_4a_variables, Create_Random_Missing_Data(Data=Ability_sim_data, Variables=Ability_4a_variables, number_of_missing = n_missing));
  
  #6a- SRPPER11_CaPS, SRPPER18_CaPS, SRPPER23_CaPS, SRPPER46_CaPS, SRPPER15_CaPS, SRPPER28r1 <-
  Ability_6a_variables<-c("SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS", "SRPPER15_CaPS", "SRPPER28r1")
  ##Raw_Score_Conversion_RSSS
  Ability_6a_RSSS_scoring<-RSSS_scoring_wrapper(Ability_IPAR, Ability_6a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Ability_sim_data, Variables=Ability_6a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Ability_6a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Ability_IPAR, Ability_6a_variables, Create_Random_Missing_Data(Data=Ability_sim_data, Variables=Ability_6a_variables, number_of_missing = n_missing));
  
  #8a- SRPPER11_CaPS, SRPPER18_CaPS, SRPPER23_CaPS, SRPPER46_CaPS, SRPPER15_CaPS, SRPPER28r1, SRPPER14r1, SRPPER26_CaPS 
  Ability_8a_variables<-c("SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS", "SRPPER15_CaPS", "SRPPER28r1", "SRPPER14r1", "SRPPER26_CaPS")
  ##Raw_Score_Conversion_RSSS
  Ability_8a_RSSS_scoring<-RSSS_scoring_wrapper(Ability_IPAR, Ability_8a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Ability_sim_data, Variables=Ability_8a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Ability_8a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Ability_IPAR, Ability_8a_variables, Create_Random_Missing_Data(Data=Ability_sim_data, Variables=Ability_8a_variables, number_of_missing = n_missing));
  
  if(n_missing==0){
    #Merge scored sim data for correlations 
    Scored_data$Sim[["Ability"]]=list(
      "True TScore"=Ability_sim_Tscore,
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Ability_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Ability_6a_ThetaSEeap_scoring$TScore, "SF_8a"=Ability_8a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Ability_4a_ThetaSEeap_scoring$SE, "SF_6a"=Ability_6a_ThetaSEeap_scoring$SE,  "SF_8a"=Ability_8a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
      "Bland_Altman"=list(
        "SF_4a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
        "SF_6a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
        "SF_8b"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$TScore, "SF_6a"=Ability_6a_RSSS_scoring$TScore, "SF_8a"=Ability_8a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$SE, "SF_6a"=Ability_6a_RSSS_scoring$SE, "SF_8a"=Ability_8a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Ability_6a_RSSS_scoring$raw_sum_score, "SF_8a"=Ability_8a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$SD, "SF_6a"=Ability_6a_RSSS_scoring$SD, "SF_8a"=Ability_8a_RSSS_scoring$SD)))
    
    Scored_data$Missing[["Ability"]]<-list()}
  
  else{
    #Merge scored sim data for correlations 
    Scored_data$Missing$Ability[[as.character(n_missing)]]=list(
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Ability_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Ability_6a_ThetaSEeap_scoring$TScore, "SF_8a"=Ability_8a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Ability_4a_ThetaSEeap_scoring$SE, "SF_6a"=Ability_6a_ThetaSEeap_scoring$SE,  "SF_8a"=Ability_8a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$TScore, "SF_6a"=Ability_6a_RSSS_scoring$TScore, "SF_8a"=Ability_8a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$SE, "SF_6a"=Ability_6a_RSSS_scoring$SE, "SF_8a"=Ability_8a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Ability_sim_Tscore, Ability_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Ability_sim_Tscore, Ability_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Ability_sim_Tscore, Ability_8a_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Ability_sim_Tscore, Ability_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Ability_6a_RSSS_scoring$raw_sum_score, "SF_8a"=Ability_8a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Ability_4a_RSSS_scoring$SD, "SF_6a"=Ability_6a_RSSS_scoring$SD, "SF_8a"=Ability_8a_RSSS_scoring$SD)))
  }
}

Short_Form_Variables[["Ability"]]=list("SF_4a"=Ability_4a_variables , "SF_6a"=Ability_6a_variables, "SF_8a"=Ability_8a_variables)
rm(list=ls()[which(startsWith(ls(), "Ability")==TRUE)])

########################################################################################################################
###Pain Interference####################################################################################################
########################################################################################################################

Pain_Interference_IPAR<-Item_Parameters$Pain_Interference
Pain_Interference_Response_Probabilities<-Item_Most_Probable_Response(Pain_Interference_IPAR[,3:ncol(Pain_Interference_IPAR)])
Pain_Interference_sim_data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Pain_Interference_Response_Probabilities)
Pain_Interference_sim_Tscore<-Pain_Interference_sim_data$Tscore
Pain_Interference_sim_data<-Pain_Interference_sim_data[,which(names(Pain_Interference_sim_data)!="Tscore")]

for(n_missing in 0:3){
  #4a- PAININ9, PAININ22, PAININ31, PAININ34
  Pain_Interference_4a_variables<-c("PAININ9", "PAININ22", "PAININ31", "PAININ34")
  ##Raw_Score_Conversion_RSSS
  Pain_Interference_4a_RSSS_scoring<-RSSS_scoring_wrapper(Pain_Interference_IPAR, Pain_Interference_4a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Pain_Interference_sim_data, Variables=Pain_Interference_4a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Pain_Interference_4a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Pain_Interference_IPAR, Pain_Interference_4a_variables, Create_Random_Missing_Data(Data=Pain_Interference_sim_data, Variables=Pain_Interference_4a_variables, number_of_missing = n_missing));
  
  #6a- PAININ9, PAININ22, PAININ31, PAININ34, PAININ12, PAININ36
  Pain_Interference_6a_variables<-c("PAININ9", "PAININ22", "PAININ31", "PAININ34", "PAININ12", "PAININ36")
  ##Raw_Score_Conversion_RSSS
  Pain_Interference_6a_RSSS_scoring<-RSSS_scoring_wrapper(Pain_Interference_IPAR, Pain_Interference_6a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Pain_Interference_sim_data, Variables=Pain_Interference_6a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Pain_Interference_6a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Pain_Interference_IPAR, Pain_Interference_6a_variables, Create_Random_Missing_Data(Data=Pain_Interference_sim_data, Variables=Pain_Interference_6a_variables, number_of_missing = n_missing));
  
  #6b- PAININ3, PAININ8, PAININ9, PAININ10, PAININ14, PAININ26
  Pain_Interference_6b_variables<-c("PAININ3", "PAININ8", "PAININ9", "PAININ10", "PAININ14", "PAININ26")
  ##Raw_Score_Conversion_RSSS
  Pain_Interference_6b_RSSS_scoring<-RSSS_scoring_wrapper(Pain_Interference_IPAR, Pain_Interference_6b_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Pain_Interference_sim_data, Variables=Pain_Interference_6b_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Pain_Interference_6b_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Pain_Interference_IPAR, Pain_Interference_6b_variables, Create_Random_Missing_Data(Data=Pain_Interference_sim_data, Variables=Pain_Interference_6b_variables, number_of_missing = n_missing));
  
  #8a- PAININ9, PAININ22, PAININ31, PAININ34, PAININ12, PAININ36, PAININ3, PAININ13
  Pain_Interference_8a_variables<-c("PAININ9", "PAININ22", "PAININ31", "PAININ34", "PAININ12", "PAININ36", "PAININ3", "PAININ13")
  ##Raw_Score_Conversion_RSSS
  Pain_Interference_8a_RSSS_scoring<-RSSS_scoring_wrapper(Pain_Interference_IPAR, Pain_Interference_8a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Pain_Interference_sim_data, Variables=Pain_Interference_8a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Pain_Interference_8a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Pain_Interference_IPAR, Pain_Interference_8a_variables, Create_Random_Missing_Data(Data=Pain_Interference_sim_data, Variables=Pain_Interference_8a_variables, number_of_missing = n_missing));
  
  if(n_missing==0){
    #Merge scored sim data for correlations 
    Scored_data$Sim[["Pain_Interference"]]=list(
      "True TScore"=Pain_Interference_sim_Tscore,
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Pain_Interference_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Pain_Interference_6a_ThetaSEeap_scoring$TScore, "SF_6b"=Pain_Interference_6b_ThetaSEeap_scoring$TScore, "SF_8a"=Pain_Interference_8a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Pain_Interference_4a_ThetaSEeap_scoring$SE, "SF_6a"=Pain_Interference_6a_ThetaSEeap_scoring$SE, "SF_6b"=Pain_Interference_6b_ThetaSEeap_scoring$SE,  "SF_8a"=Pain_Interference_8a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6b_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$TScore, "SF_6a"=Pain_Interference_6a_RSSS_scoring$TScore, "SF_6b"=Pain_Interference_6b_RSSS_scoring$TScore, "SF_8a"=Pain_Interference_8a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$SE, "SF_6a"=Pain_Interference_6a_RSSS_scoring$SE, "SF_6b"=Pain_Interference_6b_RSSS_scoring$SE, "SF_8a"=Pain_Interference_8a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6b_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6b_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Pain_Interference_6a_RSSS_scoring$raw_sum_score, "SF_6b"=Pain_Interference_6b_RSSS_scoring$raw_sum_score, "SF_8a"=Pain_Interference_8a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$SD, "SF_6a"=Pain_Interference_6a_RSSS_scoring$SD, "SF_6b"=Pain_Interference_6b_RSSS_scoring$SD, "SF_8a"=Pain_Interference_8a_RSSS_scoring$SD)))
    
    Scored_data$Missing[["Pain_Interference"]]<-list()}
  
  else{
    #Merge scored sim data for correlations 
    Scored_data$Missing$Pain_Interference[[as.character(n_missing)]]=list(
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Pain_Interference_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Pain_Interference_6a_ThetaSEeap_scoring$TScore, "SF_6b"=Pain_Interference_6b_ThetaSEeap_scoring$TScore, "SF_8a"=Pain_Interference_8a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Pain_Interference_4a_ThetaSEeap_scoring$SE, "SF_6a"=Pain_Interference_6a_ThetaSEeap_scoring$SE, "SF_6b"=Pain_Interference_6b_ThetaSEeap_scoring$SE,  "SF_8a"=Pain_Interference_8a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_4a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6b_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_8a_ThetaSEeap_scoring$TScore, S1_is_true_model = TRUE))),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$TScore, "SF_6a"=Pain_Interference_6a_RSSS_scoring$TScore, "SF_6b"=Pain_Interference_6b_RSSS_scoring$TScore, "SF_8a"=Pain_Interference_8a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$SE, "SF_6a"=Pain_Interference_6a_RSSS_scoring$SE, "SF_6b"=Pain_Interference_6b_RSSS_scoring$SE, "SF_8a"=Pain_Interference_8a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6b_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6b"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Pain_Interference_sim_Tscore, Pain_Interference_8a_RSSS_scoring$TScore)$`RMSE by Tscore`), 
        "Bland_Altman"=list(
          "SF_4a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_4a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6a"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6a_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_6b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_6b_RSSS_scoring$TScore, S1_is_true_model = TRUE),
          "SF_8b"=Bland_Altman_Raw(Pain_Interference_sim_Tscore, Pain_Interference_8a_RSSS_scoring$TScore, S1_is_true_model = TRUE))),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Pain_Interference_6a_RSSS_scoring$raw_sum_score, "SF_6b"=Pain_Interference_6b_RSSS_scoring$raw_sum_score, "SF_8a"=Pain_Interference_8a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Pain_Interference_4a_RSSS_scoring$SD, "SF_6a"=Pain_Interference_6a_RSSS_scoring$SD, "SF_6b"=Pain_Interference_6b_RSSS_scoring$SD, "SF_8a"=Pain_Interference_8a_RSSS_scoring$SD)))}
}

Short_Form_Variables[["Pain_Interference"]]=list("SF_4a"=Pain_Interference_4a_variables , "SF_6a"=Pain_Interference_6a_variables, "SF_6b"=Pain_Interference_6b_variables, "SF_8a"=Pain_Interference_8a_variables)
rm(list=ls()[which(startsWith(ls(), "Pain_Interference")==TRUE)])

########################################################################################################################
###Fatigue##############################################################################################################
########################################################################################################################
#Fatigue Items to Reverse Score
#FATIMP40, FATEXP24, FATEXP31, FATEXP42, FATEXP44, FATEXP54, AN5, AN7

Fatigue_IPAR<-Item_Parameters$Fatigue
Fatigue_Response_Probabilities<-Item_Most_Probable_Response(Fatigue_IPAR[,2:ncol(Fatigue_IPAR)])
Fatigue_sim_data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Fatigue_Response_Probabilities)
Fatigue_sim_Tscore<-Fatigue_sim_data$Tscore
Fatigue_sim_data<-Fatigue_sim_data[,which(names(Fatigue_sim_data)!="Tscore")]

for(n_missing in 0:3){
  #4a- HI7, AN3, FATEXP41, FATEXP40
  Fatigue_4a_variables<-c("HI7", "AN3", "FATEXP41", "FATEXP40")
  ##Raw_Score_Conversion_RSSS
  Fatigue_4a_RSSS_scoring<-RSSS_scoring_wrapper(Fatigue_IPAR, Fatigue_4a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Fatigue_sim_data, Variables=Fatigue_4a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Fatigue_4a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Fatigue_IPAR, Fatigue_4a_variables, Create_Random_Missing_Data(Data=Fatigue_sim_data, Variables=Fatigue_4a_variables, number_of_missing = n_missing));
  
  #6a- HI7, AN3, FATEXP41, FATEXP40, FATEXP35, FATIMP49
  Fatigue_6a_variables<-c("HI7", "AN3", "FATEXP41", "FATEXP40", "FATEXP35", "FATIMP49")
  ##Raw_Score_Conversion_RSSS
  Fatigue_6a_RSSS_scoring<-RSSS_scoring_wrapper(Fatigue_IPAR, Fatigue_6a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Fatigue_sim_data, Variables=Fatigue_6a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Fatigue_6a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Fatigue_IPAR, Fatigue_6a_variables, Create_Random_Missing_Data(Data=Fatigue_sim_data, Variables=Fatigue_6a_variables, number_of_missing = n_missing));
  
  #7a- FATEXP20, FATEXP5, FATEXP18, FATIMP33, FATIMP30, FATIMP21, FATIMP40
  Fatigue_7a_variables<-c("FATEXP20", "FATEXP5", "FATEXP18", "FATIMP33", "FATIMP30", "FATIMP21", "FATIMP40")
  ##Raw_Score_Conversion_RSSS
  Fatigue_7a_RSSS_scoring<-RSSS_scoring_wrapper(Fatigue_IPAR, Fatigue_7a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Fatigue_sim_data, Variables=Fatigue_7a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Fatigue_7a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Fatigue_IPAR, Fatigue_7a_variables, Create_Random_Missing_Data(Data=Fatigue_sim_data, Variables=Fatigue_7a_variables, number_of_missing = n_missing));
  
  #8a- HI7, AN3, FATEXP41, FATEXP40, FATEXP35, FATIMP49, FATIMP3, FATIMP16
  Fatigue_8a_variables<-c("HI7", "AN3", "FATEXP41", "FATEXP40", "FATEXP35", "FATIMP49", "FATIMP3", "FATIMP16")
  ##Raw_Score_Conversion_RSSS
  Fatigue_8a_RSSS_scoring<-RSSS_scoring_wrapper(Fatigue_IPAR, Fatigue_8a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Fatigue_sim_data, Variables=Fatigue_8a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Fatigue_8a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Fatigue_IPAR, Fatigue_8a_variables, Create_Random_Missing_Data(Data=Fatigue_sim_data, Variables=Fatigue_8a_variables, number_of_missing = n_missing));
  
  #13a- HI7, HI12, AN1, AN2, AN3, AN4, AN5, AN7, AN8, AN12, AN14, AN15, AN16 
  Fatigue_13a_variables<-c("HI7", "HI12", "AN1", "AN2", "AN3", "AN4", "AN5", "AN7", "AN8", "AN12", "AN14", "AN15", "AN16")
  ##Raw_Score_Conversion_RSSS
  Fatigue_13a_RSSS_scoring<-RSSS_scoring_wrapper(Fatigue_IPAR, Fatigue_13a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Fatigue_sim_data, Variables=Fatigue_13a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Fatigue_13a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Fatigue_IPAR, Fatigue_13a_variables, Create_Random_Missing_Data(Data=Fatigue_sim_data, Variables=Fatigue_13a_variables, number_of_missing = n_missing));
  
  if(n_missing==0){
    #Merge scored sim data for correlations 
    Scored_data$Sim[["Fatigue"]]=list(
      "True TScore"=Fatigue_sim_Tscore,
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Fatigue_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Fatigue_6a_ThetaSEeap_scoring$TScore, "SF_7a"=Fatigue_7a_ThetaSEeap_scoring$TScore, "SF_8a"=Fatigue_8a_ThetaSEeap_scoring$TScore, "SF_13a"=Fatigue_13a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Fatigue_4a_ThetaSEeap_scoring$SE, "SF_6a"=Fatigue_6a_ThetaSEeap_scoring$SE, "SF_7a"=Fatigue_7a_ThetaSEeap_scoring$SE,  "SF_8a"=Fatigue_8a_ThetaSEeap_scoring$SE, "SF_13a"=Fatigue_13a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`)),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$TScore, "SF_6a"=Fatigue_6a_RSSS_scoring$TScore, "SF_7a"=Fatigue_7a_RSSS_scoring$TScore, "SF_8a"=Fatigue_8a_RSSS_scoring$TScore, "SF_13a"=Fatigue_13a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$SE, "SF_6a"=Fatigue_6a_RSSS_scoring$SE, "SF_7a"=Fatigue_7a_RSSS_scoring$SE, "SF_8a"=Fatigue_8a_RSSS_scoring$SE, "SF_13a"=Fatigue_13a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_RSSS_scoring$TScore)$`RMSE by Tscore`)),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Fatigue_6a_RSSS_scoring$raw_sum_score, "SF_7a"=Fatigue_7a_RSSS_scoring$raw_sum_score, "SF_8a"=Fatigue_8a_RSSS_scoring$raw_sum_score, "SF_13a"=Fatigue_13a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$SD, "SF_6a"=Fatigue_6a_RSSS_scoring$SD, "SF_7a"=Fatigue_7a_RSSS_scoring$SD, "SF_8a"=Fatigue_8a_RSSS_scoring$SD, "SF_13a"=Fatigue_13a_RSSS_scoring$SD)))
    
    Scored_data$Missing[["Fatigue"]]<-list()}
  
  
  else{
    #Merge scored sim data for correlations
    Scored_data$Missing$Fatigue[[as.character(n_missing)]]=list(
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Fatigue_4a_ThetaSEeap_scoring$TScore, "SF_6a"=Fatigue_6a_ThetaSEeap_scoring$TScore, "SF_7a"=Fatigue_7a_ThetaSEeap_scoring$TScore, "SF_8a"=Fatigue_8a_ThetaSEeap_scoring$TScore, "SF_13a"=Fatigue_13a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Fatigue_4a_ThetaSEeap_scoring$SE, "SF_6a"=Fatigue_6a_ThetaSEeap_scoring$SE, "SF_7a"=Fatigue_7a_ThetaSEeap_scoring$SE,  "SF_8a"=Fatigue_8a_ThetaSEeap_scoring$SE, "SF_13a"=Fatigue_13a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`)),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$TScore, "SF_6a"=Fatigue_6a_RSSS_scoring$TScore, "SF_7a"=Fatigue_7a_RSSS_scoring$TScore, "SF_8a"=Fatigue_8a_RSSS_scoring$TScore, "SF_13a"=Fatigue_13a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$SE, "SF_6a"=Fatigue_6a_RSSS_scoring$SE, "SF_7a"=Fatigue_7a_RSSS_scoring$SE, "SF_8a"=Fatigue_8a_RSSS_scoring$SE, "SF_13a"=Fatigue_13a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_RSSS_scoring$TScore)$`Scale level RMSE`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Fatigue_sim_Tscore, Fatigue_4a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_6a"=RMSE(Fatigue_sim_Tscore, Fatigue_6a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_7a"=RMSE(Fatigue_sim_Tscore, Fatigue_7a_RSSS_scoring$TScore)$`RMSE by Tscore`, 
          "SF_8a"=RMSE(Fatigue_sim_Tscore, Fatigue_8a_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_13a"=RMSE(Fatigue_sim_Tscore, Fatigue_13a_RSSS_scoring$TScore)$`RMSE by Tscore`)),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$raw_sum_score, "SF_6a"=Fatigue_6a_RSSS_scoring$raw_sum_score, "SF_7a"=Fatigue_7a_RSSS_scoring$raw_sum_score, "SF_8a"=Fatigue_8a_RSSS_scoring$raw_sum_score, "SF_13a"=Fatigue_13a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Fatigue_4a_RSSS_scoring$SD, "SF_6a"=Fatigue_6a_RSSS_scoring$SD, "SF_7a"=Fatigue_7a_RSSS_scoring$SD, "SF_8a"=Fatigue_8a_RSSS_scoring$SD, "SF_13a"=Fatigue_13a_RSSS_scoring$SD)))
  }}

Short_Form_Variables[["Fatigue"]]=list("SF_4a"=Fatigue_4a_variables , "SF_6a"=Fatigue_6a_variables, "SF_7a"=Fatigue_7a_variables, "SF_8a"=Fatigue_8a_variables, "SF_13a"=Fatigue_13a_variables)
rm(list=ls()[which(startsWith(ls(), "Fatigue")==TRUE)])

########################################################################################################################
###Physical Function####################################################################################################
########################################################################################################################

Physical_Function_IPAR<-Item_Parameters$Physical_Function
Physical_Function_Response_Probabilities<-Item_Most_Probable_Response(Physical_Function_IPAR[,2:ncol(Physical_Function_IPAR)])
Physical_Function_sim_data<-Generate_Random_Data_from_Response_Probabilities(sim_data_length, Physical_Function_Response_Probabilities)
Physical_Function_sim_Tscore<-Physical_Function_sim_data$Tscore
Physical_Function_sim_data<-Physical_Function_sim_data[,which(names(Physical_Function_sim_data)!="Tscore")]

for(n_missing in 0:3){
  
  #4a- PFA11, PFA21, PFA23, PFA53
  Physical_Function_4a_variables<-c("PFA11", "PFA21", "PFA23", "PFA53")
  ##Raw_Score_Conversion_RSSS
  Physical_Function_4a_RSSS_scoring<-RSSS_scoring_wrapper(Physical_Function_IPAR, Physical_Function_4a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Physical_Function_sim_data, Variables=Physical_Function_4a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Physical_Function_4a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Physical_Function_IPAR, Physical_Function_4a_variables, Create_Random_Missing_Data(Data=Physical_Function_sim_data, Variables=Physical_Function_4a_variables, number_of_missing = n_missing));
  
  #6b- PFA11, PFA21, PFA23, PFA53, PFC12, PFB1
  Physical_Function_6b_variables<-c("PFA11", "PFA21", "PFA23", "PFA53", "PFC12", "PFB1")
  ##Raw_Score_Conversion_RSSS
  Physical_Function_6b_RSSS_scoring<-RSSS_scoring_wrapper(Physical_Function_IPAR, Physical_Function_6b_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Physical_Function_sim_data, Variables=Physical_Function_6b_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Physical_Function_6b_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Physical_Function_IPAR, Physical_Function_6b_variables, Create_Random_Missing_Data(Data=Physical_Function_sim_data, Variables=Physical_Function_6b_variables, number_of_missing = n_missing));
  
  
  #8b- PFA11, PFA21, PFA23, PFA53, PFC12, PFB1, PFA5, PFA4
  Physical_Function_8b_variables<-c("PFA11", "PFA21", "PFA23", "PFA53", "PFC12", "PFB1", "PFA5", "PFA4")
  ##Raw_Score_Conversion_RSSS
  Physical_Function_8b_RSSS_scoring<-RSSS_scoring_wrapper(Physical_Function_IPAR, Physical_Function_8b_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Physical_Function_sim_data, Variables=Physical_Function_8b_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Physical_Function_8b_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Physical_Function_IPAR, Physical_Function_8b_variables, Create_Random_Missing_Data(Data=Physical_Function_sim_data, Variables=Physical_Function_8b_variables, number_of_missing = n_missing));
  
  #10a- PFA1, PFC36r1, PFC37, PFA5, PFA3, PFA11, PFA16r1, PFB26, PFA55, PFC45r1
  Physical_Function_10a_variables<-c("PFA1", "PFC36r1", "PFC37", "PFA5", "PFA3", "PFA11", "PFA16r1", "PFB26", "PFA55", "PFC45r1")
  ##Raw_Score_Conversion_RSSS
  Physical_Function_10a_RSSS_scoring<-RSSS_scoring_wrapper(Physical_Function_IPAR, Physical_Function_10a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Physical_Function_sim_data, Variables=Physical_Function_10a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Physical_Function_10a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Physical_Function_IPAR, Physical_Function_10a_variables, Create_Random_Missing_Data(Data=Physical_Function_sim_data, Variables=Physical_Function_10a_variables, number_of_missing = n_missing));
  
  #10b- PFA11, PFA56, PFA21, PFA53, PFA9, PFB28r1, PFA1, PFA6, PFB3, PFB44
  Physical_Function_10b_variables<-c("PFA11", "PFA56", "PFA21", "PFA53", "PFA9", "PFB28r1", "PFA1", "PFA6", "PFB3", "PFB44")
  ##Raw_Score_Conversion_RSSS
  Physical_Function_10b_RSSS_scoring<-RSSS_scoring_wrapper(Physical_Function_IPAR, Physical_Function_10b_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Physical_Function_sim_data, Variables=Physical_Function_10b_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Physical_Function_10b_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Physical_Function_IPAR, Physical_Function_10b_variables, Create_Random_Missing_Data(Data=Physical_Function_sim_data, Variables=Physical_Function_10b_variables, number_of_missing = n_missing));
  
  #20a- PFA11, PFA12, PFA16r1, PFA34, PFA38, PFA51, PFA55, PFA56, PFB19r1, PFB22, PFB24, PFB26, PFC45r1, PFC46, PFA1, PFA3, PFA5, PFC12, PFC36r1, PFC37
  Physical_Function_20a_variables<-c("PFA11", "PFA12", "PFA16r1", "PFA34", "PFA38", "PFA51", "PFA55", "PFA56", "PFB19r1", "PFB22", "PFB24", "PFB26", "PFC45r1", "PFC46", "PFA1", "PFA3", "PFA5", "PFC12", "PFC36r1", "PFC37")
  ##Raw_Score_Conversion_RSSS
  Physical_Function_20a_RSSS_scoring<-RSSS_scoring_wrapper(Physical_Function_IPAR, Physical_Function_20a_variables, Create_Random_Missing_Data_With_Mean_Replacement(Data=Physical_Function_sim_data, Variables=Physical_Function_20a_variables, number_of_missing = n_missing))
  ##Pattern_Scoring_ThetaSEeap
  Physical_Function_20a_ThetaSEeap_scoring<-ThetaSEeap_wrapper(Physical_Function_IPAR, Physical_Function_20a_variables, Create_Random_Missing_Data(Data=Physical_Function_sim_data, Variables=Physical_Function_20a_variables, number_of_missing = n_missing));
  
  if(n_missing==0){
    #Merge scored sim data for correlations 
    
    Scored_data$Sim[["Physical_Function"]]=list(
      "True TScore"= Physical_Function_sim_Tscore,
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Physical_Function_4a_ThetaSEeap_scoring$TScore, "SF_6b"=Physical_Function_6b_ThetaSEeap_scoring$TScore, "SF_8b"=Physical_Function_8b_ThetaSEeap_scoring$TScore, "SF_10a"=Physical_Function_10a_ThetaSEeap_scoring$TScore, "SF_10b"=Physical_Function_10b_ThetaSEeap_scoring$TScore, "SF_20a"=Physical_Function_20a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Physical_Function_4a_ThetaSEeap_scoring$SE, "SF_6b"=Physical_Function_6b_ThetaSEeap_scoring$SE, "SF_8b"=Physical_Function_8b_ThetaSEeap_scoring$SE,  "SF_10a"=Physical_Function_10a_ThetaSEeap_scoring$SE, "SF_10b"=Physical_Function_10b_ThetaSEeap_scoring$SE, "SF_20a"=Physical_Function_20a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`)),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$TScore, "SF_6b"=Physical_Function_6b_RSSS_scoring$TScore, "SF_8b"=Physical_Function_8b_RSSS_scoring$TScore, "SF_10a"=Physical_Function_10a_RSSS_scoring$TScore, "SF_10b"=Physical_Function_10b_RSSS_scoring$TScore, "SF_20a"=Physical_Function_20a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$SE, "SF_6b"=Physical_Function_6b_RSSS_scoring$SE, "SF_8b"=Physical_Function_8b_RSSS_scoring$SE, "SF_10a"=Physical_Function_10a_RSSS_scoring$SE, "SF_10b"=Physical_Function_10b_RSSS_scoring$SE, "SF_20a"=Physical_Function_20a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_RSSS_scoring$TScore)$`RMSE by Tscore`)),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$raw_sum_score, "SF_6b"=Physical_Function_6b_RSSS_scoring$raw_sum_score, "SF_8b"=Physical_Function_8b_RSSS_scoring$raw_sum_score, "SF_10a"=Physical_Function_10a_RSSS_scoring$raw_sum_score, "SF_10b"=Physical_Function_10b_RSSS_scoring$raw_sum_score, "SF_20a"=Physical_Function_20a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$SD, "SF_6b"=Physical_Function_6b_RSSS_scoring$SD, "SF_8b"=Physical_Function_8b_RSSS_scoring$SD, "SF_10a"=Physical_Function_10a_RSSS_scoring$SD, "SF_10b"=Physical_Function_10b_RSSS_scoring$SD, "SF_20a"=Physical_Function_20a_RSSS_scoring$SD)))
    
    Scored_data$Missing[["Physical_Function"]]<-list()}
  
  else{
    #Merge scored sim data for correlations
    Scored_data$Missing$Physical_Function[[as.character(n_missing)]]=list(
      "ThetaSEeap_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Physical_Function_4a_ThetaSEeap_scoring$TScore, "SF_6b"=Physical_Function_6b_ThetaSEeap_scoring$TScore, "SF_8b"=Physical_Function_8b_ThetaSEeap_scoring$TScore, "SF_10a"=Physical_Function_10a_ThetaSEeap_scoring$TScore, "SF_10b"=Physical_Function_10b_ThetaSEeap_scoring$TScore, "SF_20a"=Physical_Function_20a_ThetaSEeap_scoring$TScore),
        "SE"=data.frame("SF_4a"=Physical_Function_4a_ThetaSEeap_scoring$SE, "SF_6b"=Physical_Function_6b_ThetaSEeap_scoring$SE, "SF_8b"=Physical_Function_8b_ThetaSEeap_scoring$SE,  "SF_10a"=Physical_Function_10a_ThetaSEeap_scoring$SE, "SF_10b"=Physical_Function_10b_ThetaSEeap_scoring$SE, "SF_20a"=Physical_Function_20a_ThetaSEeap_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_ThetaSEeap_scoring$TScore)$`Scale level RMSE`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_ThetaSEeap_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_ThetaSEeap_scoring$TScore)$`RMSE by Tscore`)),
      "RSSS_scoring_dfs"=list(
        "TScore"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$TScore, "SF_6b"=Physical_Function_6b_RSSS_scoring$TScore, "SF_8b"=Physical_Function_8b_RSSS_scoring$TScore, "SF_10a"=Physical_Function_10a_RSSS_scoring$TScore, "SF_10b"=Physical_Function_10b_RSSS_scoring$TScore, "SF_20a"=Physical_Function_20a_RSSS_scoring$TScore),
        "SE"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$SE, "SF_6b"=Physical_Function_6b_RSSS_scoring$SE, "SF_8b"=Physical_Function_8b_RSSS_scoring$SE, "SF_10a"=Physical_Function_10a_RSSS_scoring$SE, "SF_10b"=Physical_Function_10b_RSSS_scoring$SE, "SF_20a"=Physical_Function_20a_RSSS_scoring$SE),
        "RMSE"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_RSSS_scoring$TScore)$`Scale level RMSE`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_RSSS_scoring$TScore)$`Scale level RMSE`),
        "RMSE by Tscore"=list(
          "SF_4a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_4a_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_6b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_6b_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_8b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_8b_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_10a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10a_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_10b"=RMSE(Physical_Function_sim_Tscore, Physical_Function_10b_RSSS_scoring$TScore)$`RMSE by Tscore`,
          "SF_20a"=RMSE(Physical_Function_sim_Tscore, Physical_Function_20a_RSSS_scoring$TScore)$`RMSE by Tscore`)),
      "Descriptive_stats_dfs"=list(
        "Raw_sum_score"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$raw_sum_score, "SF_6b"=Physical_Function_6b_RSSS_scoring$raw_sum_score, "SF_8b"=Physical_Function_8b_RSSS_scoring$raw_sum_score, "SF_10a"=Physical_Function_10a_RSSS_scoring$raw_sum_score, "SF_10b"=Physical_Function_10b_RSSS_scoring$raw_sum_score, "SF_20a"=Physical_Function_20a_RSSS_scoring$raw_sum_score),
        "SD"=data.frame("SF_4a"=Physical_Function_4a_RSSS_scoring$SD, "SF_6b"=Physical_Function_6b_RSSS_scoring$SD, "SF_8b"=Physical_Function_8b_RSSS_scoring$SD, "SF_10a"=Physical_Function_10a_RSSS_scoring$SD, "SF_10b"=Physical_Function_10b_RSSS_scoring$SD, "SF_20a"=Physical_Function_20a_RSSS_scoring$SD)))    
    
  }}

Short_Form_Variables[["Physical_Function"]]=list("SF_4a"=Physical_Function_4a_variables , "SF_6b"=Physical_Function_6b_variables, "SF_8b"=Physical_Function_8b_variables, "SF_10a"=Physical_Function_10a_variables, "SF_10b"=Physical_Function_10b_variables, "SF_20a"=Physical_Function_20a_variables)
rm(list=ls()[which(startsWith(ls(), "Physical_Function")==TRUE)])

#########################################################

