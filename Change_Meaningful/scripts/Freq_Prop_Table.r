Freq_Prop_Table <- function(Data)
{
  Freq_Tbl<-table(Data, useNA="always")
  Prop_Tbl<-round(prop.table(Freq_Tbl),2)*100
  Freq_Prop_Tbl<-Freq_Tbl
  for(i in 1:length(Freq_Prop_Tbl))
  {
    Freq_Prop_Tbl[i]<-paste0(Freq_Tbl[i], "(",Prop_Tbl[i],"%)")
  }
  
  return(Freq_Prop_Tbl)
}