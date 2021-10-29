#set working directory to read in basic scripts
setwd("C:\\Users\\Robert Chapman\\Desktop\\scripts\\")
#load in the script to generate item probabilities
source("Item_Most_Probable_Response.r")

#Create a dataframe with the 6 columns for the (1) slope, (4) thresholds, (1)number of response categories
IPAR<-data.frame(matrix(nrow=0,ncol=6))
names(IPAR)<-c("a","cb1","cb2","cb3","cb4","NCAT")

#populate the dataframe with calibration statistics
#first item
IPAR[1,"a"]<-2.5
IPAR[1,"cb1"]<- -2.5
IPAR[1,"cb2"]<- -1.25
IPAR[1,"cb3"]<- 0
IPAR[1,"cb4"]<- 2.25
IPAR[1,"NCAT"]<- 5
rownames(IPAR)[1]<-"Robs_Item_1"

#second item
IPAR[2,"a"]<-1.25
IPAR[2,"cb1"]<- -0.75
IPAR[2,"cb2"]<- -0.5
IPAR[2,"cb3"]<- 0
IPAR[2,"cb4"]<- 0.15
IPAR[2,"NCAT"]<- 5
rownames(IPAR)[2]<-"Robs_Item_2"

#second item
IPAR[3,"a"]<- 0.75
IPAR[3,"cb1"]<- 0.125
IPAR[3,"cb2"]<- 0.375
IPAR[3,"cb3"]<- 0.5
IPAR[3,"cb4"]<- 0.75
IPAR[3,"NCAT"]<- 5
rownames(IPAR)[3]<-"Robs_Item_3"

#run the script to generate item probabilities
Item_Probabilities<-Item_Most_Probable_Response(IPAR)

#Item_Probabilites is a "list" object, so you have to interact with it in a different way:
#plot response 1 in Red
plot(c(10:90),Item_Probabilities$Robs_Item_1$`1`,col="red" ,type="l", lty=1, xlab="T-score", ylab="probabilities")

#add response 2 in green
lines(c(10:90),Item_Probabilities$Robs_Item_1$`2`,col="green" ,type="l", lty=1)
#add response 3 in purple
lines(c(10:90),Item_Probabilities$Robs_Item_1$`3`,col="purple" ,type="l", lty=1)
#add response 4 in orange
lines(c(10:90),Item_Probabilities$Robs_Item_1$`4`,col="orange" ,type="l", lty=1)
#add response 5 in yellow
lines(c(10:90),Item_Probabilities$Robs_Item_1$`5`,col="yellow" ,type="l", lty=1)


