
rsss<-function(ipar,model=1,minTheta=-4.0,maxTheta=4.0,inc=0.01,prior.mean=0.0,prior.sd=1.0,D=1.0,maxCat=5,minScore=1,Tscore=T){
	NCAT<-ipar[,"NCAT"]
	DISC<-ipar[,"a"]
	CB<-ipar[paste("cb",1:(maxCat-1),sep="")]

	ni<-dim(ipar)[1] # number of items

	theta<-seq(minTheta,maxTheta,by=inc) # populate theta vector
	nq<-length(theta) # number of quadrature points
	
	pp<-array(0,c(nq,ni,maxCat))

	if (model==1) {
		for (i in 1:ni) {
			ps<-matrix(0,nq,NCAT[i]+1);
			ps[,1]<-1;
			ps[,NCAT[i]+1]<-0;
			for (k in 1:(NCAT[i]-1)) {
				ps[,k+1]<-1/(1+exp(-D*DISC[i]*(theta-CB[i,k])));
			}
			pp[,i,1]<-1-ps[,1];
			pp[,i,NCAT[i]]<-ps[,NCAT[i]];
			for (k in 1:NCAT[i]) {
				pp[,i,k]=ps[,k]-ps[,k+1];
			}
		}
	} else if (model==2) {
		for (i in 1:ni) {
			cb<-unlist(CB[i,]);
			cb<-c(0,cb);
			zz<-matrix(0,nq,NCAT[i]);
			sdsum<-0;
			den<-rep(0,nq);

			for (k in 1:NCAT[i]) {
				sdsum<-sdsum+cb[k];
				zz[,k]<-exp(D*DISC[i]*(k*theta-sdsum));
				den<-den+zz[,k];
			}
			for (k in 1:NCAT[i]) {
				pp[,i,k]<-zz[,k]/den;
			}
		}
	}

	min.Raw.Score<-0 #minimum obtainable raw score
	max.Raw.Score<-sum(ipar[,"NCAT"])-ni #maximum obtainable raw score

	nScore<-max.Raw.Score-min.Raw.Score+1 #number of score points
	TCCinv<-numeric(nScore) #initialize TCC scoring table
	Raw.Score<-min.Raw.Score:max.Raw.Score #raw scores

	LH<-matrix(0,nq,nScore) #initializing distribution of summed scores

	ncat<-ipar[1,"NCAT"]
	maxScore<-0
	LH[,1:ncat]<-pp[,1,1:ncat]
	idx<-ncat

	for (i in 2:ni) {
		ncat<-ipar[i,"NCAT"] #number of categories for item i
		maxScore<-ncat-1 #maximum score for item i
		score<-0:maxScore #score values for item i
		prob<-pp[,i,1:ncat] #category probabilities for item i

		pLH<-matrix(0,nq,nScore) #place holder for LH

		for (k in 1:ncat) {
			for (h in 1:idx) {
				sco<-Raw.Score[h]+score[k]
				position<-which(Raw.Score==sco)
				pLH[,position]<-pLH[,position]+LH[,h]*prob[,k]
			}
		}
		idx<-idx+maxScore
		LH<-pLH
	}

	Scale.Score<-numeric(nScore) #score table for EAP
	SE<-numeric(nScore) #SE for EAP

	prior<-dnorm((theta-prior.mean)/prior.sd)
	posterior<-LH*prior #posterior distribution
	den<-colSums(posterior)
	den<-matrix(rep(den,rep(nq,nScore)),nq,nScore)
	posterior<-posterior/den

	for (j in 1:nScore) {
		Scale.Score[j]<-sum(posterior[,j]*theta)/sum(posterior[,j]) #EAP
		SE[j]<-sqrt(sum(posterior[,j]*(theta-Scale.Score[j])^2)/sum(posterior[,j])) #EAP
	}

	if (minScore==1) Raw.Score<-Raw.Score+ni

	if (Tscore) {
		Scale.Score=round(Scale.Score*10+50,1)
		SE=round(SE*10,1)
	}


	rsss.table<-data.frame(Raw=Raw.Score,Scale=Scale.Score,SE)

	return(rsss.table)
}


RSSS_scoring_wrapper<-function(IPAR, Variables, Data)
{
  RSSS_table<-rsss(IPAR[Variables,2:ncol(IPAR)])
  
  RSSS_scoring<-data.frame(matrix(NA,ncol=3,nrow=nrow(Data)))
  names(RSSS_scoring)<-c("raw_sum_score","Tscore","SE")
  RSSS_scoring$raw_sum_score<-round(apply(Data[,Variables], 1, sum, na.rm = FALSE))
  RSSS_scoring$SD<-apply(Data[,Variables], 1, sd)
  
  for(i in 1:nrow(RSSS_table))
  {
    matched_scores=which(RSSS_scoring$raw_sum_score==RSSS_table$Raw[i])
    RSSS_scoring$Tscore[matched_scores]<-RSSS_table$Scale[i]
    RSSS_scoring$SE[matched_scores]<-RSSS_table$SE[i]
  }
  
  return(RSSS_scoring)
  
}

