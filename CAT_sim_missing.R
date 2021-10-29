
CAT_sim_missing<-function(ipar, responses, missing=F)
{
  ################################################################################
  # 5. simulate a CAT based on the item parameters and the full-bank responses
  ################################################################################ 
  
  maxCat <- 5 # maximum number of categories
  minTheta <- -4.0    #lower bound of theta -- USED IN MAX CATEGORY DETERMINATION AND IN CAT SIMULATION
  maxTheta <- 4.0     #larger bound of theta -- USED IN MAX CATEGORY DETERMINATION AND IN CAT SIMULATION
  RespInc <- 0.25    #spacing between two bounds
  D <- 1
  
  
  CATmethod <- 3 # 3=default (MPWI); 1=MFI, 2=MLWI, 3=MPWI, 4=MEI, 5=MEPV, 6=MEPWI 7:11=do not use herein
  maxNI <- 12 # maximum number of items
  minNI <- 12 # minimum number of items
  maxSE <- .3 # standard error stopping rule
  EAP_Mean <- 0 # 0=default; mean for EAP scoring
  EAP_SD <- 1 # 1=default; sd for EAP scoring
  ScoreInc <- 0.1 # 0.1=default, uses minTheta and maxTheta from above
  StartItem <- 0 # 0=default, starts at prior mean, any other value assumes that row in the parameter file
  
  ni<-dim(ipar)[1]                               #num of item
  theta <- seq(minTheta ,maxTheta ,by=RespInc)   #theta
  nq.maxresp <- length(theta)                    #num of quadratures
  
  
  #drops the recoded score to response from the ipar
  ipar$Response_to_Score <- "1=1;2=2;3=3;4=4;5=5"
  ipar$Score_to_Response <- "1=1;2=2;3=3;4=4;5=5"
  item.par <- ipar[,1:(ncol(ipar)-1)]
  model<-1
  D<-1.0
  simulateTheta<-F
  nSimulee<-0
  popMean<-0
  popSD<-1
  eapFullLength<-T
  ScoreInc <- 0.1 # 0.1=default
  inc=ScoreInc
  topN=1
  CATmethod <- 3 # 3=default (MPWI); 1=MFI, 2=MLWI, 3=MPWI, 4=MEI, 5=MEPV, 6=MEPWI 7:11=do not use herein
  selection.method=CATmethod
  interim.Theta=1
  se.method=1
  StartItem <- 0 # 0=default, starts at prior mean, any other value assumes that row in the parameter file
  first.item.selection=ifelse(StartItem==0,1,3)
  first.at.theta=0
  first.item=StartItem
  show.Theta.Audit.Trail=F
  plot.usage=F
  plot.info=F
  plot.prob=F
  add.final.theta=F
  bank.diagnosis=F
  prior.dist=1
  prior.mean=EAP_Mean
  prior.sd=EAP_SD
  file.items.used=""
  file.theta.history=""
  file.se.history=""
  file.final.theta.se=""
  file.other.thetas=""
  file.likelihood.dist=""
  file.posterior.dist=""
  file.matrix.info=""
  file.full.length.theta=""
  file.selected.item.resp=""
  #max.prob.resp becomes resp data
  
  skipped_items=c()
  
  #response data is the max probabilities file with 1s, etc. & response by theta
  resp.data=as.matrix(responses)
  colnames(resp.data) <- paste0("R",1:ni)
  theta<-seq(minTheta,maxTheta,ScoreInc)
  #nq= # of score increments, 81
  nq=length(theta);
  if (first.item.selection==2 && first.at.theta>=minTheta && first.at.theta<=maxTheta) {
    start.theta<-first.at.theta;
  } else start.theta<-prior.mean;
  if (prior.dist==1) {
    prior<-dnorm((theta-prior.mean)/prior.sd);
  } else if (prior.dist==2) {
    prior<-exp((theta-prior.mean)/prior.sd)/(1+exp((theta-prior.mean)/prior.sd))^2;
  } else prior<-dnorm(theta);
  
  #each examinee = level of 1s, 2s, 3s & theta quadrature
  nExaminees<-dim(resp.data)[1];
  items.used<-matrix(NA,nExaminees,maxNI);
  selected.item.resp<-matrix(NA,nExaminees,maxNI);
  ni.administered<-numeric(nExaminees);
  theta.CAT<-rep(NA,nExaminees);
  sem.CAT<-rep(NA,nExaminees);
  theta.history<-matrix(NA,nExaminees,maxNI);
  se.history<-matrix(NA,nExaminees,maxNI);
  #how is this used?
  posterior.matrix<-matrix(NA,nExaminees,nq);
  #and how is this used?
  LH.matrix<-matrix(NA,nExaminees,nq);
  #item parameters
  NCAT<-item.par[,"NCAT"];
  DISC<-item.par[,"a"];
  CB<-item.par[paste("cb",1:(maxCat-1),sep="")];
  
  missing_items<-c()
  
  prep.prob.info<-function(){
    #make an array with nq=score increments, ni=number of items & number of categories & maxCat= number of response categories
    pp<-array(0,c(nq,ni,maxCat));
    #make a matrix with 0's with score increments and number of items
    matrix.info<-matrix(0,nq,ni);
    #we only used model=1, so not sure what model=2 is
    if (model==1) {
      #cycle through items
      for (i in 1:ni) {
        #ps matrix is a bunch of zeros with score increments as rows and an item's response categories as row
        ps<-matrix(0,nq,NCAT[i]+1);
        ps[,1]<-1;
        ps[,NCAT[i]+1]<-0;
        #what is in this matrix?
        for (k in 1:(NCAT[i]-1)) {
          ps[,k+1]<-1/(1+exp(-D*DISC[i]*(theta-CB[i,k])));
        }
        #fill in info matrix
        for (k in 1:NCAT[i]) {
          pp[,i,k]<-ps[,k]-ps[,k+1];
          matrix.info[,i]<-matrix.info[,i]+(D*DISC[i]*(ps[,k]*(1-ps[,k])-ps[,k+1]*(1-ps[,k+1])))^2/pp[,i,k];
        }
      }
    } 
    #unused matrix
    else if (model==2) {
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
        AX<-rep(0,nq); BX<-rep(0,nq);
        for (k in 1:NCAT[i]) {
          pp[,i,k]<-zz[,k]/den;
          AX<-AX+k^2*pp[,i,k];
          BX<-BX+k*pp[,i,k];
        }
        matrix.info[,i]<-D^2*DISC[i]^2*(AX-BX^2);
      }
    }
    list(pp=pp,matrix.info=matrix.info);
  }
  #matrix.prob.info is the matrix with prep probabilities & info objects
  matrix.prob.info<-prep.prob.info();
  pp<-matrix.prob.info$pp;
  matrix.info<-matrix.prob.info$matrix.info;
  calcInfo<-function(th) {
    info<-numeric(ni);
    if (model==1) {
      for (i in 1:ni) {
        if (items.available[i]==TRUE) {
          ps<-numeric(NCAT[i]+1);
          ps[1]<-1;
          ps[NCAT[i]+1]<-0;
          for (k in 1:(NCAT[i]-1)) {
            ps[k+1]<-1/(1+exp(-D*DISC[i]*(th-CB[i,k])));
          }
          prob<-numeric(NCAT[i]);
          for (k in 1:NCAT[i]) {
            prob[k]<-ps[k]-ps[k+1];
            info[i]<-info[i]+(D*DISC[i]*(ps[k]*(1-ps[k])-ps[k+1]*(1-ps[k+1])))^2/prob[k];
          }
        }
      }
    } else if (model==2) {
      for (i in 1:ni) {
        if (items.available[i]==TRUE) {
          zz<-numeric(NCAT[i]);
          sdsum<-0; den<-0;
          cb<-unlist(CB[i,]); cb<-c(0,cb);
          for (k in 1:(NCAT[i])) {
            sdsum<-sdsum+cb[k];
            zz[k]<-exp(D*DISC[i]*(k*th-sdsum));
            den<-den+zz[k];
          }
          AX<-0; BX<-0;
          prob<-numeric(NCAT[i]);
          for (k in 1:NCAT[i]) {
            prob[k]<-zz[k]/den;
            AX<-AX+k^2*prob[k];
            BX<-BX+k*prob[k];
          }
          info[i]<-D^2*DISC[i]^2*(AX-BX^2);
        }
      }
    }
    return(info);
  }
  
  calc.PW.info<-function(pos) {
    info<-numeric(ni);
    info<-apply(matrix.info*pos,2,sum);
    info[items.available==FALSE]<-0;
    return(info);
  }
  
  select.maxInfo<-function() {
    item.selected<-info.index[1];
    #print("first");
    return (item.selected);}
  
  select.maxInfo_missing<-function() {
    skipped_item<-NA
    print("second");
    #info.index<-info.index[!info.index %in% skipped_items]
    #info.index<-unique(c(skipped_items, info.index))
    if((ni.given+1) %in% missing_items){
      item.selected<-info.index[2]
      skipped_item<-info.index[1]
      print(paste((ni.given+1),rownames(ipar)[item.selected],"skipped", rownames(ipar)[info.index[1]]))
      print(skipped_item)
      print(info.index)}
    
    else{
      item.selected<-info.index[1];
      print(paste((ni.given+1),rownames(ipar)[item.selected]));}
    
    return (c(item.selected,skipped_item));}
  
  calcSE<-function(examinee,ngiven,th) {
    info<-0;
    if (model==1) {
      for (i in 1:ngiven) {
        itm<-items.used[examinee,i];
        ps<-numeric(NCAT[itm]+1);
        ps[1]<-1;
        ps[NCAT[itm]+1]<-0;
        for (k in 1:(NCAT[itm]-1)) {
          ps[k+1]<-1/(1+exp(-D*DISC[itm]*(th-CB[itm,k])));
        }
        prob<-numeric(NCAT[itm]);
        for (k in 1:NCAT[itm]) {
          prob[k]<-ps[k]-ps[k+1];
          info<-info+(D*DISC[itm]*(ps[k]*(1-ps[k])-ps[k+1]*(1-ps[k+1])))^2/prob[k];
        }
      }
    } else if (model==2) {
      info<-0;
      for (i in 1:ngiven) {
        itm<-items.used[examinee,i];
        cb<-unlist(CB[itm,]);
        cb<-c(0,cb);
        zz<-numeric(NCAT[itm]);
        sdsum<-0; den<-0;
        
        for (k in 1:NCAT[itm]) {
          sdsum<-sdsum+cb[k];
          zz[k]<-exp(D*DISC[itm]*(k*th-sdsum));
          den<-den+zz[k];
        }
        AX<-0; BX<-0;
        prob<-numeric(NCAT[itm]);
        for (k in 1:NCAT[itm]) {
          prob[k]<-zz[k]/den;
          AX<-AX+k^2*prob[k];
          BX<-BX+k*prob[k];
        }
        info<-info+D^2*DISC[itm]^2*(AX-BX^2);
      }
    }
    SEM<-1/sqrt(info);
    return(SEM);
  }
  calcEAP<-function (examinee,ngiven) {
    LH<-rep(1,nq);
    for (i in 1:ngiven) {
      item<-items.used[examinee,i];
      resp<-resp.data[examinee,paste("R",item,sep="")];
      prob<-pp[,item,resp];
      LH<-LH*prob;
    }
    posterior<-prior*LH;
    EAP<-sum(posterior*theta)/sum(posterior);
    if (se.method==1) {
      SEM<-sqrt(sum(posterior*(theta-EAP)^2)/sum(posterior));
    } else if (se.method==2) {
      SEM<-calcSE(examinee,ngiven,EAP);
    }
    return(list(THETA=EAP,SEM=SEM,LH=LH,posterior=posterior));
  }
  
  
  #### below is the actual CAT...
  #below is the PROMIS Standard MPWI method

  if (selection.method==3) {
    #through each examinee
    for (j in 1:nExaminees) {
      critMet<-FALSE;
      #vector of items to be administered, based on number of items - all initially set to "TRUE"
      items.available<-rep(TRUE,ni);
      #if in response data there's a missing item response, don't administer by declaring it "FALSE"
      items.available[is.na(resp.data[j,paste("R",1:ni,sep="")])]<-FALSE;
      #either administer the CAT stopping rule max number of items or items available 
      max.to.administer<-ifelse(sum(items.available)<=maxNI,sum(items.available),maxNI);
      ni.given<-0;
      #unclear what first.item.selection == 4 means, but maybe is to allow for existing prior?
      if (first.item.selection==4) theta.current<-ext.theta$theta[j]
      #set theta at prior, theta=0
      else theta.current<-start.theta
      #set the posterior distribution equal to prior (density of normal dist of theta, 0m sd=10)
      posterior<-prior;
      while (critMet==FALSE && ni.given<max.to.administer) {
        #using the posterior, assign a information value to each item that can be administered
        array.info<-calc.PW.info(posterior);
        ni.available<-sum(array.info>0); 
        #sort items based on information
        info.index<-rev(order(array.info));
        item.selected<-select.maxInfo();
        #below is item selection based on existing prior or defined start items
        if (ni.given==0) {
          if (first.item.selection==3 && first.item>=1 && first.item<=ni) {
            if (items.available[first.item]==TRUE) {
              item.selected<-first.item
            }
          } else if (first.item.selection==2 || first.item.selection==4) {
            array.info<-calcInfo(theta.current);
            info.index<-rev(order(array.info));
            item.selected<-select.maxInfo();
          }
        }
        #create a new dataframe based on the response data
        resp<-resp.data[j,paste("R",item.selected,sep="")];
        #below logs the responses
        prob<-pp[,item.selected,resp];
        posterior<-posterior*prob;
        ni.given<-ni.given+1;
        items.used[j,ni.given]<-item.selected;
        items.available[item.selected]<-FALSE;
        selected.item.resp[j,ni.given]<-resp.data[j,paste("R",item.selected,sep="")];
        estimates<-calcEAP(j,ni.given);
        theta.history[j,ni.given]<-estimates$THETA;
        se.history[j,ni.given]<-estimates$SEM;
        theta.current<-estimates$THETA;
        if (ni.given>=max.to.administer || (estimates$SEM<=maxSE && ni.given>=minNI)) {
          critMet<-TRUE;
          theta.CAT[j]<-estimates$THETA;
          sem.CAT[j]<-estimates$SEM;
          LH.matrix[j,]<-estimates$LH;
          posterior.matrix[j,]<-estimates$posterior;
          ni.administered[j]<-ni.given;
        }
      }
      if (show.Theta.Audit.Trail) plot.theta.audit.trail();
    }
  }

  ################################################################################
  # 6. prepare output
  ################################################################################

  #CAT_items<-apply(se.history, 1, function(x){list(1:ifelse(x[12]<0.3, (min(which(x[4:12]<0.3))+3),12))})
CAT_items<-lapply(1:dim(se.history)[1], function(x){1:ifelse(se.history[x,12]<0.3, (min(which(se.history[x,4:12]<0.3))+3),12)})

return(
  as.data.frame(
   t(mapply(function(x,y){
      c("Tscore"=(theta.history[x,y]*10+50),
        "SE"=(se.history[x,y]*10),
        "NI"=y)},
           1:nrow(se.history),
           sapply(CAT_items, max)))))
  
  #final.theta.se<-data.frame(
  #  "Theta"=theta.history[max(CAT_items)],
  #  "SE"=se.history[max(CAT_items)],
  #  "NI"=max(CAT_items));
    #row.names = rownames(ipar)[items.used[CAT_items]]);
  ###problem^ missing only works for a single response- it'll have to redone to accoutn for multiple  
  ################################################################################
  ### Missing items CAT below ###
  ################################################################################
  
  #Redo the CAT, after initial completion, with items missing
  if(missing!=F){
    #### below is the actual CAT...
    #below is the PROMIS Standard MPWI method
    
    missing_items=sample(CAT_items, missing)
    
    if (selection.method==3) {
      #through each examinee
      for (j in 1:nExaminees) {
        critMet<-FALSE;
        #vector of items to be administered, based on number of items - all initially set to "TRUE"
        items.available<-rep(TRUE,ni);
        #if in response data there's a missing item response, don't administer by declaring it "FALSE"
        items.available[is.na(resp.data[j,paste("R",1:ni,sep="")])]<-FALSE;
        #either administer the CAT stopping rule max number of items or items available 
        max.to.administer<-ifelse(sum(items.available)<=maxNI,sum(items.available),maxNI);
        ni.given<-0;
        #unclear what first.item.selection == 4 means, but maybe is to allow for existing prior?
        if (first.item.selection==4) theta.current<-ext.theta$theta[j]
        #set theta at prior, theta=0
        else theta.current<-start.theta
        #set the posterior distribution equal to prior (density of normal dist of theta, 0m sd=10)
        posterior<-prior;
        while (critMet==FALSE && ni.given<max.to.administer) {
          #using the posterior, assign a information value to each item that can be administered
          array.info<-calc.PW.info(posterior);
          ni.available<-sum(array.info>0); 
          #sort items based on information
          info.index<-rev(order(array.info));
          info.index<-info.index[!info.index %in% skipped_items]
          item.selected<-select.maxInfo_missing()[1];
          skipped_items<-c(skipped_items,select.maxInfo_missing()[2]);
          #below is item selection based on existing prior or defined start items
          if (ni.given==0) {
            if (first.item.selection==3 && first.item>=1 && first.item<=ni) {
              if (items.available[first.item]==TRUE) {
                item.selected<-first.item
              }
            } else if (first.item.selection==2 || first.item.selection==4) {
              array.info<-calcInfo(theta.current);
              info.index<-rev(order(array.info));
              item.selected<-select.maxInfo_missing()[1];
              skipped_items<-c(skipped_items,select.maxInfo_missing()[2]);
            }
          }
          #create a new dataframe based on the response data
          resp<-resp.data[j,paste("R",item.selected,sep="")];
          #below logs the responses
          prob<-pp[,item.selected,resp];
          posterior<-posterior*prob;
          ni.given<-ni.given+1;
          items.used[j,ni.given]<-item.selected;
          items.available[item.selected]<-FALSE;
          selected.item.resp[j,ni.given]<-resp.data[j,paste("R",item.selected,sep="")];
          estimates<-calcEAP(j,ni.given);
          theta.history[j,ni.given]<-estimates$THETA;
          se.history[j,ni.given]<-estimates$SEM;
          theta.current<-estimates$THETA;
          if (ni.given>=max.to.administer || (estimates$SEM<=maxSE && ni.given>=minNI)) {
            critMet<-TRUE;
            theta.CAT[j]<-estimates$THETA;
            sem.CAT[j]<-estimates$SEM;
            LH.matrix[j,]<-estimates$LH;
            posterior.matrix[j,]<-estimates$posterior;
            ni.administered[j]<-ni.given;
          }
        }
        if (show.Theta.Audit.Trail) plot.theta.audit.trail();
      }
    }
    
    final.theta.se<-data.frame(
      "Theta"=theta.history[CAT_items],
      "SE"=se.history[CAT_items],
      row.names = rownames(ipar)[items.used[CAT_items]]);
    print(sort(missing_items))
    
    
  }
  
  #return(final.theta.se)
  
}

#######################################################################################

#ipar2<-Item_Parameters$Depression

#responses<-read.csv("temp.csv", row.names=1)

#CAT_sim_missing(ipar, responses)
#CAT_sim_missing(ipar, responses,3)

