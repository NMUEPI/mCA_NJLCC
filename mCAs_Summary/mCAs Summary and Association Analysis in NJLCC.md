# **mCAs Summary and Association Analysis in NJLCC**

    NJLCC <- read.csv("NJLCC_mCA_upload.csv")

## **Counts and prevalence of mosaic chromosomal alterations (mCAs) by event type and genomic location among Nanjing Lung Cancer Cohort (NJLCC) participants**

    func1=function(type){
       sumdata=data.frame(Copy_number_change=type,No.of.events_Overall_populations=paste(sum(NJLCC[,get(type)])," (",sum(NJLCC[,get(paste0(type,"g",sep=""))]),", ",paste0(sprintf("%1.2f",sum(NJLCC[,get(paste0(type,"g",sep=""))])/nrow(NJLCC)*100),"%"),")",sep=""),Proportion_Overall_populations=paste0(sprintf("%1.2f",sum(NJLCC[,get(type)])/sum(NJLCC[,All])*100),"%",sep=""),No.of.events_case=paste(sum(NJLCC[LC==1,get(type)])," (",sum(NJLCC[LC==1,get(paste0(type,"g",sep=""))]),", ",paste0(sprintf("%1.2f",sum(NJLCC[LC==1,get(paste0(type,"g",sep=""))])/nrow(NJLCC[LC==1,])*100),"%"),")",sep=""),Proportion_case=paste0(sprintf("%1.2f",sum(NJLCC[LC==1,get(type)])/sum(NJLCC[LC==1,All])*100),"%",sep=""),No.of.events_control=paste(sum(NJLCC[LC==0,get(type)])," (",sum(NJLCC[LC==0,get(paste0(type,"g",sep=""))]),", ",paste0(sprintf("%1.2f",sum(NJLCC[LC==0,get(paste0(type,"g",sep=""))])/nrow(NJLCC[LC==0,])*100),"%"),")",sep=""),Proportion_control=paste0(sprintf("%1.2f",sum(NJLCC[LC==0,get(type)])/sum(NJLCC[LC==0,All])*100),"%",sep=""))

    }
    type <- c("All","gain","loss","loh","undermine")
    cl <- makeCluster(5,type="FORK")
    res <- parLapply(cl,type,func1)
    sumdata <- do.call('rbind',res)
    stopCluster(cl)
    # write.csv(sumdata,"table1.csv",row.names=F)

## **Association Analysis**

    NJLCC$Allg=ifelse(NJLCC$All==0,0,1)#mCA
    NJLCC$lohg=ifelse(NJLCC$loh==0,0,1)#CN_LOH
    NJLCC$gaing=ifelse(NJLCC$gain==0,0,1)#Gain
    NJLCC$lossg=ifelse(NJLCC$loss==0,0,1)#Loss
    NJLCC$undermineg=ifelse(NJLCC$undermine==0,0,1)#Undetermined

    NJLCC$lossg_r=NJLCC$lossg
    NJLCC[NJLCC$lossg!=1 & NJLCC$Allg==1,]$lossg_r=NA

    NJLCC$lohg_r=NJLCC$lohg
    NJLCC[NJLCC$lohg!=1 & NJLCC$Allg==1,]$lohg_r=NA

    NJLCC$gaing_r=NJLCC$gaing
    NJLCC[NJLCC$gaing!=1 & NJLCC$Allg==1,]$gaing_r=NA

    func2=function(type){
     tmod0<-glm(LC~get(type),data=NJLCC,fam=binomial)
     asso_res<-data.frame(stra=type,Incident_lung_cancer=paste(nrow(NJLCC[LC==1 & get(type)!=0,]),"/",nrow(NJLCC[LC==1 & !is.na(get(type)),])," (",scales::percent(nrow(NJLCC[LC==1 & get(type)!=0,])/nrow(NJLCC[LC==1 & !is.na(get(type)),]),0.01),")",sep=""),Incident_without_lung_cancer=paste(nrow(NJLCC[LC==0 & get(type)!=0,]),"/",nrow(NJLCC[LC==0 & !is.na(get(type)),])," (",scales::percent(nrow(NJLCC[LC==0 & get(type)!=0,])/nrow(NJLCC[LC==0 & !is.na(get(type)),]),0.01),")",sep=""),OR_CI_glm0=paste(round(exp(summary(tmod0)$coef[2,1]),2)," ","(",round(exp(confint(tmod0)[2,1]),2),",",round(exp(confint(tmod0)[2,2]),2),")",sep=""),P_glm0=summary(tmod0)$coef[2,4])
     return(asso_res)
    }
    type <- c("Allg","lohg_r","gaing_r","lossg_r")
    cl <- makeCluster(5,type="FORK")
    res <- parLapply(cl,type,func2)
    assodata <- do.call('rbind',res)
    stopCluster(cl)
    # write.csv(assodata,"asso_mca.csv",row.names=F)

