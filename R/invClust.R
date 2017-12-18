invClust<-function(roi,wh=1,geno,annot,SNPtagg="n",SNPsel=1:ncol(geno),method=1,dim=1,pc=0,ngroups=1, ...)
{
  if(is.character(roi))
    roidat<-read.table(roi,header=TRUE,as.is=TRUE)

  if(is.data.frame(roi))
    roidat<-roi  

  tagg<-FALSE
  if(SNPsel[1]=="tagg")
  {
    ssel<-strsplit(strsplit(roidat[wh,6],";")[[1]],",")
    SNPsel<-which(colnames(geno)%in%unlist(ssel))
    snpsel<-lapply(ssel, function(x) x[x%in%colnames(geno)])

    tagg<-TRUE
  }

  if(method==1)
  {
    invcallGenos<-invClustEM(roi=roi,wh=wh ,geno=geno[,SNPsel],annot=annot[SNPsel,],dim=dim,pc=pc,ngroups=ngroups,...)
    invs<-apply(invcallGenos["genotypes"],1,function(x) which(max(x)==x)[1])-1
    unc1<-1-apply(invcallGenos["genotypes"],1,function(x) max(x))

  
    if(SNPtagg=="y")
    {    
      ssel<-strsplit(strsplit(roidat[wh,6],";")[[1]],",")
      snpsel<-lapply(ssel, function(x) x[x%in%colnames(geno)])


      if(is.na(snpsel)[1])
      { 
        warning("No reference SNPs available in ROI file retuning cumputed clusters with no reference")
        return(invcallGenos)
      }

      invs<-apply(invcallGenos["genotypes"],1,function(x) which(max(x)==x)[1])-1

      genoinvMat<-matrix(as.raw(invs+1),ncol=1)
      rownames(genoinvMat)<-rownames(geno)
      colnames(genoinvMat)<-c("genotype1")
      inversionGenotypes<-new("SnpMatrix",genoinvMat)
  
      LD<-sapply(1:length(snpsel), function(y) mean(ld(inversionGenotypes,geno[,snpsel[[y]]],stats="R")))

      #relabel with first group of SNPs tagging inverted allele.
      if(LD[1]<0)    
      invs<- -invs+2                  

      attr(invcallGenos,"genoRef")<-list(invs=invs,unc=unc1)
    }
       
    out<-invcallGenos
    out
  }

  if(method==2)
  {  
     if(!tagg)
     {
       invcallGenos<-invClustEM(roi=roi,wh=wh ,geno=geno[,SNPsel],annot=annot[SNPsel,],pc=pc,dim=2)

       invs<-apply(invcallGenos["genotypes"],1,function(x) which(max(x)==x)[1])-1
       unc1<-1-apply(invcallGenos["genotypes"],1,function(x) max(x))

       sel2<-invs==2
       sel1<-invs==1
       sel0<-invs==0


       inv_0<-invClustEM(roi=roi,wh=wh ,geno=geno[sel0,SNPsel],annot=annot[SNPsel,],pc=0,dim=2)
       i0<-apply(inv_0["genotypes"],1,function(x) which(max(x)==x)[1])
       i0unc<-1-apply(inv_0["genotypes"],1,function(x) max(x))

       suppressWarnings(inv_11<-Mclust(invcallGenos$datin$y[sel1,1],G=c(1,2),modelNames="E"))
       suppressWarnings(inv_12<-Mclust(invcallGenos$datin$y[sel1,2],G=c(1,2),modelNames="E"))
      
       i1<-inv_11$classification
       i1unc<-inv_11$uncertainty

       BD<-(inv_11$BIC[2]-inv_11$BIC[1])<(inv_12$BIC[2]-inv_12$BIC[1])
       if(inv_11$BIC[1])
       {
         i1<-inv_12$classification
         i1unc<-inv_12$uncertainty
       }

       inv_2<-invClustEM(roi=roi,wh=wh ,geno=geno[sel2,SNPsel],annot=annot[SNPsel,],pc=0,dim=2)
       i2<-apply(inv_2["genotypes"],1,function(x) which(max(x)==x)[1])
       i2unc<-1-apply(inv_2["genotypes"],1,function(x) max(x))

       ind<-as.integer(c(length(table(i0)),length(table(i1)),length(table(i2))))
       #at least 2 subclusters in heterrozygous and 3 clusters in one of the homozygous
       suitable<-sum(c(2,3)%in%ind)

       BIC0diff<-inv_0$EMestimate$bicDiff
       BIC2diff<-inv_2$EMestimate$bicDiff

       if(suitable==2)
       {
                           
         if(BIC2diff>BIC0diff)
         {
           i0<-rep(1,length(i0))
           ind<-c(1,2,3)
           i0unc<-Mclust(invcallGenos$datin$y[sel0,1],G=c(1),modelNames="E")$uncertainty
         }else{ 
           i2<-rep(1,length(i2))
           ind<-c(3,2,1)
           i2unc<-Mclust(invcallGenos$datin$y[sel2,1],G=c(1),modelNames="E")$uncertainty
         } 

         unc2<-invs
         unc2[sel0]<-i0unc
         unc2[sel1]<-i1unc
         unc2[sel2]<-i2unc

         if(ind[1]==3)
         {
           sel0temp<-sel2
           sel2<-sel0
           sel0<-sel0temp

           i0temp<-i2
           i2<-i0
           i0<-i0temp

         }

         whmin1<-tapply(invcallGenos$datin$y[sel1,2],as.factor(i1),mean) 
         whmin1<-order(whmin1)
         sel11<-i1==whmin1[1] 
         sel12<-i1==whmin1[2] 


         whmin2<-tapply(invcallGenos$datin$y[sel2,2],as.factor(i2),mean) 
         whmin2<-order(whmin2)
         sel21<-i2==whmin2[1] 
         sel22<-i2==whmin2[2] 
         sel23<-i2==whmin2[3] 
  
         inv6<-invs
         inv6[sel0]<-1
         inv6[sel1][sel11]<-4
         inv6[sel1][sel12]<-2
         inv6[sel2][sel21]<-6
         inv6[sel2][sel22]<-5
         inv6[sel2][sel23]<-3

         unc<-apply(cbind(unc1,unc2),1,max)

       }else{
        warning("method 2 do not seem suitable for this inversion returning method 1\n")
        out<-invcallGenos
     
        return(out)
       }


     }else{ #call the genotypes by clusterign with tagg snps

         snp1<-snpsel[[1]]
         invcallGenos1<-invClustEM(roi=roi,wh=wh ,geno=geno[,snp1],annot=annot[annot$name%in%snp1,],pc=pc,dim=2,...)
         invs1<-apply(invcallGenos1["genotypes"],1,function(x) which(max(x)==x)[1])-1

         snp2<-snpsel[[2]]
         invcallGenos2<-invClustEM(roi=roi,wh=wh ,geno=geno[,snp2],annot=annot[annot$name%in%snp2,],pc=pc,,dim=2,...)
         invs2<-apply(invcallGenos2["genotypes"],1,function(x) which(max(x)==x)[1])-1

         snp3<-snpsel[[3]]
         invcallGenos3<-invClustEM(roi=roi,wh=wh ,geno=geno[,snp3],annot=annot[annot$name%in%snp2,],pc=pc,dim=2,...)
         invs3<-apply(invcallGenos3["genotypes"],1,function(x) which(max(x)==x)[1])-1

         invcallGenos<-invClustEM(roi=roi,wh=wh ,geno=geno,annot=annot,dim=2,pc=pc,...)
         invs<-apply(invcallGenos["genotypes"],1,function(x) which(max(x)==x)[1])-1
         
         gn<-paste(invs1,invs2,invs3,sep="|")

         ind<-c(1,3,7,9)
         tb<-table(invs1,invs2) 
         sel3<-switch(which(tb[ind]==0),
                      invs1==2&invs2==2, 
                      invs1==0&invs2==2, 
                      invs1==2&invs2==0,
                      invs1==1&invs2==1) 

         tb<-table(invs2,invs3) 
         sel1<-switch(which(tb[ind]==0),
                      invs2==2&invs3==2, 
                      invs2==0&invs3==2, 
                      invs2==2&invs3==0,
                      invs2==1&invs3==1) 


         tb<-table(invs1,invs3) 
         sel2<-switch(which(tb[ind]==0),
                      invs1==2&invs3==2, 
                      invs1==0&invs3==2, 
                      invs1==2&invs3==0,
                      invs1==1&invs3==1) 



         inv6<-rep(NA,length(invs))
         inv6[sel3]<-6
         inv6[sel1]<-1
         inv6[sel2]<-3

         inv6[invs1==1 & invs2==1]<-2
         inv6[invs1==1 & invs3==1]<-4
         inv6[invs2==1 & invs3==1]<-5

         repcall<-sel3 + sel1 + sel2 + (invs1==1 & invs2==1)+ (invs1==1 & invs3==1) + (invs2==1 & invs3==1)

         inv6[repcall!=1]<-NA 
         unc<- as.numeric(repcall==1)

     }
  
    if(SNPtagg=="y" & !tagg)
    {    
      ssel<-strsplit(strsplit(roidat[wh,6],";")[[1]],",")
      snpsel<-lapply(ssel, function(x) x[x%in%colnames(geno)])



      if(is.na(snpsel)[1])
      { 
        out<-invcallGenos

        warning("No reference SNPs available in ROI file retuning cumputed clusters with no reference")
        return(out)
      }

      ginv<-inv6
      genotype1<-rep(0,length(ginv))
      genotype1[ginv%in%c(1)]<-2
      genotype1[ginv%in%c(2,4)]<-1

      genotype2<-rep(0,length(ginv))
      genotype2[ginv%in%c(3)]<-2
      genotype2[ginv%in%c(5,2)]<-1

      genotype3<-rep(0,length(ginv))
      genotype3[ginv%in%c(6)]<-2
      genotype3[ginv%in%c(5,4)]<-1


      genoinvMat<-matrix(as.raw(c(genotype1,genotype2,genotype3)+1),ncol=3)
      rownames(genoinvMat)<-rownames(geno)
      colnames(genoinvMat)<-c("genotype1","genotype2","genotype3")
      inversionGenotypes<-new("SnpMatrix",genoinvMat)
  
      LD<-lapply(1:3, function(y)
      {    
        LD1<-ld(inversionGenotypes,geno[,snpsel[[y]]],stats="R.squared")
        rmna<-colSums(apply(LD1,2,is.na))<1
        tb<-apply(matrix(LD1[,rmna],nrow=3),2,function(x) order(-x)[1])   
        names(table(tb))[order(-table(tb))][1]
      }) 
     
      ldfac<-factor(unlist(LD))
      levels(ldfac)<-c("1","3","6")
      ldfac<-as.numeric(as.character(ldfac))


      #relabel homozygous
      inv6New<-inv6
      inv6New[inv6==ldfac[1]]<-6 
      inv6New[inv6==ldfac[2]]<-1 
      inv6New[inv6==ldfac[3]]<-3 


      cl3<-which(unlist(LD)==6) #which classified group correponds to tagg group 3
      inv6New[inv6==cl3]<-6 

      #relabel heterozygous
      mn1<-apply(invcallGenos$datin$y[inv6New==1,],2,mean)
      mn3<-apply(invcallGenos$datin$y[inv6New==3,],2,mean)
      mn6<-apply(invcallGenos$datin$y[inv6New==6,],2,mean)

      gr2<-apply(invcallGenos$datin$y[inv6==2,],2,mean)  
      gr4<-apply(invcallGenos$datin$y[inv6==4,],2,mean)  
      gr5<-apply(invcallGenos$datin$y[inv6==5,],2,mean)  
      matgr<-cbind(gr2,gr4,gr5)

      mm<-(mn1+mn3)/2-matgr
      dist<-apply(mm,2,function(x) sqrt(sum(x^2)))
      selind<-which(dist==min(dist))
      indhet<-c(2,4,5)
      inv6New[inv6==indhet[selind]]<-2

      selothers<-which(indhet%in%indhet[-selind])
      indhet<-indhet[selothers]

      mm<-(mn1+mn6)/2-matgr[,selothers]
      dist<-apply(mm,2,function(x) sqrt(sum(x^2)))
      selind<-which(dist==min(dist))
      inv6New[inv6==indhet[selind]]<-4

      selind<-which(dist==max(dist))
      inv6New[inv6==indhet[selind]]<-5

      inv6<-inv6New
    }

    out<-invcallGenos
    attr(out,"genoRef")<-list(invs=inv6,unc=unc)
    
  }
   
   out
}


invClustEM<-function(roi,wh,geno,annot,dim,pc,ngroups,...)
{

  #read candidate inversions from file
  if(file.exists(as.character(roi)[1]))
  {
    cat("\n file <", roi, "> found, reading data from file. \n")

    sd<-read.table(roi,header=TRUE,as.is=TRUE)
  }else{
   if(length(roi)==4)
   { 
      cat("\n analizing region chr", roi[1,1],":",as.character(roi[1,4]),"\n")
      sd<-roi
   }else{
      stop("\n unsupported roi format!")
   }
  } 
 
  #read chromosome 
  chr<-sd[wh,1]

  #select break-points
  bp<-unlist(sd[wh,2:3])
 
  snpsInv<-annot[,"chromosome"]==chr & annot[,"position"]>bp[1] & annot[,"position"]<bp[2] 
  genoInv<-geno[,snpsInv]
  snpsum<-col.summary(genoInv)

  #select SNPs
  use<- snpsum$MAF >0.0001 #& snpsum$z.HWE^2 <16
  use[is.na(use)]<-FALSE 
  genoDat <- matrix(as.numeric(genoInv[, use]), ncol=sum(use, na.rm=TRUE))

  #analize regions with more than 10 SNPs
  if(NCOL(genoDat)<10)
  {
    warning("  less than 10 SNPs in region!")
    #res<-NA
    #res
  }

    cat("  computing clustering with",NCOL(genoDat),"SNPs \n")


    #simple imputation of missings 
    impt<-function(x)
    {
      ans<-x
      ans[x==0]<-sample(x[!(x==0)],sum((x==0)),replace=TRUE) 
      ans
    }


    #do multidimentional scaling on imputed data
    genoDatIMP<-sapply(1:NCOL(genoDat),function(x)impt(as.vector(genoDat[,x])))
    genoDatIMP<- -(genoDatIMP-3)
    rownames(genoDatIMP)<-1:NROW(genoDat)

    d <- dist(genoDatIMP)
    fit <- cmdscale(d,eig=TRUE, k=5)
    varexpl<-sum((fit$eig[1:2])^2)/sum((fit$eig)^2)  

    yy<-fit$points/max(fit$points)
    y<-yy[,1:2]    
 


    if(identical(pc,0)) #NULL initial conditions for pca (x)
    {   
        #classify without population stratification
        x<-rep(0,NROW(y))
        ngroups<-1
        muGroup0<-rep(0,ngroups)
        sigmaGroup0<-rep(0,ngroups)
        piGroup0<-rep(1,ngroups)/ngroups

    }else{             #initial conditions for pca (x)
        #classify with population stratification

        x<-pc

        if((pc=="y2")[1])
         x<-y[,2]

        if(length(x)!=NROW(y))  
           stop("different number of subjects between genotypes and pc components!")           
       
        #infer muGroup initial conditions with Mclust  
        #remove warnings of Mclust
#        op<-options()
#        options(warn=-1)   
#          kk<-Mclust(x,1:ngroups)
#        optopns<-op

#        ngroups<-kk$G
#        cat("  clustering",ngroups,"subpopulations \n")


        cl<-kmeans(x,ngroups)$cluster
        fcl<-as.factor(cl)
        muGroup0<-tapply(x,fcl,mean)

      
        #initalize in small variance relative to range per number of groups
#        rng<-max(x)-min(x)

        sigmaGroup0<- tapply(x,fcl,sd)/3 #rep(rng/ngroups/30,ngroups)
        piGroup0<-rep(1,ngroups)/ngroups


    }

 
    if(dim==1) #one MDS component
    {
       y<-y[,1]
       yy<-y
       varexpl<-sum((fit$eig[1])^2)/sum((fit$eig)^2)  
        

       sg<-sign(skewness(y))

        qnt<-quantile(y,0.25)

        if(sg<0)
          qnt<-quantile(y,0.75)
         

        mu10<-rep(qnt,ngroups)
        mu20<-rep(qnt+sg*sd(y)*2,ngroups)

        sigma0<-rep(sd(y)/2,ngroups)

#       rng<-max(y)-min(y)
#       mu10<-rep(min(y)+rng/10,ngroups)
#       mu20<-rep(max(y)-rng/10,ngroups)
#       sigma0<-rep(rng/10,ngroups)


        p0<-rep(0.9,ngroups)
        q0<-rep(0.1,ngroups)

        cat("computing model near null inversion frequency \n")
        EMestimate<-EM1D(y,mu10,mu20,sigma0,p0,q0,x,muGroup0,sigmaGroup0,piGroup0,...)      
        like<-EMestimate$likelihood

        p0<-rep(0.5,ngroups)
        q0<-rep(0.5,ngroups)

        cat("computing model near 50% inversion frequency \n")
        EMestimate2<-EM1D(y,mu10,mu20,sigma0,p0,q0,x,muGroup0,sigmaGroup0,piGroup0,...)      
        like2<-EMestimate2$likelihood

        if(like2[length(like2)]>like[length(like)])
          EMestimate<-EMestimate2

       #test whether one component is more likely
       lk<-sum(log(exp(-(((y-mean(y))/sd(y))^2)/2)/sd(y)/sqrt(2*pi)))
 
       bicME<-(4)*log(length(y))-2*like[length(like)]
       bicOneComp<-(2)*log(length(y))-2*(lk)

    
       if(bicOneComp<bicME & c(pc!="y2")[1])
       {    
          npd<-exp(-(((EMestimate$ysam-mean(y))/sd(y))^2)/2)/sd(y)/sqrt(2*pi)
     
          thetanew<-c(mean(y),0,0,sd(y),1,0,0,1,0)
          EMestimate<-list(thetanew=thetanew,theta=c(mu10,sigma0),likelihood=lk,classInd=cbind(rep(1,length(y)),rep(0,length(y)),rep(0,length(y))),groupInd=rep(1,length(y)),npd=npd,ysam=EMestimate$ysam)
       }

       EMestimate$bicDiff<-bicOneComp-bicME

       mu<-EMestimate$thetanew[1:3]
       w<-c(EMestimate$thetanew[5]^2,EMestimate$thetanew[6]^2,2*EMestimate$thetanew[5]*EMestimate$thetanew[6])
       sds<-EMestimate$thetanew[4]    
 
       #overlap integral
       npd1<-w[1]*exp(-(((EMestimate$ysam-mu[1])/sds)^2)/2)/sds/sqrt(2*pi)
       npd2<-w[2]*exp(-(((EMestimate$ysam-mu[2])/sds)^2)/2)/sds/sqrt(2*pi)
       npd3<-w[3]*exp(-(((EMestimate$ysam-mu[3])/sds)^2)/2)/sds/sqrt(2*pi)
       dy<-diff(EMestimate$ysam)[1]

       nn<-sqrt(sum((npd1+npd2+npd3)^2*dy))   
       EMestimate$quality<-(nn-sqrt(sum((npd1*npd1+npd2*npd2+npd3*npd3)*dy)))/nn



     }else{    #two MDS components

        ##initial conditions for MDS (y)

        if(ngroups==1)
        {
          sg<-sign(skewness(y[,1]))

          qnt<-quantile(y[,1],0.25)

         if(sg<0)
           qnt<-quantile(y[,1],0.75)         

          mu10<-cbind(rep(qnt,ngroups),0)
          mu20<-cbind(rep(qnt+sg*sd(y[,1])*2,ngroups),0)


          p0<-rep(0.9,ngroups)
          q0<-rep(0.1,ngroups)

        }else{

          mu10<-cbind(rep(min(y[,1]),ngroups),rep(min(y[,2]),ngroups))
          mu20<-cbind(rep(max(y[,1]),ngroups),rep(max(y[,2]),ngroups))
 
          mu10<-cbind(rep(min(y[,1]),ngroups),0)
          mu20<-cbind(rep(max(y[,1]),ngroups),0)
       
          p0<-rep(0.8,ngroups)
          q0<-rep(0.2,ngroups)

        }

        sigma0<-lapply(1:ngroups, function(x) diag(2))

       
        cat("computing model near null inversion frequency \n")
        EMestimate<-EM(y,mu10,mu20,sigma0,p0,q0,x,muGroup0,sigmaGroup0,piGroup0,...)
        like<-EMestimate$likelihood        


       #test in a sigle subpopulation model whether one component 
       #or a solution near 50% inversion frequency are more likely
       if(ngroups==1)
       { 
         p0<-rep(0.5,ngroups)
         q0<-rep(0.5,ngroups)

         cat("computing model near 50% inversion frequency \n")
         EMestimate2<-EM(y,mu10,mu20,sigma0,p0,q0,x,muGroup0,sigmaGroup0,piGroup0,...)
         like2<-EMestimate2$likelihood

         if(like2[length(like2)]>like[length(like)])
           EMestimate<-EMestimate2

         #compute the one component pdf
         ysam1<-seq(min(y[,1])-2*sd(y[,1]),max(y[,1])+2*sd(y[,1]),length.out=500)
         ysam2<-seq(min(y[,2])-2*sd(y[,2]),max(y[,2])+2*sd(y[,2]),length.out=500)

         pd<-array(0, dim=c(length(ysam1),length(ysam2)))
         normConst<-0
         ysam<-cbind(rep(ysam1,length(ysam2),),rep(ysam2,each=length(ysam1)))

         mm<-apply(y,2,mean)
         sigma<-cov(y)
         invss<-ginv(sigma)

         yy1<-ysam[,1]-mm[1]
         yy2<-ysam[,2]-mm[2]
         ym<-cbind(yy1,yy2)
         pd<-matrix(exp(-sapply(1:NROW(ym), function(x) ym[x,]%*%(invss%*%ym[x,]))/2)/sqrt(det(sigma))/(2*pi), ncol=length(ysam1))
         dy1<-mean(diff(ysam1))
         dy2<-mean(diff(ysam2))

         normConst<-sum(pd[-1,-1]*dy1*dy2)
         npd<-pd[-1,-1]/normConst

         #compute likelihood of the data
         yy1<-y[,1]-mm[1]
         yy2<-y[,2]-mm[2]
         ym<-cbind(yy1,yy2)

         rr<-exp(-sapply(1:NROW(ym), function(x) ym[x,]%*%(invss%*%ym[x,]))/2)/sqrt(det(sigma))/(2*pi)/normConst

         lk<-sum(log(rr))

         bicME<-(8)*log(length(y))-2*like[length(like)]
         bicOneComp<-(5)*log(length(y))-2*(lk)

         bicME<-(8)*2-2*like[length(like)]
         bicOneComp<-(5)*2-2*(lk)

invs<-apply(EMestimate$classInd,1,function(x) which(max(x)==x)[1])-1

 bt<-sapply(names(table(invs)), function(ii)
 {
   ycl<-y[invs==ii,]
   mn<-apply(ycl,2,mean)
   sum(sapply(1:nrow(ycl), function(x) (ycl[x,]-mn)%*%(ycl[x,]-mn)))/(length(ycl)-1)
 })

 mn<-apply(y,2,mean)
 tot<-sum(sapply(1:nrow(y), function(x) (y[x,]-mn)%*%(y[x,]-mn)))/(length(y)-1)


#         if(bicOneComp<bicME & c(pc!="y2")[1])
    if(1-mean(bt)/tot<0.45 & c(pc!="y2")[1])
         {  

              thetanew<-c(mm,0,0,sigma,1,0,0,1,0)
              EMestimate<-list(thetanew=thetanew,theta=c(mm,sigma),likelihood=lk,classInd=cbind(rep(1,nrow(y)),rep(0,nrow(y)),rep(0,nrow(y))),groupInd=rep(1,nrow(y)),npd=npd,ysam1=ysam1[-1],ysam2=ysam2[-1],xsam=EMestimate$xsam)
           }

          EMestimate$bicDiff<-bicOneComp-bicME
      }

    }


  res<-list(EMestimate=EMestimate,datin=list(x=x,y=yy,ids=rownames(geno),varexpl=varexpl))
  class(res)<-"invClust"
  res

    
}

print.invClust<-function(x,...)
{
  cat("Inversion genotype clustering \n")
  cat("-object of class invClust- \n")
  cat("  fields: $EMestimate: mixture model parameters \n") 
  cat("          $datin: fitted data \n")
  cat("  subjects:", NROW(x$datin$y), "\n")
  cat("  groups fitted:", NCOL(x$EMestimate$groupInd), "\n")
  cat("  overall inversion allele frequency:",mean(x$EMestimate$classInd[,3]+x$EMestimate$classInd[,2]*2)/2, "\n")   
  cat("  variance explained by ",NCOL(x$datin$y)," MDS componet(s):", x$datin$varexpl, "\n") 
}

"[.invClust"<-function(x,i,j)
{
  if(i=="genotypes")
  {
    ret<-round(x$EMestimate$classInd[,c(1,3,2)],3)
    colnames(ret)<-c("NI/NI","NI/I","I/I")
    rownames(ret)<-x$datin$ids
    return(ret)
   }

  if(i=="groups")
  { ret<-round(x$EMestimate$groupInd,3)
    rownames(ret)<-x$datin$ids
    return(ret)
  }
}


invGenotypes <- function(x) 
{
   if (!inherits(x, "invClust"))
     stop("x is not of class 'invClust'")

   temp <- x["genotypes"]
   geno <- apply(temp, 1, which.max)
   p <- x$EMestimate$thetanew$pnew
   q <- x$EMestimate$thetanew$qnew

   if (p>=q)
   {
     out <- factor(geno) 
     levels(out)<-c("NI/NI", "NI/I", "I/I")[1:length(levels(out))]
   }
   
   if (q>p)
   {
     out <- factor(geno) 
     levels(out)<-c("I/I", "NI/I", "NI/NI")[1:length(levels(out))]
   }

   if("genoRef"%in%names(attributes(x)) )
   {
     out<-attr(x,"genoRef")$invs
     attr(out,"unc")<-attr(x,"genoRef")$unc
   } 
   out
}



invQuality <- function(x,...) {
  
  dim <- ncol(x$datin$y)
  theta <- x$EMestimate$thetanew

  if(sum(x$datin$x)!=0) {
     stop(" not implemented for analyisis with population stratification")
  } 

  if (dim == 1) {
   mu <- c(theta$mu1new[1,1], theta$mu3new[1,1], theta$mu2new[1,1])
   sd <- rep(theta$sigmanew, 3)
   w <- c(theta$pnew^2 , 2*theta$pnew*theta$qnew, theta$qnew^2)
   ans <- getQuality(mu, sd, w, ...)
  }

  else if (dim >= 2) {
   mu <-rbind(theta$mu1new, theta$mu3new, theta$mu2new)
   sigma <-theta$sigmanew[[1]]
   w <- c(theta$pnew^2 , 2*theta$pnew*theta$qnew, theta$qnew^2)

   ans <- getQuality2D(mu, sigma, w, ...)
  }

 ans
}



getQuality2D <- function(mu, sigma, w, iter=10000) {
  J <- nrow(mu)
  p <- NULL
  for (j in 1:J) {
      X <- rmvnorm(n = iter, mean = mu[j,], sigma = sigma)
      Y <- sapply(1:J, function(s) w[s] * dmvnorm(X, mu[s,], sigma))
      p <- c(p, mean(apply(Y, 1, which.max) == j))
  }
  out <- sum(p * w)
  out
}

getQuality <- function(mu, sds, w, iter=10000) {
        J <- length(mu)
        p <- NULL
        for (j in 1:J) {
            X <- rnorm(iter, mu[j], sds[j])
            Y <- sapply(1:J, function(s) w[s] * dnorm(X, mu[s],
                sds[s]))
            p <- c(p, mean(apply(Y, 1, which.max) == j))
        }
        out <- sum(p * w)
        out
}

invFreq <- function(x) {
 if ( inherits(x, "try-error"))
  ans <- NA
 else
  ans <- mean(x$EMestimate$classInd[, 
        3] + x$EMestimate$classInd[, 2] * 2)/2
 ans
}



plot.invClust<-function(x,y,wh="yy",...)
{

  pch<-1

  if("genoRef"%in%names(attributes(x))) 
  {
     invgenos<-attr(x,"genoRef")$invs
     pch<-as.character(invgenos)
  } 

  if(NCOL(x$EMestimate$groupInd)==1)
  {
    if(NCOL(x$datin$y)==1)
    {
      if(wh=="yy")
      { 
        d<-hist(x$datin$y,br=50,plot=FALSE);
        plot(x$EMestimate$ysam,x$EMestimate$npd,col="red",type="l",ylim=c(0,max(c(d$density,x$EMestimate$npd))),xlab="First MDS component",ylab="Frequency");
        hist(x$datin$y,freq=FALSE,br=50,add=TRUE)
       }

      if(wh=="like")
      { 
        like<-x$EMestimate$likelihood
        plot(like,xlab="iteration",ylab="loglikelihood")
        lines(1:length(like),like,col="blue",lty=2)
      }

    }else{
      contour(x$EMestimate$ysam1,x$EMestimate$ysam2,x$EMestimate$npd,drawlabels=FALSE,nlevels=100,xlab="First MDS component", ylab="Second MDS component")
      cc<-round(x["genotypes"][,2]+x["genotypes"][,3]*2)+2
      points(x$datin$y[,1],x$datin$y[,2],col=cc,pch=pch)
    #  legend("topleft",legend=c("NI/NI","NI/I","I/I"),col=c(2,3,4),pch=pch) 
    }
  }else{ #stratification was computed
    if(NCOL(x$datin$y)==1)
    {
      contour(x$EMestimate$xsam,x$EMestimate$ysam,t(x$EMestimate$npd),drawlabels=FALSE,nlevels=100,xlab="genome wide pca", ylab="First MDS component")
      #cc<-round(x["genotypes"][,2]+x["genotypes"][,3]*2)+2
      points(x$datin$x,x$datin$y,col=cc,pch=pch)
    #  legend("topleft",legend=c("NI/NI","NI/I","I/I"),col=c(2,3,4),pch=pch) 
    }else{
      if(wh=="yy")
      {
        npdmarg3<-matrix(rep(0,dim(x$EMestimate$npd)[1]*dim(x$EMestimate$npd)[2]),ncol=dim(x$EMestimate$npd)[2])
        for(jj in 1:dim(x$EMestimate$npd)[3])
        {
          npdmarg3<-npdmarg3+x$EMestimate$npd[,,jj]
        }
        contour(x$EMestimate$ysam1,x$EMestimate$ysam2,npdmarg3,drawlabels=FALSE,nlevels=100,xlab="First MDS component", ylab="Second MDS component")
        cc<-round(x["genotypes"][,2]+x["genotypes"][,3]*2)+2
        points(x$datin$y[,1],x$datin$y[,2],col=cc,pch=pch)
     #   legend("topleft",legend=c("NI/NI","NI/I","I/I"),col=c(2,3,4),pch=pch)   
      }

      if(wh=="xy")
      {

        xrng<-c(min(x$datin$x,x$EMestimate$xsam),max(x$datin$x,x$EMestimate$xsam))
        yrng<-c(min(x$datin$y,x$EMestimate$ysam1),max(x$datin$y,x$EMestimate$ysam1))

        npdmarg3<-matrix(rep(0,dim(x$EMestimate$npd)[1]*dim(x$EMestimate$npd)[3]),ncol=dim(x$EMestimate$npd)[3])
        for(jj in 1:dim(x$EMestimate$npd)[2])
        {
          npdmarg3<-npdmarg3+x$EMestimate$npd[,jj,]
        }
        contour(x$EMestimate$xsam,x$EMestimate$ysam1,t(npdmarg3),drawlabels=FALSE,nlevels=100,xlab="First PCA component", ylab="First MDS component",xlim=xrng,ylim=yrng)
        cc<-round(x["genotypes"][,2]+x["genotypes"][,3]*2)+2
        points(x$datin$x,x$datin$y[,1],col=cc,pch=pch)
      #  legend("topleft",legend=c("NI/NI","NI/I","I/I"),col=c(2,3,4),pch=pch)   
      }

    }
  }
}

EM1D<-function(y,mu10,mu20,sigma0,p0,q0,x,muGroup0,sigmaGroup0,piGroup0,tol=10^(-5),it=1000)
{
  
  mu1<-mu10
  mu2<-mu20 
  sigma<-sigma0
  p<-p0
  q<-q0
  mu3<-(mu10+mu20)/2

  muGroup<-muGroup0
  sigmaGroup<-sigmaGroup0
  piGroup<-piGroup0

  lk<-c() 
  i<-0
  dd<-Inf
  while(i<it & tol<dd)
  {
    i<-i+1
    theta<-rbind(mu1,mu2,mu3,sigma,p,q,muGroup,piGroup,sigmaGroup)
    N<-length(y)

    pnew<-c()
    qnew<-c()
    piGroupNew<-c()
    mu1new<-c()
    mu2new<-c()
    mu3new<-c()
    muGroupNew<-c()
    sigmanew<-c()
    sigmaGroupnew<-c()


    #compute partition function given observation i (sum over all groups and genotypes)
    Z<-rep(0,length(y))
    for(gg in 1:length(muGroup0))
    {
      if(sigmaGroup[gg]==0)
      { 
        fx<-1
      }else{
        fx<-exp(-(((x-muGroup[gg])/sigmaGroup[gg])^2)/2)/sigmaGroup[gg]
      }

      f1Group<-piGroup[gg]*p[gg]^2*exp(-(((y-mu1[gg])/sigma[gg])^2)/2)/sigma[gg]*fx
      f2Group<-piGroup[gg]*q[gg]^2*exp(-(((y-mu2[gg])/sigma[gg])^2)/2)/sigma[gg]*fx
      f3Group<-piGroup[gg]*(2*p[gg]*q[gg])*exp(-(((y-(mu1[gg]+mu2[gg])/2)/sigma[gg])^2)/2)/sigma[gg]*fx 
      Z<-Z+f1Group+f2Group+f3Group
    }

    for(gg in 1:length(muGroup0))
    {

      #subpopulation likelihoods
      if(sigmaGroup[gg]==0)
      { 
        fx<-1
      }else{
        fx<-exp(-(((x-muGroup[gg])/sigmaGroup[gg])^2)/2)/sigmaGroup[gg]
      }

      f1Group<-piGroup[gg]*p[gg]^2*exp(-(((y-mu1[gg])/sigma[gg])^2)/2)/sigma[gg]*fx
      f2Group<-piGroup[gg]*q[gg]^2*exp(-(((y-mu2[gg])/sigma[gg])^2)/2)/sigma[gg]*fx
      f3Group<-piGroup[gg]*(2*p[gg]*q[gg])*exp(-(((y-(mu1[gg]+mu2[gg])/2)/sigma[gg])^2)/2)/sigma[gg]*fx

      #responsabilities
      w1<-f1Group/Z
      w2<-f2Group/Z
      w3<-f3Group/Z

      #new adimixture frequencies
      pnew<-c(pnew,sum(2*w1+w3)/sum(w1+w2+w3)/2)
      qnew<-c(qnew,sum(2*w2+w3)/sum(w1+w2+w3)/2)

      #new group admixture
      piGroupNew<-c(piGroupNew,sum(w1+w2+w3)/N)  

      #new means
      muanew<-(sum(y*(w1-w2))*sum(w2-w1)+sum(y*(w1+w2+w3))*sum(w1+w2))/(sum(w1-w2)*sum(w2-w1)+sum(w1+w2+w3)*sum(w1+w2))
      mubnew<-sum((y-muanew)*(w1-w2))/sum(w1+w2)

      mu1newg<-muanew+mubnew
      mu2newg<-muanew-mubnew
      mu3newg<-muanew 

      mu1new<-c(mu1new,mu1newg)
      mu2new<-c(mu2new,mu2newg)
      mu3new<-c(mu3new,mu3newg)


      muGroupNewg<-sum(x*(w1+w2+w3))/sum(w1+w2+w3)
      muGroupNew<-c(muGroupNew,muGroupNewg)
        

      #new variance
      sigmanew<-c(sigmanew,sqrt((sum(w1*(y-mu1newg)^2)+sum(w2*(y-mu2newg)^2)+sum(w3*(y-mu3newg)^2))/sum(w1+w2+w3)))
      sigmaGroupnew<-c(sigmaGroupnew,sqrt((sum(w1*(x-muGroupNewg)^2)+sum(w2*(x-muGroupNewg)^2)+sum(w3*(x-muGroupNewg)^2))/sum(w1+w2+w3)))

    }

    thetanew<-rbind(mu1new,mu2new,mu3new,sigmanew,pnew,qnew,muGroupNew,piGroupNew,sigmaGroupnew)
    
    #compute log-likelihood of the data
    #normalizing constant
    ysam<-seq(min(y)-max(sigma),max(y)+max(sigma),length.out=1000)
    xsam<-seq(min(x)-max(sigmaGroup),max(x)+max(sigmaGroup),length.out=1000)
   
    dx<-mean(diff(xsam))
    dy<-mean(diff(ysam))

    #trick for reducing to one dimension
    if(dx==0)
    {
       dx<-1
       xsam<-c(0,0) 
    }  


    classInd<-matrix(0,nrow=length(y),ncol=3)
    pd<-matrix(0, nrow=length(ysam),ncol=length(xsam))
    groupInd<-c()
    NG<-0
    normConst<-0
    for(gg in 1:length(muGroup0))
    {
      thetaGroup<-thetanew[,gg] 
      
      #adjusted model  
      for(jj in 1:length(xsam))
      { 
         if(thetaGroup[9]==0)
         { 
            fx<-1
         }else{
            fx<-exp(-(((xsam[jj]-thetaGroup[7])/thetaGroup[9])^2)/2)/thetaGroup[9]
         }

        f1<-thetaGroup[8]*thetaGroup[5]^2*exp(-(((ysam-thetaGroup[1])/thetaGroup[4])^2)/2)/thetaGroup[4]*fx
        f2<-thetaGroup[8]*thetaGroup[6]^2*exp(-(((ysam-thetaGroup[2])/thetaGroup[4])^2)/2)/thetaGroup[4]*fx
        f3<-thetaGroup[8]*(2*thetaGroup[5]*thetaGroup[6])*exp(-(((ysam-thetaGroup[3])/thetaGroup[4])^2)/2)/thetaGroup[4]*fx
 
      
        pd[,jj]<-pd[,jj]+f1+f2+f3
      }
      

     
      #observed likelihood
      if(thetaGroup[9]==0)
      { 
         fx<-1
       }else{
         fx<-exp(-(((x-thetaGroup[7])/thetaGroup[9])^2)/2)/thetaGroup[9]
       }

      f1<-thetaGroup[8]*thetaGroup[5]^2*exp(-(((y-thetaGroup[1])/thetaGroup[4])^2)/2)/thetaGroup[4]*fx
      f2<-thetaGroup[8]*thetaGroup[6]^2*exp(-(((y-thetaGroup[2])/thetaGroup[4])^2)/2)/thetaGroup[4]*fx
      f3<-thetaGroup[8]*(2*thetaGroup[5]*thetaGroup[6])*exp(-(((y-thetaGroup[3])/thetaGroup[4])^2)/2)/thetaGroup[4]*fx

      #classification, posteriors for inversion genotypes and groups
      classInd<-classInd+cbind(f1,f2,f3)

      groupInd<-cbind(groupInd,f1+f2+f3)

    }
    
    normConst<-sum(pd[-1,-1]*dx*dy)
    npd<-pd[-1,-1]/normConst

    lk<-c(lk,sum(log(rowSums(classInd)/normConst)))

    classInd<-classInd/rowSums(classInd)
    groupInd<-groupInd/rowSums(groupInd)


    #compute convergence
    dd<-sum((theta-thetanew)^2)
    #cat("convergence:",dd,"likelihood:",pd,"\n")

    cat("doing iteration:", i, "tol:", dd, "\n")
    p<-pnew
    q<-qnew
    piGroup<-piGroupNew
    mu1<-mu1new
    mu2<-mu2new
    mu3<-mu3new
    muGroup<-muGroupNew
    sigma<-sigmanew
    sigmaGroup<-sigmaGroupnew

   }
    cat("tolerance of:", dd, "achieved with:",i,"iterations \n")
    res<-list(thetanew=thetanew,theta=theta,likelihood=lk,classInd=classInd,groupInd=groupInd,npd=npd,ysam=ysam[-1],xsam=xsam[-1])
}



EM<-function(y,mu10,mu20,sigma0,p0,q0,x,muGroup0,sigmaGroup0,piGroup0,tol=10^(-5),it=1000)
{ 
  mu1<-mu10
  mu2<-mu20 
  sigma<-sigma0
  p<-p0
  q<-q0
  mu3<-(mu10+mu20)/2

  muGroup<-muGroup0
  sigmaGroup<-sigmaGroup0
  piGroup<-piGroup0

  lk<-c() 
  i<-0
  dd<-Inf
  DD<-c()

  while(i<it & tol<dd)
  {
    i<-i+1
    theta<-list(mu1,mu2,mu3,sigma,p,q,muGroup,piGroup,sigmaGroup)
    N<-NROW(y)

    pnew<-c()
    qnew<-c()
    piGroupNew<-c()
    mu1new<-c()
    mu2new<-c()
    mu3new<-c()
    muGroupNew<-c()
    sigmanew<-list()
    sigmaGroupnew<-c()
    count<-0

    expYY<-function(yy,mm,invss,ss,gg)
    {
       yy1<-yy[,1]-mm[gg,1]
       yy2<-yy[,2]-mm[gg,2]
       ym<-cbind(yy1,yy2)
       rr<-exp(-sapply(1:NROW(ym), function(x) ym[x,]%*%(invss%*%ym[x,]))/2)/sqrt(det(ss))
    }


    #compute partition function given observation i (sum over all groups and genotypes)
    Z<-rep(0,NROW(y))
    for(gg in 1:length(muGroup0))
    {
  
      if(sigmaGroup[gg]==0)
      { 
        fx<-1
      }else{
        fx<-exp(-(((x-muGroup[gg])/sigmaGroup[gg])^2)/2)/sigmaGroup[gg]
      }

      invss<-ginv(sigma[[gg]])
      ss<-sigma[[gg]]
      ExpYY1<-expYY(y,mu1,invss,ss,gg)
      f1Group<-piGroup[gg]*p[gg]^2*ExpYY1*fx

      ExpYY2<-expYY(y,mu2,invss,ss,gg)
      f2Group<-piGroup[gg]*q[gg]^2*ExpYY2*fx

      ExpYY3<-expYY(y,mu3,invss,ss,gg)
      f3Group<-piGroup[gg]*(2*p[gg]*q[gg])*ExpYY3*fx 
      Z<-Z+f1Group+f2Group+f3Group
    }


    for(gg in 1:length(muGroup0))
    {
      count<-count+1


      if(sigmaGroup[gg]==0)
      { 
        fx<-1
      }else{
        fx<-exp(-(((x-muGroup[gg])/sigmaGroup[gg])^2)/2)/sigmaGroup[gg]
      }

      #subpopulation likelihoods
      invss<-ginv(sigma[[gg]])
      ss<-sigma[[gg]]
      ExpYY1<-expYY(y,mu1,invss,ss,gg)
      f1Group<-piGroup[gg]*p[gg]^2*ExpYY1*fx

      ExpYY2<-expYY(y,mu2,invss,ss,gg)
      f2Group<-piGroup[gg]*q[gg]^2*ExpYY2*fx

      ExpYY3<-expYY(y,mu3,invss,ss,gg)
      f3Group<-piGroup[gg]*(2*p[gg]*q[gg])*ExpYY3*fx

      #responsabilities
      w1<-exp(log(f1Group)-log(Z))
      w2<-exp(log(f2Group)-log(Z))
      w3<-exp(log(f3Group)-log(Z))

      #new adimixture frequencies
      pnew<-c(pnew,sum(2*w1+w3)/sum(w1+w2+w3)/2)
      qnew<-c(qnew,sum(2*w2+w3)/sum(w1+w2+w3)/2)

      #new group admixture
      piGroupNew<-c(piGroupNew,sum(w1+w2+w3)/N)  

      #new means
      muanew<-sapply(1:2, function(x)(sum(y[,x]*(w1-w2))*sum(w2-w1)+sum(y[,x]*(w1+w2+w3))*sum(w1+w2))/(sum(w1-w2)*sum(w2-w1)+sum(w1+w2+w3)*sum(w1+w2)))
      mubnew<-sapply(1:2, function(x)sum((y[,x]-muanew[x])*(w1-w2))/sum(w1+w2))

      mu1newg<-muanew+mubnew
      mu2newg<-muanew-mubnew
      mu3newg<-muanew 

      mu1new<-rbind(mu1new,matrix(mu1newg,nrow=1))
      mu2new<-rbind(mu2new,matrix(mu2newg,nrow=1))
      mu3new<-rbind(mu3new,matrix(mu3newg,nrow=1))


      muGroupNewg<-sum(x*(w1+w2+w3))/sum(w1+w2+w3)
      muGroupNew<-c(muGroupNew,muGroupNewg)
        

      #new variance  
      mm1<-sapply(1:NROW(y),function(x) as.vector(w1[x]*t(t((y[x,]-mu1newg)))%*%(y[x,]-mu1newg)))
      mm1<-matrix(rowSums(mm1),ncol=2)    

      mm2<-sapply(1:NROW(y),function(x) as.vector(w2[x]*t(t((y[x,]-mu2newg)))%*%(y[x,]-mu2newg)))
      mm2<-matrix(rowSums(mm2),ncol=2)    

      mm3<-sapply(1:NROW(y),function(x) as.vector(w3[x]*t(t((y[x,]-mu3newg)))%*%(y[x,]-mu3newg)))
      mm3<-matrix(rowSums(mm3),ncol=2) 

      sigmanew[[count]]<-(mm1+mm2+mm3)/sum(w1+w2+w3)
  
      sigmaGroupnew<-c(sigmaGroupnew,sqrt((sum(w1*(x-muGroupNewg)^2)+sum(w2*(x-muGroupNewg)^2)+sum(w3*(x-muGroupNewg)^2))/sum(w1+w2+w3)))

    }

    thetanew<-list(mu1new=mu1new,mu2new=mu2new,mu3new=mu3new,sigmanew=sigmanew,
    pnew=pnew,qnew=qnew,muGroupNew=muGroupNew,piGroupNew=piGroupNew,sigmaGroupnew=sigmaGroupnew)
    
    #compute log-likelihood of the data
    #normalizing constant

    mx1<-max(sapply(1:length(muGroup0), function(g) sigmanew[[g]][1,1]))
    mn1<-min(sapply(1:length(muGroup0), function(g) sigmanew[[g]][1,1]))
    mx2<-max(sapply(1:length(muGroup0), function(g) sigmanew[[g]][2,2]))
    mn2<-min(sapply(1:length(muGroup0), function(g) sigmanew[[g]][2,2]))

    ysam1<-seq(min(y[,1])-mn1,max(y[,1])+mx1,length.out=50)
    ysam2<-seq(min(y[,2])-mn2,max(y[,2])+mx2,length.out=50)
    ysam<-cbind(rep(ysam1,length(ysam2),),rep(ysam2,each=length(ysam1)))

    xsam<-seq(min(x)-max(sigmaGroupnew),max(x)+max(sigmaGroupnew),length.out=50)
   
    dx<-mean(diff(xsam))
    dy<-c(mean(diff(ysam1)),mean(diff(ysam2)))

    #trick for reducing to one dimension
    if(dx<10^(-5))
    {
       dx<-1
       xsam<-c(0,0) 
    }  

    if(sum(dy<10^(-5))>1)
    {
      dy<-c(1,1)
      ysam1<-c(0,0)
      ysam2<-ysam1
      ysam<-cbind(rep(ysam1,length(ysam2),),rep(ysam2,each=length(ysam1)))
    }

    classInd<-matrix(0,nrow=NROW(y),ncol=3)
    pd<-array(0, dim=c(length(ysam1),length(ysam2),length(xsam)))
    groupInd<-c()
    NG<-0
    normConst<-0

    for(gg in 1:length(muGroup0))
    {
      
      #adjusted model  
      for(jj in 1:length(xsam))
      { 

        if(sigmaGroupnew[gg]==0)
        { 
          fx<-1
        }else{
          fx<-exp(-(((xsam[jj]-muGroupNew[gg])/sigmaGroupnew[gg])^2)/2)/sigmaGroupnew[gg]
        }


        invss<-ginv(sigmanew[[gg]])
        ss<-sigmanew[[gg]]
        ExpYY1<-expYY(ysam,mu1new,invss,ss,gg)
        f1Group<-piGroupNew[gg]*pnew[gg]^2*ExpYY1*fx
        
        ExpYY2<-expYY(ysam,mu2new,invss,ss,gg)
        f2Group<-piGroupNew[gg]*qnew[gg]^2*ExpYY2*fx

        ExpYY3<-expYY(ysam,mu3new,invss,ss,gg)
        f3Group<-piGroupNew[gg]*(2*pnew[gg]*qnew[gg])*ExpYY3*fx

        pdy1y2<-matrix(f1Group+f2Group+f3Group, ncol=length(ysam1))
        pd[,,jj]<-pd[,,jj]+pdy1y2
      }
      
      
        #observed likelihood
        if(sigmaGroupnew[gg]==0)
        { 
          fx<-1
        }else{
          fx<-exp(-(((x-muGroupNew[gg])/sigmaGroupnew[gg])^2)/2)/sigmaGroupnew[gg]
        }

        invss<-ginv(sigmanew[[gg]])
        ss<-sigmanew[[gg]]

        ExpYY1<-expYY(y,mu1new,invss,ss,gg)
        f1Group<-piGroupNew[gg]*pnew[gg]^2*ExpYY1*fx

        ExpYY2<-expYY(y,mu2new,invss,ss,gg)
        f2Group<-piGroupNew[gg]*qnew[gg]^2*ExpYY2*fx

        ExpYY3<-expYY(y,mu3new,invss,ss,gg)
        f3Group<-piGroupNew[gg]*(2*pnew[gg]*qnew[gg])*ExpYY3*fx


        #classification, posteriors for inversion genotypes and groups
        classInd<-classInd+cbind(f1Group,f2Group,f3Group)
        groupInd<-cbind(groupInd,f1Group+f2Group+f3Group)
    } 
    
    normConst<-sum(pd[-1,-1,-1]*dx*dy[1]*dy[2])
    npd<-pd[-1,-1,-1]/normConst

    lk<-c(lk,sum(log(rowSums(classInd)/normConst)))

    classInd<-classInd/rowSums(classInd)
    groupInd<-groupInd/rowSums(groupInd)


    #compute convergence
    d1<-sum(((p-pnew))^2)
    d2<-sum(((q-qnew))^2)
    d3<-sum(((piGroup-piGroupNew))^2)
    d4<-sum(((mu1-mu1new))^2)
    d5<-sum(((mu2-mu2new))^2)
    d6<-sum(((mu3-mu3new))^2)
    d7<-sum(((muGroup-muGroupNew))^2)
     
    d8<-0
    for(gg in 1:length(muGroup0))
      d8<-d8+sum(((sigma[[gg]]-sigmanew[[gg]]))^2)

    d9<-sum(((sigmaGroup-sigmaGroupnew))^2)
    dd<-max(sqrt(c(d1,d2,d3,d4,d5,d6,d7,d8,d9)))

    #cat("convergence:",dd,"likelihood:",pd,"\n")

    cat("doing iteration:", i, "tol:", dd, "\n")
    #DD<-rbind(DD,sqrt(c(d1,d2,d3,d4,d5,d6,d7,d8,d9)))
    #write.table(DD,file="Diagnostic.txt",quote=FALSE,row.name=FALSE,col.name=FALSE)
    p<-pnew
    q<-qnew
    piGroup<-piGroupNew
    mu1<-mu1new
    mu2<-mu2new
    mu3<-mu3new
    muGroup<-muGroupNew
    sigma<-sigmanew
    sigmaGroup<-sigmaGroupnew

   }
    cat("tolerance of:", dd, "achieved with:",i,"iterations \n")
    res<-list(thetanew=thetanew,theta=theta,likelihood=lk,classInd=classInd,groupInd=groupInd,npd=npd,ysam1=ysam1[-1],ysam2=ysam2[-1],xsam=xsam[-1])
}



selectSNPsLD<-function(callinv,geno,annot,method,selsubpop=1:nrow(geno),LDthr=0.6)
{

  seloutSNPs<-annot$name!="."
  nms<-annot[seloutSNPs,2]
  annotRed<-annot[seloutSNPs,]

 

  if(method==1)
  { 
    if(!class(callinv)=="invClust")
      stop("method 1 only implemented for callinv of class invClust")

    ginv<-apply(callinv["genotypes"],1,function(x) which(max(x)==x)[1])-1
    genotype1<-ginv[selsubpop]
    

    genoinvMat<-matrix(as.raw(genotype1+1),ncol=1)
    rownames(genoinvMat)<-rownames(geno)[selsubpop]
    colnames(genoinvMat)<-"genotype1"
    inversionGenotypes<-new("SnpMatrix",genoinvMat)

    LD<-ld(inversionGenotypes,geno[selsubpop,seloutSNPs],stats=c("R.squared","R"))
 
    snps1<-annotRed[LD$R.squared[1,]>LDthr & LD$R >0,] 
    snps2<-annotRed[LD$R.squared[1,]>LDthr & LD$R <0,] 
   

    SNPLDnew<-list(snps1[,2],snps2[,2])
    allele<-c(rep("I",length(SNPLDnew[[1]])),rep("N",length(SNPLDnew[[2]]))   )

    LDselSNPs<-c(LD$R.squared[LD$R.squared[1,]>LDthr & LD$R >0],LD$R.squared[LD$R.squared[1,]>LDthr & LD$R <0])

    SNPtagg<-data.frame(rs=unlist(SNPLDnew),allele=allele,LD=round(LDselSNPs,3))
    SNPtagg<-SNPtagg[order(-SNPtagg$LD),]
    #SNPtagg<-SNPtagg[grep("rs",SNPtagg$rs),]
    SNPtaggstr<-paste(unlist(lapply(SNPLDnew,paste,collapse=",")),collapse=";")

    write.table(SNPtaggstr,file="SNPtaggstr.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(file="SNPtagg.txt",SNPtagg,quote=FALSE,row.names=FALSE)
    out<-SNPtaggstr
  }
  
  if(method==2)
  { 
    ginv<-attr(callinv,"genoRef")$invs[selsubpop]
    genotype1<-rep(0,length(ginv))
    genotype1[ginv%in%c(1)]<-2
    genotype1[ginv%in%c(2,4)]<-1

    genotype2<-rep(0,length(ginv))
    genotype2[ginv%in%c(3)]<-2
    genotype2[ginv%in%c(5,2)]<-1

    genotype3<-rep(0,length(ginv))
    genotype3[ginv%in%c(6)]<-2
    genotype3[ginv%in%c(5,4)]<-1


    genoinvMat<-matrix(as.raw(c(genotype1,genotype2,genotype3)+1),ncol=3)
    rownames(genoinvMat)<-rownames(geno)[selsubpop]
    colnames(genoinvMat)<-c("genotype1","genotype2","genotype3")
    inversionGenotypes<-new("SnpMatrix",genoinvMat)
  

    LD<-ld(inversionGenotypes,geno[selsubpop,seloutSNPs],stats="R.squared")
    snps1<-annotRed[LD[1,]>LDthr,]
    snps2<-annotRed[LD[2,]>LDthr,]
    snps3<-annotRed[LD[3,]>LDthr,]

    SNPLDnew<-list(snps1[,2],snps2[,2],snps3[,2])
    allele<-c(rep("N",length(SNPLDnew[[1]])),rep("I",length(SNPLDnew[[2]])) ,rep("B",length(SNPLDnew[[3]]))  )

    LDselSNPs<-c(LD[1,LD[1,]>LDthr],LD[2,LD[2,]>LDthr],LD[3,LD[3,]>LDthr])

    SNPtagg<-data.frame(rs=unlist(SNPLDnew),allele=allele,LD=round(LDselSNPs,3))
    SNPtagg<-SNPtagg[order(-SNPtagg$LD),]
    #SNPtagg<-SNPtagg[grep("rs",SNPtagg$rs),]
    SNPtaggstr<-paste(unlist(lapply(SNPLDnew,paste,collapse=",")),collapse=";")

    write.table(SNPtaggstr,file="SNPtaggstr.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(file="SNPtagg.txt",SNPtagg,quote=FALSE,row.names=FALSE)
    out<-SNPtaggstr
  }
  out
}



