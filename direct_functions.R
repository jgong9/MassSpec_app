make.formula<-function(block.name,block.instances){
  n<-length(block.name)
  myname<-NULL
  for (j in (1:n)){
    thing.to.add<-paste(block.name[j],block.instances[,j],sep="")
    thing.to.add[block.instances[,j]==0]<-""
    myname<-paste(myname,thing.to.add,sep="")
  }
  return(as.vector(as.character(myname)))
}

findbest<-function(feature.ID,m.z,offset,ppm.thresh,block.name,block.mass,block.max,second.dim=NULL,numiter=3){
  #offset<-1.007
  # m.z are the m/z values
  # offset is the bias correction, such that an anchor chemical is correctly identified 
  # block.name is the vector of names for the blocks. Usually atoms.
  # block.mass is the exact mass of each block
  # block.max is the largest block size considered
  make.vec<-function(n){return(0:n)}
  big.matrix<-data.matrix(expand.grid(lapply(block.max,make.vec)))
  colnames(big.matrix)<-block.name
  masses<-as.vector(big.matrix%*%block.mass)
  # correct bias and get residual for SD estimation
  for (k in (1:numiter)){
   m.z.corrected<-m.z-offset
   #masses.outer<-outer(masses,m.z.corrected,FUN="-")
   massdiff.outer<-abs(log(outer(masses,m.z.corrected,FUN="/")))*1e6
   myindex<-apply(massdiff.outer,2,which.min)
  # first get an estimate of bias and error
   #myindex<-apply(abs(masses.outer),2,which.min)
   block.instances<-big.matrix[myindex,]
   masses.predicted<-masses[myindex]
   ppm.diff<-round(abs(log(masses.predicted/m.z.corrected))*1e6,2)
   which.fail<-which(ppm.diff>ppm.thresh)
   
   if (length(which.fail)>0){bias<-mean(m.z.corrected[-which.fail]-masses.predicted[-which.fail])} else {
     bias<-mean(m.z.corrected-masses.predicted)}# don't include failures in the calc
   offset<-offset+bias
  }
  # print(offset)
  #our.formula<-make.formula(block.name,block.instances)
  #table(our.formula==their.formula)
  residual<-masses.predicted-m.z.corrected
  sigma<-sd(residual)
  # plot(m.z.corrected,residual)
  #prior<-1/nrow(big.matrix)
  #like<-dnorm(masses.outer,0,sigma)
  #posterior<-prior*like  #posterior[is.na(posterior)]<-0
  #posterior<-t(t(posterior)/colSums(posterior))
  #myindex<-apply(posterior,2,which.max)
  block.instances<-big.matrix[myindex,]
  our.formula<-make.formula(block.name,block.instances)
  
  #table(our.formula==their.formula)
  
  #maxposterior<-apply(posterior,2,max)
  masses.predicted<-masses[myindex]
  
  #plot(m.z,residual)
  
  #their.formula<-gsub(" ","",data.original$FormulaCH2)
  #table(our.formula==their.formula)
  
  #finaldf<-data.frame(feature.ID,m.z,m.z.corrected,residual,our.formula,masses.predicted,second.dim)
  # print(length(feature.ID))
  # print(length(m.z))
  # print(length(m.z.corrected))
  # print(length(residual))
  # print(length(our.formula))
  # print(length(masses.predicted))
  # print(length(second.dim))
  # 
  mycolornum<-(colSums(massdiff.outer<ppm.thresh))
  mycolornum[mycolornum>1]<-2
  
  if (is.null(second.dim)){second.dim<-rep(NA,length(feature.ID))}
  is.fail<-rep(FALSE,length(feature.ID))
  if (length(which.fail)>0){is.fail[which.fail]<-TRUE}
  finaldf<-data.frame(feature.ID,m.z,m.z.corrected,residual,our.formula,masses.predicted,second.dim,is.fail,mycolornum,ppm.diff)
  # fig<-plot_ly(data = finaldf, x = ~m.z.corrected, y = ~second.dim,text = ~paste("ID:",feature.ID,
  #                             "<br>predicted formula:",our.formula,
  #                             "<br>obs. corrected mass:",round(m.z.corrected,6),
  #                             "<br>predicted mass:",round(masses.predicted,6),
  #                "<br>ppm diff:",ppm.diff),
  #              color=~mycolornum,colors=c("red","black","light blue"))
  # print(fig)
  return(finaldf)
}


findbest2<-function(feature.ID,m.z,offset,ppm.thresh,block.name,block.mass,block.max,second.dim=NULL,numiter=3){
  #offset<-1.007
  # m.z are the m/z values
  # offset is the bias correction, such that an anchor chemical is correctly identified 
  # block.name is the vector of names for the blocks. Usually atoms.
  # block.mass is the exact mass of each block
  # block.max is the largest block size considered
  make.vec<-function(n){return(0:n)}
  big.matrix<-data.matrix(expand.grid(lapply(block.max,make.vec)))
  colnames(big.matrix)<-block.name
  exact.masses<-as.vector(big.matrix%*%block.mass)
  # correct bias and get residual for SD estimation
  numkeep<-3
  #numiter<-2
  #offset<- 1.0071074 # positive ion mode
  #ppm.thresh<-10
  n<-length(m.z)
  resid<-matrix(NA,n,numkeep)
  myindex<-matrix(NA,n,numkeep)
  resid.matrix<-matrix(NA,n,numkeep)
  ppmdiff.matrix<-matrix(NA,n,numkeep)
  masses.predicted<-rep(NA,n)
  
  for (k in (1:numiter)){
    m.z.corrected<-m.z-offset
    for (i in (1:length(m.z))){
      # print(i)
      tempresid<-m.z.corrected[i]-exact.masses
      myindex[i,]<-order(abs(tempresid))[1:numkeep]
      resid.matrix[i,]<-tempresid[myindex[i,]]
      ppmdiff.matrix[i,]<-abs(log(exact.masses[myindex[i,]]/m.z.corrected[i]))*1e6
    }
    masses.predicted<-exact.masses[myindex[,1]]
    ppm.diff<-round(ppmdiff.matrix,2)
    
    ## Joonho's correction
    which.fail<-which(ppm.diff[,1]>ppm.thresh)
    ##
    
    # print(sum(ppm.diff[,1]<ppm.thresh))
    resid<-resid.matrix[,1]
    bias<-mean(resid[ppm.diff[,1]<ppm.thresh])
    # print(bias)
    offset<-offset+bias
  }
  # print(offset)
  sigma<-sd(resid)
  block.instances<-big.matrix[myindex[,1],]
  our.formula<-make.formula(block.name,block.instances)
  masses.predicted<-exact.masses[myindex]
  mycolornum<-(rowSums(ppm.diff<ppm.thresh))
  mycolornum[mycolornum>1]<-2
  
  if (is.null(second.dim)){second.dim<-rep(NA,length(feature.ID))}
  is.fail<-rep(FALSE,length(feature.ID))
  if (length(which.fail)>0){is.fail[which.fail]<-TRUE}
  finaldf<-data.frame(feature.ID,m.z,m.z.corrected,
                      residual = resid,
                      our.formula,masses.predicted,
                      second.dim,
                      is.fail,
                      mycolornum,
                      ppm.diff=ppm.diff[,1]
                      )
  # fig<-plot_ly(data = finaldf, x = ~m.z.corrected, y = ~second.dim,text = ~paste("ID:",feature.ID,
  #                                                                                "<br>predicted formula:",our.formula,
  #                                                                                "<br>obs. corrected mass:",round(m.z.corrected,6),
  #                                                                                "<br>predicted mass:",round(masses.predicted,6),
  #                                                                                "<br>ppm diff:",ppm.diff),
  #              color=~mycolornum)
  # print(fig)
  return(finaldf)
}

prepare.data<-function(filename){
  mydata<-read.csv(filename)
  feature.ID<-mydata[,1]
  m.z<-mydata[,2]
  if (ncol(mydata)>2){second.dim<-mydata[,3]}else{second.dim<-NULL}
  return(data.frame(feature.ID,m.z,second.dim))
}



