###Charles Doss

###methods for getting information by calculating second/cross moments analytically
### analogously to method for first moments.

getBDinform.PO <- function(partialData, Lhat,Mhat,beta.immig,
                           delta=1e-3, n=1024,r=4, prec.tol=1e-12,
                           prec.fail.stop=TRUE){
  means1 <- all.cond.mean.PO(data=partialData,lambda=Lhat,mu=Mhat, nu=beta.immig*Lhat,
                             delta=delta,n=n,r=r, prec.tol=prec.tol,
                             prec.fail.stop=prec.fail.stop)
  means2 <- all.cond.mean2.PO(data=partialData,lambda=Lhat,mu=Mhat, nu=beta.immig*Lhat, 
                              delta=delta, n=n,r=r, prec.tol=prec.tol,
                             prec.fail.stop=prec.fail.stop)
  ENplus <- means1[1]; ENminus <- means1[2]; EHoldtime <- means1[3];
  ENplusSq <- means2[1]; ENminusSq <- means2[2]; EHoldtimeSq <- means2[3]; 
  ENplusNminus <- means2[4]; ENplusHoldtime <- means2[5]; ENminusHoldtime <- means2[6];
  getBDinform.PO.manual(ENplus=ENplus, ENminus=ENminus, EHoldtime=EHoldtime,
                        ENplusSq=ENplusSq, ENminusSq=ENminusSq,
                        EHoldtimeSq=EHoldtimeSq,
                        ENplusNminus=ENplusNminus,
                        ENplusHoldtime=ENplusHoldtime,
                        ENminusHoldtime=ENminusHoldtime,
                        L=Lhat, m=Mhat, beta.immig= beta.immig, 
                        T=getT(partialData));
}


###############################


