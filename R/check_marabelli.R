marabellimeta(penet ,  RR_studies=T, RR ,OR_studies=T, OR , ages=20:85, zero_studies=T, zero_OR ,pl=T ,ylim=c(0,1),xlim=c(20,85),cex.lab=1,cex.axis=1,lwd=1)
{



  ##Define the total likelihood function for monica method for penetrance
  lik=function(lambda, kappa,S2,S3,S4,ages_penet,y.pm.all,Sigma.pm.all,y.pm.RR,sigma.pm.RR ,A.RR, V.RR, A.lo.RR, A.hi.RR, lambda0, kappa0, A0.RR, V0.RR, A0.lo.RR, A0.hi.RR,y.pm.OR,sigma.pm.OR, A.OR, V.OR, A.lo.OR, A.hi.OR,A0.OR, V0.OR, A0.lo.OR, A0.hi.OR){
    #Define a function for logit transformation
    logit=function(x) {log(x/(1-x))}
    LS2=1
    for(i in 1:S2)
    {
      Arg=c(c(1-exp(-((ages_penet[[i]]-15)/lambda)^kappa)))
      Arg=logit(Arg)
      m=y.pm.all[[i]]
      cmg=Sigma.pm.all[[i]]
      L1=dmvnorm(Arg, m, cmg)
      LS2=LS2*L1

    }

    LS3=1
    for(i in 1:S3)
    {
      L=f.RR(lambda, kappa, y.pm.RR[i],sigma.pm.RR[i] ,A.RR[i], V.RR[i], A.lo.RR[i], A.hi.RR[i], lambda0, kappa0, A0.RR[i], V0.RR[i], A0.lo.RR[i], A0.hi.RR[i])
      LS3=LS3*L

    }

    LS4=1
    for(i in 1:S4)
    {


      L=f.OR(lambda, kappa,y.pm.OR[i],sigma.pm.OR[i], A.OR[i], V.OR[i], A.lo.OR[i], A.hi.OR[i], lambda0, kappa0, A0.OR[i], V0.OR[i], A0.lo.OR[i], A0.hi.OR[i])



      LS4=LS4*L

    }

    final_lik=LS2*LS3*LS4
    return(final_lik)

  }

  #Define a function to calculate the likelihood for papers providing the RR or SIR
  f.RR=function(lambda, kappa, y.pm.RR,sigma.pm.RR ,A.RR, V.RR, A.lo.RR, A.hi.RR, lambda0, kappa0, A0.RR, V0.RR, A0.lo.RR, A0.hi.RR){
    #Use the i.weib function to calculate the integrals in the RR formula (F1=numerator; F2=denominator)

    sigma_m=function(mu,stdev){
      mu=log(mu)
      ss=(exp(mu))^2*stdev

      #Calculate the standard error
      sqrt(ss)
    }






    se.RR=sigma_m(mu=y.pm.RR, stdev=sigma.pm.RR)


    A.RR=A.RR-15
    A.lo.RR= A.lo.RR-15
    A.hi.RR=A.hi.RR-15
    A0.RR=A0.RR-15
    A0.lo.RR=A0.lo.RR-15
    A0.hi.RR=A0.hi.RR-15

    F1=i.weib(lambda, kappa, A.RR, V.RR, A.lo.RR, A.hi.RR)
    F2=i.weib(lambda0, kappa0,A0.RR, V0.RR, A0.lo.RR, A0.hi.RR)
    #Define the likelihood for a paper providing the RR
    if(F1==0 &F2==0){est=0}else{
      est=F1/F2}

    dnorm(est, y.pm.RR, se.RR)
  }

  #Define a function to calculate the likelihood for papers providing the OR
  f.OR=function(lambda, kappa,y.pm.OR,sigma.pm.OR, A.OR, V.OR, A.lo.OR, A.hi.OR, lambda0, kappa0, A0.OR, V0.OR, A0.lo.OR, A0.hi.OR){

    sigma_m=function(mu,stdev){
      #bring back to original scale
      mu=log(mu)
      ss=(exp(mu))^2*stdev
      #Calculate the standard error
      sqrt(ss)
    }



    se.OR=sigma_m(mu=y.pm.OR, stdev=sigma.pm.OR)

    A.OR=A.OR-15
    A.lo.OR= A.lo.OR-15
    A.hi.OR=A.hi.OR-15
    A0.OR=A0.OR-15
    A0.lo.OR=A0.lo.OR-15
    A0.hi.OR=A0.hi.OR-15

    #Define the likelihoods for the numerator of the OR (LN for a/b), and for the denominator of the OR (LD for c/d)
    N1=i.weib(lambda, kappa, A.OR, V.OR, A.lo.OR, A.hi.OR)
    N2=i.weib(lambda0, kappa0, A.OR, V.OR, A.lo.OR, A.hi.OR)
    if(N1==0 &N2==0){LN=0}else{
      LN=N1/N2}

    D1=i.esp(lambda, kappa, A0.OR, V0.OR, A0.lo.OR, A0.hi.OR)
    D2=i.esp(lambda0, kappa0, A0.OR,V0.OR ,A0.lo.OR, A0.hi.OR)
    if(D1==0 &D2==0){LD=0}else{
      LD=D1/D2}

    # #Define the likelihood for a paper providing the OR

    if(LN==0 &LD==0){est=0}else{
      est=LN/LD}
    dnorm(est,y.pm.OR , se.OR)
  }

  ###############################ready penetrance studies#####################
  #Define a function for logit transformation
  logit=function(x) {log(x/(1-x))}

  sigma=function(mu, ICmin, ICmax){
    # pass through the logit scale to make the CIs more symmetric)
    logit=function(x) {log(x/(1-x))}
    logmu=logit(mu)
    logICmin=logit(ICmin)
    logICmax=logit(ICmax)
    m=max(abs(logmu-logICmin), abs(logICmax-logmu))
    #conf=100-((100-95)/2)
    #z_value=qnorm(conf, mean = 0, sd = 1)
    s=(m/1.96)
    return(s)
  }

   penet_sim=function(S1,penet)
   {
     penetall=list() #list to store penetrance values
     cov_mat=list() #list object to store matrices
     ages_penet=list()
     #calculate sigma for each penetrance values and CI

     penet$penet_sigma=mapply(sigma,penet$penetrance,penet$penet_ci_lower,penet$penet_ci_upper)

     n_cov=10000000 #large number of values


     for (i in 1:S1)
     {


       penet_sub=subset(penet,penet$study_label==i) #subset penetrances specific to a study
       penet_sub=penet_sub[order(penet_sub$ages_penet),] #order by age

       #check the monotonocity of the penetrance values and its confidence intervals

       #arranging Penetrance in each study by ages and then check monotonocity within each study

       penet_monotonocity=data.frame(index_penet = all(sign(penet_sub$penetrance -lag(penet_sub$penetrance, default = first(penet_sub$penetrance))) >= 0),
                                     index_ci_lower = all(sign(penet_sub$penet_ci_lower - lag(penet_sub$penet_ci_lower, default = first(penet_sub$penet_ci_lower))) >= 0),
                                     index_ci_upper = all(sign(penet_sub$penet_ci_upper - lag(penet_sub$penet_ci_upper, default = first(penet_sub$penet_ci_upper))) >= 0))


       if(any(penet_monotonocity$index_penet==F|penet_monotonocity$index_ci_lower==F|penet_monotonocity$index_ci_upper==F)==T)
       {
         {stop("Penetrance values and their corresponding CI should be monotonic with age")}
       }

       func1 <- function(x, y) rnorm(n_cov, mean = x, sd = y)
       L_list <- mapply(func1,logit(penet_sub$penetrance), penet_sub$penet_sigma) #for each age generate large number of normal random numbers with corresponding meand and SD



       viable_rows=as.data.frame(get_viable_matrix(L_list)) #rcpp function to check monotonocity of the generated numbers across age
       colnames(viable_rows)=c('status')
       L_list=as.data.frame(L_list)
       mat2=L_list[viable_rows$status==1,]  #extract only monotonic rows



       cmg=cov(mat2)

       diag(cmg)=c(penet_sub$penet_sigma)^2  #variance covariance matrix of the penetrance

       penetall=list.append(penetall,logit(penet_sub$penetrance)) #append penetrance values to the list
       cov_mat=list.append(cov_mat,assign(paste("cmg",i,sep = "_"),cmg)) #append cov matrices of each study to a list
       ages_penet=list.append(ages_penet,penet_sub$ages_penet)
     }
     names(cov_mat)=paste("cmg",1:S1,sep = "_")

     penet_list=list('penet'=penetall,'cov_mat'=cov_mat,'ages_penet'=ages_penet)

     return(penet_list)
   }





   ############################################## Functions for Data preparation Relative Risk studies####################################
   #################################### #################################### #################################### ####################################



   #function to calculate variance of log RR/OR

   sigmalog=function(mu, ICmin, ICmax){
     mu=log(mu)  #log of R est
     ICmin=log(ICmin) #min CI limit
     ICmax=log(ICmax) #max CI limit
     #Choose the bigger interval, to be more conservative
     m=max(abs(mu-ICmin), abs(ICmax-mu))
     #Calculate the var in log scale
     s=(m/1.96)^2

   }

   #calculating variance of the provided RR estimates

   if(RR_studies){ #start if

     RR$RR_sigma=mapply(sigmalog,RR$R.est,RR$RR.ci.lower,RR$RR.ci.upper)




     #check if each age column is presented .If provided replace missing values with standard values A,Ao=63, V,V0=14.00726, A.lo,A0.lo=20, A.hi,A0.hi=95.

     #A-age of onset for carriers
     if ('A' %in% colnames(RR))
     {
       RR$A[which(is.na(RR$A)==T)]=63
     }else{
       message('There is no cloumn named A in RR. Replacing A for all studies by  63')
       RR$A=63
     }

     #A0-age of onset for non-carriers
     if ('A0' %in% colnames(RR))
     {
       RR$A0[which(is.na(RR$A0)==T)]=63
     }else{
       message('There is no cloumn named A0 in RR. Replacing A0 for all studies by  63')
       RR$A0=63
     }


     #V-sd for age of onset for carriers
     if ('V' %in% colnames(RR))
     {
       RR$V[which(is.na(RR$V)==T)]=14.00726
     }else{
       message('There is no cloumn named V in RR. Replacing V for all studies by  14.00726')
       RR$V=14.00726
     }

     #V0-sd for age of onset for non-carriers
     if ('V0' %in% colnames(RR))
     {
       RR$V0[which(is.na(RR$V0)==T)]=14.00726
     }else{
       message('There is no cloumn named V0 in RR. Replacing V0 for all studies by  14.00726')
       RR$V0=14.00726
     }


     #A.lo-lowest age observed for of carrier age of onset
     if ('A.lo' %in% colnames(RR))
     {
       RR$A.lo[which(is.na(RR$A.lo)==T)]=20
     }else{
       message('There is no cloumn named A.lo in RR. Replacing A.lo for all studies by  20')
       RR$A.lo=20
     }

     #A.hi-highest age observed for of carrier age of onset
     if ('A.hi' %in% colnames(RR))
     {
       RR$A.hi[which(is.na(RR$A.hi)==T)]=95
     }else{
       message('There is no cloumn named A.hi in RR. Replacing A.hi for all studies by  95')
       RR$A.hi=95
     }

     #A0.hi-highest age observed for of carrier age of onset
     if ('A0.hi' %in% colnames(RR))
     {
       RR$A0.hi[which(is.na(RR$A0.hi)==T)]=95
     }else{
       message('There is no cloumn named A0.hi in RR. Replacing A0.hi for all studies by  95')
       RR$A0.hi=95
     }

     #A0.hi-highest age observed for of carrier age of onset
     if ('A0.lo' %in% colnames(RR))
     {
       RR$A0.lo[which(is.na(RR$A0.lo)==T)]=20
     }else{
       message('There is no cloumn named A0.lo in RR. Replacing A0.lo for all studies by  20')
       RR$A0.lo=20
     }

   } #end if
   ############################################## Functions for Data preparation Odds Ratio studies####################################
   #################################### #################################### #################################### ####################################
   #calculating variance of the provided OR estimates

   if(OR_studies){#start if

     OR$OR_sigma=mapply(sigmalog,OR$OR.est,OR$OR.ci.lower,OR$OR.ci.upper)




     #check if each age column is presented .If provided replace missing values with standard values A=63, V,=14.00726, A.lo=20, A.hi=95 and age distribution of cases =controls

     #A-age of onset for carriers
     if ('A' %in% colnames(OR))
     {
       OR$A[which(is.na(OR$A)==T)]=63
     }else{
       message('There is no cloumn named A in OR. Replacing A for all studies by  63')
       OR$A=63
     }

     #A0-age of onset for non-caORiers
     if ('A0' %in% colnames(OR))
     {
       OR$A0[which(is.na(OR$A0)==T)]=OR$A[which(is.na(OR$A0)==T)]
     }else{
       message('There is no cloumn named A0 in OR. Replacing A0 for all studies by  A')
       OR$A0=OR$A
     }


     #V-sd for age of onset for caORiers
     if ('V' %in% colnames(OR))
     {
       OR$V[which(is.na(OR$V)==T)]=14.00726
     }else{
       message('There is no cloumn named V in OR. Replacing V for all studies by  14.00726')
       OR$V=14.00726
     }

     #V0-sd for age of onset for non-caORiers
     if ('V0' %in% colnames(OR))
     {
       OR$V0[which(is.na(OR$V0)==T)]=OR$V[which(is.na(OR$V0)==T)]
     }else{
       message('There is no cloumn named V0 in OR. Replacing V0 for all studies by  V')
       OR$V0=OR$V
     }


     #A.lo-lowest age observed for of caORier age of onset
     if ('A.lo' %in% colnames(OR))
     {
       OR$A.lo[which(is.na(OR$A.lo)==T)]=20
     }else{
       message('There is no cloumn named A.lo in OR. Replacing A.lo for all studies by  20')
       OR$A.lo=20
     }

     #A.hi-highest age observed for of caORier age of onset
     if ('A.hi' %in% colnames(OR))
     {
       OR$A.hi[which(is.na(OR$A.hi)==T)]=95
     }else{
       message('There is no cloumn named A.hi in OR. Replacing A.hi for all studies by  95')
       OR$A.hi=95
     }

     #A0.hi-highest age observed for of caORier age of onset
     if ('A0.hi' %in% colnames(OR))
     {
       OR$A0.hi[which(is.na(OR$A0.hi)==T)]=OR$A.hi[which(is.na(OR$A0.hi)==T)]
     }else{
       message('There is no cloumn named A0.hi in OR. Replacing A0.hi for all studies by  A.hi')
       OR$A0.hi=OR$A.hi
     }

     #A0.hi-highest age observed for of caORier age of onset
     if ('A0.lo' %in% colnames(OR))
     {
       OR$A0.lo[which(is.na(OR$A0.lo)==T)]=OR$A.lo[which(is.na(OR$A0.lo)==T)]
     }else{
       message('There is no cloumn named A0.lo in OR. Replacing A0.lo for all studies by  A.lo')
       OR$A0.lo=OR$A.lo
     }

   }#end if


   ## if the studies include zero studies do the same##################################
   if(zero_studies)
   { #start if

     #check if each age column is presented .If provided replace missing values with standard values A=63, V,=14.00726, A.lo=20, A.hi=95 and age distribution of cases =controls


     if ('A' %in% colnames(zero_OR))
     {
       zero_OR$A[which(is.na(zero_OR$A)==T)]=63
     }else{
       message('There is no cloumn named A in zero_OR. Replacing A for all studies by  63')
       zero_OR$A=63
     }


     if ('A0' %in% colnames(zero_OR))
     {
       zero_OR$A0[which(is.na(zero_OR$A0)==T)]=zero_OR$A[which(is.na(zero_OR$A0)==T)]
     }else{
       message('There is no cloumn named A0 in zero_OR. Replacing A0 for all studies by  A')
       zero_OR$A0=zero_OR$A
     }



     if ('V' %in% colnames(zero_OR))
     {
       zero_OR$V[which(is.na(zero_OR$V)==T)]=14.00726
     }else{
       message('There is no cloumn named V in zero_OR. Replacing V for all studies by  14.00726')
       zero_OR$V=14.00726
     }


     if ('V0' %in% colnames(zero_OR))
     {
       zero_OR$V0[which(is.na(zero_OR$V0)==T)]=zero_OR$V[which(is.na(zero_OR$V0)==T)]
     }else{
       message('There is no cloumn named V0 in zero_OR. Replacing V0 for all studies by  V')
       zero_OR$V0=zero_OR$V
     }



     if ('A.lo' %in% colnames(zero_OR))
     {
       zero_OR$A.lo[which(is.na(zero_OR$A.lo)==T)]=20
     }else{
       message('There is no cloumn named A.lo in zero_OR. Replacing A.lo for all studies by  20')
       zero_OR$A.lo=20
     }

     #A.hi-highest age observed fzero_OR of cazero_ORier age of onset
     if ('A.hi' %in% colnames(zero_OR))
     {
       zero_OR$A.hi[which(is.na(zero_OR$A.hi)==T)]=95
     }else{
       message('There is no cloumn named A.hi in zero_OR. Replacing A.hi for all studies by  95')
       zero_OR$A.hi=95
     }

     #A0.hi-highest age observed fzero_OR of cazero_ORier age of onset
     if ('A0.hi' %in% colnames(zero_OR))
     {
       zero_OR$A0.hi[which(is.na(zero_OR$A0.hi)==T)]=zero_OR$A.hi[which(is.na(zero_OR$A0.hi)==T)]
     }else{
       message('There is no cloumn named A0.hi in zero_OR. Replacing A0.hi for all studies by  A.hi')
       zero_OR$A0.hi=95
     }

     #A0.hi-highest age observed fzero_OR of cazero_ORier age of onset
     if ('A0.lo' %in% colnames(zero_OR))
     {
       zero_OR$A0.lo[which(is.na(zero_OR$A0.lo)==T)]=zero_OR$A.lo[which(is.na(zero_OR$A0.lo)==T)]
     }else{
       message('There is no cloumn named A0.lo in zero_OR. Replacing A0.lo for all studies by  A.lo')
       zero_OR$A0.lo=zero_OR$A.lo
     }

     zero_OR$OR.est=((zero_OR$carrier.cases+0.5)/(zero_OR$non_carrier.cases+0.5))/((0.5)/(zero_OR$non_carrier.controls+0.5))
     zero_OR$OR_sigma=1/(zero_OR$carrier.cases+0.5)+1/(zero_OR$non_carrier.cases+0.5)+1/(0.5)+1/(zero_OR$non_carrier.controls+0.5)
   }

   #end if

   ################################################preparing Inputs for the Bayesian method######################################
   S1=length(unique(penet$study_label)) # number of studies reporting age specific penetrance


   if(RR_studies){
     S2=nrow(RR)}else{
       S2=0
     } # number of studies reporting RR/SIR



   if(zero_studies & OR_studies)
   {S3=nrow(OR)+nrow(zero_OR)}else if(zero_studies & !OR_studies)
   {S3=nrow(zero_OR)}else if(!zero_studies & OR_studies)
   {S3=nrow(OR)}else
   {S3=0}

   # number of studies reporting OR

   S=S1+S2+S3 #s- Total number of studies

   #Studies reporting age-specific penetrance
   penet_list=penet_sim(S1,penet)
   y.pm.all=unname(penet_list$penet) # penetrance values as input to the Bayesian method
   Sigma.pm.all=unname(penet_list$cov_mat)
   ages_penet=unname(penet_list$ages_penet)





   #Studies reporting RR

   if(RR_studies){
     y.pm.RR=RR$R.est #reported RR
     sigma.pm.RR=RR$RR_sigma # variance of RR
     A.RR=RR$A #mean age of onset for carriers
     V.RR=RR$V # SD of age of onset for carriers
     A.lo.RR=RR$A.lo #lowest age of onset for carriers
     A.hi.RR=RR$A.hi #highest age of onset for carriers
     A0.RR=RR$A0 # mean age of onset non-carriers
     V0.RR=RR$V0 # SD of age of onset non-carriers
     A0.lo.RR=RR$A0.lo #lowest age of onset for non-carriers
     A0.hi.RR=RR$A0.hi #highest age of onset for non-carriers
   }



   #Studies reporting OR
   # if there are zero studies get both sets if not just consider OR
   if(zero_studies & OR_studies){
     y.pm.OR=c(OR$OR.est,zero_OR$OR.est) #reported OR
     sigma.pm.OR=c(OR$OR_sigma,zero_OR$OR_sigma) # variance of OR
     A.OR=c(OR$A,zero_OR$A) #mean age of onset for cases
     V.OR=c(OR$V,zero_OR$V) # SD of age of onset for cases
     A.lo.OR=c(OR$A.lo,zero_OR$A.lo) #lowest age of onset for cases
     A.hi.OR=c(OR$A.hi,zero_OR$A.hi) #highest age of onset for cases
     A0.OR=c(OR$A0,zero_OR$A0) # mean age at inclusion controls
     V0.OR=c(OR$V0,zero_OR$V0) # SD of age at inclusion controls
     A0.lo.OR=c(OR$A0.lo,zero_OR$A0.lo) #lowest age at inclusion controls
     A0.hi.OR=c(OR$A0.hi,zero_OR$A0.hi) #highest age at inclusion controls


   }else if(!zero_studies & OR_studies){

     y.pm.OR=OR$OR.est #reported OR
     sigma.pm.OR=OR$OR_sigma # variance of OR
     A.OR=OR$A #mean age of onset for cases
     V.OR=OR$V # SD of age of onset for cases
     A.lo.OR=OR$A.lo #lowest age of onset for cases
     A.hi.OR=OR$A.hi #highest age of onset for cases
     A0.OR=OR$A0 # mean age at inclusion controls
     V0.OR=OR$V0 # SD of age at inclusion controls
     A0.lo.OR=OR$A0.lo #lowest age at inclusion controls
     A0.hi.OR=OR$A0.hi #highest age at inclusion controls

   }else if(zero_studies & !OR_studies){

     y.pm.OR=zero_OR$OR.est #reported OR
     sigma.pm.OR=zero_OR$OR_sigma # variance of OR
     A.OR=zero_OR$A #mean age of onset for cases
     V.OR=zero_OR$V # SD of age of onset for cases
     A.lo.OR=zero_OR$A.lo #lowest age of onset for cases
     A.hi.OR=zero_OR$A.hi #highest age of onset for cases
     A0.OR=zero_OR$A0 # mean age at inclusion controls
     V0.OR=zero_OR$V0 # SD of age at inclusion controls
     A0.lo.OR=zero_OR$A0.lo #lowest age at inclusion controls
     A0.hi.OR=zero_OR$A0.hi #highest age at inclusion controls
   }




   #Define the i.weib function to calculate the integrals of the numerator of the OR, a/b
   i.weib=function(lambda, kappa, A, V, A.lo, A.hi){
     f.a = function(a, lambda, kappa, A, V, A.lo, A.hi){
       dweibull(a, kappa, lambda)*dnorm(a, A, V)
     }
     I1=integrate(f.a, A.lo, A.hi, lambda=lambda, kappa=kappa, A=A, V=V)$value
     I1/(pnorm(A.hi, A, V)-pnorm(A.lo, A, V))
   }

   #Define the i.esp function to calculate the integrals of the denominator of the OR, c/d
   i.esp=function(lambda, kappa, A, V, A.lo, A.hi){
     f.a = function(a, lambda, kappa, A, V, A.lo, A.hi){
       exp(-(a/lambda)^kappa)*dnorm(a, A, V)
     }
     I1=integrate(f.a, A.lo, A.hi, lambda=lambda, kappa=kappa, A=A, V=V)$value
     I1/(pnorm(A.hi, A, V)-pnorm(A.lo, A, V))
   }
   #
   #
   kappa0_m=2.25881
   lambda0_m=159.714
   #
   # ##Manual estimation of kappa and lambda max
   # # #
   m=0.05
   kappa=seq(1.5, 8, m)
   l=length(kappa)
   #Create a vector of n lambda_s
   n=300
   lambda=seq(51, 200, length=n)
   #Create the Mx matrix which shows the maximum likelihood associated with the different kappa_s and lambda_s
   Mx=matrix(, l, 3)
   #Create an empty matrix and calculate the maximum likelihood associated with lambda

   M=matrix(, n, 2)
   for (i in 1:l){
     for (j in 1:n){

       M[j,1]=lik(lambda[j], kappa[i],S2,S3,S4,ages_penet,y.pm.all,Sigma.pm.all,y.pm.RR,sigma.pm.RR, A.RR, V.RR, A.lo.RR, A.hi.RR, lambda0_m, kappa0_m, A0.RR, V0.RR, A0.lo.RR, A0.hi.RR,y.pm.OR,sigma.pm.OR, A.OR, V.OR, A.lo.OR, A.hi.OR,A0.OR, V0.OR, A0.lo.OR, A0.hi.OR)
       M[j,2]=lambda[j]
     }
     #Sort the M matrix based on lik values (first column)
     ML=M[order(M[,1]),]
     #Populate Mx matrix with max lik and lambda associated with kappa
     Mx[i, 3]=ML[n,1]
     Mx[i, 2]=ML[n,2]
     Mx[i, 1]=kappa[i]
   }
   # #Sort the Mx matrix based on the lik values (3rd column) (ascending)
   ordMx=Mx[order(Mx[,3]),]

   #Consider the last values, which correspond to kappa and lambda associated with lik max
   kappaM=ordMx[l,1]
   lambdaM=ordMx[l,2]


   kap.est<-kappaM
   lam.est<-lambdaM


   mylist=list('kappa'=kap.est,'lambda'=lam.est)

   return(mylist)

















}
