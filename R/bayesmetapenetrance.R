
###function to carry out bayes meta analysis
BayesMetaPenetrance<-function( penet ,  RR_studies=T, RR ,OR_studies=T, OR , ages=seq(40,80,10), zero_studies=T, zero_OR ,CrI=F,pl=T ,ylim=c(0,1),xlim=c(20,85))
{

 if(max(ages)>85)
 {
   {stop(paste('Maximum possible age is 85'))}
 }


  #checking column_names of Penet
  colnames_penetrance=c("penetrance","penet_ci_lower","penet_ci_upper","ages_penet","study_label")

  if(any(!(colnames(penet) %in% colnames_penetrance))==T)
  {
    index=which(!(colnames(penet) %in% colnames_penetrance))
    {stop(paste0('In penet, there is a column name mismatch in ', colnames_penetrance[index] ,sep = "\n"  ))}
  }


  if(length(colnames(penet))!=length(colnames_penetrance))
  {
    missing_penet=setdiff(colnames_penetrance,colnames(penet))
    {stop(paste0('In penet, there is no column called ', missing_penet ,sep = "\n"  ))}
  }

  # checking whethaer all penetrance values and their CI are positive values
  if(any(penet$penetrance<0)==T)
  {
    {stop("All the penetrance values must be > 0")}
  }

  if(any(penet$penet_ci_lower<0|penet$penet_ci_upper<0)==T)
  {
    {stop("All the confidence Limits of the penetrance values must be > 0")}
  }

  # checking whether all penetrance values <1 and their CI are <=1
  if(any(penet$penetrance>=1)==T)
  {
    {stop("All the penetrance values must be <1")}
  }

  if(any(penet$penet_ci_lower>=1|penet$penet_ci_upper>=1)==T)
  {
    {stop("All the confidence Limits of the penetrance values must be <= 1")}
  }



  #for each of the penet value and CI the lower CI< penetarance < upper CI
  if(any(mapply(function(x,y,z) x<y & y<z,x=penet$penet_ci_lower,y=penet$penetrance,z=penet$penet_ci_upper)==F)==T)
  {
    {stop("Penetrance values should be in between their confidence Interval")}
  }


  ###################################################################################################################################
  #checks for RR studies if available

  if(RR_studies){ #start if


   if(!('R.est' %in% colnames(RR))==T)
    {
      {stop(paste0(' There is no column named R.est in RR' ))}
    }

  if(!('RR.ci.lower' %in% colnames(RR))==T)
  {
    {stop(paste0('There is no column named RR.ci.lower in RR' ))}
  }

  if(!('RR.ci.upper' %in% colnames(RR))==T)
  {
    {stop(paste0('There is no column named RR.ci.upper in RR' ))}
  }




  if(any(RR$R.est<0|RR$RR.ci.lower<0 |RR$RR.ci.upper<0)==T)
  {
    {stop("All the Relative Risk estimates and their confidence limits must be > 0")}
  }

  if(any(mapply(function(x,y,z) x<y & y<z,x=RR$RR.ci.lower,y=RR$R.est,z=RR$RR.ci.upper)==F)==T)
  {
    {stop("Relative Risk values should be in between their confidence Interval")}
  }

  } #end if

  ###################################################################################################################################
  #checks for OR studies  if available
 if(OR_studies){ #start if
  if(!('OR.est' %in% colnames(OR))==T)
  {
    {stop(paste0(' There is no column named OR.est in OR' ))}
  }

  if(!('OR.ci.lower' %in% colnames(OR))==T)
  {
    {stop(paste0('There is no column named OR.ci.lower in OR' ))}
  }

  if(!('OR.ci.upper' %in% colnames(OR))==T)
  {
    {stop(paste0('There is no column named OR.ci.upper in OR' ))}
  }


  if(any(OR$OR.est<0|OR$OR.ci.lower<0 |OR$OR.ci.upper<0)==T)
  {
    {stop("All the Odds ratio estimates and their confidence limits must be > 0")}
  }


  if(any(mapply(function(x,y,z) x<y & y<z,x=OR$OR.ci.lower,y=OR$OR.est,z=OR$OR.ci.upper)==F)==T)
  {
    {stop("Odds ratio values should be in between their confidence Interval")}
  }


 }# end if
  ###################################################################################################################################
  #checks for Zero .OR studies   if available
  if(zero_studies)
  {#start if
    if(!('carrier.cases' %in% colnames(zero_OR))==T)
    {
      {stop(paste0(' There is no column named carrier.cases in zero_OR' ))}
    }

    if(!('non_carrier.cases' %in% colnames(zero_OR))==T)
    {
      {stop(paste0(' There is no column named non_carrier.cases in zero_OR' ))}
    }

    if(!('non_carrier.controls' %in% colnames(zero_OR))==T)
    {
      {stop(paste0(' There is no column named non_carrier.controls in zero_OR' ))}
    }

    if(any(zero_OR$carrier.cases %% 1 != 0 |zero_OR$non_carrier.cases %% 1 != 0  |zero_OR$non_carrier.controls %% 1 != 0 )==T)
    {
      {stop("All values in columns carrier.cases,non_carrier.cases,non_carrier.controls must be integers")}
    }

    if(any(zero_OR$carrier.cases<=0|zero_OR$non_carrier.cases<=0 |zero_OR$non_carrier.controls<=0)==T)
    {
      {stop("All values in columns carrier.cases,non_carrier.cases,non_carrier.controls must be positive")}
    }

  } # end if
###########################################################functions required for Bayesian Meta analysis##############################################################################
  #log 500 function for proposal distributions
  log500=function(x){logb(x, base=500)}

  ####update kappa of studies reporting age-specific penetramce

  update.kappa.pm<-function(y.pm, sigma.pm, ages_penet, kappa.pm, lambda.pm, a, b) {
    mu.pm= pweibull(ages_penet, kappa.pm, lambda.pm) #compute penetrance at each age at current kappa,lambda value

    mu.pm<-logit(mu.pm) # pass in logit scale

    L.pm<-dmvnorm(y.pm, mu.pm, sigma.pm) # compute multivariate normal likelihood -current kappa


    if(L.pm==0) L.pm=10e-300 # making sure that we dont get a 0 for dmvnorm(y.pm, mu.pm, sigma.pm)
    g.pm<-log(L.pm)+log(dgamma(kappa.pm, shape=a, scale=b)) # likelihood* prior of kappa (target posterior given current kappa in log scale)

    kappa.new<-rgamma(1, shape=kappa.pm^2/abs(log500(kappa.pm+0.01)), scale=abs(log500(kappa.pm+0.01))/kappa.pm) # draw a new kappa value with the proposal distribution




    mu.new= pweibull(ages_penet, kappa.new, lambda.pm) #penetrance at new kappa

    mu.new=logit(mu.new)

    L.new<-dmvnorm(y.pm, mu.new, sigma.pm) # compute multivariate normal likelihood -new kappa

    if(L.new==0) L.new=10e-300
    g.new<-log(L.new)+log(dgamma(kappa.new, shape=a, scale=b)) # likelihood* prior (target posterior given new kappa in log scale)

    accept.prob<-exp(g.new-g.pm+log(dgamma(kappa.pm, shape=kappa.new^2/abs(log500(kappa.new+0.01)), scale=abs(log500(kappa.new+0.01))/kappa.new))-log(dgamma(kappa.new, shape=kappa.pm^2/abs(log500(kappa.pm+0.01)), scale=abs(log500(kappa.pm+0.01))/kappa.pm))) # metropolis step acceptance probability in log scale

    if (!is.na(accept.prob)) {
      if(accept.prob > 1) {
        update = 1        #if acceptance probability >1 accept the drawn new kappa as current kappa
      } else {
        update = rbinom(1, size=1, prob=accept.prob)  #else draw a binomial random variable with acceptance probability . If it is 1 update else do not update
      }
      if (update==1) kappa.pm<-kappa.new
    }


    return(kappa.pm)
  }


  update.lambda.pm<-function(y.pm, sigma.pm, ages_penet, kappa.pm, lambda.pm, c, d) {
    mu.pm= pweibull(ages_penet, kappa.pm, lambda.pm) #compute penetrance at each age at current kappa,lambda value

    mu.pm<-logit(mu.pm) # pass in logit scale

    L.pm<-dmvnorm(y.pm, mu.pm, sigma.pm) # compute multivariate normal likelihood -current lambda


    if(L.pm==0) L.pm=10e-300 # making sure that we dont get a 0 for dmvnorm(y.pm, mu.pm, sigma.pm)

    g.pm<-log(L.pm)+log(dgamma(lambda.pm, shape=c, scale=d)) # likelihood* prior of lambda  (target posterior given current lambda in log scale)

    lambda.new<-rgamma(1, shape=lambda.pm^1.6, scale=1/lambda.pm^0.6) # draw a new lambda value with the proposal distribution

    mu.new= pweibull(ages_penet, kappa.pm, lambda.new)  #penetrance at new lambda

    mu.new=logit(mu.new)

    L.new<-dmvnorm(y.pm, mu.new, sigma.pm) #multivariate nomal likelihood for new lmabda
    if(L.new==0) L.new=10e-300

    g.new<-log(L.new)+log(dgamma(lambda.new, shape=c, scale=d)) # likelihood* prior (target posterior given new lambda in log scale)

    accept.prob<-exp(g.new-g.pm+log(dgamma(lambda.pm, shape=lambda.new^1.6, scale=1/lambda.new^0.6))-log(dgamma(lambda.new, shape=lambda.pm^1.6, scale=1/lambda.pm^0.6)))

    if (!is.na(accept.prob)) {
      if (accept.prob > 1) {
        update = 1
      } else {
        update = rbinom(1, size=1, prob=accept.prob)
      }
      if (update==1) lambda.pm<-lambda.new
    }


    return(lambda.pm)
  }


#######################Functions to update kappa and  lambda of Relative Risk and SIR estimates#######################

  #Define the i.weib function to calculate the integrals in the RR formula (Marabelli et al., 2016)
  i.weib=function(lambda, kappa, A, V, A.lo, A.hi){
    f.a = function(a, lambda, kappa, A, V, A.lo, A.hi){
      dweibull(a, kappa, lambda)*dnorm(a, A, V)
    }
    I1=integrate(f.a, A.lo, A.hi, lambda=lambda, kappa=kappa, A=A, V=V)$value
    I1/(pnorm(A.hi, A, V)-pnorm(A.lo, A, V))
  }


  update.kappa.RR<-function(y.pm.RR, sigma.pm.RR, kappa.pm, lambda.pm, a, b,A, V, A.lo, A.hi, lambda0, kappa0, A0, V0, A0.lo, A0.hi) {


    F1=i.weib(lambda.pm, kappa.pm, A, V, A.lo, A.hi)  # numerator for the RR equations -current kappa,lambda
    F2=i.weib(lambda0, kappa0, A0, V0, A0.lo, A0.hi) # denominator of the RR equation

    mu.pm= log(F1)-log(F2)  #RR equation in log Scale

    L.pm<-dnorm(log(y.pm.RR), mu.pm, sqrt(sigma.pm.RR)) # normal likelihood-current kappa

    if(L.pm==0) L.pm=10e-300 # making sure that we dont get a 0 for dmvnorm(y.pm, mu.pm, sigma.pm)

    g.pm<-log(L.pm)+log(dgamma(kappa.pm, shape=a, scale=b)) #target posterior -current kappa


    kappa.new<-rgamma(1, shape=kappa.pm^2/abs(log500(kappa.pm+2000)), scale=abs(log500(kappa.pm+2000))/kappa.pm) # generate new kappa value-

    F1.new=i.weib(lambda.pm, kappa.new, A, V, A.lo, A.hi) #numerator -new kappa value (Note denominator does not change)
    mu.new= log(F1.new)-log(F2) # RR given new kappa



    L.new<-dnorm(log(y.pm.RR), mu.new, sqrt(sigma.pm.RR)) #likelihood -new kappa
    if(L.new==0) L.new=10e-300

    g.new<-log(L.new)+log(dgamma(kappa.new, shape=a, scale=b)) #posterior -new kappa

    accept.prob<-exp(g.new-g.pm+log(dgamma(kappa.pm, shape=kappa.new^2/abs(log500(kappa.new+2000)), scale=abs(log500(kappa.new+2000))/kappa.new))-log(dgamma(kappa.new, shape=kappa.pm^2/abs(log500(kappa.pm+2000)), scale=abs(log500(kappa.pm+2000))/kappa.pm)))


    if (!is.na(accept.prob)) {
      if(accept.prob > 1) {
        update = 1
      } else {
        update = rbinom(1, size=1, prob=accept.prob)
      }
      if (update==1) kappa.pm<-kappa.new
    }


    return(kappa.pm)
  }
  ##############################################
  update.lambda.RR<-function(y.pm.RR, sigma.pm.RR, kappa.pm, lambda.pm, c, d,A, V, A.lo, A.hi, lambda0, kappa0, A0, V0, A0.lo, A0.hi) {
    F1=i.weib(lambda.pm, kappa.pm, A, V, A.lo, A.hi) # numerator for the RR equations -current kappa,lambda
    F2=i.weib(lambda0, kappa0, A0, V0, A0.lo, A0.hi) # denominator of the RR equation

    mu.pm= log(F1)-log(F2) #RR equation in log Scale

    L.pm<-dnorm(log(y.pm.RR), mu.pm, sqrt(sigma.pm.RR)) # normal likelihood-current lambda
    if(L.pm==0) L.pm=10e-300

    g.pm<-log(L.pm)+log(dgamma(lambda.pm, shape=c, scale=d)) #target posterior -current lambda

    lambda.new<-rgamma(1, shape=lambda.pm^1.1, scale=1/lambda.pm^0.1) # generate new lambda value-

    F1.new=i.weib(lambda.new, kappa.pm, A, V, A.lo, A.hi) #numerator -new lambda value (Note denominator does not change)
    mu.new= log(F1.new)-log(F2) # RR given new lambda


    L.new<-dnorm(log(y.pm.RR), mu.new, sqrt(sigma.pm.RR)) #likelihood -new lambda
    if(L.new==0) L.new=10e-300
    g.new<-log(L.new)+log(dgamma(lambda.new, shape=c, scale=d)) #posterior -new lambda

    accept.prob<-exp(g.new-g.pm+log(dgamma(lambda.pm, shape=lambda.new^1.1, scale=1/lambda.new^0.1))-log(dgamma(lambda.new, shape=lambda.pm^1.1, scale=1/lambda.pm^0.1)))


    if (!is.na(accept.prob)) {
      if (accept.prob > 1) {
        update = 1
      } else {
        update = rbinom(1, size=1, prob=accept.prob)
      }
      if (update==1) lambda.pm<-lambda.new
    }


    return(lambda.pm)
  }

  #######################Functions to update kappa and  lambda of Odds ratio estimates#######################




  #Define the i.weib function to calculate the integrals of the numerator of the OR, a/b  (Marabelli et al., 2016)

  i.weib=function(lambda, kappa, A, V, A.lo, A.hi){
    f.a = function(a, lambda, kappa, A, V, A.lo, A.hi){
      dweibull(a, kappa, lambda)*dnorm(a, A, V)
    }
    I1=integrate(f.a, A.lo, A.hi, lambda=lambda, kappa=kappa, A=A, V=V)$value
    I1/(pnorm(A.hi, A, V)-pnorm(A.lo, A, V))
  }

  #Define the i.esp function to calculate the integrals of the denominator of the OR, c/d (Marabelli et al., 2016)
  i.esp=function(lambda, kappa, A, V, A.lo, A.hi){
    f.a = function(a, lambda, kappa, A, V, A.lo, A.hi){
      exp(-(a/lambda)^kappa)*dnorm(a, A, V)
    }
    I1=integrate(f.a, A.lo, A.hi, lambda=lambda, kappa=kappa, A=A, V=V)$value
    I1/(pnorm(A.hi, A, V)-pnorm(A.lo, A, V))
  }







  update.kappa.OR<-function(y.pm.OR, sigma.pm.OR, kappa.pm, lambda.pm, a, b,A, V, A.lo, A.hi, lambda0, kappa0, A0, V0, A0.lo, A0.hi) {
    #Use the i.weib  and i.esp function to calculate the integrals in the OR formula  (LN for a/b), and for the denominator of the OR (LD for c/d)
    N1=i.weib(lambda.pm, kappa.pm, A, V, A.lo, A.hi)
    N2=i.weib(lambda0, kappa0, A, V, A.lo, A.hi)
    LN=log(N1)-log(N2)
    D1=i.esp(lambda.pm, kappa.pm, A0, V0, A0.lo, A0.hi)
    D2=i.esp(lambda0, kappa0,  A0, V0, A0.lo, A0.hi)
    LD=log(D1)-log(D2)

    mu.pm= LN-LD # OR in log scale   current kappa

    L.pm<-dnorm(log(y.pm.OR), mu.pm, sqrt(sigma.pm.OR)) #likelihood- current kappa
    if(L.pm==0) L.pm=10e-300 # making sure that we dont get a 0 for dmvnorm(y.pm, mu.pm, sigma.pm)

    g.pm<-log(L.pm)+log(dgamma(kappa.pm, shape=a, scale=b)) #posterior -currentkappa

    kappa.new<-rgamma(1, shape=kappa.pm^2/abs(log500(kappa.pm+200000)), scale=abs(log500(kappa.pm+200000))/kappa.pm) #geberate new kappa


    N1.new=i.weib(lambda.pm, kappa.new, A, V, A.lo, A.hi)
    D1.new=i.esp(lambda.pm, kappa.new,A0, V0, A0.lo, A0.hi)
    LN.new=log(N1.new)-log(N2)
    LD.new=log(D1.new)-log(D2)

    mu.new= LN.new-LD.new  #new OR

    L.new<-dnorm(log(y.pm.OR), mu.new, sqrt(sigma.pm.OR)) #likelihood-new kappa
    if(L.new==0) L.new=10e-300

    g.new<-log(L.new)+log(dgamma(kappa.new, shape=a, scale=b))

    accept.prob<-exp(g.new-g.pm+log(dgamma(kappa.pm, shape=kappa.new^2/abs(log500(kappa.new+200000)), scale=abs(log500(kappa.new+200000))/kappa.new))-log(dgamma(kappa.new, shape=kappa.pm^2/abs(log500(kappa.pm+200000)), scale=abs(log500(kappa.pm+200000))/kappa.pm))) #posterior new kappa


    if (!is.na(accept.prob)) {
      if(accept.prob > 1) {
        update = 1
      } else {
        update = rbinom(1, size=1, prob=accept.prob)
      }
      if (update==1) kappa.pm<-kappa.new
    }

    return(kappa.pm)
  }


  ########################################################
  update.lambda.OR<-function(y.pm.OR, sigma.pm.OR, kappa.pm, lambda.pm, c, d,A, V, A.lo, A.hi, lambda0, kappa0, A0, V0, A0.lo, A0.hi) {
    #Use the i.weib  and i.esp function to calculate the integrals in the OR formula  (LN for a/b), and for the denominator of the OR (LD for c/d)
    N1=i.weib(lambda.pm, kappa.pm, A, V, A.lo, A.hi)
    N2=i.weib(lambda0, kappa0, A, V, A.lo, A.hi)
    LN=log(N1)-log(N2)
    D1=i.esp(lambda.pm, kappa.pm, A0, V0, A0.lo, A0.hi)
    D2=i.esp(lambda0, kappa0, A0, V0, A0.lo, A0.hi)
    LD=log(D1)-log(D2)

    mu.pm= LN-LD # OR in log scale   current lambda

    L.pm<-dnorm(log(y.pm.OR), mu.pm, sqrt(sigma.pm.OR))#likelihood- current lambdA
    if(L.pm==0) L.pm=10e-300

    g.pm<-log(L.pm)+log(dgamma(lambda.pm, shape=c, scale=d)) #posterior- current lambda

    lambda.new<-rgamma(1, shape=lambda.pm^0.8, scale=1/lambda.pm^-0.2) #generate new lambda


    N1.new=i.weib(lambda.new, kappa.pm, A, V, A.lo, A.hi)
    D1.new=i.esp(lambda.new, kappa.pm, A0, V0, A0.lo, A0.hi)
    LN.new=log(N1.new)-log(N2)
    LD.new=log(D1.new)-log(D2)

    mu.new= LN.new-LD.new #OR new lambda

    L.new<-dnorm(log(y.pm.OR), mu.new, sqrt(sigma.pm.OR)) #likelihood -new lambda
    if(L.new==0) L.new=10e-300

    g.new<-log(L.new)+log(dgamma(lambda.new, shape=c, scale=d)) #posterior new lambda

    accept.prob<-exp(g.new-g.pm+log(dgamma(lambda.pm, shape=lambda.new^0.8, scale=1/lambda.new^-0.2))-log(dgamma(lambda.new, shape=lambda.pm^0.8, scale=1/lambda.pm^-0.2)))

    if (!is.na(accept.prob)) {
      if (accept.prob > 1) {
        update = 1
      } else {
        update = rbinom(1, size=1, prob=accept.prob)
      }
      if (update==1) lambda.pm<-lambda.new
    }
    return(lambda.pm)
  }




  #functions to update a,b c and  d

  update.a<-function(kappa.all, a, b) {
    g.pm<-0
    for (i in 1:S){
      g.pm<-g.pm+log(dgamma(kappa.all[i], shape=a, scale=b)) #posterior for current a
    }
    #
    a.new<-runif(1, max(la, a-consta), min(ua, a+consta)) #generate new a


    g.new<-0
    for (i in 1:S){
      g.new<-g.new+log(dgamma(kappa.all[i], shape=a.new, scale=b)) #posterior for new a
    }

    accept.prob<-exp(g.new-g.pm+log(dunif(a, max(la, a.new-consta), min(ua, a.new+consta)))-log(dunif(a.new, max(la, a-consta), min(ua, a+consta))))

    if (accept.prob > 1) {
      update = 1
    } else {
      update = rbinom(1, size=1, prob=accept.prob)
    }
    if (update==1) a<-a.new

    return(a)
  }


  ### Function to update b (uses kappa_s from all studies)

  update.b<-function(kappa.all, a, b) {
    g.pm<-0
    for (i in 1:S){
      g.pm<-g.pm+log(dgamma(kappa.all[i], shape=a, scale=b)) #posterior for current b
    }

    b.new<-runif(1, max(lb, b-constb), min(ub, b+constb)) #generate new b


    g.new<-0
    for (i in 1:S){
      g.new<-g.new+log(dgamma(kappa.all[i], shape=a, scale=b.new)) #posterior for new b
    }

    accept.prob<-exp(g.new-g.pm+log(dunif(b, max(lb, b.new-constb), min(ub, b.new+constb)))-log(dunif(b.new, max(lb, b-constb), min(ub, b+constb))))

    if (accept.prob > 1) {
      update = 1
    } else {
      update = rbinom(1, size=1, prob=accept.prob)
    }
    if (update==1) b<-b.new

    return(b)
  }





  ### Function to update c (uses lambda_s from all studies)

  update.c<-function(lambda.all, c, d) {
    g.pm<-0
    for (i in 1:S){
      g.pm<-g.pm+log(dgamma(lambda.all[i], shape=c, scale=d)) #posterior for current c
    }

    c.new<-runif(1, max(lc, c-constc), min(uc, c+constc))  #generate new c


    g.new<-0
    for (i in 1:S){
      g.new<-g.new+log(dgamma(lambda.all[i], shape=c.new, scale=d)) #posterior for new c
    }

    accept.prob<-exp(g.new-g.pm+log(dunif(c, max(lc, c.new-constc), min(uc, c.new+constc)))-log(dunif(c.new, max(lc, c-constc), min(uc, c+constc))))

    if (accept.prob > 1) {
      update = 1
    } else {
      update = rbinom(1, size=1, prob=accept.prob)
    }
    if (update==1) c<-c.new

    return(c)
  }

  ### Function to update d (uses lambda_s from all studies)

  update.d<-function(lambda.all, c, d) {
    g.pm<-0
    for (i in 1:S){
      g.pm<-g.pm+log(dgamma(lambda.all[i], shape=c, scale=d)) #posterior for current d0
    }

    d.new<-runif(1, max(ld, d-constd), min(ud, d+constd)) #generate new d


    g.new<-0
    for (i in 1:S){
      g.new<-g.new+log(dgamma(lambda.all[i], shape=c, scale=d.new))
    }

    accept.prob<-exp(g.new-g.pm+log(dunif(d, max(ld, d.new-constd), min(ud, d.new+constd)))-log(dunif(d.new, max(ld, d-constd), min(ud, d+constd)))) #posterior for new d

    if (accept.prob > 1) {
      update = 1
    } else {
      update = rbinom(1, size=1, prob=accept.prob)
    }
    if (update==1) d<-d.new

    return(d)
  }

  ###function to carry out bayes meta analysis
  main.Bayes<-function(S1, S2, S3, y.pm.all, Sigma.pm.all, ages_penet,y.pm.RR,sigma.pm.RR,y.pm.OR,sigma.pm.OR, A.RR, V.RR, A.lo.RR, A.hi.RR,A.OR, V.OR, A.lo.OR, A.hi.OR,A0.RR, V0.RR, A0.lo.RR, A0.hi.RR,A0.OR, V0.OR, A0.lo.OR, A0.hi.OR,kappa.all, lambda.all, a, b, c, d, n.iter, n.burn)
  {
    #kappa.all -#initial kappa values
    #lambda.all #initial lambda values
    S<-S1+S2+S3 #total number of studies

    kappa.all.iter<-matrix(data=NA, nrow=n.iter, ncol=length(kappa.all)) #matrix to store kappa values
    lambda.all.iter<-matrix(data=NA, nrow=n.iter, ncol=length(lambda.all)) #matrix to store lambda values
    a.all.iter<-numeric(n.iter) #matrix to store a values
    b.all.iter<-numeric(n.iter)#matrix to store b values
    c.all.iter<-numeric(n.iter)#matrix to store c values
    d.all.iter<-numeric(n.iter)#matrix to store d values
    penet<-matrix(nrow=n.iter, ncol=length(ages)) #matrix to store age specific penetrance values values


    j<-1 #starting iterations


    for (j in 1:n.iter){ #start

      for (i in 1:S1){ #start updating kappa for studies reporting age specific penetrance
        kappa.all[i]<-update.kappa.pm(y.pm=y.pm.all[[i]], sigma.pm=Sigma.pm.all[[i]], ages_penet=ages_penet[[i]], kappa.pm=kappa.all[i], lambda.pm=lambda.all[i], a=a, b=b)
        if (kappa.all[i]<0) { break}


        lambda.all[i]<-update.lambda.pm(y.pm=y.pm.all[[i]], sigma.pm=Sigma.pm.all[[i]], ages_penet = ages_penet[[i]], kappa.pm=kappa.all[i], lambda.pm=lambda.all[i], c=c, d=d)
        if (lambda.all[i]<0) { break}


      }#end updating kappa for studies reporting age specific penetrance
      #

    if(S2!=0)
    {

      for (i in (S1+1):(S1+S2)){#start updating kappa for studies reporting RR


        kappa.all[i]<-update.kappa.RR(y.pm.RR=y.pm.RR[i-S1], sigma.pm.RR=sigma.pm.RR[i-S1], kappa.pm=kappa.all[i], lambda.pm=lambda.all[i], a=a, b=b,A=A.RR[i-S1], V=V.RR[i-S1], A.lo=A.lo.RR[i-S1], A.hi=A.hi.RR[i-S1], lambda0=lambda0, kappa0=kappa0, A0=A0.RR[i-S1], V0=V0.RR[i-S1], A0.lo=A0.lo.RR[i-S1], A0.hi=A0.hi.RR[i-S1])

        if (kappa.all[i]<0) { break}


        lambda.all[i]<-update.lambda.RR(y.pm.RR=y.pm.RR[i-S1], sigma.pm.RR=sigma.pm.RR[i-S1], kappa.pm=kappa.all[i], lambda.pm=lambda.all[i], c=c, d=d,A=A.RR[i-S1], V=V.RR[i-S1], A.lo=A.lo.RR[i-S1], A.hi=A.hi.RR[i-S1], lambda0=lambda0, kappa0=kappa0, A0=A0.RR[i-S1], V0=V0.RR[i-S1], A0.lo=A0.lo.RR[i-S1], A0.hi=A0.hi.RR[i-S1])
        if (lambda.all[i]<0) { break}

      }#end updating kappa for studies reporting RR

    }  #end if


   if(S3!=0)
      {

      for (i in (S1+S2+1):(S1+S2+S3)){ #start updating kappa for studies reporting OR

        kappa.all[i]<-update.kappa.OR(y.pm.OR=y.pm.OR[i-(S1+S2)], sigma.pm.OR=sigma.pm.OR[i-(S1+S2)], kappa.pm=kappa.all[i], lambda.pm=lambda.all[i], a=a, b=b,A=A.OR[i-(S1+S2)], V=V.OR[i-(S1+S2)], A.lo=A.lo.OR[i-(S1+S2)], A.hi=A.hi.OR[i-(S1+S2)], lambda0=lambda0, kappa0=kappa0, A0=A0.OR[i-(S1+S2)], V0=V0.OR[i-(S1+S2)], A0.lo=A0.lo.OR[i-(S1+S2)], A0.hi=A0.hi.OR[i-(S1+S2)])

        if (kappa.all[i]<0) { break}


        lambda.all[i]<-update.lambda.OR(y.pm.OR=y.pm.OR[i-(S1+S2)], sigma.pm.OR=sigma.pm.OR[i-(S1+S2)], kappa.pm=kappa.all[i], lambda.pm=lambda.all[i], c=c, d=d,A=A.OR[i-(S1+S2)], V=V.OR[i-(S1+S2)], A.lo=A.lo.OR[i-(S1+S2)], A.hi=A.hi.OR[i-(S1+S2)], lambda0=lambda0, kappa0=kappa0, A0=A0.OR[i-(S1+S2)], V0=V0.OR[i-(S1+S2)], A0.lo=A0.lo.OR[i-(S1+S2)], A0.hi=A0.hi.OR[i-(S1+S2)])

        if (lambda.all[i]<0) { break}


      }#end updating kappa for studies reporting OR

    }    #end if

      if ((kappa.all[i]<0) || (lambda.all[i]<0)) { break} #update a,b,c,d
      a<-update.a(kappa.all, a, b)

      b<-update.b(kappa.all, a, b)

      c<-update.c(lambda.all, c, d)

      d<-update.d(lambda.all, c, d)


      penet[j,]=pweibull(ages,a*b,c*d) # calculate age specific penetrance at provided ages

      kappa.all.iter[j,]<-kappa.all
      lambda.all.iter[j,]<-lambda.all
      a.all.iter[j]<-a
      b.all.iter[j]<-b
      c.all.iter[j]<-c
      d.all.iter[j]<-d

    } #end ietrations




    if (any(kappa.all<0) || any(lambda.all<0)) {return("error")}

    kappa.mean<-numeric(S) #vector to store mean of kappa for each study
    lambda.mean<-numeric(S) #vector to store mean of lambda for each study


    for (i in 1:S){#storing kappa,lambda mean after burn in
      kappa.mean[i]<-mean(kappa.all.iter[(n.burn+1):n.iter,i])
      lambda.mean[i]<-mean(lambda.all.iter[(n.burn+1):n.iter,i])
    }

    a.mean<-mean(a.all.iter[(n.burn+1):n.iter]) ## mean of a,b,c,d  after burn in
    b.mean<-mean(b.all.iter[(n.burn+1):n.iter])
    c.mean<-mean(c.all.iter[(n.burn+1):n.iter])
    d.mean<-mean(d.all.iter[(n.burn+1):n.iter])

    penet<-penet[(n.burn+1):(n.iter),] #consider the penetrance values after burn in


    penet.est<-apply(penet, 2, mean) # take mean at each age.These are the final penetrance values

    penet.ci<-apply(penet,2,quantile, probs=c(0.025, 0.975)) #Take 0.025 and 0.975 quantile as confidence interval limits




    result<-vector("list", length = 3) #preparing the output as a list
    result[[1]]<-ages
    result[[2]]<-penet.est
    result[[3]]<-penet.ci



    names(result)<-c('Ages',"penetrance", "penetrance_CI")

    return(result)


  }





  #############################################################################################################################
  ################ Functions for Data_preparation####################################################

  ###################age_specific_penetrance#####################

  #Define a function for logit transformation
  logit=function(x) {log(x/(1-x))}

  #Define a function to calculate the standard errors given the CIs provided by the papers for penetrance
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



  #function to generate variance covariance matrices of penetrances
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



      viable_rows=as.data.frame(.get_viable_matrix(L_list)) #rcpp function to check monotonocity of the generated numbers across age
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



  #Specification of no of iterations and burn in
  n.iter=30000
  n.burn=15000

  #kappa and lambda of approximated weibull curve for non carriers
  kappa0=3.65
  lambda0=143.2426

  kappa.all<-rep(4, S) #starting points  values for kappa and lambda
  lambda.all<-rep(101, S)


  #parameter starting values and hyper prior parameters.
  a=17.5
  b=0.2
  c=53
  d=1.67
  la=7.5 #this and other limits below are global variable
  ua=27.5
  lb=0.15
  ub=0.25
  lc=43
  uc=63
  ld=1.32
  ud=2.02
  constb=0.04 #parameters for proposal distributions of hyperparameters.
  constc=8
  constd=0.22
  consta=9




  result_correct<-main.Bayes(S1, S2, S3, y.pm.all, Sigma.pm.all, ages_penet,y.pm.RR,sigma.pm.RR,y.pm.OR,sigma.pm.OR, A.RR, V.RR, A.lo.RR, A.hi.RR,A.OR, V.OR, A.lo.OR, A.hi.OR,A0.RR, V0.RR, A0.lo.RR, A0.hi.RR,A0.OR, V0.OR, A0.lo.OR, A0.hi.OR,kappa.all, lambda.all, a, b, c, d, n.iter, n.burn)


if(pl)
{
  penetrance_lower=result_correct$penetrance_CI[1,]
  penetrance_upper=result_correct$penetrance_CI[2,]
  plot(x=ages, y= result_correct$penetrance, ylim=ylim,xlim=xlim, type='l',col='darkgreen',lwd=3,xlab='Age',ylab='Penetrance',cex.lab=1.5, cex.axis=1.2)
  legend('topleft',legend=c("Penetrance Estimate ","95% Credible Interval "),col=c('darkgreen','Orange'),lty=c(1,2),bty='n',cex = 1.2,lwd=c(3,3))
lines(x=ages, y= penetrance_lower,col='Orange',lwd=3,lty=2)
  lines(x=ages, y= penetrance_upper,col='Orange',lwd=3,lty=2)
}
ifelse(CrI==F, return(result_correct[1:2]),return(result_correct))

}


