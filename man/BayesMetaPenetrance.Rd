\name{BayesMetaPenetrance}
\alias{BayesMetaPenetrance}
\alias{BayesMetapenetrance}
\title{Estimates consensus age-specific penetrance of cancer for carriers of a pathogenic gene mutation}
\description{
\code{BayesMetaPenetrance} is used to estimate meta-analytic age-specific risk of getting cancer (penetrance) for carriers of pathogenic variants of a specific gene. Information from studies reporting different types of risk measuers (these risk measures include age-specific penetrance, relative risk (RR), standard incidence ratio (SIR), and odds ratio (OR).) on that gene are integrated to give the age-specific penetrance for carriers up to age 85. In this package, SIR is treated in the same manner as RR under the assumption of rare disease.
}
\usage{
BayesMetaPenetrance(penet,RR_studies=TRUE,RR,OR_studies=TRUE,OR,zero_studies=FALSE,zero_OR,
ages=seq(40,80,10),n.iter=30000,n.burn=15000,CrI=FALSE,pl=FALSE,ylim=c(0,1),xlim=c(40,80))
}
\arguments{
  \item{penet}{
  A data frame containing information for studies reporting age-specific penetrance with following columns:
  \tabular{rl}{
  \bold{\emph{study_label}}\tab label to identify the set of penetrance values from a single study.\cr
  \bold{\emph{penetrance}}\tab age-specific penetrance value. \cr
  \bold{\emph{penet_ci_lower}}\tab lower 95\% confidence limit. \cr
  \bold{\emph{penet_ci_upper}}\tab upper 95\% confidence limit. \cr
  \bold{\emph{ages_penet}}\tab ages at which penetrance is reported.
  }
  }
   \item{RR_studies}{
  A logical variable to indicate whether studies reporting RR or SIR are included in the meta analysis. Default is TRUE.
  }
  \item{RR}{
   If \code{RR_studies} is TRUE, a data frame containing information for all RR and SIR studies with following columns where each row corresponds to a single study
  \tabular{rl}{
  \bold{\emph{R.est}}\tab reported RR/SIR value. \cr
  \bold{\emph{RR.ci.lower}}\tab lower 95\% confidence limit. \cr
  \bold{\emph{RR.ci.upper}}\tab upper 95\% confidence limit. \cr
  \bold{\emph{A}}\tab mean age of onset for carriers. Default is 63. \cr
  \bold{\emph{V}}\tab  sd of age of onset for carriers. Default is 14.00726.  \cr
  \bold{\emph{A.lo}}\tab minimum possible age of onset for carriers. Default is 20.   \cr
  \bold{\emph{A.hi}}\tab maximum possible age of onset for carriers. Default is 95.  \cr
  \bold{\emph{A0}}\tab mean age of onset for non-carriers. Default is 63.  \cr
  \bold{\emph{V0}}\tab sd of age of onset for non-carriers. Default is 14.00726.   \cr
  \bold{\emph{A0.lo}}\tab minimum possible age of onset for non-carriers. Default is 20.   \cr
  \bold{\emph{A0.hi}}\tab maximum possible age of onset for non-carriers. Default is 95.
   }
  }
  \item{OR_studies}{
  A logical variable to indicate whether studies reporting OR are included in the meta analysis. Default is TRUE.
  }
   \item{OR}{
   If \code{OR_studies} is TRUE, a data frame containing information for all OR studies with following columns where each row corresponds to a single study
  \tabular{rl}{
  \bold{\emph{OR.est}}\tab reported OR value. \cr
  \bold{\emph{OR.ci.lower}}\tab lower 95\% confidence limit. \cr
  \bold{\emph{OR.ci.upper}}\tab upper 95\% confidence limit. \cr
  \bold{\emph{A}}\tab mean age of onset for cases. Default is 63.  \cr
  \bold{\emph{V}}\tab  sd of age of onset for cases. Default is 14.00726.\cr
  \bold{\emph{A.lo}}\tab minimum possible age of onset among cases. Default is 20.  \cr
  \bold{\emph{A.hi}}\tab maximum possible age of onset among cases. Default is 95.  \cr
  \bold{\emph{A0}}\tab mean age of controls (at inclusion in study). Default is value A.\cr
  \bold{\emph{V0}}\tab sd of  age of controls (at inclusion in study). Default is value V. \cr
  \bold{\emph{A0.lo}}\tab minimum  possible age of controls (at inclusion in study). Default is value A.lo.\cr
  \bold{\emph{A0.hi}}\tab maximum  possible age of controls (at inclusion in study). Default is value A.hi.
   }
  }
  \item{zero_studies}{
  A logical variable to indicate whether information from case control studies where no mutations were detected in controls are included in the meta analysis. These studies should not be included in \code{OR}. Default if FALSE.
  }
  \item{zero_OR}{
   If \code{zero_studies} is TRUE, a data frame containing information for such studies with following columns where each row corresponds to a single study
  \tabular{rl}{
  \bold{\emph{carrier.cases}}\tab number of carrier cases. \cr
  \bold{\emph{non_carrier.cases}}\tab  number of non-carrier cases. \cr
  \bold{\emph{non_carrier.controls}}\tab  number of non-carrier controls. \cr
 \bold{\emph{A.lo}}\tab minimum possible age of onset among cases. Default is 20.  \cr
  \bold{\emph{A.hi}}\tab maximum possible age of onset among cases. Default is 95.  \cr
  \bold{\emph{A0}}\tab mean age of controls (at inclusion in study). Default is value A.\cr
  \bold{\emph{V0}}\tab sd of  age of controls (at inclusion in study). Default is value V. \cr
  \bold{\emph{A0.lo}}\tab minimum  possible age of controls (at inclusion in study). Default is value A.lo.\cr
  \bold{\emph{A0.hi}}\tab maximum  possible age of controls (at inclusion in study). Default is value A.hi.
   }
  }
  \item{ages}{
  Ages at which penetrance values are required in the output. Default is 40 to 80 years at increments of 10.
  }
  \item{n.iter}{
  Number of MCMC iterations to run. Default is 30000.
  }
  \item{n.burn}{
  Number of MCMC iterations to burn-in. Default is 15000.
  }
  \item{CrI}{
  If CrI=TRUE, returns the 95\% credible intervals of the estimated age specific penetrance values. Default is FALSE.
  }
  \item{pl}{
  If pl=TRUE, returns a plot of the estimated age specific penetrance value vs age along with 95\% credible intervals. Default is FALSE.
  }
  \item{ylim}{
  If pl=TRUE, numeric vector of length 2, giving the y coordinates. Default is c(0,1).
  }
  \item{xlim}{
  If pl=TRUE, numeric vector of length 2, giving the x coordinates. Default is c(40,80).
  }
}
\details{
The \code{BayesMetaPenetrance} function estimates age-specific risks of developing cancer associated with mutations
in a specific gene, i.e., penetrance estimation. A meta-analysis approach based on a Bayesian hierarchical random-effects model is used to obtain penetrance estimates integrating
studies reporting different types of risk measures (e.g., penetrance, relative risk, odds
ratio, and standard incidence ratio) while accounting for associated uncertainties.

The cumulative penetrance \eqn{F_s(t|\kappa_s,\lambda_s)} at age \eqn{t} for study \eqn{s} is assumed to be given by the c.d.f. of a Weibull distribution with shape parameter \eqn{\kappa_s} and scale parameter \eqn{\lambda_s}.


The prior distributions are \eqn{\pi(\kappa_s|a, b) =} Gamma\eqn{(a, b), \;  \pi(\lambda_s|c, d) =} Gamma\eqn{(c, d),\;} where \eqn{a} and \eqn{c} are shape parameters and \eqn{b} and \eqn{d} are scale parameters.

Continuous uniform distributions are assumed for all the hyper-parameters. Specifically, \eqn{\pi(a|l_a, u_a) =} U\eqn{(l_a, u_a), \; \pi(b|l_b, u_b) =} U\eqn{(l_b, u_b), \;  \pi(c|l_c, u_c) =} U\eqn{(l_c, u_c), \;} and \eqn{ \; \pi(d|l_d, u_d) =} U\eqn{(l_d, u_d)} with \eqn{l_a, u_a,l_b, u_b,l_c, u_c,l_d, u_d} pre-specified.

Posterior distributions are obtained via a Markov chain Monte Carlo algorithm. The default values for age distributions are based on the distribution of age of onset of breast cancer for US general population (https://seer.cancer.gov/statfacts/html/breast.html). Credible interval at a given age is the \eqn{0.025^{th}} and \eqn{0.0975^{th}} quantiles of the posterior distribution of penetrance estimate at that age.
}
\value{\code{BayesMetaPenetrance} returns a list of two objects
\item{Ages}{Ages corresponding to the estimated penetrance values}
\item{penetrance}{Estimated age-specific penetrance values}
}
\references{
1. Ruberu, T. L. M, Braun, D, Parmigiani, G. and Biswas, S. Bayesian Meta-Analysis of Penetrance for Cancer Risk (2023). \href{https://doi.org/10.48550/arXiv.2304.01912}{arXiv:2304.01912} [STAT.ME].
\cr
2. Marabelli, M., Cheng, S. C., and Parmigiani, G. (2016). Penetrance of ATM Gene Mutations in Breast Cancer: A Meta-Analysis of Different Measures of Risk. Genetic Epidemiology, 40, 425-431.
}
\author{
Thanthirige Lakshika Maduwanthi Ruberu, Swati Biswas
}
\examples{

#penetrance value
set.seed(123)

penet1=c(0.134,0.230,0.400,0.571,0.706) #penetrance values from study 1
penet2=c(0.092,0.146, 0.260,0.403,0.556) #penetrance values from study 2


penet1_ci_low=c( 0.113, 0.203,0.369,0.539,0.675) #lower 95\% confidence limit of values in penet1
penet1_ci_hi=c(0.155, 0.256,0.430,0.601,0.734) #upper 95\% confidence limit of values in penet1

penet2_ci_low=c( 0.066,0.114, 0.221 , 0.359,0.508) #lower 95\% confidence limit of values in penet2
penet2_ci_hi=c(0.117,0.176,0.298,0.446,0.600) #upper 95\% confidence limit of values in penet2



ages_penet1=c(35,40,50,60,70) # ages corresponding to values in penet 1
ages_penet2=c(40,50,60,70,80) # ages corresponding to values in penet 2



study_number=c(rep(1,5),rep(2,5)) # 1 for 5 records of study 1, 2 for study 2.

penet=data.frame('penetrance'=c(penet1,penet2),
              'penet_ci_lower'=c(penet1_ci_low,penet2_ci_low),
                'penet_ci_upper'=c(penet1_ci_hi,penet2_ci_hi),
                'ages_penet'=c(ages_penet1,ages_penet2),
                'study_label'=study_number)


## RR values
R.est=c(4.71) #RR or SIR values reported by each study
RR.ci.lower=c(3.80) #corresponding 95\% lower limit
RR.ci.upper=c(5.84) #corresponding 95\% upper limit


#Age related summaries from each study

A=c(60.21)
V=c(17.34)
A.lo=c(20)
A.hi=c(95)
A0=c(67.61)
V0=c( 14.33)
A0.lo=c(20)
A0.hi=c(95)



RR=data.frame(R.est,RR.ci.lower,RR.ci.upper,A,V,A.lo,A.hi,A0,V0,A0.lo,A0.hi)
colnames(RR)=c("R.est","RR.ci.lower","RR.ci.upper","A","V","A.lo","A.hi","A0","V0","A0.lo","A0.hi")


###OR values
OR.est=c(6.25) #OR value reported by each study
OR.ci.lower=c( 1.38 ) #corresponding 95\% lower limit
OR.ci.upper=c(28.23) #corresponding 95\% lower limit

#Age related summaries from each study
A=c(66.69)
V=c(15.34)
A.lo=c(20)
A.hi=c(95)
A0=c(66.69)
V0=c(15.34)
A0.lo=c(20)
A0.hi=c(95)

OR=data.frame(OR.est,OR.ci.lower,OR.ci.upper,A,V,A.lo,A.hi,A0,V0,A0.lo,A0.hi)
colnames(OR)=c("OR.est","OR.ci.lower","OR.ci.upper","A","V","A.lo","A.hi","A0","V0","A0.lo","A0.hi")


#######studies with no mutations in controls
carrier.cases=c(1)

non_carrier.cases=c(99)

non_carrier.controls=c(100)

#Age related summaries from each study
A=c(48)
V=c(13.57)
A.lo=c(25)
A.hi=c(78)
A0=c(48)
V0=c(13.57)
A0.lo=c(25)
A0.hi=c(78)

zero_OR=data.frame(carrier.cases,non_carrier.cases,non_carrier.controls,
          A,V,A.lo,A.hi,A0,V0,A0.lo,A0.hi)
colnames(zero_OR)=c("carrier.cases","non_carrier.cases","non_carrier.controls",
          "A","V","A.lo","A.hi","A0","V0","A0.lo","A0.hi")


#The following example command uses a very small number of MCMC iterations
#to allow a quick initial check of the code but for actual analysis, it is recommended
#to use the default number of iterations and burn-in.

BayesMetaPenetrance(penet,RR_studies=TRUE,RR=RR,OR_studies=TRUE,OR=OR,zero_studies=TRUE,
zero_OR=zero_OR,ages=seq(40,80,10),n.iter=1000,n.burn=100,CrI=FALSE,pl=FALSE)
}


