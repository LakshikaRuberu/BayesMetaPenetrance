\name{bayesmeta}
\alias{bayesmeta}
\alias{BayesMeta}
\title{Estimates consensus age-specific penetrance of cancer for carriers of a gene mutation}
\description{
\code{bayesmeta} is used to estimate meta-analytic age-specific risk of getting cancer (penetrance) by combining information from studies reporting age-specific penetrance, relative risk(RR), standard incidence ratio (SIR), and odds ratio (OR). This function gives the age-specific penetrance at each specified age up to age 85.
}
\usage{
bayesmeta(penet,RR_studies=T,RR,OR_studies=T,OR,ages=20:85,zero_studies=T,zero_OR,pl=F ,ylim=c(0,1),xlim=c(min(ages),max(ages)))
}
\arguments{
  \item{penet}{
  A data frame containing information for studies reporting age-specific penetrance with following columns:
  \tabular{rl}{
  \bold{\emph{study_label}}\tab label identify the set of penetrance values from a single study. E.g, all the age specific penetrance values corresponding to first study should have 1 as the study label \cr
  \bold{\emph{penetrance}}\tab age-specific penetrance value \cr
  \bold{\emph{penet_ci_lower}}\tab lower 95\% confidence limit \cr
  \bold{\emph{penet_ci_upper}}\tab upper 95\% confidence limit \cr
  \bold{\emph{ages_penet}}\tab age at which penetrance is reported
  }
  }
   \item{RR_studies}{
  A logical variable to indicate whether studies reporting RR or SIR are included in the meta analysis. Default is TRUE.
  }
  \item{RR}{
   If \code{RR_studies} is TRUE, a data frame containing information for all RR and SIR studies with following columns where each row corresponding to a single study:
  \tabular{rl}{
  \bold{\emph{R.est}}\tab reported RR/SIR value \cr
  \bold{\emph{RR.ci.lower}}\tab lower 95\% confidence limit \cr
  \bold{\emph{RR.ci.upper}}\tab upper 95\% confidence limit \cr
  \bold{\emph{A}}\tab mean age of onset for carriers. Leave NA for studies that do not report. A deafult value of 63 will be used for NA values and if column is not provided. \cr
  \bold{\emph{V}}\tab  sd of age of onset for carriers. Leave NA for studies that do not report. A deafult value of 14.00726 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A.lo}}\tab lowest age of among for carriers.  Leave NA for studies that do not report. A deafult value of 20 will be used for NA values and if column is not provided. \cr \cr
  \bold{\emph{A.hi}}\tab highest age of among for carriers.  Leave NA for studies that do not report. A deafult value of 95 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0}}\tab mean age of onset for non-carriers. Leave NA for studies that do not report. A deafult value of 63 will be used for NA values and if column is not provided. \cr
  \bold{\emph{V0}}\tab sd of age of onset for non-carriers.  Leave NA for studies that do not report. A deafult value of 20 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0.lo}}\tab lowest age of onset for non-carriers.  Leave NA for studies that do not report. A deafult value of 20 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0.hi}}\tab highest age of onset for non-carriers.  Leave NA for studies that do not report. A deafult value of 95 will be used for NA values and if column is not provided.
   }
  }
  \item{OR_studies}{
  A logical variable to indicate whether studies reporting OR are included in the meta analysis. Default is TRUE.
  }
   \item{OR}{
   If \code{OR_studies} is TRUE, a data frame containing information for all OR studies with following columns with each row corresponding to a single study:
  \tabular{rl}{
  \bold{\emph{OR.est}}\tab reported OR value \cr
  \bold{\emph{OR.ci.lower}}\tab lower 95\% confidence limit \cr
  \bold{\emph{OR.ci.upper}}\tab upper 95\% confidence limit \cr
  \bold{\emph{A}}\tab mean age of onset for cases. Leave NA for studies that do not report. A default value of 63 will be used for NA values and if column is not provided. \cr
  \bold{\emph{V}}\tab  SD of age of onset for cases. Leave NA for studies that do not report. A default value of 14.00726 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A.lo}}\tab Lowest age of onset among cases.  Leave NA for studies that do not report. A default value of 20 will be used for NA values and if column is not provided. \cr \cr
  \bold{\emph{A.hi}}\tab Highest age of onset among cases.  Leave NA for studies that do not report. A default value of 95 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0}}\tab mean age of controls( at study inclusion). Leave NA for studies that do not report. A default value of 63 will be used for NA values and if column is not provided. \cr
  \bold{\emph{V0}}\tab SD of  age of controls( at study inclusion).  Leave NA for studies that do not report. A default value of 20 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0.lo}}\tab Lowest  age of controls( at study inclusion).  Leave NA for studies that do not report. A default value of 20 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0.hi}}\tab Highest  age of controls( at study inclusion).  Leave NA for studies that do not report. A default value of 95 will be used for NA values and if column is not provided.
   }
  }
  \item{zero_studies}{
  A logical variable to indicate whether information from case control studies where no mutations were detected in controls are included in the meta analysis. Default if TRUE.
  }
  \item{zero_OR}{
   If \code{zero_studies} is TRUE, a data frame containing information for such studies with following columns with each row corresponding to a single study:
  \tabular{rl}{
  \bold{\emph{carrier.cases}}\tab number of carrier cases. \cr
  \bold{\emph{non_carrier.cases}}\tab  number of non-carrier cases. \cr
  \bold{\emph{non_carrier.controls}}\tab  number of non-carrier controls. \cr
  \bold{\emph{A}}\tab mean age of onset for cases. Leave NA for studies that do not report. A default value of 63 will be used for NA values and if column is not provided. \cr
  \bold{\emph{V}}\tab  SD of age of onset for cases. Leave NA for studies that do not report. A default value of 14.00726 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A.lo}}\tab Lowest age of onset among cases.  Leave NA for studies that do not report. A default value of 20 will be used for NA values and if column is not provided. \cr \cr
  \bold{\emph{A.hi}}\tab Highest age of onset among cases.  Leave NA for studies that do not report. A default value of 95 will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0}}\tab mean age of controls( at study inclusion). Leave NA for studies that do not report. The corresponding value of A column will be used for NA values and if column is not provided. \cr
  \bold{\emph{V0}}\tab SD of  age of controls( at study inclusion).  Leave NA for studies that do not report. The corresponding value of V column will be used for NA values and if column is not provided. \cr
  \bold{\emph{A0.lo}}\tab Lowest  age of controls( at study inclusion).  Leave NA for studies that do not report. The corresponding value of A.lo column will be used NA values and if column is not provided. \cr
  \bold{\emph{A0.hi}}\tab Highest  age of controls( at study inclusion).  Leave NA for studies that do not report. The corresponding value of A.hi column will be used for NA values and if column is not provided.
   }
  }
  \item{ages}{
  Ages at which penetrance values are required.Default is 20:85 at increments of 1. Maximum possible age is 85
  }
  \item{pl}{
  If pl=TRUE, returns a plot of the estimated age specific penetrance value vs age with confidence intervals. Default is FALSE.
  }
  \item{ylim}{
  If pl=TRUE, numeric vector of length 2, giving the y coordinates. Default is c(0,1).
  }
  \item{xlim}{
  If pl=TRUE, numeric vector of length 2, giving the x coordinates. Default is c(min(\code{ages}),max(\code{ages})).
  }
}
\details{
The \code{bayesmeta} function estimates consensus age-specific risks of developing cancer associated with mutations
in a specific gene, i.e., penetrance estimation. A meta-analysis approach based on a Bayesian hierarchical random-effects model is used to obtain penetrance estimates integrating
studies reporting different types of risk measures (e.g., penetrance, relative risk, odds
ratio, and standard incidence ratio) while accounting for associated uncertainties.

The cumulative penetrance\eqn{F_s(t|\kappa_s,\lambda_s)} at age \eqn{t} for study \eqn{s} is assumed to be given by the c.d.f. of a Weibull distribution with shape parameter \eqn{\kappa_s} and scale parameter \eqn{\lambda_s}.


The prior distributions are \eqn{\pi(\kappa_s|a, b) =} Gamma\eqn{(a, b), \;  \pi(\lambda_s|c, d) =} Gamma\eqn{(c, d),\;} where \eqn{a} and \eqn{c} are shape parameters and \eqn{b} and \eqn{d} are scale parameters.

Continuous uniform distributions are assumed for all the hyper-parameters. Specifically, \eqn{\pi(a|l_a, u_a) =} U\eqn{(l_a, u_a), \; \pi(b|l_b, u_b) =} U\eqn{(l_b, u_b), \;  \pi(c|l_c, u_c) =} U\eqn{(l_c, u_c), \;} and \eqn{ \; \pi(d|l_d, u_d) =} U\eqn{(l_d, u_d)} with \eqn{l_a, u_a,l_a, u_a,l_a, u_a,l_a, u_a} pre-specified.

Posterior distributions are obtained via a Markov chain Monte Carlo algorithm, that employ a Metropolis-Hastings algorithm within Gibbs sampling. The default values for age distributions are based on the distribution of age of onset of breast cancer for US general population (https://seer.cancer.gov/statfacts/html/breast.html). Confidence interval at a given age is the \eqn{0.025^{th}} and \eqn{0.0975^{th}} quantiles of the posterior dstribution of penetrance estimate at that age.
}
\value{\code{bayesmeta} returns a list of three objects--- the estimated penetrance values , confidence Interval(using quantiles of the posterior distribution) and the \eqn{95\%} highest posterior density credible interval
\item{penetrance}{Estimated age-specif penetrance values}
 \item{penetrance_CI}{Confidence interval from the qunatiles}
 \item{penetrance_HPD.CI}{95\%  Highest posterior density credible interval}
}
\references{
1. Marabelli, M., Cheng, S. C., and Parmigiani, G. (2016). Penetrance of ATM Gene Mutations in Breast Cancer: A Meta-Analysis of DifferentMeasures of Risk. Genetic Epidemiology, 40, 425-431.
}
\author{
Thanthirige Lakshika Maduwanthi Ruberu
}
\examples{
#penetrance value

penet1=c(0.016,0.042,0.154,0.309,0.493,0.642)
penet2=c(0.023,0.088,0.113,0.138,0.166)

penet1_ci_low=c(0.005,0.016,0.060,0.126,0.213,0.301)
penet1_ci_hi=c(0.059,0.117,0.364,0.670,0.870,0.979)
penet2_ci_low=c(0.009,0.035,0.055,0.074,0.090)
penet2_ci_hi=c(0.059,0.214,0.227,0.255,0.294)

ages_penet1=c(35,40,50,60,70,80)
ages_penet2=c(40,50,60,70,80)

study_number=c(rep(1,6),rep(2,5))

penet=data.frame('penetrance'=c(penet1,penet2),'penet_ci_lower'=c(penet1_ci_low,penet2_ci_low),
                'penet_ci_upper'=c(penet1_ci_hi,penet2_ci_hi),
                 'ages_penet'=c(ages_penet1,ages_penet2),'study_label'=study_number)
## RR values
R.est=c(2.4,2.37,2.9,2.43,2.67)
RR.ci.lower=c(1.3,1.51,1.9,1.32,0.82)
RR.ci.upper=c(4.3,3.78,4.4,4.09,10.56)



A=c(54.35088,49.90698,52,53.7,63)
V=c(13.23212,9.797507,10.25,15.5,14.00726)
A.lo=c(28,31,38,35,20)
A.hi=c(78,68,79,97,95)
A0=c(57,63,63,63,63)
V0=c(8.618916,14.00726,14.00726,14.00726,14.00726)
A0.lo=c(40,20,20,20,20)
A0.hi=c(79,95,95,95,95)



RR=data.frame(R.est,RR.ci.lower,RR.ci.upper,A,V,A.lo,A.hi,A0,V0,A0.lo,A0.hi)
colnames(RR)=c("R.est","RR.ci.lower","RR.ci.upper","A","V","A.lo","A.hi","A0","V0","A0.lo","A0.hi")

###OR values
OR.est=c(6.9311,3.0897,0.5013,2.10,1.82,1.74,2.10,4.4032)
OR.ci.lower=c(0.8513,0.5035538,0.07008921,1.71,1.46,1.46,1,0.5135)
OR.ci.upper=c(56.4294,18.95828,3.584784,2.57,2.27,2.07,4.1,37.7539)


A=c(63,54,34.6645885,56.3529622,62.1,56.5,55.8,47.54)
V=c(14.00726,10,3.872059,14.00726,11.44,20.0038,12,11.47)
A.lo=c(20,30,20,17,21,18,20,18)
A.hi=c(95,79,40,98,94,98,95,95)
A0=c(63,54,34.6645885,54.2128249,61.22,48.5,72,46.99)
V0=c(14.00726,10,3.872059,14.00726,11.82,21.0016,7.5,12.44)
A0.lo=c(20,30,20,15,21.8,11,60,18)
A0.hi=c(95,79,40,94,94.3,95,95,95)

OR=data.frame(OR.est,OR.ci.lower,OR.ci.upper,A,V,A.lo,A.hi,A0,V0,A0.lo,A0.hi)
colnames(OR)=c("OR.est","OR.ci.lower","OR.ci.upper","A","V","A.lo","A.hi","A0","V0","A0.lo","A0.hi")

#######zero _OR
carrier.cases=c(2,1)
non_carrier.cases=c(41,99)

non_carrier.controls=c(150,100)
A=c(38.14,48)
V=c(5.4,13.57)
A.lo=c(20,25)
A.hi=c(46,78)
A0=c(41,48)
V0=c(5.4,13.57)
A0.lo=c(20,25)
A0.hi=c(46,78)

zero_OR=data.frame(carrier.cases,non_carrier.cases,non_carrier.controls,
          A,V,A.lo,A.hi,A0,V0,A0.lo,A0.hi)
colnames(zero_OR)=c("carrier.cases","non_carrier.cases","non_carrier.controls",
          "A","V","A.lo","A.hi","A0","V0","A0.lo","A0.hi")

set.seed(123)
bayesmeta(penet,RR_studies=TRUE,RR=RR,OR_studies=TRUE,OR=OR,ages=20:85,
        zero_studies=TRUE,zero_OR=zero_OR,pl=T ,ylim=c(0,1),xlim=c(min(ages),max(ages)))

}


