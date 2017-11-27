#
# Script to run the CSR Model on data from the CAS Loss Reserve Database
# Uses Stan for the MCMC run
# by Glenn Meyers
####### Most of the script is in functions - They are called at the end.
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
#
# some input
#
losstype="cpdloss"  #"incloss" if incurred loss or "cpdloss" if paid loss
#
# get packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(parallel)
library(doParallel)
scodeU = "
data{
int<lower=1> len_data;
real logprem[len_data];
real logloss[len_data];
int<lower=1,upper=10> w[len_data];
int<lower=1,upper=10> d[len_data];
}
parameters{
real r_alpha[9];
real r_beta[9];
real logelr;
real <lower=0,upper=100000> a_ig[10];
real gamma;
}
transformed parameters{
real alpha[10];
real beta[10];
real speedup[10];
real sig2[10];
real sig[10];
real mu[len_data];
alpha[1] = 0;
for (i in 2:10) alpha[i] = r_alpha[i-1];
for (i in 1:9) beta[i] = r_beta[i];
beta[10] = 0;
speedup[1] = 1;
for (i in 2:10) speedup[i] = speedup[i-1]*(1-gamma);
sig2[10] = gamma_cdf(1/a_ig[10],1,1);
for (i in 1:9) sig2[10-i] = sig2[11-i]+gamma_cdf(1/a_ig[i],1,1);
for (i in 1:10) sig[i] = sqrt(sig2[i]);
for (i in 1:len_data){
mu[i] = logprem[i]+logelr+alpha[w[i]]+beta[d[i]]*speedup[w[i]];
}
}
model {
r_alpha ~ normal(0,3.162);
r_beta ~ normal(0,3.162);
for (i in 1:10) a_ig[i] ~ inv_gamma(1,1);
logelr ~ normal(-.4,3.162);
gamma ~ normal(0,0.05);
//gamma ~ normal(0,0.0005);
for (i in 1:len_data) logloss[i] ~ normal(mu[i],sig[d[i]]);
}
generated quantities{
vector[len_data] log_lik;
for (i in 1:len_data) log_lik[i] = normal_lpdf(logloss[i]|mu[i],sig[d[i]]);
}
"
#
# function to get Schedule P triangle data given ins group and line of business
#
ins.line.data=function(g.code){
  b=subset(a,a$GRCODE==g.code)
  name=b$GRNAME
  grpcode=b$GRCODE
  w=b$AccidentYear
  d=b$DevelopmentLag
  cum_incloss=b[,6]
  cum_pdloss=b[,7]
  bulk_loss=b[,8]
  dir_premium=b[,9]
  ced_premium=b[,10]
  net_premium=b[,11]
  single=b[,12]
  posted_reserve97=b[,13]
  # get incremental paid losses - assume data is sorted by ay and lag
  inc_pdloss=numeric(0)
  for (i in unique(w)){
    s=(w==i)
    pl=c(0,cum_pdloss[s])
    ndev=length(pl)-1
    il=rep(0,ndev)
    for (j in 1:ndev){
      il[j]=pl[j+1]-pl[j]
    }
    inc_pdloss=c(inc_pdloss,il)
  }
  data.out=data.frame(grpcode,w,d,net_premium,dir_premium,ced_premium,
                      cum_pdloss,cum_incloss,bulk_loss,inc_pdloss,single,posted_reserve97)
  return(data.out)
}
#
# initialization function for scodeU
#
initU=function(chain_id){
  set.seed(123+chain_id)
  list(r_alpha=rnorm(9,0,0.2),r_beta=runif(9),a_ig=runif(10),
       logelr=runif(1,-0.75,-0.5),gamma=rnorm(1,0,0.1))
}
pars.list=c("alpha","beta","gamma","logelr","sig","log_lik")
#
# dummy data for compiling
#
data.dummy=list(len_data = 55,
                logprem  = rep(8,55),
                logloss  = rep(8,55),
                w        = c(1:10,1:9,1:8,1:7,1:6,1:5,1:4,1:3,1,2,1),
                d        = c(rep(1,10),rep(2,9),rep(3,8),rep(4,7),rep(5,6),
                             rep(6,5),rep(7,4),rep(8,3),rep(9,2),10))

#
# compile the univariate model
#
fitU = stan(model_code=scodeU,data=data.dummy,seed=123,init=initU,chains=0)
#
# set up function to run stan model and create output
#
model_function=function(grpcode){
  #
  # read and aggregate the insurer data and
  # set up training and test data frames
  #
  cdata=ins.line.data(grpcode)
  w=cdata$w-1987
  d=cdata$d
  #
  # sort the data in order of d, then w within d
  #
  o1=100*d+w
  o=order(o1)
  w=w[o]
  d=d[o]
  premium=cdata$net_premium[o]
  cpdloss=cdata$cum_pdloss[o]
  cpdloss=pmax(cpdloss,1)
  incloss=cdata$cum_incloss[o]-cdata$bulk_loss[o]
  incloss=pmax(incloss,1)
  adata=data.frame(grpcode,w,d,premium,cpdloss,incloss)
  rdata=subset(adata,(adata$w+adata$d)<12)
  if(losstype=="incloss") rloss=rdata$incloss else rloss=rdata$cpdloss
  if(losstype=="incloss") aloss=adata$incloss else aloss=adata$cpdloss
  #
  #  data for the model
  #
  data.u=list(len_data = length(rdata$w),
              logprem  = log(rdata$premium),
              logloss  = log(rloss),
              w        = rdata$w,
              d        = rdata$d)
  #
  # run the model
  #
  stan_thin=1
  stan_iter=5000
  Rhat_target=1.05
  max_Rhat=2
  while ((max_Rhat > Rhat_target)&(stan_thin<65)){
    fitU1=stan(fit = fitU, data = data.u,init=initU,
               seed = 123,iter=stan_iter,thin=stan_thin,
               chains = 4,pars=pars.list,
               control=list(adapt_delta=.9999),
               refresh=0)
    fitU1_summary=as.matrix(summary(fitU1)$summary)[1:32,c(1,3,10)]
    mrh=subset(fitU1_summary,is.na(fitU1_summary[,3])==F)
    max_Rhat=round(max(mrh[,3]),4)
    print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
    stan_thin=2*stan_thin
    stan_iter=2*stan_iter
  }
  stan_thin=stan_thin/2
  #
  # goodness of fit statistics for comparing models
  #
  loglik1=extract_log_lik(fitU1)
  loo1 <- loo(loglik1)
  #
  # extract information from stan output to process in R
  #
  b=extract(fitU1)
  alpha=b$alpha
  beta=b$beta
  gamma=b$gamma
  logelr=b$logelr
  sig=b$sig
  #
  # simulate outcomes for d=10 using parallel processing
  #
  set.seed(123)
  at.wd10=matrix(rdata$cpdloss[55],length(logelr),10)
  for (w in 2:10){
    at.wd10[,w]=rlnorm(length(logelr),log(premium[w])+logelr+alpha[,w],
                         sig[,10])
  }
  #
  # calculate loss statistics and output to data frame
  #
  Premium=subset(rdata,rdata$d==1)$premium
  ss.wd10=rep(0,10)
  ms.wd10=rep(0,10)
  #
  ms.wd10[1]=mean(at.wd10[,1])
  for (w in 2:10){
    ms.wd10[w]=mean(at.wd10[,w])
    ss.wd10[w]=sd(at.wd10[,w])
  }
  Pred.CSR=rowSums(at.wd10)
  ms.td10=mean(Pred.CSR)
  ss.td10=sd(Pred.CSR)
  CSR.Estimate=round(ms.wd10)
  CSR.SE=round(ss.wd10)
  CSR.CV=round(CSR.SE/CSR.Estimate,4)
  act=sum(subset(aloss,adata$d==10)[1:10])
  pct.CSR=sum(Pred.CSR<=act)/length(Pred.CSR)*100
  #
  # put CSR accident year statistics into a data frame
  #
  CSR.Estimate=c(CSR.Estimate,round(ms.td10))
  CSR.SE=c(CSR.SE,round(ss.td10))
  CSR.CV=c(CSR.CV,round(ss.td10/ms.td10,4))
  Premium=c(Premium,sum(Premium))
  Outcome=subset(aloss,adata$d==10)
  Outcome=c(Outcome,sum(Outcome))
  Group=rep(grpcode,11)
  CSR.Pct=c(rep(NA,10),pct.CSR)
  W=c(1:10,"Total")
  risk=data.frame(Group,W,Premium,CSR.Estimate,CSR.SE,CSR.CV,Outcome,CSR.Pct)
  Group=grpcode
  Premium=Premium[11]
  CSR.Estimate=CSR.Estimate[11]
  CSR.SE=CSR.SE[11]
  Outcome=Outcome[11]
  CSR.Pct=CSR.Pct[11]
  mean_gamma=round(mean(gamma),3)
  elpd_loo=round(loo1$elpd_loo,3)
  se_elpd_loo=round(loo1$se_elpd_loo,3)
  p_loo=round(loo1$p_loo,3)
  se_p_loo=round(loo1$se_p_loo,3)
  looic=round(loo1$looic,3)
  se_looic=round(loo1$se_looic,3)
  SumStats=data.frame(Group,Premium,CSR.Estimate,CSR.SE,Outcome,CSR.Pct,
                      mean_gamma,elpd_loo,se_elpd_loo,p_loo,se_p_loo,
                      looic,se_looic,stan_thin)
  output=list(risk=risk,Predictive_Outcome=Pred.CSR,
              ParmSummary=fitU1_summary,
              logelr=logelr,alpha=alpha,beta=beta,
              gamma=gamma,sig=sig,SumStats=SumStats)
  return(output)
}
#
# Single triangle
#
# insurer.data="The directory and file name for the CAS Loss Reserve Database"

insurer.data="http://www.casact.org/research/reserve_data/comauto_pos.csv"
g=353
#outfilename=paste("CSR CA",g,".csv")
a=read.csv(insurer.data)
co_model=model_function(g)
print(co_model$ParmSummary)
print(" ")
print(co_model$risk)
print("")
print(co_model$SumStats)
par(mfrow=c(1,1))
hist(co_model$Predictive_Outcome,xlab="Simulated Outcomes",
     main="Predictive Distribution of Outcomes",
     sub=paste("Mean =",co_model$SumStats$CSR.Estimate,
               " SE =",co_model$SumStats$CSR.SE))
#write.csv(co_model$risk,file=outfilename)
t1=Sys.time()
print(t1-t0)
#
