library(LURKVecchia)

###### Libraries needed in Vecchia ###########
library(GPvecchia)
library(Matrix)
library(fields)
### SCAD
library(ncvreg)

# For reading in Excel data 
library(readxl)
# Proper scoring rules (log-score and CRPS)
library(scoringRules)

####################################################################################################
####################################################################################################
### SECTION 6
### Simulations from the paper that test the Vecchia approximation convergence of betas ###########
####################################################################################################
####################################################################################################
# ggplot2 and reshape2 are required for the fancy plot
library(ggplot2)
library(reshape2)
##### true parameters
n=500
p=5
beta0=rep(0,p)
theta0=c(9,.3,.5,1) # variance,range,smoothness,nugget


##### simulate data
set.seed(99999)
locs=cbind(runif(n),runif(n)) # random sample on unit square
Sigma.X=exp(-rdist(sample(1:p))/5)
X=mvrnorm(n,rep(0,p),Sigma.X) # correlated predictors (in practice, should include intercept)
Sigma0=theta0[1]*Matern(rdist(locs),range=theta0[2],smoothness=theta0[3])+theta0[4]*diag(n)
epsilon=mvrnorm(1,mu=rep(0,n),Sigma0)
y=X%*%beta0+epsilon


##### convergence for increasing m

# exact transformation
prec.chol=solve(t(chol(Sigma0)))
y.tilde0=prec.chol%*%y
X.tilde0=prec.chol%*%X
beta.hat.exact=solve(crossprod(X.tilde0),crossprod(X.tilde0,y.tilde0))

## vecchia solutions
ms=seq(0,10) #,15,20) #,30,40)
beta.hat.vecchia=matrix(nr=p,nc=length(ms))
for(i.m in 1:length(ms)){
  vecchia.approx=vecchia_specify(locs,ms[i.m])
  transformed.data=transform.iid(cbind(y,X),vecchia.approx,theta0[1:3],theta0[4])
  y.tilde=transformed.data[,1]
  X.tilde=transformed.data[,-1]
  beta.hat.vecchia[,i.m]=as.numeric(solve(crossprod(X.tilde),crossprod(X.tilde,y.tilde)))
}


df <- melt(beta.hat.vecchia)
df.2<- cbind(df,"exact" = rep(beta.hat.exact,11))
df.final <- melt(df.2, id.vars = c("Var1","Var2"))

(p1 <- ggplot(data = df.final,aes(x = Var2,y = value,
                                  linetype = factor(variable),color = factor(Var1)))+
    geom_line(size=1.25)+labs(y = expression(hat(beta)), x = "m")+
    scale_color_viridis_d(option = "D",aesthetics = c("color"))+
    theme(aspect.ratio = 0.75,legend.position = "bottom")+
    scale_x_continuous(breaks = scales::pretty_breaks(5))+
    scale_linetype_manual(values=c("dashed", "solid"))
)
