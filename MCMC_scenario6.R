
for(process in 1:200){

# required libraries
library(truncnorm)
library(data.table)
library(msm)
library(ordinal)
library(rootSolve)

# ----------------------------------------------------------------
# SIMULATION OF DATA FOR ORDINAL BAYESIAN LATENT VARIABLE MODEL
# 17TH APRIL 2018
# ----------------------------------------------------------------



load("seed_data_scenario6.RData")

set.seed(process)
# FIRST EQUATION: DISEASE STATUS IS PNORM(Ui)
#uiprobvec = pnorm(uivec)

#divec = rep(0,numitems)
#for (i in 1:numitems)
#{ divec[i] = rbinom(1,1,uiprobvec[i]) }



# SECOND EQUATION: RATER CLASSIFICATIONS
probvec = rep(0,numcat)
wijvec  = NULL
for (j in 1:numraters)
{
    for (i in 1:numitems)
    {
        probvec=rep(0,numcat)
        for (k in 1:numcat)
        {
            probvec[k] = pnorm(alphavec[k+1]- (ajvec[j] + bjvec[j]*uivec[i])) -
            pnorm(alphavec[k]  - (ajvec[j] + bjvec[j]*uivec[i]))
        }
        value = rmultinom(1, size = 1, prob=probvec)
        for (q in 1:numcat)
        { if (value[q,1] == 1) { wij = q } }
        wijvec = c(wijvec,wij)
        
    }
}


cattable=table(wijvec)
round(prop.table(cattable)*50,0)

fullrater = rep(seq(1,numraters),numitems)
fullrater = sort(fullrater)
fullitem  = rep(seq(1,numitems),numraters)
fulldivec = rep(divec,numraters)


# Simulated Binary Dataset contains rater, item and classification and true disease status

simdataset=data.frame(fullrater,fullitem,wijvec,fulldivec)
colnames(simdataset) = c("rater","item","wij","di")


fit<-clmm(as.factor(wij)~ di+(1|rater)+(1|item), link = "probit", Hess=TRUE,
control = clmm.control(maxIter = 100,maxLineIter = 100),threshold = "flexible", gradTol=1e-4, data=simdataset)












K <- length(unique(simdataset$wij)) # <- Num. of categoris
m <- length(unique(simdataset$item)) # <- Num. of cases
n <- length(unique(simdataset$rater)) # <- Num. of raters

D <- simdataset$di[1:m]


w <- matrix(ncol=n,nrow=m)
ii <- 0
for(i in unique(simdataset$item)){
    jj <- 0
    ii <- ii + 1
    for(j in unique(simdataset$rater)){
        jj <- jj + 1
        w[ii,jj] <- simdataset$wij[which(simdataset$item==i & simdataset$rater==j)]
    }
}

# prior settings
gamma <- 0.1
tau <- 0.1

# Initial Values
A <- c(fit$alpha[1],fit$alpha[2],fit$alpha[3],1000)
A[K] <- Inf # set A[K] <- Inf

z0 <- rep(0, m)
z <- matrix(1, nrow=m, ncol=n)
u <- rnorm(m, 0, 1)

a <- rnorm(n, 0, 1)
#a[1] <- 0 # reference rater
b <- rexp(n, rate=1)


mu_a <- 0
tau_a <- 1
mu_b <- 0
tau_b <- 1
#mu_beta <- 0
#tau_beta <- 1
mu_u <- -0.98
tau_u <- 1
#beta <- rnorm(1,0,1)

delta= 10

##################################
## Main MCMC Function ############
##################################


MCMC <- function(m, n, z0, z, u, w, D, a, b, mu_a, tau_a, mu_b, tau_b, A, gamma = 0.1, tau = 0.1, mu_u, tau_u){
    
    # sample z_i0
    
    z0 <- (D)*rtnorm(m, lower=0, upper=Inf, mean=u, sd=1)+(1-D)*rtnorm(m, lower=-Inf, upper=0, mean=u, sd=1)
    
    # sample z_ij (a_1 <- 0)
    
    W1 <- which(w==1, arr.ind=T)
    z[W1] <- rtruncnorm(dim(W1)[1], a=-Inf, b= A[1], mean=a[W1[,2]]+b[W1[,2]]*u[W1[,1]], sd=1)
    W2 <- which(w==2, arr.ind=T)
    z[W2] <- rtruncnorm(dim(W2)[1], a=A[1], b= A[2], mean=a[W2[,2]]+b[W2[,2]]*u[W2[,1]], sd=1)
    W3 <- which(w==3, arr.ind=T)
    z[W3] <- rtruncnorm(dim(W3)[1], a=A[2], b= A[3], mean=a[W3[,2]]+b[W3[,2]]*u[W3[,1]], sd=1)
    W4 <- which(w==4, arr.ind=T)
    z[W4] <- rtruncnorm(dim(W4)[1], a=A[3], b= A[4], mean=a[W4[,2]]+b[W4[,2]]*u[W4[,1]], sd=1)
    
    
    # sample u_i (!! include z_0 !!)
    
    u <- rnorm(m,((z-matrix(a, nrow=m, ncol=n, byrow=TRUE))%*%b+z0+tau_u*mu_u)/(sum(b^2)+1+tau_u), sqrt(1/(sum(b^2)+1+tau_u)))
    model <- function(c){
        F1 <- mean(pnorm(u+c)) -  mean(D)
        c(F1 = F1)
    }
    ss <- multiroot(f = model, start = c(1))
    u <- u + ss$root
    
    
    # mu_u ~ N(0, 1/omega) : omega=1
    # mu_u <- rnorm(1, tau_u*sum(u)/(m*tau_u+1), sqrt(1/(m*tau_u+1)))
    mu_u <- 0
   #  tau_u <- rgamma(1, (0.1+m/2), (0.1+0.5*sum( (u-mu_u)^2  ))  )
     tau_u <- 1
    
    # sample a_j (a_1 : reference group)
    
    #    for(j in 2:n){
    #        a[j] <- rnorm(1,(sum(z[,j]-b[j]*u)+mu_a*tau_a)/(m+tau_a), sqrt(1/(m+tau_a)))
    #    }
#    a[4] <- 0
    for(j in 1:n){
        a[j] <- rnorm(1,(sum(z[,j]-b[j]*u)+mu_a*tau_a)/(m+tau_a), sqrt(1/(m+tau_a)))
    }
#    a <- a-mean(a)
    
    # if Num. of raters (n) is large, then use the following lines
    # a <- rnorm(n,(colSums(z-u%*%t(b))+mu_a*tau_a)/(m+tau_a), sqrt(1/(m+tau_a)))
    # a[1] <- 0
    
    # sample b_j
#    b[4] <- 1
    for(j in 1:n){
        b[j] <- rnorm(1,(sum((z[,j]-a[j])*u)+mu_b*tau_b)/(sum(u^2)+tau_b), sqrt(1/(sum(u^2)+tau_b)))
    }
    
    # beta ~ N(0, 1/0.1)
    # beta <- rnorm(1,(sum(D%*%(z-matrix(a, nrow=m, ncol=n, byrow=T)-matrix(u, nrow=m, ncol=1)%*%matrix(b, nrow=1, ncol=n)))+mu_beta*tau_beta)/(n*sum(D^2)+tau_beta), sqrt(1/(n*sum(D^2)+tau_beta)))
    
    
    # sample mu_a
    
    #mu_a <- rnorm(1, tau_a*sum(a)/(n*tau_a+tau) , sqrt(1/(n*tau_a+tau)))
    mu_a <- 0
    
    # sample tau_a
    
    tau_a <- rgamma(1, gamma+n/2, gamma+0.5*sum((a-mu_a)^2))
    
    # sample mu_b
    
    mu_b <- rnorm(1, tau_b*sum(b)/(n*tau_b+tau) , sqrt(1/(n*tau_b+tau)))
    
    # sample tau_b
    
    tau_b <- rgamma(1, gamma+n/2, gamma+0.5*sum((b-mu_b)^2))
    
    # sample mu_beta
    
#    mu_beta <- rnorm(1, tau_beta*beta/(tau_beta+tau) , sqrt(1/(tau_beta+tau)))
    
    # sample tau_beta
    
#    tau_beta <- rgamma(1, gamma+1/2, gamma+0.5*sum((beta-mu_beta)^2))
    
    
    # sample A (alpha)
    U <- NULL
    L <- NULL
    # sample A1
    L[1] <- max(-delta,A[2]-delta,max(z[which(w==1, arr.ind=T)]))
    U[1] <- min(delta,A[2],min(z[which(w==2, arr.ind=T)]))
    A[1] <- runif(1, L[1], U[1])
    # sample A2 tp Ak-2
    for(k in 2:(K-2)){
        U[k] <- min(A[k-1]+delta,A[k+1],min(z[which(w==k+1, arr.ind=T)]))
        L[k] <- max(A[k-1],A[k+1]-delta,max(z[which(w==k, arr.ind=T)]))
        A[k] <- runif(1, L[k], U[k])
    }
    U[K-1] <- min(A[K-2]+delta,A[K],min(z[which(w==K, arr.ind=T)]))
    L[K-1] <- max(A[K-2],max(z[which(w==K-1, arr.ind=T)]))
    A[K-1] <- runif(1, L[K-1], U[K-1])
    
 
        
    return(list(z0=z0,z=z,
    u=u,
    w=w,
    D=D,
    a=a,
    b=b,
    mu_a=mu_a,
    tau_a=tau_a,
    mu_b=mu_b,
    tau_b=tau_b,
#    beta=beta,
    mu_u=mu_u,
    tau_u=tau_u,
#    mu_beta=mu_beta,
#    tau_beta=tau_beta,
    A=A
#    good1_1 = good1_1,
#    good2_1 = good2_1,
#    good3_1 = good3_1,
#    good4_1 = good4_1
    ))
}


#################################
update <- list()
update[[1]] <- list(z0=z0,z=z,
u=u,
w=w,
D=D,
a=a,
b=b,
mu_a=mu_a,
tau_a=tau_a,
mu_b=mu_b,
tau_b=tau_b,
#mu_beta=mu_beta,
#tau_beta=tau_beta,
#beta=beta,
mu_u=mu_u,
tau_u=tau_u,
A=A)

for(t in 2:30000){
    update[[t]] <- MCMC(m=m,
    n=n,
    z0=update[[t-1]]$z0,
    z=update[[t-1]]$z,
    u=update[[t-1]]$u,
    w=update[[t-1]]$w,
    D=update[[t-1]]$D,
    a=update[[t-1]]$a,
    b=update[[t-1]]$b,
    mu_a=update[[t-1]]$mu_a,
    tau_a=update[[t-1]]$tau_a,
    mu_b=update[[t-1]]$mu_b,
    tau_b=update[[t-1]]$tau_b,
#    mu_beta=update[[t-1]]$mu_beta,
#    tau_beta=update[[t-1]]$tau_beta,
    A=update[[t-1]]$A,
#    beta=update[[t-1]]$beta,
    mu_u=update[[t-1]]$mu_u,
    tau_u=update[[t-1]]$tau_u)
    print(t);
}

out <- list()
for(ind in 1:3000){
    IND <- ind*5+15000
    out[[ind]] <- list(z0=update[[IND]]$z0,z=update[[IND]]$z,
    a=update[[IND]]$a,
    b=update[[IND]]$b,
    u=update[[IND]]$u,
    A=update[[IND]]$A,
    mu_a=update[[IND]]$mu_a,
    mu_b=update[[IND]]$mu_b,
    tau_a=update[[IND]]$tau_a,
    tau_b=update[[IND]]$tau_b,
#    beta=update[[IND]]$beta,
mu_a=update[[IND]]$mu_a,
    mu_b=update[[IND]]$mu_b,
#    mu_beta=update[[IND]]$mu_beta,
    tau_a=update[[IND]]$tau_a,
    tau_b=update[[IND]]$tau_b
 #   tau_beta=update[[IND]]$tau_beta,

#    good1_1=update[[IND]]$good1_1,
#    good2_1=update[[IND]]$good2_1,
#    good3_1=update[[IND]]$good3_1,
#    good4_1=update[[IND]]$good4_1
    )
}


save(out, file=paste0("/Volumes/Seagate\ Backup\ Plus\ Drive/ordinal/s6/",process,".RData"))
#save(out, file=paste0("/Volumes/SONY_32X/ordinal",process,".RData"))
#save(out, file=paste0("~/Desktop/Kerrie_3/",process,".RData"))
remove(list=c("update","out"))

}






sp1 <- sp2 <- sp3 <- se1 <- se2 <- se3 <- matrix(nrow=3000, ncol=50)
sp1.mean <- sp2.mean <- sp3.mean <- se1.mean <- se2.mean <- se3.mean <- matrix(nrow=3000, ncol=1)

for(i in 1:3000){
    sp1[i,] <- rowMeans(1-out[[i]]$sp1)
    sp2[i,] <- rowMeans(1-out[[i]]$sp2)
    sp3[i,] <- rowMeans(1-out[[i]]$sp3)
    se1[i,] <- rowMeans(1-out[[i]]$se1)
    se2[i,] <- rowMeans(1-out[[i]]$se2)
    se3[i,] <- rowMeans(1-out[[i]]$se3)
    sp1.mean[i] <- (1-out[[i]]$sp1.mean)
    sp2.mean[i] <- (1-out[[i]]$sp2.mean)
    sp3.mean[i] <- (1-out[[i]]$sp3.mean)
    se1.mean[i] <- (1-out[[i]]$se1.mean)
    se2.mean[i] <- (1-out[[i]]$se2.mean)
    se3.mean[i] <- (1-out[[i]]$se3.mean)
}

pdf("ROC1.pdf", width=6, height=6)

t<-39
plot(colMeans(sp1)[t],colMeans(se1)[t], xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate", ylab="True Positive Rate",cex=1.5)
points(colMeans(sp2)[t],colMeans(se2)[t],cex=1.5)
points(colMeans(sp3)[t],colMeans(se3)[t],cex=1.5)
lines(c(1,colMeans(sp1)[t]),c(1,colMeans(se1)[t]), lwd=1.5)
lines(c(colMeans(sp1)[t],colMeans(sp2)[t]),c(colMeans(se1)[t],colMeans(se2)[t]), lwd=1.5)
lines(c(colMeans(sp2)[t],colMeans(sp3)[t]),c(colMeans(se2)[t],colMeans(se3)[t]), lwd=1.5)
lines(c(colMeans(sp3)[t],0),c(colMeans(se3)[t],0), lwd=1.5)

t<-18
points(colMeans(sp1)[t],colMeans(se1)[t],pch=t,col=t,cex=1.5)
points(colMeans(sp2)[t],colMeans(se2)[t],pch=t,col=t,cex=1.5)
points(colMeans(sp3)[t],colMeans(se3)[t],pch=t,col=t,cex=1.5)
lines(c(1,colMeans(sp1)[t]),c(1,colMeans(se1)[t]),col=t, lwd=1.5)
lines(c(colMeans(sp1)[t],colMeans(sp2)[t]),c(colMeans(se1)[t],colMeans(se2)[t]),col=t, lwd=1.5)
lines(c(colMeans(sp2)[t],colMeans(sp3)[t]),c(colMeans(se2)[t],colMeans(se3)[t]),col=t, lwd=1.5)
lines(c(colMeans(sp3)[t],0),c(colMeans(se3)[t],0),col=t, lwd=1.5)

t<-4
points(colMeans(sp1)[t],colMeans(se1)[t],pch=t,col=t,cex=1.5)
points(colMeans(sp2)[t],colMeans(se2)[t],pch=t,col=t,cex=1.5)
points(colMeans(sp3)[t],colMeans(se3)[t],pch=t,col=t,cex=1.5)
lines(c(1,colMeans(sp1)[t]),c(1,colMeans(se1)[t]),col=t, lwd=1.5)
lines(c(colMeans(sp1)[t],colMeans(sp2)[t]),c(colMeans(se1)[t],colMeans(se2)[t]),col=t, lwd=1.5)
lines(c(colMeans(sp2)[t],colMeans(sp3)[t]),c(colMeans(se2)[t],colMeans(se3)[t]),col=t, lwd=1.5)
lines(c(colMeans(sp3)[t],0),c(colMeans(se3)[t],0),col=t, lwd=1.5)

abline(a=0,b=1, lty=2)

t<-4
points(mean(colMeans(sp1)),mean(colMeans(se1)),pch=t,col=t,cex=2)
points(mean(colMeans(sp2)),mean(colMeans(se2)),pch=t,col=t,cex=2)
points(mean(colMeans(sp3)),mean(colMeans(se3)),pch=t,col=t,cex=2)
lines(c(1,mean(colMeans(sp1))),c(1,mean(colMeans(se1))),col=t, lty=3, lwd=3)
lines(c(mean(colMeans(sp1)),mean(colMeans(sp2))),c(mean(colMeans(se1)),mean(colMeans(se2))),col=t, lty=3, lwd=3)
lines(c(mean(colMeans(sp2)),mean(colMeans(sp3))),c(mean(colMeans(se2)),mean(colMeans(se3))),col=t, lty=3, lwd=3)
lines(c(mean(colMeans(sp3)),0),c(mean(colMeans(se3)),0),col=t, lty=3, lwd=3)

legend("bottomright", legend=c("Rater 39","Rater 18","Rater 4","All"), col=c(1,2,3,4),pch=c(1,2,3,4),lty=c(1,1,1,3), lwd=c(1.5,1.5,1.5,1.5))
dev.off()

pdf("ROC2.pdf", width=6, height=6)

a <- matrix(nrow=3000, ncol=50)
b <- matrix(nrow=3000, ncol=50)
b.u <- matrix(nrow=3000, ncol=90)
U <- matrix(nrow=3000, ncol=90)
#beta <- NULL
#mu.beta <- NULL
mu.a <- NULL
#tau.beta <- NULL
tau.b <- NULL
for(i in 1:3000){
  a[i,] <- out[[i]]$a
  b[i,] <- out[[i]]$b
  b.u[i] <- mean(outer(out[[i]]$b,out[[i]]$u))
  U[i,]<- out[[i]]$u
  #    beta[i] <- out[[i]]$beta
  #    mu.beta[i] <- out[[i]]$mu_beta
  #    tau.beta[i] <- out[[i]]$tau_beta
  mu.a[i] <- mean(out[[i]]$a)
  #    tau.b[i] <- out[[i]]$tau_b
}

a.est <- mean(a)
b.est <- mean(b)
bu.est <- colMeans(b.u)
U.est <- colMeans(U)
#beta.est <- mean(beta)

a.mu.est <- mean(mu.a)

#beta.mu.est <- mean(mu.beta)


fpr <- function(h){
  temp <- mean(rowSums((pnorm(h-(a.est+b.est*U)))*(1-pnorm(U)))/rowSums(1-pnorm(U)))
}

se <- function(h){
  temp <- mean(rowSums((1-pnorm(h-(a.est+b.est*U)))*(pnorm(U)))/rowSums(pnorm(U)))
}


fpr.mu <- function(h){
  temp <- (1-pnorm(h-(a.mu.est)))
}

se.mu <- function(h){
  temp <- (1-pnorm(h-(a.mu.est)))
}


x <- sapply(seq(-10,10,length=1000), function(x) 1-fpr(x))
y <- sapply(seq(-10,10,length=1000), function(x) se(x))
lo <- loess(y~x)
plot(x,y, col="red", xlab="False Positive Rate", ylab="True Positive Rate",  type="l" , xlim=c(0,1), lwd=2)
#lines(x,predict(lo), col='red', lwd=2)
#x <- fpr.mu(seq(-4,10,length=1000))
#y <- se.mu(seq(-4,10,length=1000))
#lo <- loess(y~x)
#lines(x,predict(lo), col='green', lwd=2)
#points(x,y,col='green',pch=2)
abline(a=0,b=1, lty=2)

t<-4
points(mean(colMeans(sp1)),mean(colMeans(se1)),pch=t,col=t,cex=2)
points(mean(colMeans(sp2)),mean(colMeans(se2)),pch=t,col=t,cex=2)
points(mean(colMeans(sp3)),mean(colMeans(se3)),pch=t,col=t,cex=2)
lines(c(1,mean(colMeans(sp1))),c(1,mean(colMeans(se1))),col=t, lty=3, lwd=3)
lines(c(mean(colMeans(sp1)),mean(colMeans(sp2))),c(mean(colMeans(se1)),mean(colMeans(se2))),col=t, lty=3, lwd=3)
lines(c(mean(colMeans(sp2)),mean(colMeans(sp3))),c(mean(colMeans(se2)),mean(colMeans(se3))),col=t, lty=3, lwd=3)
lines(c(mean(colMeans(sp3)),0),c(mean(colMeans(se3)),0),col=t, lty=3, lwd=3)

t<-5
points(mean(sp1.mean),mean(se1.mean),pch=t,col=t,cex=2)
points(mean(sp2.mean),mean(se2.mean),pch=t,col=t,cex=2)
points(mean(sp3.mean),mean(se3.mean),pch=t,col=t,cex=2)
lines(c(1,mean(sp1.mean)),c(1,mean(se1.mean)),col=t, lty=3, lwd=3)
lines(c(mean(sp1.mean),mean(sp2.mean)),c(mean(se1.mean),mean(se2.mean)),col=t, lty=3, lwd=3)
lines(c(mean(sp2.mean),mean(sp3.mean)),c(mean(se2.mean),mean(se3.mean)),col=t, lty=3, lwd=3)
lines(c(mean(sp3.mean),0),c(mean(se3.mean),0),col=t, lty=3, lwd=3)

library(pROC)
roc1 <- roc(di ~ wij, simdataset, smooth=FALSE)
lines(1-roc1$specificities, roc1$sensitivities, type="l", col="black", lwd=2)


legend("bottomright", legend=c("Smoothed ROC","Empirical ROC 1","Empirical ROC 2","pROC"), col=c(2,4,5,1),pch=c(NA,4,5,NA),lty=c(1,3,3,1), lwd=c(1.5,1.5,1.5,1.5))

dev.off()



L <- NULL
for(i in 1:1000){
    temp<-out[[i]]$z - matrix(out[[i]]$a, nrow=80, ncol=50, byrow=T) - outer(out[[i]]$u,out[[i]]$b) 
    L[i] <-sum(temp^2)
}

hist(L, breaks=30, freq=F, main="Goodness-of-fit", xlab="", xlim=c(3500,4500))
lines(1:8000, dchisq(1:8000, 80*50, ncp = 0, log = FALSE), col="red")


g1 <- g2 <- g3 <- g4 <- NULL
for(i in 1:3000){
  g1[i] <- out[[i]]$good1
  g2[i] <- out[[i]]$good2
  g3[i] <- out[[i]]$good3
  g4[i] <- out[[i]]$good4
}

table(simdataset$wij)/(length(simdataset$wij))

mean(g1)
mean(g2)
mean(g3)
mean(g4)


quantile(g1, c(0.025,0.975))
quantile(g2, c(0.025,0.975))
quantile(g3, c(0.025,0.975))
quantile(g4, c(0.025,0.975))


> table(simdataset$wij)/(length(simdataset$wij))

1          2          3          4 
0.60088889 0.28422222 0.05400000 0.06088889 
> 
  > mean(g1)
[1] 0.6019884
> mean(g2)
[1] 0.2508393
> mean(g3)
[1] 0.07170268
> mean(g4)
[1] 0.07546963
> 
  > 
  > quantile(g1, c(0.025,0.975))
2.5%     97.5% 
  0.5896922 0.6139600 
> quantile(g2, c(0.025,0.975))
2.5%     97.5% 
  0.2340219 0.2666800 
> quantile(g3, c(0.025,0.975))
2.5%      97.5% 
  0.06445585 0.08027303 
> quantile(g4, c(0.025,0.975))
2.5%      97.5% 
  0.06325482 0.09015392
