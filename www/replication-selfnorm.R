###################################################
### Obtain self-normalized stats and process stats
###################################################

calc.self.stats.cont <- function(scores, ordervar){
  scorecom <- cbind(scores, ordervar)
  scorecom <- scorecom[order(ordervar),]
  scores <- as.matrix(scorecom[,-ncol(scorecom)])
  
  N <- nrow(scores)
  d <- ncol(scores)
  if(d>1)
  {
    forward <- apply(scores,2,cumsum);
    backward <- apply(scores[N:1,],2,cumsum);
  }
  if(d==1)
  {forward<-matrix(cumsum(scores),N,1);
  backward<-matrix(cumsum(scores[N:1]),N,1);}
  
  denomi<-array(0,c(d,d,N))
  x<-matrix(rep(1:N,d),N,d)/N
  numer<-t(t(forward)-t(x)*forward[N,])/sqrt(N)
  
  for(k in 1:(N-1))
  {
    x1<-matrix(rep(1:k,d),k,d)/k
    x2<-matrix(rep(1:(N-k),d),(N-k),d)/(N-k)
    
    if (k>1)
      denomi1<-t(t(forward[1:k,])-t(x1)*forward[k,])
    else
      denomi1<-t(t(t(forward[1:k,]))-t(x1)*forward[k,])
    
    if (k<(N-1))
      denomi2<-t(t(backward[1:(N-k),])-t(x2)*backward[N-k,])
    else
      denomi2<-t(t(t(backward[1:(N-k),]))-t(x2)*backward[N-k,])
    
    denomi[,,k]<-{t(denomi1)%*%denomi1+t(denomi2)%*%denomi2}/N^2
  }
  
  subG<-rep(NA,N-1)
  G <- list()
  
  for(k in 1:(N-1))
  {
    subG[k]<-numer[k,]%*%solve(denomi[,,k])%*%t(t(numer[k,]))
    
  }
  
  G <- list(max(abs(subG)), subG)
  
  return(G)
}

###################################################
### Fit models and obtain stats
###################################################
if (file.exists("statsval1.rda") & file.exists("statsval2.rda") ){
  load("statsval1.rda")
  load("statsval2.rda")
} else{
  library("mlmRev")
  library("lmerTest")
  library("merDeriv")
  data("bdf")

  
  m2 <- lmer(langPOST ~ langPRET + (1 | schoolNR), data = bdf,
             REML = FALSE)
  
  psi1 <- estfun.lmerMod(m2, level = 1)
  
  mz_gefp1  <- as.matrix(psi1[,2])
  mz_gefp2  <- as.matrix(psi1[,4])
  statsval1 <- calc.self.stats.cont(mz_gefp1, bdf$IQ.verb)
  statsval2 <-  calc.self.stats.cont(mz_gefp2, bdf$IQ.verb)
  
  save(statsval1, file = "statsval1.rda")
  save(statsval2, file = "statsval2.rda")
}

###################################################
### Produce Figure 
###################################################
par(mfrow = c(1, 2))
set.seed(1090)
plot(sort(bdf$IQ.verb)[-1], statsval1[[2]], xlab = "Verbal IQ", ylab = "Empirical fluctuation process",  type = "l", ylim = c(0, 55), 
     main = "Fixed effect coefficient")
abline(h=40.1, col="red")
plot(sort(bdf$IQ.verb)[-1], statsval2[[2]], xlab = "Verbal IQ", ylab = "Empirical fluctuation process",  type = "l", ylim = c(0, 55), 
     main = "Residual variance")
abline(h=40.1, col="red")