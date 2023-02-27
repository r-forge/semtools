## Data-generating process
dgp <- function(nobs = 960, diff = 4, nlevels = 4, 
                parms= "beta0",
                change = 0.5)
{
  ## Generates data from a multilevel linear model that has parameter change.
  ## Also generates an ordinal auxiliary
  ## variable related to the invariance (loosely called "days").
  
  stopifnot(require("mvtnorm") & require("Matrix"))
  subject = nobs/10
  
  ses <- diff
  sampsize <- nobs
  sampsize.per.nobs <- nobs %/% nlevels
  
  ## Deal with cases where nobs/nlevels is not an integer
  cognitive <- rep(NA, subject)
  cognitive[1:subject] <- rep(1:nlevels, each = subject/nlevels)
  cognitiveobs <- rep(1:nlevels, each = nobs/nlevels)
  # if (any(is.na(cognitive))) cognitive[is.na(cognitive)] <- sample((1:nlevels),
  #    sum(is.na(cognitive)))
  cognitive <- c(rep(runif(subject*change,13,16),
                     each = (nobs/subject)), 
                 rep(runif(subject*(1-change),16,18),
                     each = (nobs/subject)))
  
  
  
  
  half.level <- floor(nlevels * change)
  
  ## Define parameter vectors/matrices:
  #10.47, 10.28 251.41
  theta <- c(251.41, 10.47, 23.78, 23.78*5.7171*0.08, 5.7171, 25.59^2)
  #asymp <- c(24.443385, 4.966558, 858.854901, 125.384457, 35.112119, 329.674893)
  asymp <- c(96.42830, 18.93079, 4186.74171, 585.49750, 160.76751, 1042.52364) 
  # sqrt(c(43.98, 10.28, 7.036*10^4, 1.838*10^3,
  #           8.560*10^2, 5.957*10^3) *18)        
  # beta: 
  beta <- theta[1:2]
  
  # residual:
  residual <- theta[length(theta)]
  
  # varcov
  varcov <- matrix(NA, 2, 2)
  varcov[lower.tri(varcov,diag=TRUE)] <- theta[3:5]
  varcov <- as.matrix(forceSymmetric(varcov, "L"))
  diag(varcov) <- diag(varcov)^2
  
  ## Initialize data matrix:
  datmat <- matrix(NA, sampsize, 4)
  colnames(datmat) <- c("Reaction", "Days", "Subject", "Cog")
  
  ## Generate data for invariant individuals:
  # n.invariant <- sum(cognitive <= half.level)
  n.invariant <- round(sampsize*change, 0)
  ## Based on beta, b, residual, generate data:
  ## response variable "Reaction"
  b <- as.vector(t(rmvnorm(subject, mean = c(0, 0),
                           varcov)))
  datmat[1:n.invariant,2] <- rep(1: (nobs/subject), subject/half.level)
  ## identifying subjects "Subject"
  datmat[1:n.invariant,3] <- rep(1: (subject/half.level),
                                 each = (nobs/subject))
  Z <- matrix(0, nobs, (subject * 2))
  for (i in 1: subject) {
    Z[(1 + (i-1) * (nobs/subject)):(i * (nobs/subject)), (2*i-1)] <- 1
    Z[(1 + (i-1) * (nobs/subject)):(i * (nobs/subject)), (2*i)] <- 1:(nobs/subject)
  }
  
  datmat[1:n.invariant,1] <- rnorm(n.invariant,
                                   mean=(beta[1] +
                                           beta[2] *
                                           datmat[1:n.invariant, 2] +
                                           Z[1:n.invariant, ] %*% b),
                                   sd=sqrt(residual))
  
  
  ## Now do the same for invariant-violating levels
  n.vary <- nlevels - half.level
  tmp.n <- (nobs - n.invariant)
  
  if (parms=="beta0") {parmsvector=c(1, 0, 0, 0, 0, 0)}
  if (parms=="beta1") {parmsvector=c(0, 1, 0, 0, 0, 0)}
  if (parms=="var1") {parmsvector=c(0, 0, 1, 0, 0, 0)}
  if (parms=="cov") {parmsvector=c(0, 0, 0, 1, 0, 0)}
  if (parms=="var2") {parmsvector=c(0, 0, 0, 0, 1, 0)}
  if (parms=="residual") {parmsvector=c(0, 0, 0, 0, 0, 1)}
  if (parms=="resbeta") {parmsvector=c(0, 1, 0, 0, 0, 1)}
  
  theta.tmp <- theta + ((ses*asymp*parmsvector)/sqrt(subject/2))
  # beta, lamtheta or varcov, residual:
  tmp.beta <- theta.tmp[1:2]
  tmp.residual <- theta.tmp[length(theta.tmp)]
  tmp.varcov <- matrix(NA, 2, 2)
  tmp.varcov[lower.tri(tmp.varcov,diag=TRUE)] <- theta.tmp[3:5]
  tmp.varcov <- as.matrix(forceSymmetric(tmp.varcov, "L"))
  diag(tmp.varcov) <- diag(tmp.varcov)^2
  
  tmp.b <- as.vector(t(rmvnorm(subject, mean = c(0, 0),
                               tmp.varcov)))
  
  ## Based on beta, varcov, residual, generate data:
  ## response variable "Reaction"
  
  datmat[(n.invariant+1):nobs,2] <- rep(1: (nobs/subject), (subject/n.vary))
  ## identifying subjects "Subject"
  datmat[(n.invariant+1):nobs,3] <- rep(((subject/half.level)+1): subject, 
                                        each = (nobs/subject))
  
  datmat[(n.invariant+1):nobs,1] <- rnorm(tmp.n,
                                          mean = (tmp.beta[1] +
                                                    tmp.beta[2] *
                                                    datmat[(n.invariant+1):nobs, 2] +
                                                    Z[(n.invariant+1):nobs, ] %*% tmp.b),
                                          sd=sqrt(tmp.residual))
  
  # datmat[,3] <- as.factor(datmat[,3])
  datmat[,4] <- cognitiveobs
  list(as.data.frame(datmat), cognitive)
}

## Evaluate power simulation on a single dgp() scenario
testpower <- function(nrep = 20, size = 0.05, ordfun = NULL,
                      parnum = parnum, sim = "sim1",
                      test = test, estimator = estimator,
                      verbose = TRUE, ...)
{
  pval <- matrix(rep(NA, length(test) * nrep), ncol = length(test))
  colnames(pval) <- test
  dfun <- "dgp"
  for(i in 1:nrep) {
    ## simulate data
    d <- do.call(dfun,list(...))
    ## fit model
    if (substr(sim,4,4)==1){
      mz <- tryCatch(mzlmmfit2.lme4(d[[1]]), error=function(e) e, warning=function(w) w)
      
    }
    if(substr(sim,4,4)==2){
      mz <- tryCatch(mzfit.lavaan(d[[1]]), error=function(e) e, warning=function(w) w)
    }
    
    ## Requires estfun.ltm/estfun.ltm
    if(!inherits(mz, "error") & !inherits(mz, "warning")){
      ord_gefp <- vector("list",length(parnum))
      for (j in 1: length(parnum)){
        
        if (estimator=="MMLinfo" & substr(sim,4,4)==1) {
          ord_gefp[[j]]  <-  tryCatch(gefp(mz, fit = NULL, vcov = info_full_lme4,
                                           order.by = d[[2]],
                                           sandwich = FALSE,
                                           scores = estfun.lmerMod(level = 1),
                                           parm = eval(parse(text=parnum[j]))),
                                      error=function(e) e, warning=function(w) w)
        }
        
        if (estimator=="MMLinfo" & substr(sim,4,4)==2) {
          ord_gefp[[j]]  <-  tryCatch(gefp(mz, fit = NULL, vcov = info_full_lavaan,
                                           order.by = d[[2]],
                                           sandwich = FALSE,
                                           scores = estfun.lavaan,
                                           parm = eval(parse(text=parnum[j]))),
                                      error=function(e) e, warning=function(w) w)
        }
        if (estimator=="MMLnull") {
          ord_gefp[[j]]  <- tryCatch(gefp(mz, fit = NULL, vcov = NULL,
                                          order.by = d[[2]],
                                          sandwich = FALSE,
                                          scores = estfun.lmerMod(level = 1),
                                          parm = eval(parse(text=parnum[j]))),
                                     error=function(e) e, warning=function(w) w)
        }
        
        if(!inherits(ord_gefp[[j]], "error") & !inherits(ord_gefp[[j]], "warning")) {
          pval[i, (j-1)*3+1] <- sctest(ord_gefp[[j]], functional = maxBB)$p.value
          pval[i, (j-1)*3+2] <- sctest(ord_gefp[[j]], functional = meanL2BB)$p.value
          pval[i, (j-1)*3+3] <- sctest(ord_gefp[[j]], functional = supLM(0.1))$p.value
        }
      }
    }
  }
  rval <- colMeans(pval < size, na.rm = TRUE)
  
  return(rval)
}

## Loop over scenarios
simulation <- function(diff = seq(0, 4, by = 0.25),
                       parms=c("beta", "varcov", "residual"),
                       nobs = 240, nlevels = 10, estimator = "MMLNULL",
                       subject = 24, change = seq(0.1, 0.9, 0.1),
                       sim="sim1",
                       verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nlevels = nlevels, nobs = nobs,
                     parms = parms,
                     change = change)
  nprs <- nrow(prs)
  
  parnum <- c(1, 2, 3, 4, 5, 6)
  
  
  #tname <- c("ordmax","ordwmax","catdiff") # had suplm here
  tname <- c("dmax", "CvM", "supLM")
  
  test <- paste(rep(tname,length(parnum)), rep(parnum, each=length(tname)),sep="")
  ntest <- length(test)
  
  ## Only simulate critical values if we need it.
  if (file.exists("critvals.rda")){
    load("critvals.rda")
    ## Check to make sure it has what we need; if not,
    ## rerun critvals.
    if (!all(nlevels %in% as.numeric(names(cval)))) critvals(nlevels)
  } else {
    cval <- critvals(nlevels)
    save(cval, file="critvals.rda")
  }
  cval.conds <- as.numeric(names(cval))
  
  pow <- matrix(rep(NA,ntest*nprs),ncol=ntest)
  do.parallel <- require("parallel")
  # do.parallel=FALSE
  if (do.parallel){
    pow <- mclapply(1:nprs, function(i){
      testpower(diff = prs$diff[i], nobs = prs$nobs[i],
                parms = prs$parms[i],
                parnum = parnum, nlevels = prs$nlevels[i],
                change = prs$change[i],
                test = test, sim = sim, estimator = estimator,
                verbose = verbose,...)},
      mc.cores = round(.75*detectCores()),
      mc.preschedule = FALSE)
    pow <- t(matrix(unlist(pow), ntest, nprs))
  } else {
    pow <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    for(i in 1:nprs) {
      if(verbose) print(prs[i,])
      
      ## Find critical values we need
      ordfun <- cval[[which(cval.conds == prs$nlevels[i])]]
      
      pow[i,] <- testpower(diff = prs$diff[i], nobs = prs$nobs[i],
                           parms = prs$parms[i],parnum = parnum,
                           nlevels = prs$nlevels[i], 
                           change = prs$change[i],
                           test = test, sim = sim, estimator = estimator,
                           verbose = verbose,...)
    }
  }
  
  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(rep(tname,length(parnum)), each = nprs),
                      levels=c("dmax", "CvM", "supLM"))
  rval$pars <- factor(rep(rep(parnum,each=length(tname)),each=nprs),
                      levels=parnum)
  rval$nobs <- factor(rval$nobs)
  rval$parms <- factor(rval$parms)
  rval$power <- as.vector(pow)
  return(rval)
}

### Simulate critical values
critvals <- function(nlevels = c(4, 8, 12), verbose = TRUE, ...)
{
  do.parallel <- require("parallel")
  if (do.parallel){
    cval <- mclapply(nlevels, function(x) ordL2BB(rep(1, x)/x, nobs = 5000),
                     mc.cores = length(nlevels))
  } else {
    cval <- lapply(nlevels, function(x) ordL2BB(rep(1, x)/x, nobs = 5000))
  }
  names(cval) <- nlevels
  
  save(cval, file="critvals.rda")
}


info_full_lme4 <- function(object,...)solve(vcov.lmerMod(object, full = TRUE)*
                                              length(unique(object@frame$Subject)))


info_full_lavaan <- function(object,...)solve(vcov(object, remove.duplicated=TRUE)*
                                                nrow(object@Data@X[[1]]))

mzlmmfit2.lme4 <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
  data <- as.data.frame(data)
  stopifnot(c("Reaction", "Subject", "Days") %in% names(data))
  data$Subject <- factor(data$Subject)
  ## require package
  stopifnot(require("lme4"))
  
  ## fit model
  rval <- try(lmer(Reaction ~ Days + (Days|Subject), data,
                   REML=FALSE))
  return(rval)
}

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

testpower2 <- function(nrep = 30, size = "95%", 
                      verbose = TRUE, ...)
{
  parnum <- c(1, 2, 3, 4, 5, 6)
  
  tname <- c("ordmax") # had suplm here
  
  statsval <- matrix(rep(NA, 6* nrep), ncol = 6)
  colnames(statsval) <- paste(rep(tname, length(parnum)), rep(parnum,
                                                              each = length(tname)),sep="")
  
  dfun <- "dgp"
  for(i in 1:nrep) {
    data <- do.call(dfun,list(...))
    
    mz <- tryCatch(try(mzlmmfit2.lme4(data[[1]]), silent = TRUE),
                   error=function(e) e, warning=function(w) w)
    
    if(!is(mz, "warning") & !is(mz, "try-error")){
      psi <- estfun.lmerMod(mz, level = 1) %*% root.matrix(vcov(mz, full = TRUE))
      #psi <- estfun.lmerMod(mz, level = 1)
      mz_gefp1  <- as.matrix(psi[,1])
      mz_gefp2 <- as.matrix(psi[,2])
      mz_gefp3 <- as.matrix(psi[,3])
      mz_gefp4 <- as.matrix(psi[,4])
      mz_gefp5 <- as.matrix(psi[,5])
      mz_gefp6 <- as.matrix(psi[,6])
      
      statsval[i, 1] <- calc.self.stats.cont(mz_gefp1, data[[2]])[[1]]
      statsval[i, 2] <- calc.self.stats.cont(mz_gefp2, data[[2]])[[1]]
      statsval[i, 3] <- calc.self.stats.cont(mz_gefp3, data[[2]])[[1]]
      statsval[i, 4] <- calc.self.stats.cont(mz_gefp4, data[[2]])[[1]]
      statsval[i, 5] <- calc.self.stats.cont(mz_gefp5, data[[2]])[[1]]
      statsval[i, 6] <- calc.self.stats.cont(mz_gefp6, data[[2]])[[1]]
    }
  }
  
  
  crit <- matrix(c(29.6,
                   40.1,
                   68.6,
                   81.5,
                   103.6,
                   160.0,
                   183.8,
                   218.8,
                   318.3), 3, 3, byrow = FALSE)
  rownames(crit) <- c("90%", "95%", "99%")
  colnames(crit) <- c(1, 3, 6)
  
  rval <- colMeans(t(apply(statsval, 1, function(x) x > c(crit[which(rownames(crit)==size), 1],
                                                          crit[which(rownames(crit)==size), 1],
                                                          crit[which(rownames(crit)==size), 1],
                                                          crit[which(rownames(crit)==size), 1],
                                                          crit[which(rownames(crit)==size), 1],
                                                          crit[which(rownames(crit)==size), 1]))), na.rm = TRUE)
  if(verbose) print(rval)
  
  return(rval)
}

## Loop over scenarios
simulation2 <- function(diff = seq(0, 4, by = 0.25), parms = "beta0",
                       nlevels = 4, change = 0.5, 
                       nobs = c(120, 480, 960), 
                       verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nlevels = nlevels, nobs = nobs,
                     parms = parms,
                     change = change)
  nprs <- nrow(prs)
  
  test <- c("ordmax1", 
            "ordmax2", 
            "ordmax3",
            "ordmax4",
            "ordmax5",
            "ordmax6")
  ntest <- length(test)
  
  
  pow <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  #do.parallel <- require("parallel")
  
  do.parallel <- TRUE
  if(do.parallel){
    pow <- mclapply(1:nprs, function(i){
      testpower2(diff = prs$diff[i], nobs = prs$nobs[i],
                parms = prs$parms[i],
                nlevels = prs$nlevels[i],
                change = prs$change[i],
                verbose = verbose, ...)},
      mc.cores = round(.75*detectCores()),
      mc.preschedule = FALSE)
    pow <- t(matrix(unlist(pow), ntest, nprs))
  } else{
    for(i in 1:nprs) {
      if(verbose) print(prs[i,])
      
      pow[i,] <- testpower2(diff = prs$diff[i], nobs = prs$nobs[i],
                           parms = prs$parms[i],
                           nlevels = prs$nlevels[i],
                           change = prs$change[i],
                           verbose = verbose, ...)
    }
  }
  
  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(c("ordmax",
                            "ordmax", 
                            "ordmax", 
                            "ordmax", 
                            "ordmax",
                            "ordmax"), each = nprs),
                      levels = c("ordmax"))
  rval$pars <- factor(rep(c("1", 
                            "2", 
                            "3", 
                            "4",
                            "5",
                            "6"), each = nprs),
                      levels = c("1", "2", "3", "4", "5", "6"))
  rval$nobs <- factor(rval$nobs)
  rval$power <- as.vector(pow)
  return(rval)
}


if (FALSE){
  ###################################################
  ### Load Packages
  ###################################################
  library("lavaan")
  library("strucchange")
  library("merDeriv")
  library("lattice")
  source("simulation_replication.R")

  
  ###################################################
  ### Original Simulation
  ###################################################
  sim3 <- simulation(nobs = c(240, 480, 960),
                     nrep = 1000, subject = c(24, 48, 96), nlevels = 4,
                     change = 0.5, 
                     diff = c(0, 1, 2, 3, 4), sim = "sim1",
                     estimator = "MMLinfo",
                     parms = c("beta0", "beta1"))
  save(sim3, file = "sim3.rda")
  
  sim3$test <- factor(as.character(sim3$test), 
                      levels = c("dmax", "CvM", "maxLM"), 
                      labels = c("dmax", "CvM", "maxLM"))
  ## make labels 
  levels(sim3$pars) <-c("beta0", "beta1","var1", "cov", "var2", "residual")
  sim3$nobs <- factor(sim3$nobs)
  levels(sim3$nobs) <- paste("n=", levels(sim3$nobs), sep="")
  parlabs <-  c(expression(beta[0]), expression(beta[1]), expression(sigma[0]^2),
                expression(sigma["01"]), expression(sigma[1]^2), 
                expression(sigma[r]^2))

  ####################################################
  ### Self normalized Simulation
  ###################################################
  sim2 <- simulation2(diff = seq(0, 4, by = 1),
                     nobs = c(240, 480, 960), 
                     nrep = 1000, parms = c("beta0", "beta1"), nlevels = 4,
                     change = 0.5, size = "95%")
  
  ## make labels
  sim2$test <- factor(as.character(sim2$test), 
                      levels = c("ordmax"), 
                      labels = c("SN"))
  parlabs <- c(expression(beta[0]), expression(beta[1]), expression(sigma[0]^2),
               expression(sigma["01"]), expression(sigma[1]^2), 
               expression(sigma[r]^2))
  #levels(sim1$pars) <- c("lambda1","beta1", "var1", "var2", "residual")
  sim2$nobs <- factor(sim2$nobs)
  levels(sim2$nobs) <- paste("n=", levels(sim2$nobs), sep="")
  save(sim2, file="sim2.rda")
  
  ###################################################
  ### Output Figures and Tables
  ###################################################
  ## Figure 1 
  trellis.par.set(theme = canonical.theme(color = FALSE))
  mykey <- simpleKey(levels(sim3$test), points = TRUE, lines = TRUE)
  xyplot(power ~ diff | pars + nobs, group = ~ test, data = sim3,
         subset = (parms == "beta0" & 
                     diff %in% c(0,1,2,3,4)),
         type = "b",
         xlab = expression(paste("Violation Magnitude (", beta[0], ")")),
         ylab = "Power", key = mykey, as.table = TRUE, ylim = c(0, 1.2), 
         strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                             factor.levels, ...)
           } else {
             strip.default(which.given, factor.levels =
                             parlabs, ...)
           }
         })
  
  ## Figure 2
  trellis.par.set(theme = canonical.theme(color = FALSE))
  mykey <- simpleKey(levels(sim3$test), points = TRUE, lines = TRUE)
  xyplot(power ~ diff | pars + nobs, group = ~ test, data = sim3,
         subset = (parms == "beta1" & 
                     diff %in% c(0,1,2,3,4)),
         type = "b",
         xlab = expression(paste("Violation Magnitude (", beta[1], ")")),
         ylab = "Power", key = mykey, as.table = TRUE, ylim = c(0, 1.2), 
         strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                             factor.levels, ...)
           } else {
             strip.default(which.given, factor.levels =
                             parlabs, ...)
           }
         })
  
  ## Figure 3
  trellis.par.set(theme = canonical.theme(color = FALSE))
  mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
  xyplot(power ~ diff | pars + nobs, group = ~ test, data = sim2,
         subset = (parms == "beta0" & 
                     diff %in% c(0,1,2,3,4)),
         type = "b",
         xlab = expression(paste("Violation Magnitude (", beta[0], ")")),
         ylab = "Power", key = mykey, as.table = TRUE,
         strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                             factor.levels, ...)
           } else {
             strip.default(which.given, factor.levels =
                             parlabs, ...)
           }
         })
  
  ## Figure 4
  trellis.par.set(theme = canonical.theme(color = FALSE))
  mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
  xyplot(power ~ diff | pars + nobs, group = ~ test, data = sim2,
         subset = (parms == "beta1" & 
                     diff %in% c(0,1,2,3,4)),
         type = "b",
         xlab = expression(paste("Violation Magnitude (", beta[1], ")")),
         ylab = "Power", key = mykey, as.table = TRUE,
         strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                             factor.levels, ...)
           } else {
             strip.default(which.given, factor.levels =
                             parlabs, ...)
           }
         })
  
  ## Table 1 
  mz_sim_sub <- subset(sim2, parms %in% c("beta0"))
  mz_sim_sub$pars <- factor(mz_sim_sub$pars)
  
  levels(mz_sim_sub$pars) <- c("$\\beta_{0}$", "$\\beta_{1}$", "$\\sigma_{0}^2$",
                               "$\\sigma_{01}$", "$\\sigma_{1}^2$", 
                               "$\\sigma_{r}^2$")
  levels(mz_sim_sub$pars) <- gsub("testing", "", levels(mz_sim_sub$pars), fixed = TRUE)
  levels(mz_sim_sub$nobs) <- gsub("obs","", levels(mz_sim_sub$nobs),fixed = TRUE)
  levels(mz_sim_sub$test) <- c("SN")
  mz_tab <- round(ftable(100 * xtabs(power ~ nobs + pars + test + diff,
                                     data = mz_sim_sub, subset = diff %in% c(seq(0, 4, by = 1))),
                         col.vars = "diff"), digits = 1)
  mz_tab <- format(mz_tab, quote = FALSE)[-2, -4]
  mz_tab[1,] <- c("Observation", "Tested Parameter", "Statistic", c(format(seq(0, 4, by = 1))))
  mz_tab <- paste(apply(mz_tab, 1, paste, collapse = " & "), "\\\\")
  mz_tab[c(1, length(mz_tab))] <- paste(mz_tab[c(1, length(mz_tab))], "\\hline")
  writeLines(mz_tab)
  
  ## Table 2
  levels(mz_sim_sub$pars) <- c("$\\beta_{0}$", "$\\beta_{1}$", "$\\sigma_{0}^2$",
                               "$\\sigma_{01}$", "$\\sigma_{1}^2$", 
                               "$\\sigma_{r}^2$")
  levels(mz_sim_sub$pars) <- gsub("testing", "", levels(mz_sim_sub$pars), fixed = TRUE)
  levels(mz_sim_sub$nobs) <- gsub("obs","", levels(mz_sim_sub$nobs),fixed = TRUE)
  levels(mz_sim_sub$test) <- c("SN")
  mz_tab <- round(ftable(100 * xtabs(power ~ nobs + pars + test + diff,
                                     data = mz_sim_sub, subset = diff %in% c(seq(0, 4, by = 1))),
                         col.vars = "diff"), digits = 1)
  mz_tab <- format(mz_tab, quote = FALSE)[-2, -4]
  mz_tab[1,] <- c("Observation", "Tested Parameter", "Statistic", c(format(seq(0, 4, by = 1))))
  mz_tab <- paste(apply(mz_tab, 1, paste, collapse = " & "), "\\\\")
  mz_tab[c(1, length(mz_tab))] <- paste(mz_tab[c(1, length(mz_tab))], "\\hline")
  writeLines(mz_tab)
}