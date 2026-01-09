intconst <- function(...) {
  ## send in a simple cfa command and ordinal dataset, receive model syntax for 
  ## integer identification constraints
  ddd <- list(...)

  if (!("do.fit" %in% names(ddd))) {
    origdf <- TRUE
    ddd$do.fit <- FALSE
  } else {
    origdf <- ddd$do.fit
    ddd$do.fit <- FALSE
  }

  ## set up initial model
  LAV <- do.call("cfa", ddd)

  if (length(lavInspect(LAV, 'ordered')) != length(lavNames(LAV))) stop("ERROR: not all ovs are ordinal; code doesn't work for this situation.")

  dat <- lavInspect(LAV, 'data')

  nlevs <- apply(dat, 2, function(x) length(unique(x)))
  thetaparm <- lavInspect(LAV, 'parameterization') == "theta"

  mats <- lavInspect(LAV, 'est')
  freeloads <- which(mats$lambda != 0, arr.ind = TRUE)
  varsperlv <- tapply(freeloads[,'row'], freeloads[,'col'], function(x) rownames(mats$lambda)[x])
  varnames <- rownames(mats$lambda)

  newsyn <- NULL
  ## loadings
  for (i in 1:length(varsperlv)) {
    ovnums <- match(varsperlv[[i]], varnames)
    newsyn <- paste0(newsyn, "\n f", i, " =~ ", paste(paste0("ld", ovnums), varsperlv[[i]], sep = "*", collapse = " + "))
  }
    
  ## intercepts, residual variances and thresholds (fix intercept to 1 if binary)
  for (i in 1:length(nlevs)) {
    newsyn <- paste0(newsyn, "\n ", varnames[i], " ~ int", i, "*1", ifelse(nlevs[i] == 2, " + 0*1", ""))

    newsyn <- paste0(newsyn, "\n ", varnames[i], " ~", ifelse(thetaparm, "", "*"), "~ resid", i, "*", varnames[i])

    newsyn <- paste0(newsyn, "\n ", varnames[i], " | ", paste(paste0(letters[i], 1:(nlevs[i] - 1)), paste0("t", 1:(nlevs[i] - 1)), sep = "*", collapse = " + "))
  }

  newsyn <- paste0(newsyn, "\n")
  
  ## latent means/vars/covars
  for (i in 1:length(varsperlv)) {
    ovnums <- match(varsperlv[[i]], varnames)
    ovno2 <- ovnums[nlevs[ovnums] > 2]
    newsyn <- paste0(newsyn, "\n f", i, " ~ lvmn", i, "*1")
    newsyn <- paste0(newsyn, "\n f", i, " ~~ lvvar", i, "*f", i, ifelse(length(ovno2) == 0, paste0(" + 1*f", i), ""))
  }

  newsyn <- paste0(newsyn, "\n")

  if (ncol(mats$lambda) > 1 & !lavInspect(LAV, "options")$orthogonal) {
    for (i in 1:(ncol(mats$lambda) - 1)) {
      newsyn <- paste0(newsyn, "\n f", i, " ~~ ", paste0("f", (i + 1):ncol(mats$lambda), collapse = " + "))
    }
  }

  newsyn <- paste0(newsyn, "\n")
  
  ## intercept, loading, threshold restrictions
  for (i in 1:length(varsperlv)) {
    ovnums <- match(varsperlv[[i]], varnames)
    ovno2 <- ovnums[nlevs[ovnums] > 2]
    if (length(ovno2) > 0) {
      newsyn <- paste0(newsyn, "\n ", paste0("int", ovnums, collapse = " + "), " == 0")
    }
    newsyn <- paste0(newsyn, "\n ", paste0("ld", ovnums, collapse = " + "), " == ", length(ovnums))
    
    Kmax <- max(nlevs[ovnums])
    newsyn <- paste0(newsyn, "\n ", paste0(letters[ovnums], "1 == ", round(.5 + Kmax/nlevs[ovnums], 3), collapse = "\n "))
    ## only variables with >2 categories have upper thresholds
    if (length(ovno2) > 0) {
      newsyn <- paste0(newsyn, "\n ", paste0(letters[ovno2], nlevs[ovno2] - 1, " == ", round(.5 + Kmax * (nlevs[ovno2] - 1)/nlevs[ovno2], 3), collapse = "\n "))
    }
    newsyn <- paste0(newsyn, "\n")
  }

  newsyn
}

if (FALSE) {
  library(lavaan)
  source("integer_constraint.R")
  
  makeord <- function(Data, vars = NULL, ncat = 2){
    if(length(vars) == 0) vars <- 1:NCOL(Data)

    if(length(ncat) != 1 && length(ncat) != length(vars)) stop("bad ncat")

    ## for differing numbers of categories
    Data <- rbind(Data, NA)
    Data[nrow(Data), vars] <- ncat
  
    Data[,vars] <- apply(Data[,vars,drop = FALSE], 2,
                         function(x){
                           nc <- tail(x, 1)
                           tmpp <- (1/nc) + runif(1, -.1/nc, .1/nc)
                           brks <- c(min(x, na.rm = TRUE) - .1, seq(quantile(x, tmpp, na.rm = TRUE),
                                                                    quantile(x, 1 - tmpp, na.rm = TRUE),
                                                                    length.out = (nc-1)), max(x, na.rm = TRUE) + .1)
                           cut(x, breaks = brks, labels = FALSE)})

    Data[-nrow(Data),]
  }

  hsdat <- makeord(HolzingerSwineford1939[, paste0('x', 1:9)], ncat = 3)

  HS.model <- ' visual  =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed   =~ x7 + x8 + x9 '
  
  modsyn <- intconst(HS.model, data = hsdat, ordered = TRUE)

  cat(modsyn)

  fit <- lavaan(modsyn, data = hsdat, ordered = TRUE)
  summary(fit)

  fit2 <- cfa(HS.model, data = hsdat, ordered = TRUE)

  fitMeasures(fit)[1:4]
  fitMeasures(fit2)[1:4]

  ## theta parameterization
  modsyn <- intconst(HS.model, data = hsdat, ordered = TRUE, parameterization = "theta")
  fit3 <- lavaan(modsyn, data = hsdat, ordered = TRUE, parameterization = "theta")
  fitMeasures(fit3)[1:4]

  ## differing numbers of categories
  hsdat2 <- makeord(HolzingerSwineford1939[, paste0('x', 1:9)], ncat = rep(2:4, 3))
  modsyn <- intconst(HS.model, data = hsdat2, ordered = TRUE)
  cat(modsyn)
  
  fit <- lavaan(modsyn, data = hsdat2, ordered = TRUE)
  fit2 <- cfa(HS.model, data = hsdat2, ordered = TRUE)

  fitMeasures(fit)[1:4]
  fitMeasures(fit2)[1:4]

  ## two categories
  hsdat <- makeord(HolzingerSwineford1939[, paste0('x', 1:9)], ncat = 2)
  modsyn <- intconst(HS.model, data = hsdat, ordered = TRUE)

  fit4 <- lavaan(modsyn, data = hsdat, ordered = TRUE)
  fit5 <- cfa(HS.model, data = hsdat, ordered = TRUE)

  fitMeasures(fit4)[1:4]
  fitMeasures(fit5)[1:4]
}
