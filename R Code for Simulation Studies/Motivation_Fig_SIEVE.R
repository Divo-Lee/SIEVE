#####################################################
#SIEVE: One-stop differential expression, variability, and skewness
#analyses using RNA-Seq data
#Authors: Hongxiang Li and Tsung Fei Khang
#Email: chelsea.divo@hotmail.com
#Date: 22 May 2023
#R Codes for motivation
#Part 1: Motivational Figures (Fig. 1)
#####################################################

## R packages downloaded
# install.packages("compositions")
# BiocManager::install("polyester")
# BiocManager::install("edgeR", force = T)
# devtools::install_github("Divo-Lee/clrDV")
# install.packages("sn")
library(polyester); library(compositions)
library(edgeR); library(sn); library(SIEVE)

##############################################
### GSE123658_read_counts.gene_level.txt.gz  
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123658
### Valentim Data                            
##############################################
Data = read.table(gzfile('...GSE123658_read_counts.gene_level.txt.gz'), sep="\t", header = T, check.names = TRUE, row.names=1)
dim(Data)

# filter
CPM <- cpm(Data)
keep <- (rowMeans(CPM[,1:43]) > 0.5  & rowMeans(CPM[,44:82]) > 0.5 & apply(Data[,1:43], 1, function(x) length(x[x==0])/length(x)) < 0.85 & apply(Data[,44:82], 1, function(x) length(x[x==0])/length(x)) < 0.85)
Data <- Data[keep, ]; dim(Data)

# h = healthy volunteers
# d1 = type 1 diabetic patients
 # data_h <- Data[,1:43]
data_d1 <- Data[,44:82]; dim(data_d1)


################################
### GSE150318_counts.csv.gz    
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150318
### Kelmer Data                
################################
Data2 = read.csv(gzfile('...GSE150318_counts.csv.gz'), header = T, check.names = TRUE, row.names=1)
dim(Data2)

data_10weeks <- Data2[ , seq(1,228,2)]
data_20weeks <- Data2[ , seq(2,228,2)]
Data2 <- cbind(data_10weeks, data_20weeks)
# filter
CPM <- cpm(Data2)
keep <- (rowMeans(CPM[,1:114]) > 0.5 & rowMeans(CPM[,115:228]) > 0.5 & apply(Data2[,1:114], 1, function(x) length(x[x==0])/length(x)) < 0.85 & apply(Data2[,115:228], 1, function(x) length(x[x==0])/length(x)) < 0.85)
Data2 <- Data2[keep, ]; dim(Data2)
data_10weeks <- Data2[, 1:114]; dim(data_10weeks)
 #data_20weeks <- Data2[, 115:228]; dim(data_20weeks)


#################################
### GSE179250_gene_reads.txt.gz 
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179250
### Zhou Data (Liver)           
#################################
Data_liver <- read.table('.../GSE179250_gene_reads.txt', sep = "\t", header = T, check.names = T, row.names = 1)
dim(Data_liver)
#colnames(Data_liver)
#Data_liver[1:5, 1:7]
Data_liver <- Data_liver[, -c(1:3)]
CPM <- cpm(Data_liver)
keep <- (rowMeans(CPM) > 0.5 &
         apply(Data_liver, 1, function(k) length(k[k==0])/length(k)) < 0.85)
Data_liver <- Data_liver[keep, ]
dim(Data_liver)


#####################################################################
# modified plot.selm() in "sn" package, for the name of xlab and ylab
plot.sn <- function(x, param.type="CP", which = c(1:4), caption,
                    panel = if (add.smooth) panel.smooth else points, main = "",
                    # sub.caption = NULL,
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
                    id.n = 3, labels.id = names(x@residuals.dp),
                    cex.id = 0.75, identline = TRUE, add.smooth = getOption("add.smooth"),
                    label.pos = c(4, 2), cex.caption = 1)
{
  if(!(is(x, "selm"))) stop("object not of class 'selm'")
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  dots <- list(...)
  nmdots <- names(dots)
  p <- slot(x, "size")["p"]
  if(missing(caption))  { caption <-  if(p> 1)
    c("Residuals vs Fitted Values",
      "Residual values and fitted error distribution",
      "Q-Q plot of (scaled DP residuals)^2",
      "P-P plot of (scaled DP residuals)^2") else
        c("Boxplot of observed values",
          "Empirical values and fitted distribution",
          "Q-Q plot of (scaled DP residuals)^2",
          "P-P plot of (scaled DP residuals)^2")}
  all.par <- slot(x, "param")
  param.type <- tolower(param.type)
  param <- all.par[[param.type]]
  if(is.null(param)) { message(paste(
    "Requested param.type='", param.type, "' evaluates to NULL.", sep=""))
    if(param.type == "pseudo-cp" & x@family== "SN")
      message("Pseudo-CP makes no sense for SN family")
    if(param.type == "cp" & x@family== "SC")
      message("CP makes no sense for SC family")
    if(param.type == "cp" & x@family== "ST")
      message("CP of ST family requires nu>4")
    stop("Consider another choice of param.type (DP or pseudo-CP)")
  }
  r <- residuals(x, param.type)
  r.lab <- paste(toupper(param.type), "residuals")
  dp <- if(length(all.par$fixed) > 0) all.par$dp.complete else all.par$dp
  nu. <- switch(x@family, ST = dp[p+3], SN = Inf, SC=1)
  rs <- slot(x,"residuals.dp")/dp[p+1]
  rs2 <- rs^2
  n <- slot(x, "size")["n.obs"]
  yh <- fitted(x, param.type)
  w <- weights(x)
  if (!is.null(w)) {
    wind <- (w != 0)
    r <- r[wind]
    yh <- yh[wind]
    w <- w[wind]
    labels.id <- labels.id[wind]
  }
  else w <- rep(1,n)
  rw <- n*w/slot(x,"size")["nw.obs"]
  cex.pts <- rw * if("cex" %in% nmdots) dots$cex else par("cex")
  if (is.null(id.n))
    id.n <- 0
  else {
    id.n <- as.integer(id.n)
    if (id.n < 0 || id.n > n)
      stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
  }
  if (id.n > 0) {
    if (is.null(labels.id))
      labels.id <- paste(1:n)
    iid <- 1:id.n
    # show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    show.rs <- sort.list(rs2, decreasing = TRUE)[iid]
    # rs2.lab <- paste("(scaled DP residuals)^2")
    text.id <- function(x, y, ind, adj.x = TRUE) {
      labpos <- if (adj.x)
        label.pos[1 + as.numeric(x > mean(range(x)))]
      else 3
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
           pos = labpos, offset = 0.25)
    }
  }
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (show[1]) {
    if(all(is.na(r)) & p>1)  message(paste("CP residuals not available;",
                                           "consider param.type='DP' or 'pseudo-CP'"))
    else {
      if(p == 1){
        y <-  (x@residuals.dp + x@fitted.values.dp)
        boxplot(y, plot=TRUE, col="gray85", border="gray60")
      }
      else { # p>1
        # if (id.n > 0)
        #    ylim <- extendrange(r = ylim, f = 0.08)
        ylim <- range(r, na.rm = TRUE)
        plot(yh, r, xlab = "Fitted values", ylab = r.lab, main = main,
             ylim = ylim, type = "n")
        panel(yh, r, ...)  # previously it included 'cex=cex.pts'
        # if (one.fig) title(sub = sub.caption, ...)
        if (id.n > 0) {
          y.id <- r[show.rs]
          y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
          text.id(yh[show.rs], y.id, show.rs)
        }
        abline(h = 0, lty = 2, col = "gray")
      } }
    mtext(caption[1], 3, 0.5, cex = cex.caption) }
  if (show[2]) {
    if(all(is.na(r)) & p>1) message(
      "CP residuals not available; consider param.type='DP' or 'pseudo-CP'")
    else {
      if (p == 1){
        y <-  (x@residuals.dp + x@fitted.values.dp)
        dp0 <- dp
        xlab="CLR-transformed count"}
      else {
        y <- r
        dp0 <- as.numeric(c(dp[1]-param[1], dp[-(1:p)]))
        xlab=r.lab
      }
      h <- hist(rep(y, w), plot=FALSE)
      extr <- extendrange(x=h$breaks)
      x.pts <- seq(max(extr), min(extr), length=501)
      d.fn <- get(paste("d", tolower(x@family), sep=""), inherits = TRUE)
      pdf <- d.fn(x.pts, dp=dp0)
      plot(c(h$mids, x.pts), c(h$density, pdf), type="n", main=main,
           xlab=xlab,  ylab="Density")
      hist(rep(y, w), col="gray95", border="gray60", probability=TRUE,
           freq=FALSE, add=TRUE)
      lines(x.pts, pdf, ...)
      rug(y, ticksize=0.02, ...)
      # if (id.n > 0) {     rug(y, ticksize=0.015, ...)
      #   text(y[show.rs], 0, labels.id[show.rs], srt=90, cex=0.5, pos=1,
      #   offset=0.2) }
      mtext(caption[2], 3, 0.25, cex = cex.caption)
    }}
  if (show[3]) {
    ylim <- c(0, max(pretty(rs2)))
    q <- qf((1:n)/(n+1), 1, nu.)
    plot(q, sort(rs2), xlab="Theoretical values", ylab="Empirical values",
         ylim=ylim, type="p", main=main, ...)   # cex=cex.pts
    if(identline) abline(0, 1, lty = 2, col = "gray50")
    # if (one.fig) title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25, cex = cex.caption)
    if (id.n > 0) text.id(q[n+1-iid], rs2[show.rs], show.rs)
  }
  if (show[4]) {
    p <- (1:n)/(n+1)
    pr <- pf(sort(rs2), 1, nu.)
    plot(p, pr, xlab="Theoretical values", ylab="Empirical values",
         xlim=c(0,1), ylim=c(0,1), main=main, ...) # cex=cex.pts,
    if(identline) abline(0, 1, lty = 2, col = "gray50")
    # if (one.fig)  title(sub = sub.caption, ...)
    mtext(caption[4], 3, 0.25, cex = cex.caption)
    if(identline) abline(0, 1, lty = 2, col = "gray50")
    if (id.n > 0)  text.id(p[n+1-iid], pr[n+1-iid], show.rs)
  }
  # if (!one.fig && par("oma")[3] >= 1)
  #     mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}



### Fig. 1
##########################################
### Simulation One, Valentim Data, T1D ###
##########################################
params <- get_params(data_d1)
N.genes <- 5000
N.samples <- 200
dat0 = create_read_numbers(params$mu,params$fit, params$p0,
                           m=N.genes, n = N.samples,
                           seed=123)
row.names(dat0) <- paste('gene', 1:N.genes, sep='')

# filter
CPM <- cpm(dat0)
keep <- (rowMeans(CPM[,1:200]) > 0.5  &  apply(dat0[,1:200], 1, function(x) length(x[x==0])/length(x)) < 0.85)
dat0 <- dat0[keep, ]; dim(dat0)

# CLR-transformation
dat0[dat0 == 0] <- 1/2
clr_simu.1 <- t(clr(t(dat0)))
clr_simu.1 <- matrix(as.numeric(clr_simu.1), nrow = dim(dat0)[1], ncol = dim(dat0)[2])


# plot clr-transformed counts of gene1 of simulated data
# Fig. 1 (a)
fit_sn1 <- selm(clr_simu.1[1, ] ~ 1, family="SN") # the first gene in simulated data
plot.sn(fit_sn1, which=2, caption = NULL, main =NULL)
title(adj=0, "(a)")

# MLE
clr.SN.fit(clr_simu.1[1, ])


# we can check any other genes in simulated data,
# almost all genes fit skew-normal model well
#for (i in 2:100) {
#  fit_sn <- selm(clr_simu.1[i, ] ~ 1, family="SN")
#  plot.sn(fit_sn, which=2,
#          caption = NULL,
#          main =paste('gene', i, sep=''))
#}



### KS Test
cp_to_dp <- function(mean=NULL, sd=NULL, skewness=NULL){
  b <- sqrt(2/pi)

  if(skewness >= 0){
    r <- (2*skewness/(4-pi))^(1/3)
  } else {
    r <- -(2*(- skewness)/(4-pi))^(1/3)
  }

  alpha <- r/sqrt(2/pi - (1-2/pi)*r^2)
  delta <- alpha/sqrt(1 + alpha^2)
  omega <- sd/sqrt(1 - (b^2)*delta^2)
  xi <- mean - b*omega*delta
  return(c(xi, omega, alpha))
} # map CP to DP


##
clrSeq_result_simu1 <- clr.SN.fit(data = clr_simu.1)
clrSeq_result_simu1 <- as.data.frame(clrSeq_result_simu1)
simu1_sn_CP <- clrSeq_result_simu1[, c("mu", "sigma", "gamma")]


simu1_sn_DP <- matrix(NA, nrow = dim(simu1_sn_CP)[1], ncol = 3)
for (i in 1:dim(simu1_sn_CP)[1]) {
  simu1_sn_DP[i,] <- c(cp_to_dp(simu1_sn_CP[i,1],
                                simu1_sn_CP[i,2],
                                simu1_sn_CP[i,3]))
}

colnames(simu1_sn_DP) <- c("xi1", "omega1", "alpha1")
simu1_sn_DP <- as.data.frame(simu1_sn_DP)


# simulation 1, KS test
simu1_KS_test_pvalue <- vector()
for (i in 1:dim(simu1_sn_DP)[1]) {
  ks <- ks.test(clr_simu.1[i, ],
                "psn",
                xi = simu1_sn_DP$xi1[i],
                omega = simu1_sn_DP$omega1[i],
                alpha = simu1_sn_DP$alpha1[i])
  simu1_KS_test_pvalue[i] <- ks$p.value
}
sum(simu1_KS_test_pvalue < 0.05)
# 39/5000 = 0.0078
# skew-normal fits 99.22% of the genes well

# Fig. S1 (b)
hist(simu1_KS_test_pvalue,
     freq = F, breaks = 20,
     xlab =expression(p),
     main = NULL, border = "grey83")
abline(v = 0.05, lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(b)")



#######################################################
### Simulation Two, Kelmer data, at 10 weeks of age ###
#######################################################
# simulate dataset
params.2 <- get_params(data_10weeks)
N.genes.2 <- 5000
N.samples.2 <- 200
dat2 = create_read_numbers(params.2$mu, params.2$fit, params.2$p0,
                           m=N.genes.2, n = N.samples.2,
                           seed=1234)
row.names(dat2) <- paste('gene', 1:N.genes.2, sep='')

# filter
CPM <- cpm(dat2)
keep <- (rowMeans(CPM[,1:N.samples.2]) > 0.5  &  apply(dat2[,1:N.samples.2], 1, function(x) length(x[x==0])/length(x)) < 0.85)
dat2 <- dat2[keep, ]; dim(dat2)


# clr-transformation
dat2[dat2 == 0] <- 1/2
clr_simu.2 <- t(clr(t(dat2)))
clr_simu.2 <- matrix(as.numeric(clr_simu.2), nrow = dim(dat2)[1], ncol = dim(dat2)[2])

# plot clr-transformed counts of gene1 of simulated data
# Fig. 1 (c)
fit_sn2 <- selm(clr_simu.2[1, ] ~ 1, family="SN")
plot.sn(fit_sn2, which=2, caption = NULL,
        main =NULL, cex.id = 0)
title(adj=0, "(c)")

# MLE
clr.SN.fit(clr_simu.2[1, ])

# we can check any other genes in simulated data,
# almost all genes fit skew-normal model well
#for (i in 2:100) {
#  fit_sn2 <- selm(clr_simu.2[i, ] ~ 1, family="SN")
#  plot.sn(fit_sn2, which=2,
#          caption = NULL,
#          main =paste('gene', i, sep=''))
#}


## KS test, simulation 2
clrSeq_result_simu2 <- clr.SN.fit(data = clr_simu.2)
clrSeq_result_simu2 <- as.data.frame(clrSeq_result_simu2)
simu2_sn_CP <- clrSeq_result_simu2[, c("mu", "sigma", "gamma")]

simu2_sn_DP <- matrix(NA, nrow = dim(simu2_sn_CP)[1], ncol = 3)
for (i in 1:dim(simu2_sn_CP)[1]) {
  simu2_sn_DP[i,] <- c(cp_to_dp(simu2_sn_CP[i,1],
                                simu2_sn_CP[i,2],
                                simu2_sn_CP[i,3]))
}

colnames(simu2_sn_DP) <- c("xi2", "omega2", "alpha2")
simu2_sn_DP <- as.data.frame(simu2_sn_DP)


# KS test
simu2_KS_test_pvalue <- vector()
for (i in 1:dim(simu2_sn_DP)[1]) {
  ks <- ks.test(clr_simu.2[i, ],
                "psn",
                xi = simu2_sn_DP$xi2[i],
                omega = simu2_sn_DP$omega2[i],
                alpha = simu2_sn_DP$alpha2[i])
  simu2_KS_test_pvalue[i] <- ks$p.value
}
sum(simu2_KS_test_pvalue < 0.05)
# 57/5000 = 0.0114
# skew-normal model fits 98.86% of the genes well

# Fig. S1 (d)
hist(simu2_KS_test_pvalue,
     freq = F, breaks = 20,
     xlab =expression(p),
     main = NULL, border = "grey83")
abline(v = 0.05, lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(b)")



##########################################
### Simulation Three, Zhou data, Liver ###
##########################################
# simulate dataset
params.3 <- get_params(Data_liver)
N.genes.3<- 5000
N.samples.3 <- 200
dat3 = create_read_numbers(params.3$mu, params.3$fit, params.3$p0,
                           m=N.genes.3, n = N.samples.3,
                           seed=12345)
row.names(dat3) <- paste('gene', 1:N.genes.3, sep='')

# filter
CPM <- cpm(dat3)
keep3 <- (rowMeans(CPM[,1:N.samples.3]) > 0.5  &  apply(dat3[,1:N.samples.3], 1, function(x) length(x[x==0])/length(x)) < 0.85)
dat3 <- dat3[keep3, ]; dim(dat3)


# clr-transformation
dat3[dat3 == 0] <- 1/2
clr_simu.3 <- t(clr(t(dat3)))
clr_simu.3 <- matrix(as.numeric(clr_simu.3), nrow = dim(dat3)[1], ncol = dim(dat3)[2])

# plot clr-transformed counts of gene1 of simulated data
# Fig. 1 (e)
fit_sn3 <- selm(clr_simu.3[1, ] ~ 1, family="SN")
plot.sn(fit_sn3, which=2, caption = NULL,
        main =NULL, cex.id = 0)
title(adj=0, "(e)")

# MLE
clr.SN.fit(clr_simu.3[1, ])

# we can check any other genes in simulated data,
# almost all genes fit skew-normal model well
#for (i in 2:100) {
#  fit_sn3 <- selm(clr_simu.3[i, ] ~ 1, family="SN")
#  plot.sn(fit_sn3, which=2,
#          caption = NULL,
#          main =paste('gene', i, sep=''))
#}


## KS test, simulation 3
clrSeq_result_simu3 <- clr.SN.fit(data = clr_simu.3)
clrSeq_result_simu3 <- as.data.frame(clrSeq_result_simu3)
simu3_sn_CP <- clrSeq_result_simu3[, c("mu", "sigma", "gamma")]


simu3_sn_DP <- matrix(NA, nrow = dim(simu3_sn_CP)[1], ncol = 3)
for (i in 1:dim(simu3_sn_CP)[1]) {
  simu3_sn_DP[i,] <- c(cp_to_dp(simu3_sn_CP[i,1],
                                simu3_sn_CP[i,2],
                                simu3_sn_CP[i,3]))
}

colnames(simu3_sn_DP) <- c("xi3", "omega3", "alpha3")
simu3_sn_DP <- as.data.frame(simu3_sn_DP)


# KS test
simu3_KS_test_pvalue <- vector()
for (i in 1:dim(simu3_sn_DP)[1]) {
  ks <- ks.test(clr_simu.3[i, ],
                "psn",
                xi = simu3_sn_DP$xi3[i],
                omega = simu3_sn_DP$omega3[i],
                alpha = simu3_sn_DP$alpha3[i])
  simu3_KS_test_pvalue[i] <- ks$p.value
}
sum(simu3_KS_test_pvalue < 0.05)
# 65/5000 = 0.013
# skew-normal model fits 98.7% of the genes well

# Fig. 1 (f)
hist(simu3_KS_test_pvalue,
     freq = F, breaks = 20,
     xlab =expression(p),
     main = NULL, border = "grey83")
abline(v = 0.05, lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(f)")
###END###