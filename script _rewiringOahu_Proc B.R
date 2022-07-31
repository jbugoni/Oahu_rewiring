#Analyses for the manuscript testing rewiring & drivers of interaction frequencies on the Oahu seed dispersal network
#Last update: 2022 July 31

### Required pacakges
require(bipartite)
require(plotrix) # for plotting

#Before start, load functions available at: github.com/vanderleidebastiani/rewiring
#Such functions are needed to allow rewiring following species extinctions:

# Small change in the internal function "extinction" of the bipartite package
# The difference is that it returns rexcl and cexcl (default NULL to both) with complete row (rexcl) or column (cexcl) extinct step by step
#source("./Functions/extinction.mod.R")
extinction.mod <- function (web, participant = "both", method = "random", ext.row = NULL, ext.col = NULL) 
{
  RES <- list(rexcl = NULL, cexcl = NULL)
  partis <- c("lower", "higher", "both")
  partis.match <- pmatch(participant, partis)
  if (is.na(partis.match)) 
    stop("Choose participant: lower/higher/both.\n")
  meths <- c("random", "abundance", "degree", "external")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) 
    stop("Choose extinction method: random/abundance/degree.\n")
  nr <- NROW(web)
  nc <- NCOL(web)
  if (partis.match == 3 & meths.match == 1) 
    partis.match <- sample(2, 1)
  if (meths.match == 1) {
    rexcl <- sample(nr, 1)
    cexcl <- sample(nc, 1)
    if (partis.match == 1){
      RES$rexcl <- web[rexcl, , drop = FALSE]
      web[rexcl, ] <- 0
    }
    if (partis.match == 2) {
      RES$cexcl <- web[, cexcl, drop = FALSE]
      web[, cexcl] <- 0
    }
  }
  if (meths.match == 2) {
    web <- web[sample(1:nrow(web)), sample(1:ncol(web)), drop = FALSE]
    rseq <- order(rowSums(web))
    cseq <- order(colSums(web))
    if (partis.match == 1){
      RES$rexcl <- web[rseq[1], , drop = FALSE]
      web[rseq[1], ] <- 0
    }
    if (partis.match == 2){ 
      RES$cexcl <- web[, cseq[1], drop = FALSE]
      web[, cseq[1]] <- 0
    }
    if (partis.match == 3) {
      if (min(rowSums(web)) < min(colSums(web))) {
        RES$rexcl <- web[rseq[1], , drop = FALSE]
        web[rseq[1], ] <- 0
      }
      else {
        if (min(rowSums(web)) > min(colSums(web))) {
          RES$cexcl <- web[, cseq[1], drop = FALSE]
          web[, cseq[1]] <- 0
        }
        else {
          if (sample(2, 1) == 1) {
            RES$rexcl <- web[rseq[1], , drop = FALSE]
            web[rseq[1], ] <- 0
          }
          else {
            RES$cexcl <- web[, cseq[1], drop = FALSE]
            web[, cseq[1]] <- 0
          }
        }
      }
    }
  }
  if (meths.match == 3) {
    if (partis.match == 1) {
      sequ <- rowSums(web > 0)
      which.ex <- which(sequ == max(sequ))
      if (length(which.ex) > 1) {
        ex <- sample(which.ex, size = 1)
      }
      else {
        ex <- which.ex
      }
      RES$rexcl <- web[ex, , drop = FALSE]
      web[ex, ] <- 0
    }
    if (partis.match == 2) {
      sequ <- colSums(web > 0)
      which.ex <- which(sequ == max(sequ))
      if (length(which.ex) > 1) 
        ex <- sample(which.ex, size = 1)
      else ex <- which.ex
      RES$cexcl <- web[, ex, drop = FALSE]
      web[, ex] <- 0
    }
  }
  if (meths.match == 4) {
    rseq <- ext.row
    cseq <- ext.col
    if (partis.match == 1){ 
      RES$rexcl <- web[rseq[1], , drop = FALSE]
      web[rseq[1], ] <- 0
    }
    if (partis.match == 2) {
      RES$cexcl <- web[, cseq[1], drop = FALSE]
      web[, cseq[1]] <- 0
    }
  }
  RES$web <- web
  return(RES)
}

#source("./Functions/one.second.extinct.mod.R")
# Small change in the internal function "one.second.extinct" of the bipartite package. The change allows rewiring of interations after each step of extinction
# The arguments are the same of "one.second.extinct" function and additional argumets are described below.
# New arguments:
# rewiring - Logical argument to specify if allow rewiring (default rewiring = FALSE).
# probabilities.rewiring1 - A matrix with probabilities of rewiring, must be the same dimensions of the web (i.e. network). See section Methods in Vizentin-Bugoni et al. [2020] Meth. Ecol. Evol. for details). This matrix is required in step ii of framework (default probabilities.rewiring1 = NULL).
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods in Vizentin-Bugoni et al. [2020] Meth. Ecol. Evol. for details). This matrix is required in step iii of framework (default probabilities.rewiring2 = NULL).
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner" and "multiple.trials.each.partner". See section Methods in Vizentin-Bugoni et al. [2020] Meth. Ecol. Evol. for details); (default method.rewiring = "one.try.single.partner").
one.second.extinct.mod <- function(web, participant = "higher", method = "abun", ext.row = NULL, ext.col = NULL, 
                                   rewiring = FALSE, probabilities.rewiring1 = NULL, probabilities.rewiring2 = NULL,
                                   method.rewiring = "one.try.single.partner") {
  dead <- matrix(nrow = 0, ncol = 3)
  colnames(dead) <- c("no", "ext.lower", "ext.higher")
  m2 <- web
  i <- 1
  METHOD.REWIRING = c("one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner", "multiple.trials.each.partner")
  method.rewiring <- pmatch(method.rewiring, METHOD.REWIRING)
  if (length(method.rewiring) > 1) {
    stop("\n Only one argument is accepted in method.rewiring \n")
  }
  if (is.na(method.rewiring)) {
    stop("\n Invalid method.rewiring \n")
  }
  if(method.rewiring == 4 | method.rewiring == 5){
    keep.trying <- TRUE
  } else {
    keep.trying <- FALSE
  }
  method.rewiring <- ifelse(method.rewiring == 4, 1, ifelse(method.rewiring == 5, 2, method.rewiring))
  if(rewiring){
    if(any(web%%1!=0)){
      stop("\n If rewiring is TRUE the web must must contain only integers \n")
    }
    if(is.null(rownames(web)) | is.null(colnames(web))){
      stop("\n If rewiring is TRUE the web must must rownames and colnames\n")
    }
    if(is.null(probabilities.rewiring1) | is.null(probabilities.rewiring2)){
      stop("\n If rewiring is TRUE probabilities.rewiring1 and probabilities.rewiring1 must not be NULL\n")
    }
  }
  repeat {
    ext.temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
    if(rewiring){
      if(!is.null(ext.temp$rexcl)){
        sp.ext <- rownames(ext.temp$rexcl)
        sp.try.rewiring <- which(ext.temp$rexcl>0)
        sp.surv <- seq_len(nrow(ext.temp$web))
        sp.surv <- sp.surv[-1*which(rownames(ext.temp$web) %in% sp.ext)]
        for(jj in sp.try.rewiring){
          sp.surv.temp <- sp.surv
          go <- TRUE
          m <- 0
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$rexcl[1, jj]
          }
          while (go) {
            m <- m+1
            sp.surv.prob1 <- probabilities.rewiring1[sp.surv.temp, jj]
            sp.add <- sample(as.character(sp.surv.temp), 1, prob = sp.surv.prob1)
            sp.surv.prob2 <- probabilities.rewiring2[as.numeric(sp.add), jj]
            n.add <- rbinom(1, trials, sp.surv.prob2)
            if(n.add>0){
              ext.temp$web[as.numeric(sp.add), jj] <- ext.temp$web[as.numeric(sp.add), jj]+n.add
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              if(!keep.trying){
                go <- FALSE  
              } else{
                if((method.rewiring == 1 & n.add>0) | (method.rewiring == 2 & n.add == trials)){
                  go <- FALSE
                } else{
                  sp.surv.temp <- sp.surv.temp[-1*which(sp.surv.temp%in% sp.add)]
                  trials <- trials-n.add
                  if(length(sp.surv.temp)<1){
                    go <- FALSE
                  }
                }
              }
            } else {
              if(ext.temp$rexcl[1, jj] == m){
                go <- FALSE
              }
            }
          }
        }
      }
      if(!is.null(ext.temp$cexcl)){
        sp.ext <- colnames(ext.temp$cexcl)
        sp.try.rewiring <- which(ext.temp$cexcl>0)
        sp.surv <- seq_len(ncol(ext.temp$web))
        sp.surv <- sp.surv[-1*which(colnames(ext.temp$web) %in% sp.ext)]
        for(ii in sp.try.rewiring){
          sp.surv.temp <- sp.surv
          go <- TRUE
          m <- 0
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$cexcl[ii, 1]
          }
          while (go) {
            m <- m+1
            sp.surv.prob1 <- probabilities.rewiring1[ii, sp.surv.temp]
            sp.add <- sample(as.character(sp.surv.temp), 1, prob = sp.surv.prob1)
            sp.surv.prob2 <- probabilities.rewiring2[ii, as.numeric(sp.add)]
            n.add <- rbinom(1, trials, sp.surv.prob2)
            if(n.add>0){
              ext.temp$web[ii, as.numeric(sp.add)] <- ext.temp$web[ii, as.numeric(sp.add)]+n.add
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              if(!keep.trying){
                go <- FALSE  
              } else{
                if((method.rewiring == 1 & n.add>0) | (method.rewiring == 2 & n.add == trials)){
                  go <- FALSE
                } else{
                  sp.surv.temp <- sp.surv.temp[-1*which(sp.surv.temp%in% sp.add)]
                  trials <- trials-n.add
                  if(length(sp.surv.temp)<1){
                    go <- FALSE
                  }
                }
              }
            } else {
              if(ext.temp$cexcl[ii, 1] == m){
                go <- FALSE
              }
            }
          }
        }
      }
    }
    n <- ext.temp$web
    dead <- rbind(dead, c(i, attributes(m2 <- empty(n, count = TRUE))$empty))
    if (participant == "lower" & NROW(m2) < 2) 
      break
    if (participant == "higher" & NCOL(m2) < 2) 
      break
    if (participant == "both" & min(dim(m2)) < 2) 
      break
    if (any(dim(n) == 1)) 
      break
    if (method == "external") {
      ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > ext.col[1]] - 1
      ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > ext.row[1]] - 1
      ext.row <- ext.row[-1]
      ext.col <- ext.col[-1]
    }
    i <- i + 1
  }
  dead2 <- rbind(dead, c(NROW(dead) + 1, NROW(m2), NCOL(m2)))
  if (participant == "lower" & method == "degree") {
    if (length(table(dead[, 2])) > 1) 
      dead2[, 2] <- 1
  }
  if (nrow(dead) + 1 != nrow(dead2)) 
    stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
  if (participant == "lower") 
    supposed.length <- NROW(web)
  if (participant == "higher") 
    supposed.length <- NCOL(web)
  if (participant == "both") 
    supposed.length <- NROW(dead2)
  if (NROW(dead2) != supposed.length) {
    missing <- supposed.length - NROW(dead2)
    addit1 <- (NROW(dead2) + 1):(NROW(dead2) + missing)
    addit2n3 <- rep(0, times = missing)
    dead2 <- rbind(dead2, as.matrix(data.frame(addit1, addit2n3, addit2n3)))
  }
  out <- dead2
  class(out) <- "bipartite"
  attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
  attr(out, "exterminated")
  return(out)
}


# This way robustness function can be applied to a list of results returned by the function "one.second.extinct.mod"
#source("./Functions/calc.mean.one.second.extinct.mod.R")
# Small code extracted of the "second.extinct" function to calculate mean of multiple replicates of extinction sequence and set specific attributes to object returned
# Code extracted from the "second.extinct" function to calcultate mean of multiple replicates of extinction sequence and set specific attributes to the object returned

calc.mean.one.second.extinct.mod <- function(o){
  lengths <- sapply(o, nrow)
  z <- o[[which.max(lengths)]]
  @@ -15,4 +15,4 @@ calc.mean.one.second.extinct.mod <- function(o){
    class(out) <- "bipartite"
    attr(out, "exterminated") <- attr(o[[1]],"exterminated")
    return(out)
  } 
}

#source("./Functions/matrix.p1.R") #not used here

# Small function to facilitate the calculation of confidence interval (default conf.level = 0.95)
# Function to facilitate the calculation of confidence interval (default conf.level = 0.95)
IC <- function(X, conf.level = 0.95){
  m <- mean(X)
  s <- sd(X)
  n <- length(X)
  se <- s/sqrt(n)
  error <- qt((1-conf.level)/2, df = n-1, lower.tail = FALSE)*se
  res <- c(mean = m, sd = s, lower = m-error, upper = m+error)
  return(res)
}}

#######################################################################                                           

#Load data
#Full dataset availabe at: VINE, Hawaii, 2022, "Vizentin-Bugoni, J., J. H. Sperry, J. P. Kelley, J. T. Foster, D. R. Drake, S. B. Case, J. M. Gleditsch, A. M. Hruska, R. C. Wilcox, C. E. Tarwater. "Dataset from: Mechanisms underlying interaction frequencies and robustness in a novel seed dispersal network"", https://doi.org/10.7910/DVN/IDQQRC, Harvard Dataverse, V2, UNF:6:OewPamEPdSiLAwFdM4omYQ== [fileUNF]

#load data
net=read.table("net.txt", header=TRUE) # Load the weighted interaction network
net1 <- net[-c(42:43), ]# two plant morphotypes (UNK_101 and UNK_74) removed
net1=as.matrix(net1) #network converted to matrix

bird_abund=read.table("bird_abund.txt", header=TRUE) #bird abundances
plant_abund=read.table("plant_abund.txt", header=TRUE) #plant (fruit) abundances
bird_gape=read.table("bird_gape.txt", header=TRUE) #bird gape width
plant_seedsize=read.table("plant_seedsize.txt", header=TRUE) #seed width
bird_DS=read.table("bird_dist_spatial.txt", h=T) #bird spatial distribution
plant_DS=read.table("plant_dist_spatial.txt", h=T)#plant spatial distribution
bird_DT=read.table("bird_temp_dist.txt", h=T) #bird temporal distribution
plant_DT=read.table("plant_temp_dist.txt", h=T) #fruit temporal distribution

bird_null=read.table("bird_null.txt", header=TRUE) #A vector necessary later 
plant_null=read.table("plant_null.txt", header=TRUE) #A vector necessary later 

#Building probablistic matrices that will define pairwise probabilities of rewiring

#Relative abundances of pairwise bird and plant species
bird_abund1=bird_abund/sum(bird_abund) #bird relative abundances
bird_abund1=as.matrix(bird_abund1)# transform to matrix

plant_abund1=plant_abund/sum(plant_abund) #plant relative abundances
plant_abund1=as.matrix(plant_abund1)
plant_abund2 <- plant_abund1[-c(42:43), ]# UNK_101 and UNK_74 removed

# pairwise combination of bird and plant relative abundances
abundance <- plant_abund2%*%t(bird_abund1)
abundance                                              

# Relative abundances of plants
abundance_pl <- sweep(abundance, 2, colSums(abundance), "/")
abundance_pl

# Relative abundances of birds
abundance_b <- sweep(abundance, 1, rowSums(abundance), "/") 
abundance_b

# Morphological matching
bird_gape=as.matrix(bird_gape)
plant_seedsize = as.matrix(plant_seedsize)
plant_seedsize1 <- plant_seedsize[-c(42:43), ]# UNK_101 and UNK_74 removed
plant_seedsize1 = as.matrix(plant_seedsize1)
morphological <- 1-(as.matrix(vegdist(rbind(bird_gape,plant_seedsize1), method = "gower"))[(length(bird_gape)+1):(length(bird_gape)+length(plant_seedsize1)),(1:length(plant_seedsize1))])

# spatial overlap (Relative spatial coexistence)
plant_DS1 <- plant_DS[-c(42:43), ]# UNK_101 and UNK_74 removed
spatial <- plant_DS1%*%t(bird_DS)
spatial1 <- spatial/max(spatial) 

#temporal (phenological) overlap
bird_DT=as.matrix(bird_DT)
plant_DT=as.matrix(plant_DT)
plant_DT1 <- plant_DT[-c(42:43), ]# UNK_101 and UNK_74 removed
pheno_overlap <- plant_DT1%*%t(bird_DT)
pheno_overlap1 <- pheno_overlap/max(pheno_overlap)

## NUll models where all interactions have the same probability
one <- net1
one[] <- 1
one

############# Models combining mechanisms 

MS <- morphological1*spatial1
MP <- morphological1*pheno_overlap1
SP <- spatial1*pheno_overlap1
AM <- abundance*morphological1
AS <- abundance*spatial1
AP <- abundance*pheno_overlap1
MSP <- morphological1*spatial1*pheno_overlap1
AMS <- abundance*morphological1*spatial1
APS <- abundance*pheno_overlap1*spatial1
AMP <- abundance*morphological1*pheno_overlap1
AMSP <- abundance*morphological1*spatial1*pheno_overlap1

#all matrices must sum 1, so:
MS=MS/sum(MS)
MP=MP/sum(MP)
SP=SP/sum(SP)
MSP=MSP/sum(MSP)
AM=AM/sum(AM)
AS=AS/sum(AS)
AP=AP/sum(AP)
AMS=AMS/sum(AMS)
APS=APS/sum(APS)
AMP=AMP/sum(AMP)
AMSP=AMSP/sum(AMSP)

### Simulate secondary extinctions in networks 

## Define the number of replications
nrep <- 1000

## Define the participant to be extinct
#Options: "lower" or "higher"
participant <- "lower"
participant

## Define method to remove species
#Options: "random" or "degree"
method <- "degree"
method

## Define method of rewiring
#Options: "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner", "multiple.trials.each.partner"
method.rewiring <- "one.try.single.partner"
method.rewiring

## Define the matrices of probability based on which rewiring will take place
# To plants loss ("lower" participants)
probabilities.rewiring1 <- abundance_pl  # plant relative abundances # For random choice in this step, use "plant_null" or "bird_null"
probabilities.rewiring1 # Probability of choice of a potential partner

## Run secondary extinctions with specified parameters (results are in a list)
# The argument probabilities.rewiring2 to specify the probabilities of rewiring 
RES.without.rewiring <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = FALSE), simplify = FALSE) #this model is the beckmark to calculate effect size latter.

RES.with.rewiring.M <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = morphological1, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.P <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = pheno_overlap1, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.S <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = spatial1, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.MS <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = MS, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.MP <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = MF, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.SP <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = SF, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.MSP <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = MSF, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.A <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = A, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.AM <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AM, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.AS <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AS, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.AP <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AF, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.AMS <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AMS, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.APS <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AFS, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.AMP <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AMF, method.rewiring = method.rewiring), simplify = FALSE)

RES.with.rewiring.AMSP <- replicate(nrep, one.second.extinct.mod(net1, participant = participant, method = method, rewiring = TRUE, probabilities.rewiring1 = probabilities.rewiring1, probabilities.rewiring2 = AMSF, method.rewiring = method.rewiring), simplify = FALSE)

## Calculate robustness to species extinctions
RES.robustness.without.rewiring <- sapply(RES.without.rewiring, robustness)
RES.robustness.with.rewiring.M <- sapply(RES.with.rewiring.M, robustness)
RES.robustness.with.rewiring.P <- sapply(RES.with.rewiring.P, robustness)
RES.robustness.with.rewiring.S <- sapply(RES.with.rewiring.S, robustness)
RES.robustness.with.rewiring.MS <- sapply(RES.with.rewiring.MS, robustness)
RES.robustness.with.rewiring.MP <- sapply(RES.with.rewiring.MP, robustness)
RES.robustness.with.rewiring.SP <- sapply(RES.with.rewiring.SP, robustness)
RES.robustness.with.rewiring.MSP <- sapply(RES.with.rewiring.MSP, robustness)

RES.robustness.with.rewiring.A <- sapply(RES.with.rewiring.A, robustness)
RES.robustness.with.rewiring.AM <- sapply(RES.with.rewiring.AM, robustness)
RES.robustness.with.rewiring.AS <- sapply(RES.with.rewiring.AS, robustness)
RES.robustness.with.rewiring.AP <- sapply(RES.with.rewiring.AP, robustness)
RES.robustness.with.rewiring.AMS <- sapply(RES.with.rewiring.AMS, robustness)
RES.robustness.with.rewiring.AFS <- sapply(RES.with.rewiring.APS, robustness)
RES.robustness.with.rewiring.AMP <- sapply(RES.with.rewiring.AMP, robustness)
RES.robustness.with.rewiring.AMSP <- sapply(RES.with.rewiring.AMSP, robustness)

# Organize the results
res.robustness <- data.frame(robustness.without.rewiring = RES.robustness.without.rewiring, 
                             robustness.with.rewiring.M = RES.robustness.with.rewiring.M,
                             robustness.with.rewiring.P = RES.robustness.with.rewiring.P,
                             robustness.with.rewiring.S = RES.robustness.with.rewiring.S,
                             robustness.with.rewiring.MS = RES.robustness.with.rewiring.MS,
                             robustness.with.rewiring.MP = RES.robustness.with.rewiring.MP,
                             robustness.with.rewiring.SP = RES.robustness.with.rewiring.SP,
                             robustness.with.rewiring.MSP = RES.robustness.with.rewiring.MSP)#), # If abundance is already in probabilities.rewiring1 <- abundance_pl" or "bird_abund", then multi-mechanism models (including Abundances) are not considered. Otherwise, include the remaining models (last part, after #)

robustness.with.rewiring.A = RES.robustness.with.rewiring.A,
robustness.with.rewiring.AM = RES.robustness.with.rewiring.AM,
robustness.with.rewiring.AS = RES.robustness.with.rewiring.AS,
robustness.with.rewiring.AP = RES.robustness.with.rewiring.AP,
robustness.with.rewiring.AMS = RES.robustness.with.rewiring.AMS,
robustness.with.rewiring.APS = RES.robustness.with.rewiring.APS,
robustness.with.rewiring.AMP = RES.robustness.with.rewiring.AMP,
robustness.with.rewiring.AMSP = RES.robustness.with.rewiring.AMSP)


# Compute 95% confidence intervals
res.robustness.summary <- as.data.frame(t(sapply(res.robustness, IC)))
res.robustness.summary$rewiring <- c("0", "M", "P", "S", "MS", "MP", "SP", "MSP")#, "A", "AM", "AS","AP","AMS","APS", "AMP","AMSP") # If abundance is already in probabilities.rewiring1 <- abundance_pl" or "bird_abund", then multi-mechanism models are not considered. Otherwise, include the remaining models (last part, after #)
res.robustness.summary

# Plot results for quality control
plot(NA, xlim = c((1-0.5),(nrow(res.robustness.summary)+0.5)), ylim = c(0.0,1), 
     type = "n", las = 1, xaxt ="n", 
     ylab = "Robustness", xlab ="Mode of rewiring",
     main = paste("Participant =", participant, ";", "Method =", method))

axis(1, at = 1:nrow(res.robustness.summary), labels = res.robustness.summary$rewiring)

plotCI(1:nrow(res.robustness.summary), res.robustness.summary$mean, 
       ui = res.robustness.summary$upper, 
       li = res.robustness.summary$lower, 
       add = TRUE, pch = 19, cex = 0.7) # mean and confidence intervals
abline(h = c(res.robustness.summary$lower[1], res.robustness.summary$upper[1]), lty = 3) # lines of confidence interval to scenario without rewiring

#Now to obtain D (effect sizes), subtract mean robustness obtained in each scenario - robustness obtained without rewiring ("RES.robustness.without.rewiring") 
#Final figures using these results were produced outside R

#Repeat the same procedure for all extinction scenarios. 
#For all scenarios, see Figure 3 and the Methods section in see Vizentin-Bugoni et al. (in review - Proc B)





#END OF THE MODELING for REWIRING
###########################
##########################

#MODELING Drivers of interaction frequencies, following Vazquez et al. (2009 Ecology)

#creating equiprobable model
null=one/sum(one) # all interactions receive the same probability

#normalizing single-mechanims matrices
abundance2=abundance/sum(abundance)
morphological2=morphological1/sum(morphological1)
spatial2=spatial1/sum(spatial1)
pheno_overlap2=pheno_overlap1/sum(pheno_overlap1)

#Recombining matrices to have models including more that one mechanism
MS2<- morphological2*spatial2
MP2 <- morphological2*pheno_overlap2
SP2 <- spatial2*pheno_overlap2
AM2 <- abundance2*morphological2
AS2 <- abundance2*spatial2
AP2 <- abundance2*pheno_overlap2
MSP2 <- morphological2*spatial2*pheno_overlap2
AMS2 <- abundance2*morphological2*spatial2
APS2 <- abundance2*pheno_overlap2*spatial2
AMP2 <- abundance2*morphological2*pheno_overlap2
AMSP2 <- abundance2*morphological2*spatial2*pheno_overlap2

#all matrices must sum 1
MS2=MS2/sum(MS2)
MP2=MP2/sum(MP2)
SP2=SP2/sum(SP2)
AM2=AM2/sum(AM2)
AS2=AS2/sum(AS2)
AP2=AP2/sum(AP2)
MSP2=MSP2/sum(MSP2)
AMS2=AMS2/sum(AMS2)
APS2=APS2/sum(APS2)
AMP2=AMP2/sum(AMP2)
AMSP2=AMSP2/sum(AMSP2)

#loading mlik function from Vazquez et al. https://figshare.com/articles/dataset/Supplement_1_R_functions_used_for_the_analysis_/3531572?file=5604407
mlik<-function(o,m.p,par){
  lik<--dmultinom(o,prob=m.p,log=T) #Negative log likelihood
  aic=2*lik+2*par
  res=data.frame(lik,aic)
  return(res)
}

#likelihood analyses
aic0=mlik(net1, null, par=0) #equiprobable interactions
aicA=mlik(net1, abundance2, par=1)
aicM=mlik(net1, morphological2, par=1)
aicS=mlik(net1, spatial2, par=1)
aicP=mlik(net1, pheno_overlap2, par=1)
aicMS=mlik(net1, MS2, par=2)
aicMP=mlik(net1, MP2, par=2)
aicSP=mlik(net1, SP2, par=2) # model with the lowest AIC 
aicMSP=mlik(net1, MSP2, par=3)
aicAM=mlik(net1, AM2, par=2)
aicAS=mlik(net1, AS2, par=2)
aicAP=mlik(net1, AP2, par=2)
aicAMS=mlik(net1, AMS2, par=3)
aicAPS=mlik(net1, APS2, par=3)
aicAMP=mlik(net1, AMP2, par=3)
aicAMSP=mlik(net1, AMSP2, par=4)

#Results summarized on the Table 1 of the paper.
#END
