rm (list =ls())
setwd ("/home/n-z/sbukhar/eFDR/sc_counts/")

stopifnot(require("edgeR"))
stopifnot(require("plyr"))

#### Parameters
Brain_region = "D"
####

targets = readTargets("Targets.txt")
targets = targets[which(targets$BrainRegion == Brain_region),]

get_DGE = function(targets)
{
  Raw_DGE = readDGE(targets)
  MetaTags <- grep("^__", rownames(Raw_DGE))
  Raw_DGE = Raw_DGE[-MetaTags,]# removing meta tags
  keep <- rowSums(cpm(Raw_DGE)>1) >= 2
  print(table(keep))
  Filtered_DGE = Raw_DGE[keep,]
  Filtered_DGE = calcNormFactors(Filtered_DGE)
  return(Filtered_DGE)
}

DGE_all = get_DGE(targets)

time = as.factor(targets$TimePoint)
intruder = as.factor(targets$Treatment)
intruder = relevel(intruder, ref="Control")
group = paste(time, intruder, sep = ".")

get_fit = function(y, design)
{
  y = estimateGLMCommonDisp(y,design)
  y = estimateGLMTrendedDisp(y,design)
  y = estimateGLMTagwiseDisp(y,design)
  fit = glmFit (y, design)
  return(fit)
}
design.pairwise = model.matrix(~time+time:intruder)
fit = get_fit(DGE_all,design.pairwise)

library(doMC)
ncore = parallel::detectCores()
registerDoMC(cores = ncore)


c30_o = as.data.frame(topTags(glmLRT(fit, coef=4), n=Inf, sort.by="PValue"))
c30_o$names = rownames(c30_o)
c60_o = as.data.frame(topTags(glmLRT(fit, coef=5), n=Inf, sort.by="PValue"))
c120_o = as.data.frame(topTags(glmLRT(fit, coef=6), n=Inf, sort.by="PValue"))

setwd("/home/n-z/sbukhar/eFDR/sc_perm/")
c30_null = read.delim(paste(Brain_region,"_null30_2017-10-20.txt", sep = ""), header = T, sep = "\t")
c60_null = read.delim(paste(Brain_region,"_null60_2017-10-20.txt", sep = ""), header = T, sep = "\t")
c120_null = read.delim(paste(Brain_region,"_null120_2017-10-20.txt", sep = ""), header = T, sep = "\t")


get_eFDR = function(c30_null, c30_o)
{
  rownames(c30_null) = c30_null$names
  c30_null = c30_null[,c(-2)]
  R = dim(c30_null)[2]
  qval = c()
  qval <- foreach(i = 1:dim(c30_o)[1], .combine=c) %dopar% {
    # how many p-value are lower among the random matrix (lower is better)
    length(which(c30_null <= c30_o$PValue[i])) / R
  }
  #to avoide getting a q value as +1 in both sides.
  qval <- (qval+1)/(length(qval)+1)
  return(qval)
}

c30_o$eFDR = get_eFDR (c30_null, c30_o)
c60_o$eFDR = get_eFDR (c60_null, c60_o)
c120_o$eFDR = get_eFDR (c120_null, c120_o)

write.table(c30_o, file = paste(Brain_region, "con_30_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
write.table(c60_o, file = paste(Brain_region, "con_60_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
write.table(c120_o, file = paste(Brain_region, "con_120_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
