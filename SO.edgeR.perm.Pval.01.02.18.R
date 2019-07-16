### README for future ###
# This script computes eFDR using null p values.
# First this scripts computes observed p-alues 
### END README ###

rm(list = ls())
stopifnot(require("edgeR"))

#### Parameters
Brain_region = "T"
dir = ".tagwise"
####

setwd ("/Volumes/GoogleDrive/My Drive/social_oppurtunity/challenge.opp.TRN/counts/")

# reading targets file
targets = readTargets("sb_mrsb_opportunity_targets.txt")
# selecting brain region
targets = targets[which(targets$region == Brain_region),]

#function to return edgeR DGE object after filtering low count genes. 
getDGE = function(targets)
{
  Raw_DGE = readDGE(targets)
  #Grabing indexes of last summary lines from HTSEQ-Count
  MetaTags <- grep("^__", rownames(Raw_DGE))
  
  Raw_DGE = Raw_DGE[-MetaTags,]
  
  keep <- rowSums(cpm(Raw_DGE)>0.5) >= 5
  print(table(keep))
  Filtered_DGE = Raw_DGE[keep,]
  Filtered_DGE = calcNormFactors(Filtered_DGE)
}

DGE = getDGE(targets)
#Making a group variable for pairwise comparisons
targets$group = paste(targets$nest_condition, targets$timepoint, sep = ".")
group = as.factor(targets$group)
group = relevel(group, ref = "NOCLUTCH.9am")

design_Group = model.matrix(~group)


get_fit = function(y, design)
{
  y = estimateGLMCommonDisp(y,design, verbose=T)
  y = estimateGLMTrendedDisp(y,design, verbose=T)
  y = estimateGLMTagwiseDisp(y,design)
  fit = glmFit(y, design)
  return(fit)
}
#Using glmFIT to find negative binomial fit for each gene.
fit.Group = get_fit(DGE, design_Group)

# using LRT to rank gene for each following pairwise contrast.
Nest_o = as.data.frame(topTags(glmLRT(fit.Group, coef=6), n=Inf, sort.by="PValue"))
Eggs_o = as.data.frame(topTags(glmLRT(fit.Group, coef=2), n=Inf, sort.by="PValue"))
HEarly_o = as.data.frame(topTags(glmLRT(fit.Group, coef=5), n=Inf, sort.by="PValue"))
HMid_o = as.data.frame(topTags(glmLRT(fit.Group, contrast=c(0,0,1,0,0,0,-1,0)), n=Inf, sort.by="PValue"))
HLate_o = as.data.frame(topTags(glmLRT(fit.Group, contrast=c(0,0,0,1,0,0,0,-1)), n=Inf, sort.by="PValue"))

#setting up a writing directory for results.
setwd("/home/n-z/sbukhar/eFDR/so_perm/")

# reading most recent null distribution files. Please specify most recent null distribution files here.
Nest_null = read.delim(paste(Brain_region,"_nullNest2017-12-28.txt", sep = ""), header = T, sep = "\t")
Eggs_null = read.delim(paste(Brain_region,"_nullEggs2017-12-28.txt", sep = ""), header = T, sep = "\t")
HEarly_null = read.delim(paste(Brain_region,"_nullHEarly2017-12-28.txt", sep = ""), header = T, sep = "\t")
HMid_null = read.delim(paste(Brain_region,"_nullHMid2017-12-28.txt", sep = ""), header = T, sep = "\t")
HLate_null = read.delim(paste(Brain_region,"_nullHLate2017-12-28.txt", sep = ""), header = T, sep = "\t")

#eFDR computation function. This function takes two arguments (1) null distribution as a matrix (2) edgeR LRT contrast dataframe which has $PValue column
# This function uses parallel computing, so make sure library(doParallel) is already installed.
get_eFDR = function(c30_null, c30_o)
{
  rownames(c30_null) = c30_null$names
  c30_null = c30_null[,c(-2)]
  R = dim(c30_null)[2]
  qval = c()
  
  qval <- foreach(i = 1:dim(c30_o)[1], .combine=c) %dopar% {
    # how many p-value are lower among the random matrix
    length(which(c30_null <= c30_o$PValue[i])) / R
  }
  #to avoide getting a 0 q value adding +1 on both sides.
  qval <- (qval+1)/(length(qval)+1)
  return(qval)
}
# computing eFDR for each contrast using appropriate null distributions.
Nest_o$eFDR = get_eFDR (Nest_null, Nest_o)
Eggs_o$eFDR = get_eFDR (Eggs_null, Eggs_o)
HEarly_o$eFDR = get_eFDR (HEarly_null, HEarly_o)
HMid_o$eFDR = get_eFDR (HMid_null, HMid_o)
HLate_o$eFDR = get_eFDR (HLate_null, HLate_o)

# writing results. 
write.table(Nest_o, file = paste(Brain_region, "con_Nest_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
write.table(Eggs_o, file = paste(Brain_region, "con_Eggs_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
write.table(HEarly_o, file = paste(Brain_region, "con_HEarly_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
write.table(HMid_o, file = paste(Brain_region, "con_HMid_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
write.table(HLate_o, file = paste(Brain_region, "con_HLate_eFDR", Sys.Date(),"txt", sep = "."), quote = F, row.names = F, sep = "\t")
