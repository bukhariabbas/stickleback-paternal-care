rm (list =ls())
setwd ("/Users/syedabbasbukhari/Google Drive/social_oppurtunity/challenge.opp.TRN/")

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

#DGE_all = get_DGE(targets)

time = as.factor(targets$TimePoint)
intruder = as.factor(targets$Treatment)
intruder = relevel(intruder, ref="Control")
group = paste(time, intruder, sep = ".")

get_fit = function(y, design)
{
  y = estimateGLMCommonDisp(y,design)
  y = estimateGLMTrendedDisp(y,design)
  y = estimateGLMTagwiseDisp(y,design)
  fit = glmFit(y, design)
  return(fit)
}
design.pairwise = model.matrix(~time+time:intruder)

perm_Effect30 = c()
perm_Effect60 = c()
perm_Effect120 = c()

library(plyr)
totalperm = 500
for(i in 1:totalperm)
{
  t <- as.numeric(Sys.time())
  seed <- 1e8 * (t - floor(t))
  set.seed(seed)
  targets_perm = targets
  DGE_perm = get_DGE(targets_perm)
  DGE_perm$samples = DGE_perm$samples[sample(rownames(DGE_perm$samples)),]
  perm.fit = get_fit(DGE_perm, design.pairwise)
  c30 = as.data.frame(topTags(glmLRT(perm.fit, coef=4), n=Inf, sort.by="PValue"))
  c60 = as.data.frame(topTags(glmLRT(perm.fit, coef=5), n=Inf, sort.by="PValue"))
  c120 = as.data.frame(topTags(glmLRT(perm.fit, coef=6), n=Inf, sort.by="PValue"))
  
  if(i == 1)
  {
    perm_Effect30$iter1 = c30$PValue
    perm_Effect30$names = rownames(c30)
    perm_Effect30 = as.data.frame(perm_Effect30)
    
    perm_Effect60$iter1 = c60$PValue
    perm_Effect60$names = rownames(c60)
    perm_Effect60 = as.data.frame(perm_Effect60)
    
    perm_Effect120$iter1 = c120$PValue
    perm_Effect120$names = rownames(c120)
    perm_Effect120 = as.data.frame(perm_Effect120)
  }
  else
  {
    c1 = c()
    c1$iter = c30$PValue
    c1 = as.data.frame(c1)
    c1$names = rownames(c30)
    colnames(c1)= c(paste("iter",i, sep = "."),"names")
    perm_Effect30 = plyr::join(perm_Effect30, c1)
    
    c2 = c()
    c2$iter = c60$PValue
    c2$names = rownames(c60)
    c2 = as.data.frame(c2)
    colnames(c2)= c(paste("iter",i, sep = "."),"names")
    perm_Effect60 = plyr::join(perm_Effect60, c2)
    
    c3 = c()
    c3$iter = c120$PValue
    c3$names = rownames(c120)
    c3 = as.data.frame(c3)
    colnames(c3)= c(paste("iter",i, sep = "."),"names")
    perm_Effect120 = plyr::join(perm_Effect120, c3)
  }
  i = i + 1
}
setwd ("/home/n-z/sbukhar/eFDR/")
write.table(perm_Effect30, file = paste(Brain_region, "_null30_",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(perm_Effect60, file = paste(Brain_region, "_null60_",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(perm_Effect120, file = paste(Brain_region, "_null120_",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")

save.image(file = paste(Brain_region,"perm",Sys.Date(),"RData",sep = "."))
