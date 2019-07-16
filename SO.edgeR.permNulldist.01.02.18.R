#This scripts computes null distributions which later will be used to compute eFDRs.

rm(list = ls())
stopifnot(require("edgeR"))

#### Parameters
Brain_region = "T"
dir = ".tagwise"
####

setwd ("/Volumes/GoogleDrive/My Drive/social_oppurtunity/challenge.opp.TRN/counts/")

SC_targets = readTargets("Targets.txt")
SO_targets = readTargets("sb_mrsb_opportunity_targets.txt")
all_targets = c()
all_targets = as.data.frame(union(SC_targets$files,SO_targets$files))
colnames(all_targets) = "files"


targets = readTargets("sb_mrsb_opportunity_targets.txt")
targets = targets[which(targets$region == Brain_region),]

getDGE = function(targets)
{
  Raw_DGE = readDGE(targets)
  MetaTags <- grep("^__", rownames(Raw_DGE))
  
  Raw_DGE = Raw_DGE[-MetaTags,]
  
  keep <- rowSums(cpm(Raw_DGE)>0.5) >= 4
  Filtered_DGE = Raw_DGE[keep,]
  Filtered_DGE = calcNormFactors(Filtered_DGE)
}

DGE = getDGE(all_targets)
targets$group = paste(targets$nest_condition, targets$timepoint, sep = ".")
group = as.factor(targets$group)
group = relevel(group, ref = "NOCLUTCH.9am")

design_Group = model.matrix(~group)


get_fit = function(y, design)
{
  y = estimateGLMCommonDisp(y,design, verbose=F)
  y = estimateGLMTrendedDisp(y,design, verbose=F)
  y = estimateGLMTagwiseDisp(y,design)
  fit = glmFit(y, design)
  return(fit)
}

fit.Group = get_fit(DGE, design_Group)

perm_Nest = c()
perm_Eggs = c()
perm_HEarly = c()
perm_HMid = c()
perm_HLate = c()

library(plyr)
totalperm = 500
for(i in 1:totalperm)
{
  # using a time based seed.
  t <- as.numeric(Sys.time())
  seed <- 1e8 * (t - floor(t))
  set.seed(seed)
  targets_perm = targets
  DGE_perm = getDGE(targets_perm)
  #shuffling sample names
  DGE_perm$samples = DGE_perm$samples[sample(rownames(DGE_perm$samples)),]
  # computing glmFIT using shuffled samples and original design matrix which was used to obtain original p values.
  perm.fit = get_fit(DGE_perm, design_Group)
  #since reference is "NOCLUTCH.9am" so thats to estimate Nest, Eggs and Early simply using their coeficients. To compare
  # Mid against control at 1p and Late against control at 5pm. specific contrasts are designated.
  Nest = as.data.frame(topTags(glmLRT(perm.fit, coef=6), n=Inf, sort.by="PValue"))
  Eggs = as.data.frame(topTags(glmLRT(perm.fit, coef=2), n=Inf, sort.by="PValue"))
  HEarly = as.data.frame(topTags(glmLRT(perm.fit, coef=5), n=Inf, sort.by="PValue"))
  HMid = as.data.frame(topTags(glmLRT(perm.fit, contrast=c(0,0,1,0,0,0,-1,0)), n=Inf, sort.by="PValue"))
  HLate = as.data.frame(topTags(glmLRT(perm.fit, contrast=c(0,0,0,1,0,0,0,-1)), n=Inf, sort.by="PValue"))
 # head(HLate)
  # collecting null p values in respected dataframes.
  if(i == 1)
  {
    perm_Nest$iter1 = Nest$PValue
    perm_Nest$names = rownames(Nest)
    perm_Nest = as.data.frame(perm_Nest)
    
    perm_Eggs$iter1 = Eggs$PValue
    perm_Eggs$names = rownames(Eggs)
    perm_Eggs = as.data.frame(perm_Eggs)
    
    perm_HEarly$iter1 = HEarly$PValue
    perm_HEarly$names = rownames(HEarly)
    perm_HEarly = as.data.frame(perm_HEarly)
    
    perm_HMid$iter1 = HMid$PValue
    perm_HMid$names = rownames(HMid)
    perm_HMid = as.data.frame(perm_HMid)
    
    perm_HLate$iter1 = HLate$PValue
    perm_HLate$names = rownames(HLate)
    perm_HLate = as.data.frame(perm_HLate)
    
  }
  else
  {
    c1 = c()
    c1$iter = Nest$PValue
    c1 = as.data.frame(c1)
    c1$names = rownames(Nest)
    colnames(c1)= c(paste("iter",i, sep = "."),"names")
    perm_Nest = plyr::join(perm_Nest, c1)
    
    c2 = c()
    c2$iter = Eggs$PValue
    c2 = as.data.frame(c2)
    c2$names = rownames(Eggs)
    colnames(c2)= c(paste("iter",i, sep = "."),"names")
    perm_Eggs = plyr::join(perm_Eggs, c2)
    
    c3 = c()
    c3$iter = HEarly$PValue
    c3 = as.data.frame(c3)
    c3$names = rownames(HEarly)
    colnames(c3)= c(paste("iter",i, sep = "."),"names")
    perm_HEarly = plyr::join(perm_HEarly, c3)
    
    c4 = c()
    c4$iter = HMid$PValue
    c4 = as.data.frame(c4)
    c4$names = rownames(HMid)
    colnames(c4)= c(paste("iter",i, sep = "."),"names")
    perm_HMid = plyr::join(perm_HMid, c4)
    
    c5 = c()
    c5$iter = HLate$PValue
    c5 = as.data.frame(c5)
    c5$names = rownames(HLate)
    colnames(c5)= c(paste("iter",i, sep = "."),"names")
    perm_HLate = plyr::join(perm_HLate, c5)
  }
}

# writing null distributions to a specified directory below.
setwd ("/Users/syedabbasbukhari/Google Drive/social_oppurtunity/challenge.opp.TRN/")
### Creat Results Directory.
if(dir.exists(file.path(paste(getwd(),"/",Brain_region,"_Null_perm",totalperm,".",Sys.Date(), sep = ""))) == FALSE)
{
  dir.create(file.path(paste(getwd(),"/",Brain_region,"_Null_perm",totalperm,".",Sys.Date(), sep = "")) ,showWarnings = TRUE)
  setwd(file.path(paste(getwd(),"/",Brain_region,"_Null_perm",totalperm,".",Sys.Date(), sep = "")))
} else {
  setwd(file.path(paste(getwd(),"/",Brain_region,"_Null_perm",totalperm,".",Sys.Date(), sep = "")))
}

write.table(perm_Nest, file = paste(Brain_region, "_nullNest",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(perm_Eggs, file = paste(Brain_region, "_nullEggs",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(perm_HEarly, file = paste(Brain_region, "_nullHEarly",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(perm_HMid, file = paste(Brain_region, "_nullHMid",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(perm_HLate, file = paste(Brain_region, "_nullHLate",Sys.Date(), ".txt", sep = ""), row.names = F, quote = F, sep = "\t")

# saving session objects.
save.image(file = paste(Brain_region,"perm_SO",Sys.Date(),sep = "."))
