library(VirFinder)
args <- commandArgs(TRUE)
predResult <- VF.pred(args[1])
write.table(predResult, args[2], sep='\t')
