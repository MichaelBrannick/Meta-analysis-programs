### Leave 1 out with McNatt data
library(metafor)
library(xlsx)
####### read in McNatt data
McNattDat <- read.xlsx("McNattData.xlsx", sheetName="Sheet1")
McNattDat
####### run the meta-analysis
McNattRes1 <- rma(yi=d, vi=v, method="DL", data=McNattDat)
McNattRes1
####### sort by precision
McNattDat.sort <- McNattDat[order(McNattDat$v), ]
McNattDat.sort
####### rerun meta-analysis on sorted data
McNattRes2  <- rma(yi=d, vi=v, method="DL", data=McNattDat.sort)
####### plot forest by precision
forest(McNattRes2)
#########################################################
####### test for funnel asymmetry
regtest(McNattRes1)
###### look for outliers with residuals
rstudent(McNattRes1)
###### study in row 10 is suspicious, rerun without it
McNattRes3 <- rma(yi=d, vi=v, method="DL", data=McNattDat[-c(10), ])
McNattRes3
######  run the leave one out analysis
res2 <- leave1out(McNattRes1)
res2
######  forest plot of overall results leaving out each study
forest(res2$estimate, res2$se^2)
###### also shows impact of study #10


