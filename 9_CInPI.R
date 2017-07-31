library(metafor)
library(xlsx)
####################################################
# Function Higgins.PI needs these arguments:
# M is the overall meta-analytic mean in the output
# (find under 'estimate' in the printout).
# Tausq is the estimate tau^2 in the output.
# SEM is the estimate 'se' after the mean.
# k is the number of studies in the output.
#####################################################
Higgins.PI <- function(M, Tausq, SEM, k){
   df1 <- k-2
   SEPI <- sqrt(Tausq+SEM^2)
   PI.lb <- M-qt(.975, df=df1)*SEPI
   PI.ub <- M+qt(.975, df=df1)*SEPI
   cbind(PI.lb,PI.ub)
}
#####################################################
McNattDat <- read.xlsx("/Users/michaelbrannick/Desktop/MetaClassData/McNattData.xlsx", sheetName="Sheet1")
McNattDat
McNattRes1 <- rma(yi=d, vi=v, method="DL", data=McNattDat)
McNattRes1
confint(McNattRes1)
predict(McNattRes1)
Higgins.PI(1.0919, .4643, .1861, 17)

######################################################
dat.mcdaniel1994
McDanielRes1 <- rma(ni=ni, ri=ri, method="DL", measure="ZCOR", data=dat.mcdaniel1994)
McDanielRes1
confint(McDanielRes1)
predict(McDanielRes1)
Higgins.PI(.2368, .0268, .0164, 160)
