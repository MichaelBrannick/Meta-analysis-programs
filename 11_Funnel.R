library(metafor)
library(xlsx)
########################################
dat.mcdaniel1994
McDanielRes1 <- rma(ri=ri, ni=ni, method="DL", measure="ZCOR", data=dat.mcdaniel1994)
McDanielRes1
funnel(McDanielRes1)
predict(McDanielRes1)
########################################
McDanielRes2 <- trimfill(McDanielRes1)
McDanielRes2
funnel(McDanielRes2)
predict(McDanielRes2)
#########################################
McDanielRes3 <- rma(ri=ri, ni=ni, method="DL", measure="COR", data=dat.mcdaniel1994)
McDanielRes3
funnel(McDanielRes3)
predict(McDanielRes3)
#########################################
McDanielRes4 <- trimfill(McDanielRes3)
McDanielRes4
funnel(McDanielRes4)
predict(McDanielRes4)