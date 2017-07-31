library(metafor)
library(xlsx)
McLeodDat <- read.xlsx("/Users/michaelbrannick/Desktop/McLeod2007.xlsx", sheetName="Data")
McLeodDat
McLeodRes1 <- rma(ni=N, ri=r, method="DL", measure="ZCOR", data= McLeodDat)
McLeodRes1
McLeodRes2 <- rma(ni=N, ri=r, method="DL", measure="ZCOR", mods = ~factor(Dx), data= McLeodDat)
McLeodRes2
McLeodRes3 <- rma(ni=N, ri=r, method="DL", measure="ZCOR", mods = ~factor(Dx), knha=TRUE, data= McLeodDat)
McLeodRes3
McLeodRes4 <- rma(ni=N, ri=r, method="DL", measure="ZCOR", mods = ~factor(Dx)-1, data= McLeodDat)
McLeodRes4
McLeodRes5 <- rma(ni=N, ri=r, method="DL", measure="ZCOR", mods = ~factor(Dx)+Age, data= McLeodDat)
McLeodRes5
