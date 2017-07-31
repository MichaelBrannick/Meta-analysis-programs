library(metafor)
library(xlsx)

#############################################
# McLeod 2007 is a dataset I worked with on
# earlier videos. Effect sizes represent
# correlations between
# parenting and child depression
#############################################
McLeodDat <- read.xlsx(file.choose(), sheetName = "Data")
McLeodDat
str(McLeodDat)

McLeodRes1 <- rma(ni = N, ri = r, method = "DL",
                  measure = "ZCOR", data = McLeodDat)
McLeodRes1
#############################################
#  Unsorted forest
#############################################
forest(McLeodRes1)

#############################################
#  Sorted by effect size.
#  To refer to a variable within a dataset,
#  name the dataset, use a dollar sign, then
#  name the variable (see below).
##############################################
McLeod.ES <- McLeodDat[order(McLeodDat$r), ]
McLeod.ES
McLeodRes2 <- rma(ni=N, ri=r, method="DL", measure="ZCOR", data=McLeod.ES)
McLeodRes2
forest(McLeodRes2)

##############################################
#  Sort by precision
##############################################
McLeod.V <- McLeodDat[order(McLeodDat$v), ]
McLeod.V
McLeodRes3 <- rma(ni = N, ri = r, method = "DL",
                  measure = "ZCOR", data = McLeod.V)
forest(McLeodRes3)

###############################################
# Sort by moderator, then ES, test for moderator
# forest by moderator and ES with test results
################################################
McLeod.Dx.ES <- McLeodDat[order(McLeodDat$Dx, McLeodDat$r), ]
McLeodRes4 <- rma(ni = N, ri = r, method = "DL",
                  measure = "ZCOR", mods = ~factor(Dx),
                  data = McLeod.Dx.ES)
McLeodRes4
forest(McLeodRes4)

################################################
# Rockstuhl data
################################################
RocksDat <- read.xlsx(file.choose(), sheetName="Data")

################################################
# sorted by date
################################################
RocksDat.Yr <- RocksDat[order(RocksDat$Year), ]
RocksRes1 <- rma(ni = N, ri = r, method = "DL",
                 measure = "ZCOR", data = RocksDat.Yr)
RocksRes1
forest(RocksRes1)
RocksCU <- cumul(RocksRes1, order = order(RocksDat.Yr$Year))
forest(RocksCU)
