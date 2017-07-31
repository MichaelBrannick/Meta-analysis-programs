library(metafor)
library(xlsx)
SleepDat <- read.xlsx(file.choose(), sheetName="Mods")
SleepDat

############################  1 input
SleepDat.Grp.ES <- SleepDat[order(SleepDat$group, SleepDat$yi), ]
SleepDat.Grp.ES
SleepDatRes1 <- rma(yi = yi, vi = vi,
                    method = "DL", data = SleepDat.Grp.ES)
SleepDatRes1
forest(SleepDatRes1)
# to see the contents of the object created by the forest routine
#forestPars <- forest(SleepDatRes1)
#forestPars
#############################  2 single graph with all data

###################################
# subset data - last two groups
###################################
SleepDat.34 <- SleepDat.Grp.ES[26:48, ]
SleepDatRes.34 <- rma(yi = yi, vi = vi, method = "DL",
                      mods = ~factor(group), data = SleepDat.34)

###########################################################
#  I should have subset in the forest rather than the data,
#  but pretend we only have 2 groups.
############################################################
forest(SleepDatRes.34)

############################################################ plot with just 2 groups and 23 ES
forest(SleepDatRes.34, slab =i paste(SleepDat.34$Author, SleepDat.34$Year, sep =", "))

############################################################ add study labels
forest(SleepDatRes.34, slab = paste(SleepDat.34$Author, SleepDat.34$Year, sep = ", "), atransf = exp)

############################################################# translate to odds from OR
forest(SleepDatRes.34, slab = paste(SleepDat.34$Author, SleepDat.34$Year, sep = ", "), atransf = exp, xlab = "Odds Ratio")

############################################################## label ES scale Odds not OR
forest(SleepDatRes.34, slab = paste(SleepDat.34$Author, SleepDat.34$Year, sep = ", "), atransf = exp, xlab = "Odds Ratio", at = c(-1.5, 0, 1.5, 3))

###########################################################
forest(SleepDatRes.34, atransf = exp, xlab = "Odds Ratio", slab = paste(SleepDat.34$Author, SleepDat.34$Year2, sep = ", "), at = c(-1.5, 0, 1.5, 3))
op <- par(cex = 1, font = 2)
text(5.5, 25, "OR [95% CI]")
text(-3, 25, "Author, Date")
par(op)
# you have to play with the numbers for placement; trial and error
# start with zero, zero to be in the graph and work from there

################################################################ add column labels
forest(SleepDatRes.34, atransf = exp, xlab = "Odds Ratio", slab = paste(SleepDat.34$Author, SleepDat.34$Year2, sep = ", "), at = c(-1.5, 0, 1.5, 3))
op <- par(cex = .75, font = 2)
text(3, 20, "Insomnia")
text(3, 19, "1.54 [1.26, 1.93]" )
text(3, 6, "Obstructive Apnea")
text(3, 5, "1.92 [1.57, 2.34]")
par(op)
op <- par(cex = 1, font = 2)
text(5.5, 25, "OR [95% CI]")
text(-3, 25, "Author, Date")
par(op)
