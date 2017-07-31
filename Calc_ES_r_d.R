# Always save image of the working environment to the current directory.
#save.image('Calc_ES_r_d.RData')
# And, alwasy load the image of the previous working environment before start working on the R script.
#load('Calc_ES_r_d.RData')
# I have attached an Excel file with hypothetical data coded from 4 studies. 3 of the studies have d statistics, the 4th is a correlation coefficient.  All are missing something or other in terms of information. We want to compute a value of d and the variance of d for each of the studies. This is shown in Sheet 1.
# I would do this in Excel, as shown in Sheet 2.
# But suppose you upload the data (sheet 3) to R.  How would you calculate the required quantities in R?
##################################################
# In order to use functions you have to 
# write them, select
# them and run them.
##################################################
calc.Vd.from.d <- function(n1, n2, d){
    f <- (n1 + n2)/(n1 * n2)
    s <- (d^2)/(2 * (n1 + n2))
    f + s
}
calc.d.from.r <- function(r){
    (2*r)/sqrt(1-r^2)
}
calc.Vr.from.r <- function(r, N){
    num <- (1-r^2)^2
    den <- N-1
    num/den
}
calc.Vd.from.r <- function(r, N){
	num1 <-(1-r^2)^2
	den1 <-N-1
    num2 <- 4*(num1/den1)
    den2 <- (1-r^2)^3
    num2/den2
}
calc.Vd.from.95ci <- function(ciup, cilo){
    num <- ciup-cilo
    den <- 2*1.96
    (num/den)^2
}
calc.g <- function(n1, n2, d){
    num <- 3
    den <- (4*(n1+n2-2))-1
    (1-(num/den))*d
}
calc.Vg <- function(n1, n2, Vd){
    num <- 3
    den <- (4*(n1+n2-2))-1
    (1-(num/den))*Vd
}
##################################################
library(metafor)
library(xlsx)
##################################################
es3 <- read.xlsx('Calc_ES_r_d.xlsx',
                 sheetIndex = 3)
options(digits=3)
es3
# Calculation of Vd from Ns and ds and also from
# confidence intervals. But this is not particularly
# easier than doing that with excel.
es3$Vd <- with(es3, calc.Vd.from.d(Nexp, Nc, d))
es3$Yd <- es3$d
es3
es3$Vd[3:4] <- with(es3[3:4, ],calc.Vd.from.95ci(CILo,CIUp))
es3
es3$Vd[4] <- with(es3[4, ],calc.Vd.from.r(r, Nt))
es3$Yd[4] <- with(es3[4, ],calc.d.from.r(r) )
es3
#####################that's it#####################