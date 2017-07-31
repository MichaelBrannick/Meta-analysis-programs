####################################
# Code for video 6, method and
# measure for r and d.
####################################
# setup - install R & load libraries if not already done
####################
library(metafor)
library(xlsx)
#########################
# setup read datasets if not already done
#########################
# RocksDat <- read.xlsx("/Users/michaelbrannick/Desktop/RocksLMX_AC.xlsx",  sheetName="data")
# AppetiteDat <- read.xlsx("/Users/michaelbrannick/Desktop/AppetiteStudy.xlsx", sheetName="Sheet1")
McNattDat <- read.xlsx("/Users/michaelbrannick/Desktop/McNattData.xlsx", sheetName="Data")
######################### 
# new slide
#########################
McNattDat
McNattRes1 <- rma(yi=d, vi=v, data=McNattDat)
McNattRes1 
######################### 
# new slide
#########################
McNattRes2 <- rma(yi=d, vi=v, method="DL", data=McNattDat)
McNattRes2
McNattRes3 <- rma(yi=d, vi=v, method="FE", data=McNattDat)
McNattRes3
