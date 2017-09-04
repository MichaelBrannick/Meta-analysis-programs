#install.packages("forestplot")
library(forestplot)
library(xlsx)
initial.data <- file.choose()
Creds <- read.xlsx(initial.data, sheetName="Sheet1")
#Creds <- read.xlsx("/Users/michaelbrannick/Desktop/ReFor.xlsx", sheetName="Sheet1")
#### Check the data
Creds
tabs<-cbind(
  c("Author", "Baltes et al. (1999)", "Ben-Shakhar & Elaad (2003)", "Dwight & Donovan (2003)", 
    "Hosoda et al. (2003)", "Kierein & Gold (2000)", "Nguyen & Ryan (2008)", "Richman et al. (1999)", 
    "Stewart & Roth (2001)"),
  c("k", "36", "169", "12", "62", 
    "13", "116", "95", "12"))
#  c("r(ES,N)", "-.21", ".01", ".33", "-.03", "-.16", ".04", "-.11", ".53"))
tabs
HigsrM <-1:9
HigsrL <- HigsrM
HigsrU <- HigsrM
HigsrM[1]<- NA
HigsrL[1]<- NA
HigsrU[1]<- NA
HigsrM[2:9]<- Creds$HigM
HigsrL[2:9]<- Creds$HigL
HigsrU[2:9]<- Creds$HigU
#
HigszM <-1:9
HigszL <- HigszM
HigszU <- HigszM
HigszM[1]<- NA
HigszL[1]<- NA
HigszU[1]<- NA
HigszM[2:9]<- Creds$HigMz
HigszL[2:9]<- Creds$HigLz
HigszU[2:9]<- Creds$HigUz
#
SnHsM <-1:9
SnHsL <- SnHsM
SnHsU <- SnHsM
SnHsM[1]<- NA
SnHsL[1]<- NA
SnHsU[1]<- NA
SnHsM[2:9]<- Creds$SnHM
SnHsL[2:9]<- Creds$SnHL
SnHsU[2:9]<- Creds$SnHU

#pdf("CredInts.pdf", family="Times", height=10, width=7)
#par(mar=c(4,4,4,4))

forestplot(tabs, 
           txt_gp = fpTxtGp(ticks = gpar(cex=.8),      # numerical tick labels
                            xlab  = gpar(cex = .8),    # label you input at bottom
                            label = gpar(cex = .8)),   # stuff in 'tabs' file & legend
           legend = c("Schmidt-Hunter", "Higgins-t", "Higgins-z"),          # the legend you input
           fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI, fpDrawCircleCI), # marker (effect size) shapes      
           is.summary = c(TRUE,rep(FALSE,8)),          # to leave space
           mean = cbind(SnHsM, HigsrM, HigszM),
           lower = cbind(SnHsL, HigsrL, HigszL),
           upper = cbind(SnHsU, HigsrU, HigszU),
           clip =c(-1.1, 2.5),                         # limits of forest (bottom, top)
           lty.ci=c(1,2,3),                            # line type for wings - solid, dashed, etc.
           lwd.ci = 2,                                 # line width (thickness) for wings
           col=fpColors(box=c("blue","green", "darkred")), # colors of marks (boxes)
           vertices=TRUE,                              # markers at end of wings
           xlab="Standardized Mean Difference",        # input your label
           new_page = TRUE,
           boxsize=.2,                                 # size of markers
           #grid=TRUE
           grid=structure(c(-1,-.5, 0,.5, 1,1.5, 2, 2.5), # vertical line placements
                          gp=gpar(lty=2, lwd=2))          # vertical line type & width
           )
#dev.off()
#help(forestplot)
#
# abovle line [col=...] is for different colors for each method in the graph
# taken out for the pub version to be black only
