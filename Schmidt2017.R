######################################################################################
# Schmidt and Hunter Meta-analysis in r, the correlation coefficient.
# The user can choose artifact corrections from none (bare bones) to range restriction
# either direct or indirect and reliability for the independent and dependent variables.
# The program works for correlations that are 
# individually corrected, and also allows for assumed distributions.
# Additional features include bootstrap confidence intervals for credibility value and
# the user may specify Morris estimation rather than Schmidt-Hunter.
# ####################################################################################
# This version of the program expects to find (you must input)
#  the observed correlation, ri
#  the observed sample size ni
# You may also input the following:
#  the ratio of restricted to unrestricted SD in the IV, RRuX;
#  the range restricted reliability of the independent variable, rXXi;
#  the range restricted reliability of the dependent variable, rYYi.
#  For example (these data from the built-in dataset from the Le and Schmidt Windows program):
#
# Study	ri	Nr	RRuX	RXXi	Ryyi
# 1	  0.35	68	0.58	0.55	0.8	
# 2	  0.07	68	0.58	0.55	0.6	
# 3	  0.11	68	0.58	0.55	0.8	
# 4	  0.31	68	0.58	0.55	0.6	
# 5	  0.18	68	0.678	0.67	0.8	
# 6	  0.36	68	0.678	0.67	0.6	
# 7	  0.40	68	0.678	0.67	0.8	
# 8	  0.13	68	0.678	0.67	0.6	
# 9	  0.49	68	0.869	0.80	0.8	
# 10	0.23	68	0.869	0.80	0.6	
# 11	0.29	68	0.869	0.80	0.8	
# 12	0.44	68	0.869	0.80	0.6	
#  
# The program will compute the rest of the quantities you need.  
# 
# Studies missing data for ri and ni will be excluded from analyses by default
#
# This program is for INDIVIDUAL artifact corrections, not distribution artifacts.
#
# You may also obtain boostrap estimates for the confidence interval of the lower
#  bound of the credibility interval, and set the number of iterations for the 
#  boostrap (the default is 10000 iterations).  
# If you use the Morris model and boostrap confidence intervals, expect the program
#  to take a minute or two to run becuase both Morris and boostrap are iterative.
######################################################################################
# preparation - also install any needed supporting programs (xlsx, metafor)
rm(list=ls()) # Cleans (empties) global environment
dev.off()     # Cleans plots/graphics
cat("\014")   # Cleans out the console
######################################################################################
#library(xlsx)                                             # required to read an Excel file
library(metafor)                                          # required for Morris 
library(boot)                                             # required for bootstrapping
data.file <- file.choose()                                # find the file you want to read
setwd(dirname(data.file))                                 # change the working direct to match file
initial.data <- read.csv(data.file)  # read the file - choose the sheet
main.data <- initial.data                                 # rename
main.data                                                 # print the data
str(main.data)
#####################################################################################
#  SnH computations follow Schmidt & Hunter (2015)
#  Methods of Meta-Analysis (3rd ed.)
#####################################################################################
# Studies that are missing values for ri and ni will be excluded from analyses
#
# If you do not have data for an artifact, omit it in the function call.
# The program uses the value '1' by default for artifacts, which has the effect
# of ignorning the factor, that is, no inclusion of artifacts results in bare-bones.
#
# Input MetaModel = 'SnH' in the function call for Schmidt-Hunter;
# use 'Morris' for Morris.
#####################################################################################
Schmidt <-function(ri, ni, rxxi=rep.int(1,length(ri)), ryyi=rep.int(1,length(ri)),
                   ui=rep.int(1,length(ri)), moderator=rep.int(1,length(ri)),Bias_Correct=TRUE,
                   k_Correct=TRUE, IndirectRR=TRUE, MetaModel='SnH', PredInt=FALSE, Boots=FALSE,
                   BootIter=5000, Missing=FALSE)
{                                                       # begin function Schmidt - set up
  moderator <- factor(moderator)
  input <- data.frame(ri,ni,rxxi,ryyi,ui,moderator)     # generate dataset based on input
  initial.count <- nrow(input)                          # count number of initial studies
  # 
  #Function for imputing missing artifact values
  #  
  impute <- function(data, type) {
    for (i in which(sapply(data, is.numeric))) 
    {
      data[is.na(data[, i]), i] <- type(data[, i],  na.rm = TRUE)
    }
    return(data)}
  
  #Excludes studies with either missing ni or ri only (if user is imputing artifacts) or excludes if missing any value   
  if (Missing) {
    complete.data <- input[complete.cases(input[,c(1,2,6)]),] # exclude studies with missing ri, ni, or moderator
    
    complete.count <- nrow(complete.data)                   # count remaining studies
    missing.cases <- initial.count-complete.count           # count number of excluded studies 
    artifact.data <- data.frame(complete.data$rxxi,complete.data$ryyi,complete.data$ui)               # create dataset for imputing artifact values only
    impute.data <- impute(artifact.data,mean)               # impute artifact values based on means
    
    complete.data$rxxi <- impute.data[,1]                                # update artifact input values 
    complete.data$ryyi <- impute.data[,2]
    complete.data$ui <- impute.data[,3]
    
  } else {
    complete.data <- input[complete.cases(input[,c(1,2,3,4,5,6)]),] # exclude studies with any missing values
    complete.count <- nrow(complete.data)                           # count remaining studies
    missing.cases <- initial.count-complete.count                   # count remaining studies
  }  
  
  
  #
  # Function for prediction intervals         #
  #
  Preds <- function(rbar.a, V.rho, V.m, k)  # define the function and arguments
  {                                           # begin computations of the function
    SE.PI <- sqrt(V.rho+V.m)                  # combined variance of mean and rho
    t80.value <- qt(.900, df=(k-2))           # t value for 80 percent prediction interval
    PI80.U <- rbar.a + t80.value*SE.PI        # 80 pct prediction interval - upper
    PI80.L <- rbar.a - t80.value*SE.PI        # 80 pct prediction interval - lower
    PIs <- c(PI80.L, PI80.U)                  # group prediction interval into a vector
  }                                           # end prediction interval computations
  # how many levels of the moderator?
  #   moderator <- factor(moderator)            # make sure the moderator is categorical
  mod.levels <- length(levels(moderator))   # count the number of levels
  mod.spots <- 26*mod.levels                # 20 is the number of output rows
  mod.count <- as.numeric(complete.data$moderator)
  complete.data$mod.count <- mod.count
  Desc.data <- matrix(1:mod.spots,ncol=mod.levels) #placeholder for results output by moderator
  ########################################################a###############################
  #  Main Loop
  for(zz in 1:mod.levels){                     # run the whole thing as many times as levels
    ri <- subset(complete.data$ri, mod.count==zz)
    ni <- subset(complete.data$ni, mod.count==zz)
    rxxi <- subset(complete.data$rxxi, mod.count==zz)
    ryyi <- subset(complete.data$ryyi, mod.count==zz)
    ui <- subset(complete.data$ui, mod.count==zz)
    mod.data <- subset(complete.data, mod.count==zz)
    ri <- ifelse(ri == 0, 0.0000001, ri)        # avoid dividing by zero later on
    k <- length(ri)                             # number of effect sizes
    nR <- ri*ni                                 # weight r by N
    Nsum <- sum(ni)                             # find sum of N
    rbar <- sum(nR)/Nsum                        # find sample-weighted mean r
    bias.factor <- (2*ni-2)/(2*ni-1)             # to find the unbiased estimates S&H (2015) p. 67
    if (Bias_Correct)
      {
    ri <- ri/bias.factor                        # unbiased estimates
      }
    PIs <- c(NA, NA)                            # place holder for prediction intervals
    CR.CI.LB <- NA                              # place holders for confidence interval for lower bound
    CR.CI.UB <- NA
    REVC.CI.LB <- NA
    REVC.CI.UB <- NA
    Mbt.CI.LB <- NA
    Mbt.CI.UB <- NA
    #
    #######################################################################################
    ### Start computations for the Schmidt-Hunter model
    if(MetaModel=='SnH')  # start computations for the Schmidt-Hunter model
    {
      ##########################################
      # Direct Range Restriction
      ##########################################
      if(!IndirectRR){ # start computing for Direct Range Restriction 
        nR <- ri*ni                             # weight r by N
        Nsum <- sum(ni)                            # find sum of N
        rbar <- sum(nR)/Nsum                       # find sample-weighted mean r
        V.obs <- sum(ni*(ri-rbar)^2)/(sum(ni))     # find weighted observed variance
        if(k_Correct){V.obs <- V.obs*(k/(k-1))}
  #      rxxa <- 1-ui^2*(1-rxxi)                   # find reliability in applicant sample
        rxxa <- ((rxxi^.5/ui)/sqrt(1+rxxi*(ui^-2-1)))^2 # find reliability applicant sample (direct RR)
        Ux <- 1/ui                                # inverse of range restriction ratio
        rC <- (Ux*ri)/(sqrt(rxxa*(ryyi+Ux^2*ri^2-ri^2)))    # find corrected correlations
        A.compound <- ri/rC                      # find A, the attenuation factor
        wi <- A.compound^2*ni                    # find the weights for the meta
        rbarC.1 <- sum(wi*rC)/sum(wi)            # find the mean corrected correlations
        V.rC <- sum(wi*(rC-rbarC.1)^2)/sum(wi)   # find the variance of the corrected correlations
        if(k_Correct){V.rC <- V.rC*(k/(k-1))}   # apply k correction to variance estimate
        V.eo <- (1-rbar^2)^2/(ni-1)              # simple observed error variance for each study
        V.ec1 <- V.eo/A.compound^2               # simple error variance of corrected correlations
        err.adj <- 1/((Ux^2-1)*ri^2+1)           # range restriction error variance correction
        V.ve <- sum(wi*V.ec1*err.adj^2)/sum(wi)  # find error variance of corrected correlations
        V.rho.1 <- V.rC-V.ve                     # find the random-effects variance component
        if (V.rho.1 < 0) {V.rho.1 <- 0}          # if REVC is less than zero, set to zero
        SD.rho.1 <- sqrt(V.rho.1)                # find SD rho
        SEM <- sqrt(V.rC/k)                      # find standard error of the mean of corrected correlations
        V.m <- SEM^2                             # Sampling variance of the mean
        CI95.U.1 <- rbarC.1 + 1.96*SEM           # find the confidence interval for the mean of corrected corrs
        CI95.L.1 <- rbarC.1 - 1.96*SEM
        CR80.U.1 <- rbarC.1 + 1.28*SD.rho.1      # find the bounds of the credibility interval
        CR80.L.1 <- rbarC.1 - 1.28*SD.rho.1

        # Bootstrap CI for lower bound
        SH.DirectRR.f <- function(d, i){         # Boostrap function for credibility interval lower bound
          d2 <- d[i,]
          boot.ri <- d2$ri
          boot.ri <- ifelse(boot.ri == 0, 0.0000001, boot.ri)   # avoid dividing by zero later on
          boot.ni <- d2$ni
          boot.rxxi <- d2$rxxi
          boot.ryyi <- d2$ryyi
          boot.u <- d2$ui
          nR <- boot.ri*boot.ni
          Nsum <- sum(boot.ni)
          rbar <- sum(nR)/Nsum
          Ux <- 1/boot.u
          rxxa <- ((rxxi^.5/ui)/sqrt(1+rxxi*(ui^-2-1)))^2 #
          rC <- (Ux*ri)/(sqrt(rxxa*(ryyi+Ux^2*ri^2-ri^2)))    # find corrected correlations
          A.compound <- ri/rC                      # find A, the attenuation factor
          wi <- A.compound^2*boot.ni               # find the weights for the meta
          rbarC.1 <- sum(wi*rC)/sum(wi)            # find the mean corrected correlations
          V.rC <- sum(wi*(rC-rbarC.1)^2)/sum(wi)   # find the variance of the corrected correlations
          if(k_Correct){V.rC <- V.rC*(k/(k-1))}   # apply k correction to variance estimate
          V.eo <- (1-rbar^2)^2/(boot.ni-1)         # simple observed error variance for each study
          V.ec1 <- V.eo/A.compound^2               # simple error variance of corrected correlations
          err.adj <- 1/((Ux^2-1)*boot.ri^2+1)      # range restriction error variance correction
          V.ve <- sum(wi*V.ec1*err.adj^2)/sum(wi)  # find error variance of corrected correlations
          V.rho.1 <- V.rC-V.ve                     # find the random-effectrn(s variance component
          if (V.rho.1 < 0) {V.rho.1 <- 0}          # if REVC is less than zero, set to zero
          SD.rho.1 <- sqrt(V.rho.1)                # find SD rho
          lb <- rbarC.1 - 1.28*SD.rho.1            # find lower bound
          DirectBts <- c(rbarC.1, V.rho.1, lb)     # collect estimates to be examined for confidenc intervals
          return(DirectBts)
        } # end boostrap function for lower bound credibility
        if (Boots) # start bootstrap for direct RR for Schmidt & Hunter
        {
          SH.DirectRR.bootstrap <- boot(mod.data, SH.DirectRR.f, R = BootIter)
          SH.DirectRR.bootstrap1 <- boot.ci(SH.DirectRR.bootstrap, type="bca", index=1)
          Mbt.CI.LB <- SH.DirectRR.bootstrap1$bca[4]
          Mbt.CI.UB <- SH.DirectRR.bootstrap1$bca[5]
          SH.DirectRR.bootstrap2 <- boot.ci(SH.DirectRR.bootstrap, type="bca", index=2)
          REVC.CI.LB <- SH.DirectRR.bootstrap2$bca[4]
          REVC.CI.UB <- SH.DirectRR.bootstrap2$bca[5]
          SH.DirectRR.bootstrap3 <- boot.ci(SH.DirectRR.bootstrap, type="bca", index=3)
          CR.CI.LB <- SH.DirectRR.bootstrap3$bca[4]
          CR.CI.UB <- SH.DirectRR.bootstrap3$bca[5]
        } # end bootstrap for direct RR
       
      # Common output                          # so the names are the same regardless of model
      rbar.a <- rbarC.1
      CI95.L <- CI95.L.1
      CI95.U <- CI95.U.1
      V.rho <- V.rho.1
      V.e <-V.ve
      PctVE <- V.e/(V.e+V.rho)
      SD.rho <- SD.rho.1
      CR80.L <- CR80.L.1
      CR80.U <- CR80.U.1
      #
      PI80.U <- NA                             # If prediction intervals are not wanted
      PI80.L <- NA
      if (PredInt) {                           # If prediction intervals are requested
        PIs <- Preds(rbar.a,V.rho,V.m,k)
      }# end PredInt
      }# end Direct RR
      
      if (IndirectRR) {
        #######################################################################
        # Indirect Range Restriction
        #######################################################################
        nR <- ri*ni                                 # weight r by N
        Nsum <- sum(ni)                             # find sum of N
        rbar <- sum(nR)/Nsum                        # find sample-weighted mean r
        V.obs <- sum(ni*(ri-rbar)^2)/(sum(ni))      # find weighted observed variance
        if(k_Correct){V.obs <- V.obs*(k/(k-1))}    # k correction
        RXXa <-1-ui^2*(1-rxxi)                      # find the reliability of the unrestricted sample
        UT <- 1/sqrt((ui^2-(1-RXXa))/RXXa)              # find UT, disattenuation for indirect RR
        rdiss.rel <- ri/sqrt(rxxi*ryyi)                 # find correlation disattenuated for reliability in X and Y
        rC.2 <- 
          (rdiss.rel*UT)/(sqrt((UT^2-1)*rdiss.rel^2+1)) # find corrected correlations
        A.compound.2 <- ri/rC.2                         # find compound correction factor
        wi.2 <- ni*A.compound.2^2                       # find the weights
        rbarC.2 <- sum(wi.2*rC.2)/sum(wi.2)             # find the mean corrected correlation
        V.rC.2 <- sum(wi.2*(rC.2-rbarC.2)^2)/sum(wi.2)  # find the corrected total variance
        if(k_Correct){V.rC.2 <- V.rC.2*(k/(k-1))}   # apply k correction to variance estimate
        V.eo.2 <- (1-rbar^2)^2/(ni-1)                   # simple observed error variance for each study
        V.ec2 <- V.eo.2/A.compound.2^2                  # find the corrected error variance
        err.adj.2 <- 1/((UT^2-1)*ri^2+1)                # find the error adjustment for indirect range restriction
        V.ve.2 <- sum(wi.2*V.ec2*err.adj.2^2)/sum(wi.2) # find the error variance for the meta
        V.rho.2 <- V.rC.2-V.ve.2                        # find the random-effects variance component
        if (V.rho.2 < 0) {V.rho.2 <- 0}                 # if REVC is less than zero, set to zero
        SD.rho.2 <- sqrt(V.rho.2)                       # find SD rho
        SEM2 <- sqrt(V.rC.2/k)                          # find standard error of the mean of corrected correlations
        V.m <- SEM2^2
        CI95.U.2 <- rbarC.2 + 1.96*SEM2                 # find the confidence interval for the mean of corrected corrs
        CI95.L.2 <- rbarC.2 - 1.96*SEM2
        CR80.U.2 <- rbarC.2 + 1.28*SD.rho.2             # find the lower bounds of the credibility interval
        CR80.L.2 <- rbarC.2 - 1.28*SD.rho.2
        # Common output
        rbar.a <- rbarC.2
        CI95.L <- CI95.L.2
        CI95.U <- CI95.U.2
        V.rho <- V.rho.2
        V.e <-V.ve.2
        PctVE <- V.e/(V.e+V.rho)
        SD.rho <- SD.rho.2
        CR80.L <- CR80.L.2
        CR80.U <- CR80.U.2
        #
        #
        PI80.U <- NA
        PI80.L <- NA
        PIs <- c(PI80.L, PI80.U)
        if (PredInt) {
          PIs <- Preds(rbar.a,V.rho,V.m,k)
        }# end PredInt
        # Bootstrap CI for lower bound
        SH.IndirectRR.f <- function(d, i){
          d2 <- d[i,]
          boot.ri <- d2$ri
          boot.ri <- ifelse(boot.ri == 0, 0.0000001, boot.ri)  
          boot.ni <- d2$ni
          boot.rxxi <- d2$rxxi
          boot.ryyi <- d2$ryyi
          boot.u <- d2$ui
          nR <- boot.ri*boot.ni
          Nsum <- sum(boot.ni)
          rbar <- sum(nR)/Nsum
          RXXa <-1-boot.u^2*(1-boot.rxxi)                  # find the reliability of the unrestricted sample
          UT <- 1/sqrt((boot.u^2-(1-RXXa))/RXXa)           # find UT, disattenuation for indirect RR
          rdiss.rel <- boot.ri/sqrt(boot.rxxi*boot.ryyi)   # find correlation disattenuated for reliability in X and Y
          rC.2 <- 
            (rdiss.rel*UT)/(sqrt((UT^2-1)*rdiss.rel^2+1))  # find corrected correlations
          A.compound.2 <- boot.ri/rC.2                     # find compound correction factor
          wi.2 <- boot.ni*A.compound.2^2                   # find the weights
          rbarC.2 <- sum(wi.2*rC.2)/sum(wi.2)              # find the mean corrected correlation
          V.rC.2 <- sum(wi.2*(rC.2-rbarC.2)^2)/sum(wi.2)   # find the corrected total variance
          if(k_Correct){V.rC.2 <- V.rC.2*(k/(k-1))}       # apply k correction to variance estimate
          V.eo.2 <- (1-rbar^2)^2/(boot.ni-1)               # simple observed error variance for each study
          V.ec2 <- V.eo.2/A.compound.2^2                   # find the corrected error variance
          err.adj.2 <- 1/((UT^2-1)*boot.ri^2+1)            # find the error adjustment for indirect range restriction
          V.ve.2 <- sum(wi.2*V.ec2*err.adj.2^2)/sum(wi.2)  # find the error variance for the meta
          V.rho.2 <- V.rC.2-V.ve.2                         # find the random-effects variance component
          if (V.rho.2 < 0) {V.rho.2 <- 0}                  # if REVC is less than zero, set to zero
          SD.rho.2 <- sqrt(V.rho.2)                        # find SD rho
          lb2 <- rbarC.2 - 1.28*SD.rho.2            # find lower bound
          IndirectBts <- c(rbarC.2, V.rho.2, lb2)     # collect estimates to be examined for confidenc intervals
          return(IndirectBts)
        }# end bootstrap function for lower bound
        #Bootstrapped Statistic
        if (Boots) { #start bootstrap for indirect variance of rho
          SH.IndirectRR.bootstrap <- boot(mod.data, SH.IndirectRR.f, R = BootIter)
          SH.IndirectRR.bootstrap1 <- boot.ci(SH.IndirectRR.bootstrap, type="bca", index=1)
          Mbt.CI.LB <- SH.IndirectRR.bootstrap1$bca[4]
          Mbt.CI.UB <- SH.IndirectRR.bootstrap1$bca[5]
          SH.IndirectRR.bootstrap2 <- boot.ci(SH.IndirectRR.bootstrap, type="bca", index=2)
          REVC.CI.LB <- SH.IndirectRR.bootstrap2$bca[4]
          REVC.CI.UB <- SH.IndirectRR.bootstrap2$bca[5]
          SH.IndirectRR.bootstrap3 <- boot.ci(SH.IndirectRR.bootstrap, type="bca", index=3)
          CR.CI.LB <- SH.IndirectRR.bootstrap3$bca[4]
          CR.CI.UB <- SH.IndirectRR.bootstrap3$bca[5]
        } #end bootstrap for indirect lower bound
     }# end Indirect (SnH)
    }# End Schmidt-Hunter model computations
    ################################################################################
    if(MetaModel=='Morris') # begin Morris computations
      ################################################################################
    {
      if(!IndirectRR){ # start computing for Direct Range Restriction   
        ##########################################
        # Direct Range Restriction
        ##########################################
        nR <- ri*ni                            # weight r by N
        Nsum <- sum(ni)                        # find sum of N
        rbar <- sum(nR)/Nsum                   # find sample-weighted mean r
        V.obs <- sum(ni*(ri-rbar)^2)/(sum(ni)) # find weighted observed variance
        if(k_Correct){V.obs <- V.obs*(k/(k-1))}
        rxxa <- ((rxxi^.5/ui)/sqrt(1+rxxi*(ui^-2-1)))^2 # find reliability applicant sample (direct RR)
        Ux <- 1/ui                                # inverse of range restriction ratio
        rC <- (Ux*ri)/(sqrt(rxxa*(ryyi+Ux^2*ri^2-ri^2)))    # find corrected correlations
        A.compound <- ri/rC                      # find A, the attenuation factor
        V.eo <- (1-rbar^2)^2/(ni-1)              # simple observed error variance for each study
        V.ec1 <- V.eo/A.compound^2               # simple error variance of corrected correlations
        err.adj <- 1/((Ux^2-1)*ri^2+1)           # range restriction error variance correction
        V.ve <- V.ec1*err.adj^2                  # refined error variance of corrected correlations 
        morris.dat <- data.frame(cbind(rC,V.ve)) # collect  estimates
        morris1 <- rma(yi=rC,vi=V.ve,data=morris.dat,
                       control=list(maxiter=1000, stepadj=.5), method="REML")        # run the random-effects meta with REML
        Morris.M.rho <- morris1$b                              # output the mean
        rownames(Morris.M.rho) <- c()                          # strips value of intercept label
        colnames(Morris.M.rho) <- c("Morris.M.rho")
        Morris.CI95.L <- morris1$ci.lb                         # output the lower CI bound
        Morris.CI95.U <- morris1$ci.ub                         # output the upper CI bound
        Morris.V.rho <- morris1$tau2                           # random-effects variance component
        Morris.SD.rho <- sqrt(morris1$tau2)                    # output for Morris RHO sd
        Morris.CR80.L <- (Morris.M.rho - 1.28*Morris.SD.rho)   # Lower CR Bound
        Morris.CR80.U <- (Morris.M.rho + 1.28*Morris.SD.rho)   # Upper CR Bound
        colnames(Morris.CR80.L) <- c("Morris.CR80.L")
        
        # Common output
        rbar.a <- Morris.M.rho
        CI95.L <- Morris.CI95.L
        CI95.U <- Morris.CI95.U
        V.m <- ((morris1$ci.ub-morris1$ci.lb)/(2*1.96))^2      # compute the variance of the mean
        V.rho <- Morris.V.rho
        SD.rho <- Morris.SD.rho
        PctVE <- 100-morris1$I2
        CR80.L <- Morris.CR80.L
        CR80.U <- Morris.CR80.U
        REVC.CI.LB <- V.rho-1.96*morris1$se.tau2
        REVC.CI.UB <- V.rho+1.96*morris1$se.tau2
        #
        PI80.U <- NA
        PI80.L <- NA
        if (PredInt) {
          PIs <- Preds(rbar.a,V.rho,V.m,k)
        }# end PredInt
        # Bootstrap CI for lower bound
        Morris.DirectRR.f <- function(d, i){ # begin Morris bootstrap function
          d2 <- d[i,]
          boot.ri <- d2$ri
          boot.ri <- ifelse(boot.ri == 0, 0.0000001, boot.ri)   
          boot.ni <- d2$ni
          boot.rxxi <- d2$rxxi
          boot.ryyi <- d2$ryyi
          boot.u <- d2$ui
          nR <- boot.ri*boot.ni                         # weight r by N
          Nsum <- sum(boot.ni)                          # find sum of N
          rbar <- sum(nR)/Nsum                          # find sample-weighted mean r
          
          Ux <- 1/boot.u
          r.dis.rr <- 
            boot.ri*Ux/sqrt((Ux^2-1)*boot.ri^2+1)       # disattenuate for direct range restriction
          a1 <- boot.ri/r.dis.rr                        # find attenuation factor a1 for direct RR
          r.dis.rxx <- boot.ri/sqrt(rxxi)               # disattenuate for reliability of X (IV in selected sample)
          a2 <- boot.ri/r.dis.rxx                       # find attenuation factor a2 for reliability of X
          r.dis.ryy <- boot.ri/sqrt(boot.ryyi)          # disattenuate for reliabiliyt of Y (DV in selected sample)
          a3 <- boot.ri/r.dis.ryy                       # find attenuation factor for reliability of Y
          A.compound <- a1*a2*a3                        # find compound attenuation factor
          rC <- boot.ri/A.compound                      # find corrected correlations
          wi <- A.compound^2*boot.ni                    # find the weights for the meta
          rbarC.1 <- sum(wi*rC)/sum(wi)                 # find the mean corrected correlations
          V.rC <- sum(wi*(rC-rbarC.1)^2)/sum(wi)        # find the variance of the corrected correlations
          V.eo <- (1-rbar^2)^2/(boot.ni-1)              # simple observed error variance for each study
          V.ec1 <- V.eo/A.compound^2                    # simple error variance of corrected correlations
          err.adj <- 1/((Ux^2-1)*boot.ri^2+1)           # range restriction error variance correction
          V.ve <- V.ec1*err.adj^2                       # refined error variance of corrected correlations 
          morris.dat <- data.frame(cbind(rC,V.ve))      # collect  estimates
          morris1 <- rma(yi=rC,vi=V.ve,data=morris.dat,
                         control=list(maxiter=1000, stepadj=.5), method="REML")        # run the random-effects meta with REML
          Morris.M.rho <- morris1$b                     # output the mean
          rownames(Morris.M.rho) <- c()
          colnames(Morris.M.rho) <- c("Morris.M.rho")
          Morris.CI95.L <- morris1$ci.lb                 # output the lower CI bound
          Morris.CI95.U <- morris1$ci.ub                 # output the upper CI bound
          Morris.V.rho <- morris1$tau2                   # random-effects variance component
          Morris.SD.rho <- sqrt(morris1$tau2)            # output for Morris RHO sd
          lb3 <- Morris.M.rho - 1.28*Morris.SD.rho            # find lower bound
          MorrisDirect <- c(Morris.M.rho, Morris.V.rho, lb3)     # collect estimates to be examined for confidenc intervals
          return(MorrisDirect)
        } # end function for bootstrap
        #Bootstrapped Statistic
        if (Boots) { #begin boostrap for direct RR for Morris
          Morris.Dir.bootstrap <- boot(mod.data, Morris.DirectRR.f, R = BootIter)
          Morris.Dir.bootstrap1 <- boot.ci(Morris.Dir.bootstrap, type="bca", index=1)
          Mbt.CI.LB <- Morris.Dir.bootstrap1$bca[4]
          Mbt.CI.UB <- Morris.Dir.bootstrap1$bca[5]
          Morris.Dir.bootstrap2 <- boot.ci(Morris.Dir.bootstrap, type="bca", index=2)
          REVC.CI.LB <- Morris.Dir.bootstrap2$bca[4]
          REVC.CI.UB <- Morris.Dir.bootstrap2$bca[5]
          Morris.Dir.bootstrap3 <- boot.ci(Morris.Dir.bootstrap, type="bca", index=3)
          CR.CI.LB <- Morris.Dir.bootstrap3$bca[4]
          CR.CI.UB <- Morris.Dir.bootstrap3$bca[5]
        } #end boostrap for direct RR for Morris
      } # end Direct Range Restriction
      if(IndirectRR){ # Start computations for indirect range restriction
        #######################################################################
        # Indirect Range Restriction
        ########################################################################
        nR <- ri*ni                                   # weight r by N
        Nsum <- sum(ni)                               # find sum of N
        rbar <- sum(nR)/Nsum                          # find sample-weighted mean r
        V.obs <- sum(ni*(ri-rbar)^2)/(sum(ni))        # find weighted observed variance
        RXXa <-1-ui^2*(1-rxxi)                          # find the reliability of the unrestricted sample
        UT <- 1/sqrt((ui^2-(1-RXXa))/RXXa)              # find UT, disattenuation for indirect RR
        rdiss.rel <- ri/sqrt(rxxi*ryyi)                 # find correlation disattenuated for reliability in X and Y
        rC.2 <- 
          (rdiss.rel*UT)/(sqrt((UT^2-1)*rdiss.rel^2+1)) # find corrected correlations
        A.compound.2 <- ri/rC.2                         # find compound correction factor
        wi.2 <- ni*A.compound.2^2                       # find the weights
        rbarC.2 <- sum(wi.2*rC.2)/sum(wi.2)             # find the mean corrected correlation
        V.rC.2 <- sum(wi.2*(rC.2-rbarC.2)^2)/sum(wi.2)  # find the corrected total variance
        V.eo.2 <- (1-rbar^2)^2/(ni-1)                   # simple observed error variance for each study
        V.ec2 <- V.eo.2/A.compound.2^2                  # find the corrected error variance
        err.adj.2 <- 1/((UT^2-1)*ri^2+1)                # find the error adjustment for indirect range restriction
        V.ve <- V.ec2*err.adj.2^2                       # refined error variance of corrected correlations  
        morris.dat <- data.frame(cbind(rC.2,V.ve))      # collect  estimates
        morris1 <- rma(yi=rC.2,vi=V.ve,data=morris.dat,
                       control=list(maxiter=1000, stepadj=.5), method="REML")        # run the random-effects meta with REML
        Morris.M.rho <- morris1$b                       # output the mean
        rownames(Morris.M.rho) <- c()
        colnames(Morris.M.rho) <- c("Morris.M.rho")
        Morris.CI95.L <- morris1$ci.lb                  # output the lower CI bound
        Morris.CI95.U <- morris1$ci.ub                  # output the upper CI bound
        Morris.V.rho <- morris1$tau2
        Morris.SD.rho <- sqrt(morris1$tau2)             # output for Morris RHO sd
        Morris.CR80.L <- (Morris.M.rho - 1.28*Morris.SD.rho)   # Lower CR Bound
        Morris.CR80.U <- (Morris.M.rho + 1.28*Morris.SD.rho)   # Upper CR Bound
        colnames(Morris.CR80.L) <- c("Morris.CR80.L")
        # Common output
        rbar.a <- Morris.M.rho
        CI95.L <- Morris.CI95.L
        CI95.U <- Morris.CI95.U
        V.rho <- Morris.V.rho
        SD.rho <- Morris.SD.rho
        V.m <- ((morris1$ci.ub-morris1$ci.lb)/(2*1.96))^2 # find the variance of the mean
        PctVE <- 100-morris1$I2
        CR80.L <- Morris.CR80.L
        CR80.U <- Morris.CR80.U
        REVC.CI.LB <- V.rho-1.96*morris1$se.tau2
        REVC.CI.UB <- V.rho+1.96*morris1$se.tau2
        #
        PI80.U <- NA
        PI80.L <- NA
        if (PredInt) {
          PIs <- Preds(rbar.a,V.rho,V.m,k)
        }# end PredInt 
        # Bootstrap CI for lower bound
        Morris.IndirectRR.f <- function(d, i){ # begin Morris bootstrap function
          d2 <- d[i,]
          boot.ri <- d2$ri
          boot.ri <- ifelse(boot.ri == 0, 0.0000001, boot.ri)   # avoid dividing by zero later on
          boot.ni <- d2$ni
          boot.rxxi <- d2$rxxi
          boot.ryyi <- d2$ryyi
          boot.u <- d2$ui
          nR <- boot.ri*boot.ni
          Nsum <- sum(boot.ni)
          rbar <- sum(nR)/Nsum
          RXXa <-1-boot.u^2*(1-boot.rxxi)                 # find the reliability of the unrestricted sample
          UT <- 1/sqrt((boot.u^2-(1-RXXa))/RXXa)          # find UT, disattenuation for indirect RR
          rdiss.rel <- boot.ri/sqrt(boot.rxxi*boot.ryyi)  # find correlation disattenuated for reliability in X and Y
          rC.2 <- 
            (rdiss.rel*UT)/(sqrt((UT^2-1)*rdiss.rel^2+1)) # find corrected correlations
          A.compound.2 <- boot.ri/rC.2                    # find compound correction factor
          wi.2 <- boot.ni*A.compound.2^2                  # find the weights
          rbarC.2 <- sum(wi.2*rC.2)/sum(wi.2)             # find the mean corrected correlation
          V.rC.2 <- sum(wi.2*(rC.2-rbarC.2)^2)/sum(wi.2)  # find the corrected total variance
          V.eo.2 <- (1-rbar^2)^2/(boot.ni-1)              # simple observed error variance for each study
          V.ec2 <- V.eo.2/A.compound.2^2                  # find the corrected error variance
          err.adj.2 <- 1/((UT^2-1)*boot.ri^2+1)           # find the error adjustment for indirect range restriction
          V.ve <- V.ec2*err.adj.2^2                       # refined error variance of corrected correlations  
          morris.dat <- data.frame(cbind(rC.2,V.ve))      # collect  estimates
          morris1 <- rma(yi=rC.2,vi=V.ve,data=morris.dat,
                         control=list(maxiter=1000, stepadj=.5), method="REML")        # run the random-effects meta with REML
          Morris.M.rho <- morris1$b                       # output the mean
          rownames(Morris.M.rho) <- c()
          colnames(Morris.M.rho) <- c("Morris.M.rho")
          Morris.CI95.L <- morris1$ci.lb                  # output the lower CI bound
          Morris.CI95.U <- morris1$ci.ub                  # output the upper CI bound
          Morris.V.rho <- morris1$tau2
          Morris.SD.rho <- sqrt(morris1$tau2)             # output for Morris RHO sd
          lb4 <- Morris.M.rho - 1.28*Morris.SD.rho            # find lower bound
          MorrisIndirect <- c(Morris.M.rho, Morris.V.rho, lb4)     # collect estimates to be examined for confidenc intervals
          return(MorrisIndirect)
        } # end Morris indirect bootstrap
        if (Boots) { # begin indirect RR bootstrap for Morris
          Morris.Indir.bootstrap <- boot(mod.data, Morris.IndirectRR.f, R = BootIter)
          Morris.Indir.bootstrap1 <- boot.ci(Morris.Indir.bootstrap, type="bca", index=1)
          Mbt.CI.LB <- Morris.Indir.bootstrap1$bca[4]
          Mbt.CI.UB <- Morris.Indir.bootstrap1$bca[5]
          Morris.Indir.bootstrap2 <- boot.ci(Morris.Indir.bootstrap, type="bca", index=2)
          REVC.CI.LB <- Morris.Indir.bootstrap2$bca[4]
          REVC.CI.UB <- Morris.Indir.bootstrap2$bca[5]
          Morris.Indir.bootstrap3 <- boot.ci(Morris.Indir.bootstrap, type="bca", index=3)
          CR.CI.LB <- Morris.Indir.bootstrap3$bca[4]
          CR.CI.UB <- Morris.Indir.bootstrap3$bca[5]
        } # end indirect RR bootstrap for Morris
      } #end IndirectRR for Morris
    } # End Morris either way
    
    Results <- cbind(missing.cases, k, Nsum, rbar, rbar.a, V.obs,        # 1 - 6
                     CI95.L, CI95.U, V.rho, REVC.CI.LB, REVC.CI.UB,      # 7 - 11
                     SD.rho, PctVE,                                      # 12, 13
                     CR80.L, CR80.U, CR.CI.LB, CR.CI.UB, PIs[1], PIs[2], # 14-19
                     Mbt.CI.LB, Mbt.CI.UB, sqrt(V.obs))                 # 20, 21, 22
    colnames(Results) <- c("Missing Cases","Studies", "Total N", "BBMean","Mean", "BBVar", "95%CI_Lo", "95%CI_Hi", "Vrho", 
                           "V95lo", "V95Hi", "SD(rho)", "PctVE", "80%CR_Lo",
                           "80%CR_Hi","CRlow.CI.Low", "CRlow.CI.Hi", "80%PI_Lo", "80%PI_Hi", "BtM95Lo", "BtM95Hi","SD(r)")
    MetaM <- I(MetaModel)
    Results <-data.frame(Results, MetaM)
    main.title <- "Results of Psychometric Meta-analysis"
    second.title <- "Descriptives"
    third.title <- "Model"
    nada <- " "
    text1 <- "Number of Excluded Studies"; Output.missing <- sprintf("%5.0f", Results[1])
    text2 <- "Number of Correlations"; Outuput.k <- sprintf("%5.0f", Results[2])
    text3 <- "Total Sample Size"; Output.sumN <- sprintf("%5.0f", Results[3])
    text4 <- "Observed (bare-bones) Weighted Mean"; Output.M <-sprintf("%5.3f", Results[4])
    text5 <- "Observed (bare-bones) Weighted Variance";Output.Vobs <-sprintf("%5.3f", Results[6])
    text6 <- "Observed (bare-bones) SD-r"; Output.SDobs <- sprintf("%5.3f", Results[22])
    text7 <- "Model"; Output.Model <- sprintf("%6s", Results[23])
    text8 <- "Corrected Mean"; Output.rbar.a <- sprintf("%5.3f", Results[5])
    text9 <- "95 % CI Lower Bound Corrected Mean"; Output.CILB <- sprintf("%5.3f",Results[7])
    text10 <- "95 % CI Upper Bound Corrected Mean"; Output.CIUB <- sprintf("%5.3f",Results[8])
    text11 <- "Bootstrap 95CI Lower Bound Corrected Mean"; Output.BtCILB <- sprintf("%5.3f",Results[20])
    text12 <- "Bootstrap 95CI Upper Bound Corrected Mean"; Output.BtCIUB <- sprintf("%5.3f",Results[21])
    text13 <- "Variance of Rho"; Output.Vrho <- sprintf("%5.3f", Results[9])
    text14 <- "Bootstrap 95 % CI Lower Bound for V(Rho)"; Output.Vrho95Low <- sprintf("%5.3f", Results[10])
    text15 <- "Bootstrap 95 % CI Upper Bound for V(Rho)"; Output.Vrho95up <- sprintf("%5.3f", Results[11])
    text16 <- "Standard Deviation of Rho"; Output.SDrho <- sprintf("%5.3f",Results[12])
    text17 <- "Proportion Variance due to Artifacts"; Output.PctVE <- sprintf("%5.3f",Results[13])
    text18 <- "80 % Credibility Interval Lower Bound"; Output.CR80Lo <- sprintf("%5.3f",Results[14])
    text19 <- "80 % Credibility Interval Upper Bound"; Output.CR80Hi <- sprintf("%5.3f",Results[15])
    text20 <- "80 % Prediction Interval Lower Bound"; Output.PI80Lo <- sprintf("%5.3f",Results[18])
    text21 <- "80 % Prediction Interval Upper Bound"; Output.PI80Hi <- sprintf("%5.3f",Results[19])
    text22 <- "95 % CI.LB for Credibility Lower Bound"; Output.CRCILB <- sprintf("%5.3f",Results[16])
    text23 <- "95 % CI.UB for Credibility Lower Bound"; Output.CRCIUB <- sprintf("%5.3f",Results[17])
    Desc.labels <- c(nada, second.title, text1, text2, text3, text4,
                     text5, text6, nada, text7, text8, text9, text10, text11, text12, text13,
                     text14, text15, text16, text17, text18, text19, text20, text21,text22, text23)
    #str(Desc.labels)
    Desc.data[,zz] <- c(nada, nada, Output.missing, Outuput.k, Output.sumN, Output.M, Output.Vobs,
                        Output.SDobs, nada, Output.Model, Output.rbar.a, Output.CILB, Output.CIUB,
                        Output.BtCILB, Output.BtCIUB,
                        Output.Vrho, Output.Vrho95Low, Output.Vrho95up, Output.SDrho, 
                        Output.PctVE, Output.CR80Lo, Output.CR80Hi, Output.PI80Lo, Output.PI80Hi,
                        Output.CRCILB, Output.CRCIUB)
  } # end main loop
  
  Desc.out <- data.frame(Desc.labels,Desc.data)
  if (zz==1) levels(moderator) <- 'All studies'
  colnames(Desc.out) <- c('Meta-Analysis Results', levels(moderator))
  return(Desc.out)
} # End function Schmidt
#################################################################################
###########################################
# S&H 2015 Table 3.5
BBones <-Schmidt(ri=main.data$ri, ni=main.data$Nr, Bias_Correct = FALSE,
                MetaModel='SnH', PredInt=TRUE, Boots=TRUE, k_Correct = FALSE,
                BootIter=5000, Missing=TRUE)
BBones
Case4 <-Schmidt(ri=main.data$ri, ni=main.data$Nr, ui=main.data$RRuX,
        rxxi = main.data$RXXi, ryyi = main.data$Ryyi,
        Bias_Correct = FALSE, k_Correct = FALSE,
        IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=FALSE,
        BootIter=5000, Missing=TRUE)
Case4
Case2 <-Schmidt(ri=main.data$ri, ni=main.data$Nr,ui=main.data$RRuX,
                rxxi = main.data$RXXi, ryyi = main.data$Ryyi,
                Bias_Correct = FALSE, k_Correct = FALSE,
                IndirectRR=FALSE, MetaModel='SnH', PredInt=TRUE, Boots=FALSE,
                BootIter=5000, Missing=TRUE)
Case2
rel.Y.only <-Schmidt(ri=main.data$ri, ni=main.data$Nr,
                 ryyi = main.data$Ryyi,
                 Bias_Correct = TRUE, k_Correct = FALSE,
                IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=TRUE,
                BootIter=5000, Missing=TRUE)
rel.Y.only
###########################################

# Gonzalez-Mule
main.data$ui <- rep(.63,length(main.data$r))
Schmidt(ri=main.data$r, ni=main.data$N,ui=main.data$ui,
        rxxi = main.data$rxx, ryyi = main.data$ryy,
        IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=TRUE,
        BootIter=2000, Missing=TRUE)
###########################################

# Breuer
Schmidt(ri=main.data$r, ni=main.data$N,
        rxxi = main.data$rxx, ryyi = main.data$ryy,
        IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=TRUE,
        BootIter=2000, Missing=TRUE)
Breuer2 <- subset(main.data, main.data$Team.effectivenes =='TP' | main.data$Team.effectivenes=='S')
str(Breuer2)
Schmidt(ri=Breuer2$r, ni=Breuer2$N,
        rxxi = Breuer2$rxx, ryyi = Breuer2$ryyi, moderator = Breuer2$Team.effectiveness,
        IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=TRUE,
        BootIter=2000, Missing=FALSE)
Breuer3 <- subset(main.data,main.dat
                  a$Team.effectivenes=='S')
str(Breuer3)
Schmidt(ri=Breuer3$r, ni=Breuer3$N,
        rxxi = Breuer3$rxx, ryyi = Breuer3$ryyi, 
        IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=TRUE,
        BootIter=2000, Missing=FALSE)
###########################################

# Choi Sheet1
Schmidt(ri=main.data$r, ni=main.data$N,
        rxxi = main.data$rxx, ryyi = main.data$ryy, moderator=main.data$Pub,
        IndirectRR=TRUE, MetaModel='SnH', PredInt=TRUE, Boots=TRUE,
        BootIter=500, Missing=FALSE)
###########################################

# Rabl Sheet1
Schmidt(ri=main.data$r, ni=main.data$N,
        IndirectRR=FALSE, MetaModel='Morris', PredInt=TRUE, Boots=TRUE,
        BootIter=500, Missing=TRUE)
###########################################

#Morris Sheet1
Schmidt(ri=main.data$Validity, ni=main.data$N, ryyi = main.data$CriterionReliability,
        IndirectRR=FALSE, MetaModel='Morris', PredInt=TRUE, Boots=FALSE,
        BootIter=500, Missing=TRUE)