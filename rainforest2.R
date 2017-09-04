##########################################################################################################
##########################################################################################################
## versions of forest forest plots: conventional forest plot, thick forest plot, rainforest plot
##########################################################################################################
##########################################################################################################
library(metafor) #calculate MA
library(ggplot2) #plot graphs
library(scales) #for axis transformation to log scale 
library(gridExtra) #grid.arrange()
library(RColorBrewer) #for blues
library(plyr) #mdply()

##########################################################################################################
## ggplot2: new theme for forest plots
##########################################################################################################
blacktheme <- theme_set(theme_grey())
theme_set(blacktheme)
blacktheme <-theme_update(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = element_text(family="sans", color="black", size=15),
                          axis.text.y = element_text(family="sans", color="black", size=15, vjust=1),
                          axis.ticks.x = element_line(color="black"),
                          axis.ticks.y = element_blank())
blacktheme <- theme_set(blacktheme)

blues2 <- brewer.pal(9, "Blues")

##########################################################################################################
##########################################################################################################
## simulated MA, see Peters 2006, Bax 2009
##########################################################################################################
##########################################################################################################
set.seed(81124)
baseline <- runif(10, 0.2,0.4)  #10 studies with baseline risk between 0.2 and 0.4
OR <- rlnorm(10, meanlog=0, sdlog=0.17) #sample OR
size <- rlnorm(10, 5, 0.7) #size of study arms

treatment <- baseline/OR #risk in treatment group

#calculate cell frequencies
d <- baseline*size
b <- treatment*size
a <- size-b
c <- size-d

#calculate MA with fixed effect
meta <- rma(ai=a, bi=b, ci=c, di=d, measure="OR", method="FE")
meta

################
##make dataframe for plotting
################
dat <- data.frame(a,b,c,d) #cell frequencies
dat <- escalc(data=dat, measure="OR", ai=a, bi=b, ci=c, di=d, append=TRUE) #yi=log OR, vi=log variance
dat$lower <- (dat$yi - 1.96*sqrt(dat$vi)) #lower CI
dat$upper <- (dat$yi + 1.96*sqrt(dat$vi)) #upper CI

#empty row to force space between last study and summary
empty <- data.frame(a=c(NA, NA), b=c(NA, NA), c=c(NA,NA), d=c(NA,NA),
                    yi=c(NA,NA), vi=c(NA,NA),lower=c(NA,NA), upper=c(NA,NA))
dat <- rbind(dat, empty) # bind together
#add labels
dat$study <- factor(c(1:10, "", "Summary"), levels=c("Summary", "", 10:1), 
                    labels=c("Summary", "", rev(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "k"))))

dat$weight <- c(weights(meta), NA, NA) #add weights and NA for empty rows

################
##add info for summary diamond and effect line
################

## summary estimates for diamond and effect line
dat$sum.y <- c(1, 0.7, 1, 1.3, NA, NA, NA, NA, NA, NA, NA, NA)
dat$sum.x <- c(meta$ci.lb, meta$b, meta$ci.ub, meta$b, NA, NA, NA, NA, NA, NA, NA, NA)
dat$line.y <- c(1.3, 12.3, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
dat$line.x <- c(meta$b, meta$b, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

##########################################################################################################
## conventional forest plot
##########################################################################################################

conv <- ggplot(data=dat, aes(y=study, x=yi, xmin=lower, xmax=upper))+
  geom_vline(aes(x=1), linetype=2, xintercept = 0)+
  geom_point(aes(size=weight), shape=15, col="#08306B")+
  geom_errorbarh(height=0, col="#08306B")+
  geom_polygon(aes(x=sum.x, y=sum.y), color="#08306B", fill="#08306B")+
  geom_line(aes(y=line.y, x=line.x), linetype=3)+
  scale_x_continuous(name="",breaks=trans_breaks("exp", function(x) log(x), n=5), 
                     labels=trans_format("exp", math_format(expr=.x)))+
  blacktheme+
  theme(legend.position="none")
conv

##########################################################################################################
## thick forest plot
##########################################################################################################

thick <- ggplot(data=dat, aes(y=study, x=yi, xmin=lower, xmax=upper))+
  geom_vline(aes(x=1), linetype=2, xintercept = 0)+
  geom_errorbarh(aes(size=weight),height=0, shape=15, show_guide=FALSE, , col="#08306B")+
  geom_polygon(aes(x=sum.x, y=sum.y), color="#08306B", fill="#08306B")+
  geom_line(aes(x=line.x, y=line.y), linetype=3)+
  geom_point(shape="I", show_guide=FALSE, color="darkred", size=8) +
  blacktheme+
  scale_x_continuous(name="",breaks=trans_breaks("exp", function(x) log(x), n=5), 
                     labels=trans_format("exp", math_format(expr=.x)))
thick


##########################################################################################################
## rainforest
##########################################################################################################
################
##calculate densities
################

rain <- function(mean, sd, lower, upper) {
  x <- seq(lower, upper, length=500)
  dens <- dnorm(x, mean=mean, sd=sd, log=TRUE)
  dens0 <- dens -min(dens)
  return(data.frame(dens0, x))
}


dat2 <- data.frame("mean"=dat$yi, "sd"=sqrt(dat$vi), "lower"=dat$lower, "upper"=dat$upper)
dat2 <- dat2[1:10,]

################
##apply function over dataframe, make dataframe for plotting
################
raindat <- mdply(dat2, rain)

# add weights
raindat$weight <- rep(c(weights(meta)), each=500)

# compute values above and below, factor in weight
raindat$order <- rev(rep(c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), each=500)) #force ordering of studies
raindat$above.weight <- raindat$dens0*raindat$weight/20+raindat$order
raindat$below.weight <- raindat$order-raindat$dens0*raindat$weight/20

#summary est
#raindat$sum.y <- c(rep(NA, 4000), -10, -8.5, -10, -11.5, rep(NA, 996))
raindat$sum.y <- c(rep(NA, 4000), -11, -9.5, -11, -12.5, rep(NA, 996))
raindat$sum.x <- c(rep(NA, 4000), meta$ci.lb, meta$b, meta$ci.ub, meta$b, rep(NA, 996))
raindat$line.y <- c(-10, 45, rep(NA, 4998))
raindat$line.x <- c(meta$b, meta$b, rep(NA, 4998))


rainforest <- ggplot(raindat, aes(x=x, xend=x, y=above.weight, yend=order, alpha=weight, color=dens0))+
  geom_vline(aes(x=1), color="black", linetype=2, xintercept = 0)+
  geom_line(aes(x=line.x, y=line.y), col="black", linetype=3, alpha=1)+
  stat_identity(geom="segment")+
  stat_identity(geom="segment", aes(y=below.weight))+
  geom_point(aes(x=mean, y=order), shape="I", col="white", size=5)+
  geom_segment(aes(y=order+0.05, yend=order-0.05, x=x, xend=x), col="white")+
  scale_colour_gradientn(colours=blues2, guide="none")+
  geom_polygon(aes(x=sum.x, y=sum.y), color="#08306B", fill="#08306B", alpha=1)+
  ylim(c(-15, 45))+
  scale_y_continuous(breaks=c(-10, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45),
                     labels=c("Summary", rev(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "k"))))+
  scale_x_continuous(name="",breaks=trans_breaks("exp", function(x) log(x), n=5), 
                     labels=trans_format("exp", math_format(expr=.x)))+ ##needs library(scales)
  guides(alpha=FALSE)+
  blacktheme
rainforest


grid.arrange(thick, conv, rainforest,  nrow=1)
