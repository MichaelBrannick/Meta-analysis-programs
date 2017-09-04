# If not installed, install the package e1071
# Select a CRAN Mirror when a dialog box opens
#install.packages('e1071')

# Run this line just one in an R session
#library(e1071)

# sampling distribution of r, the correlation coefficient
samplesize <- 120

#this is the sample size for each correlation
size.r <- .5

#this is the parameter, rho
numbersamples <- 10000

# this is the number of samples in the empirical sampling distribution
# lines above are setup to run the simulation
out1 <- 1:numbersamples
for (i in 1:numbersamples){
theta <- rnorm(samplesize,0,1)
e1 <- rnorm(samplesize,0,1)
e2 <- rnorm(samplesize,0,1)
weight <- sqrt(size.r)
x <- weight*theta+e1*sqrt(1-size.r)
y <- weight*theta+e2*sqrt(1-size.r)
out1[i] <- cor(x,y)
}

# out1 has the sampled correlatins
summary(out1)
sd(out1)
hist(out1)
skewness(out1)

# You can also screen skewness with a boxplot() function.
#boxplot(out1)