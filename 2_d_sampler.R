# If not installed, install package e1071
# install.packages('e1071')

# Run this line just once and then comment out.
# library(e1071)

# sampling distribution of d, the standardized mean difference
n1 <- 15.5

#number of people in group 1
n2 <- 15

#number of people in group 2
SD1 <- 1

#within group standard deviation for group 1
SD2 <- 1

#within group standard deviation for group 2
M1 <- 15

#mean for group 1
M2 <- 14

#mean for group 2
numbersamples <- 1000

# this is the number of samples in the empirical sampling
# program computes the rest
df1 = n1-1
df2 = n2-1
SDpool <- sqrt((df1*SD1^2+df2*SD2^2)/(df1+df2))
delta <- (M1-M2)/SDpool
delta

#this is the parameter, delta, the standardized mean difference
# lines above are setup to run the simulation
out1 <- 1:numbersamples
for (i in 1:numbersamples){
sample1 <- rnorm(n1,M1,SD1)
sample2 <- rnorm(n2,M2,SD2)
Ms1 = mean(sample1)
Vs1 =var(sample1)
Ms2 = mean(sample2)
Vs2 = var(sample2)
d = (Ms1-Ms2)/sqrt((Vs1*df1+Vs2*df2)/(df1+df2))
out1[i] <- d
}

# out1 has the sampled values of d
summary(out1)
sd(out1)
hist(out1)
skewness(out1)


