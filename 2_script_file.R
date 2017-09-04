# Function to install packages is install.packages()
install.packages('metafor')


# It will ask for a download location.  
# Pick the closest one. Turkey, Greece, Austria
# would be fine. 


# Import metafor package to the working environment. 
library('metafor')


# Change "dat.bcg" with the full path that you 
# acquired. 
data("dat.bcg", package = "metafor")


# Print data set to the screen. You do not want to 
# see row names, hence 'row.names=FALSE' is in the 
# print() function. 
print(dat.bcg, row.names=FALSE)


# Model fitting with rma() function. You have to assing 
# output of the rma() function or it is printed to the 
# screen and lost. 
res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg,
measure="RR", data=dat.bcg)


# By assigning output to a name you can be able to 
# reach and see it again. 
res


# Calculating confidence intervals with confint() function.  
confint(res)


# Creating forest plots with forest() function. 
forest(res)