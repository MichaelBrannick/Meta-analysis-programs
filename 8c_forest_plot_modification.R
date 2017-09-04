# Load the required packages for the analyses. 
library(metafor)
library(xlsx)


myers2 <- read.xlsx('MyersCrowther2009.xls', sheetIndex=2)


# Practice with another data set from the same file. 
# sheetIndex= let you to select a different tab from 
# within a same Excel file. 
myers3 <- read.xlsx('MyersCrowther2009.xls', sheetIndex=3)

# You can screen variables from a data set with 
# str() function. 
str(myers2)
str(myers3)
str(myers3$pathology)

res2 <- rma(yi=d, vi=v, data=myers2, method='DL')
res2

res3 <- rma(yi=d, vi=v, data=myers3, method='DL')
res3

forest(res2)

# To modify the x and y limits of the graph 
# you have to get the xlim and ylim of the default
# graph. 
defaults <- forest(res2)
defaults$xlim
defaults$ylim


# Study names has to be created like this. 
# c() function is used to create a vector of names. 
# You can change these names as you like. 
study_names <- c('Ali (2004)', 'Veli (1999)', 
                 'Kırkdokuz (2003)', 'Elli (1950)', 
                 'Metin (1990)', 'Ali (1991)', 
                 'Feyyaz (1992)', 'Rıza (1994)', 
                 'Şifo (1996)', 'Recep (1997)', 
                 'De Niro (1999)', 'Al Pacino (2004)', 
                 'Di Caprio (2014)', 'Premuzic (2012)', 
                 'Judge (2014)', 'Mount (2001)', 
                 'Jackson (2001)', 'Sandberg (2005)', 
                 'Obama (2009)', 'Bush (2014)', 
                 'Gültaş (2004)', 'Middler (2000)')


# You can order by using: order=
forest(res2, order='fit')
forest(res2, order='obs')


# You can feed study names using: slab=
# You have to feed some kind of vector of strings
# as stated above. 
forest(res2, 
       order='obs', 
       slab=study_names[1:20])


# Add a main title to your figure. 
forest(res2, 
       order='obs', 
       main='My Ground Shaking Study', 
       slab=study_names[1:20])

# Change your horizontal dimension of your figure. 
# Use argument xlim=
forest(res2, 
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-3, 6), 
       slab=study_names[1:20])

# Change margins of your figure. 
par(mar=c(2,4,2,2))

# Add summary effect text with
# mlab=
forest(res2, 
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-3, 6), 
       mlab='RE Model for All Studies', 
       slab=study_names[1:20])

# Modify x axis with 
# xlim=
forest(res2, 
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-4, 6), 
       mlab='RE Model for All Studies', 
       slab=study_names[1:20])

# Append/add text to already printed 
# forest plot. 
forest(res2, 
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-4, 6), 
       mlab='RE Model for All Studies', 
       slab=study_names[1:20])
text(-2.5, 21.8, 'Author(s) and Year')
text(c(3.6, 4.9), 21.8, c('ESs', '95% CI'))

# Change printed x values of the figure with 
# at=
forest(res2, 
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-4, 6), 
       mlab='RE Model for All Studies', 
       at=c(-1, 0, 1.5, 3), 
       slab=study_names[1:20])
text(-2.5, 21.8, 'Author(s) and Year')
text(c(3.6, 4.9), 21.8, c('ESs', '95% CI'))

# Change the label of the x axis
forest(res2,
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-4, 6), 
       mlab='RE Model for All Studies', 
       xlab='Standardized Mean Difference', 
       at=c(-1, 0, 1.5, 3), 
       slab=study_names[1:20])
text(-2.5, 21.8, 'Author(s) and Year')
text(c(3.6, 4.9), 21.8, c('ESs', '95% CI'))

res2.patho <- rma(yi=d, 
                  vi=v, 
                  data=myers2, 
                  subset=(pathology==1),
                  method='DL')

addpoly(res2.patho, 
        mlab='Pathology Moderator', 
        cex=1)

# Change the label of the x axis
forest(res2,
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-4, 6), 
       mlab='RE Model for All Studies', 
       xlab='Standardized Mean Difference', 
       at=c(-1, 0, 1.5, 3), 
       slab=study_names[1:20])
text(-2.5, 21.8, 'Author(s) and Year')
text(c(3.6, 4.9), 21.8, c('ESs', '95% CI'))

# Change the y limits of the graph. 
forest(res2,
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-4, 6), 
       ylim=c(-1.5, 22.8), 
       mlab='RE Model for All Studies', 
       xlab='Standardized Mean Difference', 
       at=c(-1, 0, 1.5, 3), 
       slab=study_names[1:20])
text(-2.5, 21.8, 'Author(s) and Year')
text(c(3.6, 4.9), 21.8, c('ESs', '95% CI'))

# Change the label of the x axis
forest(res3,
       order='obs', 
       main='My Ground Shaking Study', 
       xlim=c(-6, 11), 
       # In order to create room for 
       # plotting moderators give 1st value of ylim
       # some negative numbers. 
       ylim=c(-4, 20), 
       mlab='RE Model', 
       xlab='Standardized Mean Difference', 
       at=c(-1, 0, 1.5, 3, 5, 7), 
       slab=study_names[1:17])

text(-3.5, 18.8, 'Author(s) and Year')
text(c(7, 9.3), 18.8, c('ESs', '95% CI'))
res3.patho1 <- rma(yi=d, vi=v, data=myers3, 
                  method='DL', 
                  subset=(pathology==1))
addpoly.rma(res3.patho1, 
            row=-2, 
            mlab='Proud Moderator 1', 
            cex=1)

res3.patho2 <- rma(yi=d, vi=v, data=myers3, 
                  method='DL', 
                  subset=(pathology==2))
addpoly.rma(res3.patho2, 
            row=-3,
            mlab='Another Moderator 2', 
            cex=1)

res3.patho3 <- rma(yi=d, vi=v, data=myers3, 
                  method='DL', 
                  subset=(pathology==3))
addpoly.rma(res3.patho3, 
            row=-4,
            mlab='Final Moderator', 
            cex=1)

