
library(metaplotr)
help("metaplotr")
help("crosshairs")
#crosshairs(x, y, xse, yse, x_lab = NULL, y_lab = NULL, main_lab = NULL,
#           confint = 0.95, mdrtr = NULL, mdrtr_lab = NULL, mdrtr_lab_pos = NULL,
#           lab_size = 14, pnt_size = 3, whis_on = TRUE, annotate = FALSE,
#           grid_dense = FALSE, bxplts = TRUE)

# Remove all variables in the .GlobalEnv, effectively clearing .GlobalEvn
rm(list = ls())

# help("FergusonBrannick2012")

# attach data frame to working environment.

attach(FergusonBrannick2012)

crosshairs(pub_z, dis_z, pub_z_se, dis_z_se)
#
# confint option can control whiskers length.
# crosshairs(pub_z, dis_z, pub_z_se, dis_z_se, confint = .95)
 crosshairs(pub_z, dis_z, pub_z_se, dis_z_se, confint = .7)
#crosshairs(pub_z, dis_z, pub_z_se, dis_z_se, confint = .3)


# whis_on option opens and closes whiskers.
crosshairs(pub_z, dis_z, pub_z_se, dis_z_se, whis_on = FALSE)


#
# Main and axes labels can be changed.
crosshairs(pub_z, dis_z, pub_z_se, dis_z_se,
           main_lab = 'Published vs. Dissertation Effect Sizes', 
           x_lab = 'Published Studides',
           y_lab = 'Dissertations')

# Annotated correlation and mean values can be added to the graph.

attach(Sweeney2015)
# help("Sweeney2015")

# add descriptive statistics to graph
crosshairs(inten_d, beh_d, inten_se, beh_se,
           main_lab = 'Sweeney (2015) Data', x_lab = 'Intentions',
           y_lab = 'Behaviors',annotate = TRUE)

# Boxplots can be hidden.
crosshairs(inten_d, beh_d, inten_se, beh_se,
           main_lab = 'Sweeney (2015) Data', x_lab = 'Intentions',
           y_lab = 'Behaviors',annotate = TRUE,
           bxplts = FALSE)


# Add moderator and label
attach(GenderDiff02)
#help("GenderDiff02")
crosshairs(men_z, women_z, men_se, women_se,
           main_lab = 'Ali et al. Psychopathology and Parental Acceptance',
           x_lab='Men', y_lab='Women', mdrtr = region, mdrtr_lab = 'Region',
           mdrtr_lab_pos = c(.1,.5))

#
attach(McLeod2007) #McLeod2007
library(metafor)
res1 <- rma(yi=z, vi=var, method = "DL", data = McLeod2007)
res2 <- blup(res1)
res2
# Assign data to x, standard error of x, y, standard error of y,
# variable name of a moderator (if any) here. Note how the names
# and values of the x variables came from the McLead2007 dataset.
# The names and values of the shrunken estimates came from 
# the output of the metafor program.
#
x1 <- McLeod2007$z
se.x1 <- sqrt(McLeod2007$var)
y1 <- res2$pred
se.y1 <- res2$se
#
crosshairs(x1, y1, se.x1, se.y1,
           main_lab = 'Effects of Empirical Bayes Estimation',
           x_lab = 'Parenting and Depression Correlations',
           y_lab = 'Shrunken Estimates',annotate = TRUE,
           whis_on = FALSE)

## End(Not run)
