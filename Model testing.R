
library(ape)

#opening tree
full_phy <- read.tree("Alldata_constrained_3fossil.tree")
full_phy
length(full_phy$tip.label)
is.ultrametric(full_phy)
plot(full_phy, show.tip.label=F)

#it is not ultrametric, see this post
#http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html
library(phytools)
library(phangorn)
## compute the NNLS ultrametric tree
full_phy<-nnls.tree(cophenetic(full_phy),full_phy,rooted=TRUE)
## check
is.ultrametric(full_phy)

#Now we are making the tree dichotomous
full_phy <- multi2di(full_phy)
is.ultrametric(full_phy)
is.rooted(full_phy)


#Loading the data
mydata <- read.csv("morpho_edited.csv", sep= ";", dec=",")
mydata


#Pruning the tree

full_phy <- drop.tip(full_phy, c("Scirpus polystachyus",
                                 "Trichophorum alpinum",
                                 "Trichophorum caespitosum",
                                 "Sumatroscirpus paniculatocorymbosus|CHC|NCBI Smith 10172",
                                 "Sumatroscirpus rupestris|VIE|NCBI Ford 15081B",
                                 "Eriophorum vaginatum"))

tips.to.remove <- setdiff(full_phy$tip.label, mydata[,1])
phy <- drop.tip(full_phy, tips.to.remove)
setdiff(phy$tip.label, mydata[,1]) 

#The tree is ready for the analyses!!!!!!


#Setting the character
states <- mydata[,4] 
names(states) = mydata[,1]

#Setting the character variation
sd <- mydata[,23]
names(sd) = mydata[,1]

#install.packages("motmot")
library(motmot)


# https://cran.r-project.org/web/packages/motmot/vignettes/motmot.vignette.html#acdc-and-early-burst

###### Maximum Likelihood

#Brownian motion
bm.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                           model = "bm",
                           profilePlot = TRUE,
                           print.warnings = TRUE)
bm.ml

#Brownian motion with error
bm.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              meserr = as.matrix(sd),
                              model = "bm",
                              profilePlot = TRUE,
                              print.warnings = TRUE)
bm.ml_sd


#Lambda
lambda.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                   model = "lambda",
                               profilePlot = TRUE,
                               print.warnings = TRUE)
lambda.ml

#Lambda with error
lambda.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                      meserr = as.matrix(sd),
                                      model = "lambda",
                                  profilePlot = TRUE,
                                  print.warnings = TRUE)
lambda.ml_sd

#Delta

delta.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              model = "delta",
                              upperBound = 50,
                              profilePlot = TRUE,
                              print.warnings = TRUE)
delta.ml

#Delta with error

delta.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              meserr = as.matrix(sd),
                              model = "delta",
                              upperBound = 50,
                              profilePlot = TRUE,
                              print.warnings = TRUE)
delta.ml_sd

#Kappa

kappa.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              model = "kappa",
                              profilePlot = TRUE,
                              print.warnings = TRUE)
kappa.ml

#Kappa with error

kappa.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                 meserr = as.matrix(sd),
                                 model = "kappa",
                                 profilePlot = TRUE,
                                 print.warnings = TRUE)
kappa.ml_sd


#ou

ou.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              model = "ou",
                              profilePlot = TRUE,
                           print.warnings = TRUE)
ou.ml

#Kappa with error

ou.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                 meserr = as.matrix(sd),
                                 model = "ou",
                                 profilePlot = TRUE,
                              print.warnings = TRUE)
ou.ml_sd


#ACDC

ACDC.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                           model = "ACDC",
                           upperBound = 5,
                           profilePlot = TRUE,
                           print.warnings = TRUE)
ACDC.ml

#ACDC with error

ACDC.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              meserr = as.matrix(sd),
                              model = "ACDC",
                              upperBound = 5,
                              profilePlot = TRUE,
                              print.warnings = TRUE)
ACDC.ml_sd

#Early Burst
#Force upperBound to be 0

EA.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                             model = "ACDC",
                             upperBound = -1e-06,
                             profilePlot = TRUE,
                           print.warnings = TRUE)
EA.ml

#EA with error

EA.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                meserr = as.matrix(sd),
                                model = "ACDC",
                                upperBound = -1e-06,
                                profilePlot = TRUE,
                              print.warnings = TRUE)
EA.ml_sd

#psi

psi.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                             model = "psi",
                             profilePlot = TRUE,
                            print.warnings = TRUE)
psi.ml

#psi with error

psi.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                meserr = as.matrix(sd),
                                model = "psi",
                                profilePlot = TRUE,
                               print.warnings = TRUE)
psi.ml_sd


#Free (each branch has a different rate of trait evolution)

free.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                            model = "free",
                            profilePlot = TRUE,
                            print.warnings = TRUE)
free.ml

#free with error

free.ml_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                               meserr = as.matrix(sd),
                               model = "free",
                               profilePlot = TRUE,
                               print.warnings = TRUE)
free.ml_sd


###################################################################################################
# For psi analysis, specify the tree as full.phy in order to get a more accurate speciation rates #
###################################################################################################
#psi

#psi

psi.full_phy.ml <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                            model = "psi",
                            full.phy = full_phy,
                            hiddenSpeciation = TRUE,
                            profilePlot = TRUE,
                            print.warnings = TRUE)
psi.full_phy.ml

#psi with error

psi.full_phy.ml_se <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                               meserr = as.matrix(sd),
                               model = "psi",
                               full.phy = full_phy,
                               hiddenSpeciation = TRUE,
                               profilePlot = TRUE,
                               print.warnings = TRUE)

psi.full_phy.ml_se

###############
# Plus lambda #
###############

#Delta + lambda

delta.ml.lambda <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              model = "delta",
                              profilePlot = TRUE,
                              lambdaEst = TRUE,
                              print.warnings = TRUE)
delta.ml.lambda

#Delta + lambda with error

delta.ml.lambda_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                 meserr = as.matrix(sd),
                                 model = "delta",
                                 profilePlot = TRUE,
                                 lambdaEst = TRUE,
                                 print.warnings = TRUE)
delta.ml.lambda_sd

#Kappa + lambda

kappa.ml.lambda <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              model = "kappa",
                              profilePlot = TRUE,
                              lambdaEst = TRUE,
                              print.warnings = TRUE)
kappa.ml.lambda

#Kappa + lambda with error

kappa.ml.lambda_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                 meserr = as.matrix(sd),
                                 model = "kappa",
                                 profilePlot = TRUE,
                                 lambdaEst = TRUE,
                                 print.warnings = TRUE)
kappa.ml.lambda_sd


#ou + lambda

ou.ml.lambda <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                           model = "ou",
                           profilePlot = TRUE,
                           lambdaEst = TRUE,
                           print.warnings = TRUE)
ou.ml.lambda

#ou + lambda with error

ou.ml.lambda_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                              meserr = as.matrix(sd),
                              model = "ou",
                              profilePlot = TRUE,
                              lambdaEst = TRUE,
                              print.warnings = TRUE)
ou.ml.lambda_sd


#ACDC + lambda

ACDC.ml.lambda <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                             model = "ACDC",
                             profilePlot = TRUE,
                             lambdaEst = TRUE,
                             print.warnings = TRUE)

ACDC.ml.lambda

#ACDC + lambda with error

ACDC.ml.lambda_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                meserr = as.matrix(sd),
                                model = "ACDC",
                                profilePlot = TRUE,
                                lambdaEst = TRUE,
                                print.warnings = TRUE)
ACDC.ml.lambda_sd


#psi + lambda

psi.ml.lambda <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                                   model = "psi",
                                   profilePlot = TRUE,
                                   lambdaEst = TRUE,
                                   print.warnings = TRUE)
                            
psi.ml.lambda

#psi+ lambda with error

psi.ml.lambda_sd <- transformPhylo.ML(phy = phy, y = as.matrix(states), 
                               meserr = as.matrix(sd),
                               model = "psi",
                               profilePlot = TRUE,
                               lambdaEst = TRUE,
                               print.warnings = TRUE)
psi.ml.lambda_sd



###### Bayesian

#Lambda
set.seed(12)
lambda.mcmc <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                           model = "lambda", 
                           mcmc.iteration = 50000, burn.in = 0.1, 
                           random.start = FALSE, sample.every = 1)
lambda.mcmc[1:4]
mcmc.plot(lambda.mcmc)

#Lambda with error
set.seed(12)
lambda.mcmc_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                              meserr = as.matrix(sd),
                              model = "lambda",
                              mcmc.iteration = 50000, burn.in = 0.1, 
                              random.start = FALSE, sample.every = 1)
lambda.mcmc_sd[1:4]
mcmc.plot(lambda.mcmc_sd)


#Delta

set.seed(12)
delta.mcmc <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                   model = "delta", 
                                   mcmc.iteration = 50000, burn.in = 0.1, 
                                   random.start = FALSE, sample.every = 1)
delta.mcmc[1:4]
mcmc.plot(delta.mcmc)

#Delta with error

set.seed(12)
delta.mcmc_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                      meserr = as.matrix(sd),
                                      model = "delta",
                                      mcmc.iteration = 50000, burn.in = 0.1, 
                                      random.start = FALSE, sample.every = 1)
delta.mcmc_sd[1:4]
mcmc.plot(delta.mcmc_sd)


#Kappa

set.seed(12)
kappa.mcmc <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                  model = "kappa", 
                                  mcmc.iteration = 50000, burn.in = 0.1, 
                                  random.start = FALSE, sample.every = 1)
kappa.mcmc[1:4]
mcmc.plot(kappa.mcmc)

#Kappa with error
set.seed(12)
kappa.mcmc_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                     meserr = as.matrix(sd),
                                     model = "kappa",
                                     mcmc.iteration = 50000, burn.in = 0.1, 
                                     random.start = FALSE, sample.every = 1)
kappa.mcmc_sd[1:4]
mcmc.plot(kappa.mcmc_sd)


#ou
set.seed(12)
ou.mcmc <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                  model = "ou", 
                                  mcmc.iteration = 150000, burn.in = 0.05, 
                                  random.start = FALSE, sample.every = 1)
ou.mcmc[1:4]
mcmc.plot(ou.mcmc)

#OU with error
set.seed(12)
ou.mcmc_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                     meserr = as.matrix(sd),
                                     model = "ou",
                                     mcmc.iteration = 150000, burn.in = 0.05, 
                                     random.start = FALSE, sample.every = 1)
ou.mcmc_sd[1:4]
mcmc.plot(ou.mcmc_sd)


#ACDC

set.seed(12)
ACDC.mcmc <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                               model = "ACDC", 
                               mcmc.iteration = 50000, burn.in = 0.1, 
                               random.start = FALSE, sample.every = 1)
ACDC.mcmc[1:4]
mcmc.plot(ACDC.mcmc)

#ACDC with error
set.seed(12)
ACDC.mcmc_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                  meserr = as.matrix(sd),
                                  model = "ACDC",
                                  mcmc.iteration = 50000, burn.in = 0.1, 
                                  random.start = FALSE, sample.every = 1)
ACDC.mcmc_sd[1:4]
mcmc.plot(ACDC.mcmc_sd)

#psi

set.seed(12)
psi.mcmc <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                 model = "psi", 
                                 mcmc.iteration = 50000, burn.in = 0.1, 
                                 random.start = FALSE, sample.every = 1)
psi.mcmc[1:4]
mcmc.plot(psi.mcmc)

#psi with error

set.seed(12)
psi.mcmc_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                    meserr = as.matrix(sd),
                                    model = "psi",
                                    mcmc.iteration = 50000, burn.in = 0.1, 
                                    random.start = FALSE, sample.every = 1)
psi.mcmc_sd[1:4]
mcmc.plot(psi.mcmc_sd)

###################################################################################################
# For psi analysis, specify the tree as full.phy in order to get a more accurate speciation rates #
###################################################################################################
#psi

set.seed(12)
psi.mcmc.full_phy <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                model = "psi",
                                full.phy = full_phy,
                                mcmc.iteration = 50000, burn.in = 0.1, 
                                random.start = FALSE, sample.every = 1)
psi.mcmc.full_phy[1:4]
mcmc.plot(psi.mcmc.full_phy)

#psi with error

set.seed(12)
psi.mcmc.full_phy_sd <- transformPhylo.MCMC(phy = phy, y = as.matrix(states), 
                                   meserr = as.matrix(sd),
                                   model = "psi",
                                   full.phy = full_phy,
                                   mcmc.iteration = 50000, burn.in = 0.1, 
                                   random.start = FALSE, sample.every = 1)
psi.mcmc.full_phy_sd[1:4]
mcmc.plot(psi.mcmc.full_phy_sd)

#acceptance rate seems to be near 0 for OU model


save.image("model testing.RData")
