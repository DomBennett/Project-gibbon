# D.J. Bennett
# 18/05/2015
# Calc lambda for all trees in 1_tree/

# LIBS
library( caper )
library( ape )
#library( coxme ) # what does this do? I found I didn't need it.
library( lme4 )
source(file.path ("0_data", "pglm3.4.r")) # Rob's updated/edited code for pagel's lambda test

# PARAMETERS
dist <- FALSE

# FUNCTIONS
init <- function (trait.data, tree) {
  # return list of trait data and phylomat with matching names
  to.drop <- tree$tip.label[!tree$tip.label %in% rownames(trait.data)]
  tree <- drop.tip (tree, tip=to.drop)
  trait.data <- trait.data[rownames(trait.data) %in% tree$tip.label, ]
  phylomat <- vcv (tree)
  res <- list ('phylomat'=phylomat, 'trait'=trait.data)
  res
}
getLambda <- function (tree, trait.data){
  # return lambda for logged trait data
  rownames (trait.data) <- trait.data$species
  data <- init (trait.data, tree)
  data$trait$logwt <- log (data$trait[ ,2]) # transform/normalise
  res <- pglmEstLambda(logwt ~ 1, data=data$trait, phylomat=data$phylomat,
                       plotit=FALSE)
  res$lambda
}
getLambdas <- function (trees, filename, obs, plotit=TRUE,
                        main.title=NULL) {
  # return vector of lambdas for multiple trees
  ntrees <- length (trees)
  trait.data <- read.csv(file.path ('0_data', filename))
  lambdas <- rep (NA, ntrees)
  for (i in 1:ntrees) {
    lambdas[i] <- getLambda (trees[[i]], trait.data)
  }
  if (plotit) {
    hist (lambdas, main=main.title, xlim=c(0,1))
    abline (v=obs, col='red')
  }
  lambdas
}

# INPUT
cat ('\nReading trees ....')
if (!dist) {
  trees <- list ()
  treefiles <- list.files ('1_tree', pattern='\\.tre')
  for (i in 1:length (treefiles)) {
    tree <- read.tree (file.path ('1_tree', treefiles[[i]]))
    if (class (tree) == 'multiPhylo') {
      # choose one in dist at random
      r <- sample (1:length (tree), 1)
      tree <- tree[[r]]
    }
    trees <- c (trees, list (tree))
  }
  treenames <- sub ('\\.tre', '', treefiles)
} else {
  # for dist rename 'Bunopithecus hoolock'
  trees <- read.nexus (file.path ('1_tree', 'dist.nex'))
  class (trees) <- 'list'
  for (i in 1:length (trees)) {
    j <- which (trees[[i]]$tip.label == 'Bunopithecus_hoolock')
    trees[[i]]$tip.label[j] <- 'Hoolock_hoolock'
  }
  class (trees) <- 'multiPhylo'
  treenames <- 1:length (trees)
}
cat ('\nDone. Read in [', length (trees), '].', sep='')

# CALCULATE LAMBDAS AND D
cat ('\nCalculating lambdas and D and plotting ....')
pdf ('2_results/results.pdf')
# BW
obs.bwlambda <- 0.9999339
bwlambdas <- getLambdas (trees, filename="spp.bodyweight.csv",
                         obs=obs.bwlambda, main.title='Body Weight')
# HRJB
obs.hrjblambda <- 0.9999339
hrjblambdas <- getLambdas (trees, filename="mean.homerange.JBHG.csv",
                           obs=obs.hrjblambda, main.title='Home Range (JB)')
# HRZ
obs.hrzjlambda <- 0.9999339
hrzjlambdas <- getLambdas (trees, filename="mean.homerange.ZJHG.csv",
                           obs=obs.hrzjlambda, main.title='Home Range (ZJ)')
# GS
obs.gslambda <- 0.9731351
gslambdas <- getLambdas (trees, filename="mean.groupsize.csv",
                          obs=obs.gslambda, main.title='Group size')
# D
ds <- rep (NA, length (trees))
MSdata <- read.csv (file.path ("0_data", "matsysbinary.csv"))
for (i in 1:length (trees)) {
  t <- trees[[i]]
  t$node.label <- NULL
  cdat <- comparative.data (t, MSdata, names.col="species")
  good.to.go <- length (unique (cdat$data$category)) > 1 &
    all (t$edge.length > 0)
  if (good.to.go) {
    res <- phylo.d (cdat, binvar=category)
    ds[i] <- res$DEstimate
  }
}
hist (ds, main='D estimates', xlim=c(-2,2))
abline (v=-1, col='red')
# plot all trees
if (!dist) {
  for (i in 1:length (trees)) {
    plot (trees[[i]], main=treenames[i])
  }
}
dev.off ()
cat ('\nDone.')

# WRITE
res <- data.frame (tree=treenames, BW=bwlambdas,
                   HRJB=hrjblambdas, HRZJ=hrzjlambdas, GSJ=gslambdas,
                   DEsts=ds)
write.csv (res, file='2_results/table.csv', row.names=FALSE)