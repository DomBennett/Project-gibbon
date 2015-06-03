# D.J. Bennett
# 2/06/2015
# Plot trees

# LIBS
library (ape)
library (MoreTreeTools)

# INPUT
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

# SORT NAMES
treenames <- c ('Bininda-Edmonds et al. 2007',
                'Chan et al. 2010',
                'Faurby and Svenning 2015',
                'Perelman et al. 2011*',
                'Springer et al. 2012*',
                'Thinh et al. 2010',
                'Thinh et al. 2010 edited (Nomascus annamensis added)')

# GET MEAN SPLIT AGE
ages <- rep (NA, length (trees))
for (i in 1:length (trees)) {
  ages[i] <- getSize (trees[[i]], 'rtt')
}
mean.age <- mean (ages[ages > 10])
print (mean.age)

# MAKE ULTRAMETRIC
for (i in 1:length (trees)) {
  if(ages[i] < 10) {
    if (!is.binary.tree (trees[[i]])) {
      trees[[i]] <- multi2di (trees[[i]])
    }
    trees[[i]] <- chronos (trees[[i]])
    trees[[i]]$edge.length <- trees[[i]]$edge.length*mean.age
  }
}
class (trees) <- 'multiPhylo'

# SORT BININDA NAMES
trees[[1]]$tip.label <- c ("Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes",
                           "Pongo_pygmaeus", "Hylobates_agilis", "Hylobates_lar", "Hylobates_moloch",
                           "Hylobates_muelleri", "Hylobates_pileatus", "Hylobates_klossii",
                           "Hoolock_hoolock", "Symphalangus_syndactylus", "Nomascus_concolor",
                           "Nomascus_leucogenys", "Nomascus_gabriellae")

# GET REF NAMES
refnames <- trees[[7]]$tip.label
for (i in 1:length (trees)) {
  to.drop <- NULL
  for (j in 1:length (trees[[i]]$tip.label)) {
    bool <- grepl (trees[[i]]$tip.label[j], refnames)
    if (any (bool)) {
      trees[[i]]$tip.label[j] <- refnames[bool]
    } else {
      to.drop <- c (to.drop, trees[[i]]$tip.label[j])
    }
  }
  trees[[i]] <- drop.tip (trees[[i]], tip=to.drop)
}

# PLOT
pdf (file.path ('3_plot', 'tree.pdf'))
for (i in 1:length (trees)) {
  tree <- drop.tip (trees[[i]],
                    tip=trees[[i]]$tip.label[!trees[[i]]$tip.label %in% refnames])
  write.tree (tree, file=file.path ('3_plot', treefiles[i]))
  plot (trees[[i]], main=treenames[i])
  axisPhylo ()
  mtext ('MYA', side=1, line=2)
}
dev.off()
