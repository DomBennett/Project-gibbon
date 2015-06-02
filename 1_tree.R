# D.J. Bennett
# 18/05/2015
# Read large numbers of trees and shrink to clade of interest

# PARAMETERS
# identify infile (from 0_data/) and outfile here!
infile <- 'perelman.nex'
outfile <- 'perelman.tre'
ape.names <- c ('Pongo_pygmaeus', 'Pongo_abelii', 'Hoolock_hoolock',
                'Hoolock_leuconedys', 'Hylobates_agilis',
                'Hylobates_albibarbis', 'Hylobates_klossii',
                'Hylobates_lar', 'Hylobates_moloch', 'Hylobates_muelleri',
                'Hylobates_pileatus', 'Nomascus_concolor',
                'Nomascus_gabriellae', 'Nomascus_hainanus',
                'Nomascus_leucogenys', 'Nomascus_nasutus',
                'Nomascus_siki', 'Symphalangus_syndactylus')

# LIBS
library (MoreTreeTools)

# FUNCTIONS
matchNames <- function (tree, names) {
  # Return tree with matching names in names
  for (name in names) {
    bool <- grepl (name, tree$tip.label)
    if (any (bool)) {
      tree$tip.label[bool] <- name
    }
  }
  tree
}

# DIRS
indir <- '0_data'
outdir <- '1_tree'

# READ
cat ('\nReading ....')
if (grepl ('\\.nex', infile)) {
  # if nexus
  trees <- read.nexus (file.path (indir, infile))
} else {
  # if newick
  trees <- read.tree (file.path (indir, infile))
}
if (class (trees) != 'multiPhylo') {
  trees <- list (trees)
  class (trees) <- 'multiPhylo'
}
cat ('\nDone.')

# EXTRACT
cat ('\nExtracting ....')
new.trees <- list ()
for (t in trees) {
  # first make sure names match between ape.names and tip.labels
  t <- matchNames (t, ape.names)
  temp.ape.names <- ape.names[ape.names %in% t$tip.label]
  if (length (temp.ape.names) <= 1) {
    cat ('\nToo few names')
    next
  }
  node <- getParent (t, tips=temp.ape.names)
  t <- extract.clade (phy=t, node=node)
  new.trees <- c (new.trees, list (t))
}
class (new.trees) <- 'multiPhylo'
cat ('\nDone.')

# WRITE
if (length (new.trees) >= 1) {
  write.tree (phy=new.trees, file=file.path (outdir, outfile))
} else {
  cat ('\nNo suitable trees found')
}