# 03/06/2015
# Dom Bennett
# Construct gibbon tree from Chan et al. 2010

library (MoreTreeTools)

# start with 2 tipped tree
tree <- compute.brlen (stree (2))
# name the most distant species
tree$tip.label <- c ('Hylobates_lar', 'Pongo_abelii')
# give their branches the total length of the tree age
tree$edge.length <- c (19.25, 19.25)

# add tips
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=1, node.age=2.90, tip.name='Hylobates_pileatus')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=2, node.age=4.17, tip.name='Hylobates_agilis')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=3, node.age=2.62, tip.name='Hylobates_muelleri')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=6, node.age=3.45, tip.name='Hylobates_klossii')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=7, node.age=2.77, tip.name='Hylobates_moloch')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=2, node.age=7.52, tip.name='Symphalangus_syndactylus')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=2, node.age=8.67, tip.name='Nomascus_concolor')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=3, node.age=2.37, tip.name='Nomascus_gabriellae')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=17, node.age=1.40, tip.name='Nomascus_siki')
plot (tree)
edgelabels ()
tree <- addTip (tree, edge=19, node.age=0.46, tip.name='Nomascus_leucogenys')

# Cool!
plot (tree)
axisPhylo()
write.tree (tree, file=file.path ('1_tree', 'chan.tre'))