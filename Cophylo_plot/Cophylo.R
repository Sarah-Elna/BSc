library(phytools)
library(rstudioapi)

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

my_tree = read.tree('concat1_rename_and_manual_trim3.fasta.treefile')
article_tree = read.tree('astral_tree_renamed.tre')

rename_file = read.csv('Rename_00_data_files.csv', sep = ';')

species_names = rename_file$Species_Name
names(species_names) = rename_file$File_Name

article_tree$tip.label = species_names[article_tree$tip.label]

missing_sp = (article_tree$tip.label[which(!article_tree$tip.label %in% my_tree$tip.label)])

article_tree = drop.tip(article_tree, missing_sp)

article_tree$edge.length[is.na(article_tree$edge.length)] <- 0.001

plot(my_tree)

plot(article_tree)

cophylo_t = cophylo(article_tree, my_tree, rotate = TRUE)

plot.cophylo(cophylo_t, cex=0.001)

pdf(paste(wd, "Test1.pdf", sep = ""), height = (15.3), width = (11))
plot(cophylo_t, fsize=0.8)
dev.off()
