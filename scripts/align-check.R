#!/usr/bin/env Rscript

library("ape")
library("tidyverse")
library("magrittr")
library("ips")
library("phangorn")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/hapCollapse.R")

# read data
tissues.df <- read_csv("../data/tissues-master.csv")

# make a label col
tissues.df %<>% mutate(label=if_else(taxonRank=="species",paste(catalogNumber,genus,specificEpithet,waterBody),paste(catalogNumber,genus,identificationQualifier,waterBody)))

# read fas
pacus.coi <- read.FASTA("../data/pacus-COI.fasta")

# check names
setdiff(names(pacus.coi), pull(tissues.df,catalogNumber))
setdiff(pull(tissues.df,catalogNumber),names(pacus.coi))

# align
pacus.coi.mat <- as.matrix(mafft(pacus.coi,exec="mafft"))

# trim after checking in geneious
pacus.coi.mat.trimmed <- pacus.coi.mat[,49:669]

# make a tree
pacus.tr <- nj(dist.dna(pacus.coi.mat.trimmed,model="raw",pairwise.deletion=TRUE))

# ladderize tree
pacus.tr <- ladderize(midpoint(pacus.tr))

# make nice labels
pacus.tr$tip.label <- tissues.df$label[match(pacus.tr$tip.label,tissues.df$catalogNumber)]

# remove neg branches
pacus.tr$edge.length[which(pacus.tr$edge.length<0)] <- 0

# write out the tree
pdf(file="../temp/nj.tree.pdf",width=30,height=60)
plot.phylo(pacus.tr,no.margin=TRUE,font=1,label.offset=0.0001)
dev.off()

# write out a nex
write.nexus.data(pacus.coi.mat.trimmed,file="../temp-local-only/pacus.trimmed.nex",format="dna",interleaved=FALSE)


# to same for haplotypes only
pacus.coi.haps <- hapCollapse(data=pacus.coi,cores=8)
pacus.coi.haps.mat <- as.matrix(mafft(pacus.coi.haps,exec="mafft"))
pacus.coi.haps.mat.trimmed <- pacus.coi.haps.mat[,49:669]
write.nexus.data(pacus.coi.haps.mat.trimmed,file="../temp-local-only/pacus.haps.trimmed.nex",format="dna",interleaved=FALSE)
