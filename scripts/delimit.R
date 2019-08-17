#!/usr/bin/env Rscript
# rm(list = ls())
library("ape")
library("phangorn")
library("tidyverse")
library("magrittr")
library("ggtree")
library("splits")
library("bGMYC")
library("spider")
library("RColorBrewer")
library("parallel")
source("functions.R")

# read dnas
dat.all.ali <- as.matrix(as.DNAbin(read.nexus.data("../temp-local-only/pacus.trimmed.nex")))
dat.haps.ali <- as.matrix(as.DNAbin(read.nexus.data("../temp-local-only/pacus.haps.trimmed.nex")))

# load a BEAST tree from the sampled trees
btr.haps <- ladderize(read.nexus(file="../data/pacus.COI.tre"))
# for ggtree
btr.haps.beast <- treeio::read.beast(file="../data/pacus.COI.tre")

# load table and filter
tissues.df <- read_csv("../data/tissues-master.csv") %>% mutate(identifier=if_else(!is.na(associatedSequences),associatedSequences,otherCatalogNumbers))

#tissues.df %<>% filter(catalogNumber %in% labels(dat.all.ali))

# get the lists of collapsed haps and all haps
dat.haps <- clean_dna(dna=dat.haps.ali)
dat.all <- clean_dna(dna=dat.all.ali)

# convert to char
dat.haps.char <- lapply(dat.haps, function(x) paste(x, collapse=""))
dat.all.char <- lapply(dat.all, function(x) paste(x, collapse=""))
dat.daughters.char <- dat.all.char[which(!names(dat.all.char) %in% names(dat.haps.char))]

# detect strings for all daughters and name
seqs.in <- mapply(FUN=function(x) which(str_detect(string=dat.haps.char, pattern=x)==TRUE), dat.daughters.char, SIMPLIFY=TRUE, USE.NAMES=FALSE)
seqs.in <- mapply(function(x) names(dat.haps.char[x]), seqs.in)
names(seqs.in) <- names(dat.daughters.char)

# turn into dataframe
pair.list <- unnest(enframe(seqs.in,name="daughter",value="mother"))

# annotate and collapse
haplotype.list <- pair.list %>% mutate(motherDrainage=pull(tissues.df,waterBody)[match(mother,pull(tissues.df,identifier))]) %>% 
    mutate(daughterDrainage=pull(tissues.df,waterBody)[match(daughter,pull(tissues.df,identifier))]) %>% 
    group_by(mother) %>% 
    summarise(nHaps=n()+1, drainagesHaplotypes=paste(unique(c(motherDrainage,daughterDrainage)),collapse=","), daughters=paste(daughter,collapse=",")) %>% 
    rename(identifier="mother")

# join with tissues
tissues.df.haps <- dplyr::left_join(tissues.df,haplotype.list) %>%
    mutate(statusHaplotype=if_else(identifier %in% names(dat.haps.char),true="mother",false="daughter")) %>%
    mutate(nHaps=if_else(statusHaplotype=="mother" & is.na(nHaps),1,nHaps),drainagesHaplotypes=if_else(statusHaplotype=="mother" & is.na(drainagesHaplotypes),waterBody,drainagesHaplotypes))

# note: a daughter can have multiple mothers if that sequence is shorter
#write_csv(tissues.df.haps,path="../temp/test.csv")
#dd <- read_csv("dups.csv")
#unique(pull(dd,code)[which(duplicated(pull(dd,code)))])


# make a label col
tissues.df.sub <- tissues.df.haps %>% mutate(lab=if_else(taxonRank=="species",paste0(genus," ",specificEpithet," (",nHaps,") ",drainagesHaplotypes),paste0(genus," ",identificationQualifier," (",nHaps,") ",drainagesHaplotypes)))
tissues.df.sub %<>% filter(statusHaplotype=="mother")
ftab <- tissues.df.sub %>% dplyr::select(identifier,lab)



## GMYC results (consensus)
# run gmyc simple
gmyc.res <- gmyc(btr.haps, method="single", quiet=FALSE)
summary(gmyc.res)
gmyc.spec <- spec.list(gmyc.res)
# make df
gmyc.df <- tibble(gr=paste0("gmyc", gmyc.spec$GMYC_spec),labels=as.character(gmyc.spec$sample_name))


## bGMYC results (consensus)
# run bGMYC 
# retains 100 samples after burnin
result.single <- bgmyc.singlephy(phylo=btr.haps, mcmc=11000, burnin=1000, thinning=100, t1=2, t2=length(btr.haps$tip.label), start=c(1,0.5,50))#changing py2 prior is the important parameter - see paper
# check results
#plot(result.single)
#cr.tab <- checkrates(result.single)
#mean(cr.tab[,3])
#plot(density(cr.tab[,3]))
#plot(cr.tab[,3], type="lines")
result.probmat <- spec.probmat(result.single)
splist <- bgmyc.point(probmat=result.probmat, ppcutoff=0.05)# 
names(splist) <- paste0(rep("bgmyc", length(splist)), seq(1:length(splist)))
bgmyc.df <- tibble(gr=rep(names(splist), sapply(splist, length)), labels=unlist(splist))


## localMinima results
# run localMinima
mat <- dist.dna(dat.all.ali, model="raw", pairwise.deletion=TRUE)
lmin <- localMinima2(mat)
lmin$localMinima[1]
clu <- tclust(mat, threshold=lmin$localMinima[1])
nams <- tissues.df$identifier[match(labels(mat), tissues.df$identifier)]
ggg <- lapply(clu, function(x) nams[x])
names(ggg) <- paste0(rep("locmin", length(ggg)), seq(1:length(ggg)))
locmin.df <- tibble(gr=rep(names(ggg), sapply(ggg, length)), labels=unlist(ggg))
# remove the dup haps
locmin.df %<>% filter(labels %in% btr.haps$tip.label)


## mPTP
# optimise ML on haps tree
ml.haps.tr <- optim.pml(pml(tree=btr.haps, data=as.phyDat(dat.haps.ali), model="TrN", k=4, shape=0.5), model="TrN", optInv=FALSE, optNni=FALSE, optGamma=TRUE, optEdge=TRUE)
tr <- midpoint(ml.haps.tr$tree)
is.binary(tr)
sort(tr$edge.length, decreasing=TRUE)
write.tree(tr, file="../temp-local-only/ml.haps.tr.nwk")
# do mPTP in terminal
# `mptp --tree_file ml.haps.tr.nwk --output_file ml.haps.tr.nwk.out --ml --single --minbr 0.0001`
# edit by hand to make the csv file
mptp.df <- read_csv(file="../temp-local-only/ml.haps.tr.nwk.out.csv")


## Process the data 
# merge all the frames
all.df <- plyr::join_all(list(gmyc.df,bgmyc.df,locmin.df,mptp.df), by="labels", type="inner")
df.groups <- data.frame(gmyc=all.df[,1], bgmyc=all.df[,3], locmin=all.df[,4], mptp=all.df[,5], stringsAsFactors=TRUE)

# concatenating into long table
df.concat <- data.frame(indivs=c(as.character(all.df$labels),as.character(all.df$labels),as.character(all.df$labels),as.character(all.df$labels)), delims=c(all.df[,1], all.df[,3], all.df[,4], all.df[,5]), stringsAsFactors=FALSE)
df.concat$delims <- factor(df.concat$delims)

# split
dlist <- by(df.concat, df.concat$delim, function(x) x$indivs)
names(dlist)

# loop 1
ff <- list()
for(i in 1:length(dlist)){#
   ff[[i]] <- which(dlist%in%dlist[i])#
   }#
sff <- lapply(ff, sort)
dff <- sff[!duplicated(sff)]

# loop 2
new.labs <- paste0("sp",rep(1:length(dff)))
dd <- vector(mode="character", length=length(dlist))
for(i in 1:length(dff)){#
dd[dff[[i]]] <- new.labs[i]#
   }#

# names
names(dlist) <- paste(gsub("[0-9]+", "", names(dlist)), dd, sep="-")
aaa.df <- data.frame(gr=rep(names(dlist), sapply(dlist, length)), labels=unlist(dlist), stringsAsFactors=FALSE)
aaa.df %>% group_by(labels)
head(aaa.df)
# join
bbb.df <- plyr::join_all(list(aaa.df[grep("bgmyc", aaa.df$gr), ], aaa.df[grep("^gmyc", aaa.df$gr), ], aaa.df[grep("locmin", aaa.df$gr), ], aaa.df[grep("mptp", aaa.df$gr), ]), by="labels", type="inner")
# check
head(bbb.df)
# assemble
ccc.df <- data.frame(mptp=sapply(strsplit(bbb.df[,5], split="-"), function(x) x[2]), #
    locmin=sapply(strsplit(bbb.df[,4], split="-"), function(x) x[2]), #
    bgmyc=sapply(strsplit(bbb.df[,1], split="-"), function(x) x[2]), #
    gmyc=sapply(strsplit(bbb.df[,3], split="-"), function(x) x[2]))
# create row names
rownames(ccc.df) <- bbb.df$labels
head(ccc.df)


### plot
# make better colours
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
set.seed(9)# 5 is good
cols <- sample(getPalette(n=length(as.character(unique(unlist(ccc.df))))))

# plot the tree
p <- ggtree(btr.haps.beast, ladderize=TRUE, color="grey50", size=0.8) + xlim(0,60)
p <- p %<+% ftab
p <- p + geom_tiplab(aes(label=lab), size=3) + geom_point2(aes(subset=!is.na(posterior) & posterior >= 0.95), color="orange", size=1.25)# color=gr
p <- gheatmap(p=p, data=ccc.df, width=0.5, offset=20) + scale_fill_manual(values=cols) + theme(legend.position="none")
ggsave(plot=p, filename="../temp/delim_all2.pdf", width=14, height=30, bg="transparent", limitsize=FALSE)
rm(p)
