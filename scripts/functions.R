#!/usr/bin/env Rscript

source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/clean_dna.R")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/hapCollapse.R")
# save session info to disk 
#si <- sessionInfo()
#sink("RsessionInfo.txt")
#print(si)
#sink()
#citation("rentrez")

# local minima 2
localMinima2 <- function(distobj){#
    den <- density(distobj)
    a <- rep(NA, length(den$y) - 2)
    for (i in 2:(length(den$y) - 1)) a[i - 1] <- den$y[i - 1] > den$y[i] & den$y[i + 1] > den$y[i]
    den$localMinima <- den$x[which(a)]
    den$data.name <- deparse(substitute(distobj))
    den$call <- paste("density.default(", den$data.name, ")", sep = "")
    #print(den$localMinima)
    invisible(den)
}#

# ADAPTED LOCALMINIMA FUNCTION FOR SPP DELIM (WITH BOOTSTRAP) with codon sampling
localMinimaBootCodon <- function(mat, block, reps, model, gamma, dip){#
    clu <- list()
    for (i in 1:reps){#
        si <- seq(from=block, to=ncol(mat), by=block)
        iboot <- sample(si, replace=TRUE)
        boot.ind <- as.vector(rbind(iboot - 2, iboot - 1, iboot))
        matb <- mat[, boot.ind]
        dis <- dist.dna(matb, model=model, pairwise.deletion=TRUE, gamma=gamma)# 
        #dis <- dist.ml(as.phyDat(matb), model="F81", exclude="pairwise", k=4, shape=0.19)#? phangorn
        thr <- localMinima2(dis)
        lm1 <- thr$localMinima[dip]
        clu[i] <- length(tclust(dis, threshold = lm1))
        print(clu)
        rm(dis,boot.ind,matb,thr,lm1); gc()
        }#
    unlist(clu)
}#

# mode function from https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux)); 
  ux[tab == max(tab)]
}

# function to remove NAs from a vector
rm_na <- function(x) {
    y <- x[!is.na(x)]
    return(y)
}

# make genbank fasta 
gb_format_fasta <- function(df, gene, product){
    organism.name <- ifelse(test=df$taxonRank=="species", yes=paste(df$genus, df$specificEpithet), no=paste(df$genus, "sp."))
    id.by <- str_replace_all(string=iconv(df$identifiedBy, to='ASCII//TRANSLIT'), pattern=" \\| ", replacement="; ")
    specimen.code <- ifelse(test=df$basisOfRecord=="MaterialSample", 
        yes=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, ":", df$otherCatalogNumbers), no=paste0(df$institutionCode, ":", df$collectionCode, ":", df$otherCatalogNumbers)), #
        no=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, ":", df$catalogNumber), no=paste0(df$institutionCode, ":", df$collectionCode, ":", df$catalogNumber)))
    date <- format(as.Date(df$eventDate, format="%Y-%m-%d"), format="%d-%b-%Y")
    country <- paste0(df$country, ": ", iconv(df$stateProvince, to='ASCII//TRANSLIT'), " State, ", iconv(df$waterBody, to='ASCII//TRANSLIT'), " River drainage")
    lat <- ifelse(test=df$decimalLatitude < 0, yes=paste(str_replace_all(df$decimalLatitude, "-", ""), "S"), no=paste(str_replace_all(df$decimalLatitude, "-", ""), "N"))
    lon <- ifelse(test=df$decimalLongitude < 0, yes=paste(str_replace_all(df$decimalLongitude, "-", ""), "W"), no=paste(str_replace_all(df$decimalLongitude, "-", ""), "E"))
    lat.lon <- paste(lat, lon)
    fas <- ifelse(test=df$basisOfRecord=="MaterialSample", #
        yes=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Bio_material=", specimen.code, "]", " ", "[location=mitochondrion] [mgcode=2]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]",  " ", "[Identified_by=", id.by, "]"), #
        no=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Specimen-voucher=", specimen.code, "]", " ", "[location=mitochondrion] [mgcode=2]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]", " ", "[Identified_by=", id.by, "]"))
    fas <- str_replace_all(string=fas, pattern=" \\[Lat_Lon=NA NA\\]| \\[Country=NA\\]| \\[Collection_date=NA\\]", replacement="")
    fas <- paste(fas, df$nucleotides, sep="\n")
        return(fas)
}

# feature table containing the locations of attributes of the sequence.
gb_features <- function(df, gene, product){
    specimen.code <- ifelse(test=df$basisOfRecord=="MaterialSample", 
        yes=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, "-", df$otherCatalogNumbers), no=paste0(df$institutionCode, "-", df$collectionCode, "-", df$otherCatalogNumbers)), #
        no=ifelse(test=is.na(df$collectionCode), yes=paste0(df$institutionCode, "-", df$catalogNumber), no=paste0(df$institutionCode, "-", df$collectionCode, "-", df$catalogNumber)))
    feature.tab <- paste0(paste0(">Feature", " ", df$otherCatalogNumbers, "_", gene),"\n", # 
        "<1", "\t", ">", nchar(df$nucleotides), "\t", "gene", "\n", #
        "\t", "\t", "\t", "gene", "\t", gene, "\n", #
        "<1", "\t", ">", nchar(df$nucleotides), "\t", "CDS", "\t", "\t", "\n", #
        "\t", "\t", "\t", "product", "\t", product, "\n", #
        "\t", "\t", "\t", "codon_start", "\t", "1")#
    return(feature.tab)
}
