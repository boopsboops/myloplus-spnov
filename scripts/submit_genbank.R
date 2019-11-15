#!/usr/bin/env Rscript
# rm(list = ls())

# load libs
library("tidyverse")
library("ape")


## load funcs to generate GB data
# make genbank fasta 
gb_format_fasta <- function(df, gene, product){
    organism.name <- ifelse(test=df$taxonRank=="species", yes=paste(df$genus, df$specificEpithet), no=paste(df$genus, "sp."))
    id.by <- str_replace_all(string=iconv(df$identifiedBy, to='ASCII//TRANSLIT'), pattern=" \\| ", replacement="; ")
    specimen.code <- paste0(df$institutionCode, ":", df$collectionCode, ":", df$catalogNumber)
    date <- format(as.Date(df$eventDate, format="%Y-%m-%d"), format="%d-%b-%Y")
    country <- paste0(df$country, ": ", iconv(df$stateProvince, to='ASCII//TRANSLIT'), ", ", iconv(df$waterBody, to='ASCII//TRANSLIT'))
    lat <- ifelse(test=df$decimalLatitude < 0, yes=paste(str_replace_all(df$decimalLatitude, "-", ""), "S"), no=paste(str_replace_all(df$decimalLatitude, "-", ""), "N"))
    lon <- ifelse(test=df$decimalLongitude < 0, yes=paste(str_replace_all(df$decimalLongitude, "-", ""), "W"), no=paste(str_replace_all(df$decimalLongitude, "-", ""), "E"))
    lat.lon <- paste(lat, lon)
    gcode <- ifelse(test=df$phylum=="Chordata", yes=2, no=5)
    fas <- ifelse(test=df$basisOfRecord=="MaterialSample", #
        yes=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Bio_material=", specimen.code, "]", " ", "[location=mitochondrion] [mgcode=", gcode, "]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]",  " ", "[Identified_by=", id.by, "]"), #
        no=paste0(">", df$otherCatalogNumbers, "_", gene, " ", "[organism=", organism.name, "]", " ", "[Specimen-voucher=", specimen.code, "]", " ", "[location=mitochondrion] [mgcode=", gcode, "]", " ", "[Collection_date=", date, "]", " ", "[Country=", country, "]", " ", "[Lat_Lon=", lat.lon, "]", " ", "[Identified_by=", id.by, "]"))
    fas <- str_replace_all(string=fas, pattern=" \\[Lat_Lon=NA NA\\]| \\[Country=NA\\]| \\[Collection_date=NA\\]", replacement="")
    fas <- paste(fas, df$nucleotides, sep="\n")
        return(fas)
}

# feature table containing the locations of attributes of the sequence.
gb_features <- function(df, gene, product){
    specimen.code <- paste0(df$institutionCode, ":", df$collectionCode, ":", df$catalogNumber)
    feature.tab <- paste0(paste0(">Feature", " ", df$otherCatalogNumbers, "_", gene),"\n", # 
        "<1", "\t", ">", nchar(df$nucleotides), "\t", "gene", "\n", #
        "\t", "\t", "\t", "gene", "\t", gene, "\n", #
        "<1", "\t", ">", nchar(df$nucleotides), "\t", "CDS", "\t", "\t", "\n", #
        "\t", "\t", "\t", "product", "\t", product, "\n", #
        "\t", "\t", "\t", "codon_start", "\t", "1")#
    return(feature.tab)
}

## load spreadsheets and DNAs

dwc.df <- read_csv(file="../data/tissues-master.csv")
pacus.fas <- read.FASTA(file="../data/pacus-COI-new.fasta")

# check that no seqs are missing from the spreadsheet
setdiff(names(pacus.fas), dwc.df$otherCatalogNumbers)

# filter spreadsheet by these names
dwc.df.share <- dwc.df %>% filter(otherCatalogNumbers %in% names(pacus.fas))
dim(dwc.df.share)


## to genbank submit 

# filter
dwc.df.up <- dwc.df.share %>% filter(is.na(associatedSequences))

# add to db
nuc.list <- sapply(as.character(pacus.fas), paste, collapse="")
#names(nuc.list) <- names(pacus.fas)

dwc.df.nuc <- dwc.df.up %>% mutate(nucleotides=nuc.list[match(otherCatalogNumbers, names(nuc.list))])


# check df
#write_csv(dwc.df.nuc, path="../temp/genbank-submit/nuc_df.csv")

# Here we write a GenBank format fasta file containing critical source modifying annotations
gb.fas <- gb_format_fasta(df=dwc.df.nuc, gene="COI")
write(gb.fas , file="../temp/genbank-submit/sequences.fsa", append=FALSE)# write out the fasta file


gb.feat <- gb_features(df=dwc.df.nuc, gene="COI", product="cytochrome oxidase subunit I")
write(gb.feat, file="../temp/genbank-submit/features.tbl", append=FALSE)# write out

# download from ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/
# submission template: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
# go to ~/Software/tbl2asn
# to run tbl2asn in terminal
#./tbl2asn -t template.sbt -i sequences.fsa -f features.tbl -a s -V vb -T
