##################################
# R source code file for cleaning TCGAOV data used in eclust package
# Git: this is on the TCGAOV repo, master branch
# Created by Sahir,  July 20, 2016
# Updated: January 19, 2017
# Notes:
#
##################################

# source("R/packages.R")
# source("R/functions.R")

# source("/home/bhatnaga/coexpression/rda/packages.R")
# source("/home/bhatnaga/coexpression/rda/functions.R")

# source(paste(Sys.getenv("PBS_O_WORKDIR"),"packages.R", sep="/"))
# source(paste(Sys.getenv("PBS_O_WORKDIR"),"functions.R", sep="/"))

library(pacman)
pacman::p_load(data.table)
pacman::p_load(TCGA2STAT)
pacman::p_load(pryr)
pacman::p_load(magrittr)
pacman::p_load(readxl)

# Get Affymetrix Human Genome U133A expression for ovarian cancer patients
u133a.ov <- getTCGA(disease = "OV", data.type = "mRNA_Array", type = "U133", clinical = TRUE)

# u133a.ov$clinical %>% dim
# u133a.ov$clinical %>% head
# u133a.ov$merged.dat %>% dim
#
# u133a.ov$merged.dat[1:10, 1:20]
# colnames(u133a.ov$dat)
# rownames(u133a.ov$dat)
# dim(u133a.ov$dat)

# we can't use u133a.ov$dat, because the ID's dont map properly to the clinical data
# we have to used the merged data to begin with.

# rows are people, columns are genes
DT_with_pheno <- as.data.table(u133a.ov$merged.dat)

# rows are genes, columns are people
x_mat <- as.matrix(t(DT_with_pheno[,-c("bcr", "status", "OS"), with = F]))
dimnames(x_mat)[[2]] <- DT_with_pheno$bcr

# rows are genes, columns are people
DT <- as.data.table(x_mat, keep.rownames = TRUE)
dim(DT)


# DT <- as.data.table(read.table("RawData/Expression/TCGA_489_UE.txt"), keep.rownames = T)
DT[,1,with = F]
setkey(DT, rn)
colnames(DT) %>% unique() %>% length()


# downloaded from "http://journals.plos.org/plosone/article/asset?unique&id=info:doi/10.1371/journal.pone.0018064.s015"
# rows are genes
# DT.C1 <- as.data.table(read_excel("/home/bhatnaga/coexpression/rda/RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C1 signature genes"))
# DT.C1 <- as.data.table(read_excel(paste(Sys.getenv("PBS_O_WORKDIR"),"RawData/expression/journal.pone.0018064.s015.XLS", sep="/"), sheet = "C1 signature genes"))
DT.C1 <- as.data.table(read_excel("~/git_repositories/eclust/data-raw/journal.pone.0018064.s015.XLS", sheet = "C1 signature genes"))
setkey(DT.C1, geneSymbol)

# DT.C2 <- as.data.table(read_excel("/home/bhatnaga/coexpression/rda/RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C2 signature genes"))
# DT.C2 <- as.data.table(read_excel(paste(Sys.getenv("PBS_O_WORKDIR"),"RawData/expression/journal.pone.0018064.s015.XLS", sep="/"), sheet = "C2 signature genes"))
DT.C2 <- as.data.table(read_excel("~/git_repositories/eclust/data-raw/journal.pone.0018064.s015.XLS", sheet = "C2 signature genes"))
setkey(DT.C2, geneSymbol)

# DT.C4 <- as.data.table(read_excel("/home/bhatnaga/coexpression/rda/RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C4 signature genes "))
# DT.C4 <- as.data.table(read_excel(paste(Sys.getenv("PBS_O_WORKDIR"),"RawData/expression/journal.pone.0018064.s015.XLS", sep="/"), sheet = "C4 signature genes "))
DT.C4 <- as.data.table(read_excel("~/git_repositories/eclust/data-raw/journal.pone.0018064.s015.XLS", sheet = "C4 signature genes "))
setkey(DT.C4, geneSymbol)

# DT.C5 <- as.data.table(read_excel("/home/bhatnaga/coexpression/rda/RawData/Expression/journal.pone.0018064.s015.XLS", sheet = "C5 signature genes"))
# DT.C5 <- as.data.table(read_excel(paste(Sys.getenv("PBS_O_WORKDIR"),"RawData/expression/journal.pone.0018064.s015.XLS", sep="/"), sheet = "C5 signature genes"))
DT.C5 <- as.data.table(read_excel("~/git_repositories/eclust/data-raw/journal.pone.0018064.s015.XLS", sheet = "C5 signature genes"))
setkey(DT.C5, geneSymbol)

DT.C1$geneSymbol %>% unique %>% length()
DT.C2$geneSymbol %>% unique %>% length()
DT.C4$geneSymbol %>% unique %>% length()
DT.C5$geneSymbol %>% unique %>% length()

DT_C1 <- DT[DT.C1]
DT_C1 <- DT_C1[complete.cases(DT_C1)]

DT_C2 <- DT[DT.C2]
DT_C2 <- DT_C2[complete.cases(DT_C2)]

DT_C4 <- DT[DT.C4]
DT_C4 <- DT_C4[complete.cases(DT_C4)]

DT_C5 <- DT[DT.C5]
DT_C5 <- DT_C5[complete.cases(DT_C5)]

indx <- grep('^TCGA', colnames(DT_C1))

for(j in indx){
  set(DT_C1, i=NULL, j=j, value=DT_C1[[j]]*DT_C1[['logFC']])
  set(DT_C2, i=NULL, j=j, value=DT_C2[[j]]*DT_C2[['logFC']])
  set(DT_C4, i=NULL, j=j, value=DT_C4[[j]]*DT_C4[['logFC']])
  set(DT_C5, i=NULL, j=j, value=DT_C5[[j]]*DT_C5[['logFC']])
}

# there are 543 subjects in the dat file
scores <- as.data.table(data.frame(C1 = DT_C1[, colSums(.SD), .SDcols = indx],
                                   C2 = DT_C2[, colSums(.SD), .SDcols = indx],
                                   C4 = DT_C4[, colSums(.SD), .SDcols = indx],
                                   C5 = DT_C5[, colSums(.SD), .SDcols = indx]),
                        keep.rownames = T)
setkey(scores, rn)
scores$rn %>% length()
scores$rn %>% unique %>% length()

scores[, `:=`(C1 = scale(C1), C2 = scale(C2), C4 = scale(C4), C5 = scale(C5))]
scores[, lapply(.SD, mean), .SDcols = 2:5]
scores[, lapply(.SD, var), .SDcols = 2:5]

scores[, subtype := which.max(.SD), by = rn, .SDcols = c("C1","C2","C4","C5")]
scores[, table(subtype)]
scores[subtype %in% 1:2, E := 0]
scores[subtype %ni% 1:2, E := 1]
scores[, table(subtype, E, useNA = "always")]

# IGL@ is named IGL. in the merged data but named IGL@ in the DT.C1 and dat file
# I will rename IGL. in the merged data to IGL@
# grep("IGL", colnames(u133a.ov$merged.dat), value = T)
# grep("IGL", rownames(u133a.ov$dat), value = T)
# grep("IGL", DT.C1$geneSymbol, value = T)

# Gene names starting with HLA have dashes in their names in the dat file
# but have '.' in the merged.dat file. I will rename these to have dahses in the  merged.dat file
# grep("HLA", colnames(u133a.ov$merged.dat), value = T)
# gsub("\\.","-",grep("HLA", colnames(u133a.ov$merged.dat), value = T))

# rows are subjects, columns are genes
# DT_with_pheno <- as.data.table(u133a.ov$merged.dat)
DT_with_pheno[,1:5,with=F]
# DT_with_pheno$bcr %>% unique() %>% length()
#
# DT_with_pheno[, "IGL.", with = F]
# setnames(DT_with_pheno, "IGL.", "IGL@")
# setnames(DT_with_pheno, grep("HLA", colnames(u133a.ov$merged.dat), value = T), gsub("\\.","-",grep("HLA", colnames(u133a.ov$merged.dat), value = T)))

c(DT_C1$rn, DT_C2$rn, DT_C4$rn, DT_C5$rn) %>% length()
c(DT_C1$rn, DT_C2$rn, DT_C4$rn, DT_C5$rn) %>% unique() %>%  length()

DT_temp <- DT_with_pheno[, c("bcr","status","OS",unique(c(DT_C1$rn, DT_C2$rn, DT_C4$rn, DT_C5$rn))), with = F]
setkey(DT_temp, bcr)

# rows are people, columns are genes
# scores[,c("rn","subtype","E"), with = F][DT_temp][,1:10,with = F][,table(E)]

DT_final <- scores[,c("rn","subtype","E"), with = F][DT_temp]

DT_final[,1:10, with = F]
DT_final <- DT_final[!is.na(OS)]

rm(DT_C1,DT_C2, DT_C4, DT_C5, DT_temp, DT_with_pheno, DT,DT.C1,DT.C2, DT.C4, DT.C5, x_mat, scores, j, indx)


DT_final[,table(subtype, E, useNA = "always")]
pryr::object_size(DT_final)

tcgaov <- copy(DT_final)
str(tcgaov)
devtools::use_data(tcgaov, overwrite = TRUE)

tcgaov[,table(subtype, E, useNA = "always")]

# devtools::run_examples()



# DT_pheno <- fread("RawData/Clinical/TCGA Clinical data with TP53 mutation class.csv")
# DT_pheno[BCRPATIENTBARCODE %in% "TCGA-61-1895"]
#
#
# load("RawData/Expression/PT.TCGA_OV.Affy_RMA_DATA.RData")
# DT2 <- as.data.table(DATA, keep.rownames = TRUE)
# DT2[,1,with = F]

# dim(DT)
# head(DT)
# DT
#
# DT.C1$geneSymbol %>% length()
# DT.C1$geneSymbol %>% unique %>% length()
#
#
#
# DT3 <- as.data.table(read.table("RawData/Expression/TCGA_489_UE.txt"), keep.rownames = T)
# setkey(DT3, "")
# head(DT3)
# dim(DT3)
# rownames(DT3)
#
# DT3[,c("rn","TCGA.61.2113.01A.01R"), with = F]
#
# DT[,c("rn","TCGA-61-2113-01A-01R-0668-01"), with = F]
#
# head(DT2)
# DT[, 1, with = F]
#
# colnames(DT)
#
# DT_manifest <- fread("RawData/Expression/file_manifest_Expression-Genes.txt")
#
# all(colnames(DT) %in% DT.manifest$`File Name`)


