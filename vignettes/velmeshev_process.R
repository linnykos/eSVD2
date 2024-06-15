rm(list=ls())

library(Seurat)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/raw_counts_mat.RData")
metadata <- read.csv("../../../../data/sns_autism/meta.txt", sep = "\t", header = T)
gene_metadata <- read.csv("../../../../data/sns_autism/genes.tsv", sep = "\t", header = F)

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
AnnotationFilter::supportedFilters(edb)
AnnotationDbi::keytypes(edb)
gene_name_vec <- AnnotationDbi::mapIds(edb,
                                       keys = gene_metadata[,1],
                                       keytype = "GENEID",
                                       column = "GENENAME",
                                       multiVals = "first")
summary_name_mat <- cbind(gene_name_vec, rownames(mat))
summary_name_mat2 <- summary_name_mat[-which(is.na(summary_name_mat[,1])),]

mat2 <- mat[-which(is.na(summary_name_mat[,1])),]
rownames(mat2) <- as.character(summary_name_mat2[,1])
mat2 <- mat2[-which(duplicated(rownames(mat2))),]

##########

sns <- Seurat::CreateSeuratObject(counts = mat2)
sns[["percent.mt"]] <- Seurat::PercentageFeatureSet(sns, pattern = "^MT-")

set.seed(10)
Seurat::DefaultAssay(sns) <- "RNA"

sns[["celltype"]] <- metadata$cluster
sns[["sample"]] <- metadata$sample
sns[["individual"]] <- metadata$individual
sns[["region"]] <- metadata$region
sns[["age"]] <- metadata$age
sns[["sex"]] <- metadata$sex
sns[["RNA.Integrity.Number"]] <- metadata$RNA.Integrity.Number
sns[["post.mortem.hours"]] <- metadata$post.mortem.interval..hours.
sns[["diagnosis"]] <- metadata$diagnosis
sns[["Capbatch"]] <- metadata$Capbatch
sns[["Seqbatch"]] <- metadata$Seqbatch

save(sns, date_of_run, session_info,
     file = "../../../../data/sns_autism/sns_formatted2.RData")
