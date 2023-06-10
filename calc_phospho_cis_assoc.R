# calculate phosphosite correlation with its gene's omics

library(jsonlite)
library(mongolite)
library(dplyr)
library(metap)
source("./src/source_func.R")
source("./src/config.R")


library(parallel)
library(doParallel)
library(foreach)
ncores <- 16

site_data <- read_all_data(cohorts, "_phospho_site_abundance_log2_reference_intensity_normalized_isoform_adjusted_Tumor.cct", F)
ids <- unique(reduce(lapply(site_data, rownames), c))



rna_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct", F)
pro_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.cct", F)
meth_data <- read_all_data(cohorts, "_methylation_gene_beta_value_Tumor.cct", F)
scnv_data <- read_all_data(cohorts, "_WES_CNV_gene_ratio_log2.cct", F)

calc_assoc_wrapper <- function(site_data, pheno_data, datatype, method) {
    res <- calc_all_pheno_cor(site_data, pheno_data, id, paste(ensg, ensg_ver, sep="."), "CON", cohorts)
    res$metap <- calc_metap_data(res)
    res$ensembl <- ensg
    res$ensembl_ver <- ensg_ver
    res$protein <- pro
    res$protein_ver <- pro_ver
    res$site <- site
    res$datatype <- datatype
    res$method <- method
    return(toJSON(res, auto_unbox = TRUE, null = "null", digits = I(2)))
    #db$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = I(2)))
}


cls <- makeCluster(ncores)
registerDoParallel(cls)

res <- foreach (i = seq_len(length(ids)), .combine = "c", .packages = c("metap", "jsonlite")) %dopar% {
    id <- ids[i]
    tmp <- strsplit(id, "|", fixed = T)[[1]]
    ensg <- tmp[1]
    pro <- tmp[2]
    site <- tmp[3]
    motif <- tmp[4]
    primary <- tmp[5]
    tmp <- strsplit(ensg, ".", fixed = T)[[1]]
    ensg <- tmp[[1]]
    ensg_ver <- as.numeric(tmp[[2]])
    tmp <- strsplit(pro, ".", fixed = T)[[1]]
    pro <- tmp[[1]]
    pro_ver <- as.numeric(tmp[[2]])

    json <- NULL
    json <- c(json, calc_assoc_wrapper(site_data, rna_data, "RNA", "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, pro_data, "protein", "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, scnv_data, "SCNV", "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, meth_data, "methylation", "spearman"))
    return(json)
}

stopCluster(cls)

db <- mongo(collection = "phospho_cis_cor", db = "pget", url = db_url)
db$drop()
db$insert(res)
db$index(add = '{"protein" : 1}')
db$disconnect()

