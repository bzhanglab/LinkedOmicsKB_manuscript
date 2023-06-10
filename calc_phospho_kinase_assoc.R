# phosphosite correlation with kinase and phosphatase
# now running in cluster

# 04282023 reran not only for LUAD, but all to add symbol and accurate sorter

library(jsonlite)
library(dplyr)
library(metap)
source("./src/config.R")
source("./src/source_func.R")

library(parallel)
library(doParallel)
library(foreach)
ncores <- 16
output_dir <- "/export/home/yuxingl/lokb/phospho_res"

inputs <- commandArgs(trailingOnly = TRUE)
start <- inputs[1]
end <- inputs[2]

manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type== "protein_coding", c("gene", "gene_name_BCM_refined")]
manifest <- unique(manifest) #19701
rownames(manifest) <- manifest$gene

site_data <- read_all_data(cohorts, "_phospho_site_abundance_log2_reference_intensity_normalized_isoform_adjusted_Tumor.cct", F)

# this has duplicate ENSP and site, due to bugs in site allocation for primary and secondary isforms
# Fix in data in the future
ids <- unique(reduce(lapply(site_data, rownames), c))
id_df <- str_split(ids, "\\|", simplify = T)
unique_idx <- !duplicated(id_df[, 1:3])
ids <- ids[unique_idx] #126568



rna_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct", F)
pro_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.cct", F)
meth_data <- read_all_data(cohorts, "_methylation_gene_beta_value_Tumor.cct", F)
scnv_data <- read_all_data(cohorts, "_WES_CNV_gene_ratio_log2.cct", F)

kinase_list <- read.table("/data/PGET_data_freeze/genelists_from_Sara/genelists/genelist_human_protein_kinase.txt", header = T, sep = "\t", check.names = FALSE, quote="")
phosphatase_list <- read.table("/data/PGET_data_freeze/genelists_from_Sara/genelists/genelist_human_phosphatase.txt", header = T, sep = "\t", check.names = FALSE, quote="")

calc_assoc_wrapper_kinase <- function(site_data, pheno_data, target_gene, target_ensg, datatype, method) {
    res <- calc_all_pheno_cor(site_data, pheno_data, id, target_ensg, "CON", cohorts)
    res$metap <- calc_metap_data(res)
    res$ensembl <- ensg
    res$ensembl_ver <- ensg_ver
    res$symbol <- sym
    res$protein <- pro
    res$protein_ver <- pro_ver
    res$site <- site
    res$kinase <- target_gene
    res$datatype <- datatype
    res$method <- method
    res <- add_sorter(res)
    # digits = NA. Presicion is manually controlled. pval: 2
    return(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    #db_kinase$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = I(2)))
}

calc_assoc_wrapper_phosphatase <- function(site_data, pheno_data, target_gene, target_ensg, datatype, method) {
    res <- calc_all_pheno_cor(site_data, pheno_data, id, target_ensg, "CON", cohorts)
    res$metap <- calc_metap_data(res)
    res$ensembl <- ensg
    res$ensembl_ver <- ensg_ver
    res$symbol <- sym
    res$protein <- pro
    res$protein_ver <- pro_ver
    res$site <- site
    res$phosphatase <- target_gene
    res$datatype <- datatype
    res$method <- method
    res <- add_sorter(res)
    return(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    #db_phosphatase$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = I(2)))
}

cls <- makeCluster(ncores)
registerDoParallel(cls)
for (i in start:end) {
    id <- ids[i]
    tmp <- strsplit(id, "|", fixed = T)[[1]]
    ensg_full <- tmp[1]
    sym <- manifest[ensg_full, "gene_name_BCM_refined"]
    if (is.null(sym)) { next }
    pro <- tmp[2]
    site <- tmp[3]
    motif <- tmp[4]
    primary <- tmp[5]
    tmp <- strsplit(ensg_full, ".", fixed = T)[[1]]
    ensg <- tmp[[1]]
    # version could be string as we added _PAR_Y for some
    ensg_ver <- tmp[[2]]
    tmp <- strsplit(pro, ".", fixed = T)[[1]]
    pro <- tmp[[1]]
    pro_ver <- tmp[[2]]

    dir.create(file.path(output_dir, id), showWarnings = FALSE)
    res_ki <- foreach (i = seq_len(nrow(kinase_list)), .combine = "c", .packages = c("metap", "jsonlite")) %dopar% {
        target_gene <- kinase_list[i, "gene_name"]
        target_ensg <- kinase_list[i, "gene"]
        json <- NULL
        json <- c(json, calc_assoc_wrapper_kinase(site_data, rna_data, target_gene, target_ensg, "RNA", "spearman"))
        json <- c(json, calc_assoc_wrapper_kinase(site_data, pro_data, target_gene, target_ensg, "protein", "spearman"))
        json <- c(json, calc_assoc_wrapper_kinase(site_data, scnv_data, target_gene, target_ensg, "SCNV", "spearman"))
        json <- c(json, calc_assoc_wrapper_kinase(site_data, meth_data, target_gene, target_ensg, "methylation", "spearman"))
        return(json)
    }
    saveRDS(res_ki, file.path(output_dir, id, "kinase_cor.rds"))

    res_ph <- foreach (i = seq_len(nrow(phosphatase_list)), .combine = "c", .packages = c("metap", "jsonlite")) %dopar% {
        target_gene <- phosphatase_list[i, "gene_name"]
        target_ensg <- phosphatase_list[i, "gene"]
        json <- NULL
        json <- c(json, calc_assoc_wrapper_phosphatase(site_data, rna_data, target_gene, target_ensg, "RNA", "spearman"))
        json <- c(json, calc_assoc_wrapper_phosphatase(site_data, pro_data, target_gene, target_ensg, "protein", "spearman"))
        json <- c(json, calc_assoc_wrapper_phosphatase(site_data, scnv_data, target_gene, target_ensg, "SCNV", "spearman"))
        json <- c(json, calc_assoc_wrapper_phosphatase(site_data, meth_data, target_gene, target_ensg, "methylation", "spearman"))
        return(json)
    }
    saveRDS(res_ph, file.path(output_dir, id, "phosphatase.rds"))
}

stopCluster(cls)
