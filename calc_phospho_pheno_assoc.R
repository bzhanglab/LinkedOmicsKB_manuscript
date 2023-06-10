library(jsonlite)
library(stringr)
#library(mongolite)
source("./src/source_func.R")
source("./src/config.R")


inputs <- commandArgs(trailingOnly = TRUE)
start <- inputs[1]
end <- inputs[2]

library(parallel)
library(doParallel)
library(foreach)
ncores <- 8
output_dir <- "/export/home/yuxingl/lokb/phospho_res"

# add gene symbol for searching
manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type== "protein_coding", c("gene", "gene_name_BCM_refined")]
manifest <- unique(manifest)
rownames(manifest) <- manifest$gene


site_data <- read_all_data(cohorts, "_phospho_site_abundance_log2_reference_intensity_normalized_isoform_adjusted_Tumor.cct", F)
ids <- unique(reduce(lapply(site_data, rownames), c))

# there is a bug causing duplicate and wrong protein ID for secondary isoform
# protein marked the same as primary isoform, practically duplicated (protein, site)
id_df <- str_split(ids, "\\|", simplify = T)
unique_idx <- !duplicated(id_df[, 1:3])
ids <- ids[unique_idx] #126568


ciber_data <- read_all_data(cohorts, "_RNAseq_cibersort_Tumor.cct", F)

# xcell
xcell_data <- read_all_data(cohorts, "_RNAseq_xCell_Tumor.cct", F)

# ESTIMATE
est_data <- read_all_data(cohorts, "_RNAseq_ESTIMATE_Tumor.cct", F)

# Chromesome instability
cin_data <- read_all_data(cohorts, "_WES_CNV_index.cct", F)
cin_data <- lapply(cin_data, function(x) {
    sub <- x["CIN_index", ]
    rownames(sub) <- "CIN_Score"
    return(sub)
})



# mutation burden
# now just use name in data as index for these two (CIN_score, TMB)
mb_data <- read_all_data(cohorts, "_TMB.cct", T)
# mb_data <- lapply(mb_data, function(x) {
    # rownames(x) <- "mutation burden"
    # return(x)
# })


# binary
bin_mut_data <- read_all_data(cohorts, "_binary_mutation_phenotype.cbt", F)

# mut sig
mut_sig_data <- read_all_data(cohorts, "_WES_mutation_signature.cct", T)

# purity
purity_wgs_data <- read_all_data(cohorts, "_WGS_purity_ploidy.cct", T)
purity_wgs_data <- lapply(purity_wgs_data, function(x) {
    x <- x["purity", , drop = F]  # remove ploidy
    rownames(x) <- "tumor_purity_wgs"
    return(x)
})
purity_wes_data <- read_all_data(cohorts, "_WES_purity_ploidy.cct", T)
purity_wes_data <- lapply(purity_wes_data, function(x) {
    x <- x["purity", , drop = F]  # remove ploidy
    rownames(x) <- "tumor_purity_wes"
    return(x)
})


hallmark_data <- read_all_data(cohorts, "_RNAseq_ssGSEA_hallmark_activity.cct", F)

# PROGENy
progeny_data <- read_all_data(cohorts, "_RNAseq_PROGENy_Tumor.cct", F)

# PTM-SEA
ptmsea_data <- read_all_data(cohorts, "_PTM_SEA.cct", F)

# Clinical; second row indicates data type
cli_data <- read_all_data(cohorts, "_clinical.tsi", F)

# now separate survival with other clinical. And has two sets of survival
surv_data <- read_all_data(cohorts, "_survival.tsi", F)

# convert to ID used for looking up meta data
# frontend still needs to urlencode due to +
# Idealy id is assigned and in data matrix
convert_id <- function(id) {
    return(gsub("[: /]", "_", id))
}

calc_assoc_wrapper <- function(site_data, pheno_data, program, method) {
    if (method == "spearman") {
        type <- "CON"
    } else if (method == "wilcox") {
        type <- "BIN"
    }
    phenotypes <- get_all_ids(pheno_data)
    results <- NULL
    for (pheno in phenotypes) {
        res <- calc_all_pheno_cor(site_data, pheno_data, id, pheno, type, cohorts)
        res$metap <- calc_metap_data(res)
        res$symbol <- sym
        res$protein <- pro
        res$protein_ver <- pro_ver
        res$site <- site
        res$phenotype <- convert_id(ifelse(length(program), paste0(program, ": ", pheno), pheno))
        res$method <- method
        res <- add_sorter(res)
        results <- c(results, toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
        #db$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    }
    return(results)

}

calc_cli_wrapper <- function(site_data, cli_data) {
    feature_df <- get_all_cli_features(cli_data)
    results <- NULL
    for (i in seq_len(nrow(feature_df))) {
        pheno <- feature_df[i, "feature"]
        type <- feature_df[i, "type"]
        if (pheno %in% c("Time", "Survival_event")) { next }
        res <- calc_all_cli_cor(site_data, cli_data, id, pheno, type, cohorts)
        res$metap <- calc_metap_data(res)
        res$symbol <- sym
        res$protein <- pro
        res$protein_ver <- pro_ver
        res$site <- site
        res$phenotype <- convert_id(paste0("clinical: ", pheno))
        if (type == "CON") {
            res$method <- "spearman"
        } else if (type == "BIN") {
            res$method <- "wilcox"
        } else if (type == "ORD") {
            res$method <- "jt"
        } else if (type == "CAT") {
            res$method <- "anova"
        }
        res <- add_sorter(res)
        results <- c(results, toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
        #db$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    }
    return(results)

}

calc_surv_wrapper <- function(site_data, cli_data, time_col, event_col, name="survival", has_type_row=F) {
    # has_type_row = T means true TSI file with 2nd row indicating type
    res <- calc_all_surv(site_data, cli_data, id, time_col, event_col, cohorts, has_type_row)
    res$metap <- calc_metap_data(res)
    res$symbol <- sym
    res$protein <- pro
    res$protein_ver <- pro_ver
    res$site <- site
    res$phenotype <- convert_id(paste0("clinical: ", name))
    res$method <- "cox"
    res <- add_sorter(res)
    return(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    #db$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
}


cls <- makeCluster(ncores)
registerDoParallel(cls)
foreach (i = start:end, .packages = c("metap", "jsonlite", "purrr")) %dopar% {
    id <- ids[i]
    tmp <- strsplit(id, "|", fixed = T)[[1]]
    ensg_full <- tmp[1]
    sym <- manifest[ensg_full, "gene_name_BCM_refined"]
    pro <- tmp[2]
    site <- tmp[3]
    if (is.null(sym)) { next }


    tmp <- strsplit(ensg_full, ".", fixed = T)[[1]]
    ensg <- tmp[[1]]
    ensg_ver <- tmp[[2]]
    tmp <- strsplit(pro, ".", fixed = T)[[1]]
    pro <- tmp[[1]]
    pro_ver <- tmp[[2]]

    dir.create(file.path(output_dir, id), showWarnings = FALSE)
    json <- NULL
    json <- c(json, calc_assoc_wrapper(site_data, ciber_data, "cibersort", "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, xcell_data, "xcell", "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, est_data, "ESTIMATE",  "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, mb_data, NULL, "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, hallmark_data, NULL, "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, progeny_data, "PROGENy", "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, ptmsea_data, "PTM-SEA", "spearman"))

    json <- c(json, calc_assoc_wrapper(site_data, cin_data, NULL, "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, bin_mut_data, NULL, "wilcox"))
    json <- c(json, calc_assoc_wrapper(site_data, mut_sig_data, NULL, "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, purity_wgs_data, NULL, "spearman"))
    json <- c(json, calc_assoc_wrapper(site_data, purity_wes_data, NULL, "spearman"))

    json <- c(json, calc_cli_wrapper(site_data, cli_data))
    json <- c(json, calc_surv_wrapper(site_data, surv_data, "OS_days", "OS_event", "overall survival", F))
    json <- c(json, calc_surv_wrapper(site_data, surv_data, "PFS_days", "PFS_event", "progression free survival", F))
    saveRDS(json, file.path(output_dir, id, "phospho_pheno.rds"))
}

stopCluster(cls)
