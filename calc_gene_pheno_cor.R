# run gene and phenotype associations, on clusters
# TODO: change to use unified phenotype file and wire it up with
# separate meta info table of phenotypes

library(jsonlite)
source("./src/source_func.R")
source("./src/config.R")

inputs <- commandArgs(trailingOnly = TRUE)
start <- inputs[1]
end <- inputs[2]

# overwrite to rerun for select cohorts
cohorts <- c("LUAD")

library(parallel)
library(doParallel)
library(foreach)
ncores <- 16

manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type== "protein_coding", c("gene", "gene_name_BCM_refined")]
manifest <- unique(manifest) #19701
ensg_fulls <- manifest$gene
rownames(manifest) <- ensg_fulls


rna_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct", F)
pro_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.cct", F)
scnv_data <- read_all_data(cohorts, "_WES_CNV_gene_ratio_log2.cct", F)

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

calc_assoc_wrapper <- function(mol_data, pheno_data, datatype, program, method) {
    if (method == "spearman") {
        type <- "CON"
    } else if (method == "wilcox") {
        type <- "BIN"
    }
    phenotypes <- get_all_ids(pheno_data)
    results <- NULL
    for (pheno in phenotypes) {
        res <- calc_all_pheno_cor(mol_data, pheno_data, id, pheno, type, cohorts)
        res$metap <- calc_metap_data(res)
        res$symbol <- sym
        res$ensembl <- ensg
        res$ensembl_ver <- ensg_ver
        res$datatype <- datatype
        res$phenotype <- convert_id(ifelse(length(program), paste0(program, ": ", pheno), pheno))
        res$method <- method
        res <- add_sorter(res)
        results <- c(results, toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    }
    return(results)
}

calc_cli_wrapper <- function(mol_data, cli_data, datatype) {
    feature_df <- get_all_cli_features(cli_data)
    results <- NULL
    for (i in seq_len(nrow(feature_df))) {
        pheno <- feature_df[i, "feature"]
        type <- feature_df[i, "type"]
        if (pheno %in% c("Time", "Survival_event")) { next }
        res <- calc_all_cli_cor(mol_data, cli_data, id, pheno, type, cohorts)
        res$metap <- calc_metap_data(res)
        res$symbol <- sym
        res$ensembl <- ensg
        res$ensembl_ver <- ensg_ver
        res$datatype <- datatype
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
    }
    return(results)
}

calc_surv_wrapper <- function(mol_data, cli_data, datatype, time_col, event_col, name="survival", has_type_row=F) {
    # has_type_row = T means true TSI file with 2nd row indicating type
    res <- calc_all_surv(mol_data, cli_data, id, time_col, event_col, cohorts, has_type_row)
    res$metap <- calc_metap_data(res)
    res$symbol <- sym
    res$ensembl <- ensg
    res$ensembl_ver <- ensg_ver
    res$datatype <- datatype
    res$phenotype <- convert_id(paste0("clinical: ", name))
    res$method <- "cox"
    res <- add_sorter(res)
    return(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
}

cls <- makeCluster(ncores)
registerDoParallel(cls)

foreach (i = start:end, .combine = "c",
    .packages = c("metap", "jsonlite", "survival", "survminer", "purrr")) %dopar% {
    id <- ensg_fulls[i]
    tmp <- strsplit(id, ".", fixed = T)[[1]]
    ensg <- tmp[[1]]
    ensg_ver <- tmp[[2]]
    sym <- manifest[id, "gene_name_BCM_refined"]
    res <- NULL

    # RNA
    res <- c(res, calc_assoc_wrapper(rna_data, ciber_data, "RNA", "cibersort", "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, xcell_data, "RNA", "xcell", "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, est_data, "RNA", "ESTIMATE", "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, mb_data, "RNA", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, hallmark_data, "RNA", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, progeny_data, "RNA", "PROGENy", "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, ptmsea_data, "RNA", "PTM-SEA", "spearman"))

    res <- c(res, calc_assoc_wrapper(rna_data, cin_data, "RNA", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, bin_mut_data, "RNA", NULL, "wilcox"))
    res <- c(res, calc_assoc_wrapper(rna_data, mut_sig_data, "RNA", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, purity_wgs_data, "RNA", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(rna_data, purity_wes_data, "RNA", NULL, "spearman"))

    res <- c(res, calc_cli_wrapper(rna_data, cli_data, "RNA"))
    res <- c(res, calc_surv_wrapper(rna_data, surv_data, "RNA", "OS_days", "OS_event", "overall survival", F))
    res <- c(res, calc_surv_wrapper(rna_data, surv_data, "RNA", "PFS_days", "PFS_event", "progression free survival", F))

    # protein
    res <- c(res, calc_assoc_wrapper(pro_data, ciber_data, "protein", "cibersort", "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, xcell_data, "protein", "xcell", "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, est_data, "protein", "ESTIMATE", "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, mb_data, "protein", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, hallmark_data, "protein", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, progeny_data, "protein", "PROGENy", "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, ptmsea_data, "protein", "PTM-SEA", "spearman"))

    res <- c(res, calc_assoc_wrapper(pro_data, cin_data, "protein", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, bin_mut_data, "protein", NULL, "wilcox"))
    res <- c(res, calc_assoc_wrapper(pro_data, mut_sig_data, "protein", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, purity_wgs_data, "protein", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(pro_data, purity_wes_data, "protein", NULL, "spearman"))

    res <- c(res, calc_cli_wrapper(pro_data, cli_data, "protein"))
    res <- c(res, calc_surv_wrapper(pro_data, surv_data, "protein", "OS_days", "OS_event", "overall survival", F))
    res <- c(res, calc_surv_wrapper(pro_data, surv_data, "protein", "PFS_days", "PFS_event", "progression free survival", F))

    # SCNV
    res <- c(res, calc_assoc_wrapper(scnv_data, ciber_data, "SCNV", "cibersort", "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, xcell_data, "SCNV", "xcell", "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, est_data, "SCNV", "ESTIMATE",  "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, mb_data, "SCNV", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, hallmark_data, "SCNV", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, progeny_data, "SCNV", "PROGENy", "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, ptmsea_data, "SCNV", "PTM-SEA", "spearman"))

    res <- c(res, calc_assoc_wrapper(scnv_data, cin_data, "SCNV", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, bin_mut_data, "SCNV", NULL, "wilcox"))
    res <- c(res, calc_assoc_wrapper(scnv_data, mut_sig_data, "SCNV", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, purity_wgs_data, "SCNV", NULL, "spearman"))
    res <- c(res, calc_assoc_wrapper(scnv_data, purity_wes_data, "SCNV", NULL, "spearman"))

    res <- c(res, calc_cli_wrapper(scnv_data, cli_data, "SCNV"))
    res <- c(res, calc_surv_wrapper(scnv_data, surv_data, "SCNV", "OS_days", "OS_event", "overall survival", F))
    res <- c(res, calc_surv_wrapper(scnv_data, surv_data, "SCNV", "PFS_days", "PFS_event", "progression free survival", F))

    dir.create(file.path(output_dir, sym), showWarnings = F)
    saveRDS(res, file.path(output_dir, sym, "gene_pheno_cor.rds"))
}

stopCluster(cls)
