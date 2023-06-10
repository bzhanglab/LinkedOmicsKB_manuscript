# This is fast, and can be done without parallelization (but should)

# original results were from Xinpei and did not have statistics
# rerun for all cohorts after update of LUAD data
library(jsonlite)
library(mongolite)
source("./src/source_func.R")
source("./src/config.R")


db <- mongo(collection = "gene_tn", db = "pget", url = db_url)
db$drop()

# add gene symbol for searching
manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type== "protein_coding", c("gene", "gene_name_BCM_refined")]
manifest <- unique(manifest)

rna_t_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct", F)
rna_n_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.cct", F)

pro_t_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.cct", F)
pro_n_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.cct", F)


calc_tn <- function(ensg_full, t_data, n_data, datatype) {
    tmp <- strsplit(ensg_full, ".", fixed = T)[[1]]
    ensg <- tmp[1]
    ensg_ver <- tmp[2]
    sym <- manifest[manifest$gene == ensg_full, "gene_name_BCM_refined"]
    res <- list(
        ensembl = ensg, ensembl_ver = ensg_ver,
        symbol = sym,
        datatype = datatype
    )

    for (cohort in cohorts) {
        t_d <- t_data[[cohort]]
        t_row <- NULL
        if (!is.null(t_d)) {
            t_row <- unlist(t_d[ensg_full, ])
        }
        n_d <- n_data[[cohort]]
        n_row <- NULL
        if (!is.null(n_d)) {
            n_row <- unlist(n_d[ensg_full, ])
        }
        if (!is.null(t_row) && !is.null(n_row) &&
                length(na.omit(t_row)) >= 20 && length(na.omit(n_row)) >= 10 &&
                !all(t_row == 0) && !all(n_row == 0)) {
            # T/N analysis
            wx_res <- wilcox.test(t_row[!is.na(t_row) & !is.infinite(t_row)], n_row[!is.na(n_row) & !is.infinite(n_row)], conf.int = T)
            sign <- unname(ifelse(wx_res$estimate >= 0, 1, -1))
            res[[cohort]] <- list(
                val = unname(wx_res$estimate),
                pval = wx_res$p.value * sign
            )
        }
    }

    res[['metap']] <- calc_metap_data(res)
    # limit significant number for cohort result in function below
    res <- add_sorter(res)
    #print(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    db$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
}

for (ensg_full in manifest$gene) {
    calc_tn(ensg_full, rna_t_data, rna_n_data, "RNA")
    calc_tn(ensg_full, pro_t_data, pro_n_data, "protein")
}

db$index(add = '{"symbol" : 1}')
# mongoDB 3.4 does not have wildcard index
db$index(add = paste0('{"sorter.', 'metap', '" : 1}'))
for (coh in cohorts) {
    db$index(add = paste0('{"sorter.', coh, '" : 1}'))
}
db$disconnect()

