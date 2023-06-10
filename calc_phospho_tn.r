library(jsonlite)
library(mongolite)

source("source_func.R")
source("config.R")

db <- mongo(collection = "phospho_tn", db = "pget", url = db_url)
#db$drop()

# add gene symbol for searching
manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type == "protein_coding", c("gene", "gene_name_BCM_refined")]
manifest <- unique(manifest)

t_data <- read_all_data(cohorts, "_phospho_site_abundance_log2_reference_intensity_normalized_isoform_adjusted_Tumor.cct")
n_data <- read_all_data(cohorts, "_phospho_site_abundance_log2_reference_intensity_normalized_isoform_adjusted_Normal.cct")

t_ids <- reduce(lapply(t_data, rownames), c)
n_ids <- reduce(lapply(n_data, rownames), c)
ids <- unique(c(t_ids, n_ids))
print(length(ids))


calc_tn <- function(id) {
    #print(id)
    tmp <- strsplit(id, "|", fixed = T)[[1]]
    ensg <- tmp[1]
    sym <- manifest[manifest$gene == ensg, "gene_name_BCM_refined"]
    if (length(sym) == 0) { return } # not coding gene
    pro <- tmp[2]
    site <- tmp[3]
    tmp <- strsplit(pro, ".", fixed = T)[[1]]
    pro <- tmp[[1]]
    pro_ver <- as.numeric(tmp[[2]])
    res <- list(
        symbol = sym,
        protein = pro,
        protein_ver = pro_ver,
        site = site
    )
    for (cohort in cohorts) {
        t_d <- t_data[[cohort]]
        t_row <- NULL
        if (!is.null(t_d)) {
            t_row <- unlist(t_d[id, ])
        }
        n_d <- n_data[[cohort]]
        n_row <- NULL
        if (!is.null(n_d)) {
            n_row <- unlist(n_d[id, ])
        }
        if (!is.null(t_row) && !is.null(n_row) &&
                length(na.omit(t_row)) >= 20 && length(na.omit(n_row)) >= 10) {
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
    #print(toJSON(res, auto_unbox = TRUE, null = "null", digits=NA))
    db$insert(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
}

for (id in ids) {
    calc_tn(id)
}
db$index(add = '{"protein" : 1}')
#db$index(add = '{"symbol" : "hashed"}')
# mongoDB 3.4 does not have wildcard index
db$index(add = paste0('{"sorter.', 'metap', '" : 1}'))
for (coh in cohorts) {
    db$index(add = paste0('{"sorter.', coh, '" : 1}'))
}
db$disconnect()
