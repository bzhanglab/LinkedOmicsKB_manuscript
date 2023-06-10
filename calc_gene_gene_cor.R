library(jsonlite)
library(mongolite)
library(parallel)

source("./src/config.R")
source("./src/source_func.R")


inputs <- commandArgs(trailingOnly = TRUE)
start <- as.numeric(inputs[1])
end <- as.numeric(inputs[2])
mc.cores <- 24


manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type== "protein_coding", c("gene", "gene_name_BCM_refined")]
manifest <- unique(manifest) #19701
ensg_fulls <- manifest$gene
rownames(manifest) <- ensg_fulls


rna_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct", F)
pro_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.cct", F)


calc_assoc_wrapper <- function(mol_data, datatype) {
    run <- function(target_ensg_full) {
        tmp <- strsplit(target_ensg_full, ".", fixed = T)[[1]]
        target_ensg <- tmp[[1]]
        target_ensg_ver <- tmp[[2]]
        target_sym <- manifest[target_ensg_full, "gene_name_BCM_refined"]
        if (target_sym == sym) { return(NA) }

        res <- calc_all_gene_gene_cor(mol_data, ensg_full, target_ensg_full, cohorts)
        res$metap <- calc_metap_data(res)
        res$datatype <- datatype
        res$symbol1 <- sym
        res$ensembl1 <- ensg
        res$ensembl_ver1 <- ensg_ver
        res$symbol2 <- target_sym
        res$ensembl2 <- target_ensg
        res$ensembl_ver2 <- target_ensg_ver
        res$method <- "spearman"
        res <- add_sorter(res)
        return(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA))
    }
    return(na.omit(simplify2array(mclapply(ensg_fulls, run, mc.cores = mc.cores))))
    #return(na.omit(sapply(ensg_fulls, run)))
}

for (i in start:end) {
    ensg_full <- ensg_fulls[i]
    tmp <- strsplit(ensg_full, ".", fixed = T)[[1]]
    ensg <- tmp[[1]]
    ensg_ver <- tmp[[2]]
    sym <- manifest[ensg_full, "gene_name_BCM_refined"]

    dir.create(file.path(output_dir, sym), showWarnings = F)

    res <- calc_assoc_wrapper(rna_data, "RNA")
    saveRDS(res, file.path(output_dir, sym, "gene_cor_RNA.rds"))

    res <- calc_assoc_wrapper(pro_data, "protein")
    saveRDS(res, file.path(output_dir, sym, "gene_cor_protein.rds"))
}
