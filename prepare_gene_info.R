# Gather basic gene information
# Calculate intra-omics correlation

library(jsonlite)
library(mongolite)
source("./src/source_func.R")
source("./src/config.R")
# separate cor computation now, just basic gene information and expression
output_base <- "/home/yuxingl/PGET/results"

manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
manifest <- manifest[manifest$Type == "protein_coding", c("gene", "transcript", "protein", "gene_name_BCM_refined", "Transcript_order",
    "MANE", "MANE_Plus_Clinical", "Primary_select", "Secondary_select")]
manifest <- unique(manifest) #19701
ensg_fulls <- manifest$gene
genes <- unique(manifest$gene_name_BCM_refined) #19701


meta <- read.delim("ensembl_description.txt", sep = "\t", stringsAsFactors = FALSE)
summary <- read.delim("entrez_summary.txt", sep = "\t", stringsAsFactors = FALSE, quote = "", header = FALSE)
colnames(summary) <- c("entrezgene_id", "summary")
meta <- left_join(meta, summary)

isoform_data <- read_all_data(cohorts, "_RNAseq_isoform_FPKM_log2_Tumor.cct", F)
rna_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct", F)
nrna_data <- read_all_data(cohorts, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.cct", F)
pro_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.cct", F)
npro_data <- read_all_data(cohorts, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.cct", F)
# meth has "   NA", the class is char
meth_data <- read_all_data(cohorts, "_methylation_gene_beta_value_Tumor.cct", F)
scnv_data <- read_all_data(cohorts, "_WES_CNV_gene_ratio_log2.cct", F)

extract_gene <- function(data, gene, cohorts) {
    # extract one row from whole data matrix
    d <- list()
    for (cohort in cohorts) {
        d[[cohort]] <- na.omit(unlist(data[[cohort]][gene, ]))
    }
    return(d)
}

for (i in seq_len(length(genes))) {
#for (i in start:end) {
#for (i in seq_len(50)) {
    sym <- genes[i]
    mani_sub <- manifest %>% filter(gene_name_BCM_refined == sym) %>% arrange(Transcript_order)
    ensg <- mani_sub[1, "gene"]
    ensg_nover <- sub("\\.\\d+$", "", ensg)
    output_dir <- file.path(output_base, sym)
    dir.create(output_dir, showWarnings = FALSE)


    # may have multiple entrez id, and larger one may be obsolete
    meta_ensg <- meta[meta$ensembl_gene_id == ensg_nover, ] %>% arrange(desc(entrezgene_id))
    for (lno in seq_len(nrow(meta_ensg))) {
        summary <- meta_ensg[lno, "summary"]
        entrez <- meta_ensg[lno, "entrezgene_id"]
        description <- sub(" \\[.*\\]$", "", meta_ensg[lno, "description"])
        if (!is.na(summary)) {
            break
        }
    }

    rep_iso <- unlist(mani_sub[mani_sub$Primary_select == "Yes", "transcript"])
    mani_sub$MANE <- mani_sub$MANE == "Yes"
    mani_sub$MANE_Plus_Clinical <- mani_sub$MANE_Plus_Clinical == "Yes"
    mani_sub$Primary_select <- mani_sub$Primary_select == "Yes"
    mani_sub$Secondary_select <- mani_sub$Secondary_select == "Yes"
    mani_sub <- mani_sub[, c("transcript", "protein", "Primary_select", "Secondary_select", "MANE", "MANE_Plus_Clinical")]

    # TODO: split ensembl version and add UniProt
    info <- list(
        symbol = sym, name = description, entrez = entrez, ensembl = ensg,
        summary = summary, isoforms = mani_sub
    )

    cat(toJSON(info, auto_unbox = TRUE), file = file.path(output_dir, "info.json"))

    # isoform plot data
    # data not used now, plot hidden
    isoform_output <- isoform_stat(cohorts, mani_sub$transcript, isoform_data)
    isoform_output$rep <- rep_iso
    cat(toJSON(isoform_output, auto_unbox = T), file = file.path(output_dir, "isoform.json"))

    ## RNA expression
    rna <- extract_gene(rna_data, ensg, cohorts)
    nrna <- extract_gene(nrna_data, ensg, cohorts)

    ## protein experssion
    pro <- extract_gene(pro_data, ensg, cohorts)
    npro <- extract_gene(npro_data, ensg, cohorts)

    ## meth
    meth <- extract_gene(meth_data, ensg, cohorts)

    scnv <- extract_gene(scnv_data, ensg, cohorts)

    ## omics cor
    cors <- list()
    for (cohort in cohorts) {
        data <- list(
            "protein" = pro[[cohort]],
            "RNA" = rna[[cohort]],
            "methylation" = meth[[cohort]],
            "SCNV" = scnv[[cohort]] #,
        )
        cors[[cohort]] <- feature_cor(data)
    }

    cat(toJSON(cors, auto_unbox = T),
        file = file.path(output_dir, "omic_cor.json")
    )

}
