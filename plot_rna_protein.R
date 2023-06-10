library(jsonlite)
library(mongolite)
source("./src/source_func.R")
source("./src/config.R")
library(dplyr)
library(ggrepel)
# plot metap of phosphosite vs protein for each phenotype and label some

output_dir <- "/home/yuxingl/PGET/data/rna_vs_protein"

db <- mongo(collection = "phenotype_info", db = "pget", url = db_url)
phenotypes <- unname(unlist(db$find('{}', fields = '{"_id": true}')))
db$disconnect()


for (pheno in phenotypes) {
    print(pheno)

    db <- mongo(collection = "gene_pheno_cor", db = "pget", url = db_url)
    rna_data <- db$find(paste0('{"phenotype": "', pheno, '", "datatype": "RNA"}'),
        fields = '{"sorter.metap": true, "ensembl": true, "symbol": true}')
    rna_data$logp <- rna_data$sorter$metap

    pro_data <- db$find(paste0('{"phenotype": "', pheno, '", "datatype": "protein"}'),
        fields = '{"sorter.metap": true, "ensembl": true, "symbol": true}')
    pro_data$logp <- pro_data$sorter$metap
    db$disconnect()


    plot_data <- pro_data %>% full_join(rna_data, by = c("symbol" = "symbol"))
    plot_data <- plot_data[!is.na(plot_data$logp.x) & !is.na(plot_data$logp.y), ]

    # x = protein, y = rna
    label_data <- plot_data %>% filter(logp.x >= 10 | logp.x <= -10) %>%
        mutate(diff = abs(logp.x) - abs(logp.y)) %>%
        filter(diff >= 5) %>%
        top_n(20, diff)
    p <- ggplot(plot_data, aes(logp.x, logp.y)) + geom_point(size = 1, colour = "grey") +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_abline(slope = 1, intercept = 0) +
        xlab("Signed log of meta p-values of proteins") +
        ylab("Signed log of meta p-values of mRNAs")
    if (nrow(label_data) > 0) {
        p <- p + geom_text_repel(data = label_data, mapping = aes(label = symbol, segment.color="grey"),
            size = 2, colour = "red", min.segment.length = 0, max.overlaps = NA)
    }
    ggsave(file.path(output_dir, paste0(pheno, ".png")), width = 4, height = 4)
}

# Tumor vs Normal
db <- mongo(collection = "gene_tn", db = "pget", url = db_url)
rna_data <- db$find(paste0('{"datatype": "RNA"}'),
    fields = '{"sorter.metap": true, "ensembl": true, "symbol": true}')
rna_data$logp <- rna_data$sorter$metap

pro_data <- db$find('{"datatype": "protein"}',
    fields = '{"sorter.metap": true, "ensembl": true, "symbol": true}')
pro_data$logp <- pro_data$sorter$metap
db$disconnect()


plot_data <- pro_data %>% full_join(rna_data, by = c("symbol" = "symbol"))
plot_data <- plot_data[!is.na(plot_data$logp.x) & !is.na(plot_data$logp.y), ]

# x = protein, y = rna
label_data <- plot_data %>% filter(logp.x >= 10 | logp.x <= -10) %>%
    mutate(diff = abs(logp.x) - abs(logp.y)) %>%
    filter(diff >= 5) %>%
    top_n(20, diff)
p <- ggplot(plot_data, aes(logp.x, logp.y)) + geom_point(size = 1, colour = "grey") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Signed log of meta p-values of proteins") +
    ylab("Signed log of meta p-values of mRNAs")
if (nrow(label_data) > 0) {
    p <- p + geom_text_repel(data = label_data, mapping = aes(label = symbol, segment.color="grey"),
        size = 2, colour = "red", min.segment.length = 0, max.overlaps = NA)
}
ggsave(file.path(output_dir, paste0("Tumor_normal_compare", ".png")), width = 4, height = 4)
