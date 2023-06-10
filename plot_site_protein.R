library(jsonlite)
library(mongolite)
source("./src/source_func.R")
source("./src/config.R")
library(dplyr)
library(ggrepel)
# plot metap of phosphosite vs protein for each phenotype and label some

output_dir <- "/home/yuxingl/PGET/data/site_vs_protein"

db <- mongo(collection = "phenotype_info", db = "pget", url = db_url)
phenotypes <- unname(unlist(db$find('{}', fields = '{"_id": true}')))
db$disconnect()


# no actual data for some mutations
db <- mongo(collection = "phospho", db = "pget", url = db_url)
id_map <- db$find('{}', fields = '{"ensembl": true, "protein": true, "site": true}')
id_map <- unique(id_map)
id_map <- id_map %>% mutate(site_id = paste(protein, site, sep = "_")) %>%
    select(site_id, ensembl)
db$disconnect()

for (pheno in phenotypes) {
    #tmp, mutations that do not have enough samples for any cohort; i.e. no results
    if (pheno %in% c("CDH1_mutation", "NDUFC2_mutation", "CBFB_mutation", "GLYR1_mutation",
        "PHGR1_mutation", "RPL5_mutation", "SOX17_mutation", "SPANXB1_mutation", "TGIF1_mutation", "TMEM60_mutation")) { next }
    print(pheno)

    db <- mongo(collection = "phospho_pheno_cor", db = "pget", url = db_url)
    site_data <- db$find(paste0('{"phenotype": "', pheno, '"}'),
        fields = '{"sorter.metap": true, "protein": true, "site": true}')
    site_data$metap <- site_data$sorter$metap
    site_data <- site_data %>% mutate(site_id = paste(protein, site, sep = "_")) %>%
        select(site_id, site, logp = metap)
    db$disconnect()


    db <- mongo(collection = "gene_pheno_cor", db = "pget", url = db_url)
    pro_data <- db$find(paste0('{"phenotype": "', pheno, '", "datatype": "protein"}'),
        fields = '{"sorter.metap": true, "ensembl": true, "symbol": true}')
    pro_data$logp <- pro_data$sorter$metap
    db$disconnect()


    pro_data <- pro_data %>% inner_join(id_map)
    pro_data <- pro_data[pro_data$site_id %in% site_data$site_id, ]

    plot_data <- site_data %>% full_join(pro_data, by = c("site_id" = "site_id")) %>%
        mutate(symbol_site = paste(symbol, site, sep = "_"))
    plot_data <- plot_data[!is.na(plot_data$logp.x) & !is.na(plot_data$logp.y), ]

    # x = site, y = protein
    label_data <- plot_data %>% filter(logp.x >= 10 | logp.x <= -10) %>%
        mutate(diff = abs(logp.x) - abs(logp.y)) %>%
        filter(diff >= 5) %>%
        top_n(20, diff)
    p <- ggplot(plot_data, aes(logp.x, logp.y)) + geom_point(size = 1, colour = "grey") +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_abline(slope = 1, intercept = 0) +
        xlab("Signed log of meta p-values of phosphosites") +
        ylab("Signed log of meta p-values of proteins")
    if (nrow(label_data) > 0) {
        p <- p + geom_text_repel(data = label_data, mapping = aes(label = symbol_site, segment.color="grey"),
            size = 2, colour = "red", min.segment.length = 0, max.overlaps = NA, max.time=5, max.iter=50000)
    }
    ggsave(file.path(output_dir, paste0(pheno, ".png")), width = 4, height = 4)
}


db <- mongo(collection = "phospho_tn", db = "pget", url = db_url)
site_data <- db$find(paste0('{}'),
    fields = '{"sorter.metap": true, "protein": true, "site": true}')
site_data$metap <- site_data$sorter$metap
site_data <- site_data %>% mutate(site_id = paste(protein, site, sep = "_")) %>%
    select(site_id, site, logp = metap)
db$disconnect()


db <- mongo(collection = "gene_tn", db = "pget", url = db_url)
pro_data <- db$find('{"datatype": "protein"}',
    fields = '{"sorter.metap": true, "ensembl": true, "symbol": true}')
pro_data$logp <- pro_data$sorter$metap
db$disconnect()


pro_data <- pro_data %>% inner_join(id_map)
pro_data <- pro_data[pro_data$site_id %in% site_data$site_id, ]

plot_data <- site_data %>% full_join(pro_data, by = c("site_id" = "site_id")) %>%
    mutate(symbol_site = paste(symbol, site, sep = "_"))
plot_data <- plot_data[!is.na(plot_data$logp.x) & !is.na(plot_data$logp.y), ]

# x = site, y = protein
label_data <- plot_data %>% filter(logp.x >= 10 | logp.x <= -10) %>%
    mutate(diff = abs(logp.x) - abs(logp.y)) %>%
    filter(diff >= 5) %>%
    top_n(20, diff)
p <- ggplot(plot_data, aes(logp.x, logp.y)) + geom_point(size = 1, colour = "grey") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Signed log of meta p-values of phosphosites") +
    ylab("Signed log of meta p-values of proteins")
if (nrow(label_data) > 0) {
    p <- p + geom_text_repel(data = label_data, mapping = aes(label = symbol_site, segment.color="grey"),
        size = 2, colour = "red", min.segment.length = 0, max.overlaps = NA)
}
ggsave(file.path(output_dir, paste0("Tumor_normal_compare", ".png")), width = 4, height = 4)
