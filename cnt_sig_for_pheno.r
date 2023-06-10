library(jsonlite)
library(mongolite)
source("./src/source_func.R")
source("./src/config.R")

# prepare the count of significant association results for plotting

db <- mongo(collection = "phenotype_info", db = "pget", url = db_url)
phenotypes <- unname(unlist(db$find('{}', fields = '{"_id": true}')))
db$disconnect()



get_sig_cnt1 <- function(pheno, datatype, cohort, positive=TRUE) {
    if (positive) {
        return(db$aggregate(paste0('[
            {"$match" : {"phenotype" : "', pheno, '", "datatype": "', datatype, '"}},
            {"$match" : {"sorter.', cohort, '" : {"$gt" : 6}}},
            {"$count" : "total"}
        ]')))
    } else {
        return(db$aggregate(paste0('[
            {"$match" : {"phenotype" : "', pheno, '", "datatype": "', datatype, '"}},
            {"$match" : {"sorter.', cohort, '" : {"$lt" : -6}}},
            {"$count" : "total"}
        ]')))
    }
}

db <- mongo(collection = "gene_pheno_cor", db = "pget", url = db_url)
db_gene_cnt <- mongo(collection = "pheno_sig_gene_cnt", db = "pget", url = db_url)
db_gene_cnt$drop()

for (pheno in phenotypes) {
    for (dt in c("protein", "RNA", "SCNV")) {
        cnts <- list(phenotype = pheno, datatype = dt, threshold = 6, positive = list(), negative = list())
        for (coh in c("metap", cohorts)) {
            res <- get_sig_cnt1(pheno, dt, coh, TRUE)
            if (nrow(res) == 0) {
                cnt <- 0
            } else {
                cnt <- res[1, 1]
            }
            cnts$positive[coh] <- cnt
            res <- get_sig_cnt1(pheno, dt, coh, FALSE)
            if (nrow(res) == 0) {
                cnt <- 0
            } else {
                cnt <- res[1, 1]
            }
            cnts$negative[coh] <- cnt
        }
        db_gene_cnt$insert(toJSON(cnts, auto_unbox = T, na = "null", null = "null"))
    }
}
db$disconnect()

# TN

db <- mongo(collection = "gene_tn", db = "pget", url = db_url)
for (dt in c("protein", "RNA")) {
    cnts <- list(phenotype = "Tumor_normal_compare", datatype = dt, threshold = 6, positive = list(), negative = list())
    for (coh in c("metap", cohorts)) {
        res <- db$aggregate(paste0('[
            {"$match" : {"datatype": "', dt, '"}},
            {"$match" : {"sorter.', coh, '" : {"$gt" : 6}}},
            {"$count" : "total"}
        ]'))
        if (nrow(res) == 0) {
            cnt <- 0
        } else {
            cnt <- res[1, 1]
        }
        cnts$positive[coh] <- cnt
        res <- db$aggregate(paste0('[
            {"$match" : {"datatype": "', dt, '"}},
            {"$match" : {"sorter.', coh, '" : {"$lt" : -6}}},
            {"$count" : "total"}
        ]'))
        if (nrow(res) == 0) {
            cnt <- 0
        } else {
            cnt <- res[1, 1]
        }
        cnts$negative[coh] <- cnt
    }
    db_gene_cnt$insert(toJSON(cnts, auto_unbox = T, na = "null", null = "null"))
}
db$disconnect()

db_gene_cnt$index(add = '{"phenotype" : 1, "datatype": 1}')
db_gene_cnt$disconnect()



get_sig_cnt2 <- function(pheno, cohort, positive=TRUE) {
    if (positive) {
        return(
            db$aggregate(paste0('[
                {"$match" : {"phenotype" : "', pheno, '"}},
                {"$match" : {"sorter.', cohort, '" : {"$gt" : 6}}},
                {"$count" : "total"}
            ]'))
        )
    } else {
        return(
            db$aggregate(paste0('[
                {"$match" : {"phenotype" : "', pheno, '"}},
                {"$match" : {"sorter.', cohort, '" : {"$lt" : -6}}},
                {"$count" : "total"}
            ]'))
        )
    }
}
db <- mongo(collection = "phospho_pheno_cor", db = "pget", url = db_url)
db_phospho_cnt <- mongo(collection = "pheno_sig_phospho_cnt", db = "pget", url = db_url)
db_phospho_cnt$drop()
for (pheno in phenotypes) {
    cnts <- list(phenotype = pheno, threshold = 6, positive = list(), negative = list())
    for (coh in c("metap", cohorts)) {
        res <- get_sig_cnt2(pheno, coh, TRUE)
        if (nrow(res) == 0) {
            cnt <- 0
        } else {
            cnt <- res[1, 1]
        }
        cnts$positive[coh] <- cnt
        res <- get_sig_cnt2(pheno, coh, FALSE)
        if (nrow(res) == 0) {
            cnt <- 0
        } else {
            cnt <- res[1, 1]
        }
        cnts$negative[coh] <- cnt
    }
    db_phospho_cnt$insert(toJSON(cnts, auto_unbox = T, na = "null", null = "null"))
}
db$disconnect()

#TN
db <- mongo(collection = "phospho_tn", db = "pget", url = db_url)
cnts <- list(phenotype = "Tumor_normal_compare", threshold = 6, positive = list(), negative = list())
for (coh in c("metap", cohorts)) {
    res <- db$aggregate(paste0('[
        {"$match" : {"sorter.', coh, '" : {"$gt" : 6}}},
        {"$count" : "total"}
    ]'))
    if (nrow(res) == 0) {
        cnt <- 0
    } else {
        cnt <- res[1, 1]
    }
    cnts$positive[coh] <- cnt
    res <- db$aggregate(paste0('[
        {"$match" : {"sorter.', coh, '" : {"$lt" : -6}}},
        {"$count" : "total"}
    ]'))
    if (nrow(res) == 0) {
        cnt <- 0
    } else {
        cnt <- res[1, 1]
    }
    cnts$negative[coh] <- cnt
}
db_phospho_cnt$insert(toJSON(cnts, auto_unbox = T, na = "null", null = "null"))
db$disconnect()


db_phospho_cnt$index(add = '{"phenotype" : 1}')
db_phospho_cnt$disconnect()

