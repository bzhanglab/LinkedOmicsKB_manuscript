library(dplyr)
library(metap)
library(purrr)
#library(SAGx)  # SAGx removed from Bioconductor
library(survival)
library(survminer)

read_data_table <- function(path, transpose = F) {
    if (!file.exists(path)) {
        warning(paste("File not exist:", path))
        return(NULL)
    }
    d <- read.table(path, header = T, row.names = 1, sep = "\t", check.names = FALSE, na.strings=c("NA", "   NA"), quote="")
    if (transpose) d <- t(d)
    return(d)
}
read_all_data <- function(cohorts, suffix, transpose=F) {
    data <- list()
    for (cohort in cohorts) {
        data[[cohort]] <- read_data_table(file.path(data_dir, cohort, paste0(cohort, suffix)), transpose)
    }
    return(data)
}

convert_p_sorter <- function(p) {
    # convert p-values to signed logP directly used for indexing and sorting
    # Has to prepare for NULL, as we want them in the middle
    # Using addField on the fly cannot use Index
    if (is.null(p)) { return(0) }
    sign <- sign(p)
    mlogp <- -log10(abs(p))
    signif(sign * mlogp, 4)
}

# Add sorter and also limit digits here for association results
# Now manage precision here, Use digits=NA for toJSON
add_sorter <- function(d, update_cohorts = NULL) {
    if (!"sorter" %in% names(d)) {
        d$sorter <- list()
    }
    d$sorter$metap <- convert_p_sorter(d$metap)
    if (!is.null(d$metap)) {
        d$metap <- signif(d$metap, 2)
    }
    # if not given, use all cohorts in global environment
    if (is.null(update_cohorts)) {
        update_cohorts <- cohorts
    }
    for (coh in update_cohorts) {
        # sorter always exists
        d$sorter[[coh]] <- convert_p_sorter(d[[coh]][["pval"]])
        if (!is.null(d[[coh]][["pval"]])) {
            d[[coh]][["pval"]] <- signif(d[[coh]][["pval"]], 2)
        }
        if (!is.null(d[[coh]][["val"]])) {
            d[[coh]][["val"]] <- signif(d[[coh]][["val"]], 2)
        }
    }
    d
}

calc_cor <- function(x, y, test = FALSE) {
    x <- na.omit(unlist(x[!is.infinite(x)]))
    y <- na.omit(unlist(y[!is.infinite(y)]))
    overlap <- intersect(names(x), names(y))
    if (length(overlap) < 10 || var(x[overlap], na.rm = TRUE) == 0 || var(y[overlap], na.rm = TRUE) == 0) {
        # require sample num >= 10
        return(NULL)
    }
    if (test) {
        res <- cor.test(x[overlap], y[overlap], use = "na.or.complete", method = "spearman")
        if (is.na(res$estimate)) res <- NA
    } else {
        res <- cor(x[overlap], y[overlap], use = "na.or.complete", method = "spearman")
    }
    return(res)
}

# originally used metap package v1.4. New version supports log input and ouput
# which may help simplify handling extreme cases
calc_metap <- function(pval) {
    pvals0 <- pval[!is.na(pval)]
    if (length(pvals0) > 0) {
        major_sign <- sign(sum(sign(pvals0)))
        if (major_sign == 0)  { major_sign <- 1 }
        pvals <- two2one(abs(pvals0), two = NULL, pvals0 * major_sign < 0)
        if (is.vector(pvals) && length(pvals) > 1) {
            if (median(pvals) > 0.5) {
                # equivalent once all 1-sided, but if most are close to 1 and have many digits
                # sumz output mey lose precision and become 1
                # e.g. pvals = c(0.9999999, 0.8999232, 0.9883498, 0.984223)
                pvals <- 1 - pvals
            }
            pvals[pvals == 1] <- 1 - .Machine$double.neg.eps # or it will be omitted sliently
            pvals[pvals == 0] <- .Machine$double.eps
            p <- as.numeric(sumz(pvals)[["p"]])
            # in extreme cases, it could still be 0 without warning or error
            if (p == 0) { p <- .Machine$double.xmin }
            if (p < 0.5) {
                metap <- 2 * p * major_sign
            } else {
                metap <- 2 * (1 - p) * -1 * major_sign ## add sign to metap
            }
        } else {
            metap <- pvals0[1]
        }
        return(metap)
    }
    return(NULL)
}

# only one row, and necessary data
calc_metap_data <- function(data) {
    data <- data[lapply(data, length) > 0]
    # no data for all cohorts, usually the data like protein is missing
    if (length(data) == 0 ) { return(NULL) }

    pvals <- na.omit(sapply(cohorts, function(coh) {
        return(ifelse(is.null(data[[coh]]$pval), NA, data[[coh]]$pval))
    }))
    return(calc_metap(pvals))
}

# get all ids (phenotypes, genes, sites in all cohorts)
get_all_ids <- function(data) {
    unique(reduce(lapply(data, rownames), c))
}

get_all_cli_features <- function(data) {
    unique(reduce(lapply(data, function(df) {
        return(data.frame(feature = colnames(df), type = unlist(df[1, ])))
    }), rbind.data.frame))
}

calc_all_gene_gene_cor <- function(mol_data, id1, id2, cohorts) {
    calc_all_pheno_cor(mol_data, mol_data, id1, id2, "CON", cohorts)
}

calc_all_pheno_cor <- function(mol_data, pheno_data, id, phenotype, type, cohorts) {
    cors <- list()
    for (cohort in cohorts) {
        coh_data <- pheno_data[[cohort]]
        if (is.null(coh_data) || !phenotype %in% rownames(coh_data)) {
            cors[[cohort]] <- NULL
        } else {
            cors[[cohort]] <- calc_cli_cor(
                unlist(mol_data[[cohort]][id, ]),
                unlist(pheno_data[[cohort]][phenotype, ]),
                type
            )
        }
    }
    return(cors)
}

calc_con_cor <- function(mol_vec, pheno_vec) {
    res <- calc_cor(mol_vec, pheno_vec, T)
    if (!is.null(res) && !is.na(res)) {
        pval <- ifelse(res$p.value == 0, .Machine$double.eps, res$p.value)
        return(list(
            val = unname(res$estimate),
            pval = unname(ifelse(res$estimate == 0, pval, pval * sign(res$estimate)))
        ))
    } else {
        return(NULL)
    }
}

add_count_factor <- function(x) {
    # includes sorting to decide order of factors
    x <- droplevels(as.factor(na.omit(x)))
    for (lvl in levels(x)) {
        count <- sum(x == lvl, na.rm = TRUE)
        levels(x)[levels(x) == lvl] <- paste0(lvl, " (n=", count, ")")
    }
    return(x)
}

calc_all_cli_cor <- function(mol_data, cli_data, id, phenotype, type, cohorts) {
    cors <- list()
    for (cohort in cohorts) {
        vec <- unlist(mol_data[[cohort]][id, ])
        coh_data <- cli_data[[cohort]]
        if (is.null(coh_data) || !phenotype %in% colnames(coh_data)) {
            cors[[cohort]] <- NULL
        } else {
            cli_vec <- coh_data[2:nrow(coh_data), phenotype]
            names(cli_vec) <- rownames(coh_data)[2:nrow(coh_data)]
            cors[[cohort]] <- calc_cli_cor(vec, cli_vec, type)
        }
    }
    return(cors)
}

#' Calculate association based on data type
#'
#' @param vec numeric omic vector with sample names
#' @param cli_vec phenotype data vector with sample names
#' @param type data type of cli_vec
#' @return NULL or a list of statistics and signed p-value
calc_cli_cor <- function(vec, cli_vec, type) {
    # sample on row
    vec <- vec[!is.infinite(vec) & !is.na(vec)]
    if (type == "CON") {
        # TODO age has >= 90
        return(calc_con_cor(vec, sapply(cli_vec, as.numeric)))
    } else if (type == "BIN") {
        overlap <- intersect(names(vec), names(cli_vec))
        val <- add_count_factor(cli_vec[overlap])
        if (length(na.omit(unique(val))) != 2 || length(overlap) < 20 || any(table(val) < 10) || var(vec[overlap], na.rm = TRUE) == 0) {
            # e.g. UCEC is all female
            return(NULL)
        }
        vec <- vec[overlap]
        old_factor <- as.factor(cli_vec[overlap])
        if ("0" %in% old_factor && "1" %in% old_factor) {
            # levels are alphabetically sorted
            # if it is 0 and 1, put 1 first, i.e. mutation
            # NOTE: not really tested
            x <- vec[cli_vec[overlap] == 1]
            y <- vec[cli_vec[overlap] == 0]
        } else {
            x <- vec[val == levels(val)[1]]
            y <- vec[val != levels(val)[1]]
        }
        res <- wilcox.test(x, y, conf.int = T)
        sign <- unname(ifelse(res$estimate >= 0, 1, -1))
        # now use wilcox test
        # d <- data.frame(con = vec[overlap], bin = val, stringsAsFactors = T)
        # res <- t.test(con ~ bin, d)
        pval <- ifelse(res$p.value == 0, .Machine$double.eps, res$p.value)
        return(list(
            val = unname(res$estimate),
            pval = ifelse(res$estimate == 0, pval, pval * sign(res$estimate))
        ))
    } else if (type == "CAT") {
        # NOTE, DO NOT actually have unordered categorical values now
        overlap <- intersect(names(vec), names(cli_vec))
        if (length(overlap) < 20 || var(vec[overlap], na.rm = TRUE) == 0) { return(NULL) }
        val <- add_count_factor(cli_vec[overlap])
        d <- data.frame(con = vec[overlap], cat = val, stringsAsFactors = T)
        res <- anova(lm(con ~ cat, d))
        return(list(val = unname(res$`F value`[[1]]), pval = res$`Pr(>F)`[[1]]))
    } else if (type == "ORD") {
        grp <- cli_vec
        grp[grp == "pTX" | grp == "pNX"] <- NA  # for stage
        grp <- grp[!is.na(grp)]
        overlap <- intersect(names(vec), names(grp))
        if (length(overlap) < 20 || var(vec[overlap], na.rm = TRUE) == 0) { return(NULL) }
        A <- as.matrix(vec[overlap])
        grp <- as.ordered(grp[overlap])
        res <- JT.test(data = A, class = grp, alternative = "two-sided")
        if (is.na(res$medians[1, "rank correlation"])) { return(NULL) }
        pval <- ifelse(res$p.value == 0, .Machine$double.eps, res$p.value)
        return(list(val = unname(res$medians[1, "rank correlation"]),
            pval = ifelse(res$medians[1, "rank correlation"] == 0,
                pval, sign(unname(res$medians[1, "rank correlation"])) * pval))
        )
    }
}

calc_all_surv <- function(mol_data, cli_data, id, time_col, event_col, cohorts, has_type_row) {
    cors <- list()
    for (cohort in cohorts) {
        cors[[cohort]] <- calc_surv(
            unlist(mol_data[[cohort]][id, ]),
            cli_data[[cohort]],
            time_col, event_col, has_type_row
        )
    }
    return(cors)
}

#' survival analysis
#'
#' survival analysis between an omic vector and survival data frame
#' @param vec numeric omic vector with sample names
#' @param cli_data data frame containing survival data
#' @param time_col column name of time data in cli_data
#' @param event_col column name of event data in cli_data
#' @param has_type_row whether cli_data has a second row of data types as in tsi files
#' @return NULL or a list of hazard and signed p-value
calc_surv <- function(vec, cli_data, time_col, event_col, has_type_row=F) {
    start_idx <- ifelse(has_type_row, 2, 1)
    vec <- na.omit(vec[!is.infinite(vec)])
    survdata <- cli_data[start_idx:nrow(cli_data), c(time_col, event_col)]
    if (!ncol(survdata)) { return(NULL) }
    overlap <- intersect(names(vec), rownames(survdata))
    if (length(overlap) < 10 || var(vec[overlap], na.rm = TRUE) == 0) {
        return(NULL)
    }

    vec <- vec[overlap]
    survdata <- survdata[overlap, ]
    event_no <- sum(survdata[[event_col]], na.rm = TRUE)
    group <- factor(vec > median(vec, na.rm = TRUE), c("FALSE", "TRUE"))
    if (any(table(group) < 10) || event_no < 5) {  # if most are zero, group could still be unbalanced
        return(NULL)
    }
    survobj <- list(time=as.numeric(survdata[[time_col]]), status=as.numeric(survdata[[event_col]]), x=vec, group=group)
    # may not converge
    tryCatch({
        res <- summary(coxph(Surv(time, status) ~ x, data = survobj))
        hazard <- signif(res$coefficients[1, 2], 3)
        pval <- signif(res$sctest["pvalue"], 3) * sign(hazard)
    }, error = function(e) {
        hazard <<- NA
        pval <<- NA
    })
    if (is.na(hazard) || is.na(pval)) { return(NULL) }
    if (hazard == 0 && pval == 0) { return(NULL) } # also not converged

    return(list(val = hazard, pval = pval))
}

isoform_stat <- function(cohorts, transcripts, isoform_data) {
    transcripts <- unlist(transcripts)
    all_data <- NULL
    ts_all <- list()
    expr_data <- list()
    for (cohort in cohorts) {
        # isoform <- read.table(file.path("data", paste0(cohort, "_linear_isoform_fpkm_tumor_normal_raw_log2(x+1)_BCM.txt")), header = T, row.names = 1, sep = "\t")
        isoform <- isoform_data[[cohort]]
        isoform <- isoform[transcripts, , drop = FALSE]
        data <- NULL
        for (ts in transcripts) {
            expression <- unlist(isoform[ts, ], use.names = F)
            ts_all[[ts]] <- c(ts_all[[ts]], expression)
            qt <- quantile(expression, c(0.25, 0.5, 0.75))
            data <- c(data, list(list(median = qt[[2]], q1 = qt[[1]], q3 = qt[[3]], cohort = cohort, transcript = ts)))
        }
        all_data <- c(all_data, data)
    }
    hline <- NULL
    for (ts in names(ts_all)) {
        med <- median(ts_all[[ts]])
        hline <- c(hline, list(list(transcript = ts, median = med)))
    }
    list(vline = all_data, hline = hline)
}

feature_cor <- function(data, entrez=NULL) {
    features <- names(data)
    cors <- list()
    for (i in 1:(length(features) - 1)) {
        for (j in (i + 1):length(features)) {
            f1 <- features[i]
            f2 <- features[j]
            res <- calc_cor(data[[f1]], data[[f2]], test = TRUE)
            if (!is.null(res)) {
                d <- list("x" = f1, "y" = f2, "val" = res$estimate, "pval" = res$p.value)
                cors <- c(cors, list(d))
                d <- list("x" = f2, "y" = f1, "val" = res$estimate, "pval" = res$p.value)
                cors <- c(cors, list(d))
            }
        }
    }
    cors
}




#### Copied from SAGx package, which is delisted from BioConductor
# The Lehmann Nonparametrics: Statistical Methods Based on Ranks p. 233 #
# Assumes that groups are coded in increasing numerical order #
# The suggestions posted on R list by Christopher Andrews (SUNY Buffalo, Department of Biostatistics) are gratefully acknowledged

# new idea strsplit(as.character(tet(data = A)), split = "=")[[2]]

JT.test <- function(data, class, labs = NULL, alternative = c("two-sided", "decreasing", "increasing"), ties = FALSE) {
    # data <- Response;class <- Treatment;alternative = "two-sided"
    # decreasing means that the null hypothesis states that the trend is decreasing in higher class
    alternative <- match.arg(alternative)
    calls <- strsplit(as.character(match.call()), split = "=")
    if (!is.ordered(class)) {
        class <- as.ordered(class)
        cat("class was not an ordered factor.  Redefined to be one.\n")
    }
    if (is(data, "exprSet")) stop("Pls update data to ExpressionSet, see Help")
    if (is(data, "ExpressionSet")) {
        pDataX <- pData(data)
        class <- pDataX[, paste(class)]
        data <- exprs(data)
    }
    if (is.factor(data)) {
        data <- as.numeric(data)
        data <- as.matrix(data)
        factor.ind <- TRUE
    } else {
        factor.ind <- FALSE
        data <- as.matrix(data)
    }
    if (!is.null(dim(data))) n.obs <- ncol(data) else n.obs <- length(data)
    if (!(n.obs == length(class))) {
        data <- t(data)
        n.obs <- ncol(data)
    }
    class.tab <- unique(class)
    class.tab <- class.tab[order(class.tab)]
    if (min(dim(data)) == 1) {
        sums <- 0
        upper <- length(class.tab) - 1
        for (i in 1:upper) {
            for (j in seq(i + 1, upper + 1)) {
                x <- t(as.matrix(data[, class == class.tab[i]]))
                y <- t(as.matrix(data[, class == class.tab[j]]))
                n.x <- ncol(x)
                ranked.x <- rank(c(x, y))
                r2 <- sum(ranked.x[1:n.x])
                sums <- sums + r2 - n.x * (n.x + 1) / 2
            }
        }
    } else {
        sums <- 0
        upper <- length(class.tab) - 1
        for (i in 1:upper) {
            for (j in seq(i + 1, upper + 1)) {
                x <- as.matrix(data[, class == class.tab[i]])
                y <- as.matrix(data[, class == class.tab[j]])
                n.x <- ncol(x)
                ranked.x <- t(as.matrix(apply(cbind(x, y), 1, rank)))
                r2 <- rowSums(as.matrix(ranked.x[, 1:n.x]))
                sums <- sums + r2 - n.x * (n.x + 1) / 2
            }
        }
    }
    ni <- table(class)
    EH <- (n.obs^2 - sum(ni^2)) / 4
    if (ties == FALSE) {
        STDH <- sqrt((n.obs^2 * (2 * n.obs + 3) - sum(ni^2 * (2 * ni + 3))) / 72)
    } else {
        dj <- list(apply(data, 1, table))
        term1 <- sapply(dj, function(x) sum(x * (x - 1) * (2 * x + 5)))
        term2 <- sapply(dj, function(x) sum(x * (x - 1) * (x - 2)))
        term3 <- sapply(dj, function(x) sum(x * (x - 1)))
        STDH <- sqrt(
            (n.obs * (n.obs - 1) * (2 * n.obs + 5) - sum(ni * (ni - 1) * (2 * ni + 5)) - term1) / 72 +
                sum(ni * (ni - 1) * (ni - 2) * term2) / (36 * n.obs * (n.obs - 1) * (n.obs - 2)) +
                sum(ni * (ni - 1)) * term3 / (8 * n.obs * (n.obs - 1))
        )
    }

    # continuity correction remains

    ps <- pnorm((EH - sums) / STDH)
    S1 <- (EH - sums) / EH # 09Feb2008 added, see reference Flandre and O'Quigley
    pvalues <- switch(alternative,
        "two-sided" = 2 * pmin(ps, 1 - ps),
        "decreasing" = ps,
        "increasing" = 1 - ps
    )

    # pvalues <- 2*pmin(pnorm((sums-EH)/STDH),1-pnorm((sums-EH)/STDH))

    # utres <- t(apply(data, 1, function(x) c(2*min(pnorm((sum.stat(x,class)-EH)/STDH),1-pnorm((sum.stat(x,class)-EH)/STDH)),tapply(x,class,median),cor(rank(x),rank(class) ) )))

    medians <- t(apply(data, 1, function(x) c(tapply(x, class, median), cor(rank(x), rank(class)))))
    if (factor.ind) medians <- class.tab[medians, -ncol(medians)]
    # utres <- data.frame(pvalues, medians)
    if (is.null(labs)) level <- levels(class) else level <- labs
    alternative <- paste(alternative, paste(level, collapse = switch(alternative,
        two.sided = ", ",
        decreasing = " > ",
        increasing = " < "
    )), sep = ": ")
    data.name <- paste(calls[[2]], "by", calls[[3]])
    ifelse(is.null(labs), colnames(medians)[-ncol(medians)] <- class.tab, colnames(medians)[-ncol(medians)] <- labs)
    colnames(medians)[ncol(medians)] <- "rank correlation"
    rownames(medians) <- rownames(data)
    utres <- list(
        statistic = sums, parameter = NULL, p.value = pvalues, method = "Jonckheere-Terpstra",
        null.value = NULL, alternative = alternative, medians = medians, S1 = S1, data.name = data.name
    )
    class(utres) <- c("JT-test", "htest")
    names(utres$statistic) <- "J"
    return(utres)
}