## format functions

##' This function formats the result of enrichment test (ET)
##' based on Fisher's exact test
##'
##' @title formats the result of enrichment test (ET)
##'
##' @param list a list of the result of enrichment test (ncol = 6)
##' @return data.frame
##' @export
##' @examples
##' data(example.list)
##' res <- formatET(example.list)
##'
##' @author Atsushi Fukushima
formatET <- function(list) {
    if (!is.list(list)) stop("list must be a list.")
    if (length(list) < 6) stop("The number of list must be >6.")
    df <- data.frame(matrix(unlist(list), ncol = 6, byrow = TRUE))
    names(df) <- c("TFs", "m", "n", "N", "k", "p.value")
    df$TFs <- as.character(df$TFs)
    df$m <- as.numeric(as.character(df$m))
    df$n <- as.numeric(as.character(df$n))
    df$N <- as.numeric(as.character(df$N))
    df$k <- as.numeric(as.character(df$k))
    df$p.value <- as.numeric(as.character(df$p.value))
    df.new <- df[order(df$p.value), ]
    return(df.new)
}

##' This function formats the result of Random Control (RC)
##' with random simulation based on Fisher's exact test
##'
##' @title formats the result of Random Control (RC)
##'
##' @param df a data frame of ET including E-values, FDR-BH, and Q-values
##' @param list a list of the result of RC (ncol = 2)
##' @param nRep the number of random selections (negative control)
##' @return data.frame
##' @export
##' @examples
##' data(example.df)
##' data(example.list)
##' nRep <- 100
##' res <- formatRC(example.df, example.list, nRep)
##'
##' @author Atsushi Fukushima
formatRC <- function(df, list, nRep) {
    if (!is.data.frame(df))
        stop("df must be a data.frame.")
    if (!is.list(list))
        stop("list must be a list.")
    if (length(list) < 2)
        stop("The number of list must be >2.")
    rand.df <-
        data.frame(matrix(unlist(list), ncol = 2, byrow = TRUE))
    names(rand.df) <- c("TFs", "sigCount")
    rand.df$RC <-
        as.numeric(as.character(rand.df$sigCount)) / nRep * 100
    rand.df <- rand.df[,-2]
    df.new <- merge(df, rand.df)
    df.new <- df.new[order(df.new$p.value), ]
    return(df.new)
}

##' This function extracts information about transcription factor (TF) and
##' target gene (TG) with TFactS Catalogue (v2).
##'
##' @title extracts transcription factor (TF) and target gene (TG) information
##'
##' @param DEGs a character vector of DEGs (differentially expressed genes)
##' @param catalog a data frame of TFactS catalog (ver. 2)
##' @param TF.col the name of the column that contains the TF names
##' @param TG.col the name of the column that contains the TG names
##' @return list
##' @export
##' @examples
##' data(DEGs)
##' data(catalog)
##'
##' res <- extractTFTG(DEGs, catalog)
##' head(res$TFs)
##'
##' @author Atsushi Fukushima
extractTFTG <- function(DEGs,
                        catalog,
                        TF.col = "TF..OFFICIAL_TF_CODING_GENE_NAME.",
                        TG.col = "Target.gene..OFFICIAL_GENE_NAME.") {
    if (!is.character(DEGs))
        stop("DEGs must be a character.")
    if (!length(DEGs) > 0)
        stop("There is no DEGs.")
    if (!is.data.frame(catalog))
        stop("catalog must be a data.rame.")
    if (!nrow(catalog) > 0)
        stop("There is no catalogue.")

    TFs <-
        unique(as.character(catalog[, grep(TF.col, colnames(catalog))]))
    overlapped <- catalog[, grep(TG.col, colnames(catalog))] %in% DEGs

    if (sum(overlapped) == 0) {
        stop("There is no overlap between DEGs and TF.TG.")
    }
    else {
        subcatalog <- catalog[overlapped,]
        TFs <- unique(as.character(
            subcatalog[, grep(TF.col, colnames(subcatalog))]))
    }

    all.targets <- unique(
        as.character(catalog[, grep(TG.col, colnames(catalog))])
        )
    res <- list(TFs = TFs, all.targets = all.targets)
    return(res)
}
