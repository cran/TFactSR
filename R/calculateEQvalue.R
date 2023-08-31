## calculate E-value, Q-value, and FDR-BH

##' This function calculates E-value based on .
##'
##' @title calculates E-value
##'
##' @param df a data frame containng p-values
##' @param TFs a character vector of transcription factor
##' @return data.frame
##' @references Essaghir A et al. Nucleic Acids Res. 2010 Jun;38(11):e120.
##' @export
##' @examples
##' data(DEGs)
##' data(catalog)
##'
##' tftg <- extractTFTG(DEGs, catalog)
##' TFs <- tftg$TFs
##'
##' p.value <- runif(10)/(1:10)
##' df <- data.frame(p.value = p.value)
##' res <- calculateEvalue(df, TFs)
##'
##' @author Atsushi Fukushima
calculateEvalue <- function(df, TFs) {
    if (!is.data.frame(df))
        stop("df must be a data.frame.")
    if (!is.character(TFs))
        stop("TFs must be a character.")
    if (length(TFs) < 1)
        stop("The length of TFs must be >0.")
    e.value <- df$p.value * length(TFs)
    df.new <- cbind(df, e.value)
    return(df.new)
}


##' This function calculates Q-value based on Storey.
##'
##' @title calculates Q-value
##'
##' @param df a data frame containng p-values
##' @param lambda a vector of the lambda values utilized to obtain pi0.lambda
##' @return data.frame
##' @references Storey JD, The Annals of Statistics 31:2013-2035 (2003)
##' @export
##' @examples
##' data(example.df)
##' p.value <- example.df$p.value
##' df <- data.frame(p.value = p.value)
##' res <- calculateQvalue(df)
##'
##' @author Atsushi Fukushima
calculateQvalue <- function(df, lambda = seq(0.05, 0.5, 0.01)) {
    if (!is.data.frame(df))
        stop("catalog must be a data.rame.")
    qobj <-
        qvalue::qvalue(p = df$p.value,
                    lambda = lambda,
                    pfdr = TRUE)
    df.new <- cbind(df, q.value = qobj$qvalues)
    return(df.new)
}


##' This function calculates FDR based on BH.
##'
##' @title calculates FDR by Benjamini and Hochberg method
##'
##' @param df a data frame containng p-values
##' @return data.frame
##' @references Benjamini Y and Hochberg Y, J Roy Stat Soc B 57: 289?300 (1995)
##' @export
##' @examples
##' p.value <- runif(10)/(1:10)
##' df <- data.frame(p.value = p.value)
##' res <- calculateFDRBH(df)
##'
##' @author Atsushi Fukushima
## FDR control (B-H)
calculateFDRBH <- function(df) {
    if (!is.data.frame(df))
        stop("catalog must be a data.rame.")
    FDR.BH <- stats::p.adjust(df$p.value, method = "BH")
    df.new <- cbind(df, FDR.BH)
    return(df.new)
}
