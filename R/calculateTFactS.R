##' This function calculates TFactS
##'
##' @title calculates TFactS
##'
##' @param DEGs a character vector of DEGs (differentially expressed genes)
##' @param catalog a data frame of TFactS catalog (ver. 2)
##' @param TFs a character vector of transcription factor
##' @param all.targets a character vector of all target genes
##' @param Q.value logical. If it is TRUE, Q.value by Storey method.
##' @param lambda1 a vector of the lambda values utilized to obtain pi0.lambda
##' @param lambda2 a user-specified threshold of E-value (default: 0.05)
##' @param nRep number of random selections (default: 100)
##' @param TF.col the name of the column that contains the TF names
##' @return data.frame
##' @references Essaghir A et al. Nucleic Acids Res. 2010 Jun;38(11):e120.
##' @export
##' @examples
##' data(DEGs)
##' data(catalog)
##'
##' tftg <- extractTFTG(DEGs, catalog)
##' TFs <- tftg$TFs
##' all.targets <- tftg$all.targets
##'
##' res <- calculateTFactS(DEGs, catalog, TFs, all.targets)
##'
##' @author Atsushi Fukushima
calculateTFactS <-
    function(DEGs,
            catalog,
            TFs,
            all.targets,
            Q.value = FALSE,
            lambda1 = seq(0.05, 0.5, 0.01),
            lambda2 = 0.05,
            nRep = 100,
            TF.col = "TF..OFFICIAL_TF_CODING_GENE_NAME.") {
        if (!is.character(DEGs))
            stop("DEGs must be a character.")
        if (!is.data.frame(catalog))
            stop("catalog must be a data.rame.")
        if (!is.character(all.targets))
            stop("all.targets must be a character vector.")
        if (!is.logical(Q.value))
            stop("Q.value must be a logical object.")
        if (!is.character(TFs))
            stop("TFs must be a character.")
        if (!length(TFs) > 0)
            stop("The number of TFs must be >0.")
        if (!length(grep(TF.col, colnames(catalog))) > 0)
            stop("TF.col must specify one column of the catalogue.")

        ## enrichment test
        df <- calculateEnrichmentTest(DEGs, catalog, TFs, TF.col = TF.col)
        ## calculate E-values
        df <- calculateEvalue(df, TFs)
        ## calculate FDR control by BH
        df <- calculateFDRBH(df)
        ## calculate Q-value by Storey
        if (Q.value) {
            df <- calculateQvalue(df, lambda = lambda1)
        }
        ## calculate Random Control (RC)
        df <- calculateRC(
            df,
            DEGs,
            catalog,
            TFs,
            all.targets,
            TF.col = TF.col,
            lambda = lambda2,
            nRep = nRep
        )
        return(df)
    }
