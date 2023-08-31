## calculate enrichment analysis based on Fisher's exact test

##' This function performs enrichment test (ET) based on Fisher's exact test
##'
##' @title performs enrichment analysis
##'
##' @param DEGs a character vector of DEGs (differentially expressed genes)
##' @param catalog a data frame of TFactS catalog (ver. 2)
##' @param TFs a character vector of transcription factor
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
##'
##' res <- calculateEnrichmentTest(DEGs, catalog, TFs)
##'
##' @author Atsushi Fukushima
calculateEnrichmentTest <-
    function(DEGs, catalog, TFs,
            TF.col = "TF..OFFICIAL_TF_CODING_GENE_NAME.") {
        if (!is.character(DEGs))
            stop("DEGs must be a character.")
        if (!is.data.frame(catalog))
            stop("catalog must be a data.rame.")
        if (!is.character(TFs))
            stop("TFs must be a character.")

        res <- apply(data.frame(TFs), 1, function(i) {
            TF.targets <-
                as.character(catalog[which(catalog[[TF.col]] == i), 1])
            ##' Definition:
            ##' m is the number of target genes annotated for the TF
            ##' under consideration
            ##' n is the number of query genes
            ##' N is the number of regulations in the catalog
            ##' k is the number of query genes that are annotated as
            ##' regulated by TF (i.e., the intersection between the
            ##' query and the TF signature)
            m <- length(TF.targets)
            n <- length(DEGs)
            N <- dim(catalog)[1]
            k <- length(intersect(TF.targets, DEGs))
            ##
            table <- matrix(c(k, n - k, m - k, N - m - n + k), nrow = 2)
            testRes <- stats::fisher.test(table, alternative = 'g')
            res.tmp <- list(
                i,
                m = m,
                n = n,
                N = N,
                k = k,
                p.value = as.numeric(testRes$p.value)
            )
            return(res.tmp)
        })

        df <- formatET(res)
        return(df)
    }

##' This function calculates Random Control (RC)
##'
##' @title calculates Random Control (RC)
##'
##' @param df a data frame containng p-values
##' @param DEGs a character vector of DEGs (differentially expressed genes)
##' @param catalog a data frame of TFactS catalog (ver. 2)
##' @param TFs a character vector of transcription factor
##' @param all.targets a character vector of all target genes
##' @param TF.col the name of the column that contains the TF names
##' @param lambda a user-specified threshold of E-value (default: 0.05)
##' @param nRep number of random selections (default: 100)
##' @return data.frame
##' @references Essaghir A et al. Nucleic Acids Res. 2010 Jun;38(11):e120.
##' @export
##' @examples
##' data(example.df)
##' data(catalog)
##' data(DEGs)
##'
##' tftg <- extractTFTG(DEGs, catalog)
##' TFs <- tftg$TFs
##' all.targets <- tftg$all.targets
##'
##' res <- calculateRC(example.df, DEGs, catalog, TFs, all.targets)
##'
##' @author Atsushi Fukushima
calculateRC <-
    function(df,
            DEGs,
            catalog,
            TFs,
            all.targets,
            TF.col = "TF..OFFICIAL_TF_CODING_GENE_NAME.",
            lambda = 0.05,
            nRep = 100) {
        if (!is.character(DEGs))
            stop("DEGs must be a character.")
        if (!is.data.frame(catalog))
            stop("catalog must be a data.rame.")
        if (!is.character(TFs))
            stop("TFs must be a character.")

        randRes <- apply(data.frame(TFs), 1, function(i) {
            TF.targets <- as.character(
                catalog[which(catalog[[TF.col]] == i), 1])
            ## random sampling
            sigCount <- 0
            for (j in seq_len(nRep)) {
                TF.targets <-
                    sample(all.targets, length(TF.targets), replace = TRUE)
                ##' Definition:
                ##' m is the number of target genes annotated for the TF
                ##' under consideration
                ##' n is the number of query genes
                ##' N is the number of regulations in the catalog
                ##' k is the number of query genes that are annotated as
                ##' regulated by TF (i.e., the intersection between the
                ##' query and the TF signature)
                m <- length(TF.targets)
                n <- length(DEGs)
                N <- dim(catalog)[1]
                k <- length(intersect(TF.targets, DEGs))
                table <- matrix(c(k, n - k, m - k, N - m - n + k), nrow = 2)
                testRes <- stats::fisher.test(table, alternative = 'greater')
                e.value <- as.numeric(testRes$p.value) * length(TFs)
                if (e.value <= lambda)
                    sigCount <- sigCount + 1
            }
            res.tmp <- list(i, sigCount = sigCount)
            return(res.tmp)
        })

        ## RC results table
        df.new <- formatRC(df, randRes, nRep)
        return(df.new)
    }


##' This function calculates Random Control (RC)
##'
##' @title calculates Random Control (RC) fastly?
##'
##' @param df a data frame containng p-values
##' @param DEGs a character vector of DEGs (differentially expressed genes)
##' @param catalog a data frame of TFactS catalog (ver. 2)
##' @param TFs a character vector of transcription factor
##' @param all.targets a character vector of all target genes
##' @param TF.col the name of the column that contains the TF names
##' @param lambda a user-specified threshold of E-value (default: 0.05)
##' @param nRep number of random selections (default: 100)
##' @return data.frame
##' @references Essaghir A et al. Nucleic Acids Res. 2010 Jun;38(11):e120.
##' @export
##' @examples
##' data(example.df)
##' data(catalog)
##' data(DEGs)
##'
##' tftg <- extractTFTG(DEGs, catalog)
##' TFs <- tftg$TFs
##' all.targets <- tftg$all.targets
##'
##' res <- FASTcalculateRC(example.df, DEGs, catalog, TFs, all.targets)
##'
##' @author Atsushi Fukushima
FASTcalculateRC <- function(
            df,
            DEGs,
            catalog,
            TFs,
            all.targets,
            TF.col = "TF..OFFICIAL_TF_CODING_GENE_NAME.",
            lambda = 0.05,
            nRep = 100) {
        if (!is.character(DEGs))
            stop("DEGs must be a character.")
        if (!is.data.frame(catalog))
            stop("catalog must be a data.rame.")
        if (!is.character(TFs))
            stop("TFs must be a character.")

        randRes <- apply(data.frame(TFs), 1, function(i) {
            TF.targets <- as.character(
                catalog[which(catalog[[TF.col]] == i), 1])
            nRep.sampling <- function() {
                sampled <- sample(all.targets, length(TF.targets),
                                  replace = TRUE)
                return(list(sampled))
            }
            ## random sampling
            sampled.TF.targets <- replicate(nRep, nRep.sampling())

            res.rc <- lapply(sampled.TF.targets, function(TF.targets) {
                m <- length(TF.targets)
                n <- length(DEGs)
                N <- dim(catalog)[1]
                k <- length(intersect(TF.targets, DEGs))
                table <- matrix(c(k, n - k, m - k, N - m - n + k), nrow = 2)
                testRes <- stats::fisher.test(table, alternative = 'greater')
                e.value <- as.numeric(testRes$p.value) * length(TFs)
                if (e.value <= lambda)
                    return(TRUE)
                else
                    return(FALSE)
            })

            sigCount <- sum(unlist(res.rc))

            res.tmp <- list(i, sigCount = sigCount)
            return(res.tmp)
        })
        ## RC results table
        df.new <- formatRC(df, randRes, nRep)
        return(df.new)
    }
