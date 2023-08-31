## ----the original TFactS------------------------------------------------------
library(TFactSR)
data(DEGs)
data(catalog)

tftg <- extractTFTG(DEGs, catalog)
TFs <- tftg$TFs
all.targets <- tftg$all.targets

res <- calculateTFactS(DEGs, catalog, TFs, all.targets)
head(res)

## ----Arabidopsis--------------------------------------------------------------
data(AtCatalog)
data(GenesUp_SH1H)

d <- extractTFTG(GenesUp_SH1H, AtCatalog,
                     TF.col = "TF",
                     TG.col = "target.genes")

res <- calculateTFactS(GenesUp_SH1H, AtCatalog, d$TFs, d$all.targets, TF.col = "TF")
head(res)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

