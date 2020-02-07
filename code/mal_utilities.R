## Utility functions for mal067project

library(tidyverse)
library(kableExtra)
library(limma)


#'
#'
runGSEA <- function(mal_dat, form, coef, btm=btm_ind, fdrCut=fdr_cut) {

  ## voom the data with the linear model
  design_vac <- model.matrix(formula(form), mal_dat)
  mal_voom_sm <- voom(mal_dat, design=design_vac)

  ## run camera
  cam_vac <- as_tibble(camera(mal_voom_sm, btm, design=design_vac, contrast=coef),
                       rownames="geneset")

  if(sum(cam_vac$FDR <= fdrCut) >0) {
      cat(knitr::knit_print(cam_vac %>%
        filter(FDR < fdrCut) %>%
        mutate(FDR = signif(FDR, 4),
               PValue = signif(PValue, 4)) %>%
        kable(align=rep("l", ncol(cam_vac))) %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                      full_width=FALSE, position="left") %>%
        scroll_box(height="600px")))
  } else {
     cat("No significant gene-sets found")
   }

  return(cam_vac)
} ## runGSEA
