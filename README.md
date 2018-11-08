# collab_katrina_selbach_glut1

The purpose of this repository is to provide the scripts that are used to produce the figures related to the analysis of 
disease-causing mutations affecting short linear motif content in transmembrane proteins and also the analysis of gain/loss of 
short linear motif mediated protein interactions based on the peptide array analysis as documented in the paper Meyer et al, Cell, 2018. 

See paper: https://www.sciencedirect.com/science/article/pii/S0092867418310353

## Requirements 


## Vignettes:

### motif_gains_TM_proteins.humsavar.Rmd

This vignette shows how to reproduce Figures 7B and 7C from Meyer et al. 

To render this vignette type:

> Rscript ./vignettes/render.vignette.R ./vignettes/motif_gains_TM_proteins.humsavar.Rmd ./data

### motif_gains_TM_proteins.clinvar.Rmd

This vignette shows how to reproduce Figures S6B and S6C from Meyer et al. 

To render this vignette type:

> Rscript ./vignettes/render.vignette.R ./vignettes/motif_gains_TM_proteins.clinvar.Rmd ./data
