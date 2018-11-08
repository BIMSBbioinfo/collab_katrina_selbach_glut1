# collab_katrina_selbach_glut1

The purpose of this repository is to provide the scripts that are used to produce the figures related to the analysis of 
disease-causing mutations affecting short linear motif content in transmembrane proteins and also the analysis of gain/loss of 
short linear motif mediated protein interactions based on the peptide array analysis as documented in the paper Meyer et al, Cell, 2018. 

See paper: https://www.sciencedirect.com/science/article/pii/S0092867418310353

## Requirements 

The vignettes are implemented in R and makes use of various CRAN and Bioconductor packages. 
It also depends on the disorder prediction tool IUPred. 

### slimR package 

> install.packages('devtools')
> devtools::install_github('BIMSBbioinfo/slimR')

### CRAN packages 

> install.packages(c('BiocManager', 'rmarkdown', 'knitr', 'data.table', 'ggplot2', 'ggrepel', 'pbapply', 'stringi'))

### Bioconductor packages 

> BiocManager::install(c('Biostrings', 'biomaRt', 'rtracklayer', 'GenomicRanges'))

### IUPred Disorder Predictor

IUPred source code can be dowloaded from here: http://iupred.enzim.hu/Downloads.php .
After unpacking the source code, cd to the src directory. Compile the code with "cc iupred.c -o iupred"

## Vignettes:

### motif_gains_TM_proteins.humsavar.Rmd

This vignette shows how to reproduce Figures 7B and 7C from Meyer et al. 

To render this vignette type:

> Rscript ./vignettes/render.vignette.R ./vignettes/motif_gains_TM_proteins.humsavar.Rmd ./data

### motif_gains_TM_proteins.clinvar.Rmd

This vignette shows how to reproduce Figures S6B and S6C from Meyer et al. 

To render this vignette type:

> Rscript ./vignettes/render.vignette.R ./vignettes/motif_gains_TM_proteins.clinvar.Rmd ./data
