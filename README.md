The purpose of this repository is to provide the scripts that are used to produce the figures related to the analysis of 
disease-causing mutations affecting short linear motif content in transmembrane proteins and also the analysis of gain/loss of 
short linear motif mediated protein interactions based on the peptide array analysis as documented in the paper Meyer et al, Cell, 2018. 

See paper: https://www.sciencedirect.com/science/article/pii/S0092867418310353

# Vignettes:

## motif_gains_transmembrane_proteins

The Rmarkdown scripts in this folder can be run to reproduce figures 7B, 7C, S6B, and S6C from the paper. 

### Requirements

The vignettes are implemented in R and makes use of various CRAN and Bioconductor packages.
It also depends on the disorder prediction tool IUPred.

#### slimR package

> install.packages('devtools')
> devtools::install_github('BIMSBbioinfo/slimR')

#### CRAN packages

> install.packages(c('BiocManager', 'rmarkdown', 'knitr', 'data.table', 'ggplot2', 'ggrepel', 'pbapply', 'stringi'))

#### Bioconductor packages

> BiocManager::install(c('Biostrings', 'biomaRt', 'rtracklayer', 'GenomicRanges'))

#### IUPred Disorder Predictor

IUPred source code can be dowloaded from here: http://iupred.enzim.hu/Downloads.php .
After unpacking the source code, cd to the src directory. Compile the code with "cc iupred.c -o iupred"


### motif_gains_TM_proteins.humsavar.Rmd

This vignette shows how to reproduce Figures 7B and 7C from Meyer et al. 

To render this vignette type:

> Rscript ./vignettes/motif_gains_transmembrane_proteins/render.vignette.R \\
	  ./vignettes/motif_gains_transmembrane_proteins/motif_gains_TM_proteins.humsavar.Rmd \\ 
	  ./data

### motif_gains_TM_proteins.clinvar.Rmd

This vignette shows how to reproduce Figures S6B and S6C from Meyer et al. 

To render this vignette type:

> Rscript ./vignettes/motif_gains_transmembrane_proteins/render.vignette.R \\
	  ./vignettes/motif_gains_transmembrane_proteins/motif_gains_TM_proteins.clinvar.Rmd \\
	  ./data


## peptideArrayAnalysis

The scripts in this folder reproduce the analysis of peptide array pull-down experiment analysis results, in particular Figure S2B and Data S1. 

### Required R packages

The required R packages can be installed via:

> install.packages(c('cowplot', 'data.table', 'DT', 'ggplot2', 'ggnetwork', 'intergraph', 'ggsignif', 'rmarkdown', 'DT'))

### Description
The script `preprocess_peptideArray_table.R` preprocesses the peptide pull-down results table, 
which is at `./data/20170522_Neuroarray_results.tsv`. 

The second script `findSLiMDomainPairs.R` looks for motif gains/losses in mutant peptides with respect
to the wild-type peptides and associates theses changes to PFAM domains in the detected proteins as gained/lost 
interaction partners of the peptides. 

The rmarkdown script `peptideArray_manuscript_figures.Rmd` reproduces the figure S2B and supplementary data file Data S1 from the paper. 

The R script `render.vignette.R` is used to run the Rmarkdown script. 

### Step by step: how to run the peptide array analysis

Assuming the current directory as the top-level source directory. 

1. Preprocess the peptide pull-down result table
> Rscript vignettes/peptideArrayAnalysis/preprocess_peptideArray_table.R ./data

2. Associate slims to PFAM domains 

> Rscript vignettes/peptideArrayAnalysis/findSLiMDomainPairs.R ./data 

3. Render the manuscript figures/tables 

> Rscript ./vignettes/peptideArrayAnalysis/render.vignette.R \\
          ./vignettes/peptideArrayAnalysis/peptideArray_manuscript_figures.Rmd \\
          ./data


The output is a pdf file named `network_data.clustering_goterms.pdf` and an html file named `peptideArray_manuscript_figures.html`.

