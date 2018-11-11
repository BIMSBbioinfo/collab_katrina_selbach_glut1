library('data.table')
library('parallel')

datadir <- '/data/local/buyar/collaborations/katrina/collab_katrina_selbach_glut1/data'

#get uniprot accessions and corresponding gene names from the fasta file of human uniprot sequences
fasta <- Biostrings::readAAStringSet(file.path(datadir, 'uniprot.9606.Nov_8_2018.fasta.gz'))

uniprot2geneName <- data.frame('uniprotAccession' = sub("^sp\\|(.+?)\\|.+$", "\\1", names(fasta)),
                               'geneName' = sub("^.+ GN=(.+?) .+$", "\\1", names(fasta)),
                               stringsAsFactors = FALSE)

#import peptide-array results table from Katrina
neuroarrayResultsFile <- file.path(datadir, '20170522_Neuroarray_results.tsv')
dt <- fread(neuroarrayResultsFile)

#reviewed uniprot accesssions
uniReviewed <- uniprot2geneName$uniprotAccession

#only keep review uniprot accessions
cl <- parallel::makeCluster(10)
clusterExport(cl = cl, varlist = c('dt', 'uniReviewed'))
match <- do.call(c, parLapply(cl = cl, X = 1:nrow(dt),
          fun = function(i) {
            ids <- unique(gsub(replacement = '',
                               pattern = '-[0-9]*$',
                               x = unlist(strsplit(dt[i,]$Majority.protein.IDs.prey, split = ';'))))
            if(sum(ids %in% uniReviewed) > 0) {
              return(paste(ids[ids %in% uniReviewed], collapse = ';'))
            } else {
              return(paste(ids, collapse = ';'))
            }
          }))
stopCluster(cl)

dt$uniprotMatch <- match
dt <- dt[uniprotMatch != ''] #remove rows with unmapped ids
dt$genotype <- gsub(pattern = '[0-9]+$', replacement = '', x = dt$ExperimentID)

#further process dt to assign unique ids to each peptide (wt or mutant)
dt$PeptideIDext <- ''
dt[genotype == 'wt']$PeptideIDext <- paste(dt[genotype == 'wt']$PeptideID,
                                                  dt[genotype == 'wt']$wild.type.sequence,
                                                  'wt', sep = '_')
dt[genotype == 'mut']$PeptideIDext <- paste(dt[genotype == 'mut']$PeptideID,
                                                   dt[genotype == 'mut']$mutant.sequence,
                                                   'mut', sep = '_')
#Define maxium SILAC ratio field (we have 2 replicates for each, so multiplying median
#-actually it is the mean value- by 2 and substracting minimum value, gives the maximum value)
dt$`Maximum.SILAC.ratio.Wt/Mut` <- dt$`Median.SILAC.ratio.Wt/Mut` * 2 - dt$`Minimum.SILAC.ratio.Wt/Mut`

#create peptide ids replacing uniprot ids with gene names
dt$PeptideUniprotID <- gsub('_.*$', '', dt$PeptideID)
dt$PeptideVariant <- gsub('^.*_', '', dt$PeptideID)
dt$PeptideGeneName <- uni2gene[match(dt$PeptideUniprotID, uni2gene$UNIPROTKB),]$`ENTRY-NAME`
dt[is.na(PeptideGeneName)]$PeptideGeneName <- dt[is.na(PeptideGeneName)]$PeptideUniprotID
dt$PeptideIDwithGeneName <- paste(dt$PeptideGeneName, dt$PeptideVariant, sep = '_')
#remove ARHGAP36_cntrl_1 peptides from the analysis
dt <- dt[PeptideID != 'ARHGAP36_cntrl_1']
dt$interactionID <- c(1:nrow(dt))

outFile <- gsub(x = neuroarrayResultsFile, pattern = '.tsv$', replacement = '.preprocessed.tsv')
write.table(dt, file = outFile, quote = FALSE, sep = '\t', row.names = FALSE)
