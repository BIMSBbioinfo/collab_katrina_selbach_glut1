library(data.table)
library(parallel)

args = commandArgs(trailingOnly=TRUE)

#set progres bar options
pbapply::pboptions(type = 'timer')

datadir <- args[1] #path to the `data` folder

#1. read peptide array interaction data
message("importing processed peptide array table")
arrayData <- fread(file.path(datadir, '20170522_Neuroarray_results.preprocessed.tsv'))

#2. get a list of all uniprot accessions found as interaction partners of the peptides in the array data
uniprotAccessions <- unique(arrayData$uniprotMatch) #some of them are semi-column separated.
uniprotAccessions <- unique(unlist(lapply(uniprotAccessions, function(x) {strsplit(x = x, split = ';')})))

#3. now associate each uniprot accession with PFAM domains (which pfam domains are found in which sequence?)
# pfam domain predictions (based on scanning PFAM hmms) on human protein sequences can be downloaded from here:
# 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//releases/Pfam30.0/proteomes/9606.tsv.gz'

message("importing PFAM domains")
uniprot2pfam <- as.data.table(slimR::getPFAM(organism = 9606, pfam_version = "Pfam30.0"))
uniprot2pfam <- uniprot2pfam[uniprot2pfam$seqnames %in% uniprotAccessions,]
uniprot2pfam <- unique(subset(uniprot2pfam, select = c('seqnames', 'pfam_acc', 'pfam_name', 'clan')))
colnames(uniprot2pfam) <- c('uniprotAccession', 'PFAM', 'NAME', 'CLAN')

message("importing ELM -> PFAM associations")
#download pfam clan information
pfamClans <- slimR::getPFAMClans()

# download a table of ELM classes and their cognate PFAM domains
elm2pfam <- slimR::getELMdomainInteractions()

#5. For each peptide (WT or MUT) in the arrayData table,
# find out if the interaction partners of the peptide contain
# a domain that can bind to a motif found in the peptide

#' @param arrayData peptide array data
#' @param elm2pfam Table of associations between ELM classes and PFAM classes (see getELMdomainInteractions function)
#' @param uniprot2pfam Table of uniprot accessions and pfam domains contained in each protein sequence.
#' @param nodeN Number of cpus to use for parallel operationss
associatePeptides2PfamDomains <- function(arrayData,
                                          uniprot2pfam,
                                          elm2pfam,
                                          nodeN = 10) {
  #Collect all peptides in a list as WT and MUT peptides and search for SLIM pattern matches in each peptide,
  pepSeqs <- unique(arrayData[,c('PeptideID', 'wild.type.sequence', 'mutant.sequence')])
  wtPeptides <- as.list(pepSeqs$wild.type.sequence)
  names(wtPeptides) <- pepSeqs$PeptideID

  mutPeptides <- as.list(pepSeqs$mutant.sequence)
  names(mutPeptides) <- pepSeqs$PeptideID

  #scan the peptide sequences to find out which slim classes match the given sequences
  message("scanning WT peptides for slims")
  wtSlims <- pbapply::pblapply(X = wtPeptides,
                       FUN = function(x) {
                         #avoid motifs at N' or C terminal because
                         #we don't consider the whole protein sequences
                         #but just a peptide
                         x <- paste0('XXX',x,'XXX')
                         slims = slimR::searchSLiMs(x, slimR::motifRegex)
                         return(unique(slims$SLiM))
                         })
  message("scanning mutant peptides for slims")
  mutSlims <- pbapply::pblapply(X = mutPeptides,
                               FUN = function(x) {
                                 #avoid motifs at N' or C terminal because
                                 #we don't consider the whole protein sequences
                                 #but just a peptide
                                 x <- paste0('XXX',x,'XXX')
                                 slims = slimR::searchSLiMs(x, slimR::motifRegex)
                                 return(unique(slims$SLiM))
                               })


  message("Looking through detected interactions in arrayData")
  #now for each interaction found in the arrayData table, find out if the
  #interaction partner contains a domain. If the interaction partner contains a
  #domain, find out if the domain can bind to any of the slims available in the
  #peptide sequences for both WT and MUT peptides.
  cl <- makeCluster(nodeN)
  clusterExport(cl = cl,
                varlist = c('arrayData', 'uniprot2pfam', 'elm2pfam', 'wtSlims', 'mutSlims'),
                envir = environment())
  df <- do.call(rbind, pbapply::pblapply(cl = cl, X = 1:nrow(arrayData), function(i){
    require('data.table')
    pepId <- arrayData[i]$PeptideID
    #one partner may be identified with multiple uniprot accessions separated by a ';'.
    interactor <- unlist(strsplit(x = arrayData[i]$uniprotMatch, split = ';'))
    domains <- unique(unlist(strsplit(x = uniprot2pfam[uniprot2pfam$uniprotAccession %in% interactor,]$PFAM,
                                      split = ';')))

    do.call(rbind, lapply(X = domains,
                          FUN = function(d) {
                            domainTargets <- as.character(elm2pfam[elm2pfam$Interaction_Domain_Id == d,]$ELM_identifier)
                            domainTargetsInPeptideWT <- intersect(domainTargets, wtSlims[[pepId]])
                            domainTargetsInPeptideMUT <- intersect(domainTargets, mutSlims[[pepId]])
                            lostTargets <- setdiff(domainTargetsInPeptideWT, domainTargetsInPeptideMUT)
                            gainedTargets <- setdiff(domainTargetsInPeptideMUT, domainTargetsInPeptideWT)
                            data.frame('interactionID' = i,
                                       'PeptideID' = pepId,
                                       'interactor' = paste0(interactor, collapse = ';'),
                                       'domain' = d,
                                       'domainTargetsInPeptideWT' = paste0(domainTargetsInPeptideWT, collapse = ';'),
                                       'domainTargetsInPeptideMUT' = paste0(domainTargetsInPeptideMUT, collapse = ';'),
                                       'lostTargets' = paste0(lostTargets, collapse = ';'),
                                       'gainedTargets' = paste0(gainedTargets, collapse = ';'),
                                       stringsAsFactors = FALSE)
                          }
                          ))
    }))
  stopCluster(cl = cl)
  return(data.table(df))
}

#get slim-domain interactions for wt and mutant peptides
#and find out which peptides have gained/lost slim-domain interactions
#that can be explained by gained slims and existing cognate domains in the interaction partneres
#of the mutant form of the peptides.
message("Associating slims to PFAM domains")
slimDomainInteractions <- associatePeptides2PfamDomains(arrayData = arrayData,
                                                        uniprot2pfam = uniprot2pfam,
                                                        elm2pfam = elm2pfam,
                                                        nodeN = 10)
slimDomainInteractions$domainName <- uniprot2pfam[match(slimDomainInteractions$domain, uniprot2pfam$PFAM),]$NAME
slimDomainInteractions$clan <- uniprot2pfam[match(slimDomainInteractions$domain, uniprot2pfam$PFAM),]$CLAN
slimDomainInteractions$clanName <- as.character(pfamClans[match(slimDomainInteractions$clan, pfamClans$Accession),]$ID)

message("Saving peptideArray_slimDomainInteractions.RDS at ",datadir)
saveRDS(object = slimDomainInteractions, file = file.path(datadir, 'peptideArray_slimDomainInteractions.RDS'))
