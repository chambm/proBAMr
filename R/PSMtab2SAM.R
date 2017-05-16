#' Generate SAM files from confident peptide-spectrum-matches (PSMs).
#'
#' @title Generate SAM files from PSMs.
#'
#' @param passedPSM a data.frame of PSMs (created by one of proBAMr's read functions)
#' @param exon_anno a data.frame of exon annotations
#' @param proteinseq a data.frame of protein ids and protein sequences
#' @param procodingseq a data.frame of coding sequences for each protein
#' @param snvproseq a data.frame of protein sequence changes due to simple nucleotide variations
#' @param snvprocoding a data.frame of coding changes due to simple nucleotide variations (SNPs and INDELs)
#' @param junpep a data.frame of peptide sequences coded for by novel splice junctions
#' @param junpepcoding a data.frame of coding sequences coded for by novel splice junctions
#' @param jun_anno a data.frame of genome annotation information for novel splice junctions
#' @param outfile the filepath of the output SAM file
#' @param show_progress if TRUE, long-running parts of the function will report progress information
#' @param ... additional arguments
#'
#' @return invisible(NULL)
#' @author Xiaojing Wang
#' 
#' @importFrom stringr str_split fixed
#' @importFrom fastmatch fmatch
#' @import AhoCorasickTrie
#' @export
#' 
#' @examples
#' load(system.file("extdata/GENCODE", "exon_anno.RData", package="proBAMr"))
#' load(system.file("extdata/GENCODE", "proseq.RData", package="proBAMr"))
#' load(system.file("extdata/GENCODE", "procodingseq.RData", package="proBAMr"))
#' passedPSM <- readPepXmlTab(system.file("extdata/psm/passedPSM.tab", package="proBAMr"))
#' PSMtab2SAM(passedPSM, exon, proteinseq,  procodingseq, outfile=paste0(tempdir(), '/test.sam'))
#' SAM <- read.table(paste0(tempdir(), '/test.sam'))
#'
PSMtab2SAM <- function(passedPSM, exon_anno,
                       proteinseq, procodingseq,
                       snvproseq = NULL, snvprocoding = NULL,
                       junpep = NULL, junpepcoding = NULL, jun_anno = NULL,
                       outfile = "output.sam",
                       show_progress = TRUE, ...)
{
  old <- options(stringsAsFactors = FALSE)
  on.exit(options(old), add = TRUE)
  
  setDT(passedPSM, key="PeptideSequence") # convert to data.table for fast grouping
  # group by peptide so consecutive PSMs for the same peptide do not need to be remapped
  
  # abbreviate nativeId
  
  abbreviatedNativeId = stringi::stri_replace_last(gsub("\\s*\\S+=(\\S+)", "\\1.", passedPSM$NativeID),
                                                   "", fixed=".")
  passedPSM$SpectrumUniqueId = paste(passedPSM$Source, abbreviatedNativeId)
  passedPSM[, peptide_count:=length(unique(PeptideSequence)), SpectrumUniqueId] # peptides per spectrum
  passedPSM[, psm_count:=length(unique(SpectrumUniqueId)), PeptideSequence] # PSMs per peptide 
  passedPSM = as.data.frame(passedPSM) # back to data.frame because [.data.frame is more efficient

  # create data.tables for the annotation tables, keyed by ProteinAccessions
  proteinseqDT = data.table(proteinseq); setkey(proteinseqDT, "pro_name")
  procodingseqDT = data.table(procodingseq); setkey(procodingseqDT, "pro_name")
  
  # add DNA complement columns if they are missing
  if (is.null(procodingseq$dna_complement))
  {
    message("Optimizing procodingseq for fast reverse complement access...")
    procodingseqDT$dna_complement = sapply(procodingseqDT$coding, function(x) toString(reverseComplement(DNAString(x))))
    #save("procodingseq", file="procodingseq.RData")
  }

  # if there are no SNV annotations, use same structure as proteinseq but with 0 rows
  if (!is.null(snvproseq))
  {
    snvproseqDT = data.table(snvproseq); setkey(snvproseqDT, "pro_name")
    snvprocodingDT = data.table(snvprocoding); setkey(snvprocodingDT, "pro_name")
  }
  else
  {
    snvproseqDT = proteinseqDT[0,]
    snvprocodingDT = procodingseqDT[0,]
  }
  
  # some coding sequences have NA pro_name??
  procodingseqDT = unique(procodingseqDT[!is.na(pro_name), .SD])
  snvprocodingDT = unique(snvprocodingDT[!is.na(pro_name), .SD])

  # create joined tables with both coding and translated sequences for reference and variant proteins
  refproteinDT = proteinseqDT[procodingseqDT, .(pro_name, tx_name, tx_id, peptide, coding, dna_complement)]
  snvproteinDT = snvproseqDT[snvprocodingDT, .(pro_name, tx_name, tx_id, peptide, coding, dna_complement)]
  setkey(refproteinDT, tx_name, physical=F) # secondary index used for joining refprotein to snvprotein
  
  PEP <- passedPSM$PeptideSequence
  Spectrumid = passedPSM$SpectrumUniqueId
  
  normalProteins = proteinseqDT$peptide; names(normalProteins) = proteinseqDT$pro_name
  snvProteins = snvproteinDT$peptide; names(snvProteins) = snvproteinDT$pro_name
  junProteins = junpep$peptide; names(junProteins) = junpep$pro_name
  
  # Aho-Corasick trie maps all peptides to all proteins in a single sweep
  message("Mapping peptides to proteins...")
  feedbackInterval = 0L
  if (show_progress) feedbackInterval = 1000
  idxList = AhoCorasickSearchList(PEP, list(normal=normalProteins, snv=snvProteins, jun=junProteins),
                                  groupByKeyword=T
                                  #, iterationFeedback=feedbackInterval
                                  )
  idxNormal = lapply(idxList$normal, function(x) unlist(sapply(x, "[", "Text"), use.names=F))
  idxSNV = lapply(idxList$snv, function(x) unlist(sapply(x, "[", "Text"), use.names=F))
  idxJun = lapply(idxList$jun, function(x) unlist(sapply(x, "[", "Text"), use.names=F))
  
  normalProteinNames = data.table(names=names(normalProteins), key="names")
  snvProteinNames = data.table(names=names(snvProteins), key="names")
  if (is.null(names(junProteins))) junProteinNames = NULL
  else junProteinNames = data.table(names=names(junProteins), key="names")
  
  # AhoCorasickSearchList returns keyword (protein) names; convert those to indexes
  message("Calculating protein indexes...")
  idxNormal = lapply(idxNormal, function(x) unique(refproteinDT[x, which=T]))
  idxSNV = lapply(idxSNV, function(x) unique(snvproteinDT[x, which=T]))
  idxJun = lapply(idxJun, function(x) unique(junProteinNames[x, which=T]))
  
  message("Optimizing exon_anno structure for fast access...")
  tx_by_name = .exonAnnoByName(exon_anno)
  exon_anno_names = names(tx_by_name)
  attributes(exon_anno_names)$.match.hash = NULL # I got some crashes with fmatch() before doing this

  # loop body that processes a set of PSMs mapping to the same peptide; the
  # peptide is mapped to the genome, then that mapping info is written to a connection
  # in SAM format for each PSM
  processPSM = function(start, end, con="")
  {
    peptide <- PEP[start]
    peptide_ref <- peptide
    
    if(passedPSM$Rank[i] == 1) pri <- TRUE else pri <- FALSE
    idx = idxNormal[[peptide]]
    
    if(length(idx) > 0){
      pep_g <- 'N'
      res <- .Mapping4Peptide(idx, peptide, refproteinDT, tx_by_name, exon_anno_names, primary=pri)
    }else{
      idx_v = idxSNV[[peptide]]
      if(length(idx_v) > 0){
        pep_g <- 'V'
        res <- .Mapping4Peptide(idx_v, peptide, snvproteinDT, tx_by_name, exon_anno_names, primary=pri)
        peptide_ref <- .SAAVpeptideRef(idx_v, peptide, refproteinDT, snvproteinDT)
      }else{
        idx_j = idxJun[[peptide]]
        if(length(idx_j) > 0){
          pep_g <- 'J'
          res <- .Mapping4NovJunPeptide(idx_j, peptide, junpep, junpepcoding, jun_anno, primary=pri)
        }else{
          res <- .UnknownPeptide(primary=pri)
          if (passedPSM$IsDecoy[i])
          {
            pep_g <- 'D'
            res[[1]][1] = ifelse(pri, 0x400, 0x100+0x400)
          }else{
            pep_g <- 'U'
            res[[1]][1] = ifelse(pri, 0x4, 0x100+0x4)
          }
        }
      }
    }
    
    #if (length(res) > 0 && length(res[[1]]) != 13)
    #{
    #  print(start)
    #}
    
    NH <- length(res)
    XP <- peptide
    XR <- peptide_ref
    XG <- pep_g
    
    for (i in start:end)
    {
      QNAME <- Spectrumid[i]
      XL <- passedPSM$peptide_count[i]
      
      #XF <- paste('XF:f:', round(passedPSM[i, XFcolumn], digits=4), sep='')
      XC <- passedPSM$Charge[i]
      XS <- round(passedPSM$SearchEngineScore[i], digits=4)
      XQ <- format(passedPSM$FdrScore[i], scientific=TRUE)
      #XA <- paste('XA:Z:', annoted, sep='')
      
      XM <-  ifelse(is.na(passedPSM$PeptideModifications[i]), "-", passedPSM$PeptideModifications[i])
      XN <- passedPSM$MissedCleavages[i]
      XT <- passedPSM$NumSpecificTermini[i]
      
      numRows = length(res)
      while (numRows > 0)
      {
        cat(QNAME,
            '\t', res[[numRows]][[1]],
            '\t', res[[numRows]][[2]],
            '\t', res[[numRows]][[3]],
            '\t', res[[numRows]][[4]],
            '\t', res[[numRows]][[5]], 'M', # trailing M is implied in CIGAR
            '\t', res[[numRows]][[6]],
            '\t', res[[numRows]][[7]],
            '\t', res[[numRows]][[8]],
            '\t', substr(res[[numRows]][[9]], res[[numRows]][[10]], res[[numRows]][[11]]),
            '\t', res[[numRows]][[12]],
            '\t', res[[numRows]][[13]],
            '\tNH:i:', NH,
            '\tXL:i:', XL,
            '\tXP:Z:', XP,
            '\tXR:Z:', XR,
            '\tXC:i:', XC,
            '\tXS:f:', XS,
            '\tXQ:f:', XQ,
            '\tXM:Z:', XM,
            '\tXN:i:', XN,
            '\tXT:i:', XT,
            '\tXG:Z:', XG,
            '\n',
            file=con,
            sep="")
        numRows = numRows-1
      }
    }
    invisible(NULL)
  }
  
  message("Mapping peptides to genome...")
  con = file(outfile, "a+b")
  if (show_progress) pb = txtProgressBar(style=3, min=1, max=nrow(passedPSM))
  i=1
  psmCount = nrow(passedPSM)
  while (i < psmCount)
  {
    peptidePsmCount = passedPSM$psm_count[i]
    processPSM(i, min(psmCount, i+peptidePsmCount-1), con)
    i = min(psmCount, i+peptidePsmCount)
    if (show_progress) setTxtProgressBar(pb, i)
  }
  if (show_progress) close(pb)
  close(con)

  invisible(NULL)
}


.Mapping4Peptide <- function(idx, peptide, proteinDT,
                             exon_anno_by_tx_name, exon_anno_names, primary=TRUE, ...)
{
  numTx = length(idx)
  res = vector('list', numTx)
  
  pro_names = proteinDT$pro_name
  pro_peptides = proteinDT$peptide
  tx_names = proteinDT$tx_name
  codings = proteinDT$coding
  dna_complements = proteinDT$dna_complement
  pep_len <- nchar(peptide)
  
  while(numTx > 0)
  {
    pro_name = pro_names[[idx[[numTx]]]]
    pro_seq = pro_peptides[[idx[[numTx]]]]
    tx_name = tx_names[[idx[[numTx]]]]
    sta_pos <- regexpr(peptide, pro_seq, fixed=TRUE)
    end_pos <- sta_pos + pep_len - 1
    
    coding_seq = codings[[idx[[numTx]]]]
    code_s <- (sta_pos-1) * 3 + 1
    code_e <- end_pos * 3
    
    #################remove '^N' from the coding location
    # how many 'N' in the coding sequence and the position
    Ncoding_rmN <- gsub('N', '', coding_seq)
    numX <- nchar(coding_seq) - nchar(Ncoding_rmN)
    
    code_s <- code_s - numX
    code_e <- code_e - numX

    exon_anno = exon_anno_by_tx_name[[fmatch(tx_name, exon_anno_names)]]
    
    res[[numTx]] = .Coding2GenomeMapping(code_s, code_e,
                                         exon_anno,
                                         codings[[idx[[numTx]]]],
                                         dna_complements[[idx[[numTx]]]],
                                         pep_len, primary=primary)
    numTx = numTx-1
  }
  res[!duplicated(lapply(res, '[', c(2, 3)))]
}


.Coding2GenomeMapping <- function(x, y, exon_anno, d, c, pep_len, primary=TRUE, ...)
{
  if(is.null(exon_anno)) {
    .proteinUnannotated(x, y, exon_anno, d, c, primary=primary)
  }else{
    if(((y-x)+1 != 3*pep_len) || (y > max(exon_anno$cds_end, na.rm = TRUE))){
      #if(toString(translate(DNAString(m))) != peptide){
      .peptideUnannotated(x, y, exon_anno, d, c, primary=primary)
    }else{
      .MapCoding2genome(x, y, exon_anno, d, c, primary=primary)
    }
  }
}


###Key function to convert coding position to genomic position
.MapCoding2genome <- function(c_sta, c_end, exon_anno, dna, dna_complement, primary=TRUE, ...){
  
  # use a simple fast loop to find the exon the peptide starts and ends in
  eacs = exon_anno$cds_start
  eace = exon_anno$cds_end
  idxs = vector('integer', length(eacs))
  idxe = vector('integer', length(eace))
  j=0; k=0
  for(i in 1:length(eacs))
  {
    if (eacs[[i]] <= c_sta && c_sta <= eace[[i]])
    {
      j = j+1
      idxs[[j]] = i
    }
    
    if (eacs[[i]] <= c_end && c_end <= eace[[i]])
    {
      k = k+1
      idxe[[k]] = i
    }
  }
  idxs = idxs[1:j]
  idxe = idxe[1:k]

  len <- as.integer(c_end - c_sta + 1)
  RNAME <- as.character(exon_anno$chromosome_name[1])
  MAPQ <- 255L
  RNEXT <- '*'
  PNEXT <- 0L
  TLEN <- 0L
  QUAL <- '*'
  annoted <- 0
  XA <- 'XA:Z:0'
  CS <- as.integer(c_sta)
  CE <- as.integer(c_end)
  
  if(exon_anno$strand[[1]] == '+'){
    POS <- as.integer(exon_anno$cds_chr_start[[idxs]] + c_sta - exon_anno$cds_start[[idxs]])
    SEQ <- dna
    FLAG <- ifelse(primary==TRUE, 0x00, 0x00+0x100)
    
  }else{
    POS <- as.integer(exon_anno$cds_chr_start[[idxe]] + exon_anno$cds_end[[idxe]] - c_end)
    SEQ <- dna_complement
    FLAG <- ifelse(primary==TRUE, 0x10, 0x10+0x100)
  }
  
  # when generating the CIGAR string, leave off the final 'M';
  # this way, for single-exon peptides (the majority) no memory needs to be allocated
  if(idxe == idxs){
    CIGAR <- len
  }else{
    if(exon_anno$strand[[1]] == '+'){
      #insert <- exon_anno[idxe, 'cds_chr_start'] - exon_anno[idxs, 'cds_chr_end']- 1
      part1 <- exon_anno$cds_end[[idxs]] - c_sta + 1
      part2 <- c_end - exon_anno$cds_start[[idxe]] + 1
      
      insert <- unlist(lapply(1:(idxe-idxs), function(x)
        paste(exon_anno$cds_chr_start[[idxs+x]]-exon_anno$cds_chr_end[[idxs+x-1]]- 1,
              'N', sep='')))
      if(idxe-idxs >1){
        innerexon <- unlist(lapply(1:(idxe-idxs-1), function(x)
          paste(exon_anno$cds_chr_end[[idxs+x]] -
                  exon_anno$cds_chr_start[[idxs+x]]+1, 'M', sep='')))
      }else{ innerexon <-''}
      
      #     ifelse(idxe-idxs >1, unlist(lapply(1:(idxe-idxs-1), function(x)
      #     paste(exon_anno[idxs+x, 'cds_chr_end']-exon_anno[idxs+x, 'cds_chr_start']+1,
      #             'M', sep=''))), '')
      midpattern <- rep(NA, length(insert)+length(innerexon))
      midpattern[seq(1, length(insert)+length(innerexon), by=2)] <- insert
      midpattern[seq(2, length(insert)+length(innerexon), by=2)] <- innerexon
      midpattern <- paste(midpattern, collapse='')
      
    }else{
      #insert <- exon_anno[idxs, 'cds_chr_start'] - exon_anno[idxe, 'cds_chr_end']- 1
      part1 <- c_end- exon_anno$cds_start[idxe] + 1
      part2 <- exon_anno$cds_end[idxs] - c_sta + 1
      
      insert <- unlist(lapply(1:(idxe-idxs), function(x)
        paste(exon_anno$cds_chr_start[[idxe-x]]-exon_anno$cds_chr_end[[idxe-x+1]]- 1,
              'N', sep='')))
      if(idxe-idxs >1){
        innerexon <- unlist(lapply(1:(idxe-idxs-1), function(x)
          paste(exon_anno$cds_chr_end[[idxe-x]] -
                  exon_anno$cds_chr_start[[idxe-x]]+1, 'M', sep='')))
      }else{ innerexon <-''}
      
      midpattern <- rep(NA, length(insert)+length(innerexon))
      midpattern[seq(1, length(insert)+length(innerexon), by=2)] <- insert
      midpattern[seq(2, length(insert)+length(innerexon), by=2)] <- innerexon
      midpattern <- paste(midpattern, collapse='')
      
    }
    
    CIGAR <- paste0(part1, 'M', midpattern, part2)
  }
  
  tmp <- list(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, CS, CE, QUAL, XA)
  tmp
}

###Peptide could map to a specific protein, but the protein
###don't have genomic annotation
.proteinUnannotated <- function(c_sta, c_end, exon_anno, dna, dna_complement, primary=TRUE, ...)
{
  RNAME <- '*'
  MAPQ <- 255L
  RNEXT <- '*'
  PNEXT <- 0L
  TLEN <- 0L
  QUAL <- '*'
  POS <- 0L
  SEQ <- '*'
  CS <- 0L
  CE <- 0L
  CIGAR <- '*'
  annoted <- 2
  XA <- 'XA:Z:2'
  #FLAG <- 0x4
  FLAG <- ifelse(primary==TRUE, 0x4L, 0x4L+0x100L)
  
  tmp <- list(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, CS, CE, QUAL, XA)
  tmp
}


###Peptide could map to a specific protein, but the corresponding coding
###sequence start or end position are outside the genomic annotation
.peptideUnannotated <- function(c_sta, c_end, exon_anno, dna, dna_complement, primary=TRUE, ...)
{
  
  #RNAME <- as.character(exon_anno[1, 'chromosome_name'])
  RNAME <- '*'
  MAPQ <- 255L
  RNEXT <- '*'
  PNEXT <- 0L
  TLEN <- 0L
  QUAL <- '*'
  POS <- 0L
  SEQ <- '*'
  CS <- 0L
  CE <- 0L
  CIGAR <- '*'
  annoted <- 1
  XA <- 'XA:Z:1'
  #FLAG <- 0x4
  FLAG <- ifelse(primary==TRUE, 0x4L, 0x4L+0x100L)
  #if(exon_anno[1, 'strand'] == '+'){
  #    FLAG <- ifelse(primary==T, 0x00, 0x00+0x100)
  #}else{
  #    FLAG <- ifelse(primary==T, 0x10, 0x10+0x100)
  #}
  
  tmp <- list(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, CS, CE, QUAL, XA)
  tmp
}


.SAAVpeptideRef <- function(idx_v, peptide, refproteinDT, snvproteinDT)
{
  dontCheck(snvpro <- snvproteinDT[idx_v, .(pro_name, tx_name, peptide)])
  dontCheck(refpeptide <- refproteinDT[list(tx_name=snvpro$tx_name), peptide, on="tx_name"])
  sta_pos <- unlist(lapply(snvpro$peptide, function(x) regexpr(peptide, x, fixed=TRUE)))
  pep_len <- nchar(peptide)
  end_pos <- sta_pos + pep_len - 1
  peptide_ref <- unique(substring(refpeptide, sta_pos, end_pos))
  peptide_ref
}

.Mapping4NovJunPeptide <- function(idx_j, peptide, junpep, junpepcoding,
                                   jun_anno, primary=TRUE, ...)
{
  pro <- junpep[idx_j, ]
  sta_pos <- unlist(lapply(pro[, 'peptide'], function(x)
    regexpr(peptide, x, fixed=TRUE)))
  pep_len <- nchar(peptide)
  end_pos <- sta_pos + pep_len - 1
  
  coding <- junpepcoding[match(substr(pro[, 'pro_name'], 1,
                                      nchar(pro[, 'pro_name'])-5),
                               junpepcoding[, 'pro_name']), ]
  orf <- substring(pro[, 'pro_name'],  nchar(pro[, 'pro_name']))
  if(orf=='1') mvf <- 0
  if(orf=='2') mvf <- 1
  if(orf=='3') mvf <- 2
  code_s <- (sta_pos-1) * 3 + 1 + mvf
  code_e <- end_pos * 3 + mvf
  codingseq <- substring(coding[, 'coding'], code_s, code_e)
  
  exonp <- lapply(substr(pro[, 'pro_name'], 1,
                         nchar(pro[, 'pro_name'])-5), function(x)
                           jun_anno[jun_anno[, 'pro_name']==x, ])
  
  
  #if(passedPSM[i, 'hit_rank'] == 1) pri <- TRUE else pri <- FALSE
  
  res <- mapply(function(x, y, z, m)
    .Coding2GenomeMapping(x, y, z, m, pep_len, primary=primary),
    code_s, code_e, exonp, codingseq)
  
  res <-  unique(data.frame(t(res)))
  res
}


#################################################################
# Decide whether the coding sequence is consistent with protein
# sequence, if not, use tblastn to find the coding start/end postion
#################################################################
# .CorrectCodingPosTblastn <- function(...){
#   #if(FALSE %in% lapply(translate(DNAStringSet(codingseq)),
#   #                        function(x) toString(x) == peptide)){
#   
#   pro <- proteinseq[idx, ]
#   write(paste('>peptide', '\n', peptide, sep=''),
#         file=paste(tmpfolder, '/pep.fasta', sep=''))
#   
#   write(paste(pro[, 'pro_name'], collapse='\n'),
#         file=paste(tmpfolder, '/seqid', sep=''))
#   
#   cmd <- paste('tblastn -db ', blastdb,
#                ' -query ', paste(tmpfolder, '/pep.fasta', sep=''),
#                ' -outfmt 7 -seqidlist ',
#                file=paste(tmpfolder, '/seqid', sep=''), sep='')
#   
#   myPipe <- pipe(cmd)
#   results <- try(read.table(myPipe), TRUE)
#   if(!class(results)=='try-error') {
#     colnames( results ) <- c( "QueryID",  "SubjectID",
#                               "Perc.Ident",
#                               "Alignment.Length", "Mismatches", "Gap.Openings",
#                               "Q.start",
#                               "Q.end", "S.start", "S.end", "E", "Bits" )
#     
#     code_s <- results[match(pro[, 'pro_name'],
#                             results[, 'SubjectID']), 'S.start']
#     code_e <- results[match(pro[, 'pro_name'],
#                             results[, 'SubjectID']), 'S.end']
#     codingseq <- substring(coding[, 'coding'], code_s, code_e)
#   }
#   
#   
#   #            }
# }

### For some reason, peptide cannot matach to any protein sequence
.UnknownPeptide <- function(primary=TRUE,...){
  
  RNAME <- '*'
  MAPQ <- 255L
  RNEXT <- '*'
  PNEXT <- 0L
  TLEN <- 0L
  QUAL <- '*'
  POS <- 0L
  CS <- 0L
  CE <- 0L
  SEQ <- '*'
  CIGAR <- '*'
  annoted <- 2
  XA <- 'XA:Z:2'
  #FLAG <- 0x4
  FLAG <- ifelse(primary==TRUE, 0x4L, 0x4L+0x100L)
  
  tmp <- list(list(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, CS, CE, QUAL, XA))
  tmp
}

.exonAnnoByName <- function(exon)
{
  setDT(exon)
  
  # reduce exon annotation down to one tx_id per tx_name;
  # drop the rows without cds_chr_start information
  setkey(exon, tx_name, tx_id)
  firstTxId <- exon[, .(tx_id=min(tx_id)), tx_name]
  setkey(firstTxId, tx_name, tx_id)
  uniqueTranscripts <- exon[firstTxId, .(tx_id, tx_name, strand, chromosome_name,
                                         cds_chr_start, cds_chr_end,
                                         cds_start, cds_end)][!is.na(cds_chr_start)]
  
  tx_names = unique(uniqueTranscripts$tx_name)
  
  # create a list with separate data.frames for each tx_name so that we don't have to constantly
  # create slices 
  tx_by_name = lapply(tx_names, function(x) uniqueTranscripts[x])
  
  names(tx_by_name) = tx_names
  
  tx_by_name
}
