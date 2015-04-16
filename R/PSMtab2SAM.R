##' Generate SAM files from confident peptide-spectrum-matches (PSMs). 
##' 
##' @title Generate SAM files from PSMs.
##' @param passedPSM a data frame of PSMs passed FDR.
##' @param XScolumn specify the column which represents the matching score.
##' @param exon_anno a dataframe of exon annotations.
##' @param proteinseq a dataframe containing  protein ids and protein 
##' sequences.
##' @param procodingseq a data frame cotaining coding sequence for each 
##' protein.
##' @param ... additional arguments
##' @return a dataframe containing 
##' @author Xiaojing Wang
##' @examples
##' 
##' load(system.file("extdata/GENCODE", "exon_anno.RData", package="proBAMr"))
##' load(system.file("extdata/GENCODE", "proseq.RData", package="proBAMr"))
##' load(system.file("extdata/GENCODE", "procodingseq.RData", 
##'     package="proBAMr"))
##' options(stringsAsFactors=FALSE)
##' passedPSM <- read.table(system.file("extdata", "passedPSM.tab", 
##'     package="proBAMr"), sep='\t', header=TRUE)
##' SAM <- PSMtab2SAM(passedPSM, XScolumn='mvh', exon, proteinseq, 
##'     procodingseq)
##' write.table(SAM, file=paste(tempdir(), '/test.sam', sep=''), 
##'             sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
##' 

PSMtab2SAM <- function(passedPSM, XScolumn='mvh', exon_anno, proteinseq, 
    procodingseq, ...)
{
    options(stringsAsFactors=FALSE)
    #passedPSM
    PEP <- passedPSM[, 'peptide']
    Spectrumid <- apply(passedPSM, 1, function(x){
        paste(strsplit(x['spectrum'], '\\.')[[1]][1], 
        strsplit(strsplit(x['spectrumNativeID'], ' ')[[1]][1], '=')[[1]][2], 
        strsplit(strsplit(x['spectrumNativeID'], ' ')[[1]][2], '=')[[1]][2], 
        strsplit(strsplit(x['spectrumNativeID'], ' ')[[1]][3], '=')[[1]][2], 
        sep='.')
        })
    #PEP_SEQ <- formatPep(spectra[, 'Sequence'])

    SAM <- c()

    spectrumcount <- table(Spectrumid)

    for(i in 1:dim(passedPSM)[1]){
        #print(i)
        peptide <- PEP[i]
        QNAME <- Spectrumid[i]
        idx <- grep(peptide, proteinseq[, 'peptide'], fixed=TRUE)
        if(length(idx) == 0){
            RNAME <- '*' 
            MAPQ <- 255
            RNEXT <- '*'
            PNEXT <- 0
            TLEN <- 0
            QUAL <- '*'
            POS <- 0
            SEQ <- '*'
            CIGAR <- '*'            
            FLAG <- 0x4 
            annoted <- '?'
            XA <- paste('XA:Z:', annoted, sep='')
            res <- c(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, 
                    as.character(SEQ), QUAL, XA)
            res <-  unique(data.frame(t(res))) 
        }else{
            pro <- proteinseq[idx, ]
            sta_pos <- unlist(lapply(pro[, 'peptide'], function(x) 
                        regexpr(peptide, x, fixed=TRUE)))
            pep_len <- nchar(peptide)
            end_pos <- sta_pos + pep_len - 1
            
            coding <- procodingseq[match(pro[, 'pro_name'], 
                                        procodingseq[, 'pro_name']), ]
            code_s <- (sta_pos-1) * 3 + 1
            code_e <- end_pos * 3
            codingseq <- substring(coding[, 'coding'], code_s, code_e)
            
            
            
            exonp <- lapply(pro[, 'tx_name'], function(x)
                        exon_anno[exon_anno[, 'tx_name']==x, ])
            
            exonp <- lapply(exonp, function(x){
                            if(length(unique(x[, 'tx_id'])) > 1){
                                x[grep(x[1, 'tx_id'], x[, 'tx_id'], 
                                fixed=TRUE), ]
                            }else x
                        })
                       
            if(passedPSM[i, 'hit_rank'] == 1) pri <- TRUE else pri <- FALSE
                
            res <- mapply(function(x, y, z, m) 
                    if(dim(z)[1] == 0){
                        .proteinUnannotated(x, y, z, m, primary=pri)
                    }else{
                        if((nchar(m) != 3*pep_len) | (y > max(z[, 'cds_end'], 
                            na.rm = TRUE))){
                        #if(toString(translate(DNAString(m))) != peptide){
                            .peptideUnannotated(x, y, z, m, primary=pri)
                        }else{
                            .MapCoding2genome(x, y, z, m, primary=pri)
                        }
                    }, 
                    code_s, code_e, exonp, codingseq) 
                    
            res <-  unique(data.frame(t(res)))
            
        }
        XL <- paste('XL:i:', as.numeric(spectrumcount[QNAME]), sep='')
        NH <- paste('NH:i:', dim(res)[1], sep='')
        XP <- paste('XP:Z:', peptide, sep='')
        #XF <- paste('XF:f:', round(passedPSM[i, XFcolumn], digits=4), sep='')
        XC <- paste('XC:i:', passedPSM[i, 'assumed_charge'], sep='')
        XS <- paste('XS:f:', round(as.numeric(passedPSM[i, XScolumn]), 
                        digits=4), sep='')
        #XA <- paste('XA:Z:', annoted, sep='')
        XN <- paste('XN:i:', passedPSM[i, 'num_missed_cleavages'], sep='')
        XT <- paste('XT:i:', passedPSM[i, 'NTT'], sep='')
        
        if(length(idx) > 0){
            pep_g <- 'N'
        }
        
        XG <- paste('XG:Z:', pep_g, sep='')
        
        XM <-  ifelse(is.na(passedPSM[i, 'modification']), paste('XM:Z:-'),
                    paste('XM:Z:', passedPSM[i, 'modification'], sep=''))
        
        res <- cbind(QNAME, res, NH, XL, XP, XC, XS, XM, XN, XT, XG)    
        SAM <- rbind(SAM, res)
    }

    SAM

}



.proteinUnannotated <-function(c_sta, c_end, exon_anno, cseq, primary=TRUE, 
...)
{ 
    RNAME <- '*' 
    MAPQ <- 255
    RNEXT <- '*'
    PNEXT <- 0
    TLEN <- 0
    QUAL <- '*'
    POS <- 0
    SEQ <- '*'
    CIGAR <- '*'
    annoted <- 2
    XA <- paste('XA:Z:', annoted, sep='')
    FLAG <- 0x4 
    
    tmp <- c(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, 
                as.character(SEQ), QUAL, XA)
    tmp
}



.peptideUnannotated <- function(c_sta, c_end, exon_anno, cseq, primary=TRUE, 
...)
{    
    #RNAME <- as.character(exon_anno[1, 'chromosome_name']) 
    RNAME <- '*'
    MAPQ <- 255
    RNEXT <- '*'
    PNEXT <- 0
    TLEN <- 0
    QUAL <- '*'
    POS <- 0
    SEQ <- '*'
    CIGAR <- '*'
    annoted <- 1
    XA <- paste('XA:Z:', annoted, sep='')
    FLAG <- 0x4 
    #if(exon_anno[1, 'strand'] == '+'){        
    #    FLAG <- ifelse(primary==TRUE, 0x00, 0x00+0x100) 
    #}else{
    #    FLAG <- ifelse(primary==TRUE, 0x10, 0x10+0x100)   
    #}
    
    tmp <- c(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, 
                as.character(SEQ), QUAL, XA)
    tmp
}




.MapCoding2genome <- function(c_sta, c_end, exon_anno, cseq, primary=TRUE, 
...)
{
    idxs <- intersect(which(exon_anno[, 'cds_start'] <= c_sta), 
                        which(exon_anno[, 'cds_end'] >= c_sta))
    idxe <- intersect(which(exon_anno[, 'cds_start'] <= c_end), 
                        which(exon_anno[, 'cds_end'] >= c_end))              
    len <- c_end - c_sta + 1 
    RNAME <- as.character(exon_anno[1, 'chromosome_name']) 
    MAPQ <- 255
    RNEXT <- '*'
    PNEXT <- 0
    TLEN <- 0
    QUAL <- '*'
    annoted <- 0
    XA <- paste('XA:Z:', annoted, sep='')
    
    if(exon_anno[1, 'strand'] == '+'){
            POS <- exon_anno[idxs, 'cds_chr_start'] + c_sta - 
                    exon_anno[idxs, 'cds_start']
            SEQ <- DNAString(cseq)
            FLAG <- ifelse(primary==TRUE, 0x00, 0x00+0x100) 
            
     }else{
            POS <- exon_anno[idxe, 'cds_chr_start'] +  
                    exon_anno[idxe, 'cds_end'] - c_end
            SEQ <- reverseComplement(DNAString(cseq))
            FLAG <- ifelse(primary==TRUE, 0x10, 0x10 + 0x100) 
     }
    
    if(idxe == idxs){
        CIGAR <- paste(len, 'M', sep='')        
    }else{
         if(exon_anno[1, 'strand'] == '+'){
            #insert <- exon_anno[idxe, 'cds_chr_start'] - exon_anno[idxs, 
            #                        'cds_chr_end']- 1
            part1 <- exon_anno[idxs, 'cds_end'] - c_sta + 1 
            part2 <- c_end - exon_anno[idxe, 'cds_start'] + 1  
            
            insert <- unlist(lapply(1:(idxe - idxs), function(x) 
                paste(exon_anno[idxs + x, 'cds_chr_start'] - 
                    exon_anno[idxs+x-1, 'cds_chr_end']- 1, 'N', sep='')))
            if(idxe-idxs >1){
                innerexon <- unlist(lapply(1:(idxe-idxs-1), function(x) 
                        paste(exon_anno[idxs + x, 'cds_chr_end'] - 
                        exon_anno[idxs + x, 'cds_chr_start'] + 1, 'M', 
                        sep=''))) 
            }else{ innerexon <- ''}
            
            #ifelse(idxe-idxs >1, unlist(lapply(1:(idxe-idxs-1), function(x) 
            #paste(exon_anno[idxs+x, 'cds_chr_end'] - 
            #             exon_anno[idxs+x, 'cds_chr_start']+1, 
            #             'M', sep=''))), '')
            midpattern <- rep(NA, length(insert) + length(innerexon))
            midpattern[seq(1, length(insert) + length(innerexon), 
                                                        by=2)] <- insert
            midpattern[seq(2, length(insert) + length(innerexon), 
                                                        by=2)] <- innerexon
            midpattern <- paste(midpattern, collapse='') 
            
        }else{
            #insert <- exon_anno[idxs, 'cds_chr_start'] - 
            #            exon_anno[idxe, 'cds_chr_end']- 1
            part1 <- c_end- exon_anno[idxe, 'cds_start'] + 1 
            part2 <- exon_anno[idxs, 'cds_end'] - c_sta + 1  
            
            insert <- unlist(lapply(1:(idxe-idxs), function(x) 
                paste(exon_anno[idxe - x, 'cds_chr_start'] - 
                        exon_anno[idxe-x + 1, 'cds_chr_end']- 1, 
                        'N', sep='')))
            if(idxe-idxs >1){
                innerexon <- unlist(lapply(1:(idxe-idxs-1), function(x) 
                        paste(exon_anno[idxe-x, 'cds_chr_end'] - 
                        exon_anno[idxe-x, 'cds_chr_start']+1, 'M', sep='')))
            }else{ innerexon <-''}
        
            midpattern <- rep(NA, length(insert)+length(innerexon))
            midpattern[seq(1, length(insert) + length(innerexon), 
                                                        by=2)] <- insert
            midpattern[seq(2, length(insert) + length(innerexon), 
                                                        by=2)] <- innerexon
            midpattern <- paste(midpattern, collapse='') 

        }

        CIGAR <- paste(part1, 'M', midpattern, part2, 'M', sep='')    
    }

    tmp <- c(FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, 
                as.character(SEQ), QUAL, XA)
    tmp
}






