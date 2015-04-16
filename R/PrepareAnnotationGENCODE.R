##' prepare the annotation from GENCODE. Download GTF and FASTA files from 
##' GENCODE ftp first. Read introduction for more information.
##' 
##' @title prepare annotation from GENCODE
##' @param gtfFile specify GTF file location.
##' @param CDSfasta path to the fasta file of coding sequence.
##' @param pepfasta path to the fasta file of protein sequence.
##' @param annotation_path specify a folder to store all the annotations.
##' @param dbsnp specify a snp dataset to be used for the SNP annotation, 
##' default is NULL. (e.g. "snp135")
##' @param splice_matrix whether generate a known exon splice matrix from the 
##' annotation. this is not necessary if you don't want to analyse junction 
##' results, default is FALSE. 
##' @param COSMIC whether to download COSMIC data, default is FALSE.
##' @param ... additional arguments
##' @return several .RData files containing annotations needed for further 
##' analysis.
##' @author Xiaojing Wang
##' @examples
##' 
##' gtfFile <- system.file("extdata", "test.gtf", package="proBAMr")
##' CDSfasta <- system.file("extdata", "coding_seq.fasta", package="proBAMr") 
##' pepfasta <- system.file("extdata", "pro_seq.fasta", package="proBAMr") 
##' annotation_path <- tempdir()
##' PrepareAnnotationGENCODE(gtfFile, CDSfasta, pepfasta, 
##'                 annotation_path, dbsnp=NULL, 
##'                 splice_matrix=FALSE, COSMIC=FALSE)  
##' 


PrepareAnnotationGENCODE <- function(gtfFile, CDSfasta, pepfasta, 
                annotation_path, dbsnp=NULL, 
                splice_matrix=FALSE, COSMIC=FALSE, 
 ...) {
    options(stringsAsFactors=FALSE)
    if (!is.null(dbsnp)) {
        session  <- browserSession()
        dbsnps <- trackNames(session)[grep('snp', trackNames(session), 
        fixed=TRUE)]
        dbsnp <- pmatch(dbsnp, dbsnps)
        if (is.na(dbsnp)) 
                stop("invalid dbsnp name for specified genome")
        if (dbsnp == -1) 
                stop("ambiguous dbsnp name")
    }

    options(stringsAsFactors=FALSE)
    message("Build TranscriptDB object (txdb.sqlite) ... ", appendLF=TRUE)    
    txdb <- makeTxDbFromGFF(file=gtfFile,
        format="gtf",
        exonRankAttributeName=NA,
        dataSource="Gencode v19",
        species="Homo sapiens")


    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")


    message(paste("Prepare gene/transcript/protein id mapping information", 
            "(ids.RData) ... "), appendLF=FALSE)
    
    gff <- read.delim(gtfFile, header=FALSE, comment.char="#")
    trans_row <- gff[gff[, 'V3']=='transcript', ]
    annotation_V9 <- lapply(trans_row[, 'V9'], function(x) 
                                    strsplit(x, '; ', fixed=TRUE)[[1]][1:9])
    tmp <- lapply(annotation_V9, function(x)
                            strsplit(x, ' ', fixed=TRUE))  
    tmp2 <- lapply(tmp, function(x)
                            do.call(rbind.data.frame, x)[, 2]) 
    trans_anno <- do.call(rbind.data.frame, tmp2)

    gencode_names <- c('gene_id', 'transcript_id', 'gene_type', 
                            'gene_status', 'gene_name', 'transcript_type', 
                            'transcript_status', 'transcript_name', 
                            'level')
    colnames(trans_anno) <- gencode_names
    trans_anno[, 'level']<- gsub(';', '', trans_anno[, 'level'], fixed=TRUE)

    ids <- cbind(trans_anno, trans_anno[, 'transcript_id'])
    colnames(ids) <- c('gene_name', 'tx_name', 'gene_type', 
                            'gene_status', 'description', 'transcript_type', 
                            'transcript_status', 'transcript_name', 
                            'level', 'pro_name')

    save(ids, file=paste(annotation_path, '/ids.RData', sep=''))
    packageStartupMessage(" done")

    #tr_coding <- subset(ids,pro_name!="")
    #tr_noncoding <- subset(ids,pro_name == "")
    #txdb_coding <- makeTranscriptDbFromUCSC(genome=genome, 
    #               tablename=tablename, 
    #               transcript_ids=tr_coding[, "tx_name"] )
    #saveDb(txdb_coding, 
    #       file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))

    #txdb_noncoding <- makeTranscriptDbFromUCSC(genome=genome, 
    #        tablename=tablename, transcript_ids=tr_noncoding[, "tx_name"] )
    #saveDb(txdb_noncoding, 
    #        file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))

    message("Prepare exon annotation information (exon_anno.RData) ... ", 
            appendLF=FALSE)

    transGrange <- transcripts(txdb)
    tr <- IRanges::as.data.frame(transGrange)
    cdsByTx <- cdsBy(txdb, "tx", use.names=FALSE)
    exonByTx <- exonsBy(txdb, "tx", use.names=FALSE)
    fiveutrByTx <- fiveUTRsByTranscript(txdb, use.names=FALSE)
    threeutrByTx <- threeUTRsByTranscript(txdb, use.names=FALSE)

    #####################################################
    # get error in most recent version of IRanges package
    # previous: rownames after unlist 1.1 1.2 1.3
    # now: .1 .2 .3 ... were removed, so there are some rows with same rownames
    ######################################################
    #cdss <-  IRanges::as.data.frame(IRanges::unlist(cdsByTx))     
    #exons <- IRanges::as.data.frame(IRanges::unlist(exonByTx))
    #fiveutrs <- IRanges::as.data.frame(IRanges::unlist(fiveutrByTx))
    #threeutrs <- IRanges::as.data.frame(IRanges::unlist(threeutrByTx))

    cdss <-  IRanges::as.data.frame(cdsByTx)
    exons <- IRanges::as.data.frame(exonByTx)
    fiveutrs <- IRanges::as.data.frame(fiveutrByTx)
    threeutrs <- IRanges::as.data.frame(threeutrByTx)

    #txid <- matrix(unlist(strsplit(rownames(exons), '\\.')), ncol = 2, 
    #    byrow =TRUE)[, 1]
    #txid <- gsub('=','\\.', txid)
    exon_p <- data.frame(txid=exons[, "group_name"], chr=exons[, "seqnames"], 
                exon_s=exons[, "start"], exon_e=exons[, "end"], 
                exon_rank=exons[, "exon_rank"])
    exon2tr <-  merge(exon_p, tr, by.y="tx_id", by.x="txid")
    exon2tr <- exon2tr[, -which(names(exon2tr) %in% c("seqnames", "width"))]

    #txid <- matrix(unlist(strsplit(rownames(cdss), '\\.')), ncol = 2, 
    #       byrow =TRUE)[, 1]
    #txid <- gsub('=','\\.',txid)
    cds_p <- data.frame(txid=cdss[, "group_name"], cds_s=cdss[, "start"], 
                cds_e=cdss[, "end"], exon_rank=cdss[, "exon_rank"], 
                width=cdss[, "width"])
    ttt <- split(cds_p, cds_p$txid)

        cds_p_new_list <-lapply(ttt, function(x){
        #len <- x[,'cds_e']-x[,'cds_s']+1
        #cum <- cumsum(len)
        cum <- cumsum(x[, 'width'])
        rdis <- cbind(c(1, cum[1:length(cum)-1]+1), cum)
        colnames(rdis) <- c('cds_start', 'cds_end')
        tmp <- cbind(x, rdis)
        tmp
        })


    cds_p_new <- do.call(rbind, cds_p_new_list)
    cds_p_new <- cds_p_new[, -which(names(cds_p_new) %in% c("width"))]

    #for(i in 1:length(ttt)) {
    #    print(i)
    #    ttt1 <- ttt[[i]]
    #    len <- ttt1[,'cds_e']-ttt1[,'cds_s']+1
    #    cum <- cumsum(len)
    #    rdis <- cbind(c(1,cum[1:length(cum)-1]+1),cum)
    #    colnames(rdis) <- c('cds_start','cds_end')
    #    tmp <- cbind(ttt1,rdis)
    #    cds_p_new <- rbind(cds_p_new,tmp)
    #}

    cds2exon <- merge(exon2tr, cds_p_new, by.x=c("txid", "exon_rank"), 
                        by.y=c("txid", "exon_rank"), all.x = TRUE)
    #txid <- matrix(unlist(strsplit(rownames(fiveutrs), '\\.')), ncol = 2, 
    #               byrow=TRUE)[, 1]
    #txid <- gsub('=','\\.', txid)
    fiveutr_p <- data.frame(txid=fiveutrs[, "group_name"], 
                fiveutr_s=fiveutrs[, "start"], 
                fiveutr_e=fiveutrs[, "end"], 
                exon_rank=fiveutrs[, "exon_rank"])
    fiveutr2exon <- merge(cds2exon, fiveutr_p, by.x=c("txid", "exon_rank"), 
                by.y =c("txid", "exon_rank"), all.x = TRUE)

    #txid <- matrix(unlist(strsplit(rownames(threeutrs),'\\.')), ncol = 2, 
    #         byrow =TRUE)[, 1]
    #txid <- gsub('=','\\.', txid)
    threeutr_p <- data.frame(txid=threeutrs[, "group_name"], 
                threeutr_s=threeutrs[, "start"], 
                threeutr_e=threeutrs[, "end"], 
                exon_rank=threeutrs[, "exon_rank"])
    threeutr2exon <- merge(fiveutr2exon, threeutr_p, 
                by.x=c("txid", "exon_rank"),
                by.y=c("txid", "exon_rank"), all.x = TRUE)

    #exon <- merge(threeutr2exon,ids,by.x=c("tx_name"), by.y="tx_name", 
    #         all.x = TRUE)
    exon <- threeutr2exon[order(threeutr2exon$txid, threeutr2exon$exon_rank), ]
    
    colnames(exon) <- c("tx_id","rank", "chromosome_name", "exon_chrom_start", 
    "exon_chrom_end", "start_position", "end_position", "strand", "tx_name", 
    "cds_chr_start", "cds_chr_end", "cds_start", "cds_end", "5_utr_start", 
    "5_utr_end", "3_utr_start", "3_utr_end")

    pro_name <- ids[match(exon[, 'tx_name'], ids[, 'tx_name']), 'pro_name']
    gene_name <- ids[match(exon[, 'tx_name'], ids[, 'tx_name']), 'gene_name']
    exon <- cbind(exon, pro_name, gene_name)

    save(exon, file=paste(annotation_path, '/exon_anno.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein sequence (proseq.RData) ... ", appendLF=FALSE)
    pro_seqs <- readAAStringSet(pepfasta, format= 'fasta')
    pro_name_v <- names(pro_seqs)
    pro_name <- unlist(lapply(pro_name_v, function(x) 
                            strsplit(x, '|', fixed=TRUE)[[1]][1]))
    tx_name <- pro_name
    proteinseq <- as.data.frame(pro_seqs)
    proteinseq <- cbind(proteinseq, pro_name_v, pro_name, tx_name)
    colnames(proteinseq) <- c("peptide", "pro_name_v", "pro_name", "tx_name")
    #proteinseq <- subset(proteinseq, tx_name %in% refGene[, 'name'])
    save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein coding sequence (procodingseq.RData)... ", 
                appendLF=FALSE)
    trans_seqs <- readDNAStringSet(CDSfasta, format= 'fasta')
    cds_range <- unlist(lapply(names(trans_seqs), function(x) 
                        strsplit(strsplit(x, '|CDS:', fixed=TRUE)[[1]][2], 
                        '|', fixed=TRUE)[[1]][[1]]))
    
    tx_name <- unlist(lapply(names(trans_seqs), function(x) 
                                strsplit(x, '|', fixed=TRUE)[[1]][1]))
     
    pro_name <- tx_name
    cds_sta <- unlist(lapply(cds_range, function(x) 
                            strsplit(x, '-', fixed=TRUE)[[1]][1]))
    cds_end <- unlist(lapply(cds_range, function(x) 
                            strsplit(x, '-', fixed=TRUE)[[1]][2]))  
                            
    cds_seqs <- BStringSet(subseq(trans_seqs, start=as.integer(cds_sta), 
                                    end=as.integer(cds_end)))
    
    procodingseq <- as.data.frame(cds_seqs)
    procodingseq <- cbind(procodingseq, names(trans_seqs), pro_name, 
            tx_name, '-')
    colnames(procodingseq) <- c("coding", "tx_name_full", "pro_name", 
                "tx_name", "tx_id")
    #procodingseq <- subset(procodingseq, tx_name %in% refGene[, 'name'])
    save(procodingseq, file=paste(annotation_path, '/procodingseq.RData', 
            sep=''))
    
    packageStartupMessage(" done")
    
    
    if (!is.null(dbsnp)) {
        message("Prepare dbSNP information (dbsnpinCoding.RData) ... ", 
            appendLF=FALSE)
        if(length(dbsnps) == 1&&dbsnps == 'snp128'){
            dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp], 
                table='snp128')
        }else{
            dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp], 
                    table=paste(dbsnps[dbsnp], 'CodingDbSnp', sep=''))
        }
        snpCodingTab <- getTable(dbsnp_query)
        #snpCodingTab[, 'chrom'] <- gsub('chr', '', snpCodingTab[, 'chrom'])
        chrlist <- paste('chr', c(seq(1:22),'X','Y'), sep='')
        snpCoding <- snpCodingTab[which(snpCodingTab$chrom %in% chrlist), 
            c('chrom', 'chromStart', 'chromEnd', 'name', 'alleleCount', 
            'alleles')]
        #snpCoding <- subset(snpCodingTab, chrom %in% chrlist, 
        #        select=c(chrom:name, alleleCount, alleles))
        snpCoding <- unique(snpCoding)
        #snpCoding[, 'chrom'] <- gsub('chrM', 'MT', snpCoding[, 'chrom'])
        #
        
        #save(snpCoding,file=paste(annotation_path,'/snpcoding.RData',sep=''))
        snpCoding <- GRanges(seqnames=snpCoding[, 'chrom'], 
                    ranges=IRanges(start=snpCoding[, 'chromStart'], 
                    end=snpCoding[, 'chromEnd']), strand='*', 
                    rsid=snpCoding[, 'name'], 
                    alleleCount=snpCoding[, 'alleleCount'], 
                    alleles=snpCoding[, 'alleles'])
        
        #seqlevels(snpCoding)
        
        #if(TRUE%in% grep('chr',seqlevels(snpCoding)) > 0 ) {
        #    rchar <- sub('chr','',seqlevels(snpCoding))
        #    names(rchar) <- seqlevels(snpCoding)
        #    snpCoding <- renameSeqlevels(snpCoding, rchar) }
        #if('M'%in%seqlevels(snpCoding)){
        #snpCoding <- renameSeqlevels(snpCoding, c( M='MT'))
        #}
        #chrlist <- paste(c(seq(1:22),'X','Y'))
        transGrange_snp <- transGrange
        #transGrange_snp <- keepSeqlevels(transGrange_snp, snpCoding)
        #snpCoding <- keepSeqlevels(snpCoding, transGrange_snp)
        
        #snpCoding <- keepSeqlevels(snpCoding, transGrange)
        
        dbsnpinCoding <- subsetByOverlaps(snpCoding,transGrange_snp)
        
        save(dbsnpinCoding,file=paste(annotation_path, '/dbsnpinCoding.RData', 
            sep=''))
        packageStartupMessage(" done")
    
    }

    if (COSMIC) {
        #cosmic <- trackNames(session)[grep('cosmic',trackNames(session), 
        #        fixed=TRUE)]
        message("Prepare COSMIC information (cosmic.RData) ... ", 
            appendLF=FALSE)
        
        cosmic_query <- ucscTableQuery(session, 'cosmic', table='cosmic')
        cosmicTab <- getTable(cosmic_query)
        cosmic <- GRanges(seqnames=cosmicTab[, 'chrom'], 
                ranges=IRanges(start=cosmicTab[, 'chromStart'], 
                end=cosmicTab[, 'chromEnd']), 
                strand = '*', cosid=cosmicTab[,'name'])    

        #cosmic <- keepSeqlevels(cosmic,transGrange)
        cosmic <- subsetByOverlaps(cosmic, transGrange)

        save(cosmic,file=paste(annotation_path, '/cosmic.RData', sep=''))
        packageStartupMessage(" done")        
    }
    if(splice_matrix){
        message("Prepare exon splice information (splicemax.RData) ... ", 
                appendLF=FALSE)
        index <- which(elementLengths(exonByTx)==1)
        exonByTx_mul <- exonByTx[-index]
        exons_mul <- IRanges::as.data.frame(exonByTx_mul)
        exonslist <- split(exons_mul, exons_mul$group_name)
        #system.time( exonByTx <- exonsBy(txdb,"tx", use.names=FALSE))
        splicemax_list <- lapply(exonslist, .gen_splicmatrix)
        splicemax <- do.call(rbind, splicemax_list)
        save(splicemax, file=paste(annotation_path, '/splicemax.RData', 
            sep=''))
        packageStartupMessage(" done")        
    }

}
    
.gen_splicmatrix <- function(x, 
     ...) {           
        mystrand=x[1,'strand']
        a=x[,'exon_rank']
        b=x[,'exon_id']
        n=length(a)
        if (mystrand=='+'){
            tmp=order(a)
        
        }else{
            tmp=order(a, decreasing=TRUE)

        }
        mm=cbind(b[tmp[1:(n-1)]], b[tmp[2:n]])
        mm
    }



