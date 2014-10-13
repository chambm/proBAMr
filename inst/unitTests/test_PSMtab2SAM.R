test_PSMtab2SAM <- function() 
{
    load(system.file("extdata/GENCODE", "exon_anno.RData", package="proBAMr"))
    load(system.file("extdata/GENCODE", "proseq.RData", package="proBAMr"))
    load(system.file("extdata/GENCODE", "procodingseq.RData", 
        package="proBAMr"))
    options(stringsAsFactors=FALSE)
    passedPSM <- read.table(system.file("extdata", "passedPSM.tab", 
        package="proBAMr"), sep='\t', header=TRUE)
    load(system.file("extdata/res", "SAM.RData", package="proBAMr"))
    
    checkIdentical(PSMtab2SAM(passedPSM, XScolumn='mvh', exon, proteinseq, 
        procodingseq), SAM)    
    checkEquals(dim(SAM)[2], 21)
    checkTrue(!is.na(SAM[1, 1]))
    checkEqualsNumeric(as.integer(SAM[1, 5]), 255)    
}

