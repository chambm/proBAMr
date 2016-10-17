test_PrepareAnnotationGENCODE <- function() {
    gtfFile <- system.file("extdata", "test.gtf", package="proBAMr")
    CDSfasta <- system.file("extdata", "coding_seq.fasta", package="proBAMr") 
    pepfasta <- system.file("extdata", "pro_seq.fasta", package="proBAMr") 
    annotation_path <- tempdir()
    PrepareAnnotationGENCODE(gtfFile, CDSfasta, pepfasta, 
                annotation_path)
    load(paste(annotation_path, "/ids.RData", sep=''))
    checkEqualsNumeric(length(unique(ids[, 'gene_name'])), 6)
}


