#' @title proBAMr
#' @description Maps peptide-spectrum-matches (PSMs) to a genome. The package builds SAM file
#' from shotgun proteomics data. The package also provides function
#' to prepare annotation from GTF file.
#' @importFrom utils read.delim read.table setTxtProgressBar txtProgressBar write.table 
#' @importFrom AnnotationDbi saveDb loadDb
#' @importFrom rtracklayer browserSession ucscTableQuery tableNames getTable trackNames ucscSchema
#' @importFrom data.table data.table rbindlist setkey setDT := .SD
#' @import Biostrings IRanges GenomicRanges GenomicFeatures
"_PACKAGE"

# suppress NOTE about undefined global variables due to data.table's non-standard evaluation
utils::globalVariables(c("PeptideSequence", "SpectrumUniqueId","cds_chr_end", "cds_chr_start", "cds_end",
                         "cds_start", "chromosome_name", "coding", "dna_complement", "peptide",
                         "peptide_count", "pro_name", "psm_count", "tx_id", "tx_name", "."))