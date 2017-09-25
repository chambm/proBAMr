library(testthat)
library(customProDB)
library(proBAMr)

context("PrepareAnnotationGENCODE")
options(stringsAsFactors = FALSE)

load_annotations = function(annotation_path, envir, dbsnp=FALSE, cosmic=FALSE) {
  
  stopifnot(file.exists(paste0(annotation_path, "/exon_anno.RData")))
  stopifnot(file.exists(paste0(annotation_path, "/ids.RData")))
  stopifnot(file.exists(paste0(annotation_path, "/procodingseq.RData")))
  stopifnot(file.exists(paste0(annotation_path, "/proseq.RData")))
  
  load(paste0(annotation_path, "/exon_anno.RData"), envir=envir)
  load(paste0(annotation_path, "/ids.RData"), envir=envir)
  load(paste0(annotation_path, "/procodingseq.RData"), envir=envir)
  load(paste0(annotation_path, "/proseq.RData"), envir=envir)
  
  expect_is(envir$exon, "data.frame")
  expect_is(envir$ids, "data.frame")
  expect_is(envir$proteinseq, "data.frame")
  expect_is(envir$procodingseq, "data.frame")
  
  if (dbsnp) {
    load(paste0(annotation_path, "/dbsnpinCoding.RData"), envir=envir)
    expect_is(envir$dbsnpinCoding, "GRanges")
  }
  
  if (cosmic) {
    load(paste0(annotation_path, "/cosmic.RData"), envir=envir)
    expect_is(envir$cosmic, "GRanges")
  }
}

on.update.view=customProDB:::on.update.view
on.fail.diff=customProDB:::on.fail.diff

test_that("Preparing annotation from GENCODE", {
    gtfFile <- system.file("extdata", "test.gtf", package="proBAMr")
    CDSfasta <- system.file("extdata", "coding_seq.fasta", package="proBAMr") 
    pepfasta <- system.file("extdata", "pro_seq.fasta", package="proBAMr") 
    annotation_path <- tempdir()

    if (BiocInstaller::biocVersion() < numeric_version("3.5")) {
      expect_warning(PrepareAnnotationGENCODE(gtfFile, CDSfasta, pepfasta, annotation_path))
    } else {
      PrepareAnnotationGENCODE(gtfFile, CDSfasta, pepfasta, annotation_path)
    }
      
    env = new.env()
    load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
    customProDB:::expect_equal_to_reference(env$exon, "exon.rds", on.update=on.update.view, on.fail=on.fail.diff)
    customProDB:::expect_equal_to_reference(env$ids, "ids.rds", on.update=on.update.view, on.fail=on.fail.diff)
    customProDB:::expect_equal_to_reference(env$proteinseq, "proteinseq.rds", on.update=on.update.view, on.fail=on.fail.diff)
    customProDB:::expect_equal_to_reference(env$procodingseq, "procodingseq.rds", on.update=on.update.view, on.fail=on.fail.diff)
})
