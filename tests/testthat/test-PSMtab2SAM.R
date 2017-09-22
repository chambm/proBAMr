library(testthat)
library(customProDB)
library(proBAMr)

context("PSMtab2SAM")
options(stringsAsFactors = FALSE)


on.update.view = customProDB:::on.update.view
on.fail.diff = customProDB:::on.fail.diff


test_that("PSM2SAM works with IDPicker for hg19 annotation and variant parameters", {
  
  load(system.file('extdata/hg19/ids.RData', package="proBAMr"))
  load(system.file('extdata/hg19/procodingseq.RData', package="proBAMr"))
  load(system.file('extdata/hg19/proseq.RData', package="proBAMr"))
  load(system.file('extdata/hg19/exon_anno.RData', package="proBAMr"))
  load(system.file('extdata/hg19/dbsnpinCoding.RData', package="proBAMr"))
  load(system.file('extdata/hg19/cosmic.RData', package="proBAMr"))
  
  vcfFilepath = system.file("extdata/hg19/test.vcf", package="proBAMr")
  psmFilepath = system.file("extdata/psm/test.idpDB", package="proBAMr")
  passedPSM = readIdpDB(psmFilepath, "MyriMatch:MVH")
  
  variantAnnotation <- getVariantAnnotation(vcfFilepath,
                                            ids, exon,
                                            proteinseq, procodingseq,
                                            dbsnpinCoding, cosmic)
  snvprocoding = variantAnnotation$snvprocoding
  snvproseq = variantAnnotation$snvproseq

  samFilepath = paste0(tools::file_path_sans_ext(psmFilepath), "-hg19-idpicker-with-vcf.sam")
  file.create(samFilepath)
  
  PSMtab2SAM(passedPSM, exon,
             proteinseq, procodingseq,
             snvproseq, snvprocoding,
             outfile = samFilepath,
             show_progress = FALSE)

  # compare the resulting SAM to the one in test directory
  customProDB:::expect_equal_to_reference(read.delim(samFilepath, header=FALSE),
                                          "test-hg19-idpicker-with-vcf.sam.rds",
                                          on.update = on.update.view, on.fail = on.fail.diff)
})


test_that("PSM2SAM works with PeptideShaker for mm10 annotation", {
  
  load(system.file('extdata/mm10/ids.RData', package="proBAMr"))
  load(system.file('extdata/mm10/procodingseq.RData', package="proBAMr"))
  load(system.file('extdata/mm10/proseq.RData', package="proBAMr"))
  load(system.file('extdata/mm10/exon_anno.RData', package="proBAMr"))

  psmFilepath = system.file("extdata/psm/test.psm-report", package="proBAMr")
  passedPSM = readPeptideShakerPsmReport(psmFilepath)
  samFilepath = file.path(tempdir(), basename(sub(".psm-report", "-mm10-peptideshaker-no-vcf.sam", psmFilepath)))
  file.create(samFilepath)
  
  # first run, no VCF
  PSMtab2SAM(passedPSM, exon,
             proteinseq, procodingseq,
             outfile = samFilepath,
             show_progress = FALSE)
  
  # compare the resulting SAM to the one in test directory
  customProDB:::expect_equal_to_reference(read.delim(samFilepath, header=FALSE),
                                          "test-mm10-peptideshaker-no-vcf.sam.rds",
                                          on.update = on.update.view, on.fail = on.fail.diff)
  
  vcfFilepaths = system.file(c("extdata/mm10/test-sample1.vcf.gz",
                               "extdata/mm10/test-sample2.vcf.gz"),
                             package="proBAMr")
  variantAnnotation <- getVariantAnnotation(vcfFilepaths,
                                            ids, exon,
                                            proteinseq, procodingseq,
                                            dbsnpinCoding = NULL,
                                            cosmic = NULL)
                 
  varprocoding = unique(rbind(variantAnnotation$snvprocoding, variantAnnotation$indelprocoding))
  varproseq = unique(rbind(variantAnnotation$snvproseq, variantAnnotation$indelproseq))

  samFilepath = file.path(tempdir(), basename(sub(".psm-report", "-mm10-peptideshaker-with-vcf.sam", psmFilepath)))
  file.create(samFilepath)

  # second run, with VCF
  PSMtab2SAM(passedPSM, exon,
             proteinseq, procodingseq,
             varprocoding, varproseq,
             outfile = samFilepath,
             show_progress = FALSE)
  
  # compare the resulting SAM to the one in test directory
  customProDB:::expect_equal_to_reference(read.delim(samFilepath, header=FALSE),
                                          "test-mm10-peptideshaker-with-vcf.sam.rds",
                                          on.update = on.update.view, on.fail = on.fail.diff)
})


test_that("PSM2SAM works with pepXMLTab for GENCODE annotation and no variant parameters ", {
  
  load(system.file('extdata/GENCODE/ids.RData', package="proBAMr"))
  load(system.file('extdata/GENCODE/procodingseq.RData', package="proBAMr"))
  load(system.file('extdata/GENCODE/proseq.RData', package="proBAMr"))
  load(system.file('extdata/GENCODE/exon_anno.RData', package="proBAMr"))
  
  psmFilepath = system.file("extdata/psm/passedPSM.tab", package="proBAMr")
  passedPSM = readPepXmlTab(psmFilepath, decoyPattern=NULL, searchEngineScore="mvh")
  samFilepath = paste0(tools::file_path_sans_ext(psmFilepath), "-gencode-pepxmltab-no-vcf.sam")
  file.create(samFilepath)
  
  load(system.file("extdata/res/SAM.RData", package="proBAMr"))
  
  PSMtab2SAM(passedPSM, exon,
             proteinseq, procodingseq,
             outfile = samFilepath,
             show_progress = FALSE)
  
  # compare the resulting SAM to the one in test directory
  customProDB:::expect_equal_to_reference(read.delim(samFilepath, header=FALSE),
                                          "test-gencode-pepxmltab-no-vcf.sam.rds",
                                          on.update = on.update.view, on.fail = on.fail.diff)
})
