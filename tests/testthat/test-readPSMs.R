library(testthat)
library(proBAMr)

context("readPSMs")

on.update.view = customProDB:::on.update.view
on.fail.diff = customProDB:::on.fail.diff

test_that("Reading PSMs from an idpDB", {
  result = readIdpDB(system.file("extdata/psm/test.idpDB", package="proBAMr"), "myrimatch:mvh")
  expect_is(result$PeptideSequence, "character")
  expect_is(result$PeptideModifications, "character")
  customProDB:::expect_equal_to_reference(result, "readIdpDB.rds", on.update=on.update.view, on.fail=on.fail.diff)
})

test_that("Reading PSMs from a PeptideShaker report", {
  result = readPeptideShakerPsmReport(system.file("extdata/psm/test.psm-report", package="proBAMr"))
  expect_is(result$PeptideSequence, "character")
  expect_is(result$PeptideModifications, "character")
  customProDB:::expect_equal_to_reference(result, "readPeptideShakerPsmReport.rds", on.update=on.update.view, on.fail=on.fail.diff)
})

test_that("Reading PSMs from a pepXMLTab file", {
  result = readPepXmlTab(system.file("extdata/psm/passedPSM.tab", package="proBAMr"), searchEngineScore="mvh", decoyPattern=NULL)
  expect_is(result$PeptideSequence, "character")
  expect_is(result$PeptideModifications, "character")
  customProDB:::expect_equal_to_reference(result, "readPepXmlTab.rds", on.update=on.update.view, on.fail=on.fail.diff)
})
