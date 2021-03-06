% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readPSMs.R
\name{readPepXmlTab}
\alias{readPepXmlTab}
\title{Read PSMs from a pepXMLTab file}
\usage{
readPepXmlTab(pepXmlTabFilepath, searchEngineScore = "mvh",
  decoyPattern = NULL)
}
\arguments{
\item{pepXmlTabFilepath}{the path to the pepXMLTab file to read}

\item{searchEngineScore}{the name of the score column to use as the search engine score, e.g. "mvh"}

\item{decoyPattern}{a regular expression pattern which matches decoy protein accessions
(but not target proteins), or NULL (the default) if all matches should be considered targets}
}
\value{
a data.frame with one row per peptide-spectrum-match and the following columns:
\enumerate{
\item Source (of the spectrum)
\item NativeID (of the spectrum)
\item ObservedNeutralMass (of the precursor)
\item TheoreticalNeutralMass (of the peptide)
\item MassError (monoisotopic)
\item Charge (of the precursor, as identified)
\item ScanTimeInSeconds (relative to start of the run)
\item PeptideSequence (no modifications)
\item PeptideModifications (e.g. '-17.026549@(,57.021464@3,15.994915@7', but the exact format depends on the PSM source)
\item ProteinAccessions (comma separated)
\item ProteinCount (number of associated proteins)
\item IsDecoy (TRUE or FALSE)
\item NumSpecificTermini (2, 1, or 0)
\item MissedCleavages
\item Rank (1 being the best hit)
\item FdrScore (or 1 - confidence)
\item SearchEngineScore (specific to a given search engine, e.g. 'myrimatch:mvh')
}
}
\description{
Read PSMs from a pepXMLTab file
}
\examples{
pepXmlTabFilepath = system.file("extdata/psm/passedPSM.tab", package="proBAMr")
result = readPepXmlTab(pepXmlTabFilepath)
head(result)
}
\seealso{
Other functions to read PSMs: \code{\link{readIdpDB}},
  \code{\link{readPeptideShakerPsmReport}}
}
