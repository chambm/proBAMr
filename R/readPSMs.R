
#' Read PSMs from an IDPicker DB
#'
#' @param idpDbFilepath the path to the idpDB file to read
#' @param scoreName the name of the PSM score (e.g. 'myrimatch:mvh') to extract
#'
#' @return a data.frame with one row per peptide-spectrum-match and the following columns:
#' \enumerate{
#' \item Source (of the spectrum)
#' \item NativeID (of the spectrum)
#' \item ObservedNeutralMass (of the precursor)
#' \item TheoreticalNeutralMass (of the peptide)
#' \item MassError (monoisotopic)
#' \item Charge (of the precursor, as identified)
#' \item ScanTimeInSeconds (relative to start of the run)
#' \item PeptideSequence (no modifications)
#' \item PeptideModifications (e.g. '-17.026549@(,57.021464@3,15.994915@7', but the exact format depends on the PSM source)
#' \item ProteinAccessions (comma separated)
#' \item ProteinCount (number of associated proteins)
#' \item IsDecoy (TRUE or FALSE)
#' \item NumSpecificTermini (2, 1, or 0)
#' \item MissedCleavages
#' \item Rank (1 being the best hit)
#' \item FdrScore (or 1 - confidence)
#' \item SearchEngineScore (specific to a given search engine, e.g. 'myrimatch:mvh')
#' }
#' 
#' @family functions to read PSMs
#' @export
#'
#' @examples
#' idpDbFilepath = system.file("extdata/psm/test.idpDB", package="proBAMr")
#' result = readIdpDB(idpDbFilepath, "myrimatch:mvh")
#' head(result)

readIdpDB <- function(idpDbFilepath, scoreName)
{
  old <- options(stringsAsFactors = FALSE)
  on.exit(options(old), add = TRUE)
  
  stopifnot(file.exists(idpDbFilepath))
  
  drv <- RSQLite::SQLite()
  con <- RSQLite::dbConnect(drv, idpDbFilepath)
  
  res <- RSQLite::dbSendQuery(con, paste0("SELECT Id FROM PeptideSpectrumMatchScoreName WHERE lower(Name)='", stringr::str_to_lower(scoreName), "'"))
  scoreId = RSQLite::fetch(res, n=1)["Id"]
  RSQLite::dbClearResult(res)
  
  sql = paste0("SELECT ss.Name AS Source",
               ", s.NativeID",
               ", psm.ObservedNeutralMass",
               ", psm.ObservedNeutralMass - (psm.MonoisotopicMassError - ROUND(psm.MonoisotopicMassError) * 1.0026) AS TheoreticalNeutralMass",
               ", psm.MonoisotopicMassError - ROUND(psm.MonoisotopicMassError) * 1.0026 AS MassError",
               ", psm.Charge",
               ", s.ScanTimeInSeconds",
               ", IFNULL(SUBSTR(pd.Sequence, pi.Offset+1, pi.Length), pep.DecoySequence) AS PeptideSequence",
               ", GROUP_CONCAT(DISTINCT mod.MonoMassDelta || '@' || pm.Offset) AS PeptideModifications",
               ", GROUP_CONCAT(DISTINCT pro.Accession || IFNULL('@' || (pi.Offset+1), '')) AS ProteinAccessions",
               ", COUNT(DISTINCT pro.Id) AS ProteinCount",
               ", IsDecoy AS IsDecoy",
               ", CTerminusIsSpecific+NTerminusIsSpecific AS NumSpecificTermini",
               ", MissedCleavages",
               ", psm.Rank",
               ", psm.QValue AS FdrScore",
               ", psmScore.Value AS SearchEngineScore",
               " FROM PeptideSpectrumMatch psm ",
               "JOIN Spectrum s ON psm.Spectrum=s.Id ",
               "JOIN SpectrumSource ss ON s.Source=ss.Id ",
               "JOIN PeptideInstance pi ON psm.Peptide=pi.Peptide ",
               "JOIN Protein pro ON pi.Protein=pro.Id ",
               "JOIN Peptide pep ON pi.Peptide=pep.Id ",
               "JOIN PeptideSpectrumMatchScore psmScore ON psmScore.PsmId=psm.Id AND ScoreNameId=", scoreId,
               " LEFT JOIN ProteinData pd ON pro.Id=pd.Id ",
               "LEFT JOIN PeptideModification pm ON psm.Id=pm.PeptideSpectrumMatch ",
               "LEFT JOIN Modification mod ON pm.Modification=mod.Id ",
               "GROUP BY psm.Id"
  )
  #print(sql)
  
  res <- RSQLite::dbSendQuery(con, sql)
  passedPSM <- RSQLite::fetch(res, n=-1)
  RSQLite::dbClearResult(res)
  RSQLite::dbDisconnect(con)

  passedPSM
}


#' Read PSMs from a PeptideShaker PSM report
#'
#' @param psmReportFilepath the path to the PSM report to read
#'
#' @inherit readIdpDB return
#' @family functions to read PSMs
#' @export
#'
#' @examples
#' psmReportFilepath = system.file("extdata/psm/test.psm-report", package="proBAMr")
#' result = readPeptideShakerPsmReport(psmReportFilepath)
#' head(result)
readPeptideShakerPsmReport <- function(psmReportFilepath)
{
  old <- options(stringsAsFactors = FALSE)
  on.exit(options(old), add = TRUE)
  
  stopifnot(file.exists(psmReportFilepath))

  # Source columns:
  # Index
  # Protein(s)
  # Sequence
  # AAs Before
  # AAs After
  # Position
  # Modified Sequence
  # Variable Modifications
  # Fixed Modifications
  # Spectrum File
  # Spectrum Title
  # Spectrum Scan Number
  # RT
  # m/z
  # Measured Charge
  # Identification Charge
  # Theoretical Mass
  # Isotope Number
  # Precursor m/z Error [ppm]
  # Localization Confidence
  # Probabilistic PTM score
  # D-score
  # Confidence [%]
  # Validation

  # Result columns:
  # Source NativeID ObservedNeutralMass TheoreticalNeutralMass MassError Charge ScanTimeInSeconds
  # PeptideSequence PeptideModifications ProteinAccessions ProteinCount IsDecoy
  # NumSpecificTermini MissedCleavages Rank FdrScore SearchEngineScore 
  result = read.delim(psmReportFilepath, stringsAsFactors = FALSE)
  
  if (sub(".*File:\"([^\"]+)\".*", "\\1", result$Spectrum.Title[1]) == result$Spectrum.Title[1]) {
    spectrumSources = tools::file_path_sans_ext(result$Spectrum.File)
  } else {
    spectrumSources = tools::file_path_sans_ext(sub(".*File:\"([^\"]+)\".*", "\\1", result$Spectrum.Title))
  }
  
  if (sub(".*NativeID:\"([^\"]+)\"", "\\1", result$Spectrum.Title[1]) == result$Spectrum.Title[1]) {
    spectrumNativeIDs = sub(".*?\\.(\\d+)\\.\\d+.*", "\\1", result$Spectrum.Title)
  } else {
    spectrumNativeIDs = sub(".*NativeID:\"([^\"]+)\"", "\\1", result$Spectrum.Title)
  }
  
  identifiedCharges = as.numeric(sub("\\+", "", result$Identification.Charge))
  
  peptideMods = mapply(function(x, y)
  {
    if (nzchar(x) && nzchar(y))
      paste(x, y, sep=", ")
    else if (nzchar(x))
      x
    else y
  }, result$Variable.Modifications, result$Fixed.Modifications, USE.NAMES = FALSE)
  
  data.frame(Source = spectrumSources,
             NativeID = spectrumNativeIDs,
             ObservedNeutralMass = (result$m.z * identifiedCharges) - (identifiedCharges * 1.007276),
             TheoreticalNeutralMass = result$Theoretical.Mass,
             MassError = result$Theoretical.Mass * result$Precursor.m.z.Error..ppm. / 1e6,
             Charge = identifiedCharges,
             ScanTimeInSeconds = result$RT,
             PeptideSequence = result$Sequence,
             PeptideModifications = peptideMods,
             ProteinAccessions = result$Protein.s.,
             ProteinCount = stringr::str_count(result$Protein.s., ", ")+1,
             IsDecoy = FALSE,
             NumSpecificTermini = NA,
             MissedCleavages = NA,
             Rank = 1,
             FdrScore = (100 - result$Confidence....)/100,
             SearchEngineScore = NA)
}

#' Read PSMs from a pepXMLTab file
#'
#' @param pepXmlTabFilepath the path to the pepXMLTab file to read
#' @param searchEngineScore the name of the score column to use as the search engine score, e.g. "mvh"
#' @param decoyPattern a regular expression pattern which matches decoy protein accessions
#' (but not target proteins), or NULL (the default) if all matches should be considered targets
#'
#' @inherit readIdpDB return
#' @family functions to read PSMs
#' @export
#'
#' @examples
#' pepXmlTabFilepath = system.file("extdata/psm/passedPSM.tab", package="proBAMr")
#' result = readPepXmlTab(pepXmlTabFilepath)
#' head(result)

readPepXmlTab <- function(pepXmlTabFilepath, searchEngineScore = "mvh", decoyPattern = NULL)
{
  old <- options(stringsAsFactors = FALSE)
  on.exit(options(old), add = TRUE)
  
  stopifnot(file.exists(pepXmlTabFilepath))
  passedPSM = read.delim(pepXmlTabFilepath, row.names=NULL)
  stopifnot(searchEngineScore %in% colnames(passedPSM))
  
  spectrumSource = sub("(.*?)\\.\\d+\\.\\d+.\\d+$", "\\1", passedPSM$spectrum)
  
  if (!is.null(decoyPattern))
    isDecoy = grepl(passedPSM$protein, decoyPattern)
  else
    isDecoy = FALSE
  
  data.frame(Source = spectrumSource,
             NativeID = passedPSM$spectrumNativeID,
             ObservedNeutralMass = passedPSM$precursor_neutral_mass,
             TheoreticalNeutralMass = passedPSM$calc_neutral_pep_mass,
             MassError = passedPSM$massdiff,
             Charge = passedPSM$assumed_charge,
             ScanTimeInSeconds = passedPSM$retention_time_sec,
             PeptideSequence = passedPSM$peptide,
             PeptideModifications = passedPSM$modification,
             ProteinAccessions = passedPSM$protein,
             ProteinCount = passedPSM$num_tot_proteins,
             IsDecoy = isDecoy,
             NumSpecificTermini = passedPSM$num_tol_term,
             MissedCleavages = passedPSM$num_missed_cleavages,
             Rank = passedPSM$hit_rank,
             FdrScore = NA,
             SearchEngineScore = passedPSM[, searchEngineScore])
}
