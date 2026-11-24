\name{Hoverfly}
\alias{Hoverfly}
\docType{data}
\title{
  Hoverfly species incidence data
}
\description{
  This dataset contains hoverfly species incidence data collected across
  multiple sites under two treatments (Control and Enhanced).
  The data are organized by site, treatment, patch, and sampling unit
  for meta-analysis of taxonomic, phylogenetic, and functional diversity.
}
\usage{data(Hoverfly)}
\format{
  \code{Hoverfly} is a data.frame with the following structure:

  \itemize{
    \item \code{Site}: Identifier of the study site.
    \item \code{Treatment}: Treatment category (Control or Enhanced).
    \item \code{Patch}: Identifier of the sampling patch (assemblage)
    within each site and treatment.
    \item \code{Sampling_unit}: Identifier of the sampling unit within
    each patch.
    \item Species columns: All remaining columns correspond to hoverfly
    species, with incidence recorded as 0 (non-detection) or
    1 (detection) for each sampling unit.
  }
}
\details{
  Each row represents one sampling unit within a particular patch,
  site, and treatment combination.

  Multiple sampling units are recorded within each patch.
  For incidence-based diversity analysis, the sampling units within
  a patch provide the species detection history, while patches
  represent assemblages used in the meta-analysis.

  The \code{Sampling_unit} column is metadata and should not be
  treated as a species column.
}
\keyword{datasets}

