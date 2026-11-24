\name{Spider}
\alias{Spider}
\docType{data}
\title{
  Spider species abundance data
}
\description{
  This dataset contains spider species abundance data collected across
  multiple sites under two treatments (Control and Enhanced).
  The data are organized by site, treatment, and patch for the
  meta-analysis of taxonomic, phylogenetic, and functional diversity.
}
\usage{data(Spider)}
\format{
  \code{Spider} is a data.frame with the following structure:

  \itemize{
    \item \code{Site}: Identifier of the study site.
    \item \code{Treatment}: Treatment category (Control or Enhanced).
    \item \code{Patch}: Identifier of the sampling patch within each
    site and treatment combination.
    \item Species columns: All remaining columns correspond to spider
    species, with non-negative integer values representing species
    abundances within each patch.
  }
}
\details{
  Each row represents a sampling patch within a particular site
  and treatment combination.

  The \code{Patch} column is metadata and is excluded from the
  species abundance matrix during analysis.
}
\keyword{datasets}
