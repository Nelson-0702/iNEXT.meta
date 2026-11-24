#'
#' Data information and alpha-coverage summary for reference samples
#' \code{DataInfobeta3Dmeta()} provides basic data information for each study/site and treatment combination, including the pooled assemblage (gamma reference sample) and the joint assemblage (alpha reference sample). The function also summarizes the estimated sample coverage of the alpha reference samples at twice the reference sample size across all study/site and treatment combinations.
#' 
#' @param data For \code{datatype = "abundance"}, a data.frame with site and treatment as the first two columns, optionally followed by a \code{Patch} column, and then species abundance columns. Each row represents an assemblage. For \code{datatype = "incidence_raw"}, a data.frame with site,
#' treatment, and patch as the first three columns, optionally followed by a \code{Sampling_unit} column, and then species incidence columns (0 or 1). Each row represents a sampling unit within a patch.
#' For incidence data, patches with no observed species are removed before diversity estimation; an error is raised if any site-treatment group is lost.
#' @param diversity Selection of the diversity type: \code{diversity = "TD"} for taxonomic diversity, \code{diversity = "PD"} for phylogenetic diversity, and \code{diversity = "FD"} for functional diversity.
#' @param datatype Type of input data: individual-based abundance data (\code{datatype = "abundance"}) or sampling-unit-based incidence data (\code{datatype = "incidence_raw"}). For incidence raw data, species incidence is recorded as 0 (non-detection) or 1 (detection).
#' @param PDtree Required when \code{diversity = "PD"}. A phylogenetic tree (an R \code{phylo} object) spanning all observed species across the input data. The tree can be imported from a file in Newick format.
#' @param PDreftime Argument only for \code{diversity = "PD"} and \code{PDtype = "meanPD"}. A numerical value specifying a fixed reference time for phylogenetic diversity. Default is \code{PDreftime = NULL}, in which case the root age of the phylogenetic tree spanning the pooled observed species is used.
#' @param PDtype Argument only for \code{diversity = "PD"}. Specifies the phylogenetic diversity measure: \code{PDtype = "meanPD"} for the effective number of equally divergent lineages at a specified reference time, or \code{PDtype = "AUC"} for an overall measure obtained by integrating and averaging mean phylogenetic diversity over all reference times from zero to the root age of the phylogenetic tree. Default is \code{PDtype = "AUC"}.
#' @param FDdistM Required when \code{diversity = "FD"}. A species pairwise distance matrix covering all species observed across the input data.
#' @param FDtype Argument only for \code{diversity = "FD"}. Specifies the functional diversity measure: \code{FDtype = "tau_value"} for FD at a specified threshold of functional distinctiveness, or \code{FDtype = "AUC"} for an overall FD measure obtained by integrating the tau-profile over threshold values from 0 to 1. Default is \code{FDtype = "AUC"}.
#' @param FDtau Argument only for \code{diversity = "FD"} and \code{FDtype = "tau_value"}. A numeric value between 0 and 1 specifying the threshold of functional distinctiveness used to compute FD. If \code{FDtau = NULL} (default), the threshold is set to the mean pairwise functional distance between two individuals randomly selected from the pooled data across the input data (i.e., Rao\'s quadratic entropy).
#'
#' @importFrom dplyr group_by summarise filter select across pick everything ungroup reframe pull
#' @importFrom tidyr unnest
#' @importFrom iNEXT.beta3D DataInfobeta3D
#' @importFrom magrittr %>%
#' @importFrom stringr str_to_title
#' @importFrom utils globalVariables
#'
#'
#'
#'
#' @return A list containing two components:
#' \itemize{
#'   \item \code{DataInfo}: a data.frame containing basic data information for each study/site and treatment combination. Both the pooled assemblage (gamma reference sample) and joint assemblage (alpha reference sample) are included.
#'
#'   \item \code{Summary_Coverage}: a data.frame summarizing the sample coverage of the alpha reference samples at twice the reference sample size. The five reported values are the minimum, 25th percentile, median, 75th percentile, and maximum coverage, denoted by \code{C_Min}, \code{C_25quantile}, \code{C_50quantile}, \code{C_75quantile}, and \code{C_Max}, respectively.
#' }
#'
#' For abundance data, alpha coverage in \code{Summary_Coverage} is based on \code{SC(2n)} of the joint assemblage. For incidence data, alpha coverage is based on \code{SC(2T)} of the joint assemblage.
#'
#' For abundance data, basic information in \code{DataInfo} shared by TD, PD and FD includes study/site name (\code{Site}), treatment name
#' (\code{Treatment}), pooled/joint assemblage (\code{Assemblage}), sample size (\code{n}), observed species richness (\code{S.obs}), sample coverage estimate
#' of the reference sample (\code{SC(n)}), and sample coverage estimate for twice the reference sample size (\code{SC(2n)}). Other additional information is given below.
#'
#' \itemize{
#'   \item TD: the first five species abundance frequency counts in the reference sample (\code{f1-f5}).
#'   \item PD (\code{PDtype = "meanPD"}): the observed phylogenetic diversity (\code{PD.obs}), representing the total observed branch length; the observed mean phylogenetic diversity (\code{meanPD.obs}), calculated as \code{PD.obs / Reftime}; the numbers of singletons (\code{f1*}) and doubletons (\code{f2*}) in the node/branch abundance set; their corresponding total branch lengths (\code{g1} and \code{g2}); and the reference time (\code{Reftime}). If \code{PDreftime = NULL}, the root age of the phylogenetic tree spanning the pooled observed species is used.
#'   \item PD (\code{PDtype = "AUC"}): the integration end time (\code{Tmax}) for AUC-based phylogenetic diversity, corresponding to the root age of the phylogenetic tree spanning the pooled observed species. The integration starts at zero.
#'   \item FD (\code{FDtype = "AUC"}): the minimum distance (\code{dmin}) and maximum distance (\code{dmax}) among all non-diagonal elements in the distance matrix, and the mean distance between any two individuals randomly selected from the dataset (\code{dmean}).
#'   \item FD (\code{FDtype = "tau_value"}): the number of singletons (\code{a1*}) and doubletons (\code{a2*}) among the functionally indistinct set at the specified threshold level \code{Tau}, as well as the total contribution of singletons (\code{h1}) and doubletons (\code{h2}) at the specified threshold level \code{Tau}.
#' }
#'
#'
#' For incidence data, the basic information for TD in \code{DataInfo} includes
#' study/site name (\code{Site}), treatment name (\code{Treatment}),
#' pooled/joint assemblage (\code{Assemblage}), number of sampling units
#' (\code{T}), total number of incidences (\code{U}), observed species richness
#' (\code{S.obs}), sample coverage estimate of the reference sample
#' (\code{SC(T)}), sample coverage estimate for twice the reference sample size
#' (\code{SC(2T)}), and the first five species incidence frequency counts
#' (\code{Q1-Q5}) in the reference sample. For PD and FD, the output is
#' similar to that for abundance data.
#'
#'
#'
#' @examples
#' ## (Data Information) Taxonomic diversity for abundance data
#' # See Examples 1 and 2 in the iNEXT.meta vignette for the Spider TD analyses.
#' 
#' data("Spider")
#' info_TD_abu <- DataInfobeta3Dmeta(data = Spider,
#'                                          diversity = 'TD',
#'                                          datatype = 'abundance')
#' info_TD_abu$DataInfo
#' info_TD_abu$Summary_Coverage
#'
#'
#' ## (Data Information) Taxonomic diversity for incidence data
#' # See Examples 3 and 4 in the iNEXT.meta vignette for the Hoverfly TD analyses.
#' 
#' data("Hoverfly")
#' info_TD_inc <- DataInfobeta3Dmeta(data = Hoverfly,
#'                                          diversity = 'TD',
#'                                          datatype = 'incidence_raw')
#' info_TD_inc$DataInfo
#' info_TD_inc$Summary_Coverage
#'
#' \donttest{
#' ## (Data Information) Mean phylogenetic diversity for abundance data
#' # See the PD section of the iNEXT.meta vignette for an overview of meanPD and AUC.
#' 
#' data("Spider")
#' data("Spider_tree")
#' info_meanPD_abu <- DataInfobeta3Dmeta(data = Spider,
#'                                   diversity = 'PD', datatype = 'abundance',
#'                                   PDtree = Spider_tree, PDtype = "meanPD")
#' info_meanPD_abu$DataInfo
#' info_meanPD_abu$Summary_Coverage
#'
#'
#' ## (Data Information) Phylogenetic diversity (AUC) for abundance data
#' # See Example 5 in the iNEXT.meta vignette for the Spider PD AUC analysis.
#' 
#' data("Spider")
#' data("Spider_tree")
#' info_PDAUC_abu <- DataInfobeta3Dmeta(data = Spider,
#'                                   diversity = 'PD', datatype = 'abundance',
#'                                   PDtree = Spider_tree, PDtype = "AUC")
#' info_PDAUC_abu$DataInfo
#' info_PDAUC_abu$Summary_Coverage
#'
#' ## (Data Information) Mean phylogenetic diversity for incidence data
#' # See the PD section of the iNEXT.meta vignette for the available PD measures.
#' 
#' data("Hoverfly")
#' data("Hoverfly_tree")
#' info_meanPD_inc <- DataInfobeta3Dmeta(data = Hoverfly,
#'                                   diversity = 'PD', datatype = 'incidence_raw',
#'                                   PDtree = Hoverfly_tree, PDtype = "meanPD")
#' info_meanPD_inc$DataInfo
#' info_meanPD_inc$Summary_Coverage
#'
#' ## (Data Information) Phylogenetic diversity (AUC) for incidence data
#' # See Example 5 in the iNEXT.meta vignette for a PD AUC example using Spider.
#' 
#' data("Hoverfly")
#' data("Hoverfly_tree")
#' info_PDAUC_inc <- DataInfobeta3Dmeta(data = Hoverfly,
#'                                   diversity = 'PD', datatype = 'incidence_raw',
#'                                   PDtree = Hoverfly_tree, PDtype = "AUC")
#' info_PDAUC_inc$DataInfo
#' info_PDAUC_inc$Summary_Coverage
#'
#' ## (Data Information) Functional diversity for abundance data under a specified threshold level
#' # See the FD section of the iNEXT.meta vignette for tau_value and AUC definitions.
#' 
#' data("Spider")
#' data("Spider_distM")
#' info_FDtau_abu <- DataInfobeta3Dmeta(data = Spider,
#'                                      diversity = 'FD', datatype = 'abundance',
#'                                      FDdistM = Spider_distM, FDtype = "tau_value", FDtau = NULL)
#' info_FDtau_abu$DataInfo
#' info_FDtau_abu$Summary_Coverage
#'
#'
#' ## (Data Information) Functional diversity for abundance data when all
#' ## threshold levels from 0 to 1 are considered
#' # See Example 6 in the iNEXT.meta vignette for the Spider FD AUC analysis.
#' 
#' data("Spider")
#' data("Spider_distM")
#' info_FDAUC_abu <- DataInfobeta3Dmeta(data = Spider,
#'                                      diversity = 'FD', datatype = 'abundance',
#'                                      FDdistM = Spider_distM, FDtype = "AUC")
#' info_FDAUC_abu$DataInfo
#' info_FDAUC_abu$Summary_Coverage
#'
#'
#' ## (Data Information) Functional diversity for incidence data under a specified threshold level
#' # See the FD section of the iNEXT.meta vignette for tau_value and AUC definitions.
#' 
#' data("Hoverfly")
#' data("Hoverfly_distM")
#' info_FDtau_inc <- DataInfobeta3Dmeta(data = Hoverfly,
#'                                      diversity = 'FD', datatype = 'incidence_raw',
#'                                      FDdistM = Hoverfly_distM, FDtype = "tau_value", FDtau = NULL)
#'
#' info_FDtau_inc$DataInfo
#' info_FDtau_inc$Summary_Coverage
#'
#'
#' ## (Data Information) Functional diversity for incidence data when all
#' ## threshold levels from 0 to 1 are considered
#' # See Example 6 in the iNEXT.meta vignette for an FD AUC example using Spider.
#' 
#' data("Hoverfly")
#' data("Hoverfly_distM")
#' info_FDAUC_inc <- DataInfobeta3Dmeta(data = Hoverfly,
#'                                      diversity = 'FD', datatype = 'incidence_raw',
#'                                      FDdistM = Hoverfly_distM, FDtype = "AUC")
#'                                      
#' info_FDAUC_inc$DataInfo
#' info_FDAUC_inc$Summary_Coverage
#' }
#' @export

DataInfobeta3Dmeta <- function(data, diversity = "TD", datatype = "abundance",
                               PDtree = NULL, PDreftime = NULL, PDtype = "AUC",
                               FDdistM = NULL, FDtype = "AUC", FDtau = NULL){
  
  if (datatype == "abundance") {
    data <- dplyr::select(
      data, -dplyr::any_of(c("Patch", "Sampling_unit"))
    )
  } else if (datatype == "incidence_raw") {
    data <- dplyr::select(
      data, -dplyr::any_of("Sampling_unit")
    )
  } else {
    stop("Invalid datatype.", call. = FALSE)
  }
  
  # Remove all-zero patches and detect site-treatment groups lost after filtering
  if (datatype == "incidence_raw") {
    group_cols <- names(data)[1:2]
    original_groups <- dplyr::distinct(data, dplyr::across(dplyr::all_of(group_cols)))
    species_cols <- names(data)[-(1:3)]
    data <- data |>
      dplyr::group_by(dplyr::across(1:3)) |>
      dplyr::filter(any(rowSums(as.matrix(dplyr::pick(dplyr::all_of(species_cols)))) > 0)) |>
      dplyr::ungroup()
    remaining_groups <- dplyr::distinct(data, dplyr::across(dplyr::all_of(group_cols)))
    missing_groups <- dplyr::anti_join(original_groups, remaining_groups, by = group_cols)
    if (nrow(missing_groups) > 0)
      stop("Some site-treatment combinations have no non-empty patches.", call. = FALSE)
  }
  
  if (diversity == "PD" && !PDtype %in% c("meanPD", "AUC")) {
    stop(
      "Incorrect phylogenetic diversity type. Please use 'meanPD' or 'AUC'.",
      call. = FALSE
    )
  }
  
  if (diversity == "FD" && !FDtype %in% c("tau_value", "AUC"))
    stop("Incorrect functional diversity type. Please use 'tau_value' or 'AUC'.", call. = FALSE)
  
  if(datatype == "incidence_raw" && diversity == "FD"){
    
    species_names <- names(data)[-(1:3)]
    
    nT <- data |>
      dplyr::group_by(dplyr::across(1:2)) |>
      dplyr::summarise(
        nT = dplyr::n_distinct(.data[[names(data)[3]]]),
        .groups = "drop"
      )
    
    if(any(nT$nT < 8)){
      
      message(
        "For functional diversity, the number of sampling units in your data is fewer than eight, ",
        "which may be insufficient for making reliable inferences using incidence-based models.\n",
        "In this case, species incidence frequencies are treated as abundances, and all outputs ",
        "are generated using abundance-based models.\n",
        "Accordingly, sample coverage estimation is also based on abundance data."
      )
      
      data <- data |>
        dplyr::group_by(dplyr::across(1:2)) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(species_names),
            ~ sum(.x, na.rm = TRUE)
          ),
          .groups = "drop"
        )
      
      datatype <- "abundance"
    }
  }
  
  if (datatype == "abundance"){
    
    output <- data |>
      group_by(across(1:2)) |>
      summarise(
        mat = list(
          t(pick(everything())) %>%
            DataInfobeta3D(data = ., diversity, datatype = "abundance", PDtype,
                           PDtree, PDreftime, FDdistM, FDtype, FDtau) |>
            filter(Assemblage %in% c("Pooled assemblage", "Joint assemblage"))
        ),
        .groups = "drop"
      ) |>
      unnest(mat) |>
      select(-Dataset)
    
    C_vector <- output |>
      filter(Assemblage == "Joint assemblage") |>
      pull(`SC(2n)`)
    
  } else if (datatype == "incidence_raw") {
    
    output <- data |>
      group_by(across(1:3)) |>
      summarise(
        mat = list(t(pick(everything()))),
        .groups = "drop"
      ) |>
      group_by(across(1:2)) |>
      summarise(
        mat = list(mat),
        .groups = "keep"
      ) %>%
      reframe(
        mat %>%
          DataInfobeta3D(data = ., diversity, datatype = "incidence_raw", PDtype, 
                         PDtree , PDreftime, FDdistM, FDtype, FDtau) |>
          filter(Assemblage %in% c("Pooled assemblage", "Joint assemblage"))) |>
      select(-Dataset)
    
    C_vector <- output |>
      filter(Assemblage == "Joint assemblage") |>
      pull(`SC(2T)`)
  }else{
    stop('Please specify datatype = "abundance" or "incidence_raw".', call. = FALSE)
  }
  
  summary_C <- stats::quantile(C_vector, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm = TRUE)
  
  summary_C <- as.data.frame(t(summary_C))
  
  colnames(summary_C) <- c("C_Min", "C_25quantile", "C_50quantile", "C_75quantile", "C_Max")
  
  return(list(DataInfo = output, Summary_Coverage = summary_C))
}

#' Meta-analysis of standardized 3D diversity between two treatments
#'
#' \code{iNEXTbeta3Dmeta()} performs a meta-analysis of treatment effects on coverage-standardized taxonomic, phylogenetic, or functional gamma, alpha, and beta diversity, as well as four dissimilarity indices, across multiple studies/sites. Treatment effects are quantified using either differences or response ratios and combined across studies/sites using a fixed-effect meta-analysis model.
#'
#'
#' @param data For \code{datatype = "abundance"}, a data.frame with site and treatment as the first two columns, optionally followed by a \code{Patch} column, and then species abundance columns. Each row represents an assemblage. For \code{datatype = "incidence_raw"}, a data.frame with site,
#' treatment, and patch as the first three columns, optionally followed by a \code{Sampling_unit} column, and then species incidence columns (0 or 1). Each row represents a sampling unit within a patch.
#' For incidence data, patches with no observed species are removed before diversity estimation; an error is raised if any site-treatment group is lost.
#'
#' @param diversity Selection of the diversity type: \code{diversity = "TD"} for taxonomic diversity, \code{diversity = "PD"} for phylogenetic diversity, and \code{diversity = "FD"} for functional diversity.
#'
#' @param order.q A numerical value or vector specifying the diversity order. Default is \code{order.q = 0}.
#'
#' @param effect Effect measure used to compare the two treatments. \code{effect = "difference"} calculates the diversity estimate of the first treatment minus that of the second treatment. \code{effect = "ratio"} calculates the diversity estimate of the first treatment divided by that of the second treatment. For response ratios, meta-analysis is performed on the log-ratio scale, and the pooled estimate and confidence interval are back-transformed to the original ratio scale. Default is \code{effect = "ratio"}.
#'
#' @param datatype Type of input data: individual-based abundance data (\code{datatype = "abundance"}) or sampling-unit-based incidence data (\code{datatype = "incidence_raw"}). For incidence raw data, species incidence is recorded as 0 (non-detection) or 1 (detection).
#'
#' @param level A single numeric value specifying the target sample coverage level used to compute standardized diversity and dissimilarity estimates. The specified coverage level must be between 0 and 1. If \code{level = NULL}, the smallest value greater than 0.5 among the five alpha-coverage summary statistics returned by \code{DataInfobeta3Dmeta()} is used; if all five are at most 0.5, their maximum is used. A user-specified level at most 0.5 is replaced by the default level.
#'
#' @param nboot An integer of at least 2 specifying the number of bootstrap replications used to assess sampling uncertainty. Unlike \code{iNEXTbeta3D()}, which allows \code{nboot = 0}, this function requires at least two replications to obtain variance estimates for meta-analysis. Default is \code{nboot = 10}. For more accurate uncertainty estimates, a larger number of replications (e.g., 100 or 200) is recommended.
#'
#' @param treatment_order A character vector specifying the order of the two treatments. For \code{effect = "difference"}, the effect size is the first
#' treatment minus the second treatment. For \code{effect = "ratio"}, the effect size is the first treatment divided by the second treatment.
#'
#' @param conf A positive number smaller than 1 specifying the confidence level. Default is 0.95.
#'
#' @param PDtree Required when \code{diversity = "PD"}. A phylogenetic tree (an R \code{phylo} object) spanning all observed species across the input data. The tree can be imported from a file in Newick format.
#'
#' @param PDreftime Argument only for \code{diversity = "PD"} and \code{PDtype = "meanPD"}. A numerical value specifying a fixed reference time for phylogenetic diversity. Default is \code{PDreftime = NULL}, in which case the root age of the phylogenetic tree spanning the pooled observed species is used.
#'
#' @param PDtype Argument only for \code{diversity = "PD"}. Specifies the phylogenetic diversity measure: \code{PDtype = "meanPD"} for the effective number of equally divergent lineages at a specified reference time, or \code{PDtype = "AUC"} for an overall measure obtained by integrating and averaging mean phylogenetic diversity over all reference times from zero to the root age of the phylogenetic tree. Default is \code{PDtype = "AUC"}.
#'
#' @param PDcut_number Argument only for \code{diversity = "PD"} and \code{PDtype = "AUC"}. An integer specifying the number of equally spaced integration points used to compute the AUC over the interval from zero to the root age of the phylogenetic tree. Default is \code{PDcut_number = 30}, and the minimum is 3.
#'
#' @param FDdistM Required when \code{diversity = "FD"}. A species pairwise distance matrix covering all species observed across the input data.
#'
#' @param FDtype Argument only for \code{diversity = "FD"}. Specifies the functional diversity measure: \code{FDtype = "tau_value"} for FD at a specified threshold of functional distinctiveness, or \code{FDtype = "AUC"} for an overall FD measure obtained by integrating the tau-profile over threshold values from 0 to 1. Default is \code{FDtype = "AUC"}.
#'
#' @param FDtau Argument only for \code{diversity = "FD"} and \code{FDtype = "tau_value"}. A numeric value between 0 and 1 specifying the threshold of functional distinctiveness used to compute FD. If \code{FDtau = NULL} (default), the underlying \code{iNEXTbeta3D()} function computes the mean pairwise functional distance between two individuals randomly selected from the pooled data within each treatment (i.e., Rao\'s quadratic entropy). To compare treatments at a common threshold, specify \code{FDtau} explicitly.
#'
#' @param FDcut_number Argument only for \code{diversity = "FD"} and \code{FDtype = "AUC"}. An integer specifying the number of equally spaced integration points used to compute the AUC over the interval from 0 to 1. Default is \code{FDcut_number = 30}, and the minimum is 3.
#'
#'
#'
#' @import iNEXT.3D
#' @importFrom dplyr group_by group_split summarise filter select mutate ungroup rename rename_with bind_rows across everything
#' @importFrom purrr map map2_dfr map_dfr map_df imap
#' @importFrom stringr str_to_title
#' @importFrom tibble tibble
#' @importFrom stats qnorm
#' @importFrom metafor rma
#' @importFrom iNEXT.beta3D iNEXTbeta3D
#'
#' @return
#' \code{iNEXTbeta3Dmeta()} returns a list containing seven summary data frames for gamma, alpha, and beta diversity and four dissimilarity measures (\code{1-C}, \code{1-U}, \code{1-V}, and \code{1-S}), together with the target standardized sample coverage (\code{SC}).
#' Each summary data frame contains site-specific treatment effects and a pooled fixed-effect meta-analysis estimate. The output columns are described below:
#' \describe{
#'   \item{Site}{The name of the study/site. The row labeled \code{Meta-analysis} contains the pooled estimate.}
#'   \item{Ratio / Difference}{The estimated treatment effect. For \code{effect = "ratio"}, the estimate is the ratio of the first treatment to the second treatment; for \code{effect = "difference"}, it is the first treatment minus the second treatment.}
#'   \item{SE}{The standard error of the estimated treatment effect on the analysis scale (log-ratio for response ratios and original scale for differences).}
#'   \item{LCL, UCL}{The lower and upper confidence limits at the specified confidence level. For response ratios, limits are back-transformed to the original ratio scale.}
#'   \item{Order.q}{The diversity order q.}
#'   \item{Diversity}{The diversity type: \code{"TD"}, \code{"meanPD"}, \code{"PD_AUC"}, \code{"FD_tau"}, or \code{"FD_AUC"}, as applicable.}
#'   \item{Treatment columns}{The standardized diversity or dissimilarity estimates for the two treatments. Column names follow \code{treatment_order}; pooled meta-analysis values are not provided for these columns.}
#'   \item{Weight}{The percentage inverse-variance weight assigned to each site in the fixed-effect meta-analysis.}
#'   \item{Effect}{The effect measure used to compare the two treatments, either \code{"Ratio"} or \code{"Difference"}.}
#' }
#' The \code{SC} component is a numeric value representing the target standardized sample coverage.
#'
#' @examples
#'
#' ## Taxonomic diversity for abundance data
#' # See Examples 1 (ratio) and 2 (difference) in the iNEXT.meta vignette for details and forest plots.
#' # Coverage-based standardized TD
#' 
#' data("Spider")
#' TD_abu <- iNEXTbeta3Dmeta(data = Spider, diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#'
#' TD_abu
#'
#' \dontrun{
#' ## Phylogenetic diversity for abundance data
#' # See Example 5 in the iNEXT.meta vignette for details and forest plots (PD AUC).
#' # Coverage-based standardized PD
#' 
#' data("Spider")
#' data("Spider_tree")
#' PD_abu <- iNEXTbeta3Dmeta(data = Spider, diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Spider_tree, PDtype = "AUC", PDcut_number = 30)
#' PD_abu
#'
#' ## Functional diversity for abundance data
#' # See Example 6 in the iNEXT.meta vignette for details and forest plots (FD AUC).
#' # Coverage-based standardized FD
#' 
#' data("Spider")
#' data("Spider_distM")
#' FD_abu <- iNEXTbeta3Dmeta(data = Spider, diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Spider_distM, FDtype = "AUC", FDcut_number = 30)
#' FD_abu
#' }
#' ## Taxonomic diversity for incidence data
#' # See Examples 3 (ratio) and 4 (difference) in the iNEXT.meta vignette for details and forest plots.
#' # Coverage-based standardized TD
#' 
#' data("Hoverfly")
#' TD_inc <- iNEXTbeta3Dmeta(data = Hoverfly, diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#' TD_inc
#'
#' \dontrun{
#' ## Phylogenetic diversity for incidence data
#' # See Example 5 in the iNEXT.meta vignette for a PD AUC analysis using Spider abundance data.
#' # Coverage-based standardized PD
#' 
#' data("Hoverfly")
#' data("Hoverfly_tree")
#' PD_inc <- iNEXTbeta3Dmeta(data = Hoverfly, diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Hoverfly_tree, PDtype = "AUC", PDcut_number = 30)
#' PD_inc
#'
#' ## Functional diversity for incidence data
#' # See Example 6 in the iNEXT.meta vignette for an FD AUC analysis using Spider abundance data.
#' # Coverage-based standardized FD
#' 
#' data("Hoverfly")
#' data("Hoverfly_distM")
#' FD_inc <- iNEXTbeta3Dmeta(data = Hoverfly, diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Hoverfly_distM, FDtype = "AUC", FDcut_number = 30)
#' FD_inc
#' }
#' @export

iNEXTbeta3Dmeta <- function(data, diversity = "TD", order.q = 0, effect = "ratio",
                            datatype = "incidence_raw", level = NULL, nboot = 10,
                            treatment_order, conf = 0.95, PDtree = NULL,
                            PDreftime = NULL, PDtype = "AUC", PDcut_number = 30,
                            FDdistM = NULL, FDtype = "AUC", FDtau = NULL,
                            FDcut_number = 30) 
{
  
  # Metadata are not species. Patch is the incidence assemblage ID.
  if (datatype == "abundance") {
    data <- dplyr::select(data, -dplyr::any_of(c("Patch", "Sampling_unit")))
  } else if (datatype == "incidence_raw") {
    data <- dplyr::select(data, -dplyr::any_of("Sampling_unit"))
  } else stop("Invalid datatype.", call. = FALSE)
  
  # Remove all-zero patches and detect site-treatment groups lost after filtering
  if (datatype == "incidence_raw") {
    group_cols <- names(data)[1:2]
    original_groups <- dplyr::distinct(data, dplyr::across(dplyr::all_of(group_cols)))
    species_cols <- names(data)[-(1:3)]
    data <- data |>
      dplyr::group_by(dplyr::across(1:3)) |>
      dplyr::filter(any(rowSums(as.matrix(dplyr::pick(dplyr::all_of(species_cols)))) > 0)) |>
      dplyr::ungroup()
    remaining_groups <- dplyr::distinct(data, dplyr::across(dplyr::all_of(group_cols)))
    missing_groups <- dplyr::anti_join(original_groups, remaining_groups, by = group_cols)
    if (nrow(missing_groups) > 0)
      stop("Some site-treatment combinations have no non-empty patches.", call. = FALSE)
  }
  
  if (diversity == "PD" && !PDtype %in% c("meanPD", "AUC"))
    stop("Incorrect phylogenetic diversity type. Please use 'meanPD' or 'AUC'.", call. = FALSE)
  if (diversity == "FD" && !FDtype %in% c("tau_value", "AUC"))
    stop("Incorrect functional diversity type. Please use 'tau_value' or 'AUC'.", call. = FALSE)
  if (diversity == "PD" && PDtype == "AUC" &&
      (!is.numeric(PDcut_number) || length(PDcut_number) != 1L || is.na(PDcut_number) ||
       !is.finite(PDcut_number) || PDcut_number < 3 || PDcut_number %% 1 != 0))
    stop("PDcut_number must be an integer of at least 3.", call. = FALSE)
  if (diversity == "FD" && FDtype == "AUC" &&
      (!is.numeric(FDcut_number) || length(FDcut_number) != 1L || is.na(FDcut_number) ||
       !is.finite(FDcut_number) || FDcut_number < 3 || FDcut_number %% 1 != 0))
    stop("FDcut_number must be an integer of at least 3.", call. = FALSE)
  if (!is.numeric(nboot) || length(nboot) != 1L || is.na(nboot) ||
      !is.finite(nboot) || nboot < 2 || nboot %% 1 != 0)
    stop("nboot must be an integer greater than one.", call. = FALSE)
  if (!is.null(level) && (!is.numeric(level) || length(level) != 1L ||
                          is.na(level) || !is.finite(level) || level <= 0 || level > 1))
    stop("level must be a single value between 0 and 1.", call. = FALSE)
  if (length(treatment_order) != 2 ||
      !all(treatment_order %in% unique(data[[2]])))
    stop("Invalid treatment_order.", call. = FALSE)
  if (any(!vapply(split(data[[2]], data[[1]]),
                  function(z) all(treatment_order %in% z), logical(1))))
    stop("Each site must contain both treatments.", call. = FALSE)
  
  # Preserve the original FD fallback for fewer than eight incidence patches.
  if (datatype == "incidence_raw" && diversity == "FD") {
    nT <- dplyr::summarise(dplyr::group_by(data, dplyr::across(1:2)),
                           nT = dplyr::n_distinct(.data[[names(data)[3]]]), .groups = "drop")
    if (any(nT$nT < 8)) {
      message("FD: fewer than eight patches; incidence frequencies are treated as abundances.")
      species <- names(data)[-(1:3)]
      data <- dplyr::summarise(dplyr::group_by(data, dplyr::across(1:2)),
                               dplyr::across(dplyr::all_of(species), ~ sum(.x, na.rm = TRUE)),
                               .groups = "drop")
      datatype <- "abundance"
    }
  }
  
  info <- DataInfobeta3Dmeta(data, diversity, datatype, PDtree, PDreftime,
                             PDtype, FDdistM, FDtype, FDtau)
  
  C_available <- as.numeric(
    unlist(info$Summary_Coverage, use.names = FALSE)
  )
  
  if (is.null(level) || level <= 0.5) {
    level <- if (any(C_available > 0.5, na.rm = TRUE))
      min(C_available[C_available > 0.5], na.rm = TRUE) else
        max(C_available, na.rm = TRUE)
  }
  
  get_sites <- function(treatment) {
    
    x <- data[data[[2]] == treatment, , drop = FALSE]
    
    split(x, x[[1]]) |>
      lapply(function(s) {
        
        if (datatype == "abundance") {
          return(t(dplyr::select(s, -(1:2))))
        }
        
        # Each Patch represents an incidence assemblage
        patch_list <- split(s, s[[3]]) |>
          lapply(function(p) {
            t(dplyr::select(p, -(1:3)))
          })
        
        return(patch_list)
      })
  }
  
  calc <- function(treatment) {
    res <- iNEXT.beta3D::iNEXTbeta3D(
      get_sites(treatment), diversity = diversity, q = order.q,
      datatype = datatype, level = level, nboot = nboot,
      PDtree = PDtree, PDreftime = PDreftime, PDtype = PDtype,
      PDcut_number = PDcut_number, FDdistM = FDdistM, FDtype = FDtype,
      FDtau = FDtau, FDcut_number = FDcut_number)
    purrr::map(res, ~ purrr::imap(.x, function(z, key) {
      z <- dplyr::mutate(z, Type = stringr::str_to_title(key))
      dplyr::rename_with(z, ~ ifelse(.x %in% c("Gamma", "Alpha", "Beta", "Dissimilarity"),
                                     "qD", .x))
    })) |> purrr::map_dfr(dplyr::bind_rows)
  }
  a <- calc(treatment_order[1]); b <- calc(treatment_order[2])
  types <- c("Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", "1-S")
  
  output <- lapply(types, function(type) {
    x <- merge(dplyr::filter(a, Type == type), dplyr::filter(b, Type == type),
               by = c("Dataset", "Order.q"), suffixes = c("_T1", "_T2"))
    if (!nrow(x)) return(data.frame())
    if (effect == "ratio") {
      x$Effect <- log(x$qD_T1 / x$qD_T2)
      x$vi <- log1p((x$s.e._T1 / x$qD_T1)^2 + (x$s.e._T2 / x$qD_T2)^2)
    } else if (effect == "difference") {
      x$Effect <- x$qD_T1 - x$qD_T2
      x$vi <- x$s.e._T1^2 + x$s.e._T2^2
    } else stop("Invalid effect.", call. = FALSE)
    x <- x[is.finite(x$Effect) & is.finite(x$vi) & x$vi > 0, , drop = FALSE]
    if (!nrow(x)) return(data.frame())
    
    purrr::map_dfr(split(x, x$Order.q), function(z) {
      fit <- metafor::rma(yi = z$Effect, vi = z$vi, method = "FE",
                          level = conf * 100)
      weight <- (1 / z$vi) / sum(1 / z$vi) * 100
      zscore <- stats::qnorm(1 - (1 - conf) / 2)
      eff <- c(z$Effect, as.numeric(fit$beta))
      se <- c(sqrt(z$vi), as.numeric(fit$se))
      lower <- c(z$Effect - zscore * sqrt(z$vi), as.numeric(fit$ci.lb))
      upper <- c(z$Effect + zscore * sqrt(z$vi), as.numeric(fit$ci.ub))
      is_ratio <- identical(effect, "ratio")
      if (is_ratio) {
        estimates <- exp(eff)
        lower <- exp(lower); upper <- exp(upper)
      } else estimates <- eff
      ans <- data.frame(
        Site = c(as.character(z$Dataset), "Meta-analysis"),
        Estimate = estimates, SE = se, LCL = lower, UCL = upper,
        Order.q = rep(z$Order.q[1], nrow(z) + 1L),
        Diversity = rep(as.character(z$Diversity_T1[1]), nrow(z) + 1L),
        T1 = c(z$qD_T1, NA_real_), T2 = c(z$qD_T2, NA_real_),
        Weight = c(weight, 100),Effect = if (is_ratio) "Ratio" else "Difference",
        check.names = FALSE)
      names(ans)[2] <- if (is_ratio) "Ratio" else "Difference"
      names(ans)[8:9] <- treatment_order
      ans
    })
  })
  names(output) <- paste0("Summary_", types)
  output$SC <- level
  output
}

#' Forest plot for meta-analysis of standardized 3D diversity
#'
#' \code{ggiNEXTmeta()} generates a forest plot for the meta-analysis results obtained from \code{iNEXTbeta3Dmeta()}. The function automatically identifies whether the effect measure is a difference or a response ratio. For differences, the null effect is zero and a linear scale is used. For response ratios, the null effect is one and a logarithmic scale is used.
#' 
#' 
#' @param output The output of the \code{iNEXTbeta3Dmeta} function. The sample coverage used for standardization is automatically extracted from the output and displayed as a percentage in the forest plot title.
#' @param order.q A numerical value specifying the diversity order to be displayed in the forest plot. The specified value must be available in the output of \code{iNEXTbeta3Dmeta()}.
#' @param num_round A numerical value that the values show on the plot are rounded to the specified value of decimal places.
#' @param range An optional numeric vector of length two specifying the lower and upper limits for clipping confidence intervals. For response-ratio plots, both values must be positive. Default is \code{NULL}.
#' @param type A character string specifying the diversity or dissimilarity measure to be displayed: \code{"Gamma"}, \code{"Alpha"}, \code{"Beta"}, \code{"1-C"}, \code{"1-U"}, \code{"1-V"}, or \code{"1-S"}.
#'
#' @importFrom dplyr filter
#' @importFrom forestplot forestplot fp_set_style fp_add_header fp_append_row fp_decorate_graph fp_set_zebra_style fp_add_lines fpTxtGp
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom grid grid.text gpar
#' @importFrom gridExtra grid.arrange
#'
#'
#'
#' @return A forest plot displaying the treatment effect for each study/site and the pooled meta-analysis estimate.
#' Depending on the output from \code{iNEXTbeta3Dmeta}, the effect measure is displayed as either a difference or a response ratio.
#'
#'
#' @examples
#' ## Taxonomic diversity for abundance data
#' # See Examples 1 (ratio) and 2 (difference) in the iNEXT.meta vignette for details and forest plots.
#' # Coverage-based standardized TD
#' 
#' data("Spider")
#' TD_abu <- iNEXTbeta3Dmeta(data = Spider, diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#' ggiNEXTmeta(TD_abu, order.q = 0, num_round = 3,
#'             range = NULL, type = "Gamma")
#'
#' \dontrun{
#' ## Phylogenetic diversity for abundance data
#' # See Example 5 in the iNEXT.meta vignette for details and forest plots (PD AUC).
#' # Coverage-based standardized PD
#' 
#' data("Spider")
#' data("Spider_tree")
#' PD_abu <- iNEXTbeta3Dmeta(data = Spider, diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Spider_tree, PDtype = "AUC", PDcut_number = 30)
#' ggiNEXTmeta(PD_abu, order.q = 0, num_round = 3,
#'             range = NULL, type = "Gamma")
#'
#'
#' ## Functional diversity for abundance data
#' # See Example 6 in the iNEXT.meta vignette for details and forest plots (FD AUC).
#' # Coverage-based standardized FD
#' 
#' data("Spider")
#' data("Spider_distM")
#' FD_abu <- iNEXTbeta3Dmeta(data = Spider, diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Spider_distM, FDtype = "AUC", FDcut_number = 30)
#' ggiNEXTmeta(FD_abu, order.q = 0, num_round = 3,
#'             range = NULL, type = "Gamma")
#'
#' 
#' ## Taxonomic diversity for incidence data
#' # See Examples 3 (ratio) and 4 (difference) in the iNEXT.meta vignette for details and forest plots.
#' # Coverage-based standardized TD
#' 
#' data("Hoverfly")
#' TD_inc <- iNEXTbeta3Dmeta(data = Hoverfly, diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#' ggiNEXTmeta(TD_inc, order.q = 0, num_round = 3,
#'             range = NULL, type = "Gamma")
#'
#'
#' ## Phylogenetic diversity for incidence data
#' # See Example 5 in the iNEXT.meta vignette for a PD AUC analysis using Spider abundance data.
#' # Coverage-based standardized PD
#' 
#' data("Hoverfly")
#' data("Hoverfly_tree")
#' PD_inc <- iNEXTbeta3Dmeta(data = Hoverfly, diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Hoverfly_tree, PDtype = "AUC", PDcut_number = 30)
#' ggiNEXTmeta(PD_inc, order.q = 0, num_round = 3,
#'             range = NULL, type = "Gamma")
#'
#'
#' ## Functional diversity for incidence data
#' # See Example 6 in the iNEXT.meta vignette for an FD AUC analysis using Spider abundance data.
#' # Coverage-based standardized FD
#' 
#' data("Hoverfly")
#' data("Hoverfly_distM")
#' FD_inc <- iNEXTbeta3Dmeta(data = Hoverfly, diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Hoverfly_distM, FDtype = "AUC", FDcut_number = 30)
#' ggiNEXTmeta(FD_inc, order.q = 0, num_round = 3,
#'             range = NULL, type = "Gamma")
#'
#' }
#' @export

ggiNEXTmeta <- function(output, order.q, num_round = 3, range = NULL, type = NULL)
{
  
  level <- output$SC
  
  if(is.null(type)){
    stop('Please specify type = "Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", or "1-S".', call. = FALSE)
  }
  
  output <- output[paste("Summary", type, sep = "_")]
  meta_output <- output[[1]]
  
  meta_output <- meta_output |>
    dplyr::filter(.data$Order.q == .env$order.q)
  
  if(nrow(meta_output) == 0){
    stop(paste0("order.q = ", order.q, " is not available in the output."), call. = FALSE)
  }
  
  if("Ratio" %in% names(meta_output)){
    effect_type <- "Ratio"
    effect_value <- meta_output$Ratio
    xlog_value <- TRUE
    zero_value <- 1
  }else if("Difference" %in% names(meta_output)){
    effect_type <- "Difference"
    effect_value <- meta_output$Difference
    xlog_value <- FALSE
    zero_value <- 0
  }else{
    stop('The output must contain either a "Ratio" or "Difference" column.', call. = FALSE)
  }
  
  if(!is.null(range)){
    if(length(range) != 2) stop('"range" must contain two numerical values.', call. = FALSE)
    if(effect_type == "Ratio" && any(range <= 0)){
      stop('For ratio plots, both values in "range" must be greater than 0.', call. = FALSE)
    }
  }
  
  # meta_note_q <- paste0(
  #   meta_output[nrow(meta_output), 1], " (Q = ", meta_note$Q_val,
  #   ", df = ", meta_note$df_val, ", p = ", meta_note$p_val,
  #   "; I2 = ", meta_note$I2_val, ", tau2 = ", meta_note$tau2_val, ")"
  # )
  
  name <- colnames(meta_output)[1]
  name_treat <- colnames(meta_output)[c(8, 9)]
  ll <- nrow(meta_output)
  name_site <- meta_output[, 1][-ll]
  diversity <- meta_output$Diversity[1]
  order <- meta_output$Order.q[1]
  
  if(is.null(level)){
    plot_title <- paste0(
      "[", effect_type, "] ", diversity, " ", type
    )
  }else{
    level_show <- if(level <= 1) round(level * 100, 2) else level
    
    plot_title <- paste0(
      "[", effect_type, "] ", diversity, " ", type,
      " (SC = ", level_show, "%)"
    )
  }
  
  forestplot_q <- tibble::tibble(
    mean = effect_value[-ll],
    lower = meta_output$LCL[-ll],
    upper = meta_output$UCL[-ll],
    study = name_site,
    q_T1 = round(meta_output[, 8][-ll], num_round),
    q_T2 = round(meta_output[, 9][-ll], num_round),
    Effect = round(effect_value[-ll], num_round),
    LCL = round(meta_output$LCL[-ll], num_round),
    UCL = round(meta_output$UCL[-ll], num_round),
    Weight = paste0(round(meta_output$Weight[-ll], 2), "%")
  )
  
  if(is.null(range)){
    plot_range <- c(meta_output$LCL, meta_output$UCL, zero_value)
    plot_range <- plot_range[is.finite(plot_range)]
    
    if(effect_type == "Ratio"){
      plot_range <- plot_range[plot_range > 0]
    }
    
    plot_range <- base::range(plot_range)
  }else{
    plot_range <- range
  }
  
  fig <- forestplot_q |>
    forestplot::forestplot(
      labeltext = c(study, q_T1, q_T2, Effect, LCL, UCL, Weight),
      title = plot_title,
      clip = plot_range,
      zero = zero_value,
      xlog = xlog_value
    ) |>
    forestplot::fp_set_style(
      txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1))
    ) |>
    forestplot::fp_add_header(
      study = c("", name),
      q_T1 = c(paste0(diversity, " (q = ", order, ")"), name_treat[1]),
      q_T2 = c(paste0(diversity, " (q = ", order, ")"), name_treat[2]),
      Effect = c("", effect_type),
      LCL = c("", "LCL"),
      UCL = c("", "UCL"),
      Weight = c("", "Weight")
    ) |>
    forestplot::fp_append_row(
      mean = effect_value[ll],
      lower = meta_output$LCL[ll],
      upper = meta_output$UCL[ll],
      study = "Meta analysis",
      Effect = round(effect_value[ll], num_round),
      LCL = round(meta_output$LCL[ll], num_round),
      UCL = round(meta_output$UCL[ll], num_round),
      is.summary = TRUE
    ) |>
    forestplot::fp_decorate_graph(graph.pos = 4) |>
    forestplot::fp_set_zebra_style("#EFEFEF") |>
    forestplot::fp_add_lines()
  
  print(fig)
  # grid::grid.text(meta_note_q, x = 0.5, y = 0, just = "bottom", check.overlap = TRUE)
}