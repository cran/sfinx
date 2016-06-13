# Straightforward Filtering INdeX (SFINX) This the package for the use of SFINX as implemented in the Web site interface at http://sfinx.ugent.be/ .
# Affinity purification-mass spectrometry is one of the most common techniques for the analysis of protein-protein interactions, but inferring bona fide
# interactions from the resulting data sets remains notoriously difficult.  We introduce SFINX, a Straightforward Filtering INdeX that identifies
# true-positive protein interactions in a fast, user-friendly, and highly accurate way.


#' The TIP49 dataset of protein interactions (AP-MS).
#'
#' A strictly numeric input matrix with unique proteins as rownames
#' and unique projects as colnames. The cells of the matrix are filled with
#' the associated \strong{peptide counts}. Cells that have no associated
#' peptide counts are filled with a zero. This specific example dataset is
#' derived from the publication of Sardiu et al. (see below). It contains
#' complexes involved in chromatin remodeling and consists of 35 bait-specific
#' projects and 35 negative controls.
#'
#' @format A matrix with 1581 rows (proteins) and 70 variables (projects).
#'
#' @source M.E. Sardiu, Y. Cai, J. Jin, S.K. Swanson, R.C. Conaway, J.W.
#' Conaway, L. Florens, M.P. Washburn Probabilistic assembly of human protein
#' interaction networks from label-free quantitative proteomics. Proc. Natl.
#' Acad. Sci. USA, 105 (2008), pp. 1454 to 1459.
"DataInputExampleFile"

#' A vector with proteins of interest (baits) for the TIP49 dataset.
#'
#' A character vector with all the bait proteins of interest. These proteins
#' are all present as the exact same rownames in the DataInputExampleFile, and
#' these are also the original bait proteins that were used in the original
#' publication of Sardiu et al. (see below).
#'
#' @format A character vector with 26 bait protein entries of interest.
#'
#' @source M.E. Sardiu, Y. Cai, J. Jin, S.K. Swanson, R.C. Conaway, J.W.
#' Conaway, L. Florens, M.P. Washburn Probabilistic assembly of human protein
#' interaction networks from label-free quantitative proteomics. Proc. Natl.
#' Acad. Sci. USA, 105 (2008), pp. 1454 to 1459.
"BaitIdentityExampleFile"

# Function for rowSums that can also handle vector-like input to eliminate unwanted behaviour or error messages
rowSumsNewKT <- function(InputStuff) {
    if (is.array(InputStuff)) {
        return(rowSums(InputStuff))
    } else {
        return(InputStuff)
    }
}



# Function for colSums that can also handle vector-like input to eliminate unwanted behaviour or error messages
colSumsNewKT <- function(InputStuff) {
    if (is.array(InputStuff)) {
        return(colSums(InputStuff))
    } else {
        return(InputStuff)
    }
}


#' SFINX (Straightforward Filtering INdeX).
#'
#' @description \code{sfinx} identifies the true-positive protein interactions
#' in affinity purification - mass spectrometry data sets and in similar
#' co-complex interactomics data sets. It is highly accurate, fast and
#' independent of external data input.
#'
#' It is also available via the Web interface at \url{http://sfinx.ugent.be},
#' which has extra analysis and visualization features.
#'
#' @details For most standard applications of \code{sfinx}, the arguments
#' \code{InputData} and \code{BaitVector} should be sufficient. Any
#' optimization of the other parameters is discouraged and should be explicitly
#' reported upon communication of the results.
#'
#' @param InputData A strictly numeric matrix with unique proteins as rownames
#' and unique projects as colnames. The cells of the matrix are filled with
#' the associated \strong{peptide counts}. Cells that have no associated
#' peptide counts have to be filled with a zero.
#' @param BaitVector A character vector with all the bait proteins of
#' interest. These proteins should all be present as the exact same
#' rownames in \code{InputData}. \code{sfinx} will control this, and it will
#' report possible deviations.
#' @param BackgroundRatio Advanced. A natural number equal or bigger than 2,
#' that specifies the maximal ratio of total considered projects over the
#' amount of bait projects. If this parameter equals for example 5, it will
#' take into account 4 times as much non-bait projects as it uses bait
#' projects. \code{sfinx} will preferably first select the non-bait projects
#' with most peptide counts as negative controls.
#' @param BackgroundIdentity Deprecated. A character string or character
#' vector describing the background projects. \code{"automatic"} is the
#' advised default entry. However, all extra or alternative entries will be
#' matched to the column headers and taken into account when possible.
#' @param BaitInfluence Advanced. A logical. When TRUE, \code{sfinx} uses only
#' the non-bait negative control projects with the biggest amount of data for
#' the calculation of the background, but no negative control projects
#' associated with other baits in the analysis. When FALSE, \code{sfinx} uses
#' both.
#' @param ConstantLimit Advanced. A logical. When TRUE, an internal cut-off is
#' used that is a simplified constant for the actual complete calculation of
#' the binomial equivalent. This is the version of \code{sfinx} that was used
#' in the article (Titeca et al., J. Proteome Res., 2016). When FALSE, the
#' complete calculation of the binomial equivalent is done. Some datasets with
#' many highly abundant proteins can benefit from having this parameter FALSE.
#' @param FWERType Advanced. A character string that equals \code{"B"},
#' \code{"HolmB"} or \code{"Sidak"}. \code{"B"} gives the Bonferroni correction
#' for the family wise error rate (FWER), \code{"HolmB"} gives Holm-Bonferroni
#' correction, and \code{"Sidak"} gives Sidak correction. However, note that
#' these options will only very rarely yield different output.
#'
#' @return \code{sfinx} returns a list with two elements. The first element of
#' the list contains a dataframe with the true-positive protein interactions
#' that were identified by \code{sfinx} in \code{InputData} for the proteins of
#' interest in \code{BaitVector}. The second element of the list contains a
#' string with comments about the output and the underlying data.
#'
#' @examples
#' sfinx(DataInputExampleFile, BaitIdentityExampleFile)
#'
#' sfinx(InputData = DataInputExampleFile, BaitVector =
#' BaitIdentityExampleFile, ConstantLimit = TRUE, FWERType = "Sidak")
#' 
#' 
#' @importFrom stats dbinom
#' @importFrom stats pbinom
#'
#' @export

# The SFINX Algorithm itself
sfinx <- function(InputData, BaitVector, BackgroundRatio = 5, BackgroundIdentity = "automatic", BaitInfluence = FALSE, ConstantLimit = FALSE, FWERType = "B") {

    #### Input & conversion datamatrix input
    TotalDataMatrixSW <- as.matrix(InputData)

    if (mode(TotalDataMatrixSW) != "numeric") {
        stop("Your data are not purely numeric, or something is wrong with the header or rownames. Please check your data.")
    }
    if (is.null(rownames(TotalDataMatrixSW)) | is.null(colnames(TotalDataMatrixSW))) {
        stop("Your input data lack row names or column names. Please check your data.")
    }

    #### Check-up for the presence of the baits in the matrix

    if (!any(rownames(TotalDataMatrixSW) %in% BaitVector)) {
        stop("None of the bait proteins is present in the basic data matrix. Please check whether the names of the baits are appropriate for this dataset.")
    }

    BaitIndices <- which(rownames(TotalDataMatrixSW) %in% BaitVector)
    WorkBaits <- rownames(TotalDataMatrixSW)[BaitIndices]

    #### Check whether BaitInfluence and ConstantLimit are logical

    if (!all(is.logical(BaitInfluence), is.logical(ConstantLimit))) {
        stop("The parameters BaitInfluence and ConstantLimit can only be TRUE or FALSE, no other types of entry are allowed. Please check these parameters.")
    }


    #### Check whether FWERType is in one of the allowed types: 'B', 'HolmB', 'Sidak'


    if (!(FWERType %in% c("B", "HolmB", "Sidak"))) {
        stop("The parameter FWERType can only be B, HolmB or Sidak. Please check this parameter.")
    }


    #### BackgroundRatio has to be above 2 and natural
    if (!is.numeric(BackgroundRatio)) {
        stop("The parameter BackgroundRatio is not numeric. Please correct this parameter.")
    }

    if (BackgroundRatio < 2) {
        stop("The parameter BackgroundRatio has a value below 2 or is not a natural number. Please correct this parameter.")
    }

    #### Set-up for determination of background
    BackgroundIdentityDeterminedAutomatically <- c()

    if (any("automatic" %in% BackgroundIdentity, is.null(BackgroundIdentity))) {
        BackgroundIdentityDeterminedAutomatically <- as.numeric(which(colSumsNewKT(TotalDataMatrixSW[BaitIndices, ]) == 0))
        if ("automatic" %in% BackgroundIdentity) {
            BackgroundIdentity <- BackgroundIdentity[-match("automatic", BackgroundIdentity)]
        }
        if (length(BackgroundIdentity) == 0) {
            BackgroundIdentity <- NULL
        }
    }

    if (!is.null(BackgroundIdentity)) {
        BackgroundIdentity <- match(BackgroundIdentity, colnames(TotalDataMatrixSW))
    }

    if (anyNA(BackgroundIdentity)) {
        stop("The BackgroundIdentity contains elements that are not present as column headings in the input data. Please evaluate the InputData and BackgroundIdentity parameters.")
    }


    BackgroundIdentity <- unique(c(BackgroundIdentity, BackgroundIdentityDeterminedAutomatically))  #automatic == non-bait background; rest == user-input overlapping with bait background

    #### Evaluation of the baits that are (not) used --> user output
    if (length(unique(BaitVector)) == length(BaitIndices)) {
        PrintOutput <- "All baits were found as possible preys."
    } else {
        PrintOutput <- paste("Not all baits were found as possible preys. The following were not found:", paste(setdiff(BaitVector, WorkBaits), collapse = " "),
            ".", sep = "")
    }

    # Control for baits with single peptide counts
    if (length(BaitIndices) > 1) {
        BaitsWithOnly1PeptideCount <- WorkBaits[rowSumsNewKT(TotalDataMatrixSW[BaitIndices, ]) < 2]
    } else {
        BaitsWithOnly1PeptideCount <- WorkBaits[sum(TotalDataMatrixSW[BaitIndices, ]) < 2]
    }


    if (length(BaitsWithOnly1PeptideCount) != 0) {
        PrintOutput <- paste(PrintOutput, " Some baits had only one peptide count, and were not used: ", paste(BaitsWithOnly1PeptideCount, collapse = " "), ". Please, do more experiments for this bait.",
            sep = "")

        WorkBaits <- WorkBaits[rowSumsNewKT(TotalDataMatrixSW[BaitIndices, ]) >= 2]
        BaitIndices <- BaitIndices[rowSumsNewKT(TotalDataMatrixSW[BaitIndices, ]) >= 2]
    }



    #### The core of the algorithm
    OutputListOfSFINX <- lapply(1:length(BaitIndices), function(x) {

        # Bait projects position and amount
        BaitLogicals <- TotalDataMatrixSW[BaitIndices[x], ] != 0
        BaitFrac <- sum(BaitLogicals)

        # Dynamic BaitSom and CutOff1 determination
        if (BaitFrac == 1) {
            BaitSom <- TotalDataMatrixSW[, TotalDataMatrixSW[BaitIndices[x], ] != 0]
        } else {
            BaitSom <- rowSumsNewKT(TotalDataMatrixSW[, TotalDataMatrixSW[BaitIndices[x], ] != 0])
        }

        # Calculation of the NegSom dynamically depending on BackgroundRatio and BackgroundIdentity Negatives from the non-bait projects
        NegativeNonBaitProjects <- BaitLogicals[setdiff(BackgroundIdentity, which(BaitLogicals))]
        NegativeNonBaitProjectsBeforeTheReduction <- TotalDataMatrixSW[, names(NegativeNonBaitProjects[NegativeNonBaitProjects == FALSE])]

        if (length(NegativeNonBaitProjects[NegativeNonBaitProjects == FALSE]) > 1) {
            NegativeNonBaitProjectsBeforeTheReduction <- NegativeNonBaitProjectsBeforeTheReduction[, order(colSumsNewKT(NegativeNonBaitProjectsBeforeTheReduction),
                decreasing = TRUE)]
            NegFrac <- dim(NegativeNonBaitProjectsBeforeTheReduction)[2]
        } else {
            if (length(NegativeNonBaitProjects[NegativeNonBaitProjects == FALSE]) == 1) {
                NegFrac <- 1
            } else {
                NegFrac <- 0
            }
        }

        # Negatives from the bait projects
        NegativeTrueBaitProjects <- BaitLogicals[setdiff(1:length(BaitLogicals), BackgroundIdentity)]
        NegativeTrueBaitProjectsBeforeTheReduction <- TotalDataMatrixSW[, names(NegativeTrueBaitProjects[NegativeTrueBaitProjects == FALSE])]

        if (length(NegativeTrueBaitProjects[NegativeTrueBaitProjects == FALSE]) > 1) {
            NegativeTrueBaitProjectsBeforeTheReduction <- NegativeTrueBaitProjectsBeforeTheReduction[, order(colSumsNewKT(NegativeTrueBaitProjectsBeforeTheReduction),
                decreasing = TRUE)]
            NegFracBPs <- dim(NegativeTrueBaitProjectsBeforeTheReduction)[2]
        } else {
            if (length(NegativeTrueBaitProjects[NegativeTrueBaitProjects == FALSE]) == 1) {
                NegFracBPs <- 1
            } else {
                NegFracBPs <- 0
            }
        }

        # Potential reduction of negatives to calculate the NegSom and TotalFrac
        if (NegFrac > BaitFrac * (BackgroundRatio - 1)) {
            if (BaitInfluence) {
                NegSom <- rowSumsNewKT(NegativeNonBaitProjectsBeforeTheReduction[, 1:(BaitFrac * (BackgroundRatio - 1))])
                TotalFrac <- BaitFrac + BaitFrac * (BackgroundRatio - 1)
            } else {
                NegSom <- rowSumsNewKT(cbind(NegativeNonBaitProjectsBeforeTheReduction[, 1:(BaitFrac * (BackgroundRatio - 1))], NegativeTrueBaitProjectsBeforeTheReduction))
                TotalFrac <- BaitFrac + BaitFrac * (BackgroundRatio - 1) + NegFracBPs
            }
        } else {
            NegSom <- rowSumsNewKT(cbind(NegativeNonBaitProjectsBeforeTheReduction, NegativeTrueBaitProjectsBeforeTheReduction))
            TotalFrac <- BaitFrac + NegFrac + NegFracBPs
        }

        # TotaleSom calculation binomial input
        TotaleSom <- BaitSom + NegSom

        # Binomial calculation
        BinomialSFINXDistributionData <- ifelse(TotaleSom > 1, dbinom(floor(BaitSom), ceiling(TotaleSom), sum(BaitSom)/sum(TotaleSom), log = FALSE), 1)

        if (ConstantLimit) {


            BinomialSFINXCytoscapeOutput <- sort(BinomialSFINXDistributionData[BinomialSFINXDistributionData < (sum(BaitSom)/sum(TotaleSom))^(BaitFrac)], decreasing = TRUE)
            # lower = better
        } else {
            BinomialSFINXCutOffs <- ifelse(ceiling(TotaleSom) >= BaitFrac, (dbinom(floor((TotaleSom + BaitFrac)/2), ceiling(TotaleSom), sum(BaitSom)/sum(TotaleSom),
                log = FALSE) + dbinom(ceiling((TotaleSom + BaitFrac)/2), ceiling(TotaleSom), sum(BaitSom)/sum(TotaleSom), log = FALSE))/2, 0)

            BinomialSFINXCytoscapeOutput <- sort(BinomialSFINXDistributionData[BinomialSFINXDistributionData <= BinomialSFINXCutOffs], decreasing = FALSE)  #BinomialSFINXCytoscapeOutput[[i]]; BinomialSFINXDistributionData[[i]]; BinomialSFINXCutOffs[[i]]]

        }

        # Begin output and pbinom calculation

        if (length(BinomialSFINXCytoscapeOutput) == 0) {
            BinomialInternalSFINXCytoscape <- NULL


            BaitFracCollectorList <- c(BaitFrac, NegFrac, NegFracBPs, 0)

        } else {
            BinomialInternalSFINXCytoscape <- data.frame(WorkBaits[x], as.numeric(BinomialSFINXCytoscapeOutput), names(BinomialSFINXCytoscapeOutput), as.numeric(BinomialSFINXCytoscapeOutput) +
                pbinom(BaitSom[names(BinomialSFINXCytoscapeOutput)], TotaleSom[names(BinomialSFINXCytoscapeOutput)], sum(BaitSom)/sum(TotaleSom), lower.tail = FALSE,
                  log.p = FALSE), stringsAsFactors = FALSE)

            colnames(BinomialInternalSFINXCytoscape) <- c("Baits", "Score", "Preys", "pValue")


            if (FWERType == "B") {
                # Bonferroni correction because we do not want any false positives at all
                BinomialInternalSFINXCytoscape <- BinomialInternalSFINXCytoscape[BinomialInternalSFINXCytoscape[, 4] < 0.05/(length(BinomialSFINXDistributionData) +
                  dim(BinomialInternalSFINXCytoscape)[1]), ]

            }

            if (FWERType == "HolmB") {

                OrderedBinomialInternalSFINXCytoscape <- BinomialInternalSFINXCytoscape[order(BinomialInternalSFINXCytoscape[, 4], decreasing = FALSE), ]

                n <- 1

                for (n in 1:dim(OrderedBinomialInternalSFINXCytoscape)[1]) {
                  if (OrderedBinomialInternalSFINXCytoscape[n, 4] > (0.05/(length(BinomialSFINXDistributionData) + dim(OrderedBinomialInternalSFINXCytoscape)[1] +
                    1 - n))) {
                    break
                  }
                }

                if (n == 1) {
                  BinomialInternalSFINXCytoscape <- NULL
                } else {
                  BinomialInternalSFINXCytoscape <- OrderedBinomialInternalSFINXCytoscape[1:(n - 1), ]
                }
            }

            if (FWERType == "Sidak") {

                BinomialInternalSFINXCytoscape <- BinomialInternalSFINXCytoscape[BinomialInternalSFINXCytoscape[, 4] < (1 - (1 - 0.05)^(1/(length(BinomialSFINXDistributionData) +
                  dim(BinomialInternalSFINXCytoscape)[1]))), ]

            }

            # Back-end check-up-matrix
            BaitFracCollectorList <- c(BaitFrac, NegFrac, NegFracBPs, dim(BinomialInternalSFINXCytoscape)[1])
        }

        list(BinomialInternalSFINXCytoscape, BaitFracCollectorList)


    })



    # Putting together back-end check-up-matrix
    BaitFracCollectorMatrix <- do.call("rbind", lapply(OutputListOfSFINX, function(x) {
        x[[2]]
    }))


    # Control for baits with not enough data
    if (sum(BaitFracCollectorMatrix[, 4] == 0) != 0) {
        PrintOutput <- paste(PrintOutput, " Some baits yielded no interactions because of shortage of underlying data:", paste(WorkBaits[BaitFracCollectorMatrix[,
            4] == 0], collapse = " "), ". Please, try to gather more data for these baits if you are interested in them.", sep = "")
    }

    # Control for baits with not enough negative controls
    DangerBaits <- WorkBaits[BaitFracCollectorMatrix[, 1]/dim(TotalDataMatrixSW)[2] > 0.475]

    if (length(DangerBaits) != 0) {
        PrintOutput <- paste(PrintOutput, " Some baits yielded interactions with a lower confidence, as there are not enough negative controls for them:", paste(DangerBaits,
            collapse = " "), ". Please, use more negative controls.", sep = "")
    }

    # Putting together the output matrix
    StandardBinomialOutputDataFrame <- do.call("rbind", lapply(OutputListOfSFINX, function(x) {
        x[[1]]
    }))
    names(StandardBinomialOutputDataFrame) <- c("Baits", "Scores", "Preys", "pValues")

    return(list(StandardBinomialOutputDataFrame, PrintOutput))
}


# Variables: names; functions: verbs

# Comment lines: explain the why Make names more useful

# Define default value input .onLoad & .onUnload ==> zzz.R
