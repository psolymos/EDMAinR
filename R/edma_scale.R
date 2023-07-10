#' Scale an EDMA data object
#'
#' This function implements the landmark scaling procedures described in
#' Lele and Cole (1996) which are used to rescale the landmarks for
#' specimens in which the variance-covariance matrices to two populations
#' are unequal.
#'
#' @param x an EDMA data object of class 'edma_data'.
#' @param scale_by string specifying the type of scaling. Valid options are
#' "constant", "endpoints", "geometric mean", "maximum", "median", "sneath".
#' See below for details.
#' @param L1 string specifying first landmark to use if scale_by = "endpoints"
#' @param L2 string specifying second landmark to use if scale_by = "endpoints"
#' @param scale_constant numeric specifying the scaling constant to use for
#' scale_by = "constant"
#'
#' @details
#' \code{constant} Interlandmark distances are scaled by a numeric constant
#' that is applied to all specimens.
#'
#' \code{endpoints} Interlandmark distances are scaled by the distance between
#' a pair of landmarks (\code{L1} and \code{L2}) for each specimen.
#'
#' \code{geometric mean} Interlandmark distances are scaled by the geometric
#' mean of all pairwise distances for each specimen.
#'
#' \code{maximum} Interlandmark distances are scaled by the maximum of all
#' pairwise distances for each specimen.
#'
#' \code{median} Interlandmark distances are scaled by the median of all
#' pairwise distances for each specimen.
#'
#' \code{sneath}  Interlandmark distances are scaled using the method described
#' by Sneath (1967), which uses the square-root of the sum of the squared
#' distances of each landmark to the centroid.
#'
#' @return object of class 'edma_data', with landmarks scaled according to
#' scale_by parameter. See details for details of scaling procedures. The
#' object x are appended the scaling method string and the values used for
#' scaling. The latter can be useful for comparisons, e.g., of geometric
#' means.
#'
#' @export
#'
#' @references
#' Lele, S., and T. M. Cole III. 1996. A new test for shape
#' differences when variance-covariance matrices are unequal.
#' J. Hum. Evol. 31:193-212.
#'
#' Sneath, P. H. A. 1967. Trend-surface analysis of transformation grids.
#' J. Zool. 151:65-122. Wiley.
#'
#' @examples
#' # Following example in Lele and Cole (1996)
#' X <- matrix(c(0, 0, 2, 3, 4, 1), byrow = TRUE, ncol = 2)
#' Y <- matrix(c(0, 0, 3, 3, 3, 0), byrow = TRUE, ncol = 2)
#'
#' XY <- as.edma_data(array(dim = c(3, 2, 2),
#'                          data = cbind(X, Y)))
#' XY_const <- edma_scale(XY, scale_by = "constant", scale_constant = 2)
#' XY_const$data
#' XY_const$scaling_method
#' XY_const$scaling_values
#'
#' XY_endpt <- edma_scale(XY, scale_by = "endpoints", L1 = "L1", L2 = "L3")
#' XY_endpt$data
#' XY_endpt$scaling_method
#' XY_endpt$scaling_values
#'
#' XY_geomean <- edma_scale(XY, scale_by = "geometric mean")
#' XY_geomean$data
#' XY_geomean$scaling_method
#' XY_geomean$scaling_values
#'
#' XY_max <- edma_scale(XY, scale_by = "maximum")
#' XY_max$data
#' XY_max$scaling_method
#' XY_max$scaling_values
#'
#' XY_median <- edma_scale(XY, scale_by = "median")
#' XY_median$data
#' XY_median$scaling_method
#' XY_median$scaling_values
#'
#' XY_sneath <- edma_scale(XY, scale_by = "sneath")
#' XY_sneath$data
#' XY_sneath$scaling_method
#' XY_sneath$scaling_values
#'
#' @name edma_scale
#'
edma_scale <- function(x,
                       scale_by,
                       L1 = NULL,
                       L2 = NULL,
                       scale_constant = NULL) {

    # Check that x is edma_data class object
    if (!inherits(x, "edma_data")) stop("Input must be edma_fit object.")

    # Calculate pairwise distances for all landmarks
    # Several scaling methods use these distances
    pairwise_dists <- function(x) {
        lapply(x$data,
               FUN = function(M) {
                   as.numeric(dist(M, method = "euclidean"))
               })
    }

    scaling <- match.arg(scale_by, c("constant", "endpoints",
                                     "geometric mean", "maximum",
                                     "median", "sneath"))

    if (scale_by == "constant") {
        # Check that scale_constant is supplied
        if (is.null(scale_constant)) {
            stop("scale_constant must be supplied if scaling by 'constant'.")
            }

        # Check that scale_constant is numeric
        if (!inherits(scale_constant, "numeric")) {
            stop("scale_constant must be class 'numeric'")
        }

        # List of constant values the length of the landmark list
        scaling_values <- as.list(rep(scale_constant,
                                      times = length(x$data)))
        names(scaling_values) <- names(x$data)
    }

    if (scale_by == "endpoints") {
        # Check that endpoints are supplied and can be found in the data
        if (is.null(L1) | is.null(L2)) {
            stop("L1 and L2 must be supplied if scaling by 'endpoints'.")}
        if (!(L1 %in% dimnames(x$data[[1]])[[1]])) {
            stop("L1 not found in landmark names.")
        }
        if (!(L2 %in% dimnames(x$data[[1]])[[1]])) {
            stop("L2 not found in landmark names.")
        }

        # Find distance between L1 and L2 for each landmark set
        scaling_values <-
            lapply(x$data,
                   FUN = function(M) {
                       # Coordinates for L1 and L2
                       L1_coords <- M[L1, ]
                       L2_coords <- M[L2, ]
                       as.numeric(dist(rbind(L1_coords, L2_coords),
                            method = "euclidean"))
                   })
    }

    if (scale_by == "geometric mean") {
        # Calculate geometric mean of each set of distances
        # Use log scale to avoid overflow for large N
        scaling_values <-
            lapply(pairwise_dists(x),
                   FUN = function(D){exp(sum(log(D)) / length(D))})
    }

    if (scale_by == "maximum") {
        scaling_values <-
            lapply(pairwise_dists(x),
                   FUN = function(D){max(D)})
    }

    if (scale_by == "median") {
        scaling_values <-
            lapply(pairwise_dists(x),
                   FUN = function(D){median(D)})
    }

    if (scale_by == "sneath") {
        scaling_values <-
            lapply(x$data,
                   FUN = function(M){
                       # Calculate centroid
                       centroid <- colMeans(M)

                       # Distance from landmarks to centroid
                       d <- apply(M,
                                  MARGIN = 1,
                                  FUN = function(r) {
                                      dist(rbind(r, centroid),
                                           method = "euclidean")
                                  })

                       # Return square root of the summed squared distances
                       sqrt(sum(d ^ 2))
                   })
    }

    # Center and scale x$data for each set of landmarks
    # Use a loop to access both x$data and scaling_values
    for (ii in seq_len(length(scaling_values))) {
        # Center
        x$data[[ii]] <-
            t(apply(x$data[[ii]],
                    MARGIN = 1,
                    FUN = function(r) {r - colMeans(x$data[[ii]])}))
        # Scale
        x$data[[ii]] <- x$data[[ii]] / scaling_values[[ii]]
    }

    # Store the scaling method and scaling values
    x$scaling_method <- scale_by
    x$scaling_values <- scaling_values

    return(x)
}
