#' Estimates a phylogenetic tree from allele count data using The Coalescent Model
#'
#' Estimates a phylogenetic tree from a matrix of allele counts using root distance method under The Coalescent Model.
#'
#' @param mat_allele_count A numeric matrix of allele counts where rows
#' represent populations or samples and columns represent alleles or loci.
#'    Input as allele counts, each should be non-negative integer \eqn{\le} \code{fixed.n.at.tips}.
#' @param theta A positive numeric scalar representing the scaled mutation rate (\eqn{4N\mu}).
#' Used to scale genetic distances when estimating branch lengths. Defaults to \code{1}.
#' @param fixed.n.at.tips A positive integer specifying the fixed sample size assumed at each tip of the tree.
#' Used to correct for sampling effects when computing distances. Defaults to \code{4}.
#' @param newick_br_length_digits A non-negative integer controlling the number of decimal places used when
#' formatting branch lengths in the Newick string output.
#' Defaults to \code{3}.
#' @param n.cores Optional integer specifying the number of CPU cores used
#' for parallel computation. Defaults to \code{detectCores() - 1}.
#'
#'
#' @returns A named list with two elements:
#'   \describe{
#'     \item{unlabeled_tree}{A phylogenetic tree of class \code{"phylo"} with
#'       estimated branch lengths but no tip labels assigned.}
#'     \item{labeled_tree}{A phylogenetic tree of class \code{"phylo"} with
#'       estimated branch lengths and tip labels derived from the row names
#'       of \code{mat_allele_count}.}
#'   }
#'
#' @export
#'
#' @references
#' Peng J, Rajeevan H, Kubatko L, RoyChoudhury A (2021).
#' A fast likelihood approach for estimation of large phylogenies
#' from continuous trait data.
#' \emph{Molecular Phylogenetics and Evolution}, \bold{161}, 107142.
#' \doi{10.1016/j.ympev.2021.107142}

#' @examples
#' # Load built-in example dataset (9 taxa x 2000 loci)
#' data(Human_Allele_Count_Data)
#'
#' # Inspect dimensions: rows = taxa, columns = loci
#' dim(Human_Allele_Count_Data)
#'
#' # Preview first few rows and columns
#' Human_Allele_Count_Data[1:3, 1:5]
#'
#' # Check for missing data
#' anyNA(Human_Allele_Count_Data)
#'
#' \dontrun{
#' # Estimate phylogenetic tree using the Coalescent Model
#' tree <- RD(mat_allele_count = Human_Allele_Count_Data, n.cores = NULL)
#'
#' # Summarize the result
#' print(tree)
#'
#' # Plot the labeled phylogenetic tree
#' plot(tree$labeled_tree)
#' }


RD <- function(mat_allele_count,
                            theta = 1,
                            fixed.n.at.tips = 4,
                            newick_br_length_digits = 3,
                            n.cores = NULL) {

  start_time <- Sys.time()

  # Validate input type
  if (!is.matrix(mat_allele_count) &&
      !is.data.frame(mat_allele_count)) {
    stop("mat_allele_count must be a matrix or data frame")
  }

  # Convert to matrix for consistent handling
  mat_allele_count <- as.matrix(mat_allele_count)

  # Check minimum number of rows (at least 2: one reference + one sample)
  if (nrow(mat_allele_count) < 2) {
    stop("mat_allele_count must have at least 2 rows (P+1 where P >= 1)")
  }

  # Check at least one locus is present
  if (ncol(mat_allele_count) < 1) {
    stop("mat_allele_count must have at least 1 column (locus)")
  }

  # Validate allele count values are numeric
  if (!is.numeric(mat_allele_count)) {
    stop("mat_allele_count must contain numeric values")
  }

  # Validate no negative counts
  if (any(mat_allele_count < 0, na.rm = TRUE)) {
    stop("All allele counts must be non-negative")
  }

  # Validate counts do not exceed theoretical maximum
  if (any(mat_allele_count > fixed.n.at.tips, na.rm = TRUE)) {
    stop(paste("All allele counts must be <= fixed.n.at.tips (",
               fixed.n.at.tips, ")", sep = ""))
  }

  # Validate counts are integers (with tolerance for floating point errors)
  # Using a small tolerance to handle potential floating point issues
  if (any(abs(mat_allele_count - round(mat_allele_count)) > 1e-8, na.rm = TRUE)) {
    stop("All allele counts must be integers")
  }

  # Warn about missing values but allow processing (user should be aware)
  if (any(is.na(mat_allele_count))) {
    warning("mat_allele_count contains NA values")
  }

  # Validate theta parameter
  if (!is.numeric(theta) || theta <= 0) {
    stop("theta must be a positive numeric value")
  }

  # Validate fixed.n.at.tips parameter
  if (!is.numeric(fixed.n.at.tips) || fixed.n.at.tips < 1) {
    stop("fixed.n.at.tips must be a positive integer")
  }

  # Ensure fixed.n.at.tips is an integer (not fractional)
  if (fixed.n.at.tips != round(fixed.n.at.tips)) {
    stop("fixed.n.at.tips must be an integer")
  }

  # Validate branch length digit specification
  if (!is.numeric(newick_br_length_digits) ||
      newick_br_length_digits < 0) {
    stop("newick_br_length_digits must be a non-negative integer")
  }

  # Number of taxa (excluding reference population)
  P = nrow(mat_allele_count) - 1
  # Number of loci
  L = ncol(mat_allele_count)

  cat(paste("Starting tree estimation with", P, "taxa and", L, "loci\n"))
  cat("===========================================================\n\n")

  # STEP 1: ESTIMATE TAU MATRIX
  # Tau matrix represents pairwise distances between populations
  # Based on the coalescent model with given theta parameter
  cat(paste("STEP 1: Estimate tau matrix\n"))
  mat_est_tau_x <- estimate_tau_x_paralell(mat_allele_count, P, theta, fixed.n.at.tips, n.cores)

  # Check if estimation was successful
  if (is.null(mat_est_tau_x)) {
    warning("Tau estimation failed")
  }

  # STEP 2: IDENTIFY ZERO RDD (ROUSSEEUW-DANIEL DISTANCES)
  # Zero RDD indicates pairs of taxa that are likely closely related
  # These form the basis for identifying splits in the tree
  cat(paste("STEP 2: Identify zero RDD\n"))
  array_zero_ID <- identify_zero_rdd(mat_est_tau_x)

  # STEP 3: CONVERT ZERO RDD TO SPLITS
  # Convert pairwise relationships into bipartition splits that define tree topology
  cat(paste("STEP 3: Convert zero RDD to splits\n"))
  split_set_indicator <- zero_rdd_to_splits(array_zero_ID, P)

  # STEP 4: RECONSTRUCT TREE FROM SPLITS
  # Build phylogenetic tree topology from the set of splits
  # Branch lengths are formatted with specified digit precision
  cat(paste("STEP 4: Construct tree from splits\n"))
  est_tree <- reconstruct_tree_from_splits(split_set_indicator, P, newick_br_length_digits)

  # STEP 5: LABEL TREE TIPS
  # Assign original sample/population names to tree tips
  # Uses row names from mat_allele_count if available, otherwise generates numeric labels
  est_tree_labeled <- label_tree_tips(tree = est_tree, mat_allele_count = mat_allele_count)

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  cat("\n")
  cat("===========================================================\n")
  cat(sprintf("TREE ESTIMATION COMPLETE. (%.2f minutes)\n", as.numeric(elapsed)))

  return(list(unlabeled_tree = est_tree, labeled_tree = est_tree_labeled))
}
