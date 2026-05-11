#' Reconstruct a phylogenetic tree from inferred splits
#'
#' Reconstructs a rooted phylogenetic tree using a set of inferred splits
#' (bipartitions) among taxa. Each row of the split indicator matrix defines
#' a grouping of taxa into two sides of a split, which is converted into a
#' small tree. These trees are then combined into a final consensus tree
#' using a supertree method.
#'
#' @param new_split_set_indicator A matrix describing inferred splits.
#' Each row represents one split and each column corresponds to a taxon.
#' Entries are coded as:
#' \itemize{
#' \item \code{1} – taxon belongs to the first side of the split
#' \item \code{2} – taxon belongs to the second side of the split
#' \item \code{0} – taxon not involved in the split
#' }
#' @param P Integer. Number of taxa excluding the outgroup.
#' @param newick_br_length_digits Integer specifying the number of digits
#' used when writing branch lengths in Newick format.
#'
#' @return A rooted phylogenetic tree object (class \code{phylo}) representing
#' the reconstructed tree topology.
#'
#' @details
#' The procedure works as follows:
#' \enumerate{
#' \item Each split is converted into a small tree represented in Newick format.
#' \item These trees are stored as a \code{multiPhylo} object.
#' \item A consensus tree is reconstructed using a supertree algorithm.
#' }
#'
#' The outgroup is assumed to be indexed as \code{P + 1}.
#'
#' @importFrom ape read.tree
#' @importFrom phangorn superTree
#'
#' @noRd
#'

reconstruct_tree_from_splits <- function(new_split_set_indicator,
                                         P,
                                         newick_br_length_digits) {
  label1 <- 1:P

  # We will include content of "edited.topology_to_newick_vA5" here.

  ## please check if the last row in the freq matrix is Species1

  ## 2023-02-23 replace "species1" as "P+1" in function multiphylo

  multiphylo <- function(base_tree) {
    for (i in 1:dim(base_tree)[1]) {
      ind3 <- which(base_tree[i, ] == 0)
      ind4 <- which(base_tree[i, ] == 1)
      ind5 <- which(base_tree[i, ] == 2)
      if (length(ind3) == 0) {
        grp1 <- label1[ind4]
        if (length(ind4) > 1) {
          grp1 <- paste("(", paste(grp1, collapse = ","), ")", sep = "")
        }
        grp2 <- label1[ind5]
        if (length(ind5) > 1) {
          grp2 <- paste("(", paste(grp2, collapse = ","), ")", sep = "")
        }
        grp <- paste("(", P + 1, ",(", grp1, ",", grp2, "));", sep = "")
      } else{
        grp1 <- label1[ind4]
        if (length(ind4) > 1) {
          grp1 <- paste("(", paste(grp1, collapse = ","), ")", sep = "")
        }
        grp2 <- label1[ind5]
        if (length(ind5) > 1) {
          grp2 <- paste("(", paste(grp2, collapse = ","), ")", sep = "")
        }
        grp3 <- label1[ind3]
        if (length(ind3) > 1) {
          grp3 <- paste("(", paste(grp3, collapse = ","), ")", sep = "")
        }
        grp <- paste("(", P + 1, ",(", grp1, ",", grp2, ",", grp3, "));", sep =
                       "")
      }
      tree <- read.tree(text = grp)
      mp <- list(tree)
      class(mp) <- "multiPhylo"
      if (i == 1)
        mp2 <- mp
      else
        mp2 <- c(mp2, mp)
    }
    return(mp2)
  }


  mp2 <- multiphylo(new_split_set_indicator)


  est_tree_rdm <- superTree(mp2, rooted = TRUE)

  return(est_tree_rdm)

}




#' Label tree tips using population names
#'
#' Replaces numeric tip labels in a phylogenetic tree with the
#' corresponding population names from the allele count data.
#'
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param mat_allele_count Matrix or data frame of allele counts where
#' rows correspond to populations and row names contain population labels.
#'
#' @return A tree object with updated tip labels.
#'
#' @noRd
#'
label_tree_tips <- function(tree, mat_allele_count) {
  # Create mapping between numeric tip labels and population names
  popnames <- data.frame(
    number = seq_len(nrow(mat_allele_count)),
    label  = rownames(mat_allele_count)
  )

  # Convert tree tip labels to numeric indices
  tip_index <- as.numeric(tree$tip.label)

  # Map indices to original population names
  mapped_labels <- popnames$label[match(tip_index, popnames$number)]

  # Warn if some tips could not be mapped
  if (any(is.na(mapped_labels))) {
    warning("Some tip labels could not be mapped to original names.")
  }

  # Assign mapped labels to the tree
  labeled_tree <- tree
  labeled_tree$tip.label <- mapped_labels

  return(labeled_tree)
}
