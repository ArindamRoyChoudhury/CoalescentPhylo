#' Identify zero rooted distance differences (RDD)
#'
#' Determines which rooted distance differences (RDD) are zero based on
#' the estimated internal branch lengths between taxa. For each triplet
#' of taxa \eqn{(i, j, k)}, the function compares the three pairwise
#' internal branch lengths and identifies the sister pair (the pair with
#' the largest \code{tau_x}). The corresponding RDD for the remaining
#' taxon is then set to zero.
#'
#' @param mat_est_tau_x A symmetric matrix of estimated internal branch
#' lengths (\code{tau_x}) between taxa. The matrix dimension is \code{P × P},
#' where \code{P} is the number of taxa excluding the outgroup.
#'
#' @return A 3-dimensional array of size \code{P × P × P}. Entries equal to
#' zero indicate that the corresponding rooted distance difference (RDD)
#' is inferred to be zero for that triplet configuration.
#'
#' @details
#' For each ordered triplet \eqn{(i, j, k)} with \eqn{i < j < k}, the three
#' pairwise internal branch lengths \code{(i,j)}, \code{(i,k)}, and
#' \code{(j,k)} are compared. The pair with the largest value is treated
#' as the sister taxa. If the maximum branch length is positive, the
#' corresponding RDD entry is set to zero.
#'
#' A progress bar is displayed during computation to track processing of
#' taxa.
#'
#' @import utils
#'
#' @noRd

identify_zero_rdd <- function(mat_est_tau_x){
  P <- nrow(mat_est_tau_x)
  ## P is the number of taxa other than the outgroup

  H <- P * (P - 1) * (P - 2) / 6
  ## H is the smallest number of  RD differences that have to be zero in order for it
  ## to agree to a tree-topology


  n_list_ID <- P * (P - 1) * (P - 2) / 2

  array_zero_ID <- array(dim = c(P, P, P), 1)

  # Create progress bar inside the function
  pb <- txtProgressBar(min = 0, max = P, style = 3)

  for (i in 1:P){
    setTxtProgressBar(pb, i)
    for (j in 1:P){
      for (k in 1:P)
      {
        if ((i < j) && (j < k)) {
          # Get the three tau_x values for the triplet
          dists <- c(mat_est_tau_x[i, j], mat_est_tau_x[i, k], mat_est_tau_x[j, k])
          max_idx <- which.max(dists)

          # Only set a zero if there is a clear "winner" (shared branch > 0)
          if (dists[max_idx] > 0) {
            if (max_idx == 1) {
              # i,j are sister
              array_zero_ID[k, i, j] <- 0
              array_zero_ID[k, j, i] <- 0
            } else if (max_idx == 2) {
              # i,k are sister
              array_zero_ID[j, i, k] <- 0
              array_zero_ID[j, k, i] <- 0
            } else if (max_idx == 3) {
              # j,k are sister
              array_zero_ID[i, j, k] <- 0
              array_zero_ID[i, k, j] <- 0
            }
          }
        }


      }
    }
  }
  close(pb)

  return(array_zero_ID)

}


#' Convert zero-RDD information to phylogenetic splits
#'
#' Constructs phylogenetic splits (bipartitions of taxa) from the array of
#' rooted distance differences (RDD) identified as zero. Each triplet
#' configuration that satisfies the zero-RDD condition provides information
#' about the grouping of taxa, which is used to infer candidate splits.
#'
#' The function iteratively builds and merges split sets based on the
#' triplet relationships, removes duplicate splits, and finally adds
#' end-bound splits that are implied by taxa pairs that never appear in
#' the same inferred split.
#'
#' @param array_zero_ID A 3-dimensional array of size \code{P × P × P}
#' indicating which rooted distance differences are zero. Entries equal
#' to zero correspond to triplet configurations supporting a specific
#' sister relationship.
#' @param P Integer. Number of taxa excluding the outgroup.
#'
#' @return A matrix where each row represents a split (bipartition) of the
#' taxa. Columns correspond to taxa, and entries indicate group membership:
#' \itemize{
#' \item \code{1} – taxon belongs to the first side of the split
#' \item \code{2} – taxon belongs to the second side of the split
#' \item \code{0} – taxon not involved in that split
#' }
#'
#' @details
#' The algorithm proceeds in several steps:
#' \enumerate{
#' \item Identify candidate splits implied by zero-RDD triplets.
#' \item Merge compatible split definitions inferred from different triplets.
#' \item Remove duplicate or equivalent splits.
#' \item Add end-bound splits implied by taxa pairs that never participate
#'       in a detected internal split.
#' }
#'
#' Progress bars are displayed during the major processing stages.
#'
#' @noRd
#'
zero_rdd_to_splits <- function(array_zero_ID, P){

  increment <- 2 * P

  max_split <- increment

  split_set_indicator <- matrix(nrow = max_split, ncol = P)
  ## i_split_set_indicator <- 0
  ## max no of splits = P-1 (i.e. for a binary tree), but we may not identify some of them as multiple
  ## so, maximum number of splits we could identify would be maximum assigned, i.e. H
  ## max no of taxa taking part in a split = P (i.e. for the first split)
  ## if a taxon is not in a split, then it has a default value zero in that split
  ## that is why the array is being populated by zeros.


  have_split <- matrix(nrow = P, ncol = P, 0)
  ## this matrix will keep track whether two taxa have split or not

  n_split_set <- 0
  # This will be incremented as we find new splits

  pb <- txtProgressBar(min = 0, max = P, style = 3)

  for (i in 1:P){
    setTxtProgressBar(pb, i)
    for (j in 1:P){
      for (k in 1:P)
      {
        if (((i != j) &&
             (i != k)) && ((j < k) && (array_zero_ID[i, j, k] == 0)))
        {
          i_split <- 0
          repeat {
            ## i_split is counting the number of already-identified splits compared
            i_split <- i_split + 1


            if (i_split > n_split_set)
            {
              n_split_set <- i_split

              if (n_split_set > max_split)
              {
                max_split <- max_split + increment
                new_array <- matrix(nrow = max_split, ncol = P)
                new_array[(1:(max_split - increment)), ] <- split_set_indicator
                split_set_indicator <- new_array
                rm(new_array)
              }

              this_split_set_indicator <- integer(P)
              this_split_set_indicator[i] <- 1
              this_split_set_indicator[j] <- 2
              this_split_set_indicator[k] <- 2

              split_set_indicator[n_split_set, ] <- this_split_set_indicator

              break
            }
            ## that is, if no other equivalent split is found, define (i,j,k) as a new split; then get out

            if ((split_set_indicator[i_split, i] > 0) &&
                ((
                  (split_set_indicator[i_split, j] > 0)  &&
                  (split_set_indicator[i_split, j] !=  split_set_indicator[i_split, i])
                )))
            {
              split_set_indicator[i_split, k] <- split_set_indicator[i_split, j]

              break
            }
            ## that is, if i is part of split i_split, and j is also part but in a different group than i
            ##(i.e. the split for the (i,j,k) id value is equivalent to split i_split)
            ## then k is part of the j's group
            ## then, get out of the loop
            else if ((split_set_indicator[i_split, i] > 0) &&
                     ((
                       (split_set_indicator[i_split, k] > 0)  &&
                       (split_set_indicator[i_split, k] !=  split_set_indicator[i_split, i])
                     )))
            {
              split_set_indicator[i_split, j] <- split_set_indicator[i_split, k]

              break
            }
            ## that is, if i is part of split i_split, and k is also part but in a different group than i
            ##(i.e. the split for the (i,j,k) id value is equivalent to split i_split)
            ## then j is part of the k's group
            ## then, get out of the loop

          }
        }

        if (((i != j) &&
             (i != k)) && ((j < k) && (array_zero_ID[i, j, k] == 0)))
        {
          have_split[i, j] <- 1
          have_split[j, i] <- 1
          have_split[i, k] <- 1
          have_split[k, i] <- 1
        }
      }
    }
  }


  split_set_indicator <- matrix(ncol = P, split_set_indicator[(1:n_split_set), ])

  ## so far we have all non-end-bound splits; but, there might be repetitions

  increment <- P - 1
  max_split <- increment


  next_split_set_indicator <- matrix(nrow = max_split, ncol = P)


  next_n_split_set <- 0


  merge_record <- NULL

  pb <- txtProgressBar(min = 0, max = n_split_set, style = 3)

  for (i_split in 1:n_split_set)
  {
    j_split <- 0
    setTxtProgressBar(pb, i_split)
    repeat {

      j_split <- j_split + 1

      if (j_split > next_n_split_set)
      {
        next_n_split_set <- next_n_split_set + 1


        if (next_n_split_set > max_split)
        {
          max_split <- max_split + increment
          next_array <- matrix(nrow = max_split, ncol = P)
          next_array[(1:(max_split - increment)), ] <- next_split_set_indicator
          next_split_set_indicator <- next_array
          rm(next_array)
        }



        next_split_set_indicator[next_n_split_set, ] <- split_set_indicator[i_split, ]

        break
      }

      if (((length(intersect(
        which(next_split_set_indicator[j_split, ] == 1),
        which(split_set_indicator[i_split, ] == 1)
      )) > 0) &&
      (length(intersect(
        which(next_split_set_indicator[j_split, ] == 2),
        which(split_set_indicator[i_split, ] == 2)
      )) > 0))
      &&
      ((length(intersect(
        which(next_split_set_indicator[j_split, ] == 1),
        which(split_set_indicator[i_split, ] == 2)
      )) == 0) &&
      (length(intersect(
        which(next_split_set_indicator[j_split, ] == 2),
        which(split_set_indicator[i_split, ] == 1)
      )) == 0)))
      {

        next_split_set_indicator[j_split, (union(
          which(next_split_set_indicator[j_split, ] == 1),
          which(split_set_indicator[i_split, ] == 1)
        ))] <- 1
        next_split_set_indicator[j_split, (union(
          which(next_split_set_indicator[j_split, ] == 2),
          which(split_set_indicator[i_split, ] == 2)
        ))] <- 2

        break
      }
      else if (((length(intersect(
        which(next_split_set_indicator[j_split, ] == 1),
        which(split_set_indicator[i_split, ] == 2)
      )) > 0) &&
      (length(intersect(
        which(next_split_set_indicator[j_split, ] == 2),
        which(split_set_indicator[i_split, ] == 1)
      )) > 0))
      &&
      ((length(intersect(
        which(next_split_set_indicator[j_split, ] == 1),
        which(split_set_indicator[i_split, ] == 1)
      )) == 0) &&
      (length(intersect(
        which(next_split_set_indicator[j_split, ] == 2),
        which(split_set_indicator[i_split, ] == 2)
      )) == 0)))
      {

        next_split_set_indicator[j_split, (union(
          which(next_split_set_indicator[j_split, ] == 1),
          which(split_set_indicator[i_split, ] == 2)
        ))] <- 1
        next_split_set_indicator[j_split, (union(
          which(next_split_set_indicator[j_split, ] == 2),
          which(split_set_indicator[i_split, ] == 1)
        ))] <- 2

        break
      }
    }
  }
  close(pb)

  new_split_set_indicator <- matrix(nrow = max_split, ncol = P)


  new_n_split_set <- 0

  pb <- txtProgressBar(min = 0, max = next_n_split_set, style = 3)
  for (i_split in 1:next_n_split_set)
  {
    setTxtProgressBar(pb, i_split)

    j_split <- 0

    repeat {
      j_split <- j_split + 1

      if (j_split > new_n_split_set)
      {
        new_n_split_set <- new_n_split_set + 1


        if (new_n_split_set > max_split)
        {
          max_split <- max_split + increment
          new_array <- matrix(nrow = max_split, ncol = P)
          new_array[(1:(max_split - increment)), ] <- new_split_set_indicator
          new_split_set_indicator <- new_array
          rm(new_array)
        }

        new_split_set_indicator[new_n_split_set, ] <- next_split_set_indicator[i_split, ]

        break
      }

      if (setequal(which(new_split_set_indicator[j_split, ] == 1),
                   which(next_split_set_indicator[i_split, ] == 1)))
      {
        new_split_set_indicator[j_split, (union(
          which(new_split_set_indicator[j_split, ] == 2),
          which(next_split_set_indicator[i_split, ] == 2)
        ))] <- 2

        break
      }
      else if (setequal(which(new_split_set_indicator[j_split, ] == 2),
                        which(next_split_set_indicator[i_split, ] == 2)))
      {
        new_split_set_indicator[j_split, (union(
          which(new_split_set_indicator[j_split, ] == 1),
          which(next_split_set_indicator[i_split, ] == 1)
        ))] <- 1

        break
      }
      else if (setequal(which(new_split_set_indicator[j_split, ] == 1),
                        which(next_split_set_indicator[i_split, ] == 2)))
      {
        new_split_set_indicator[j_split, (union(
          which(new_split_set_indicator[j_split, ] == 2),
          which(next_split_set_indicator[i_split, ] == 1)
        ))] <- 2

        break
      }
      else if (setequal(which(new_split_set_indicator[j_split, ] == 2),
                        which(next_split_set_indicator[i_split, ] == 1)))
      {
        new_split_set_indicator[j_split, (union(
          which(new_split_set_indicator[j_split, ] == 1),
          which(next_split_set_indicator[i_split, ] == 2)
        ))] <- 1

        break
      }

    }
  }
  close(pb)

  ## done with elimination

  ## so far we have figured out all the non-end-bound splits; next we need to figure out all the end-bound splits using have_split

  for (j in 1:P)
    for (k in 1:P)
      if ((j < k) && (have_split[j, k] == 0))
      {
        new_n_split_set <- new_n_split_set + 1



        if (new_n_split_set > max_split)
        {
          max_split <- max_split + increment
          new_array <- matrix(nrow = max_split, ncol = P)
          new_array[(1:(max_split - increment)), ] <- new_split_set_indicator
          new_split_set_indicator <- new_array
          rm(new_array)
        }


        this_split_set_indicator <- integer(P)
        this_split_set_indicator[j] <- 1
        this_split_set_indicator[k] <- 2

        new_split_set_indicator[new_n_split_set, ] <- this_split_set_indicator
      }

  new_split_set_indicator <- new_split_set_indicator[(1:new_n_split_set), ]
  return(new_split_set_indicator)

}

