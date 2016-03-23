
compute_expectations <- function(counts_mat, 
                                 by = c("all","annotation"), 
                                 norm = FALSE,
                                 annotation = NULL){
  
  by = match.arg(by)
  
  counts_info = counts_summary(counts_mat)
  
  if (by == "all"){
    if (norm){
      expectation = rowSums(counts_mat / 
                              matrix(counts_info$fragments_per_sample, 
                                     nrow = counts_info$npeak, 
                                     ncol = counts_info$nsample,
                                     byrow = TRUE))
    } else{
      expectation = counts_info$fragments_per_peak / counts_info$total_fragments
    } 
  } else if (by == "annotation"){
    anno = as.factor(annotation)
    n_anno = length(levels(anno))
    mat = matrix(nrow = counts_info$npeak, ncol = n_anno)
    if (norm){
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(counts_mat[,ix] / 
                            matrix(counts_info$fragments_per_sample[ix], 
                                   nrow = counts_info$npeak, 
                                   ncol = length(ix),
                                   byrow = TRUE))
      }
    } else{
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(counts_mat[,ix]) / sum(counts_mat[,ix])
      }
    }
    expectation = rowMeans(mat)
  } 
  return(expectation)
}



#' compute_deviations
#' 
#' Computes deviations across sets of annotations
#' @param peak_indices list of indices representing different sets of peaks
#' @param counts_mat matrix with counts
#' @param background_peaks background peaks matrix
#' @param expectation expectations computed using \code{\link{compute_expectations}}
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  matrix with deviation scores
#' @seealso  \code{\link{compute_variability}}, \code{\link{plot_variability}}
#' @export
compute_deviations <- function(counts_mat, 
                               background_peaks,
                               peak_indices = NULL, 
                               expectation = NULL,
                               norm = TRUE){
  
  stopifnot(inherits(counts_mat,"Matrix") || inherits(counts_mat,"matrix"))
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  
  counts_info <- counts_summary(counts_mat)
      
  if (is.null(peak_indices)){
    peak_indices <- lapply(1:counts_info$npeak, function(x) x)
  } else if (!is.list(peak_indices) && is.vector(peak_indices)){
    peak_indices = list(peak_indices)
  }
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  } else{
    stopifnot(length(expectation) != nrow(counts_mat))
  }
  
  stopifnot(inherits(peak_indices,"list"))
  

  # check that indices fall within appropriate bounds
  tmp <- unlist(peak_indices, use.names =F)
  if (is.null(tmp) ||
        !(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_info$npeak ||
        min(tmp) < 1){
    stop("peak_indices are not valid")
  }

  if(is.null(names(peak_indices))){
    names(peak_indices) = 1:length(peak_indices)
  }
  # remove sets of length 0
  peak_indices <- peak_indices[which(sapply(peak_indices,length)>0)]
  
  results <- BiocParallel::bplapply(peak_indices,
                                      compute_deviations_single_wrapper,
                                      counts_mat, 
                                      background_peaks,
                                      expectation,
                                      counts_info,
                                      norm = norm)
  
  out <- list()
  out$z <- t(vapply(results, function(x) x[["z"]], rep(0,counts_info$nsample))) 
  out$fc <- t(vapply(results, function(x) x[["fc"]], rep(0,counts_info$nsample))) 
  
  return(out)
}

# Helper Functions -------------------------------------------------------------

compute_deviations_single_wrapper <- function(peak_set, 
                                              counts_mat, 
                                              background_peaks,
                                              expectation,
                                              counts_info = NULL,
                                              intermediate_results=FALSE,
                                              norm = TRUE){
  if (is.null(counts_info)){
    counts_info = counts_summary(counts_mat)
  }
  ### Determine if any sets have too low expected counts
  expected_totals <- sum(expectation[peak_set]) * counts_info$fragments_per_sample
  fail_filter <- which(expected_totals < 5)  
  ###
  if (inherits(counts_mat,"dgCMatrix")){
    if (length(peak_set) == 1){
      out <- compute_deviations_single_peak_sparse(peak_set[1] - 1, counts_mat, background_peaks, expectation, counts_info, intermediate_results, norm)  
    } else{
      out <- compute_deviations_single_sparse(peak_set - 1, counts_mat, background_peaks, expectation, counts_info, intermediate_results, norm)        
    }
  } else{
    out <- compute_deviations_single_dense(peak_set - 1, counts_mat, background_peaks, expectation, counts_info, intermediate_results, norm)       
  }
  sample_names = colnames(counts_mat)
  if (length(fail_filter) > 0){
      out[["z"]][fail_filter] <- NA  
    }  
  names(out[["z"]]) = sample_names
  names(out[["fc"]]) = sample_names
  if (intermediate_results){
    names(out[["observed_deviations"]]) = sample_names
    colnames(out[["sampled_deviations"]]) = sample_names
    names(out[["observed"]]) = sample_names
    colnames(out[["sampled"]]) = sample_names
  } 
  return(out)
}


compute_deviations_single_slow <- function(peak_set, 
                                counts_mat, 
                                background_peaks,
                                expectation,
                                counts_info = NULL,
                                intermediate_results = FALSE){
  
  if (is.null(counts_info)){
    counts_info = counts_summary(counts_mat)
  }
  
  tf_count <- length(peak_set)
  
  if (tf_count == 1){
    observed <- (as.vector(counts_mat[peak_set,]) - (expectation[peak_set] * counts_info$fragments_per_sample)) /
      sqrt(expectation[peak_set] * counts_info$fragments_per_sample)
    names(observed) <- colnames(counts_mat)
    sampled_counts <-  t((as.matrix(counts_mat[background_peaks[peak_set,],]) -
                            outer(expectation[background_peaks[peak_set,]],counts_info$fragments_per_sample)) /    
                           outer(sqrt(expectation[background_peaks[peak_set,]]),
                                 sqrt(counts_info$fragments_per_sample)))
  }
  else {
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                           dims = c(1, counts_info$npeak))
    
    observed <- as.vector(tf_vec %*% ((counts_mat - outer(expectation, counts_info$fragments_per_sample)) / 
                                        outer(sqrt(expectation), sqrt(counts_info$fragments_per_sample))))
    names(observed) <- colnames(counts_mat)
    
    sampled_counts <- sample_background_peaks(counts_mat = counts_mat,
                                              background_peaks = background_peaks,
                                              peak_set = peak_set,
                                              counts_info = counts_info,
                                              expectation = expectation)
  } 
  
  mean_sampled_deviation <- colMeans(sampled_counts)
  sd_sampled_deviation <- apply(sampled_counts, 2, sd)
  
  res <- (observed - mean_sampled_deviation) / sd_sampled_deviation
  
  if (intermediate_results){
    out = list(deviations = res, observed = observed,
               sampled = sampled_counts)
    return(out)
  } else{
    return(res)
  }
}

sample_background_peaks <- function(counts_mat, peak_set, background_peaks, counts_info, expectation){
  
  niterations = ncol(background_peaks)
  sample_mat = sparseMatrix(j = as.vector(background_peaks[peak_set,1:niterations]), 
                            i = rep(1:niterations, each = length(peak_set)), 
                            x=1, 
                            dims = c(ncol(background_peaks), counts_info$npeak))
  
  sampled_counts =  as.matrix(sample_mat %*% ((counts_mat - outer(expectation, counts_info$fragments_per_sample)) / 
                                                  outer(sqrt(expectation), sqrt(counts_info$fragments_per_sample))))
  
  return(sampled_counts)
}