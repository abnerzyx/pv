#' Feature selection
#'
#' Logistic regression-based feature selection approach.
#'
#' @inheritParams PV
#' @param genotype [character] The prefix of PLINK binary files (bed/bim/fam).
#' @param phenotype [character] A space- or tab-delimited file to specify an alternate phenotype for the logistic regression analysis using the "\code{--pheno}" flag in plink. This file must have a header row. The first and second columns of the phenotype file must be "FID" and "IID", the case/control phenotype in column 3 (1 = unaffected (control), 2 = affected (case)), and covariates in remaining columns.  See the PLINK 1.9 documentation for details (\url{https://www.cog-genomics.org/plink/1.9/}).
#' @return \code{feature.selection} return a list containing the results of logistic regression analysis derived from PLINK (via \code{\link{plink.lr}}), the indices and names of selected features.
#' \item{lr.result}{The output of \code{\link{plink.lr}}}
#' \item{index}{A vector of indices of the selected features.}
#' \item{name}{A vector of names of the selected features.}
#' @seealso {\code{\link{plink.lr}}}
#' @export
#' @examples
#' input.dir <- system.file("extdata", package="pv")
#' output.dir <- system.file("extdata", package="pv")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' feature.selection.result <- feature.selection(input.dir = input.dir,
#' output.dir = output.dir,
#' genotype = "train",
#' phenotype = "train.phenotypes.txt",
#' covar.number = c(2, 3),
#' plink.path = path2plink,
#' topK = 10,
#' verbose = TRUE)
#' }
feature.selection <- function(input.dir, output.dir, genotype, phenotype, covar.number = NULL, plink.path = NULL, topK = 10, P.value = NULL, candidate.SNPs = NULL, verbose = TRUE){
  if(missing(input.dir)){
    input.dir <- getwd()
  }
  if(missing(output.dir)){
    output.dir <- getwd()
  }

  # run logistic regression analysis
  lr.result <- plink.lr(input.dir, output.dir, genotype, phenotype, covar.number, plink.path, verbose)

  if(is.null(candidate.SNPs)){
    p.val <- lr.result[which(lr.result$TEST == "ADD"), ]$P
    p.val[is.na(p.val)] <- 1

    if(is.null(P.value)){
      fea.index <- order(p.val)[1:topK]
    }else{
      fea.index <- which(p.val < P.value)
      if(length(fea.index) <1){
        stop("No SNPs pass this P-value threshold: ", P.value, ".")
      }
    }
    fea.name <- lr.result$SNP[which(lr.result$TEST == "ADD")[fea.index]]
  }else{
    genotype.path <- file.path(input.dir, genotype)
    genotype.path <- gsub('\\\\', '/', genotype.path)
    # Check data available
    plinkQC:::checkFormat(genotype.path)

    bim <- data.table::fread(paste0(genotype.path, ".bim"), data.table = FALSE)
    fea.index <- match(candidate.SNPs, bim$V2)
    if(any(is.na(fea.index))){
      stop("The specified candidate SNP (", paste(candidate.SNPs[is.na(fea.index)], sep="",collapse = ", "), ") not found in the input genotype data!")
    }
    fea.name <- candidate.SNPs
  }
  return(list(lr.result = lr.result,
              index = fea.index,
              name = fea.name))
}
