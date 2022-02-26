#' Logistic regression analysis
#'
#' A Function to perform logistic regression analysis given a case/control phenotype using PLINK 1.9.
#'
#' @inheritParams PV
#' @param genotype [character] The prefix of PLINK binary files (bed/bim/fam).
#' @param phenotype [character] To specify an alternate phenotype for the logistic regression analysis using the "\code{--pheno}" flag in plink. This file must have a header row. The first and second columns of the phenotype file must be "FID" and "IID", the case/control phenotype in column 3 (1 = unaffected (control), 2 = affected (case)), and covariates in remaining columns. If \code{phenotype} is NULL, the original fam file must contain a phenotype in column 6. See the PLINK 1.9 documentation for details (\url{https://www.cog-genomics.org/plink/1.9/}).
#' @param covar.number [vector] To specify a subset of column numbers of covariates to load from \code{phenotype} file using the "\code{--covar-number}" flag in plink. If NULL (the default), the logistic regression model without covariate adjustment.
#' @return \code{plink.lr} returns a data frame with the results of logistic regression analysis derived from PLINK.
#' @export
#' @examples
#' input.dir <- system.file("extdata", package="pv")
#' output.dir <- system.file("extdata", package="pv")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' lr.result <- plink.lr(input.dir = input.dir,
#' output.dir = output.dir,
#' genotype = "train",
#' phenotype = "train.phenotypes.txt",
#' covar.number = c(2, 3),
#' plink.path = path2plink,
#' verbose = TRUE)
#' }
plink.lr <- function(input.dir, output.dir, genotype, phenotype = NULL, covar.number = NULL, plink.path = NULL, verbose = TRUE){
  if(missing(input.dir)){
    input.dir <- getwd()
  }
  if(missing(output.dir)){
    output.dir <- getwd()
  }

  genotype.path <- file.path(input.dir, genotype)
  genotype.path <- gsub('\\\\', '/', genotype.path)

  output.path <- file.path(output.dir, genotype)
  output.path <- gsub('\\\\', '/', output.path)

  phenotype.path <- file.path(input.dir, phenotype)### updated
  phenotype.path <- gsub('\\\\', '/', phenotype.path)

  # Check PLINK software access
  plink.path <- plinkQC::checkPlink(plink.path)

  # Check data available
  plinkQC:::checkFormat(genotype.path)

  if(is.null(phenotype)){
    sys::exec_wait(plink.path,
                   args = c("--noweb", "--bfile", genotype.path, "--logistic", "beta", "--allow-no-sex", "--out", output.path),
                   std_out = verbose,
                   std_err = verbose)
  }else{
    if(is.null(covar.number)){
      sys::exec_wait(plink.path,
                     args = c("--noweb", "--bfile", genotype.path, "--logistic", "beta", "--pheno", phenotype.path, "--mpheno", 1, "--allow-no-sex", "--out", output.path),
                     std_out = verbose,
                     std_err = verbose)
    }else{
      covar.number.range <- length(read.table(phenotype.path, nrows = 1)) - 2
      if(any(covar.number > covar.number.range) | any(covar.number < 2)){
        stop("The covariant number ranges from 2 to ",covar.number.range, ", but not all specified covariant numbers [", paste(covar.number, collapse = ",") ,"] are in this range.")
      }
      sys::exec_wait(plink.path,
                     args = c("--noweb", "--bfile", genotype.path, "--logistic", "beta", "--pheno", phenotype.path, "--mpheno", 1, "--covar", phenotype.path, "--covar-number", paste(covar.number, collapse = ","), "--allow-no-sex", "--out", output.path),
                     std_out = verbose,
                     std_err = verbose)
    }
  }

  assoc.logstic <- data.table::fread(paste0(output.path, ".assoc.logistic"), data.table = FALSE)
  return(assoc.logstic)
}
