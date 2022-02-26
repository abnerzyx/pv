#' Prism Vote
#'
#' Perform Prism Vote method to build a prediction model on the training dataset, and make predictions on the test dataset.
#'
#' @param input.dir [character] The full absolute path to the directory containing the training and test dataset. If \code{input.dir} is missing, the current working directory obtained by \code{getwd()} is used.
#' @param output.dir [character] The full absolute path where the result will be written to. If \code{output.dir} is missing, the current working directory obtained by \code{getwd()} is used.
#' @param train.genotype [character] The prefix of PLINK binary files (bed/bim/fam) of the training dataset.
#' @param train.phenotype [character] A space- or tab-delimited file to specify an alternate phenotype of the training dataset for the logistic regression analysis using the "\code{--pheno}" flag in plink. This file must have a header row. The first and second columns of the phenotype file must be "FID" and "IID", the case/control phenotype in column 3 (1 = unaffected (control), 2 = affected (case)), and covariates in remaining columns. See the PLINK 1.9 documentation for details (\url{https://www.cog-genomics.org/plink/1.9/}).
#' @param test.genotype [character] The prefix of PLINK binary files (bed/bim/fam) of the test dataset.
#' @param test.phenotype [character] A space- or tab-delimited file to specify an alternate phenotype of the training dataset for the logistic regression analysis using the "\code{--pheno}" flag in plink. This file must have a header row. The first and second columns of the phenotype file must be "FID" and "IID", the case/control phenotype in column 3 (1 = unaffected (control), 2 = affected (case)), and covariates in remaining columns.  See the PLINK 1.9 documentation for details (\url{https://www.cog-genomics.org/plink/1.9/}).
#' @param covar.number.PV [vector] Used in PV model to specify a subset of column numbers of covariates to load from \code{phenotype} file using the "\code{--covar-number}" flag in plink (via \code{\link{plink.lr}}). If NULL (the default), the logistic regression model without covariate adjustment.
#' @param covar.number.LR [vector] Used in LR and PRS model to specify a subset of column numbers of covariates to load from \code{phenotype} file using the "\code{--covar-number}" flag in plink (via \code{\link{plink.lr}}). If NULL (the default), the logistic regression model without covariate adjustment.
#' @param PCA.separate [logical] If TURE, the principal components are calculated from the training dataset and then project the test dataset onto those principal components. If FALSE, the principal components are calculated from the combined data of the training and test dataset. The default value is FALSE.
#' @param PCs.count [numeric] To specify the number of top principal components that should be extracted. The default value is 10.
#' @param stratum.count [numeric] To specify the number of strata, as default the sample size of each stratum is N/stratum.count, N is the sample size. The default value is 2.
#' @param plink.path [character] The full absolute path to the PLINK executable file. The executable to run is path/to/plink.exe if you are on a Windows operating system, for Unix-like operating system this is path/to/plink. If \code{plink.path} is NULL, the PLINK PATH should be added as a system environment variable.
#' @param P.value [double] To specify the genome-wide significance P-value threshold to select the significant SNPs to build a prediction model. The default value is NULL. This value is ignored when candidate.SNPs is not NULL. When left NULL (the default), the topK or candidate.SNPs will be used. The P-value of each SNP is calculated from logistic regression analysis using PLINK 1.9 (via \code{\link{plink.lr}}).
#' @param topK [numeric] To specify the top K significant SNPs to build a prediction model. For a fair comparison, the number of the top-ranked SNPs from entire sample (for LR and PRS model) equals to the number of the unique union set of the selected SNPs from each stratum in PV. The default value is 10. This value is ignored when P.value or candidate.SNPs is not NULL.
#' @param candidate.SNPs [vector] A character vector of SNP name, used to specify the candidate SNPs to build a prediction model, ignores \code{P.value} and \code{topK}. The default value is NULL. Should match the names of SNPs in the provided PLINK binary files.
#' @param CS [logical] If TRUE, the softmax of cosine similarity will be used to calculate the probability that the samples belong to each stratum. If FALSE, the squared distance of a subject to a cluster center empirically follows a chi-squared distribution will be used. The default value is FALSE.
#' @param verbose [logical] If TRUE, the PLINK log, error, and warning information are printed to standard out. The default value is TRUE.
#' @return \code{PV} returns a list with the following components:
#' \item{stratification.result}{The output of \code{\link{stratification}}}
#' \item{feature.selection.result}{The output of \code{\link{feature.selection}}}
#' \item{LR.result}{The results of logistic regression model. A list with i) predict, the output of \code{\link{LR.model}}. ii) performance, a list containing the AUC and accuracy (Acc) value.}
#' \item{LR.PV.result}{The results of the logistic regression model under the Prism Vote framework. A list with i) predict, the output of \code{\link{LR.PV.model}}. ii) performance, a list containing the AUC and accuracy (Acc) value.}
#' \item{PRS.result}{The results of the logistic regression model based on the polygenic risk score. A list with i).predict, the output of \code{\link{LR.model}}. ii) performance, a list containing the AUC and accuracy (Acc) value.}
#' \item{PRS.PV.result}{The results of the logistic regression model based on the polygenic risk score under the Prism Vote framework. A list with i) predict, the output of \code{\link{LR.PV.model}}. ii) performance, a list containing the AUC and accuracy (Acc) value.}
#' @seealso {\code{\link{PCA}}, \code{\link{stratification}}, \code{\link{LR.model}}, \code{\link{PV.model}}, \code{\link{PRS}}, \code{\link{feature.selection}}}
#' @export
#' @examples
#' input.dir <- system.file("extdata", package="pv")
#' output.dir <- system.file("extdata", package="pv")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' pv.result <- PV(input.dir = input.dir,
#' output.dir = output.dir,
#' train.genotype = "train",
#' train.phenotype = "train.phenotypes.txt",
#' test.genotype = "test",
#' test.phenotype = "test.phenotypes.txt",
#' covar.number.PV = c(2,3),
#' covar.number.LR = c(2,3),
#' PCA.separate = FALSE,
#' PCs.count = 10,
#' stratum.count = 2,
#' plink.path = path2plink,
#' P.value = NULL,
#' topK = 10,
#' CS = FALSE,
#' candidate.SNPs = NULL,
#' verbose = TRUE)
#' }
PV <- function(input.dir, output.dir, train.genotype, train.phenotype, test.genotype, test.phenotype, covar.number.PV = NULL, covar.number.LR = NULL, PCA.separate = FALSE, PCs.count = 10, stratum.count = 2, plink.path = NULL, P.value = NULL, topK = 10, candidate.SNPs = NULL, CS = FALSE, verbose = TRUE){
  if(missing(input.dir)){
    input.dir <- getwd()
  }
  if(missing(output.dir)){
    output.dir <- getwd()
  }
  train.genotype.path <- file.path(input.dir, train.genotype)
  train.genotype.path <- gsub('\\\\', '/', train.genotype.path)

  train.phenotype.path <- file.path(input.dir, train.phenotype)
  train.phenotype.path <- gsub('\\\\', '/', train.phenotype.path)

  train.output.path <- file.path(output.dir, train.genotype)
  train.output.path <- gsub('\\\\', '/', train.output.path)

  test.genotype.path <- file.path(input.dir, test.genotype)
  test.genotype.path <- gsub('\\\\', '/', test.genotype.path)

  test.phenotype.path <- file.path(input.dir, test.phenotype)
  test.phenotype.path <- gsub('\\\\', '/', test.phenotype.path)

  test.output.path <- file.path(output.dir, test.genotype)
  test.output.path <- gsub('\\\\', '/', test.output.path)

  # check sample missing
  train.missing.flag <- checkMissingness(input.dir,
                                         output.dir,
                                         train.genotype,
                                         plink.path,
                                         verbose)

  test.missing.flag <- checkMissingness(input.dir,
                                        output.dir,
                                        test.genotype,
                                        plink.path,
                                        verbose)

  # stratification
  stratification.result <- stratification(input.dir,
                                          output.dir,
                                          train.genotype,
                                          test.genotype,
                                          stratum.count,
                                          PCA.separate,
                                          PCs.count,
                                          plink.path,
                                          CS,
                                          verbose)

  # generate the genotype and phenotype data of each stratum.
  data.split(input.dir,
             output.dir,
             train.genotype,
             train.phenotype,
             stratification.result$train.stratum,
             plink.path,
             verbose)

  # Feature selection
  selected.SNPs <- c()
  # 1. LR.PV feature selection
  feature.selection.result.PV <- list()
  for (i in 1:stratum.count) {
    feature.selection.result.PV[[i]] <- feature.selection(output.dir,
                                                       output.dir,
                                                       paste0(train.genotype, ".stratum.",i),
                                                       paste0(train.genotype, ".stratum.",i,".phenotype.txt"),
                                                       covar.number.PV,
                                                       plink.path,
                                                       topK,
                                                       P.value,
                                                       candidate.SNPs,
                                                       verbose)

    # updated
    write.table(feature.selection.result.PV[[i]]$name, file = paste0(train.output.path, ".PV_selected_SNPs_stratum.", i, ".txt"), col.names = F,row.names = F,sep = "\n",quote = F)

    selected.SNPs <- c(selected.SNPs, feature.selection.result.PV[[i]]$name)

    # convert the selected genotype data to bfile(bed/bim/fam) and 0/1/2 (.raw).
    # if(external){
    #   # external model
    #   sys::exec_wait(plink.path.temp,
    #                  args = c("--noweb", "--bfile", train.genotype.path,  "--extract", paste0(train.output.path, ".PV_selected_SNPs_stratum.", i, ".txt"), "--make-bed", "--recodeA", "--out", paste0(train.output.path, ".selected.stratum.",i)),
    #                  std_out = verbose,
    #                  std_err = verbose)
    # }

      # separate model
    plink.path.temp <- plinkQC::checkPlink(plink.path)
    sys::exec_wait(plink.path.temp,
                   args = c("--noweb", "--bfile", train.genotype.path, "--keep", paste0(train.output.path, ".stratum.",i,".sampleID.txt"), "--extract", paste0(train.output.path,".PV_selected_SNPs_stratum.", i, ".txt"), "--make-bed", "--recodeA", "--out", paste0(train.output.path, ".selected.stratum.",i)),
                   std_out = verbose,
                   std_err = verbose)

  }

  # 2. LR feature selection
  PV.selected.SNPs <- unique(selected.SNPs)
  #the number of the top-ranked SNPs from entire sample equals to the number of the unique union set of the selected SNPs from each stratum in PV.
  feature.selection.result.LR <- feature.selection(input.dir,
                                                   output.dir,
                                                   train.genotype,
                                                   train.phenotype,
                                                   covar.number.LR,
                                                   plink.path,
                                                   length(PV.selected.SNPs),
                                                   P.value,
                                                   candidate.SNPs,
                                                   verbose)

  write.table(feature.selection.result.LR$name, file = paste0(train.output.path, ".LR_selected_SNPs.txt"), col.names = F,row.names = F,sep = "\n",quote = F)

  selected.SNPs <- c(selected.SNPs, feature.selection.result.LR$name)
  write.table(selected.SNPs, file = paste0(train.output.path, ".LR_and_PV_selected_SNPs.txt"), col.names = F,row.names = F,sep = "\n",quote = F)

  sys::exec_wait(plink.path.temp,
                 args = c("--noweb", "--bfile", train.genotype.path, "--extract", paste0(train.output.path, ".LR_selected_SNPs.txt"), "--make-bed", "--recodeA", "--out", paste0(train.output.path, ".selected")),
                 std_out = verbose,
                 std_err = verbose)

  sys::exec_wait(plink.path.temp,
                 args = c("--noweb", "--bfile", test.genotype.path, "--extract", paste0(train.output.path, ".LR_and_PV_selected_SNPs.txt"), "--make-bed", "--recodeA", "--out", paste0(test.output.path, ".selected")),
                 std_out = verbose,
                 std_err = verbose)

  # Prediction
  # 1. LR model
  train.data.LR <- data.table::fread(paste0(train.output.path, ".selected.raw"), data.table = FALSE)[, -c(1,2,3,4,5,6), drop = FALSE]
  train.data.bim <- data.table::fread(paste0(train.output.path, ".selected.bim"), data.table = F)
  colnames(train.data.LR) <- train.data.bim$V2
  if(train.missing.flag){
    train.data.LR <- data.frame(apply(train.data.LR, 2, impute), check.names = FALSE)
  }
  train.pheno <- data.table::fread(train.phenotype.path, data.table = FALSE)
  rownames(train.pheno) <- train.pheno[, 2]
  train.data.LR$Y <- train.pheno[, 3]-1

  test.data <- data.table::fread(paste0(test.output.path, ".selected.raw"), data.table = FALSE)[, -c(1,2,3,4,5,6)]
  test.data.bim <- data.table::fread(paste0(test.output.path, ".selected.bim"), data.table = F)
  colnames(test.data) <- test.data.bim$V2
  if(test.missing.flag){
    test.data <- data.frame(apply(test.data, 2, impute), check.names = FALSE)
  }
  test.pheno <- data.table::fread(test.phenotype.path, data.table = FALSE)
  test.data.LR <- test.data[, feature.selection.result.LR$name, drop = FALSE]
  test.data.LR$Y <- test.pheno[,3]-1

  if(!is.null(covar.number.LR)){
    train.data.LR <- cbind(train.data.LR, train.pheno[, covar.number.LR + 2, drop = FALSE])
    test.data.LR <- cbind(test.data.LR, test.pheno[, covar.number.LR + 2, drop = FALSE])
  }

  LR.pred <- LR.model(train.data.LR, test.data.LR)
  LR.performance <- model.evaluation(test.data.LR$Y, LR.pred)

  # 2. LR + PV
  train.data.list <- list()
  test.data.list <- list()
  for (i in 1:stratum.count) {
    train.data.stratum <- data.table::fread(paste0(train.output.path, ".selected.stratum.", i, ".raw"), data.table = FALSE)[, -c(1,2,3,4,5,6), drop = FALSE]
    train.data.bim <- data.table::fread(paste0(train.output.path, ".selected.stratum.", i, ".bim"), data.table = F)
    colnames(train.data.stratum) <- train.data.bim$V2

    if(train.missing.flag){
      train.data.stratum <- data.frame(apply(train.data.stratum, 2, impute), check.names = FALSE)
    }

    stratum.fid <- data.table::fread(paste0(train.output.path, ".selected.stratum.", i, ".fam"), data.table = FALSE)

    train.data.stratum$Y <- train.pheno[as.character(stratum.fid[, 2]), 3]-1

    test.data.stratum <- test.data[, feature.selection.result.PV[[i]]$name, drop = FALSE]
    test.data.stratum$Y <- test.pheno[, 3]-1

    if(!is.null(covar.number.PV)){
      train.data.stratum <- cbind(train.data.stratum, train.pheno[as.character(stratum.fid[, 2]), covar.number.PV + 2, drop = FALSE])
      test.data.stratum <- cbind(test.data.stratum, test.pheno[, covar.number.PV + 2, drop = FALSE])
    }

    train.data.list[[i]] <- train.data.stratum
    test.data.list[[i]] <- test.data.stratum
  }

  LR.PV.pred <- LR.PV.model(train.data.list, test.data.list, stratification.result$test.prob.to.stratum)
  LR.PV.performance <- model.evaluation(test.data.list[[1]]$Y, LR.PV.pred$pred.agg)

  # 3. PRS
  beta <- feature.selection.result.LR$lr.result$BETA[which(feature.selection.result.LR$lr.result$TEST == "ADD")[feature.selection.result.LR$index]]
  beta[is.na(beta)] <- 0

  train.data.LR <- data.table::fread(paste0(train.output.path, ".selected.raw"), data.table = FALSE)[, -c(1,2,3,4,5,6), drop = FALSE]
  train.data.bim <- data.table::fread(paste0(train.output.path, ".selected.bim"), data.table = F)
  colnames(train.data.LR) <- train.data.bim$V2

  if(train.missing.flag){
    train.data.LR <- data.frame(apply(train.data.LR, 2, impute), check.names = FALSE)
  }

  train.prs <- PRS(train.data.LR, beta)
  train.data.PRS <- data.frame(PRS = train.prs, Y = train.pheno[, 3]-1, check.names = FALSE)

  test.prs <- PRS(test.data[, feature.selection.result.LR$name, drop = FALSE], beta)
  test.data.PRS <- data.frame(PRS = test.prs, Y = test.pheno[,3]-1)

  if(!is.null(covar.number.LR)){
    train.data.PRS <- cbind(train.data.PRS, train.pheno[, covar.number.LR + 2, drop = FALSE])
    test.data.PRS <- cbind(test.data.PRS, test.pheno[, covar.number.LR + 2, drop = FALSE])
  }

  PRS.pred <- LR.model(train.data.PRS, test.data.PRS)
  PRS.performance <- model.evaluation(test.data.PRS$Y, PRS.pred)


  # 4. PRS + PV
  train.data.PRS.list <- list()
  test.data.PRS.list <- list()
  for (i in 1:stratum.count) {
    beta <- feature.selection.result.PV[[i]]$lr.result$BETA[which(feature.selection.result.PV[[i]]$lr.result$TEST == "ADD")[feature.selection.result.PV[[i]]$index]]
    beta[is.na(beta)] <- 0

    train.data.stratum <- data.table::fread(paste0(train.output.path, ".selected.stratum.", i, ".raw"), data.table = FALSE)[, -c(1,2,3,4,5,6), drop = FALSE]
    train.data.bim <- data.table::fread(paste0(train.output.path, ".selected.stratum.", i, ".bim"), data.table = F)
    colnames(train.data.stratum) <- train.data.bim$V2

    if(train.missing.flag){
      train.data.stratum <- data.frame(apply(train.data.stratum, 2, impute), check.names = FALSE)
    }

    train.prs <- PRS(train.data.stratum, beta)

    stratum.fid <- data.table::fread(paste0(train.output.path, ".selected.stratum.", i, ".fam"), data.table = FALSE)

    train.data.stratum.PRS <- data.frame(PRS = train.prs, Y = train.pheno[as.character(stratum.fid[, 2]), 3]-1, check.names = FALSE)

    test.prs <- PRS(test.data[, feature.selection.result.PV[[i]]$name, drop = FALSE], beta)
    test.data.stratum.PRS <- data.frame(PRS = test.prs, Y = test.pheno[, 3]-1, check.names = FALSE)

    if(!is.null(covar.number.PV)){
      train.data.stratum.PRS <- cbind(train.data.stratum.PRS, train.pheno[as.character(stratum.fid[, 2]), covar.number.PV + 2, drop = FALSE])
      test.data.stratum.PRS <- cbind(test.data.stratum.PRS, test.pheno[, covar.number.PV + 2, drop = FALSE])
    }
    train.data.PRS.list[[i]] <- train.data.stratum.PRS
    test.data.PRS.list[[i]] <- test.data.stratum.PRS
  }

  PRS.PV.pred <- LR.PV.model(train.data.PRS.list, test.data.PRS.list, stratification.result$test.prob.to.stratum)
  PRS.PV.performance <- model.evaluation(test.data.PRS.list[[1]]$Y, PRS.PV.pred$pred.agg)


  return(list(stratification.result = stratification.result,
              feature.selection.result = list(feature.selection.result.PV,
                                              feature.selection.result.LR),
              LR.result = list(predict = LR.pred,
                               performance = LR.performance),
              LR.PV.result = list(predict = LR.PV.pred,
                                  performance = LR.PV.performance),
              PRS.result = list(predict = PRS.pred,
                                performance = PRS.performance),
              PRS.PV.result = list(predict = PRS.PV.pred,
                                   performance = PRS.PV.performance)))
}
