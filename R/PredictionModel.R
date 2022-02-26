#' Logistic regression-based Prism Vote model
#'
#' A Function to build a logistic regression prediction model on the training dataset, and make predictions on the test dataset with the framework of Prism Vote.
#'
#' @param train.data.list [list] A list containing the training data [data.frame] of each stratum. The stratum data with rows corresponding to individuals and columns to features, and must contain a column named Y providing the case/control phenotype (0 = unaffected (control), 1 = affected (case)).
#' @param test.data.list [list] A list containing the test data [data.frame]. It must contain a column named Y providing the case/control phenotype (0 = unaffected (control), 1 = affected (case)).
#' @param test.prob.to.stratum [list] Output of \code{\link{stratification}}. A list providing the probability that the test samples belong to each stratum.
#' @return \code{LR.PV.model} return a list containing the predicted probability from each stratum and the aggregated prediction.
#' \item{pred.stratum}{A list of the predicted probability from each stratum.}
#' \item{pred.agg}{A vector of the aggregated prediction.}
#' @export
#' @examples
#' input.dir <- system.file("data", package="pv")
#' output.dir <- system.file("data", package="pv")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' stratum.count <- 2
#' covar.number <-  c(2, 3)
#'
#' stratification.result <- stratification(input.dir = input.dir,
#' output.dir = input.dir,
#' train.genotype = "train",
#' test.genotype = "test",
#' stratum.count = stratum.count,
#' PCA.separate = FALSE,
#' PCs.count = 10,
#' plink.path = path2plink,
#' verbose = TRUE)
#'
#' feature.selection.result <- list()
#' for (i in 1:stratum.count) {
#'   feature.selection.result[[i]] <- feature.selection(input.dir = input.dir,
#'   output.dir = output.dir,
#'   genotype = paste0("train.stratum.",i),
#'   phenotype = paste0("train.stratum.",i,".phenotype.txt"),
#'   covar.number = covar.number,
#'   plink.path = path2plink,
#'   topK = 10,
#'   verbose = TRUE)
#' }
#'
#' train.genotype.path <- file.path(input.dir, "train.raw")
#' train.genotype.path <- gsub('\\\\', '/', train.genotype.path)
#' train.data <- data.table::fread(train.genotype.path, data.table = FALSE)[, -c(1,2,3,4,5,6)]
#' colnames(train.data) <- unlist(purrr::map(colnames(train.data),function(x) {substr(x,1, nchar(x)-2)}))
#' train.pheno <- data.table::fread(train.phenotype.path, data.table = FALSE)
#'
#' test.genotype.path <- file.path(input.dir, "test.raw")
#' test.genotype.path <- gsub('\\\\', '/', test.genotype.path)
#' test.data <- data.table::fread(test.genotype.path, data.table = FALSE)[, -c(1,2,3,4,5,6)]
#' colnames(test.data) <- unlist(purrr::map(colnames(test.data),function(x) {substr(x,1, nchar(x)-2)}))
#' test.pheno <- data.table::fread(test.phenotype.path, data.table = FALSE)
#'
#' train.data.list <- list()
#' test.data.list <- list()
#' for (i in 1:stratum.count) {
#'   train.data.stratum <- train.data[stratification.result$train.stratum.index[[i]],
#'   feature.selection.result[[i]]$index, drop = FALSE]
#'   train.data.stratum$Y <- train.pheno[stratification.result$train.stratum.index[[i]], 3]-1
#'
#'   test.data.stratum <- test.data[, feature.selection.result[[i]]$index, drop = FALSE]
#'   test.data.stratum$Y <- test.pheno[, 3]-1
#
#'   if(!is.null(covar.number){
#'     train.data.stratum <- cbind(train.data.stratum,
#'     train.pheno[stratification.result$train.stratum.index[[i]], covar.number + 2, drop = FALSE])
#'     test.data.stratum <- cbind(test.data.stratum, test.pheno[, covar.number + 2, drop = FALSE])
#'   }
#'
#'   train.data.list[[i]] <- train.data.stratum
#'   test.data.list[[i]] <- test.data.stratum
#' }
#'
#' LR.PV.pred <- LR.PV.model(train.data.list, test.data.list, stratification.result$test.prob.to.stratum)
#' }

LR.PV.model <- function(train.data.list, test.data.list, test.prob.to.stratum){
  pred.stratum <- list()
  pred.agg <- c(0)
  for(i in 1:length(train.data.list)){
    lr.fit <- glm(Y~., data = train.data.list[[i]], family = binomial)
    pred.stratum[[i]] <- predict(lr.fit, test.data.list[[i]], type = "response")
    pred.agg <- pred.agg + pred.stratum[[i]] * test.prob.to.stratum[[i]]
  }
  return(list(pred.stratum = pred.stratum,
              pred.agg = pred.agg))
}

#' Logistic regression model
#'
#' A Function to build a logistic regression prediction model on the training dataset, and make predictions on the test dataset.
#'
#' @param train.data [data.frame] The training dataset. It must contain a column named Y providing the case/control phenotype (0 = unaffected (control), 1 = affected (case)).
#' @param test.data [data.frame] The test dataset. It must contain a column named Y providing the case/control phenotype (0 = unaffected (control), 1 = affected (case)).
#' @return \code{LR.model} return a vector containing the predicted probability.
#' @export
#' @examples
#' input.dir <- system.file("extdata", package="PrismVote")
#' output.dir <- system.file("extdata", package="PrismVote")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' covar.number <-  c(2, 3)
#'
#' feature.selection.result <- feature.selection(input.dir = input.dir,
#' output.dir = output.dir,
#' genotype = "train",
#' phenotype = "train.phenotypes.txt",
#' covar.number = covar.number,
#' plink.path = path2plink,
#' topK = 10,
#' verbose = TRUE)
#'
#' train.genotype.path <- file.path(input.dir, "train.raw")
#' train.genotype.path <- gsub('\\\\', '/', train.genotype.path)
#'
#' train.phenotype.path <- file.path(input.dir, "train.phenotypes.txt")
#' train.phenotype.path <- gsub('\\\\', '/', train.phenotype.path)
#'
#' train.data <- data.table::fread(train.genotype.path, data.table = FALSE)[, -c(1,2,3,4,5,6)]
#' colnames(train.data) <- unlist(purrr::map(colnames(train.data),function(x) {substr(x,1, nchar(x)-2)}))
#' train.pheno <- data.table::fread(train.phenotype.path, data.table = FALSE)
#' train.data$Y <- train.pheno[, 3]-1
#' train.data <- cbind(train.data[, feature.selection.result$index, drop = FALSE],
#' train.pheno[, covar.number + 2, drop = FALSE])
#'
#' test.genotype.path <- file.path(input.dir, "test.raw")
#' test.genotype.path <- gsub('\\\\', '/', test.genotype.path)
#'
#' test.phenotype.path <- file.path(input.dir, "test.phenotypes.txt")
#' test.phenotype.path <- gsub('\\\\', '/', test.phenotype.path)
#'
#' test.data <- data.table::fread(test.genotype.path, data.table = FALSE)[, -c(1,2,3,4,5,6)]
#' colnames(test.data) <- unlist(purrr::map(colnames(test.data),function(x) {substr(x,1, nchar(x)-2)}))
#' test.pheno <- data.table::fread(test.phenotype.path, data.table = FALSE)
#' test.data$Y <- test.pheno[, 3]-1
#' test.data <- cbind(test.data[, feature.selection.result$index, drop = FALSE],
#' test.pheno[, covar.number + 2, drop = FALSE])
#'
#' LR.pred <- LR.model(train.data, test.data)
#' }

LR.model <- function(train.data, test.data){
  lr.fit <- glm(Y~., data = train.data, family = binomial)
  pred<- predict(lr.fit, test.data, type = "response")
  return(pred)
}
