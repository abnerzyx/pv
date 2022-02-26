#' Data stratification
#'
#' A Function to stratify samples to several strata using top principal components.
#'
#' @inheritParams PV
#' @param stratum.count [numeric] To specify the number of strata, as default the sample size of each stratum is N/stratum.count, N is the sample size.
#' @return \code{stratification} returns a list containing the following components:
#' \item{train.stratum}{A vector containing the stratum number each training sample belongs to.}
#' \item{train.stratum.index}{A list containing the index of training samples belonging to each stratum.}
#' \item{stratum.center}{A vector containing the center of each stratum.}
#' \item{train.distance.to.center}{A list containing the distance between the training samples and the center of each stratum.}
#' \item{test.distance.to.center}{A list containing the distance between the training samples and the center of each stratum.}
#' \item{train.prob.to.g}{A list containing the probability of the training samples having variable g under the hypothesis that the sample belongs to each stratum.}
#' \item{test.prob.to.g}{A list containing the probability of the test samples having variable g under the hypothesis that the sample belongs to each stratum.}
#' \item{train.prob.to.stratum}{A list containing the probability that the training samples belong to each stratum.}
#' \item{test.prob.to.stratum}{A list containing the probability that the test samples belong to each stratum.}
#' @seealso {\code{\link{PCA}}}
#' @export
#' @examples
#' input.dir <- system.file("extdata", package="pv")
#' output.dir <- system.file("extdata", package="pv")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' stratification.result <- stratification(input.dir = input.dir,
#' output.dir = input.dir,
#' train.genotype = "train",
#' test.genotype = "test",
#' stratum.count = 2,
#' PCA.separate = FALSE,
#' PCs.count = 10,
#' plink.path = path2plink,
#' CS = FALSE,
#' verbose = TRUE)
#' }
stratification <- function(input.dir, output.dir, train.genotype, test.genotype, stratum.count = 2 ,PCA.separate = FALSE, PCs.count = 10, plink.path = NULL, CS = FALSE, verbose = TRUE){
  if(missing(input.dir)){
    input.dir <- getwd()
  }
  if(missing(output.dir)){
    output.dir <- getwd()
  }

  PCA.result <- PCA(input.dir, output.dir, train.genotype, test.genotype , PCA.separate , PCs.count, plink.path , verbose)

  eigenvalue <- PCA.result$eigenvalue
  train.eigenvector <-  PCA.result$train.eigenvector
  test.eigenvector <-  PCA.result$test.eigenvector

  # calculate g value
  g.train <- cal.g(eigenvalue, train.eigenvector)
  g.test <- cal.g(eigenvalue, test.eigenvector)

  # Stratification of the training samples.
  train.stratum <- infotheo::discretize(g.train, "equalfreq", stratum.count)$X

  # The center of each stratum
  train.stratum.index <- list()
  stratum.center <- list()

  if(CS == TRUE){
    # cosine similarity
    stratum.eigenvector <- list()
    for (k in 1:stratum.count) {
      train.stratum.index[[k]] <- which(train.stratum == k)
      stratum.center[[k]] <- colMeans(train.eigenvector[train.stratum.index[[k]],])
      stratum.eigenvector[[k]] <- train.eigenvector[train.stratum.index[[k]],]
    }

    train.bayes.result <- prob.to.stratum.cs(stratum.center, stratum.eigenvector, train.eigenvector)
    test.bayes.result <- prob.to.stratum.cs(stratum.center, stratum.eigenvector, test.eigenvector)

  }else{
    stratum.g <- list()
    for (k in 1:stratum.count) {
      train.stratum.index[[k]] <- which(train.stratum == k)
      stratum.center[[k]] <- mean(g.train[train.stratum.index[[k]]])
      stratum.g[[k]] <- g.train[train.stratum.index[[k]]]
    }

    train.bayes.result <- prob.to.stratum(stratum.center, stratum.g, g.train)
    test.bayes.result <- prob.to.stratum(stratum.center, stratum.g, g.test)
  }

  return(list(train.stratum = train.stratum,
              train.stratum.index = train.stratum.index,
              stratum.center = stratum.center,
              train.distance.to.center = train.bayes.result$distance.to.center,
              test.distance.to.center = test.bayes.result$distance.to.center,
              train.prob.to.g = train.bayes.result$prob.to.g,
              test.prob.to.g = test.bayes.result$prob.to.g,
              train.prob.to.stratum = train.bayes.result$prob.to.stratum,
              test.prob.to.stratum = test.bayes.result$prob.to.stratum))
}

