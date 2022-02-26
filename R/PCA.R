#' Principal component analysis
#'
#' A Function to run principal component analysis on input dataset using PLINK 1.9.
#'
#' @param input.dir [character] The full absolute path to the directory containing the training and test dataset. If \code{input.dir} is missing, the current working directory obtained by \code{getwd()} is used.
#' @param output.dir [character] The full absolute path where the result will be written to. If \code{output.dir} is missing, the current working directory obtained by \code{getwd()} is used.
#' @param train.genotype [character] The prefix of PLINK binary files (bed/bim/fam) of the training dataset.
#' @param test.genotype [character] The prefix of PLINK binary files (bed/bim/fam) of the test dataset.
#' @param PCA.separate [logical] If TURE, the principal components are calculated from the training dataset and then project the test dataset onto those principal components. If FALSE, the principal components are calculated from the combined data of the training and test dataset. The default value is FALSE.
#' @param PCs.count [numeric] To specify the number of top principal components that should be extracted. The default value is 10.
#' @param plink.path [character] The full absolute path to the PLINK executable file. The executable to run is path/to/plink.exe if you are on Windows operating system, for Unix-like operating system this is path/to/plink. If \code{plink.path} is NULL, the PLINK PATH should be added as a system environment variable.
#' @param verbose [logical] If TRUE, the PLINK log, error, and warning information are printed to standard out. The default value is TRUE.
#' @return \code{PCA} returns a list containing eigenvalues and eigenvectors of the training and test dataset:
#' \item{eigenvalue}{A vector containing the top eigenvalues according to \code{PCs.count} specified.}
#' \item{train.eigenvector}{A data frame containing the eigenvectors of the training dataset.}
#' \item{test.eigenvector}{A data frame containing the eigenvectors of the test dataset.}
#' @export
#' @examples
#' input.dir <- system.file("extdata", package="pv")
#' output.dir <- system.file("extdata", package="pv")
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' pca.result <- PCA(input.dir = input.dir,
#' output.dir = output.dir,
#' train.genotype = "train",
#' test.genotype = "test",
#' PCA.separate = FALSE,
#' PCs.count = 10,
#' plink.path = path2plink,
#' verbose = TRUE)
#' }
PCA <- function(input.dir, output.dir, train.genotype, test.genotype, PCA.separate = FALSE, PCs.count = 10, plink.path = NULL, verbose = TRUE){
  if(missing(input.dir)){
    input.dir <- getwd()
  }
  if(missing(output.dir)){
    output.dir <- getwd()
  }

  train.genotype.path <- file.path(input.dir, train.genotype)
  train.genotype.path <- gsub('\\\\', '/', train.genotype.path)
  test.genotype.path <- file.path(input.dir, test.genotype)
  test.genotype.path <- gsub('\\\\', '/', test.genotype.path)

  output.path <- file.path(output.dir, train.genotype)
  output.path <- gsub('\\\\', '/', output.path)

  # Check PLINK software access
  plink.path <- plinkQC::checkPlink(plink.path)

  # Check data available
  plinkQC:::checkFormat(train.genotype.path)
  plinkQC:::checkFormat(test.genotype.path)

  # Run principal component analysis
  if(PCA.separate){
    sys::exec_wait(plink.path,
                   args = c("--noweb", "--bfile", train.genotype.path, "--pca", "var-wts", PCs.count, "--out", output.path),
                   std_out = verbose,
                   std_err = verbose)
    eigenvalue <- read.table(paste0(output.path, '.eigenval'))$V1
    eigenvec.var <- read.table(paste0(output.path, '.eigenvec.var'), header=F, col.names=c("CHROM", "RS","A1", "A2", paste("PC", 1:PCs.count, sep="")))
    write.table(subset(eigenvec.var, select=c("RS","A1",paste("PC",1:PCs.count, sep=""))), file=paste0(output.path, "_pca_score_file.txt"), quote=F)

    train.eigenvector <- list()
    test.eigenvector <- list()
    for (i in 1:PCs.count) {
      sys::exec_wait(plink.path,
                     args = c("--noweb", "--bfile", train.genotype.path, "--score", paste0(output.path, "_pca_score_file.txt"), 2, 3, i + 3, "header", "--out", paste0(output.path,"_pca", i,'_score_train')),
                     std_out = verbose,
                     std_err = verbose)
      sys::exec_wait(plink.path,
                     args = c("--noweb", "--bfile", test.genotype.path, "--score", paste0(output.path, "_pca_score_file.txt"), 2, 3, i + 3, "header", "--out", paste0(output.path,"_pca", i,'_score_test')),
                     std_out = verbose,
                     std_err = verbose)

      train.eigenvector[[i]] <- read.table(paste0(output.path,"_pca", i,'_score_train.profile'), header=T)$SCORE
      test.eigenvector[[i]] <- read.table(paste0(output.path,"_pca", i,'_score_test.profile'), header=T)$SCORE
      train.eigenvector <- data.frame(train.eigenvector)
      test.eigenvector <- data.frame(test.eigenvector)
    }
  }else{
    sys::exec_wait(plink.path,
                   args = c("--noweb", "--bfile", train.genotype.path, "--bmerge", paste0(test.genotype.path, ".bed"), paste0(test.genotype.path, ".bim"), paste0(test.genotype.path, ".fam"), "--make-bed", "--out", paste0(train.genotype.path, ".merge")),
                   std_out = verbose,
                   std_err = verbose)
    sys::exec_wait(plink.path,
                   args = c("--noweb", "--bfile", paste0(train.genotype.path, ".merge"), "--pca", PCs.count, "--out", paste0(train.genotype.path, ".merge")),
                   std_out = verbose,
                   std_err = verbose)

    eigenvalue <- read.table(paste0(paste0(train.genotype.path, ".merge"), '.eigenval'))$V1
    eigenvector <- read.table(paste0(paste0(train.genotype.path, ".merge"), '.eigenvec'),row.names = 2)[,-1]

    train.fam <- read.table(paste0(train.genotype.path, ".fam"), header = F)
    test.fam <- read.table(paste0(test.genotype.path, ".fam"), header = F)
    merge.fam <- read.table(paste0(train.genotype.path, ".merge.fam"), header = F)
    train.eigenvector <- eigenvector[which(merge.fam$V2 %in% train.fam$V2),]
    test.eigenvector <- eigenvector[which(merge.fam$V2 %in% test.fam$V2),]
  }
  return(list(eigenvalue = eigenvalue,
              train.eigenvector = train.eigenvector,
              test.eigenvector = test.eigenvector))
}
