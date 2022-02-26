prob.to.stratum <- function(stratum.center, stratum.g, g){
  distance.to.center <- list()
  prob.to.g <- list()
  for(i in 1:length(stratum.center)){
    distance.to.center[[i]] <- (g - stratum.center[[i]])^2/var(stratum.g[[i]])
    prob.to.g[[i]] <- pchisq(distance.to.center[[i]], df = 1, lower.tail = F)
  }
  prob.to.g.df <- data.frame(prob.to.g)
  prob.to.stratum <- as.list(data.frame(t(apply(prob.to.g.df, 1, function(x) x/sum(x)))))
  names(prob.to.stratum) <- NULL
  return(list(distance.to.center = distance.to.center,
              prob.to.g = prob.to.g,
              prob.to.stratum = prob.to.stratum))
}

# cosine similarity
prob.to.stratum.cs <- function(stratum.center, stratum.eigenvector, eigenvector){
  distance.to.center <- list()
  prob.to.g <- list()
  for(i in 1:length(stratum.center)){
    distance.to.center[[i]] <- apply(eigenvector, 1, function(x) lsa::cosine(x,stratum.center[[i]]))
    p.temp <- apply(stratum.eigenvector[[i]], 1, function(x) lsa::cosine(x,stratum.center[[i]]))
    prob.to.g[[i]] <- exp(distance.to.center[[i]]) / (sum(exp(p.temp)) + exp(distance.to.center[[i]]))

  }
  prob.to.g.df <- data.frame(prob.to.g)
  prob.to.stratum <- as.list(data.frame(t(apply(prob.to.g.df, 1, function(x) x/sum(x)))))
  names(prob.to.stratum) <- NULL
  return(list(distance.to.center = distance.to.center,
              prob.to.g = prob.to.g,
              prob.to.stratum = prob.to.stratum))
}

cal.g <- function(eigenvalue, eigenvector){
  a <- eigenvalue / sum(eigenvalue)
  g <- apply(eigenvector, 1, function(x) sum(x * a))
  return(g)
}

checkMissingness <- function(input.dir, output.dir, data.name, plink.path = NULL, verbose = TRUE){
  if(missing(input.dir)){
    input.dir <- getwd()
  }
  if(missing(output.dir)){
    output.dir <- getwd()
  }

  missing.flag <- 1
  data.path <- file.path(input.dir, data.name)
  data.path <- gsub('\\\\', '/', data.path)

  output.path <- file.path(output.dir, data.name)
  output.path <- gsub('\\\\', '/', output.path)

  # Check PLINK software access
  plink.path <- plinkQC::checkPlink(plink.path)

  # Check data available
  plinkQC:::checkFormat(data.path)

  sys::exec_wait(plink.path,
                 args = c("--noweb", "--bfile", data.path, "--missing", "--out", output.path),
                 std_out = verbose,
                 std_err = verbose)
  sample.missing <- data.table::fread(paste0(output.path, ".imiss"), data.table = FALSE)
  SNP.missing <- data.table::fread(paste0(output.path, ".lmiss"), data.table = FALSE)
  if(!all(sample.missing$N_MISS == 0)){
    warning(data.path, ".bed with missing value! \nIID: ", paste(sample.missing$IID[which(sample.missing$N_MISS != 0)],sep="",collapse = ", "), "\nSNP: ", paste(SNP.missing$SNP[which(SNP.missing$N_MISS != 0)],sep="",collapse = ", "),"\nThe missing SNP values are imputed according to its probability distribution.")
  }else{
    missing.flag <- 0
  }
  return(missing.flag)
}

impute <- function(genotype, seed = 1234){
  geno.0 <- sum(genotype == 0, na.rm = TRUE)
  geno.1 <- sum(genotype == 1, na.rm = TRUE)
  geno.2 <- sum(genotype == 2, na.rm = TRUE)
  na.index <- which(is.na(genotype))
  if(length(na.index)>0){
    set.seed(seed)
    genotype[na.index] <- sample(c(0:2), length(na.index), prob = c(geno.0, geno.1, geno.2)/sum(c(geno.0, geno.1, geno.2)), replace = TRUE)
  }
  return(genotype)
}

data.split <- function(input.dir, output.dir, genotype, phenotype, stratum.number, plink.path = NULL, verbose = TRUE){
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

  phenotype.path <- file.path(input.dir, phenotype)###### updated "output.dir" to "input.dir"
  phenotype.path <- gsub('\\\\', '/', phenotype.path)
  # Check PLINK software access
  plink.path <- plinkQC::checkPlink(plink.path)

  fam <- data.table::fread(paste0(genotype.path, ".fam"), data.table = FALSE)
  phenotype <- data.table::fread(phenotype.path, data.table = FALSE)

  for(i in 1:length(unique(stratum.number))){
    FID_IID <- fam[which(stratum.number == i), c(1,2)]
    colnames(FID_IID) <- c("FID", "IID")
    write.table(FID_IID,file = paste0(output.path, ".stratum.",i,".sampleID.txt"),row.names = FALSE, col.names = TRUE, quote = FALSE)
    FID_IID$Y <- phenotype[which(stratum.number == i), 3]
    FID_IID <- cbind(FID_IID, phenotype[which(stratum.number == i), -c(1:3)])
    write.table(FID_IID,file = paste0(output.path, ".stratum.",i,".phenotype.txt"),row.names = FALSE, col.names = TRUE, quote = FALSE)
    sys::exec_wait(plink.path,
                   args = c("--noweb", "--bfile", genotype.path, "--keep", paste0(output.path, ".stratum.",i,".sampleID.txt"), "--recodeA", "--make-bed", "--out", paste0(output.path, ".stratum.",i)),
                   std_out = verbose,
                   std_err = verbose)
  }
}

model.evaluation <- function(response, predicted.probability){
  AUC <- pROC::auc(response, predicted.probability)
  Acc <- sum(response == ifelse(predicted.probability >= 0.5, 1, 0))/length(response)
  return(list(AUC = c(AUC),
              Acc = Acc))
}
