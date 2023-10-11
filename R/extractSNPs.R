#' extract_snps_from_fasta
#'
#' Function to extract SNP positions from fasta alignment. If you require non-SNP positions, use extract_positions_from_fasta() instead.
#'
#' @importFrom Rcpp sourceCpp
#'
#' @param aln_path path to multi fasta alignment
#' @param gap_freq sites with 'N' allele (ambiguous/gap) frequency >gap_freq will be dropped (default = 0.15)
#' @param maf_freq sites with second most common allele (minor allele) frequency <maf_freq will be dropped (default = 0.01)
#' @param is_N_minor_allele specify whether 'N' allele (ambiguous/gap) should be considered as minor allele (default = F). Warning! Can lead to poor quality SNPs
#' @param save_path path to save fasta (default = out.fa). POS and SEQS files will also be created. Set to NULL to avoid saving. WARNING: Files will be overwritten without checking!
#' @param return_to_environment return snps R matrix back to environment (default = T)
#'
#' @return If return_to_environment=T, R character matrix with requested positions. colnames: positions, rownames: sequence names
#' @useDynLib FastaR
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "test.fa", package = "FastaR")
#' snp.dat <- extract_snps_from_fasta(aln_path)
#' }
#' @export
extract_snps_from_fasta <- function(aln_path, gap_freq = 0.15, maf_freq = 0.01, is_N_minor_allele = F, save_path = "out.fa", return_to_environment = T){
  aln_path = normalizePath(aln_path) # C safety (~ character causes crash)
  if(!file.exists(aln_path)) stop(paste("Can't locate file", aln_path))
  if(is_N_minor_allele == T) {
    warning("Warning! Ambiguous characters are considered minor allele, might result in lower quality SNPs")
    filter = 1
  } else {
    filter = 0
  }

  snp.param <- .extractAlnParam(aln_path, filter, gap_freq, maf_freq)

  if(snp.param$seq.length==-1) stop("Error! sequences are of different lengths!")
  if(snp.param$num.seqs==0) stop("File does not contain any sequences!")
  if(snp.param$num.snps==0) stop("File does not contain any SNPs")

  snps <- .getSNPs(aln_path, snp.param$num.seqs, snp.param$num.snps, snp.param$pos)
  snps <- matrix(snps, nrow = snp.param$num.seqs, length(snp.param$pos), byrow = T)
  colnames(snps) = snp.param$pos
  rownames(snps) = snp.param$seq.names

  if(!is.null(save_path)) FastaR::write_fasta(snps, save_path)

  if(return_to_environment) return(snps)
}


#' extract_positions_from_fasta
#'
#' Function to extract positions from fasta alignment.
#'
#' @importFrom Rcpp sourceCpp
#'
#' @param aln_path path to multi fasta alignment
#' @param pos vector of positions to extract. These positions can be SNP or non-SNP, no filtering will be performed. Set pos = 0 to return full MSA in matrix format to R environment.
#' WARNING: might hang in systems with insufficient RAM. Use at your own risk.
#' @param save_path path to save fasta (default = out.fa). POS and SEQS files will also be created. Set to NULL to avoid saving. WARNING: Files will be overwritten without checking!
#' @param return_to_environment return snps R matrix back to environment (default = T)
#'
#' @return If return_to_environment=T, R character matrix with requested positions. colnames: positions, rownames: sequence names
#' @useDynLib FastaR
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "test.fa", package = "FastaR")
#' snp.dat <- extract_positions_from_fasta(aln_path, c(1,2,3))
#' }
#' @export
extract_positions_from_fasta <- function(aln_path, pos, save_path = "out.fa", return_to_environment = T){
  # Check inputs
  aln_path = normalizePath(aln_path) # C safety (~ character causes crash)
  if(!file.exists(aln_path)) stop(paste("Can't locate file", aln_path))

  # sanity checks
  if(!is.numeric(pos)) stop("Error! pos must be numeric!")
  if(any(duplicated(pos))) stop("Error! pos cannot contain duplicates!")
  pos = pos[order(pos)] # order the positions

  snp.param <- .extractAlnParam2(aln_path) # This is only to extract seq.length
  if(snp.param$seq.length==-1) stop("Error! sequences are of different lengths!")
  if(snp.param$num.seqs==0) stop("File does not contain any sequences!")
  if(pos[length(pos)] > snp.param$seq.length) stop("Error! Requested position(s) are outside sequence length")

  if(length(pos) == 1){
    if(pos == 0){ # requets for full MSA, provide warning, this might crash machines with low RAM
      cat("Returning full MSA to environment! WARNING: Might crash the system if the MSA is very large!\n")
      snps <- .getSNPs(aln_path, snp.param$num.seqs, snp.param$seq.length, 1:snp.param$seq.length)
      snps <- matrix(snps, nrow = snp.param$num.seqs, length(pos), byrow = T)
      colnames(snps) = pos
      rownames(snps) = snp.param$seq.names
      return(snps)}
  }

  snps <- .getSNPs(aln_path, snp.param$num.seqs, length(pos), pos)
  snps <- matrix(snps, nrow = snp.param$num.seqs, length(pos), byrow = T)
  colnames(snps) = pos
  rownames(snps) = snp.param$seq.names

  if(!is.null(save_path)) FastaR::write_fasta(snps, save_path)

  if(return_to_environment) return(snps)
}

#' extract_snps_from_fasta_sparse
#'
#' Function to extract SNPs from a multi-fasta alignment and return to R in sparse matrix format. This format to perform fast computations such as cross products,
#' see https://github.com/Sudaraka88/LDWeaver/blob/main/R/BacGWES.R for a use case.
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix sparseMatrix t
#'
#' @param aln_path path to multi fasta alignment
#' @param gap_freq sites with 'N' allele (ambiguous/gap) frequency >gap_freq will be dropped (default = 0.15)
#' @param maf_freq sites with second most common allele (minor allele) frequency <maf_freq will be dropped (default = 0.01)
#' @param is_N_minor_allele specify whether 'N' allele (ambiguous/gap) should be considered as minor allele (default = F). Warning! Can lead to poor quality SNPs
#'
#' @return R list with SNPs in sparse format and additional parameters for fast computations
#' @useDynLib FastaR
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "test.fa", package = "FastaR")
#' snp.dat <- extract_snps_from_fasta_sparse(aln_path)
#' }
#' @export
extract_snps_from_fasta_sparse <- function(aln_path, gap_freq = 0.15, maf_freq = 0.01, is_N_minor_allele = F){
  # Check inputs
  aln_path = normalizePath(aln_path) # C safety (~ character causes crash)
  if(!file.exists(aln_path)) stop(paste("Can't locate file", aln_path))


  if(is_N_minor_allele == T) {
    warning("Warning! Ambiguous characters are considered minor allele, might result in lower quality SNPs")
    filter = 1
  } else {
    filter = 0
  }

  snp.param <- .extractAlnParam(aln_path, filter, gap_freq, maf_freq)

  if(snp.param$seq.length==-1) stop("Error! sequences are of different lengths!")
  if(snp.param$num.seqs==0) stop("File does not contain any sequences!")
  if(snp.param$num.snps==0) stop("File does not contain any SNPs")

  snp.data <- .extractSNPs(aln_path, snp.param$num.seqs, snp.param$num.snps, snp.param$pos)
  seq.names <- gsub("^>","",snp.data$seq.names); snp.data$seq.names = NULL
  uqe = apply(snp.data$ACGTN_table>0, 1, function(x) as.numeric(x>0)); snp.data$ACGTN_table = NULL

  snp.matrix_A <- Matrix::sparseMatrix(i=snp.data$i_A,
                                       j=snp.data$j_A,
                                       x=as.logical(snp.data$x_A),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL

  snp.matrix_C <- Matrix::sparseMatrix(i=snp.data$i_C,
                                       j=snp.data$j_C,
                                       x=as.logical(snp.data$x_C),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL

  snp.matrix_G <- Matrix::sparseMatrix(i=snp.data$i_G,
                                       j=snp.data$j_G,
                                       x=as.logical(snp.data$x_G),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL

  snp.matrix_T <- Matrix::sparseMatrix(i=snp.data$i_T,
                                       j=snp.data$j_T,
                                       x=as.logical(snp.data$x_T),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL

  snp.matrix_N <- Matrix::sparseMatrix(i=snp.data$i_N,
                                       j=snp.data$j_N,
                                       x=as.logical(snp.data$x_N),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data = NULL

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A), snp.matrix_C=Matrix::t(snp.matrix_C),
              snp.matrix_G=Matrix::t(snp.matrix_G), snp.matrix_T=Matrix::t(snp.matrix_T),
              snp.matrix_N=Matrix::t(snp.matrix_N), g = snp.param$seq.length, nsnp = snp.param$num.snps,
              nseq = snp.param$num.seqs, seq.names = seq.names, r = rowSums(uqe), uqe = uqe, POS = snp.param$pos))
}



#' write_fasta
#'
#' Function to write a sequence matrix to fasta
#'
#' @importFrom utils write.table
#'
#' @param snps snps in matrix format
#' @param fa_path path to save fasta
#'
#'
#' @examples
#' \dontrun{
#' aln_path <- write_fasta(snps, "snps.fa")
#' }
#' @export
write_fasta <- function(snps, fa_path){
  seq_path = paste(fa_path, ".seqs", sep = "")
  pos_path = paste(fa_path, ".pos", sep = "")

  ## clean up if existing and send a warning
  if(file.exists(fa_path)) {unlink(fa_path); cat("WARNING! Overwriting", fa_path, "\n")}
  if(file.exists(seq_path)) {unlink(seq_path); cat("WARNING! Overwriting", seq_path, "\n")}
  if(file.exists(pos_path)) {unlink(pos_path); cat("WARNING! Overwriting", pos_path, "\n")}

  for(i in 1:nrow(snps)){
    write.table(x = paste(">",rownames(snps)[i], sep = ""), file = fa_path, append = T, quote = F, row.names = F, col.names = F)
    write.table(x = paste(snps[i,], collapse = ""), file = fa_path, append = T, quote = F, row.names = F, col.names = F)
  }
  write.table(rownames(snps), seq_path, col.names = F, row.names = F, quote = F)
  write.table(colnames(snps), pos_path, col.names = F, row.names = F, quote = F)
}
