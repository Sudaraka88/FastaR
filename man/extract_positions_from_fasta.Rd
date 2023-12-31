% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractSNPs.R
\name{extract_positions_from_fasta}
\alias{extract_positions_from_fasta}
\title{extract_positions_from_fasta}
\usage{
extract_positions_from_fasta(
  aln_path,
  pos,
  save_path = "out.fa",
  return_to_environment = T
)
}
\arguments{
\item{aln_path}{path to multi fasta alignment}

\item{pos}{vector of positions to extract. These positions can be SNP or non-SNP, no filtering will be performed. Set pos = 0 to return full MSA in matrix format to R environment.
WARNING: might hang in systems with insufficient RAM. Use at your own risk.}

\item{save_path}{path to save fasta (default = out.fa). POS and SEQS files will also be created. Set to NULL to avoid saving. WARNING: Files will be overwritten without checking!}

\item{return_to_environment}{return snps R matrix back to environment (default = T)}
}
\value{
If return_to_environment=T, R character matrix with requested positions. colnames: positions, rownames: sequence names
}
\description{
Function to extract nucleotides at chosen positions from a fasta alignment.
}
\examples{
\dontrun{
aln_path <- system.file("extdata", "test.fa", package = "FastaR")
snp.dat <- extract_positions_from_fasta(aln_path, c(1,2,3))
}
}
