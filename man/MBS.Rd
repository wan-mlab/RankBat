% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step4_MB_subtyping.R
\name{MBS}
\alias{MBS}
\title{medulloblastoma (MB) subtyping}
\usage{
MBS(data, MB_RANK_GP, sampannote)
}
\arguments{
\item{data}{data contains the samples information, row is gene, column is sample}

\item{MB_RANK_GP}{MB_RANK_GP is the subtype-specific gene ratios selected from LaSelect() function}

\item{sampannote}{sampannote is the annotation of data}
}
\value{
a global object myMat
}
\description{
this function is used to subtype given samples. Based on gene ratios we selected from previous steps, for each sample, we predict its subtype.
}
\examples{
myMat<- (GSE21140, MB_RANK_GP, sampAnnote_GSE21140)
}
