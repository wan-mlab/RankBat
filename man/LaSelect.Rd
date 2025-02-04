% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step3_LASSO_featureselection.R
\name{LaSelect}
\alias{LaSelect}
\title{extract important gene ratio;}
\usage{
LaSelect(data, sampAnnote, all_rank_t_genes, all_reversed_gp_genes)
}
\arguments{
\item{data}{data which you want to extract gene ratios from;}

\item{sampAnnote}{sampAnnot which is the annotation file for data;}

\item{all_rank_t_genes}{all_rank_t_genes which is the differentially ranked genes selected in previous RCA() function.}

\item{all_reversed_gp_genes}{all_reversed_gp_genes which is the differentially ranked genes ratio selected in previous RRA() function.}
}
\value{
a global object MB_RANK_GP
}
\description{
this function is used to extract important gene ratios from dataset for each subtype based on genes ratios selected from RRA() function.
}
\examples{
MB_RANK_GP<-LaSelect(GSE85217, sampAnnote_GSE85217, all_rank_t_genes,all_reversed_gp_genes)
}
