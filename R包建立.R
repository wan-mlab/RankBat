#

# package in need ---------------------------------------------------------

library('devtools')
library('knitr')
devtools::has_devel(debug=TRUE)
pkgbuild::check_build_tools(debug = TRUE)


# 1.create package --------------------------------------------------------

setwd('E:/')
library('devtools')
create_package('RankBat')


#之后将去建立的R.Proj环境中去进行后续的步骤

#> getwd()
#[1] "E:/R-Studio/RStudio/new"
# > setwd('E:/')
# > getwd()
# [1] "E:/"
# > setwd('E:/RankBat/')
# > devtools::install()
# ── R CMD build ──────────────────────────────────────────────────────────────────────────────────────
# ✔  checking for file 'E:\RankBat/DESCRIPTION' ...
# ─  preparing 'RankBat':
#   ✔  checking DESCRIPTION meta-information ...
# ─  checking for LF line-endings in source and make files and shell scripts
# ─  checking for empty or unneeded directories
# ─  building 'RankBat_0.0.0.9000.tar.gz'
#
# Running "E:/R/R-4.3.1/bin/x64/Rcmd.exe" INSTALL \
# "C:\Users\86178\AppData\Local\Temp\Rtmpa4UTbS/RankBat_0.0.0.9000.tar.gz" --install-tests
# * installing to library 'E:/R/R-4.3.1/library'
# * installing *source* package 'RankBat' ...
# ** using staged installation
# ** R
# ** byte-compile and prepare package for lazy loading
# ** help
# *** installing help indices
# ** building package indices
# ** testing if installed package can be loaded from temporary location
# ** testing if installed package can be loaded from final location
# ** testing if installed package keeps a record of temporary installation path
# * DONE (RankBat)
