% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SBGNview.obj.fun.R
\name{print.SBGNview}
\alias{print.SBGNview}
\title{Generate image files}
\usage{
\method{print}{SBGNview}(x, ...)
}
\arguments{
\item{x}{A SBGNview class object.}

\item{...}{Other parameters passed to print.}
}
\value{
None
}
\description{
A wrapper to run function \code{\link{renderSbgn}} for all pathways in a SBGNview object and generate image files.
}
\examples{
### Use simulated data. Please see vignettes for more examples.
### Run `browseVignettes(package = "SBGNview")`
data('pathways.info','sbgn.xmls')
SBGNview.obj = SBGNview(simulate.data = TRUE,
                        sbgn.dir = './',
                        input.sbgn = 'P00001',
                        output.file = './test.local.file',
                        output.formats = c('pdf'),
                        min.gene.value = -1,
                        max.gene.value = 1)
 print(SBGNview.obj)
}
