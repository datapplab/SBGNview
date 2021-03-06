% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.utilities.R
\name{findPathways}
\alias{findPathways}
\title{Retrieve pathways by keywords}
\usage{
findPathways(
  keywords = NULL,
  keywords.logic = "or",
  keyword.type = "pathway.name",
  org = NULL,
  SBGNview.data.folder = "./SBGNview.tmp.data"
)
}
\arguments{
\item{keywords}{A character string or vector. The search is case in-sensitive.}

\item{keywords.logic}{A character string. Options are 'and' or 'or' (default). This will tell the function if the search require 'all' or 'any' of the keywords to be present. It only makes difference when keyword.type is 'pathway.name'.}

\item{keyword.type}{A character string. Either 'pathway.name' (default) or one of the ID types in \code{data('mapped.ids')}}

\item{org}{A character string. The KEGG species code.}

\item{SBGNview.data.folder}{A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files.}
}
\value{
A dataframe. Contains information of pathways found.
}
\description{
This function searches for pathways by input keywords.
}
\details{
If 'keyword.type' is 'pathway.name' (default), this function will search for the presence of any keyword in the pathway.name column of data(pathways.info). The search is case in-sensitive. If 'keyword.type' is one of the identifier types and 'keywords' are corresponding identifiers, this function will return pathways that include nodes mapped to input identifiers.
}
\examples{
data(pathways.info)
input.pathways <- findPathways('Adrenaline and noradrenaline biosynthesis')
}
