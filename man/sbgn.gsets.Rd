% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.utilities.R
\name{sbgn.gsets}
\alias{sbgn.gsets}
\title{Retrieve gene list or compound list from collected databases}
\usage{
sbgn.gsets(
  id.type = "ENTREZID",
  mol.type = "gene",
  species = "hsa",
  database = "pathwayCommons",
  output.pathway.name = TRUE,
  combine.duplicated.set = TRUE,
  truncate.name.length = 50,
  SBGNview.data.folder = "./SBGNview.tmp.data"
)
}
\arguments{
\item{id.type}{A character string. Default: "ENTREZID". The ID type of output gene list. One of the supported types in \code{data('mapped.ids')}}

\item{mol.type}{A character string. Default: "gene". One of 'gene' or 'cpd'}

\item{species}{A character string. Default: "hsa". The three letter species code used by KEGG. E.g. 'hsa','mmu'}

\item{database}{A character string. Default: "pathwayCommons". The database where gene list will be extracted. Acceptable values: 'MetaCyc', 'pathwayCommons', 'MetaCrop'. The value is case in-sensitive.}

\item{output.pathway.name}{Logical. Default: T. If set to 'TRUE', the names of returned list are in the format: 'pathway.id::pathway.name'. If set to 'FALSE', the format is 'pahtway.id'}

\item{combine.duplicated.set}{Logical. Default: T. Some pathways have the same geneset. If this parameter is set to 'TRUE', the output list will combine pathways that have the same gene set. The name in the list will be pathway names concatinated with '||'}

\item{truncate.name.length}{Integer. Default: 50. The pathway names will be truncated to at most that length.}

\item{SBGNview.data.folder}{A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files.}
}
\value{
A list. Each element is a genelist of a pathway.
}
\description{
Retrieve gene list or compound list from collected databases
}
\examples{
data(pathways.info)
mol.list <- sbgn.gsets(id.type = 'ENTREZID',
                       species = 'hsa',
                       database = 'pathwayCommons')
  
}
