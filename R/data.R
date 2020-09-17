#' IDs mappable by SBGNview
#' @details ID types mappable to glyph IDs in SBGNhub SBGN-ML file collection.
#' @format A list with two elements: A vector of mappable gene ID types (\code{mapped.ids[["gene"]]}), a vector of mappable compound ID types (\code{mapped.ids[["cpd"]]}).
"mapped.ids"

#' Information of collected pathways
#' @details The information of pre-collected pathways and their SBGN-ML file information, such as pathway ID, name, source database, glyph ID types and URLs to the original pathway webpages.
#' @format A data.frame
"pathways.info"

#' Number of pathways collected
#' @details The number of pathways in SBGNhub from each source database. It is calculated from data 'pathways.info'.
#' @format A data.frame
"pathways.stat"

#' Demo SBGNview object
#' @details An example SBGNview object
#' @format A SBGNview object
"SBGNview.obj"

