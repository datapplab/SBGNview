
#########################################################################################################
#' IDs mappable by SBGNview
#' @details ID types mappable to glyph IDs in SBGNhub SBGN-ML file collection.
#' @format A list with two elements: A vector of mappable gene ID types (\code{mapped.ids[["gene"]]}), a vector of mappable compound ID types (\code{mapped.ids[["cpd"]]}).
#' @usage data(mapped.ids)
"mapped.ids"

#########################################################################################################
#' Information of collected pathways
#' @details The information of pre-collected pathways and their SBGN-ML file information, such as pathway ID, name, source database, glyph ID types and URLs to the original pathway webpages.
#' @format A data.frame
#' @usage data(pathways.info)
"pathways.info"

#########################################################################################################
#' Number of pathways collected
#' @details The number of pathways in SBGNhub from each source database. It is calculated from data 'pathways.info'.
#' @format A data.frame
#' @usage data(pathways.stats)
"pathways.stats"

#########################################################################################################
#' Mapping tables available in SBGNhub
#' @details A collection of ID mapping table files available in SBGNhub. 
#' If a matching mapping table is available, it will be downloaded from \href{https://github.com/datapplab/SBGNhub/tree/master/data/id.mapping.unique.pair.name}{SBGNhub} github page.
#' @format A matrix
#' @usage data(SBGNhub.id.mapping.tables)
"SBGNhub.id.mapping.tables"

#########################################################################################################
