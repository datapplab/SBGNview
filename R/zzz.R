.onLoad <- function(libname, pkgname){
  pkgNames = rownames(installed.packages())
  
  if("SBGNview" %in% pkgNames){
    data("pathways.info", package = "SBGNview")
    data("mapped.ids", package = "SBGNview")
  }
}
