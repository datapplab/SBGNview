.onLoad <- function(libname, pkgname){
  pkgNames = rownames(installed.packages())
  
  if("SBGNview" %in% pkgNames){
    data("pathways.info", "mapped.ids")
  }
}
