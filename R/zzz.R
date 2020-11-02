
#########################################################################################################
# Load R data objects when SBGNview package is loaded and print startup message
.onLoad <- function(libname, pkgname){
  pkgNames = rownames(installed.packages())
  
  if("SBGNview" %in% pkgNames){
    data("pathways.info", package = "SBGNview")
    data("mapped.ids", package = "SBGNview")
  }
  
  message = "\n##############################################################################\nSBGNview is an open source R software package distributed under GNU General Public License version 3 (GPLv3). Details of GPLv3 are available at http://www.gnu.org/licenses/gpl-3.0.html. For more information about SBGNview, visit https://bioconductor.org/packages/release/bioc/html/SBGNview.html.\n\nTo get started, run browseVignettes(\"SBGNview\") to view vignettes for SBGNview. Use ls(\"package:SBGNview\") to get list of avaliable functions and use help(\"<functionName>\") to access a function's documentation.\n##############################################################################\n"
  
  packageStartupMessage(message)
  
}

#########################################################################################################
