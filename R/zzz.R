
#########################################################################################################
# Load R data objects when SBGNview package is loaded and print start up message
.onLoad <- function(libname, pkgname) {
  
  data("pathways.info", package = "SBGNview")
  data("mapped.ids", package = "SBGNview")
  data("SBGNhub.id.mapping.tables", package = "SBGNview")
  
  disclaimer <- "\n##############################################################################\nSBGNview is an open source software package distributed under GNU General Public License version 3 (GPLv3). Details of GPLv3 is available at http://www.gnu.org/licenses/gpl-3.0.html. \n\nUsers are required to formally cite the SBGNview paper and Pathview paper (not just mention them) in publications or products. For details, do 'citation(\"SBGNview\");  citation(\"Pathview\")' within R.\n##############################################################################"
  disclaimer <- pathview::wordwrap(disclaimer, 80)
  packageStartupMessage(disclaimer)
  
}

#########################################################################################################
