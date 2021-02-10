
#########################################################################################################
#' Generic function to modify SBGN graph features
#' 
#' This binary operator ('+') modifies a SBGNview object(first argument) using a function (second argument)
#' 
#' @param SBGNview.obj An object of class SBGNview
#' @param fn A character string. The name of any funtion that modifies and returns a SBGNview object. Some functions are available in SBGNview package: \code{\link{highlightPath}}, \code{\link{highlightNodes}}.
#' @return The returned value of *fn*
#' @examples 
#' \dontrun{
#' obj.new <- SBGNview.obj + highlightArcs(class = 'production',
#'                                         color = 'red') 
#' }
#' @export

"+.SBGNview" <- function(SBGNview.obj, fn) {
    fn(SBGNview.obj)
}

#########################################################################################################
#' Generate image files
#' 
#' A wrapper to run function \code{\link{renderSbgn}} for all pathways in a SBGNview object and generate image files.
#' 
#' @param x A SBGNview class object. 
#' @param ... Other parameters passed to print.
#' @return None
#' @examples 
#' ### Use simulated data. Please see vignettes for more examples.
#' ### Run `browseVignettes(package = "SBGNview")`
#' data('pathways.info','sbgn.xmls')
#' SBGNview.obj = SBGNview(simulate.data = TRUE,
#'                         sbgn.dir = './',
#'                         input.sbgn = 'P00001',
#'                         output.file = './test.local.file',
#'                         output.formats = c('pdf'),
#'                         min.gene.value = -1,
#'                         max.gene.value = 1)
#'  print(SBGNview.obj)
#' @export

"print.SBGNview" <- function(x, ...) {
    
    output.file <- NULL
    #SBGNview.obj <- x
    SBGNview.obj <- merge.entity.specific.parameters.list(x)
    glyphs.arcs.list <- SBGNview.obj$data
    sbgns <- names(glyphs.arcs.list)
    for (s in seq_len(length.out = length(glyphs.arcs.list))) {
        # for each sbgn file
        data.this.sbgn <- glyphs.arcs.list[[s]]
        sbgn.parameters.list <- data.this.sbgn$render.sbgn.parameters.list
        if (!is.null(output.file)) {
            sbgn.parameters.list$output.file <- gsub(SBGNview.obj$output.file, output.file, 
                sbgn.parameters.list$output.file)
        }
        
        # if input sbgn-ml file doesn't exist while printing: download file if in ID in pathways.info, else user's file
        if(!file.exists(sbgn.parameters.list$input.sbgn)) {
            if(sbgns[s] %in% pathways.info[, "pathway.id"]) {
                message(sbgn.parameters.list$input.sbgn, " not in current working directory")
                message("Downloading SBGN-ML file for pathway.id: ", sbgns[s])
                downloadSbgnFile(pathway.id = sbgns[s])
            } else {
                stop(sbgn.parameters.list$input.sbgn, " not in current working directory.\nPlease make sure SBGN-ML file is in current working directory")
            }
        }
            
        tp <- renderSbgn(input.sbgn = sbgn.parameters.list$input.sbgn, output.file = sbgn.parameters.list$output.file, 
            arcs.info = sbgn.parameters.list$arcs.info, compartment.layer.info = sbgn.parameters.list$compartment.layer.info, 
            user.data = sbgn.parameters.list$user.data, output.formats = SBGNview.obj$output.formats, 
            sbgn.id.attr = sbgn.parameters.list$sbgn.id.attr, glyphs.user = data.this.sbgn$glyphs.list, 
            arcs.user = data.this.sbgn$arcs.list, pathway.name = sbgn.parameters.list$pathway.name, 
            global.parameters.list = data.this.sbgn$global.parameters.list)
        message("Image files written: ", sbgn.parameters.list$output.file, "\n")
    }
    return(invisible())
}

#########################################################################################################
# function takes an SBGNview object and check the glyph objects parameters.list slot.
# The entity specific parameters.list can be empty or partial if user has customized the glyph 
# The function below will return a SBGNview object where the glyph and arc objects' parameters.list slot
# will be merged with global.parameters.list to create a full list to be used by downstream functions
merge.entity.specific.parameters.list <- function(obj) {
    
    for(i in seq_along(obj$data)) { # for multiple input.sbgn IDs
        
        global.parameters.list <- obj$data[[i]]$global.parameters.list
        glyphs.list <- obj$data[[i]]$glyphs.list
        #arcs.list <- obj$data[[i]]$arcs.list 
        
        for(j in seq_along(glyphs.list)) { # check each glyph's parameter.list slot
            if(length(glyphs.list[[j]]@parameters.list) == 0) {
                glyphs.list[[j]]@parameters.list <- global.parameters.list
            } else {
                partial.params.list <- glyphs.list[[j]]@parameters.list
                diff.list.names <- c(setdiff(names(global.parameters.list), names(partial.params.list)))
                # merged parameters.list
                glyphs.list[[j]]@parameters.list <- c(partial.params.list, global.parameters.list[diff.list.names]) 
            }
        } # end for loop for glyphs.list
        
        obj$data[[i]]$glyphs.list <- glyphs.list
        
    } # end main for loop
    
    return(obj)
}

#########################################################################################################
#' Retrieve output file information from a SBGNview object
#' 
#' @param obj An SBGNview class object.
#' @details This function prints the output file path recorded in a SBGNview object. 
#' @return A string. The output file information. 
#' @examples 
#' \dontrun{
#' fileOutput(SBGNview.obj) 
#' }
#' @export

fileOutput <- function(obj) {
    obj$output.file
}

#########################################################################################################
#' Set or change output file information for a SBGNview object
#' 
#' @param obj A SBGNview class object
#' @param value No need to provide 
#' @details This function sets the output file path recorded in a SBGNview object.
#' @return A SBGNview class object
#' @examples 
#' \dontrun{
#' outputFile(SBGNview.obj) <- "test.output.file"
#' }
#' @export
"outputFile<-" <- function(obj, value) {
    
    glyphs.arcs.list <- obj$data
    sbgns <- names(glyphs.arcs.list)
    for (s in seq_len(length.out = length(glyphs.arcs.list))) {
        # for each sbgn file
        data.this.sbgn <- glyphs.arcs.list[[s]]
        sbgn.parameters.list <- data.this.sbgn$render.sbgn.parameters.list
        new.output.file <- gsub(obj$output.file, value, sbgn.parameters.list$output.file)
        obj$data[[s]]$render.sbgn.parameters.list$output.file <- new.output.file
    }
    obj$output.file <- value
    return(obj)
}

#########################################################################################################
