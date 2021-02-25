
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

# get glyphs and arcs from object and use using plot.glyph, plot.arc functions 
# and other svg parameters from object and generate output svg and other files
# reconstruct entity (glyph/arc) specific parameters list
"print.SBGNview" <- function(x, ...) {
    object <- x
    for(i in seq_along(object$data)) {
        
        input.sbgn <- object$data[[i]]
        # reconstruct entity parameters list
        global.parameters.list <- input.sbgn$global.parameters.list
        glyphs <- lapply(input.sbgn$glyphs.list, FUN = get.entity.parameter.list, 
                         global.parameters.list = global.parameters.list)
        arcs <- lapply(input.sbgn$arcs.list, FUN = get.entity.parameter.list, 
                       global.parameters.list = global.parameters.list)
        svg.nodes.compartment <- ""
        svg.nodes.complex <- ""
        svg.nodes <- ""
        svg.arcs <- ""
        svg.cardinality <- ""
        svg.ports <- ""
        
        # plot glyphs
        for(glyph in glyphs) {
           
            # there is no plot.glyph() method signature for port class so skip. 
            # we plot the port if glyph@port slot contains svg code for plotting port
            if(is(glyph, "port")) next 
            
            # plot different classes of glyphs and store in respective svg variables
            if (is(glyph, "compartment.sbgn")) {
                svg.nodes.compartment <- paste(svg.nodes.compartment, plot.glyph(glyph), sep = "\n")
            } else if (is(glyph, "complex.sbgn")) {
                svg.nodes.complex <- paste(svg.nodes.complex, plot.glyph(glyph), sep = "\n")
            } else if (is(glyph, "cardinality.sbgn") | is(glyph, "stoichiometry.sbgn")) {
                if(input.sbgn$global.parameters.list$if.plot.cardinality) {
                    svg.cardinality <- paste(svg.cardinality, plot.glyph(glyph), sep = "\n")
                }
            } else {
                # plot all other glyphs
                svg.nodes <- paste(svg.nodes, plot.glyph(glyph), sep = "\n")
            }
            
            # if glyph contains port, get port svg
            if(!identical(glyph@svg.port, character(0))) {
                svg.ports <- paste(svg.ports, glyph@svg.port, sep = "\n")
            } 
            # if glyph contains clone, plot clone
            if(!length(glyph@clone) == 0) {
                # reconstruct clone glyph parameters.list slot
                clone.glyphs <- lapply(glyph@clone , FUN = get.entity.parameter.list, 
                                       global.parameters.list = global.parameters.list)
                for(i in seq_along(clone.glyphs)) {
                    svg.nodes <- paste(svg.nodes, plot.glyph(clone.glyphs[[i]]), sep = "\n") 
                }
            }
        }
        # plot arcs
        for(arc in arcs) {
            svg.arcs <- paste(svg.arcs, plot.arc(arc), sep = "\n")
        }
        
        col.panel.svg <- input.sbgn$printing.info$col.panel.svg
        pathway.name.svg <- input.sbgn$printing.info$pathway.name.svg
        stamp.svg <- input.sbgn$printing.info$stamp.svg
        
        svg.head <- sprintf(svg.header, input.sbgn$svg.dim.x, input.sbgn$svg.dim.y)
        # order of svg code matters: bigger elements (compartment and complex) first
        out.svg <- paste(svg.head, svg.nodes.compartment, svg.nodes.complex, svg.nodes,
                         svg.arcs, svg.cardinality, svg.ports, col.panel.svg, pathway.name.svg,
                         stamp.svg, svg.end, sep = "\n")
        
        Encoding(out.svg) <- "native.enc"  # This is necessary. Some node labels have special symbols that need native encoding
        
        output.file <- input.sbgn$render.sbgn.parameters.list$output.file
        output.svg.file <- paste(output.file, ".svg", sep = "")
        output.formats <- object$output.formats
        write(out.svg, output.svg.file)
        if ("pdf" %in%  output.formats) {
            rsvg::rsvg_pdf(output.svg.file, paste(output.file, ".pdf", sep = ""))
        }
        if ("png" %in% output.formats) {
            rsvg::rsvg_png(output.svg.file, paste(output.file, ".png", sep = ""))
        }
        if ("ps" %in% output.formats) {
            rsvg::rsvg_ps(output.svg.file, paste(output.file, ".ps", sep = ""))
        }
        
        message("Image files written: ", output.file)
        
    } # end main for loop
    
    return(invisible())
}

#########################################################################################################
# function to reconstruct entity (glyph/arc) parameter list. 
# if @parameters.list slot is empty, assign global.parameters.list
# else @parameters.list contains some parameters specified by user 
# which need to be merged with the global.parameters.list. 
# This function used as FUN argument in lapply in print.SBGNview function 
# for input glyphs and arcs list
get.entity.parameter.list <- function(entity, global.parameters.list) {
    if(length(entity@parameters.list) == 0) {
        entity@parameters.list <- global.parameters.list
    } else {
        partial.params.list <- entity@parameters.list
        diff.list.names <- c(setdiff(names(global.parameters.list), names(partial.params.list)))
        # merged parameters.list
        entity@parameters.list <- c(partial.params.list, global.parameters.list[diff.list.names]) 
    }
    return(entity)
}

#########################################################################################################
#' Retrieve output file information from a SBGNview object
#' 
#' @param obj A SBGNview class object.
#' @details This function prints the output file path recorded in a SBGNview object. 
#' @return A string. The output file information. 
#' @examples 
#' \dontrun{
#' outputFile(SBGNview.obj) 
#' }
#' @export
outputFile <- function(obj) {
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
#### old version of print function that calls renderSbgn()
#### updated function uses parsed data in SBGNview object to generate output
# "print.SBGNview" <- function(x, ...) {
#     
#     output.file <- NULL
#     #SBGNview.obj <- x
#     SBGNview.obj <- merge.entity.specific.parameters.list(x)
#     glyphs.arcs.list <- SBGNview.obj$data
#     sbgns <- names(glyphs.arcs.list)
#     for (s in seq_len(length.out = length(glyphs.arcs.list))) {
#         # for each sbgn file
#         data.this.sbgn <- glyphs.arcs.list[[s]]
#         sbgn.parameters.list <- data.this.sbgn$render.sbgn.parameters.list
#         if (!is.null(output.file)) {
#             sbgn.parameters.list$output.file <- gsub(SBGNview.obj$output.file, output.file, 
#                                                      sbgn.parameters.list$output.file)
#         }
#         
#         # if input sbgn-ml file doesn't exist while printing: download file if in ID in pathways.info, else user's file
#         if(!file.exists(sbgn.parameters.list$input.sbgn)) {
#             if(sbgns[s] %in% pathways.info[, "pathway.id"]) {
#                 message(sbgn.parameters.list$input.sbgn, " not in current working directory")
#                 message("Downloading SBGN-ML file for pathway.id: ", sbgns[s])
#                 downloadSbgnFile(pathway.id = sbgns[s])
#             } else {
#                 stop(sbgn.parameters.list$input.sbgn, " not in current working directory.\nPlease make sure SBGN-ML file is in current working directory")
#             }
#         }
#         
#         tp <- renderSbgn(input.sbgn = sbgn.parameters.list$input.sbgn, output.file = sbgn.parameters.list$output.file, 
#                          arcs.info = sbgn.parameters.list$arcs.info, compartment.layer.info = sbgn.parameters.list$compartment.layer.info, 
#                          user.data = sbgn.parameters.list$user.data, output.formats = SBGNview.obj$output.formats, 
#                          sbgn.id.attr = sbgn.parameters.list$sbgn.id.attr, glyphs.user = data.this.sbgn$glyphs.list, 
#                          arcs.user = data.this.sbgn$arcs.list, pathway.name = sbgn.parameters.list$pathway.name, 
#                          global.parameters.list = data.this.sbgn$global.parameters.list)
#         message("Image files written: ", sbgn.parameters.list$output.file, "\n")
#     }
#     return(invisible())
# }

#########################################################################################################
