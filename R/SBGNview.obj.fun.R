

#' Generic function to modify SBGN graph features
#' 
#' This binary operator ("+") modifies a SBGNview object(first argument) using a function (second argument)
#' 
#' @param SBGNview.obj An object of class SBGNview
#' @param fn A character string. The name of any funtion that modifies and returns a SBGNview object. Some functions are available in SBGNview package: \code{\link{highlightPath}}, \code{\link{highlightNodes}}.
#' @return The returned value of *fn*
#' @examples 
#' data(SBGNview.obj)
#' obj.new = SBGNview.obj + 
#'             highlightArcs(class = "production",color = "red") 
#' 
#' @export

"+.SBGNview" <- function(SBGNview.obj,fn){
    fn(SBGNview.obj)
}

#' A wrapper to run function \code{\link{renderSbgn}} for all pathways in a SBGNview object and generate image files.
#' @param x An object of class SBGNview
#' @param ... Other parameters passed to print.
#' @return None
#' @examples 
#' ### use simulated data. Please see vignettes for more examples
#' data("pathways.info","sbgn.xmls")
#' SBGNview.obj = SBGNview(
#'               simulate.data = TRUE
#'               ,sbgn.dir = "./"
#'               ,input.sbgn = "P00001"
#'               
#'               ,output.file = "./test.local.file" 
#'               ,output.formats = c("pdf")
#'               
#'               ,min.gene.value = -1
#'               ,max.gene.value = 1
#'             )
#'  print(SBGNview.obj)
#' @export
"print.SBGNview" = function(x,...
                            # ,output.file = NULL
                            ){
        output.file = NULL
        SBGNview.obj = x
        glyphs.arcs.list = SBGNview.obj$data
        sbgns = names(glyphs.arcs.list)
        for(s in seq_len(length.out = length(glyphs.arcs.list))){ # for each sbgn file
            data.this.sbgn = glyphs.arcs.list[[s]]
            sbgn.parameters.list = data.this.sbgn$render.sbgn.parameters.list
            if(!is.null(output.file)){
                sbgn.parameters.list$output.file = gsub(SBGNview.obj$output.file
                                                        ,output.file
                                                        ,sbgn.parameters.list$output.file )
            }
            tp = renderSbgn(
                    input.sbgn = sbgn.parameters.list$input.sbgn
                    ,output.file = sbgn.parameters.list$output.file
                    ,arcs.info = sbgn.parameters.list$arcs.info
                    ,compartment.layer.info = sbgn.parameters.list$compartment.layer.info
                    ,user.data = sbgn.parameters.list$user.data
                    ,output.formats = sbgn.parameters.list$output.formats
                    ,sbgn.id.attr = sbgn.parameters.list$sbgn.id.attr
                    
                    ,glyphs.user = data.this.sbgn$glyphs.list
                    ,arcs.user = data.this.sbgn$arcs.list
                    ,pathway.name = sbgn.parameters.list$pathway.name
                    
                    ,global.parameters.list = data.this.sbgn$global.parameters.list
                    )
            message("Image files written: ",sbgn.parameters.list$output.file,"\n")
        }
    return(invisible())
}


#' Retrieve output file information from a SBGNview object
#' @param obj A SBGNview object.
#' @details This function prints the output file path recorded in a SBGNview object. When printing the SBGNview object, output image files will be generated using the output file.
#' @return A string. The output file information. The same as parameter "output.file" when running \code{\link{SBGNview}}
#' @examples 
#' data("SBGNview.obj" )
#' outputFile(SBGNview.obj) 
#' @export
outputFile = function(obj){
    obj$outputFile
}


#' Set output file information for a SBGNview object
#' @param obj No need to provide
#' @param value No need to provide 
#' @details This function sets the output file path that can be recorded in a SBGNview object. When printing the SBGNview object, output image files will be generated using the output file.
#' @return A SBGNview object
#' @examples 
#' data("SBGNview.obj" )
#' outputFile(SBGNview.obj) = "./test.output"
#' @export
"outputFile<-" = function(obj,value){
    glyphs.arcs.list = obj$data
    sbgns = names(glyphs.arcs.list)
    for(s in seq_len(length.out = length(glyphs.arcs.list))){ # for each sbgn file
        data.this.sbgn = glyphs.arcs.list[[s]]
        sbgn.parameters.list = data.this.sbgn$render.sbgn.parameters.list
        new.output.file = gsub(obj$output.file
                                                ,value
                                                ,sbgn.parameters.list$output.file )
        obj$data[[s]]$render.sbgn.parameters.list$output.file <- new.output.file
    }
    obj$output.file <- value
    return(obj)
}

