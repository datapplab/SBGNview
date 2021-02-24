
#########################################################################################################
#' Overlay omics data on a SBGN pathway graph and output image files.
#' 
#' This function is not intended to be used directly. Use SBGNview instead. Some input arguments can be better prepared by \code{\link{SBGNview}}.
#' 
#' @param input.sbgn A character string. The path to a local SBGN-ML file.
#' @param output.file,output.formats,sbgn.id.attr These parameters are the same as the ones in \code{\link{SBGNview}}. Please see \code{\link{SBGNview}} for more details.
#' @param glyphs.user  A list, optional. Each element is a 'glyph' object. The element names are glyph IDs (attribute 'id' of XHTML element 'glyph'). Note this is not affected by parameter 'sbgn.id.attr'. The glyph elements contain glyph meta-data for plotting (e.g. text size, border width, border color etc.). Please see the \code{\link{glyph-class}} documentation for more information. By default, SBGNview will run without this argument and return a glyph list extracted from the SBGN file. User can then customize this glyph list and assign it to 'glyphs.user' in the next SBGNview run to update the graph.
#' @param arcs.user  A list, optional. Each member is an 'arc' object. The element names are arc IDs (the value of 'id' attribute in XHTML element 'arc' or 'arc.spline' in the SBGN-ML file). Some SBGN-ML files have no arc IDs, in this case SBGNview  will create an arc ID using 'source' and 'target' node IDs). The arc object contains arc meta-data for plotting (e.g. arc line width, line color etc.). Please see the \code{\link{arc-class}} documentation for more information. By default, SBGNview() will run without this argument and return an arc list. User can then customize this arc list and assign it to 'arcs.user' in the next SBGNview() run to update the arcs.
#' @param arcs.info A character string. It should be one of the following: 'parse splines', 'straight' or a string of svg code of arcs. If it is 'parse splines', this function will look for XML element 'arc.spline' in the SBGN-ML file and plot spline arcs. If it is 'straight', the function will look for element 'arc' and plot straight line arcs. If it is a string of svg code, it will write this code directly to the output svg file.
#' @param compartment.layer.info A character vector. It is a vector containing the IDs of all compartment glyphs. It determins the layer arrangement of compartments. Compartments will be drawn following their sequence in this vector. Therefore, a compartment that appears later in the vector will be on the front layer and covers the compartments that are before it in this vector. This is important. In some cases compartments have overlap. This layer information ensures that a glyph laying in the overlapped region belongs to the compartment on the top layer.
#' @param user.data A list. It holds both gene/protein data and compound data. Names are gene or compounds, each element is a numeric vector of the omics data of each molecule. 
#' @param color.panel.scale Numeric. Default: 1. It controls the relative size of color scheme panel. 
#' @param key.pos  A character string. The position of color panel. Default: 'topright'. Accepts one of 'bottomleft' , 'bottomright', 'topright', or 'topleft'. The ideal position for the color panel is 'topright' or 'bottomright'. If 'topleft' or 'bottomleft' are passed in, the key.pos location will default to 'topright'. If key.pos is set to 'none', the pathway name and color panel won't be plotted.
#' @param if.plot.cardinality Logical. Default: F. If plot cardinality glyphs.
#' @param status.node.font.scale Numeric. Default: 3. Scale the font size for status variable and unit of information nodes.
#' @param font.size.scale.complex Numeric. Default: 1.1. Scale the font size of a complex. 
#' @param font.size.scale.compartment Numeric. Default: 1.6. Scale the font size of a compartment. 
#' @param text.length.factor Numeric. Default: 2. How wide the wrapped text should be. If text is longer than the width controled by this parameter, the text is split into a new line, but only at characters in 'label.spliting.string'. Controls all glyphs.
#' @param text.length.factor.macromolecule  Numeric. Default: 2. Used to determine label text wrapping based on number of characters, font size, and node width for macromolecule glyphs. 
#' @param text.length.factor.compartment Numeric. Default: 2. Used to determine label text wrapping based on number of characters, font size, and node width for compartment glyphs.
#' @param text.length.factor.complex  Numeric. Default: 2. Used to determine label text wrapping based on number of characters, font size, and node width for complex glyphs. 
#' @param global.parameters.list List. A record of parameters fed to 'renderSbgn' for reuse. It will over-write other parameters. It is not designed to be used directly.
#' @param color.panel.n.grid Numeric. Default: 21. How many colors does the color scheme show.
#' @param col.gene.low  A character string. Default: 'green'. Color panel's color representing low gene data values. 
#' @param col.gene.high A character string. Default: 'red'. Color panel's color representing high gene data values. 
#' @param col.gene.mid A character string. Default: 'gray'. Color panel's color representing mid range gene data values.
#' @param col.cpd.low  A character string. Default: 'blue'. Color panel's color representing low compound data values. 
#' @param col.cpd.high A character string. Default: 'yellow'. Color panel's color representing high compound data values. 
#' @param col.cpd.mid A character string. Default: 'gray'. Color panel's color representing mid range compound data values. 
#' @param min.gene.value Numeric. Default: -1. Color panel's min value for gene data. Values smaller than this will have the same color as min value.
#' @param max.gene.value Numeric. Default: 1. Color panel's max value for gene data. Values greater than this will have the same color as the max value.
#' @param mid.gene.value Numeric. Default: 0. Color panel's mid value for gene data. 
#' @param min.cpd.value  Numeric. Default: -1. Color panel's min value for compound data. Values smaller than this will have the same color as min value.
#' @param max.cpd.value Numeric. Default: 1. Color panel's max value for compound data. Values greater than this will have the same color as max value.
#' @param mid.cpd.value Numeric. Default: 0. Color panel's mid value for compound data.
#' @param multimer.margin  Numeric. Default: 5. For multimers, they are represented by two partly overlapped shapes (rectangle, ellipse etc.). This parameter controls how much the two shapes overlap. 
#' @param compartment.opacity Numeric. Default: 1. Controls how transparent the compartments are.
#' @param auxiliary.opacity Numeric. Default: 1.  Controls opacity of auxiliary glyphs.
#' @param if.plot.annotation.nodes Logical. Default: F. Some SBGN-ML files have 'annotation' glyphs. By default we don't plot them.
#' @param inhibition.edge.end.shift Numeric. Default: 5. The tip of 'inhibition' arcs is a line segment. Sometimes it overlaps with target glyph's border. We can shift it along the arc to prevent the overlap.
#' @param edge.tip.size Numeric. Default: 6. Control size of edge tips. 
#' @param if.use.number.for.long.label  Logical. Default: F. If the label is too long, we can create a shorter name for it. e.g. 'macromolecule_1'.
#' @param if.write.shorter.label.mapping Logical. Default: T. If if.use.number.for.long.label is 'T', we can write the mapping between shorter name and the original label to a text file.
#' @param label.spliting.string  A character vector. Default: c(' ','-',';','/','_'). When we split text into multiple lines, these characters will be used to split label (where a new line can be created). 
#' @param complex.compartment.label.margin Numeric. Default: 8. Move the label text vertically for compartment and complex.
#' @param font.size Numeric. Default: 3. Affects font size of all types of glyphs.
#' @param font.size.scale.gene Numeric. Default: 3. Scales font size according to the node's width for large compartments. Only affect font size of 'macromolecule' glyphs.
#' @param font.size.scale.cpd Numeric. Default: 3. Scales font size according to the node's width for large compartments. Only affects font size of 'simple chemical' glyphs.
#' @param logic.node.font.scale Numeric. Default: 3. Controls the size of logical glyphs ('and', 'or', 'not' etc.).
#' @param node.width.adjust.factor Numeric. Default: 2. Change font size according to the width of its glyph. If the glyph is too large (e.g. compartment), its label may look too small. Then we can enlarge the label in proportion to the width of the glyph. It affects all types of glyphs. 
#' @param pathway.name List containing two elements: 1. pathway name 2. stamp information. See argument description in \code{\link{SBGNview}} function to change/update pathway name displayed on graph. 
#' @param pathway.name.font.size Numeric. Default: 1. When pathway names are plotted on graph, this parameter controls their font size.
#' @param if.scale.complex.font.size Logical. Default: F.  Whether to scale complex font size according to its width. If set to 'T', the 'node.width.adjust.factor.complex' argument can be used to specify the scale factor. 
#' @param if.scale.compartment.font.size  Logical. Default: F. Whether to scale compartment font size according to its width. If set to 'T', the 'node.width.adjust.factor.compartment' argument can be used to specify the scale factor.   
#' @param node.width.adjust.factor.compartment Numeric. Default: 0.02. How much the font size should change in proportion to the width of compartment. The font is scaled only if 'if.scale.compartment.font.size' is set to 'T'. To find the best scale factor which works you, start with 0.02 (default) and incrementally increase that value.  
#' @param node.width.adjust.factor.complex Numeric. Default: 0.02. How much the font size should change in proportion to the width of complex. The font is scaled only if 'if.scale.complex.font.size' is set to 'T'. To find the best scale factor which works you, start with 0.02 (default) and incrementally increase that value. 
#' @param space.between.color.panel.and.entity Numeric. Default: 100. The minimum space between color panel and any other object in the graph. The function will always try to find a location of the color panel to minimize empty space on the whole graph. This parameter controls how close it can reach a glyph.
#' @param if.plot.svg Logical. Default: T. Whether to generate svg code or only parse SBGN-ML file. This parameter is for internal use only.
#' @return A list of three elements: glyphs.list, arcs.list, global.parameters.list
#' @examples
#' \dontrun{
#' data(pathways.info)
#' SBGNview.obj <- SBGNview(simulate.data = TRUE,
#'                          sbgn.dir = './',
#'                          input.sbgn = 'P00001',
#'                          output.file = './test.local.file',
#'                          output.formats = c('pdf'),
#'                          min.gene.value = -1,
#'                          max.gene.value = 1)
#'  }
#' 
#' @export

# updated so function doesn't write files. only parses data which will be added to SBGNview object
# parsed data will be used by print.SBGNview function to write output files
renderSbgn <- function(input.sbgn, output.file, output.formats, sbgn.id.attr, 
                       glyphs.user = list(), arcs.user = list(), arcs.info = "straight", 
                       compartment.layer.info = "original", user.data = matrix("no.user.data", nrow = 1), 
                       if.plot.svg = TRUE, key.pos = "topright", color.panel.scale = 1,  # Control the relative size of color scheme panel
                       color.panel.n.grid = 21,  # how many colors doese the color scheme show
                       col.gene.low = "green", col.gene.high = "red", col.gene.mid = "gray", col.cpd.low = "blue", col.cpd.high = "yellow", 
                       col.cpd.mid = "gray", min.gene.value = -1,  # color panel min value, values smaller than this will have the min.value color
                       max.gene.value = 1, mid.gene.value = 0, min.cpd.value = -1,  # color panel min value, values smaller than this will have the min.value color
                       max.cpd.value = 1, mid.cpd.value = 0, pathway.name = "", pathway.name.font.size = 1, 
                       if.plot.cardinality = FALSE, multimer.margin = 5, compartment.opacity = 1,  # how transparent the compartments are
                       auxiliary.opacity = 1,  # opacity of auxiliary nodes
                       if.plot.annotation.nodes = FALSE,  # Some sbgn files have 'annotation' nodes. By default we don't plot them
                       inhibition.edge.end.shift = 5,  # The tip of 'inhibition' arcs is a line segment. Sometimes it overlaps with target node's border. We can shift it to prevent the overlap.
                       edge.tip.size = 6, if.use.number.for.long.label = FALSE, 
                       label.spliting.string = c(" ", ":", "-", ";", "/", "_"),  # the regular expression used to spline text to wrape labels. Can be set to '' to split by single letter. The default is space ' '. In some cases the word seperated by ' ' is too long. We can use space or '-'(i.e. '-| ')  to  split the words
                       complex.compartment.label.margin = 8,  # shift the label to the upper direction
                       if.write.shorter.label.mapping = TRUE, font.size = 3, logic.node.font.scale = 3, 
                       status.node.font.scale = 3, node.width.adjust.factor = 2,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
                       font.size.scale.gene = 3,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
                       font.size.scale.cpd = 3, # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
                       font.size.scale.complex = 1.1, font.size.scale.compartment = 1.6, if.scale.complex.font.size = FALSE, 
                       if.scale.compartment.font.size = FALSE, node.width.adjust.factor.compartment = 0.02,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
                       node.width.adjust.factor.complex = 0.02,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
                       text.length.factor = 2, text.length.factor.macromolecule = 2,  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for how wide the wrapped text should be.
                       text.length.factor.compartment = 2, text.length.factor.complex = 2, space.between.color.panel.and.entity = 100, 
                       global.parameters.list = NULL) {
    
    col.panel.params <- find.col.panel.range(user.data, max.gene.value, mid.gene.value, 
                                             min.gene.value, max.cpd.value, mid.cpd.value, min.cpd.value)
    
    if.has.gene.data <- col.panel.params$if.has.gene.data
    if.has.cpd.data <- col.panel.params$if.has.cpd.data
    max.gene.value <- col.panel.params$max.gene.value
    mid.gene.value <- col.panel.params$mid.gene.value
    min.gene.value <- col.panel.params$min.gene.value
    max.cpd.value <- col.panel.params$max.cpd.value
    mid.cpd.value <- col.panel.params$mid.cpd.value
    min.cpd.value <- col.panel.params$min.cpd.value
    
    if (is.null(global.parameters.list)) {
        global.parameters.list <- list()
        global.parameters.list$if.plot.cardinality <- if.plot.cardinality
        global.parameters.list$multimer.margin <- multimer.margin
        global.parameters.list$if.write.shorter.label.mapping <- if.write.shorter.label.mapping
        global.parameters.list$compartment.opacity <- compartment.opacity  # how transparent the compartments are
        global.parameters.list$auxiliary.opacity <- auxiliary.opacity  # opacity of auxiliary nodes
        global.parameters.list$if.plot.annotation.nodes <- if.plot.annotation.nodes  # Some sbgn files have 'annotation' nodes. By default we don't plot them
        
        # arc parameters
        global.parameters.list$inhibition.edge.end.shift <- inhibition.edge.end.shift  # The tip of 'inhibition' arcs is a line segment. Sometimes it overlaps with target node's border. We can shift it to prevent the overlap.
        global.parameters.list$edge.tip.size <- edge.tip.size
        
        # label parameters
        global.parameters.list$if.use.number.for.long.label <- if.use.number.for.long.label
        global.parameters.list$label.spliting.string <- label.spliting.string  # the regular expression used to spline text to wrape labels. Can be set to '' to split by single letter. The default is space ' '. In some cases the word seperated by ' ' is too long. We can use space or '-'(i.e. '-| ')  to  split the words
        global.parameters.list$complex.compartment.label.margin <- complex.compartment.label.margin  # shift the label to the upper direction
        global.parameters.list$font.size <- font.size
        global.parameters.list$logic.node.font.scale <- logic.node.font.scale
        global.parameters.list$status.node.font.scale <- status.node.font.scale
        global.parameters.list$font.size.scale.gene <- font.size.scale.gene  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$font.size.scale.cpd <- font.size.scale.cpd  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$font.size.scale.compartment <- font.size.scale.compartment  # 
        global.parameters.list$font.size.scale.complex <- font.size.scale.complex  # 
        
        global.parameters.list$node.width.adjust.factor <- node.width.adjust.factor  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$if.scale.compartment.font.size <- if.scale.compartment.font.size  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$if.scale.complex.font.size <- if.scale.complex.font.size  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$node.width.adjust.factor.compartment <- node.width.adjust.factor.compartment  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$node.width.adjust.factor.complex <- node.width.adjust.factor.complex  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
        global.parameters.list$text.length.factor <- text.length.factor  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
        global.parameters.list$pathway.name <- pathway.name
        global.parameters.list$pathway.name.font.size <- pathway.name.font.size
        global.parameters.list$text.length.factor.macromolecule <- text.length.factor.macromolecule  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
        global.parameters.list$text.length.factor.compartment <- text.length.factor.compartment  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
        global.parameters.list$text.length.factor.complex <- text.length.factor.complex  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
        
        global.parameters.list$key.pos <- key.pos  # ll , lr, ur, ul  # location of color panel: lower left, lower right, upper right, upper left
        global.parameters.list$color.panel.scale <- color.panel.scale  # Control the relative size of color scheme panel
        global.parameters.list$color.panel.n.grid <- color.panel.n.grid  # how many colors does the color scheme show
        global.parameters.list$col.gene.low <- col.gene.low
        global.parameters.list$col.gene.high <- col.gene.high
        global.parameters.list$col.gene.mid <- col.gene.mid
        
        global.parameters.list$col.cpd.low <- col.cpd.low
        global.parameters.list$col.cpd.high <- col.cpd.high
        global.parameters.list$col.cpd.mid <- col.cpd.mid
        global.parameters.list$min.gene.value <- min.gene.value  # color panel min value, values smaller than this will have the min.value color
        global.parameters.list$max.gene.value <- max.gene.value
        global.parameters.list$mid.gene.value <- mid.gene.value
        global.parameters.list$min.cpd.value <- min.cpd.value  # color panel min value, values smaller than this will have the min.value color
        global.parameters.list$max.cpd.value <- max.cpd.value
        global.parameters.list$mid.cpd.value <- mid.cpd.value
        global.parameters.list$space.between.color.panel.and.entity <- space.between.color.panel.and.entity
    }
    
    sbgn.xml <- read_xml(input.sbgn)
    xml_attrs(sbgn.xml) <- NULL  # Remove root node attribute. This is necessary Otherwise xml2 won't find the nodes when using xml_find_all.
    
    message("checking graph size and create margin for color panel")
    coords.range.list <- find.max.xy(sbgn.xml, arcs.info, color.panel.scale, global.parameters.list)
    max.x <- coords.range.list$max.xw
    max.y <- coords.range.list$max.yh
    min.x <- coords.range.list$min.x
    min.y <- coords.range.list$min.y
    y.margin <- coords.range.list$y.margin
    
    message("parsing ports")
    ports <- xml.to.port.glyphs(sbgn.xml, y.margin = y.margin)  # The output 'ports' is a list of port glphs
    
    message("parsing glyphs")
    parse.glyph.out.list <- parse.glyph(sbgn.xml, user.data, y.margin = y.margin, max.x = max.x, 
                                        global.parameters.list = global.parameters.list, 
                                        sbgn.id.attr = sbgn.id.attr, if.plot.svg = if.plot.svg, 
                                        glyphs.user = glyphs.user, 
                                        compartment.layer.info = compartment.layer.info, 
                                        if.plot.cardinality = if.plot.cardinality)
    glyphs <- parse.glyph.out.list$glyphs
    # find plot parameters
    min.w <- parse.glyph.out.list$min.w  # find the minimum h to set the text size
    min.w.maxNchar <- parse.glyph.out.list$min.w.maxNchar  # number of characters of the text that occupies the most narrow box
    # svg contents
    svg.ports <- parse.glyph.out.list$svg.ports
    svg.nodes <- parse.glyph.out.list$svg.nodes
    svg.nodes.complex <- parse.glyph.out.list$svg.nodes.complex
    svg.nodes.compartment <- parse.glyph.out.list$svg.nodes.compartment
    svg.cardinality <- parse.glyph.out.list$svg.cardinality  # the cadinality node are supposed to be in front of the arcs, so need to print it again at the end of the svg file
    shorter.label.mapping.list <- parse.glyph.out.list$shorter.label.mapping.list
    
    if (if.write.shorter.label.mapping & if.use.number.for.long.label & !is.vector(shorter.label.mapping.list)) {
        write.table(shorter.label.mapping.list, paste(output.file, ".shorter.label.mapping.tsv", 
                                                      sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t")
    }
    # combine glyphs and ports
    glyphs <- c(glyphs, ports)
    
    message("parsing arcs")
    arcs.result <- get.arcs(arcs.info, sbgn.xml, glyphs, if.plot.svg, y.margin, 
                            global.parameters.list, arcs.user)
    svg.arc <- arcs.result$svg.arc
    arcs.list <- arcs.result$arcs.list
    
    message("plotting color panel")
    col.panel.params <- find.col.panel.position.and.plot(y.margin, global.parameters.list, 
                                                         if.has.gene.data, if.has.cpd.data, 
                                                         parse.glyph.out.list, max.x, max.y, 
                                                         min.x, min.y)
    col.panel.svg <- col.panel.params$col.panel.svg
    col.panel.w <- col.panel.params$col.panel.w
    col.panel.y <- col.panel.params$col.panel.y
    col.panel.h <- col.panel.params$col.panel.h
    
    # add pathway.name and stamp
    stamp.svg.list <- add.stamp(col.panel.w, col.panel.y, global.parameters.list, 
                                template.text, template.text.pathway.name, min.x, 
                                max.x, max.y, y.margin)
    pathway.name.svg <- stamp.svg.list$pathway.name.svg
    stamp.svg <- stamp.svg.list$stamp.svg
    
    # generate output xml content
    # svg.dim.x = max.x+50+4*70
    # svg.dim.y = max.y+col.panel.h+50+y.margin
    svg.dim.x = max.x + col.panel.w/2 
    remove.margin <- sqrt(y.margin)#*2
    if(if.has.gene.data & if.has.cpd.data | color.panel.scale != 1) {
        remove.margin <- remove.margin + col.panel.h/2
    }
    svg.dim.y <- max.y + col.panel.h + y.margin - remove.margin
    
    # including svg code used by print function
    printing.info <- list(col.panel.svg = col.panel.svg, pathway.name.svg = pathway.name.svg,
                          stamp.svg = stamp.svg)
    
    return(list(glyphs.list = glyphs, arcs.list = arcs.list, global.parameters.list = global.parameters.list,
                coords.range.list = coords.range.list, svg.dim.x = svg.dim.x, svg.dim.y = svg.dim.y,
                printing.info = printing.info))
}

#########################################################################################################
# custom template to print label on top left corner in output files as two lines
# line 1: pathway name
# line 2: (database::id)
template.text.pathway.name <- "
<text x=\"%f\" y=\"%f\" id=\"%s\"
style=\"font-size:%fpx;text-anchor:%s;font-family:Arial,Helvetica;alignment-baseline:%s;stroke.opacity:%s;dominant-baseline:%s\"
fill=\"%s\">
<tspan x=\"%f\">%s</tspan>
<tspan x=\"%f\" dy=\"%f\">%s</tspan>
</text>
"

add.stamp <- function(col.panel.w, col.panel.y, global.parameters.list, template.text, 
                      template.text.pathway.name, min.x, max.x, max.y, y.margin) {
    
    pathway.name.font <- col.panel.w/global.parameters.list$color.panel.scale/7 * global.parameters.list$pathway.name.font.size
    pathway.name.y <- col.panel.y #- pathway.name.font
    
    pathway.name.display <- global.parameters.list$pathway.name$pathway.name.on.graph # pathwayname::databse::id
    # split pathway.name.display into 1) pathway name and 2) database :: id
    split.name.db <- regmatches(pathway.name.display, regexpr("::", pathway.name.display), invert = TRUE)
    name.of.pathway <- split.name.db[[1]][1]
    db.and.id <- paste("(", split.name.db[[1]][2], ")", sep = "")
    
    # if db.and.id = 'user.data' and name.of.pathway not in sbgn.xmls, or not in pathways.info[,'pathway.name']
    # display nothing for top left stamp, if in sbgn.xmls or pathways.info[,'pathway.name'], show first line not second line
    if(split.name.db[[1]][2] == "user.data") {
        if(name.of.pathway %in% names(sbgn.xmls) || 
           name.of.pathway %in% pathways.info[,'pathway.name']) {
            db.and.id <- ""
        } else {
            name.of.pathway <- ""
            db.and.id <- ""
        }
    } else if (split.name.db[[1]][2] == "user.named.pathway") {
        db.and.id <- ""
    }
    
    # Using template.text.pathway.name to display in two lines
    dy <- (pathway.name.font / 2) + (10*global.parameters.list$pathway.name.font.size) # space between two lines
    pathway.name.svg <- sprintf(template.text.pathway.name, min.x + 10, pathway.name.y, "pathway.name",
                                pathway.name.font, "left", "baseline", 1, "baseline", "black", 
                                min.x + 10, name.of.pathway, min.x + 10, dy, db.and.id)
    
    stamp.if.sbgnhub.h <- max.x/5/9
    stamp.y <- max(max.y, col.panel.y + col.panel.w) + stamp.if.sbgnhub.h + y.margin
    stamp.svg <- sprintf(template.text, 20, stamp.y, "stamp", 22,
                         "left", "baseline", 1, "baseline", "black", global.parameters.list$pathway.name$if.file.in.collection)
    return(list(pathway.name.svg = pathway.name.svg, stamp.svg = stamp.svg))
}


#########################################################################################################
# get information of arcs for plotting and list arc objects 
get.arcs <- function(arcs.info, sbgn.xml, glyphs, if.plot.svg, y.margin, 
                     global.parameters.list, arcs.user) {
    
    if (arcs.info == "straight") {
        print("using original edges")
        result.list <- parse.arcs(sbgn.xml, glyphs, if.plot.svg = if.plot.svg, y.margin = y.margin, 
            global.parameters.list = global.parameters.list, arcs.user = arcs.user)
        svg.arc <- result.list$svg.arc
        arcs.list <- result.list$arcs.list
    } else if (arcs.info == "parse splines") {
        message("using spline arcs")
        result.list <- parse.splines(sbgn.xml, glyphs, if.plot.svg = if.plot.svg, 
            y.margin = y.margin, global.parameters.list = global.parameters.list, 
            arcs.user = arcs.user)
        svg.arc <- result.list$svg.splines
        arcs.list <- result.list$splines.list
    } else {
        message("using spline svg")
        svg.arc <- arcs.info
        arcs.list <- list()
    }
    return(list(svg.arc = svg.arc, arcs.list = arcs.list))
}

#########################################################################################################
### old version of renderSbgn() that parses output svg code and writes files
### updated version parses all data which is added to SBGNview object

## removed if.write.files argument in updated version
# @param if.write.files Logical. Default: T. If generate image files. This parameter is for internal use only.

# renderSbgn <- function(input.sbgn, output.file, if.write.files = TRUE, output.formats, 
#                        sbgn.id.attr, glyphs.user = list(), arcs.user = list(), arcs.info = "straight", 
#                        compartment.layer.info = "original", user.data = matrix("no.user.data", nrow = 1), 
#                        if.plot.svg = TRUE, key.pos = "topright", color.panel.scale = 1,  # Control the relative size of color scheme panel
#                        color.panel.n.grid = 21,  # how many colors doese the color scheme show
#                        col.gene.low = "green", col.gene.high = "red", col.gene.mid = "gray", col.cpd.low = "blue", col.cpd.high = "yellow", 
#                        col.cpd.mid = "gray", min.gene.value = -1,  # color panel min value, values smaller than this will have the min.value color
#                        max.gene.value = 1, mid.gene.value = 0, min.cpd.value = -1,  # color panel min value, values smaller than this will have the min.value color
#                        max.cpd.value = 1, mid.cpd.value = 0, pathway.name = "", pathway.name.font.size = 1, 
#                        if.plot.cardinality = FALSE, multimer.margin = 5, compartment.opacity = 1,  # how transparent the compartments are
#                        auxiliary.opacity = 1,  # opacity of auxiliary nodes
#                        if.plot.annotation.nodes = FALSE,  # Some sbgn files have 'annotation' nodes. By default we don't plot them
#                        inhibition.edge.end.shift = 5,  # The tip of 'inhibition' arcs is a line segment. Sometimes it overlaps with target node's border. We can shift it to prevent the overlap.
#                        edge.tip.size = 6, if.use.number.for.long.label = FALSE, 
#                        label.spliting.string = c(" ", ":", "-", ";", "/", "_"),  # the regular expression used to spline text to wrape labels. Can be set to '' to split by single letter. The default is space ' '. In some cases the word seperated by ' ' is too long. We can use space or '-'(i.e. '-| ')  to  split the words
#                        complex.compartment.label.margin = 8,  # shift the label to the upper direction
#                        if.write.shorter.label.mapping = TRUE, font.size = 3, logic.node.font.scale = 3, 
#                        status.node.font.scale = 3, node.width.adjust.factor = 2,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#                        font.size.scale.gene = 3,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#                        font.size.scale.cpd = 3, # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#                        font.size.scale.complex = 1.1, font.size.scale.compartment = 1.6, if.scale.complex.font.size = FALSE, 
#                        if.scale.compartment.font.size = FALSE, node.width.adjust.factor.compartment = 0.02,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#                        node.width.adjust.factor.complex = 0.02,  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#                        text.length.factor = 2, text.length.factor.macromolecule = 2,  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for how wide the wrapped text should be.
#                        text.length.factor.compartment = 2, text.length.factor.complex = 2, space.between.color.panel.and.entity = 100, 
#                        global.parameters.list = NULL) {
#     
#     col.panel.params <- find.col.panel.range(user.data, max.gene.value, mid.gene.value, 
#                                              min.gene.value, max.cpd.value, mid.cpd.value, min.cpd.value)
#     
#     if.has.gene.data <- col.panel.params$if.has.gene.data
#     if.has.cpd.data <- col.panel.params$if.has.cpd.data
#     max.gene.value <- col.panel.params$max.gene.value
#     mid.gene.value <- col.panel.params$mid.gene.value
#     min.gene.value <- col.panel.params$min.gene.value
#     max.cpd.value <- col.panel.params$max.cpd.value
#     mid.cpd.value <- col.panel.params$mid.cpd.value
#     min.cpd.value <- col.panel.params$min.cpd.value
#     
#     if (is.null(global.parameters.list)) {
#         global.parameters.list <- list()
#         global.parameters.list$if.plot.cardinality <- if.plot.cardinality
#         global.parameters.list$multimer.margin <- multimer.margin
#         global.parameters.list$if.write.shorter.label.mapping <- if.write.shorter.label.mapping
#         global.parameters.list$compartment.opacity <- compartment.opacity  # how transparent the compartments are
#         global.parameters.list$auxiliary.opacity <- auxiliary.opacity  # opacity of auxiliary nodes
#         global.parameters.list$if.plot.annotation.nodes <- if.plot.annotation.nodes  # Some sbgn files have 'annotation' nodes. By default we don't plot them
#         
#         # arc parameters
#         global.parameters.list$inhibition.edge.end.shift <- inhibition.edge.end.shift  # The tip of 'inhibition' arcs is a line segment. Sometimes it overlaps with target node's border. We can shift it to prevent the overlap.
#         global.parameters.list$edge.tip.size <- edge.tip.size
#         
#         # label parameters
#         global.parameters.list$if.use.number.for.long.label <- if.use.number.for.long.label
#         global.parameters.list$label.spliting.string <- label.spliting.string  # the regular expression used to spline text to wrape labels. Can be set to '' to split by single letter. The default is space ' '. In some cases the word seperated by ' ' is too long. We can use space or '-'(i.e. '-| ')  to  split the words
#         global.parameters.list$complex.compartment.label.margin <- complex.compartment.label.margin  # shift the label to the upper direction
#         global.parameters.list$font.size <- font.size
#         global.parameters.list$logic.node.font.scale <- logic.node.font.scale
#         global.parameters.list$status.node.font.scale <- status.node.font.scale
#         global.parameters.list$font.size.scale.gene <- font.size.scale.gene  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$font.size.scale.cpd <- font.size.scale.cpd  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$font.size.scale.compartment <- font.size.scale.compartment  # 
#         global.parameters.list$font.size.scale.complex <- font.size.scale.complex  # 
#         
#         global.parameters.list$node.width.adjust.factor <- node.width.adjust.factor  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$if.scale.compartment.font.size <- if.scale.compartment.font.size  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$if.scale.complex.font.size <- if.scale.complex.font.size  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$node.width.adjust.factor.compartment <- node.width.adjust.factor.compartment  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$node.width.adjust.factor.complex <- node.width.adjust.factor.complex  # change font size according to the node's width, for large compartments, it is better to enlarge font size in proportion to the width
#         global.parameters.list$text.length.factor <- text.length.factor  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
#         global.parameters.list$pathway.name <- pathway.name
#         global.parameters.list$pathway.name.font.size <- pathway.name.font.size
#         global.parameters.list$text.length.factor.macromolecule <- text.length.factor.macromolecule  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
#         global.parameters.list$text.length.factor.compartment <- text.length.factor.compartment  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
#         global.parameters.list$text.length.factor.complex <- text.length.factor.complex  # calculate label length based on number of characters and font size, then the length is compared to the width of node for text wrapping. This factor is used as a scale factor for the label length and node width to match
#         
#         global.parameters.list$key.pos <- key.pos  # ll , lr, ur, ul  # location of color panel: lower left, lower right, upper right, upper left
#         global.parameters.list$color.panel.scale <- color.panel.scale  # Control the relative size of color scheme panel
#         global.parameters.list$color.panel.n.grid <- color.panel.n.grid  # how many colors does the color scheme show
#         global.parameters.list$col.gene.low <- col.gene.low
#         global.parameters.list$col.gene.high <- col.gene.high
#         global.parameters.list$col.gene.mid <- col.gene.mid
#         
#         global.parameters.list$col.cpd.low <- col.cpd.low
#         global.parameters.list$col.cpd.high <- col.cpd.high
#         global.parameters.list$col.cpd.mid <- col.cpd.mid
#         global.parameters.list$min.gene.value <- min.gene.value  # color panel min value, values smaller than this will have the min.value color
#         global.parameters.list$max.gene.value <- max.gene.value
#         global.parameters.list$mid.gene.value <- mid.gene.value
#         global.parameters.list$min.cpd.value <- min.cpd.value  # color panel min value, values smaller than this will have the min.value color
#         global.parameters.list$max.cpd.value <- max.cpd.value
#         global.parameters.list$mid.cpd.value <- mid.cpd.value
#         global.parameters.list$space.between.color.panel.and.entity <- space.between.color.panel.and.entity
#     }
#     
#     sbgn.xml <- read_xml(input.sbgn)
#     xml_attrs(sbgn.xml) <- NULL  # Remove root node attribute. This is necessary Otherwise xml2 won't find the nodes when using xml_find_all.
#     
#     message("checking graph size and create margin for color panel")
#     coords.range.list <- find.max.xy(sbgn.xml, arcs.info, color.panel.scale, global.parameters.list)
#     max.x <- coords.range.list$max.xw
#     max.y <- coords.range.list$max.yh
#     min.x <- coords.range.list$min.x
#     min.y <- coords.range.list$min.y
#     y.margin <- coords.range.list$y.margin
#     
#     message("parsing ports")
#     ports <- xml.to.port.glyphs(sbgn.xml, y.margin = y.margin)  # The output 'ports' is a list of port glphs
#     
#     message("parsing glyphs")
#     parse.glyph.out.list <- parse.glyph(sbgn.xml, user.data, y.margin = y.margin, max.x = max.x, 
#                                         global.parameters.list = global.parameters.list, 
#                                         sbgn.id.attr = sbgn.id.attr, if.plot.svg = if.plot.svg, 
#                                         glyphs.user = glyphs.user, 
#                                         compartment.layer.info = compartment.layer.info, 
#                                         if.plot.cardinality = if.plot.cardinality)
#     glyphs <- parse.glyph.out.list$glyphs
#     # find plot parameters
#     min.w <- parse.glyph.out.list$min.w  # find the minimum h to set the text size
#     min.w.maxNchar <- parse.glyph.out.list$min.w.maxNchar  # number of characters of the text that occupies the most narrow box
#     # svg contents
#     svg.ports <- parse.glyph.out.list$svg.ports
#     svg.nodes <- parse.glyph.out.list$svg.nodes
#     svg.nodes.complex <- parse.glyph.out.list$svg.nodes.complex
#     svg.nodes.compartment <- parse.glyph.out.list$svg.nodes.compartment
#     svg.cardinality <- parse.glyph.out.list$svg.cardinality  # the cadinality node are supposed to be in front of the arcs, so need to print it again at the end of the svg file
#     shorter.label.mapping.list <- parse.glyph.out.list$shorter.label.mapping.list
#     
#     if (if.write.shorter.label.mapping & if.use.number.for.long.label & !is.vector(shorter.label.mapping.list)) {
#         write.table(shorter.label.mapping.list, paste(output.file, ".shorter.label.mapping.tsv", 
#                                                       sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t")
#     }
#     # combine glyphs and ports
#     glyphs <- c(glyphs, ports)
#     
#     message("parsing arcs")
#     arcs.result <- get.arcs(arcs.info, sbgn.xml, glyphs, if.plot.svg, y.margin, 
#                             global.parameters.list, arcs.user)
#     svg.arc <- arcs.result$svg.arc
#     arcs.list <- arcs.result$arcs.list
#     
#     message("plotting color panel")
#     col.panel.params <- find.col.panel.position.and.plot(y.margin, global.parameters.list, 
#                                                          if.has.gene.data, if.has.cpd.data, 
#                                                          parse.glyph.out.list, max.x, max.y, 
#                                                          min.x, min.y)
#     col.panel.svg <- col.panel.params$col.panel.svg
#     col.panel.w <- col.panel.params$col.panel.w
#     col.panel.y <- col.panel.params$col.panel.y
#     col.panel.h <- col.panel.params$col.panel.h
#     
#     # add pathway.name and stamp
#     stamp.svg.list <- add.stamp(col.panel.w, col.panel.y, global.parameters.list, 
#                                 template.text, template.text.pathway.name, min.x, 
#                                 max.x, max.y, y.margin)
#     pathway.name.svg <- stamp.svg.list$pathway.name.svg
#     stamp.svg <- stamp.svg.list$stamp.svg
#     
#     # generate output xml content
#     # svg.dim.x = max.x+50+4*70
#     # svg.dim.y = max.y+col.panel.h+50+y.margin
#     svg.dim.x = max.x + col.panel.w/2 
#     remove.margin <- sqrt(y.margin)#*2
#     if(if.has.gene.data & if.has.cpd.data | color.panel.scale != 1) {
#         remove.margin <- remove.margin + col.panel.h/2
#     }
#     svg.dim.y = max.y + col.panel.h + y.margin - remove.margin
#     svg.header = sprintf(svg.header, svg.dim.x, svg.dim.y)
#     
#     out <- paste(svg.header, svg.nodes.compartment, svg.nodes.complex, svg.nodes,
#                  svg.arc, svg.cardinality, svg.ports, col.panel.svg, pathway.name.svg, 
#                  stamp.svg, svg.end, sep = "\n")
#     
#     Encoding(out) <- "native.enc"  # This is necessary. Some node labels have special symbols that need native encoding
#     
#     # write output file
#     if (if.write.files) {
#         output.svg.file <- paste(output.file, ".svg", sep = "")
#         write(out, output.svg.file)
#         if ("pdf" %in% output.formats) {
#             rsvg::rsvg_pdf(output.svg.file, paste(output.file, ".pdf", sep = ""))
#         }
#         if ("png" %in% output.formats) {
#             rsvg::rsvg_png(output.svg.file, paste(output.file, ".png", sep = ""))
#         }
#         if ("ps" %in% output.formats) {
#             rsvg::rsvg_ps(output.svg.file, paste(output.file, ".ps", sep = ""))
#         }
#     }
#     
#     return(list(glyphs.list = glyphs, arcs.list = arcs.list, global.parameters.list = global.parameters.list,
#                 coords.range.list = coords.range.list, svg.dim.x = svg.dim.x, svg.dim.y = svg.dim.y))
# }

#########################################################################################################
