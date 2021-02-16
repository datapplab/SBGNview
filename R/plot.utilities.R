
#########################################################################################################
link.to.sbgn.start <- "<a xlink:href=\"https://cdn.rawgit.com/sbgn/process-descriptions/b2904462d11bd8d65e9c7a1318d95d468048cb50/templates/PD_L1V1.3.svg\">"
link.to.sbgn.end <- "</a>"

svg.header <- "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<svg width=\"%f\" height=\"%f\"
xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\"
xmlns:xlink=\"http://www.w3.org/1999/xlink\"
>
<metadata id=\"metadata7\">image/svg+xml</metadata>
"
svg.end <- "
</svg>"

######################################################################################################### 
template.text <- "
<text
x=\"%f\"
y=\"%f\"
id=\"%s\"
style=\"font-size:%fpx;text-anchor:%s;font-family:Arial,Helvetica;alignment-baseline:%s;stroke.opacity:%s;dominant-baseline:%s\"
fill=\"%s\"
>
%s</text>
"

plot.text <- function(x, y, h, w, label, id, label.location, color = "black", glyph.class = "", 
                      stroke.opacity = 1, max.x = 0, if.generate.text.box = FALSE, 
                      parameters.list, glyph) {
    
    if (length(glyph@text$color) > 0) {
        color <- glyph@text$color
    }
    if (length(glyph@text$x.shift) > 0) {
        x <- x + glyph@text$x.shift
    }
    if (length(glyph@text$y.shift) > 0) {
        y <- y + glyph@text$y.shift
    }
    label <- gsub("&", "-", label)
    label <- gsub("<[^><]+>", "", label)  # when using paxtools::toSBGN to convert biopax files to SBGN files, some text in the SBGN files contain special characters :label text='<i>p</i>-chloromercuribenzoate' . the '<' and '>' will create problem in the svg file. So we need remove them.
    if (length(nchar(label)) == 0) {
        # when the label is empty, don't plot it
        print("no label")
        return("")
    }
    result.list <- break.text.into.segments(label, w, glyph.class, parameters.list = parameters.list, 
                                            max.x = max.x, glyph = glyph)
    label <- result.list$words.segments
    font.size <- result.list$font.size
    parameters.list$complex.compartment.label.margin.use <- parameters.list$complex.compartment.label.margin
    if (glyph.class == "compartment") {
        if (max.x > 0) {
            # if the graph is large, make label larger
            parameters.list$complex.compartment.label.margin.use <- parameters.list$complex.compartment.label.margin.use * 
                max(1, max.x/2000)
        }
    }
    
    label.segments <- strsplit(label, "\n")  # the label text in some of the input sbgn files are labeled with '&#xA;' in order to be seperated to different lines
    label.segments <- label.segments[[1]]
    align.shift <- font.size/3
    
    text.anchor <- "middle"
    alignment.baseline <- "middle"
    x.label.box <- x
    if (length(label.segments) > 1) {
        # if the text is labeled to be seperated into different lines
        svg <- ""
        n.lines <- length(label.segments)
        h.label.box <- n.lines * font.size
        for (i in seq_len(length.out = n.lines)) {
            label.segment <- label.segments[i]
            w.label.box <- nchar(label.segment) * font.size
            if (label.location == "center") {
                alignment.baseline <- "baseline"
                first.segment.shift <- (n.lines - 1)/2 * font.size
                start.y <- y - first.segment.shift
                y.segment <- start.y + (i - 1) * font.size
                y.segment <- y.segment + align.shift
                if (i == 1) {
                    y.label.box <- y.segment + align.shift
                }
            } else if (label.location == "top") {
                y.segment <- y - h/2 + 2 + (i - 1) * font.size
                if (glyph.class %in% c("complex", "compartment")) {
                    y.segment <- y - h/2 - n.lines * font.size + (i - 1) * font.size
                }
                text.anchor <- "middle"
                alignment.baseline <- "baseline"
                y.segment <- y.segment + align.shift * 2
                if (i == 1) {
                    y.label.box <- y.segment + align.shift
                }
            }
            svg.segment <- sprintf(template.text, x, y.segment, id, font.size, text.anchor, 
                                   alignment.baseline, stroke.opacity, alignment.baseline, color, label.segment)
            svg <- paste(svg, svg.segment, sep = "\n")
        }
    } else {
        w.label.box <- nchar(label) * font.size * 0.5
        h.label.box <- font.size
        if (glyph.class %in% c("compartment", "complex")) {
            align.shift <- 0
        }
        svg <- template.text
        if (label.location == "center") {
            alignment.baseline <- "baseline"
            svg <- sprintf(template.text, x, y + align.shift, id, font.size, text.anchor, 
                           alignment.baseline, stroke.opacity, alignment.baseline, color, label)
            y.label.box <- y + align.shift
        } else if (label.location == "top") {
            text.anchor <- "middle"
            # alignment.baseline = 'hanging'
            if (glyph.class %in% c("complex", "compartment")) {
                alignment.baseline <- "baseline"
                svg <- sprintf(template.text, x, y + align.shift - h/2 - parameters.list$complex.compartment.label.margin.use, 
                               id, font.size, text.anchor, alignment.baseline, stroke.opacity, 
                               alignment.baseline, color, label)
            } else {
                svg <- sprintf(template.text, x, y + align.shift - h/2 + 4, id, font.size, 
                               text.anchor, alignment.baseline, stroke.opacity, alignment.baseline, 
                               color, label)
            }
            y.label.box <- y - h/2 - parameters.list$complex.compartment.label.margin.use - h.label.box/3
        }
    }
    return(svg)
}

#########################################################################################################
# break long text into multiple segments
break.text.into.segments <- function(label, w, glyph.class, parameters.list, max.x = 0, glyph) {

    if.long.word <- FALSE
    text.length.factor <- parameters.list$text.length.factor.macromolecule
    if (length(glyph@text$font.size) > 0) {
        font.size <- glyph@text$font.size
    } else {
        font.size <- parameters.list$font.size

        font.size0 <- font.size
        if (length(label) == 0) {
            return(list(words.segments = "", label.margin = 2, font.size = font.size,
                        if.long.word = if.long.word, nline = 0))
        } else if (nchar(label) == 0)
        {
            return(list(words.segments = "", label.margin = 2, font.size = font.size,
                        if.long.word = if.long.word, nline = 0))
        }  # when the label is empty, don't plot it
        if (length(w) == 0) {
            w <- 70 * 2/3  # some nodes has no coordinates,
        } else if (is.na(w)) {
            w <- 70 * 2/3  # some nodes has no coordinates,
        }
        if (glyph.class == "compartment") {
            if (parameters.list$if.scale.compartment.font.size) {
                font.size <- max(glyph@shape$stroke.width * 3.5, 
                                 font.size * parameters.list$node.width.adjust.factor.compartment * w)
            } else {
                font.size <- font.size * parameters.list$node.width.adjust.factor
                font.size <- font.size * parameters.list$font.size.scale.gene *
                    parameters.list$font.size.scale.compartment
            }
            if (max.x > 0) {
                # font.size = font.size * max(1,max.x/2000)
            }
            text.length.factor <- parameters.list$text.length.factor.compartment
            
        } else if (glyph.class == "complex") {
            if (parameters.list$if.scale.complex.font.size) {
                font.size <- font.size * parameters.list$node.width.adjust.factor.complex *
                    w
            } else {
                font.size <- font.size * parameters.list$node.width.adjust.factor
                font.size <- font.size * parameters.list$font.size.scale.gene *
                    parameters.list$font.size.scale.complex
            }
            text.length.factor <- parameters.list$text.length.factor.complex
            
        } else if (glyph.class %in% c("omitted process", "uncertain process", "cardinality")) {
            # font.size = font.size * global.parameters.list$logic.node.font.scale
            font.size <- glyph@h
            
        } else if (glyph.class %in% c("tau", "or", "and", "not", "delay", "tag", "terminal")) {
            # font.size = font.size * global.parameters.list$logic.node.font.scale
            font.size <- glyph@h/2
            
        } else if (label == "SBGNhub Pathway Collection") {
            # font.size = font.size * global.parameters.list$logic.node.font.scale
            font.size <- glyph@h * 0.6
            return(list(words.segments = label, label.margin = font.size, font.size = font.size,
                        if.long.word = FALSE, nline = 1))
            
        } else if (glyph.class %in% c("state variable", "unit of information")) {
            font.size <- font.size * parameters.list$status.node.font.scale
        } else {
            font.size <- font.size * 2 * glyph@h * parameters.list$node.width.adjust.factor/70
            font.size <- font.size * parameters.list$font.size.scale.gene
        }
    }
    # sometimes the input label is already splited, then we just need to return it
    if (grepl("\n", label)) {
        label.words <- strsplit(label, "\n")[[1]]
        words.segments <- label
        nline <- length(label.words)
    } else {
        label.words <- strsplit(label, "")[[1]]
        first.word <- strsplit(label, " ")[[1]][1]
        words.segments <- ""
        if (length(label.words) == 1) {
            words.segments <- label
            label.margin <- font.size
            label.length <- (nchar(words.segments)) * font.size
            return(list(words.segments = words.segments, label.margin = label.margin,
                        font.size = font.size, if.long.word = if.long.word, nline = 1))
        }
        current.line <- ""
        nline <- 1
        for (i in seq_len(length.out = length(label.words))) {
            word.current <- label.words[i]
            word.previous <- label.words[i - 1]
            if (length(word.previous) == 0) {
                word.previous <- "currently.first.word"
            }
            current.line.to.be <- paste(current.line, word.current, sep = "")
            current.line.to.be.length <- (nchar(current.line.to.be)) * font.size

            if (current.line.to.be.length > text.length.factor * w & 
                (length(label.words) - i) > 2 & 
                (word.previous %in% parameters.list$label.spliting.string | identical(parameters.list$label.spliting.string, c("any"))) ) {
                
                # if there are less than 2 letters left ,merge them to the current line instead
                # of creating a new line with just 2 letters
                if (glyph.class == "complex") {
                    if (first.word == "Complex") {
                        words.segments <- ""
                        (break)()
                    }
                }
                if (words.segments != "") {
                    nline <- nline + 1
                }
                if (word.current %in% c(" ") | word.previous == "-" | length(grep("[a-zA-Z]",
                                                                                  word.previous)) == 0) {
                    # if this is the beginning or end of the word, don't insert '-'
                    words.segments <- paste(words.segments, word.current, sep = "\n")
                } else {
                    words.segments <- paste(paste(words.segments, "", sep = ""), word.current,
                                            sep = "\n")
                }
                current.line <- word.current
            } else {
                words.segments <- paste(words.segments, word.current, sep = "")
                current.line <- current.line.to.be
            }
        }
        words.segments <- gsub("^\n", "", words.segments)
    }
    label.margin <- (nline) * (font.size)
    if (nline > 3) {
        if.long.word <- TRUE
        label.margin <- (font.size)
    }
    return(list(words.segments = words.segments, label.margin = label.margin, 
                font.size = font.size, if.long.word = if.long.word, nline = nline))
}

######################################################################################################### 
template.path <- "
<path
d=\"%s\"
id=\"%s\"
style=\"fill:%s;fill-opacity:%f;stroke-width:%f;stroke:%s;stroke-opacity:%f\"
transform = \"rotate(%f %f %f)\"
/>
"

plot.path <- function(d, id, fill.color, stroke.width, stroke.color, angle = 0, rotate.x = 0, 
                      rotate.y = 0, fill.opacity = 0.6, stroke.opacity = 1, object) {
    
    if ("shape" %in% slotNames(object)) {
        # if we are ploting a glyph(glyph has slot 'shape')
        if (length(object@shape$fill.color) > 0) {
            fill.color <- object@shape$fill.color
        }
        if (length(object@shape$fill.opacity) > 0) {
            fill.opacity <- object@shape$fill.opacity
        }
        if (length(object@shape$stroke.width) > 0) {
            stroke.width <- object@shape$stroke.width
        }
        if (length(object@shape$stroke.opacity) > 0) {
            stroke.opacity <- object@shape$stroke.opacity
        }
        if (length(object@shape$stroke.color) > 0) {
            stroke.color <- object@shape$stroke.color
        }
        # }else{ if we are ploting an arc
    } else if ("edge" %in% slotNames(object)) {
        fill.color <- object@edge$tip.fill.color
        fill.opacity <- object@edge$tip.fill.opacity
        stroke.width <- object@edge$tip.stroke.width
        stroke.opacity <- object@edge$tip.stroke.opacity
        stroke.color <- object@edge$tip.stroke.color
    } else {
        message("Plot.path function is not implemented for this object!!! \n Please make sure the object has a slot named 'shape' or 'edge'.")
        print(object)
        stop()
    }
    svg <- sprintf(template.path, d, id, fill.color, fill.opacity, stroke.width, 
        stroke.color, stroke.opacity, angle, rotate.x, rotate.y)
    return(svg)
}

#########################################################################################################
template.rectangle <- "
<rect
width=\"%f\"
height=\"%f\"
rx=\"%f\"
ry=\"%f\"
x=\"%f\"
y=\"%f\"
id=\"%s\"
style=\"fill:%s;fill-opacity:%f;stroke:%s;stroke-width:%f;stroke-opacity:%f\" />
"

plot.rectangle <- function(x, y, h, w, id, rx = 0, ry = 0, fill.color = "white", 
                           stroke.color = "black", stroke.width = 1, fill.opacity = 0.6, 
                           stroke.opacity = 1, glyph = new("glyph")) {
    
    if (length(glyph@shape$fill.color) > 0) {
        fill.color <- glyph@shape$fill.color
    }
    if (length(glyph@shape$fill.opacity) > 0) {
        fill.opacity <- glyph@shape$fill.opacity
    }
    if (length(glyph@shape$stroke.width) > 0) {
        stroke.width <- glyph@shape$stroke.width
    }
    if (length(glyph@shape$stroke.opacity) > 0) {
        stroke.opacity <- glyph@shape$stroke.opacity
    }
    if (length(glyph@shape$stroke.color) > 0) {
        stroke.color <- glyph@shape$stroke.color
    }
    svg <- sprintf(template.rectangle, w, h, rx, ry, x, y, id, fill.color, fill.opacity, 
                   stroke.color, stroke.width, stroke.opacity)
    return(svg)
}

#########################################################################################################
# use path to plot a rectangle that has 1,2,3 or 4 rounded corner. Returns the
# content for 'd' command in 'path' the r of the rounded corners can be
# different(can be set using r1,r2,r3,r4: from top left corner, clockwise)
generate.d.for.partly.rounded.corner.rect <- function(x, y, w, h, r1, r2, r3, r4) {
    
    start.and.upperLeft <- paste("M", x, r1 + y, "Q", x, y, x + r1, y, sep = " ")
    top.and.upperRight <- paste("L", x + w - r2, y, "Q", x + w, y, x + w, y + r2, sep = " ")
    right.and.lowerRight <- paste("L", x + w, y + h - r3, "Q", x + w, y + h, 
                                  x + w - r3, y + h, sep = " ")
    bottom.and.lowerLeft <- paste("L", x + r4, y + h, "Q", x, y + h, x, y + h - r4, 
                                  "Z", sep = " ")
    d.content <- paste(start.and.upperLeft, top.and.upperRight, right.and.lowerRight, 
                       bottom.and.lowerLeft, sep = " ")
    return(d.content)
}

#########################################################################################################
template.polygon <- "
<polygon
points = \"%s\"
id=\"line\"
style=\"fill:%s; stroke:%s; stroke-opacity:%s; stroke-width:%f\"
transform = \"rotate(%f %f %f)\"
/>
"
plot.polygon <- function(object, points, fill.color, rotate.a, rotate.x, rotate.y) {
    
    if (fill.color == "black") {
        fill.color <- object@edge$tip.fill.color  # if the edge has a filled tip(e.g. production), we can change it to user defined color. If the edge type has an empty tip, we can't change it.
    }
    sprintf(template.polygon, points, fill.color, object@edge$tip.stroke.color, object@edge$tip.stroke.opacity, 
            object@edge$tip.stroke.width, rotate.a, rotate.x, rotate.y)
}

#########################################################################################################
template.ellipse <- "
<ellipse
rx=\"%f\"
ry=\"%f\"
cx=\"%f\"
cy=\"%f\"
id=\"%s\"
style=\"fill:%s;fill-opacity:%f;stroke:%s;stroke-width:%f;stroke-opacity:%f\" />
"

plot.ellipse <- function(x, y, h, w, id, rx = 0, ry = 0, fill.color = "white", 
                         fill.opacity = 1, stroke.color = "black", stroke.width = 1, 
                         stroke.opacity = 1, glyph = new("glyph")) {
    
    if (length(glyph@shape$fill.color) > 0) {
        fill.color <- glyph@shape$fill.color
    }
    if (length(glyph@shape$fill.opacity) > 0) {
        fill.opacity <- glyph@shape$fill.opacity
    }
    if (length(glyph@shape$stroke.width) > 0) {
        stroke.width <- glyph@shape$stroke.width
    }
    if (length(glyph@shape$stroke.opacity) > 0) {
        stroke.opacity <- glyph@shape$stroke.opacity
    }
    if (length(glyph@shape$stroke.color) > 0) {
        stroke.color <- glyph@shape$stroke.color
    }
    svg <- sprintf(template.ellipse, w/2, h/2, x, y, id, fill.color, fill.opacity, 
                   stroke.color, stroke.width, stroke.opacity)
    return(svg)
}

#########################################################################################################
plot.ellipse.edge.tip <- function(x, y, h, w, id, rx = 0, ry = 0, fill.color = "white", 
                                  fill.opacity = 1, stroke.color = "black", stroke.width = 1, 
                                  stroke.opacity = 1, object) {
    
    if (fill.color == "black") {
        fill.color <- object@edge$tip.fill.color  # if the edge has a filled tip(e.g. production), we can change it to user defined color. If the edge type has empty tip fill color, we can't change it.
    }
    svg <- sprintf(template.ellipse, w/2, h/2, x, y, id, fill.color, 
                   object@edge$tip.fill.opacity, object@edge$tip.stroke.color, 
                   object@edge$tip.stroke.width, object@edge$tip.stroke.opacity)
    return(svg)
}

#########################################################################################################
template.line <- "
<line
%s
id=\"%s\"
style=\"stroke:%s;stroke-width:%f;stroke-opacity:%f\" />
"

plot.line <- function(LineCoordinates, id, stroke.color = "black", stroke.opacity = 1, 
                      stroke.width = 1, object = NULL) {
    
    if (!is.null(object)) {
        stroke.color <- object@edge$line.stroke.color
        stroke.width <- object@edge$line.width
    }
    stroke.width <- 2 * stroke.width
    svg <- template.line
    svg <- sprintf(template.line, LineCoordinates, id, stroke.color, stroke.width, stroke.opacity)
    return(svg)
}

#########################################################################################################
plot.LineWings <- function(x, y, w, h, orientation, id) {
    
    leftLineCoordinates <- paste("x1=", x - w, " y1=", y, " x2=", x - w/2, " y2=", 
                                 y, "", sep = "\"")
    rightLineCoordinates <- paste("x1=", x + w/2, " y1=", y, " x2=", x + w, " y2=",
                                  y, "", sep = "\"")
    if (orientation == "vertical") {
        leftLineCoordinates <- paste("x1=", x, " y1=", y - h, " x2=", x, " y2=",
                                     y - h/2, "", sep = "\"")
        rightLineCoordinates <- paste("x1=", x, " y1=", y + h/2, " x2=", x, " y2=",
                                      y + h, "", sep = "\"")
    }
    svg.line.l <- plot.line(leftLineCoordinates, id)
    svg.line.r <- plot.line(rightLineCoordinates, id)
    svg <- paste(svg.line.l, svg.line.r, sep = "\n")
    return(svg)
}

#########################################################################################################
#' An object to store glyph information
#' 
#' @slot compartment 
#' A character string. The compartment this glyph belongs to.
#' 
#' @slot x,y,h,w
#' Numeric. The x,y location and hight, width of the glyph.
#' 
#' @slot id
#' A character string. The ID of the glyph (determined by the 'id' attribute in element 'glyph').
#' 
#' @slot label
#' A character string. label of the glyph. Extracted from the element 'label'.
#' 
#' @slot glyph.class 
#' A character string. Class of the glyph.
#' 
#' @slot user.data 
#' A numeric vector. The omics data mapped to this glyph.
#' 
#' @slot stroke.width,stroke.opacity,fill.opacity 
#' Numeric. Controls glyph border width or opacity of node fill and border.
#' 
#' @slot fill.color,stroke.color 
#' A character string. The fill color and stroke color of glyph.
#' 
#' @slot label_location 
#' A character string. One of 'center' or 'top'. If set to 'center', the glyph label text will be displayed at the center of the glyph. If set to 'top', the label will be displayed at the top of the glyph. 
#' 
#' @slot orientation 
#' A character string. One of 'left','right','up','down'. This only applies to glyphs of classes 'terminal' and 'tag'. It will change the orientation of the node. 
#' 
#' @slot parameters.list 
#' A list. The parameters.list slot is a copy of the global.parameters.list. The global.parameters.list contains the '...' parameters passed into \code{\link{SBGNview}} which will be used by \code{\link{renderSbgn}}. By default this slot for all glyph objects in the output SBGNview object will be an empty list. You can add named elements to this slot to customize each individual glyph. For information about which parameters can be set in the parameters.list slot, see the arguments in \code{\link{renderSbgn}}. A customized list of glyph objects can be passed as input into SBGNview function using the ‘glyphs.user’ argument (see \code{\link{renderSbgn}} for more information).  
#'              
#' @slot shape 
#' A list. Parameters controlling how the node is rendered. Available elements can be accessed as follow:
#'              glyph@@shape$stroke.color\cr
#'                      -Character (e.g. 'red'). The color of glyph border. \cr
#'              glyph@@shape$stroke.opacity\cr
#'                      -Numeric (e.g. 0.8). The opacity of glyph border. Must be between 0 and 1.\cr
#'              glyph@@shape$stroke.width\cr
#'                      -Numeric. The thickness of glyph border.\cr
#'              glyph@@shape$fill.opacity\cr
#'                      -Numeric. The opacity of glyph fill. Must be between 0 and 1.\cr
#'              glyph@@shape$fill.color\cr
#'                      -Numeric. The glyph fill color. Note that this will overwrite colors generated from omics data!\cr
#'              
#' @slot text
#' A list. Parameters controlling how node label is rendered. Available elements can be accessed as follow:
#'              glyph@@text$y.shift\cr
#'                      -Numeric. Move glyph label vertically by this value.\cr
#'              glyph@@text$x.shift\cr
#'                      -Numeric. Move glyph label horizontally by this value.\cr
#'              glyph@@text$color\cr
#'                      -A character string. The color of glyph label.\cr
#'              glyph@@text$font.size\cr
#'                      -A character string. The font size of glyph label.\cr
#'                      
#' @details User can modify glyph objects to change the way how they are rendered. See the "glyphs.user" argument in \code{\link{renderSbgn}}. 

setClass("glyph", slots = c(compartment = "character", x = "numeric", y = "numeric", 
            h = "numeric", w = "numeric", id = "character", label = "character", label.margin = "numeric", 
            label_location = "character", svg.port = "character", orientation = "character", 
            user.data = "vector", clone = "list", fill.color = "character", stroke.width = "numeric", 
            stroke.color = "character", stroke.opacity = "numeric", fill.opacity = "numeric", 
            parameters.list = "list", glyph.class = "character", text = "list", shape = "list"))

setClass("port", contains = "glyph")

#########################################################################################################
#' An object to store arc information
#' 
#' @slot target,source,id,arc.class 
#' A character string. Information extracted from element 'arc'.
#' 
#' @slot start.x,start.y,end.x,end.y 
#' Numeric. Information extracted from elements 'start', 'end' or 'next'.
#' 
#' @slot stroke.opacity 
#' Numeric. Controls the line of an arc (not tip of the arc).
#' 
#' @slot parameters.list 
#' A list. The parameters.list slot is a copy of the global.parameters.list. The global.parameters.list contains the '...' parameters passed into \code{\link{SBGNview}} which will be used by \code{\link{renderSbgn}}. By default this slot for all arc objects in the output SBGNview object will be an empty list. You can add named elements to this slot to customize each individual arc. For information about which parameters can be set in the parameters.list slot, see the arguments in \code{\link{renderSbgn}}. A customized list of arc objects can be passed as input into SBGNview function using the ‘arcs.user’ argument (see \code{\link{renderSbgn}} for more information). 
#'              
#' @slot edge 
#' A list. An arc in SBGN map normally consists of a line and a tip at the end of the line. This list holds variables that controls arc properties. Available elements can be accessed as follow:
#'           The following three parameters control the line:\cr
#'              arc@@edge$line.stroke.color\cr
#'              arc@@edge$line.stroke.opacity\cr
#'              arc@@edge$line.width\cr
#'                      
#'           The following five parameters control the tip:\cr
#'              arc@@edge$tip.stroke.opacity\cr
#'              arc@@edge$tip.stroke.color\cr
#'              arc@@edge$tip.stroke.width\cr
#'              arc@@edge$tip.fill.color\cr
#'              arc@@edge$tip.fill.opacity\cr
#'              
#' @details Information in an 'arc' object comes from two sources: 1. SBGN-ML file's 'arc' element ('source', 'target', coordinates etc.). 2. Parameters specified when running function \code{\link{SBGNview}}.  User can modify arc objects to change the way how they are rendered. See the "arcs.user" argument in \code{\link{renderSbgn}}.

setClass("arc", slots = c(target = "character", source = "character", id = "character", 
                          arc.class = "character", start.y = "numeric", start.x = "numeric", end.y = "numeric", 
                          end.x = "numeric", stroke.opacity = "numeric", parameters.list = "list", 
                          edge = "list"))

setClass("spline", slots = c(id = "character", spline.coords = "vector", edge = "list", 
                                  parameters.list = "list"))

#########################################################################################################
#' An object to store information of spline arcs
#' 
#' @slot target,source,id,arc.class 
#' A character string. Information extracted from element 'arc.spline'.
#' 
#' @slot start.x,start.y,end.x,end.y 
#' Numeric. Information extracted from elements 'start', 'end' or 'next'.
#' 
#' @slot stroke.opacity 
#' Numeric. Controls the line of an arc (not tip of arc).
#' 
#' @slot components 
#' A list of 'arc' and 'spline' objects. A spline arc is represented by several components: 1. The two ends of the arc are represented by straight line 'arc' objects. 2. splines connecting the ends are represented by 'spline' objects.
#' 
#' @slot parameters.list 
#' A list. The parameters.list slot is a copy of the global.parameters.list. The global.parameters.list contains the '...' parameters passed into \code{\link{SBGNview}} which will be used by \code{\link{renderSbgn}}. By default this slot for all spline.arc objects in the output SBGNview object will be an empty list. You can add named elements to this slot to customize each individual spline.arc. For information about which parameters can be set in the parameters.list slot, see the arguments in \code{\link{renderSbgn}}. A customized list of arc objects can be passed as input into SBGNview function using the ‘arcs.user’ argument (see \code{\link{renderSbgn}} for more information).
#'              
#' @slot edge 
#' A list. An arc in SBGN map normally consists of a line and a tip shape at the end of the line. This list holds variables that controls arc properties. Available elements can be accessed as follow:
#'  
#'           The following three parameters control the line:\cr
#'              arc@@edge$line.stroke.color\cr
#'              arc@@edge$line.stroke.opacity\cr
#'              arc@@edge$line.width\cr
#'                      
#'           The following five parameters control the tip:\cr
#'              arc@@edge$tip.stroke.opacity\cr
#'              arc@@edge$tip.stroke.color\cr
#'              arc@@edge$tip.stroke.width\cr
#'              arc@@edge$tip.fill.color\cr
#'              arc@@edge$tip.fill.opacity\cr
#'              
#' @details Arc information comes from two sources:1. SBGN-ML file's 'arc' element ('source', 'target', coordinates etc.). 2. Parameters specified when running \code{\link{SBGNview}}.  User can modify arc objects to change the way how they are rendered. See the "arcs.user" argument in \code{\link{renderSbgn}}.

setClass("spline.arc", slots = c(id = "character", source = "character", target = "character", 
        arc.class = "character", edge = "list", components = "list", parameters.list = "list"))

#########################################################################################################
# generic function to plot glyph and arc defined 
# by classes inheriting from "glyph" and "arc" classes

setGeneric("plot.glyph", function(object) {
    standardGeneric("plot.glyph")
})

setGeneric("plot.arc", function(object) {
    standardGeneric("plot.arc")
})

#########################################################################################################
###### Classes implementing the glyph class and defining plot.glyph generic function
#########################################################################################################
setClass("macromolecule.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("macromolecule.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    label <- object@label
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (length(object@shape$stroke.width) > 0) {
        stroke.width <- object@shape$stroke.width
    }
    stroke.color <- 1
    if (length(object@stroke.color) > 0) {
        stroke.color <- object@stroke.color
    }
    stroke.opacity <- 1
    if (length(object@stroke.opacity) > 0) {
        stroke.opacity <- object@stroke.opacity
    }
    fill.opacity <- 1
    if (length(object@fill.opacity) > 0) {
        fill.opacity <- object@fill.opacity
    }
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    # svg = paste(svg.rect,sep='\n')
    svg <- ""
    
    if (is.na(user.data[1])) {
        
    }
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        for (i in seq_len(length.out = length(user.data))) {
            fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                mol.type = "gene")
            glyph.data.i <- object
            glyph.data.i@fill.color <- fill.color
            glyph.data.i@fill.opacity <- 1
            glyph.data.i@user.data <- c("no.user.data")
            glyph.data.i@w <- w/length(user.data)
            glyph.data.i@x <- x - w/2 + (i - 1/2) * glyph.data.i@w
            glyph.data.i@label <- " "
            glyph.data.i@stroke.opacity <- 0.5
            glyph.data.i@stroke.color <- "black"
            glyph.data.i@shape$stroke.opacity <- 0.5  # if we are plotting data shape, don't plot the border
            glyph.data.i@shape$stroke.width <- 0.0 # changed 0.1 to 0.0. if we are plotting data shape, don't plot the border
            svg <- paste(svg, plot.glyph(glyph.data.i), sep = "\n")
        }
        fill.opacity <- 0
        stroke.opacity <- 1
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, rx = rx, ry = ry, stroke.width = stroke.width, fill.color = fill.color, 
        stroke.opacity = stroke.opacity, fill.opacity = fill.opacity)
    svg <- paste(svg, svg.rect, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("entity.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("entity.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    label <- object@label
    label.location <- "center"
    stroke.width <- 1
    stroke.opacity <- 1
    if (length(object@shape$stroke.width) > 0) {
        stroke.width <- object@shape$stroke.width
    }
    if (length(object@stroke.opacity) > 0) {
        stroke.opacity <- object@stroke.opacity
    }
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    # svg.rect =
    # plot.rectangle(glyph=object,x=x-w/2,y=y-h/2,h=h,w=w,id=id,rx=rx,ry=ry,stroke.width
    # = stroke.width,stroke.opacity = stroke.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    # svg = paste(svg.rect,sep='\n')
    svg <- ""
    
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        for (i in seq_len(length.out = length(user.data))) {
            fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                mol.type = "gene")
            glyph.data.i <- object
            glyph.data.i@fill.color <- fill.color
            glyph.data.i@user.data <- c("no.user.data")
            glyph.data.i@w <- w/length(user.data)
            glyph.data.i@x <- x - w/2 + (i - 1/2) * glyph.data.i@w
            glyph.data.i@label <- " "
            glyph.data.i@shape$stroke.width <- 0.5
            glyph.data.i@stroke.opacity <- 0
            svg <- paste(svg, plot.glyph(glyph.data.i), sep = "\n")
        }
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, rx = rx, ry = ry, stroke.width = stroke.width, stroke.opacity = stroke.opacity)
    svg <- paste(svg, svg.rect, svg.text, sep = "\n")
    
    return(svg)
})

#########################################################################################################
setClass("compartment.sbgn", contains = "glyph", slots = c(max.x = "numeric"))
setMethod("plot.glyph", signature("compartment.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "top"
    # ry = rx
    max.x <- object@max.x
    rx <- ry <- 20
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, id = id, 
        h = h, w = w, rx = rx, ry = ry, stroke.width = object@shape$stroke.width, 
        fill.opacity = parameters.list$compartment.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location, 
        max.x = max.x)
    svg <- paste(svg.rect, svg.text, sep = "\n")
    # svg = ''
    return(svg)
})

#########################################################################################################
setClass("complex.sbgn", contains = "glyph", slots = c(if.complex.empty = "logical"))
setMethod("plot.glyph", signature("complex.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    if (length(nchar(object@label_location)) == 0) {
        label.location <- "top"
    } else if (object@label_location != "center") {
        label.location <- "top"
    } else {
        label.location <- object@label_location
    }
    if.complex.empty <- FALSE
    if (length(object@if.complex.empty) != 0) {
        if.complex.empty <- object@if.complex.empty
    }
    if (if.complex.empty) {
        label.location <- "center"
    }
    rx <- ry <- 8
    corner.r <- rx
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    d1 <- (paste("M", x, y + corner.r, x, y + h - corner.r, x + corner.r, y + h, 
        x + w - corner.r, y + h, x + w, y + h - corner.r, x + w, y + corner.r, x + 
            w - corner.r, y, x + corner.r, y, "Z", sep = " "))
    svg.path1 <- plot.path(object = object, d = d1, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black")
    
    # label='complex' # try not to show complex name
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path1, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("complex_multimer.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("complex_multimer.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    corner.r <- parameters.list$edge.tip.size * 3
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    d1 <- (paste("M", x, y + corner.r, x, y + h - corner.r, x + corner.r, y + h, 
        x + w - corner.r, y + h, x + w, y + h - corner.r, x + w, y + corner.r, x + 
            w - corner.r, y, x + corner.r, y, "Z", sep = " "))
    svg.path1 <- plot.path(object = object, d = d1, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black", fill.opacity = 1)
    
    x <- x + parameters.list$edge.tip.size * 3/2
    y <- y + parameters.list$edge.tip.size * 3/2
    d2 <- (paste("M", x, y + corner.r, x, y + h - corner.r, x + corner.r, y + h, 
        x + w - corner.r, y + h, x + w, y + h - corner.r, x + w, y + corner.r, x + 
            w - corner.r, y, x + corner.r, y, "Z", sep = " "))
    svg.path2 <- plot.path(object = object, d = d2, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black")
    
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path2, svg.path1, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("state_variable.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("state_variable.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label <- gsub(" ", "", label)
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
        fill.opacity = parameters.list$auxiliary.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.ellipse, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("unspecified_entity.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("unspecified_entity.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    fill.opacity <- 1
    # if(length(object@fill.opacity)>0){ fill.opacity = object@fill.opacity }
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (length(object@shape$stroke.width) > 0) {
        stroke.width <- object@shape$stroke.width
    }
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    # svg.ellipse = plot.ellipse(glyph=object,x,y,h,w,id,stroke.width = stroke.width)
    # svg = paste(svg.ellipse,sep='\n')
    svg <- ""
    
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        fill.color <- color.from.value(user.data[1], global.parameters.list = parameters.list, 
            mol.type = "gene")
        svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
            fill.color = fill.color, stroke.opacity = 0)
        svg <- paste(svg, svg.ellipse, sep = "\n")
        if (length(user.data) > 1) {
            for (i in 2:length(user.data)) {
                fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                  mol.type = "gene")
                glyph.data.i <- object
                glyph.data.i@fill.color <- fill.color
                glyph.data.i@user.data <- c("no.user.data")
                glyph.data.i@x <- x - w/2 + (i - 1) * w/length(user.data)
                glyph.data.i@label <- " "
                glyph.data.i@shape$stroke.width <- 0.5
                start.x <- x - w/2 + (i - 1) * w/length(user.data)
                start.y.1 <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                start.y.2 <- y - (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                if (start.x > x) {
                  large.flag <- 0
                } else {
                  large.flag <- 1
                }
                d.1 <- (paste("M", start.x, start.y.1, "A", w/2, h/2, 1, large.flag, 
                  0, start.x, start.y.2, sep = " "))
                svg.path.1 <- plot.path(object = object, d = d.1, id = id, fill.color = fill.color, 
                  stroke.width = 1, stroke.color = "black", stroke.opacity = 0)
                
                start.x <- x - w/2 + (i) * w/length(user.data)
                start.y.1 <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                start.y.2 <- y - (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                if (start.x > x) {
                  large.flag <- 0
                } else {
                  large.flag <- 1
                }
                d.2 <- (paste("M", start.x, start.y.1, "A", w/2, h/2, 1, large.flag, 
                  0, start.x, start.y.2, sep = " "))
                svg.path.2 <- plot.path(object = object, d = d.2, id = id, fill.color = "white", 
                  stroke.width = 1, stroke.color = "black", stroke.opacity = 0)
                svg <- paste(svg, svg.path.1, svg.path.2, sep = "\n")
            }
        }
        fill.opacity <- 0
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
        fill.opacity = fill.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg, svg.ellipse, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("unit_of_information.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("unit_of_information.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label <- gsub(" ", "", label)
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width, fill.opacity = parameters.list$auxiliary.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("stoichiometry.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("stoichiometry.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("cardinality.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("cardinality.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("biological_activity.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("biological_activity.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("submap.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("submap.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "top"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg <- svg.rect
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    # svg = paste(svg.rect,svg.wings,sep='\n')
    svg <- paste(svg.rect, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("process.sbgn", contains = "glyph", slots = c(if.show.label = "logical"))
setMethod("plot.glyph", signature("process.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (!object@if.show.label) {
        label <- ""
    } else {
        label <- object@label
    }
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg <- svg.rect
    # svg.wings = plot.LineWings(x,y,w,h,orientation,id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("uncertain_process.sbgn", contains = "glyph", slots = c(if.show.label = "logical"))
setMethod("plot.glyph", signature("uncertain_process.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (!object@if.show.label) {
        label <- "?"
    } else {
        label <- object@label
    }
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("and.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("and.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- "AND"
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.ellipse, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("or.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("or.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- "OR"
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.ellipse, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("not.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("not.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- "NOT"
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.ellipse, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("tau.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("tau.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- "t"
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.ellipse, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("macromolecule_multimer.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("macromolecule_multimer.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    fill.opacity <- 1
    if (length(object@fill.opacity) > 0) {
        fill.opacity <- object@fill.opacity
    }
    label <- object@label
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    stroke.opacity <- 1
    if (length(object@stroke.opacity) > 0) {
        stroke.opacity <- object@stroke.opacity
    }
    label.location <- "center"
    rx <- min(h, w)/10
    ry <- min(h, w)/10
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    # svg.rect.up =
    # plot.rectangle(glyph=object,x=x-w/2,y=y-h/2,id=id,h=h,w=w,rx=rx,ry=ry,stroke.width
    # = stroke.width,fill.opacity = 1)
    svg.rect.down <- plot.rectangle(glyph = object, x = x - w/2 + parameters.list$multimer.margin, 
        y = y - h/2 + parameters.list$multimer.margin, id = id, h = h, w = w, 
        rx = rx, ry = ry, stroke.width = stroke.width)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    # svg = paste(svg.rect.down,svg.rect.up,sep='\n')
    svg <- paste(svg.rect.down, sep = "\n")
    # svg=''
    
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        for (i in seq_len(length.out = length(user.data))) {
            fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                mol.type = "gene")
            glyph.data.i <- new("macromolecule.sbgn")
            glyph.data.i@fill.color <- fill.color
            glyph.data.i@user.data <- c("no.user.data")
            glyph.data.i@compartment <- object@compartment
            glyph.data.i@y <- y
            glyph.data.i@h <- h
            glyph.data.i@w <- w/length(user.data)
            glyph.data.i@x <- x - w/2 + (i - 1/2) * glyph.data.i@w
            glyph.data.i@label <- " "
            glyph.data.i@id <- paste(id, "user.data", i, sep = "_")
            glyph.data.i@shape$stroke.width <- 0.5
            glyph.data.i@stroke.opacity <- 0
            glyph.data.i@fill.opacity <- 1
            svg <- paste(svg, plot.glyph(glyph.data.i), sep = "\n")
        }
        fill.opacity <- 0
    }
    svg.rect.up <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, id = id, 
        h = h, w = w, rx = rx, ry = ry, stroke.width = stroke.width, fill.opacity = fill.opacity)
    svg <- paste(svg, svg.rect.up, svg.text, sep = "\n")
    
    return(svg)
})

#########################################################################################################
setClass("outcome.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("outcome.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
        print(orientation)
    } else {
        orientation <- "horizontal"
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, fill.color = "black", 
        stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg <- svg.ellipse
    return(svg)
})

#########################################################################################################
setClass("simple_chemical_multimer.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("simple_chemical_multimer.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    w <- h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    fill.opacity <- 1
    # if(length(object@fill.opacity)>0){ fill.opacity = object@fill.opacity }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    # svg.ellipse = plot.ellipse(glyph=object,x,y,h,w,id,stroke.width =
    # stroke.width,fill.opacity = 1)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg.ellipse2 <- plot.ellipse(glyph = object, x + 5, y + 5, h, w, id)
    svg <- paste(svg.ellipse2, sep = "\n")
    
    if (!identical(user.data[1], "no.user.data")) {
        # if(user.data[1] != 'no.user.data'){
        fill.color <- color.from.value(user.data[1], global.parameters.list = parameters.list, 
            mol.type = "cpd")
        svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
            fill.color = fill.color, fill.opacity = 1)
        svg <- paste(svg, svg.ellipse, sep = "\n")
        if (length(user.data) > 1) {
            for (i in 2:length(user.data)) {
                fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                  mol.type = "cpd")
                glyph.data.i <- object
                glyph.data.i@fill.color <- fill.color
                glyph.data.i@user.data <- c("no.user.data")
                glyph.data.i@fill.opacity <- 1
                # glyph.data.i@w = w/length(user.data)
                glyph.data.i@x <- x - w/2 + (i - 1) * w/length(user.data)
                glyph.data.i@label <- " "
                glyph.data.i@shape$stroke.width <- 0.5
                start.x <- x - w/2 + (i - 1) * w/length(user.data)
                start.y.1 <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                start.y.2 <- y - (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                if (start.x > x) {
                  large.flag <- 0
                } else {
                  large.flag <- 1
                }
                d.1 <- (paste("M", start.x, start.y.1, "A", w/2, h/2, 1, large.flag, 
                  0, start.x, start.y.2, sep = " "))
                svg.path.1 <- plot.path(object = object, d = d.1, id = id, fill.color = fill.color, 
                  stroke.width = 1, stroke.color = "black", fill.opacity = 1, stroke.opacity = 0)
                
                start.x <- x - w/2 + (i) * w/length(user.data)
                start.y.1 <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                start.y.2 <- y - (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                if (start.x > x) {
                  large.flag <- 0
                } else {
                  large.flag <- 1
                }
                d.2 <- (paste("M", start.x, start.y.1, "A", w/2, h/2, 1, large.flag, 
                  0, start.x, start.y.2, sep = " "))
                svg.path.2 <- plot.path(object = object, d = d.2, id = id, fill.color = "white", 
                  stroke.width = 1, stroke.color = "black", stroke.opacity = 0)
                svg <- paste(svg, svg.path.1, svg.path.2, sep = "\n")
            }
        }
        fill.opacity <- 0
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
        fill.opacity = fill.opacity)
    svg <- paste(svg, svg.ellipse, svg.text, sep = "\n")
    
    return(svg)
})

#########################################################################################################
setClass("omitted_process.sbgn", contains = "glyph", slots = c(if.show.label = "logical"))
setMethod("plot.glyph", signature("omitted_process.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (!object@if.show.label) {
        label <- "\\\\"
    } else {
        label <- object@label
    }
    
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.rect <- plot.rectangle(glyph = object, x = x - w/2, y = y - h/2, h = h, w = w, 
        id = id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.rect, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("simple_chemical.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("simple_chemical.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    w <- h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    label <- object@label
    label.location <- "center"
    stroke.color <- 1
    if (length(object@stroke.color) > 0) {
        stroke.color <- object@stroke.color
    }
    stroke.width <- 1
    if (length(object@shape$stroke.width) > 0) {
        stroke.width <- object@shape$stroke.width
    }
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    # svg.ellipse = plot.ellipse(glyph=object,x,y,h,w,id,stroke.width =
    # stroke.width,fill.color = fill.color)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    # svg = paste(svg.ellipse,sep='\n')
    svg <- ""
    
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        fill.color <- color.from.value(user.data[1], global.parameters.list = parameters.list, 
            mol.type = "cpd")
        svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
            fill.color = fill.color, stroke.opacity = 0.5)
        svg <- paste(svg, svg.ellipse, sep = "\n")
        if (length(user.data) > 1) {
            svg.sep.data.line <- ""
            for (i in 2:length(user.data)) {
                fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                  mol.type = "cpd")
                glyph.data.i <- object
                glyph.data.i@fill.color <- fill.color
                glyph.data.i@user.data <- c("no.user.data")
                # glyph.data.i@w = w/length(user.data)
                glyph.data.i@x <- x - w/2 + (i - 1) * w/length(user.data)
                glyph.data.i@label <- " "
                glyph.data.i@shape$stroke.width <- 0.5
                start.x <- x - w/2 + (i - 1) * w/length(user.data)
                start.y.1 <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                start.y.2 <- y - (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                
                Coordinates <- paste("x1=", start.x, " y1=", start.y.1, " x2=", start.x, 
                  " y2=", start.y.2, "", sep = "\"")
                svg.sep.data.line <- paste(svg.sep.data.line, plot.line(LineCoordinates = Coordinates, 
                  id = id, stroke.width = 0.0, stroke.opacity = 0.5), sep = "\n")
                
                if (start.x > x) {
                  large.flag <- 0
                } else {
                  large.flag <- 1
                }
                d.1 <- (paste("M", start.x, start.y.1, "A", w/2, h/2, 1, large.flag, 
                  0, start.x, start.y.2, sep = " "))
                object.data <- object  # if we are ploting shape for data, don't plot border
                object.data@shape$stroke.opacity <- 0
                svg.path.1 <- plot.path(object = object.data, d = d.1, id = id, fill.color = fill.color, 
                  stroke.width = 0.1, stroke.color = "black", stroke.opacity = 0)
                
                start.x <- x - w/2 + (i) * w/length(user.data)
                start.y.1 <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                start.y.2 <- y - (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
                if (start.x > x) {
                  large.flag <- 0
                } else {
                  large.flag <- 1
                }
                d.2 <- (paste("M", start.x, start.y.1, "A", w/2, h/2, 1, large.flag, 
                  0, start.x, start.y.2, sep = " "))
                svg.path.2 <- plot.path(object = object.data, d = d.2, id = id, fill.color = "white", 
                  stroke.width = 0.1, stroke.color = "black", stroke.opacity = 0)
                svg <- paste(svg, svg.path.1, svg.path.2, sep = "\n")
            }
            svg <- paste(svg, svg.sep.data.line, sep = "\n")
        }
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width, 
        fill.opacity = 0)
    svg <- paste(svg, svg.ellipse, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("dissociation.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("dissociation.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse.o <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    svg.ellipse.i <- plot.ellipse(glyph = object, x, y, h/2, w/2, id)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg <- paste(svg.ellipse.o, svg.ellipse.i, svg.wings, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("nucleic_acid_feature_multimer.sbgn", contains = "glyph", slots = c(r4 = "numeric", r3 = "numeric"))
setMethod("plot.glyph", signature("nucleic_acid_feature_multimer.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    fill.opacity <- 1
    if (length(object@fill.opacity) > 0) {
        fill.opacity <- object@fill.opacity
    }
    r1 <- 0
    r2 <- 0
    r3 <- max(h, w)/10
    r4 <- max(h, w)/10
    if (length(object@r4) > 0) {
        r4 <- object@r4
    }
    if (length(object@r3) > 0) {
        r3 <- object@r3
    }
    
    stroke.width <- 1
    if (length(object@shape$stroke.width) > 0) {
        stroke.width <- object@shape$stroke.width
    }
    stroke.opacity <- 1
    if (length(object@stroke.opacity) > 0) {
        stroke.opacity <- object@stroke.opacity
    }
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    # svg.path =
    # plot.path(object=object,d=d,id=id,fill.color='white',stroke.width=stroke.width,stroke.color
    # = 'black',fill.opacity = 1)
    
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    
    d2 <- generate.d.for.partly.rounded.corner.rect(x = x - w/2 + parameters.list$multimer.margin, 
        y = y - h/2 + parameters.list$multimer.margin, w = w, h = h, r1 = r1, 
        r2 = r2, r3 = r3, r4 = r4)
    svg.path2 <- plot.path(object = object, d = d2, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black", stroke.opacity = 1)
    
    # svg = paste(svg.path2,svg.path,sep='\n')
    svg <- paste(svg.path2, sep = "\n")
    
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        for (i in seq_len(length.out = length(user.data))) {
            fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                mol.type = "gene")
            glyph.data.i <- new("nucleic_acid_feature.sbgn")
            glyph.data.i@fill.color <- fill.color
            glyph.data.i@user.data <- c("no.user.data")
            glyph.data.i@y <- y
            glyph.data.i@h <- h
            glyph.data.i@w <- w/length(user.data)
            glyph.data.i@x <- x - w/2 + (i - 1/2) * glyph.data.i@w
            glyph.data.i@label <- " "
            glyph.data.i@compartment <- object@compartment
            glyph.data.i@id <- paste(id, "user.data", i, sep = "_")
            glyph.data.i@shape$stroke.width <- 0.5
            glyph.data.i@stroke.opacity <- 0
            glyph.data.i@fill.opacity <- 1
            if (i == 1) {
                glyph.data.i@r3 <- 0
            } else if (i == length(user.data)) {
                glyph.data.i@r4 <- 0
            } else {
                glyph.data.i@r4 <- 0
                glyph.data.i@r3 <- 0
            }
            svg <- paste(svg, plot.glyph(glyph.data.i), sep = "\n")
        }
        fill.opacity <- 0
    }
    d <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y - h/2, w = w, 
        h = h, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = 1, stroke.color = "black", fill.opacity = fill.opacity, stroke.opacity = 1)
    svg <- paste(svg, svg.path, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("nucleic_acid_feature.sbgn", contains = "glyph", slots = c(r4 = "numeric", r3 = "numeric"))
setMethod("plot.glyph", signature("nucleic_acid_feature.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    fill.color <- "white"
    if (length(object@fill.color) > 0) {
        fill.color <- object@fill.color
    }
    fill.opacity <- 1
    if (length(object@fill.opacity) > 0) {
        fill.opacity <- object@fill.opacity
    }
    label <- object@label
    label.location <- "center"
    r1 <- 0
    r2 <- 0
    r3 <- max(h, w)/10
    r4 <- max(h, w)/10
    if (length(object@r4) > 0) {
        r4 <- object@r4
    }
    if (length(object@r3) > 0) {
        r3 <- object@r3
    }
    stroke.width <- 1
    if (length(object@shape$stroke.width) > 0) {
        stroke.width <- object@shape$stroke.width
    }
    stroke.opacity <- 1
    if (length(object@stroke.opacity) > 0) {
        stroke.opacity <- object@stroke.opacity
    }
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    
    svg <- ""
    # if(user.data[1] != 'no.user.data'){
    if (!identical(user.data[1], "no.user.data")) {
        for (i in seq_len(length.out = length(user.data))) {
            fill.color <- color.from.value(user.data[i], global.parameters.list = parameters.list, 
                mol.type = "gene")
            glyph.data.i <- object
            glyph.data.i@fill.color <- fill.color
            glyph.data.i@fill.opacity <- 1
            glyph.data.i@user.data <- c("no.user.data")
            glyph.data.i@w <- w/length(user.data)
            glyph.data.i@x <- x - w/2 + (i - 1/2) * glyph.data.i@w
            glyph.data.i@label <- " "
            glyph.data.i@shape$stroke.width <- 0.5
            glyph.data.i@stroke.opacity <- 0
            if (i == 1) {
                glyph.data.i@r3 <- 0
            } else if (i == length(user.data)) {
                glyph.data.i@r4 <- 0
            } else {
                glyph.data.i@r4 <- 0
                glyph.data.i@r3 <- 0
            }
            svg <- paste(svg, plot.glyph(glyph.data.i), sep = "\n")
        }
        fill.opacity <- 0
    }
    d <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y - h/2, w = w, 
        h = h, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = fill.color, 
        stroke.width = stroke.width, stroke.color = "black", stroke.opacity = stroke.opacity, 
        fill.opacity = fill.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg, svg.path, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("perturbation.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("perturbation.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    d <- (paste("M", x, y, x + h/3, y + h/2, x, y + h, x + w, y + h, x + w - h/3, 
        y + h/2, x + w, y, "Z", sep = " "))
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("perturbing_agent.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("perturbing_agent.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    d <- (paste("M", x, y, x + h/3, y + h/2, x, y + h, x + w, y + h, x + w - h/3, 
        y + h/2, x + w, y, "Z", sep = " "))
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    # svg.ellipse = plot.ellipse(glyph=object,180,210,5,5,id,fill.color = 'red')
    svg <- paste(svg.path, svg.text, sep = "\n")
    # svg = paste(svg.path,svg.text,svg.ellipse,sep='\n')
    return(svg)
})

#########################################################################################################
setClass("annotation.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("annotation.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w
    y <- object@y - h
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    # stroke.opacity = object@stroke.opacity
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    label.location <- "center"
    d1 <- (paste("M", x, y, x, y + h, x + w/4, y + h, x + min(h, w)/10, y + 3 * h/2, 
        x + w/2, y + h, x + w, y + h, x + w, y + h/5, x + w - h/5, y, "Z", sep = " "))
    
    svg.path1 <- plot.path(object = object, d = d1, id = id, fill.color = "white", 
        stroke.width = 1, stroke.color = "black", stroke.opacity = stroke.opacity)
    
    d2 <- (paste("M", x + w - h/5, y + h/5, x + w, y + h/5, x + w - h/5, y, "Z", 
        sep = " "))
    svg.path2 <- plot.path(object = object, d = d2, id = id, fill.color = "black", 
        stroke.width = 1, stroke.color = "black", stroke.opacity = stroke.opacity)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location, 
        stroke.opacity = stroke.opacity)
    svg <- paste(svg.path1, svg.path2, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("existence.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("existence.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, fill.color = "white", 
        stroke.width = stroke.width)
    d <- (paste("M", x, y - h/2, "a", h/2, h/2, 0, 0, 1, 0, h, sep = " "))
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "black", 
        stroke.width = 1, stroke.color = "black")
    # svg.text =
    # plot.text(glyph=object,glyph.class=object@glyph.class,global.parameters.list=global.parameters.list,x=x,y=y,h=h,w=w,id=id,label=label,label.location=label.location)
    svg <- paste(svg.ellipse, svg.path, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("clone_simple_chemical.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("clone_simple_chemical.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- "clone"
    label <- object@label
    label.location <- "center"
    start.x <- x - w/2 + w/20
    start.y <- y + (h/2/(w/2) * ((w/2)^2 - (start.x - x)^2)^0.5)
    d <- (paste("M", start.x, start.y, "A", w/2, h/2, 1, 0, 0, x + w/2 - w/20, start.y, 
        sep = " "))
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "black", 
        stroke.width = stroke.width, stroke.color = "black")
    svg <- svg.path
    return(svg)
})

#########################################################################################################
setClass("clone.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("clone.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- "clone"
    label <- object@label
    label.location <- "center"
    r1 <- 0
    r2 <- 0
    r3 <- max(h, w)/10
    r4 <- max(h, w)/10
    d1 <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y + h/2 - h/5, 
        w = w, h = h/5, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path.bigbox <- plot.path(object = object, d = d1, id = id, fill.color = "black", 
        stroke.width = stroke.width, stroke.color = "black")
    
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y + h/2 - h/10, h = h, w = w, id = id, label = label, label.location = label.location, 
        color = "white")
    svg <- paste(svg.path.bigbox, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("clone_marker.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("clone_marker.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    r1 <- 0
    r2 <- 0
    r3 <- min(h, w)/10
    r4 <- min(h, w)/10
    d1 <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y - h/2, w = w, 
        h = h, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    svg.path.bigbox <- plot.path(object = object, d = d1, id = id, fill.color = "white", 
        stroke.width = 1, stroke.color = "black")
    d2 <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y - h/2 + 3/5 * 
        h, w = w, h = 2/5 * h, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path.innerbox <- plot.path(object = object, d = d2, id = id, fill.color = "gray", 
        stroke.width = stroke.width, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y - 1/5 * h, h = h, w = w, id = id, label = label, label.location = label.location)
    white.top.line.to.cover <- paste("M", x - w/2, y - h/2, x + w/2, y - h/2, sep = " ")
    svg.path.topWhiteLine <- plot.path(object = object, d = white.top.line.to.cover, 
        id = id, fill.color = "gray", stroke.width = 5, stroke.color = "white")
    svg <- paste(svg.path.bigbox, svg.path.innerbox, svg.path.topWhiteLine, svg.text, 
        sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("terminal.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("terminal.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    orientation <- object@orientation
    label.location <- "center"
    d <- (paste("M", x + h/3, y, x, y + h/2, x + h/3, y + h, x + w - h/3, y + h, 
        x + w, y + h/2, x + w - h/3, y, "Z", sep = " "))
    if (orientation == "right") {
        d <- (paste("M", x, y, x, y + h, x + w - w/3, y + h, x + w, y + h/2, x + 
            w - w/3, y, "Z", sep = " "))
    } else if (orientation == "left") {
        d <- (paste("M", x + w/3, y, x, y + h/2, x + w/3, y + h, x + w, y + h, x + 
            w, y, "Z", sep = " "))
    } else if (orientation == "up") {
        d <- (paste("M", x + w/2, y, x, y + h/3, x, y + h, x + w, y + h, x + w, y + 
            h/3, "Z", sep = " "))
    } else if (orientation == "down") {
        d <- (paste("M", x, y, x, y + h - h/3, x + w/2, y + h, x + w, y + h - h/3, 
            x + w, y, "Z", sep = " "))
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("tag.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("tag.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    orientation <- object@orientation
    label.location <- "center"
    d <- (paste("M", x + h/3, y, x, y + h/2, x + h/3, y + h, x + w - h/3, y + h, 
        x + w, y + h/2, x + w - h/3, y, "Z", sep = " "))
    if (orientation == "right") {
        d <- (paste("M", x, y, x, y + h, x + w - w/3, y + h, x + w, y + h/2, x + 
            w - w/3, y, "Z", sep = " "))
    } else if (orientation == "left") {
        d <- (paste("M", x + w/3, y, x, y + h/2, x + w/3, y + h, x + w, y + h, x + 
            w, y, "Z", sep = " "))
    } else if (orientation == "up") {
        d <- (paste("M", x + w/2, y, x, y + h/3, x, y + h, x + w, y + h, x + w, y + 
            h/3, "Z", sep = " "))
    } else if (orientation == "down") {
        d <- (paste("M", x, y, x, y + h - h/3, x + w/2, y + h, x + w, y + h - h/3, 
            x + w, y, "Z", sep = " "))
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = stroke.width, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("phenotype.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("phenotype.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x - w/2
    y <- object@y - h/2
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    d <- (paste("M", x + h/3, y, x + 5, y + h/2, x + h/3, y + h, x + w - h/3, y + 
        h, x + w - 5, y + h/2, x + w - h/3, y, "Z", sep = " "))
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = 1, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x + w/2, y = y + h/2, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path, svg.text, sep = "\n")
    return(svg)
})

######################################################################################################### 
setClass("association.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("association.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, fill.color = "black", 
        stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg <- paste(svg.ellipse, svg.wings, sep = "\n")
    return(svg)
})

######################################################################################################### 
setClass("interaction.sbgn.arc", contains = "glyph")
setMethod("plot.glyph", signature("interaction.sbgn.arc"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    Coordinates <- paste("x1=", x - w/2, " y1=", y + h/2, " x2=", x + w/2, " y2=", 
        y - h/2, "", sep = "\"")
    svg <- svg.ellipse
    return(svg)
})

#########################################################################################################
setClass("interaction.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("interaction.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    Coordinates <- paste("x1=", x - w/2, " y1=", y + h/2, " x2=", x + w/2, " y2=", 
        y - h/2, "", sep = "\"")
    svg.line <- plot.line(Coordinates, id)
    svg <- svg.ellipse
    return(svg)
})

#########################################################################################################
setClass("delay.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("delay.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- "t"
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    label.location <- "center"
    rx <- 0
    ry <- 0
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    svg.wings <- plot.LineWings(x, y, w, h, orientation, id)
    svg.wings <- ""
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.ellipse, svg.wings, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("source_and_sink.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("source_and_sink.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    if (length(object@orientation) > 0) {
        orientation <- object@orientation
    } else {
        orientation <- "horizontal"
    }
    stroke.width <- 1
    if (object@compartment == "free.node.within.compartment") {
        stroke.width <- 3
    }
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id, stroke.width = stroke.width)
    Coordinates <- paste("x1=", x - w/2, " y1=", y + h/2, " x2=", x + w/2, " y2=", 
        y - h/2, "", sep = "\"")
    svg.line <- plot.line(Coordinates, id)
    svg <- paste(svg.ellipse, svg.line, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("location.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("location.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    h <- w
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    svg.ellipse <- plot.ellipse(glyph = object, x, y, h, w, id)
    angle <- pi/6
    Coordinates1 <- paste("x1=", x - w/2 * cos(angle), " y1=", y + h/2 * sin(angle), 
        " x2=", x + w/2 * sin(angle), " y2=", y - h/2 * cos(angle), "", sep = "\"")
    svg.line1 <- plot.line(Coordinates1, id)
    Coordinates2 <- paste("x1=", x - w/2 * (cos(angle) - sin(angle))/2, " y1=", y + 
        h/2 * (sin(angle) - cos(angle))/2, " x2=", x + w/2 * sin(pi/4), " y2=", y + 
        h/2 * cos(pi/4), "", sep = "\"")
    svg.line2 <- plot.line(Coordinates2, id)
    svg <- paste(svg.ellipse, svg.line1, svg.line2, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("state_variable.ER.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("state_variable.ER.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    r1 <- h/2
    r2 <- r1
    r3 <- r1
    r4 <- r1
    d <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y - h/2, w = w, 
        h = h, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = 1, stroke.color = "black", fill.opacity = 1)
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path, svg.text, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("variable_value.sbgn", contains = "glyph")
setMethod("plot.glyph", signature("variable_value.sbgn"), function(object) {
    w <- object@w
    h <- object@h
    x <- object@x
    y <- object@y
    parameters.list <- object@parameters.list
    id <- object@id
    user.data <- object@user.data
    label <- object@label
    label.location <- "center"
    r1 <- h/2
    r2 <- r1
    r3 <- r1
    r4 <- r1
    d <- generate.d.for.partly.rounded.corner.rect(x = x - w/2, y = y - h/2, w = w, 
        h = h, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
    svg.path <- plot.path(object = object, d = d, id = id, fill.color = "white", 
        stroke.width = 1, stroke.color = "black")
    svg.text <- plot.text(glyph = object, glyph.class = object@glyph.class, parameters.list = parameters.list, 
        x = x, y = y, h = h, w = w, id = id, label = label, label.location = label.location)
    svg <- paste(svg.path, svg.text, sep = "\n")
    return(svg)
})


#########################################################################################################
###### Classes implementing the arc class and defining plot.arc generic function
#########################################################################################################
setMethod("plot.arc", signature("spline.arc"), function(object) {
    components <- object@components
    svg <- ""
    for (i in seq_len(length.out = length(components))) {
        component <- components[[i]]
        component@edge <- object@edge
        component@parameters.list <- object@parameters.list
        if (!is(component, "spline")) {
            component@edge$line.stroke.opacity <- 0
        }
        svg <- paste(svg, plot.arc(component), sep = "\n")
    }
    return(svg)
})

#########################################################################################################
spline.svg.template <- "<path d=\"M%s %s C %s %s, %s %s,%s %s\" stroke=\"%s\" stroke-opacity=\"%s\" stroke-width=\"%s\" fill-opacity=\"0\" fill=\"transparent\"/>"
setMethod("plot.arc", signature("spline"), function(object) {
    svg <- do.call(sprintf, c(list(spline.svg.template), c(object@spline.coords, 
        object@edge$line.stroke.color, object@edge$line.stroke.opacity, object@edge$line.width)))
    return(svg)
})

#########################################################################################################
setClass("absolute_stimulation.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("absolute_stimulation.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y - parameters.list$edge.tip.size, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y + parameters.list$edge.tip.size, sep = "")
    points2 <- paste(end.x - parameters.list$edge.tip.size * 3, ",", end.y, 
        " ", end.x - parameters.list$edge.tip.size * 6, ",", end.y - parameters.list$edge.tip.size, 
        " ", end.x - parameters.list$edge.tip.size * 6, ",", end.y + parameters.list$edge.tip.size, 
        sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, fill.color = "white", 
            atan((start.y - end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
        svg.triangle2 <- plot.polygon(object = object, points2, fill.color = "white", 
            atan((start.y - end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, fill.color = "white", 
            atan((start.y - end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
        svg.triangle2 <- plot.polygon(object = object, points2, fill.color = "white", 
            atan((start.y - end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, svg.triangle2, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("absolute_inhibition.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("absolute_inhibition.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    parameters.list <- object@parameters.list
    id <- object@id
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    d <- (paste("M", end.x - 1, end.y - parameters.list$edge.tip.size * 3/2, 
        end.x - 1, end.y + parameters.list$edge.tip.size * 3/2, sep = " "))
    svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
        stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
            end.x))/pi * 180, end.x, end.y)
    if (start.x < end.x) {
        d1 <- (paste("M", end.x - parameters.list$edge.tip.size * 3/4, end.y - 
            parameters.list$edge.tip.size * 3/2, end.x - parameters.list$edge.tip.size * 
            3/4, end.y + parameters.list$edge.tip.size * 3/2, sep = " "))
        svg.endLine1 <- plot.path(object = object, d = d1, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    } else {
        d1 <- (paste("M", end.x + parameters.list$edge.tip.size * 3/4, end.y - 
            parameters.list$edge.tip.size * 3/2, end.x + parameters.list$edge.tip.size * 
            3/4, end.y + parameters.list$edge.tip.size * 3/2, sep = " "))
        svg.endLine1 <- plot.path(object = object, d = d1, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.endLine, svg.endLine1, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("necessary_stimulation.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("necessary_stimulation.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    end.line.length <- parameters.list$edge.tip.size * 3/2
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y - end.line.length, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y + end.line.length, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
        d <- (paste("M", end.x - parameters.list$edge.tip.size * 4, end.y - 
            end.line.length - parameters.list$edge.tip.size * 3/4, end.x - 
            parameters.list$edge.tip.size * 4, end.y + end.line.length + parameters.list$edge.tip.size * 
            3/4, sep = " "))
        svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
        d <- (paste("M", end.x + parameters.list$edge.tip.size * 4, end.y - 
            end.line.length - parameters.list$edge.tip.size * 3/4, end.x + 
            parameters.list$edge.tip.size * 4, end.y + end.line.length + parameters.list$edge.tip.size * 
            3/4, sep = " "))
        svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, svg.endLine, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("equivalence_arc.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("equivalence_arc.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    svg <- svg.line
    return(svg)
})

######################################################################################################### 
setClass("consumption.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("consumption.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    svg <- svg.line
    return(svg)
})

#########################################################################################################
setClass("next.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("next.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    svg <- svg.line
    return(svg)
})

#########################################################################################################
setClass("production.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("production.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    end.size <- parameters.list$edge.tip.size * 2
    points <- paste(end.x, ",", end.y, " ", end.x - end.size, ",", end.y - end.size/2, 
        " ", end.x - end.size, ",", end.y + end.size/2, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "black", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "black", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("interaction.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("interaction.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    
    points.end <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y - parameters.list$edge.tip.size, " ", end.x - parameters.list$edge.tip.size * 
        2, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 3, ",", 
        end.y + parameters.list$edge.tip.size, sep = "")
    points.start <- paste(start.x, ",", start.y, " ", start.x + parameters.list$edge.tip.size * 
        3, ",", start.y - parameters.list$edge.tip.size, " ", start.x + parameters.list$edge.tip.size * 
        2, ",", start.y, " ", start.x + parameters.list$edge.tip.size * 3, 
        ",", start.y + parameters.list$edge.tip.size, sep = "")
    if (start.x < end.x) {
        svg.triangle.end <- plot.polygon(object = object, points.end, "black", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
        svg.triangle.start <- plot.polygon(object = object, points.start, "black", 
            atan((start.y - end.y)/(start.x - end.x))/pi * 180, start.x, start.y)
    } else {
        svg.triangle.end <- plot.polygon(object = object, points.end, "black", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
        svg.triangle.start <- plot.polygon(object = object, points.start, "black", 
            atan((start.y - end.y)/(start.x - end.x))/pi * 180 + 180, start.x, start.y)
    }
    svg <- paste(svg.line, svg.triangle.end, svg.triangle.start, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("assignment.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("assignment.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y - parameters.list$edge.tip.size, " ", end.x - parameters.list$edge.tip.size * 
        2, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 3, ",", 
        end.y + parameters.list$edge.tip.size, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "black", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "black", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("catalysis.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("catalysis.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    parameters.list <- object@parameters.list
    id <- object@id
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    angle <- atan((start.y - end.y)/(start.x - end.x))
    circle.r <- parameters.list$edge.tip.size
    if (start.x < end.x) {
        circle.center.x <- end.x - circle.r * cos(angle)
        circle.center.y <- end.y - circle.r * sin(angle)
    } else {
        circle.center.x <- end.x + circle.r * cos(angle)
        circle.center.y <- end.y + circle.r * sin(angle)
    }
    svg.ellipse <- plot.ellipse.edge.tip(circle.center.x, circle.center.y, circle.r * 
        2, circle.r * 2, id, fill.opacity = 0.8, object = object)
    svg <- paste(svg.line, svg.ellipse, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("positive_influence.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("positive_influence.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y - parameters.list$edge.tip.size, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y + parameters.list$edge.tip.size, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("stimulation.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("stimulation.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y - parameters.list$edge.tip.size, " ", end.x - parameters.list$edge.tip.size * 
        3, ",", end.y + parameters.list$edge.tip.size, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("negative_influence.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("negative_influence.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    if (start.x < end.x) {
        d <- (paste("M", end.x - parameters.list$inhibition.edge.end.shift, 
            end.y - parameters.list$edge.tip.size * 3/2, end.x - parameters.list$inhibition.edge.end.shift, 
            end.y + parameters.list$edge.tip.size * 3/2, sep = " "))
        svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    } else {
        d <- (paste("M", end.x + parameters.list$inhibition.edge.end.shift, 
            end.y - parameters.list$edge.tip.size * 3/2, end.x + parameters.list$inhibition.edge.end.shift, 
            end.y + parameters.list$edge.tip.size * 3/2, sep = " "))
        svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.endLine, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("inhibition.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("inhibition.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    
    inhibition.edge.end.shift <- parameters.list$inhibition.edge.end.shift
    if (start.x < end.x) {
        d1 <- (paste("M", start.x, end.y, end.x - inhibition.edge.end.shift, end.y, 
            sep = " "))
        svg.line <- plot.path(object = object, d = d1, id = id, fill.color = "white", 
            stroke.color = "white", stroke.width = 1, atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
        
        d <- (paste("M", end.x - inhibition.edge.end.shift, end.y - parameters.list$edge.tip.size * 
            3/2, end.x - inhibition.edge.end.shift, end.y + parameters.list$edge.tip.size * 
            3/2, sep = " "))
        svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    } else {
        d1 <- (paste("M", start.x, end.y, end.x + inhibition.edge.end.shift, end.y, 
            sep = " "))
        svg.line <- plot.path(object = object, d = d1, id = id, fill.color = "white", 
            stroke.color = "white", stroke.width = 1, atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
        
        d <- (paste("M", end.x + inhibition.edge.end.shift, end.y - parameters.list$edge.tip.size * 
            3/2, end.x + inhibition.edge.end.shift, end.y + parameters.list$edge.tip.size * 
            3/2, sep = " "))
        # object@edge$tip.stroke.color = 'purple'
        svg.endLine <- plot.path(object = object, d = d, id = id, fill.color = "black", 
            stroke.width = 1, stroke.color = "black", atan((start.y - end.y)/(start.x - 
                end.x))/pi * 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.endLine, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("unknown_influence.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("unknown_influence.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3/2, ",", end.y - parameters.list$edge.tip.size * 3/2, " ", end.x - 
        parameters.list$edge.tip.size * 3, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3/2, ",", end.y + parameters.list$edge.tip.size * 3/2, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("modulation.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("modulation.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    parameters.list <- object@parameters.list
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    points <- paste(end.x, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3/2, ",", end.y - parameters.list$edge.tip.size * 3/2, " ", end.x - 
        parameters.list$edge.tip.size * 3, ",", end.y, " ", end.x - parameters.list$edge.tip.size * 
        3/2, ",", end.y + parameters.list$edge.tip.size * 3/2, sep = "")
    if (start.x < end.x) {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180, end.x, end.y)
    } else {
        svg.triangle <- plot.polygon(object = object, points, "white", atan((start.y - 
            end.y)/(start.x - end.x))/pi * 180 + 180, end.x, end.y)
    }
    svg <- paste(svg.line, svg.triangle, sep = "\n")
    return(svg)
})

#########################################################################################################
setClass("logic_arc.sbgn.arc", contains = "arc")
setMethod("plot.arc", signature("logic_arc.sbgn.arc"), function(object) {
    start.x <- object@start.x
    start.y <- object@start.y
    end.x <- object@end.x
    end.y <- object@end.y
    id <- object@id
    if (length(object@stroke.opacity) != 0) 
        stroke.opacity <- object@stroke.opacity else stroke.opacity <- 1
    LineCoordinates <- paste("x1=", start.x, " y1=", start.y, " x2=", end.x, " y2=", 
        end.y, "", sep = "\"")
    svg.line <- plot.line(object = object, LineCoordinates, id, stroke.opacity = stroke.opacity)
    svg <- svg.line
    return(svg)
})

#########################################################################################################
