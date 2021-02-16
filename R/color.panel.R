
#########################################################################################################
# changed function name from 'color.from.01.value' to 'color.from.value'. Function used in plot.utilities.R
color.from.value <- function(value, mol.type, global.parameters.list) {
    
    if (mol.type == "gene") {
        min.value <- global.parameters.list$min.gene.value
        max.value <- global.parameters.list$max.gene.value
        mid.value <- global.parameters.list$mid.gene.value
        high <- global.parameters.list$col.gene.high
        mid <- global.parameters.list$col.gene.mid
        low <- global.parameters.list$col.gene.low
    } else if (mol.type == "cpd") {
        min.value <- global.parameters.list$min.cpd.value
        max.value <- global.parameters.list$max.cpd.value
        mid.value <- global.parameters.list$mid.cpd.value
        high <- global.parameters.list$col.cpd.high
        mid <- global.parameters.list$col.cpd.mid
        low <- global.parameters.list$col.cpd.low
    }
    
    if (is.na(value) | identical(value, "no.user.data")) {
        return("white")
    }
    if (is.na(mid)) {
        low <- col2rgb(low)
        high <- col2rgb(high)
        red <- (low[1, 1] + (high[1, 1] - low[1, 1]) * (value - min.value)/(max.value - 
            min.value))/255
        blue <- (low[2, 1] + (high[2, 1] - low[2, 1]) * (value - min.value)/(max.value - 
            min.value))/255
        green <- (low[3, 1] + (high[3, 1] - low[3, 1]) * (value - min.value)/(max.value - 
            min.value))/255
    } else {
        isodd <- value%%2 == 1
        low <- col2rgb(low)
        row.names(low) <- c("red", "blue", "green")
        mid <- col2rgb(mid)
        row.names(mid) <- c("red", "blue", "green")
        high <- col2rgb(high)
        row.names(high) <- c("red", "blue", "green")
        # lower <- floor(value/2) upper <- value - lower
        if (value < mid.value) {
            red <- (low[1, 1] + (mid[1, 1] - low[1, 1]) * (max(min.value, value) - 
                min.value)/(mid.value - min.value))/255  # if the value is out of color key range, use the color as min/max value
            blue <- (low[2, 1] + (mid[2, 1] - low[2, 1]) * (max(min.value, value) - 
                min.value)/(mid.value - min.value))/255
            green <- (low[3, 1] + (mid[3, 1] - low[3, 1]) * (max(min.value, value) - 
                min.value)/(mid.value - min.value))/255
        } else {
            red <- (mid[1, 1] + (high[1, 1] - mid[1, 1]) * (min(max.value, value) - 
                mid.value)/(max.value - mid.value))/255
            blue <- (mid[2, 1] + (high[2, 1] - mid[2, 1]) * (min(max.value, value) - 
                mid.value)/(max.value - mid.value))/255
            green <- (mid[3, 1] + (high[3, 1] - mid[3, 1]) * (min(max.value, value) - 
                mid.value)/(max.value - mid.value))/255
        }
    }
    return(rgb(red, blue, green))
}

#########################################################################################################
# color.panel(100,100,10)
color.panel <- function(x, y, mol.type, col.panel.w, global.parameters.list, if.insert.sbgn.link = TRUE) {
    
    if (mol.type == "gene") {
        min.value <- global.parameters.list$min.gene.value
        max.value <- global.parameters.list$max.gene.value
        mid.value <- global.parameters.list$mid.gene.value
    } else if (mol.type == "cpd") {
        min.value <- global.parameters.list$min.cpd.value
        max.value <- global.parameters.list$max.cpd.value
        mid.value <- global.parameters.list$mid.cpd.value
    }
    n.grid <- global.parameters.list$color.panel.n.grid
    if (n.grid%%2 == 1) {
        n.grid <- n.grid - 1
    }
    col.values <- seq(max.value, min.value, by = (min.value - max.value)/n.grid)
    panel.grid.w <- col.panel.w/n.grid
    
    col.panel.h <- col.panel.w/7
    y.loc.text <- y + (2 + 1/3) * col.panel.h
    y.loc.link.to.SBGN.notation <- y.loc.text + col.panel.h * 1.2
    font.size.color.panel <- col.panel.h
    alignment.baseline <- "baseline"
    
    svg <- ""
    for (i in seq_len(length.out = length(col.values))) {
        col <- color.from.value(col.values[i], global.parameters.list = global.parameters.list, 
            mol.type = mol.type)
        x.loc.text <- x - (i - 0.5) * panel.grid.w
        svg.rect <- plot.rectangle(x = x - i * panel.grid.w, y = y, h = col.panel.h, 
            w = panel.grid.w + 1, id = "col.pan", rx = 0, ry = 0, stroke.width = 0, 
            fill.color = col, stroke.opacity = 0, fill.opacity = 1)
        svg <- paste(svg, svg.rect, sep = "\n")
    }
    
    j <- 0
    if (min.value > 100) 
        min.value <- formatC(min.value, format = "e", digits = 0)
    if (mid.value > 100) 
        mid.value <- formatC(mid.value, format = "e", digits = 0)
    if (max.value > 100) 
        max.value <- formatC(max.value, format = "e", digits = 0)
    for (x.loc in seq(x - col.panel.w - panel.grid.w, x, length.out = 11)) {
        j <- j + 1
        LineCoordinates <- paste("x1=", x.loc, " y1=", y, " x2=", x.loc, " y2=", 
            y + col.panel.h * 5/4, "", sep = "\"")
        svg.tick <- plot.line(LineCoordinates, id = "tick", stroke.width = col.panel.w/100)
        svg <- paste(svg, svg.tick, sep = "\n")
        if (j == 11) {
            svg.text.high <- sprintf(template.text, x.loc, y.loc.text, "color.panel.max.value", 
                font.size.color.panel, "middle", alignment.baseline, 1, alignment.baseline, 
                "black", as.character(max.value))
        } else if (j == 6) {
            svg.text.mid <- sprintf(template.text, x.loc, y.loc.text, "color.panel.mid.value", 
                font.size.color.panel, "middle", alignment.baseline, 1, alignment.baseline, 
                "black", as.character(mid.value))
            if (if.insert.sbgn.link) {
                svg.text.link.to.SBGN.panel <- sprintf(template.text, x.loc, y.loc.link.to.SBGN.notation, 
                  "Link to SBGN notation", font.size.color.panel, "middle", alignment.baseline, 
                  1, alignment.baseline, "blue", "Link to SBGN notation")
                svg.text.link.to.SBGN.panel <- paste(link.to.sbgn.start, svg.text.link.to.SBGN.panel, 
                  link.to.sbgn.end, sep = "\n")
            } else {
                svg.text.link.to.SBGN.panel <- ""
            }
        } else if (j == 1) {
            svg.text.low <- sprintf(template.text, x.loc, y.loc.text, "color.panel.low.value", 
                font.size.color.panel, "middle", alignment.baseline, 1, alignment.baseline, 
                "black", as.character(min.value))
        }
    }
    
    svg <- paste(svg, svg.text.high, svg.text.low, svg.text.mid, 
                 svg.text.link.to.SBGN.panel, sep = "\n")
    return(svg)
}

#########################################################################################################
find.col.panel.range <- function(user.data, max.gene.value, mid.gene.value, min.gene.value, 
                                 max.cpd.value, mid.cpd.value, min.cpd.value) {
    
    if (!"max.gene" %in% names(user.data)) {
        user.data[["max.gene"]] <- 1
        user.data[["min.gene"]] <- -1
        if.has.gene.data <- FALSE
    } else {
        if.has.gene.data <- TRUE
    }
    if (!"max.cpd" %in% names(user.data)) {
        user.data[["max.cpd"]] <- 1
        user.data[["min.cpd"]] <- -1
        if.has.cpd.data <- FALSE
    } else {
        if.has.cpd.data <- TRUE
    }
    if (is.null(max.gene.value)) 
        max.gene.value <- signif(user.data[["max.gene"]], 2) + 0.1
    if (is.null(min.gene.value)) 
        min.gene.value <- signif(user.data[["min.gene"]], 2) - 0.1
    if (is.null(max.cpd.value)) 
        max.cpd.value <- signif(user.data[["max.cpd"]], 2) + 0.1
    if (is.null(min.cpd.value)) 
        min.cpd.value <- signif(user.data[["min.cpd"]], 2) - 0.1
    if (is.null(mid.gene.value)) {
        mid.gene.value <- signif(median(c(max.gene.value, min.gene.value)), 2) + 
            0.1
    }
    if (is.null(mid.cpd.value)) {
        mid.cpd.value <- signif(median(c(max.cpd.value, min.cpd.value)), 2) + 0.1
    }
    return(list(if.has.gene.data = if.has.gene.data, if.has.cpd.data = if.has.cpd.data, 
                max.gene.value = max.gene.value, mid.gene.value = mid.gene.value, min.gene.value = min.gene.value, 
                max.cpd.value = max.cpd.value, mid.cpd.value = mid.cpd.value, min.cpd.value = min.cpd.value))
}

#########################################################################################################
find.col.panel.position.and.plot <- function(y.margin, global.parameters.list, if.has.gene.data, 
                                             if.has.cpd.data, parse.glyph.out.list, max.x, max.y, min.x, min.y) {
    
    add.scale.factor <- 14
    if(global.parameters.list$if.scale.compartment.font.size) {
        y.margin <- y.margin - (add.scale.factor * 100*global.parameters.list$node.width.adjust.factor.compartment)
    }
    if (global.parameters.list$if.scale.complex.font.size){
        y.margin <- y.margin - (add.scale.factor * 100*global.parameters.list$node.width.adjust.factor.complex)
    }
    if (global.parameters.list$font.size.scale.compartment != 1.6) { # 1.6 default for font.size.scale.compartment
        y.margin <- y.margin - (add.scale.factor*global.parameters.list$font.size.scale.compartment)
    }
    if (global.parameters.list$font.size.scale.complex != 1.1) { # 1.1 default for font.size.scale.complex
        y.margin <- y.margin - (add.scale.factor*global.parameters.list$font.size.scale.complex)
    }
    if (global.parameters.list$font.size != 3) { # 3 default for font.size
        y.margin <- y.margin - (4*global.parameters.list$font.size)
    }
    if (global.parameters.list$pathway.name.font.size != 1) { # 1 default for pathway.name.font.size
        y.margin <- y.margin - (add.scale.factor*global.parameters.list$pathway.name.font.size)
    }
    
    col.panel.w <- y.margin * 2.4 * global.parameters.list$color.panel.scale # scale color.panel
       
    if (global.parameters.list$key.pos == "none") {
        print("no color panel and pathway name will be plotted")
        col.panel.svg <- ""
        col.panel.svg.chemical <- ""
        col.panel.h <- 0
        col.panel.y <- 0
    } 
    else {
        if (all(if.has.gene.data, if.has.cpd.data)) {
            col.panel.h <- col.panel.w
        } else {
            col.panel.h <- col.panel.w/2
        }
        result.list <- find.key.pos(parse.glyph.out.list, col.panel.w, col.panel.h, 
                                    ymargin = y.margin, max.y, max.x, min.x, min.y, key.pos = global.parameters.list$key.pos, 
                                    space.between.color.panel.and.entity = global.parameters.list$space.between.color.panel.and.entity)
        col.panel.x <- result.list$col.panel.x
        col.panel.y <- result.list$col.panel.y
        
        # plot color panel
        if (if.has.gene.data & if.has.cpd.data) {
            col.panel.svg <- color.panel(x = col.panel.x, y = col.panel.y, mol.type = "gene", 
                                         col.panel.w = col.panel.w, global.parameters.list = global.parameters.list, 
                                         if.insert.sbgn.link = FALSE)
            col.panel.svg.chemical <- color.panel(x = col.panel.x, y = col.panel.y + 
                                                      col.panel.w * 0.45, mol.type = "cpd", col.panel.w = col.panel.w, 
                                                  global.parameters.list = global.parameters.list)
        } else if (if.has.cpd.data) {
            col.panel.svg <- ""
            col.panel.svg.chemical <- color.panel(x = col.panel.x, y = col.panel.y, 
                                                  mol.type = "cpd", col.panel.w = col.panel.w, global.parameters.list = global.parameters.list)
        } else if (if.has.gene.data) {
            col.panel.svg.chemical <- ""
            col.panel.svg <- color.panel(x = col.panel.x, y = col.panel.y, mol.type = "gene", 
                                         col.panel.w = col.panel.w, global.parameters.list = global.parameters.list)
        } else {
            col.panel.svg.chemical <- ""
            col.panel.svg <- ""
        }
    }
    col.panel.svg <- paste(col.panel.svg, col.panel.svg.chemical, sep = "\n")
    return(list(col.panel.svg = col.panel.svg, col.panel.w = col.panel.w, 
                col.panel.y = col.panel.y, col.panel.h = col.panel.h))
}

#########################################################################################################
# new version of find.key.pos. If key.pos equals 'topleft' or 'bottomleft', the output is messed up
# since we have stamps and lables in top/bottom left
# so only accept key.pos of 'topright' or 'bottomright'. if other, default to 'topright'
# and show message to user 
find.key.pos <- function(parse.glyph.out.list, col.panel.w, col.panel.h, ymargin, 
                         max.y, max.x, min.x, min.y, key.pos, space.between.color.panel.and.entity) {
    
    #dist.to.top <- ymargin + min.y # original
    dist.to.top <- ymargin + 10

    glyph.coors <- parse.glyph.out.list$glyph.coors
    
    if(!key.pos %in% c("bottomleft", "topleft", "bottomright", "topright", "none")) {
        stop("Not a valid input for key.pos arguemnt. Please check renderSbgn function documentation.")
    }
    if(key.pos %in% c("bottomleft", "topleft")) {
        print("key.pos cannot be 'topleft' or 'bottomleft' because it will interfere with labels in those positions. key.pos value set to 'topright'")
        key.pos <- "topright"
    }
    
    if(key.pos == "none"){
        return("no color key")
    }
    else if (key.pos == "bottomright") {
        if.nodes.within.col.panel.best.area <- glyph.coors[, "xw"] > (max.x - col.panel.w) - 
            space.between.color.panel.and.entity & glyph.coors[, "yh"] > max.y - 
            col.panel.h - space.between.color.panel.and.entity
    } else if (key.pos == "topright") {
        if.nodes.within.col.panel.best.area <- glyph.coors[, "xw"] > (max.x - col.panel.w) - 
            space.between.color.panel.and.entity & glyph.coors[, "y"] < dist.to.top + 
            col.panel.h + space.between.color.panel.and.entity
    } 
    
    if (any(if.nodes.within.col.panel.best.area)) {
        nodes.within.col.panel.best.area <- glyph.coors[if.nodes.within.col.panel.best.area, 
        ]
        if (is.vector(nodes.within.col.panel.best.area)) {
            nodes.within.col.panel.best.area <- as.matrix(t(nodes.within.col.panel.best.area))
        }
        if (key.pos == "bottomright") {
            col.panel.y <- max(nodes.within.col.panel.best.area[, "yh"]) + col.panel.h/10
            col.panel.x <- max.x
        } else if (key.pos == "topright") {
            #col.panel.y <- max(0 + dist.to.top/10, min(nodes.within.col.panel.best.area[, "y"]) - col.panel.h - col.panel.h/10)
            col.panel.y <- dist.to.top/2 + 10
            col.panel.x <- max.x
        } 
    } else {
        if (key.pos == "bottomright") {
            col.panel.y <- max.y - col.panel.h + col.panel.h/10
            col.panel.x <- max.x
        } else if (key.pos == "topright") {
            col.panel.y <- dist.to.top
            col.panel.x <- max.x
        } 
    }
    return(list(col.panel.x = col.panel.x, col.panel.y = col.panel.y))
}

#########################################################################################################
### old version of function that allows key.pos to be set to 'bottomleft', 'topleft', 'none'
### replaced with new version (above) that only allows key.pos to be set to 'topright' or 'bottomright'
### since other locations contain stamps and labels which interfere with key.pos if set to those locations.
# find.key.pos <- function(parse.glyph.out.list, col.panel.w, col.panel.h, ymargin, 
#     max.y, max.x, min.x, min.y, key.pos, space.between.color.panel.and.entity) {
#     
#     dist.to.top <- ymargin + min.y
#     glyph.coors <- parse.glyph.out.list$glyph.coors
#     if (key.pos == "none") {
#         return("no color key")
#     } else if (key.pos == "bottomright") {
#         if.nodes.within.col.panel.best.area <- glyph.coors[, "xw"] > (max.x - col.panel.w) - 
#             space.between.color.panel.and.entity & glyph.coors[, "yh"] > max.y - 
#             col.panel.h - space.between.color.panel.and.entity
#     } else if (key.pos == "topright") {
#         # if.nodes.within.col.panel.best.area =
#         # glyph.coors[,'xw']>(max.x-col.panel.w)-space.between.color.panel.and.entity &
#         # glyph.coors[,'y'] < ymargin + col.panel.h +
#         # space.between.color.panel.and.entity
#         if.nodes.within.col.panel.best.area <- glyph.coors[, "xw"] > (max.x - col.panel.w) - 
#             space.between.color.panel.and.entity & glyph.coors[, "y"] < dist.to.top + 
#             col.panel.h + space.between.color.panel.and.entity
#     } else if (key.pos == "bottomleft") {
#         if.nodes.within.col.panel.best.area <- glyph.coors[, "x"] < (col.panel.w + 
#             min.x) + space.between.color.panel.and.entity & glyph.coors[, "yh"] > 
#             max.y - col.panel.h - space.between.color.panel.and.entity
#     } else if (key.pos == "topleft") {
#         if.nodes.within.col.panel.best.area <- glyph.coors[, "x"] < (col.panel.w + 
#             min.x) + space.between.color.panel.and.entity & glyph.coors[, "y"] < 
#             dist.to.top + col.panel.h + 4
#     }
#     if (any(if.nodes.within.col.panel.best.area)) {
#         nodes.within.col.panel.best.area <- glyph.coors[if.nodes.within.col.panel.best.area, 
#             ]
#         if (is.vector(nodes.within.col.panel.best.area)) {
#             nodes.within.col.panel.best.area <- as.matrix(t(nodes.within.col.panel.best.area))
#         }
#         if (key.pos == "bottomright") {
#             col.panel.y <- max(nodes.within.col.panel.best.area[, "yh"]) + col.panel.h/10
#             col.panel.x <- max.x
#         } else if (key.pos == "topright") {
#             col.panel.y <- max(0 + dist.to.top/10, min(nodes.within.col.panel.best.area[, 
#                 "y"]) - col.panel.h - col.panel.h/10)
#             col.panel.x <- max.x
#         } else if (key.pos == "bottomleft") {
#             col.panel.y <- max(nodes.within.col.panel.best.area[, "yh"]) + col.panel.h/10
#             col.panel.x <- col.panel.w + min.x
#         } else if (key.pos == "topleft") {
#             col.panel.y <- max(0 + dist.to.top/10, min(nodes.within.col.panel.best.area[, 
#                 "y"]) - col.panel.h - col.panel.h/10)
#             col.panel.x <- col.panel.w + min.x
#         }
#     } else {
#         if (key.pos == "bottomright") {
#             col.panel.y <- max.y - col.panel.h + col.panel.h/10
#             col.panel.x <- max.x
#         } else if (key.pos == "topright") {
#             col.panel.y <- dist.to.top
#             col.panel.x <- max.x
#         } else if (key.pos == "bottomleft") {
#             col.panel.y <- max.y - col.panel.h + col.panel.h/10
#             col.panel.x <- col.panel.w + min.x
#         } else if (key.pos == "topleft") {
#             col.panel.y <- dist.to.top
#             col.panel.x <- col.panel.w + min.x
#         }
#     }
#     return(list(col.panel.x = col.panel.x, col.panel.y = col.panel.y))
# }

#########################################################################################################
