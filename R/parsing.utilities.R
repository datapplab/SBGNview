
#########################################################################################################
# data('pathways.info','sbgn.xmls','pathway.completeness.cutoff.info')
utils::globalVariables(c("pathways.info", "mapped.ids", "sbgn.xmls", 
                         "pathway.completeness.cutoff.info", "bods"))

#' @import igraph 
#' @import xml2
#' @import rsvg
#' @import pathview
#' @import Rdpack
#' @import rmarkdown
#' @import knitr
#' @import SBGNview.data
#' @import grDevices
#' @import methods
#' @importFrom stats median runif var
#' @importFrom AnnotationDbi select columns
#' @importFrom SummarizedExperiment assays
#' @import utils
#' @import httr
#' @import KEGGREST

#########################################################################################################
# Parse input SBGN-ML file from local or download the file
# get pathway database, and SBGN gene and compound ID type
parse.input.sbgn <- function(input.sbgn, output.file, show.pathway.name, sbgn.dir, 
                             sbgn.gene.id.type, sbgn.cpd.id.type, sbgn.id.attr, 
                             SBGNview.data.folder = "./SBGNview.tmp.data") {
    database <- NULL
    output.file.sbgn <- paste(output.file, "_", input.sbgn, sep = "")
    input.sbgn.full.path <- paste(sbgn.dir, input.sbgn, sep = "/")
    if (!file.exists(input.sbgn.full.path)) {
        # if the SBGN-ML file is not downloaded, in this case, input.sbgn is pathway.id
        message(input.sbgn, " does not look like a local file.\nUsing it as pathway ID in 'pathways.info'\n\n ")
        is.input <- pathways.info[, "pathway.id"] == input.sbgn
        if (!any(is.input)) {
            stop("The input ID is not included in the pre-generated SBGN-ML file database! Please check 'data(\"pathways.info\")' for supported pathway IDs!\n\n")
        } else {
            database <- pathways.info[is.input, "database"]
            sub.database <- pathways.info[is.input, "sub.database"]
            
            pathway.name <- pathways.info[is.input, "pathway.name"]
            pathway.name <- gsub("[\\&\\/\\*\\?\\<\\>\\|\\:\"]", "_", pathway.name)
            pathway.name.on.graph <- paste(pathway.name, sub.database, input.sbgn, sep = "::")
           
            if (show.pathway.name) {
                output.file.sbgn <- paste(output.file.sbgn, "_", pathway.name, sep = "")
            }
            
            sbgn.gene.id.type <- pathways.info[is.input, "macromolecule.ID.type"]
            sbgn.cpd.id.type <- pathways.info[is.input, "simple.chemical.ID.type"]
            input.sbgn.full.path <- downloadSbgnFile(pathway.id = input.sbgn, download.folder = sbgn.dir)
            message("SBGN-ML files downloaded to: ", input.sbgn.full.path)
        }
    } else {
        # if the SBGN-ML file is downloaded. In this way input.sbgn is file.name
        message(input.sbgn, " looks like a local file. Using it directly.\n ")
        if (input.sbgn %in% pathways.info[, "file.name"]) {
            # if this local SBGN-ML file is our pregenerated file instead of user's own file
            is.input <- pathways.info[, "file.name"] == input.sbgn
            if (!any(is.input)) {
                stop("The input file name is not included in the pre-generated SBGN-ML file database! Please check 'data(\"pathways.info\")' for supported pathway file names!\n\n")
            }
            database <- pathways.info[is.input, "database"]
            sbgn.gene.id.type <- pathways.info[is.input, "macromolecule.ID.type"]
            sbgn.cpd.id.type <- pathways.info[is.input, "simple.chemical.ID.type"]
            message("Using pre-generated SBGN-ML files and ID mapping tables.\n", 
                "ID mapping tables will be downloaded into folder: ", SBGNview.data.folder)
        } else {
            # if this is user's own file
            message("Using user defined pathway SBGN-ML file!!\n")
            database <- "user.data"
        }
        
        pathway.name.on.graph <- paste(input.sbgn, database, sep = "::")
    }
    if (database == "MetaCrop") {
        sbgn.id.attr <- "label"
    } else if (database == "MetaCyc") {
        sbgn.id.attr <- "id"
    } else if (database == "pathwayCommons") {
        sbgn.id.attr <- "id"
    }
    if.file.in.collection <- "Rendered by SBGNview"
    pathway.name.on.graph <- list(pathway.name.on.graph = pathway.name.on.graph, 
                                  if.file.in.collection = if.file.in.collection)
    
    return(list(input.sbgn.full.path = input.sbgn.full.path, output.file.sbgn = output.file.sbgn, 
                sbgn.gene.id.type = sbgn.gene.id.type, sbgn.cpd.id.type = sbgn.cpd.id.type, 
                pathway.name.on.graph = pathway.name.on.graph, if.file.in.collection = if.file.in.collection, 
                sbgn.id.attr = sbgn.id.attr, database = database))
}

#########################################################################################################
# parse input user data by changing ID tpype to SBGN ID type for gene and compound data
# keeps record of data ID that are already converted 
parse.omics.data <- function(gene.data, cpd.data, input.sbgn.full.path, database, 
                             user.data.recorded, gene.id.type, cpd.id.type, id.mapping.gene, 
                             id.mapping.cpd, node.sum, org, sbgn.gene.id.type, sbgn.cpd.id.type, 
                             simulate.data, SBGNview.data.folder = "./SBGNview.tmp.data") {
    
    
    if (!is.null(gene.data) & is.null(sbgn.gene.id.type)) {
        stop("Must provide 'sbgn.gene.id.type'!")
    }
    if (!is.null(cpd.data) & is.null(sbgn.cpd.id.type)) {
        stop("Must provide 'sbgn.cpd.id.type'!")
    }
    
    if (is(gene.data, "SummarizedExperiment")) {
        gene.data <- SummarizedExperiment::assays(gene.data)$counts
    }
    if (is(cpd.data, "SummarizedExperiment")) {
        cpd.data <- SummarizedExperiment::assays(cpd.data)$counts
    }
    if (simulate.data) {
        gene.data.converted <- simulate.user.data(sbgn.file = input.sbgn.full.path, n.samples = 3)
        mm=range(gene.data.converted, na.rm=TRUE)
        gene.data.converted <- list(gene.data.converted, max.gene=mm[2], min.gene=mm[1], max.cpd=mm[2], min.cpd=mm[1])
        cpd.data.converted <- gene.data.converted[1]
        
    } else {
        # convert input IDs to glyph IDs SBGNview can render multiple SBGN-ML files but
        # the input omics data only need to be processed once for each database(each
        # glyph ID type)
        if (database %in% names(user.data.recorded)) {
            print("Data ID already converted!!")
            gene.data.converted <- user.data.recorded[[database]][["gene.data"]]
            cpd.data.converted <- user.data.recorded[[database]][["cpd.data"]]
        } else {
            gene.data.converted <- gene.data
            if (!is.null(gene.data)) {
                # if user provided gene data, we need its ID type
                if (is.na(gene.id.type)) {
                    stop("No omics gene ID type specified. Please set parameter gene.id.type!")
                }
                if (is.vector(gene.data)) {
                    if (is.character(gene.data[1])) {
                        # if we need to generate count data from input data
                        gene.data <- as.matrix(table(gene.data))
                    } else {
                        gene.data <- as.matrix(gene.data)
                    }
                }
                if (gene.id.type != sbgn.gene.id.type) {
                    gene.data.converted <- changeDataId(gene.data, input.type = gene.id.type, 
                                                        output.type = sbgn.gene.id.type, sum.method = node.sum, mol.type = "gene", 
                                                        org = org, id.mapping.table = id.mapping.gene, SBGNview.data.folder = SBGNview.data.folder)
                    # if (identical(gene.data.converted, "no.id.mapped")) {
                    #     warning("no id mapped!")
                    # }
                }
                mm=range(gene.data.converted, na.rm=TRUE)
                gene.data.converted <- list(gene.data.converted, max.gene=mm[2], min.gene=mm[1])
            }
            
            cpd.data.converted <- cpd.data
            if (!is.null(cpd.data)) {
                if (is.na(cpd.id.type)) {
                    stop("No omics compound ID type specified. Please set parameter cpd.id.type!")
                } else {
                    if (is.vector(cpd.data)) {
                        if (is.character(cpd.data[1])) {
                            cpd.data.converted <- as.matrix(table(cpd.data))
                        } else {
                            cpd.data.converted <- as.matrix(cpd.data)
                        }
                    }
                    if (cpd.id.type != sbgn.cpd.id.type) {
                        cpd.data.converted <- changeDataId(cpd.data.converted, input.type = cpd.id.type, 
                                                           output.type = sbgn.cpd.id.type, sum.method = node.sum, mol.type = "cpd", 
                                                           id.mapping.table = id.mapping.cpd, SBGNview.data.folder = SBGNview.data.folder)
                        # if (cpd.data.converted == "no.id.mapped") {
                        #     warning("no id mapped!")
                        # }
                    }
                    mm=range(cpd.data.converted, na.rm=TRUE)
                    cpd.data.converted <- list(cpd.data.converted, max.cpd=mm[2], min.cpd=mm[1])
                }
            }
            user.data.recorded[[database]] <- list(gene.data=gene.data.converted, cpd.data=cpd.data.converted)            
        }
    }
    user.data <- c(gene.data.converted, cpd.data.converted)  # input data is a list, name is sbgn entity id, content is a vector of data values. Therefore the number of values can be different
    return(list(user.data = user.data, user.data.recorded = user.data.recorded))
}

#########################################################################################################
# Use users input arcs list or get spline information from sbgn file, and create objects of spline.arc class. 
# generate svg code for plotting splines
parse.splines <- function(sbgn.xml, glyphs, if.plot.svg = TRUE, y.margin = 0, 
                          global.parameters.list, arcs.user = list()) {
    
    if (length(arcs.user) == 0) {
        svg.splines <- ""
        splines.list <- list()
        arc.splines <- xml2::xml_find_all(sbgn.xml, ".//edge.spline.info")
        arc.splines <- xml2::xml_find_all(arc.splines[[length(arc.splines)]], ".//arc.spline")  # if there are more than one version of routed edges, we use the latest version
        
        #splines <- rep(NULL, times = length(arc.splines))
        splines <- list()
        for (i in seq_along(arc.splines)) {
            arc.spline <- arc.splines[[i]]
            spline.info <- xml2::xml_attrs(arc.spline)  # extract attributes of this glyph
            spline.arc <- new("spline.arc")
            spline.arc@id <- paste(spline.info["id"], spline.info["class"], sep = "==")
            # spline.arc@id = spline.info['id']
            spline.arc@source <- spline.info["source"]
            spline.arc@target <- spline.info["target"]
            spline.arc@arc.class <- spline.info["class"]
            #spline.arc@parameters.list <- global.parameters.list
            spline.arc@parameters.list <- list()
            spline.arc@edge <- list(line.stroke.color = "black", line.stroke.opacity = 1, 
                                    line.width = 2, tip.stroke.opacity = 1, tip.stroke.color = "black", 
                                    tip.stroke.width = 1, tip.fill.color = "black", tip.fill.opacity = 1)
            spline.class <- gsub(" ", "_", spline.info["class"])
            
            components <- get.spline.compo(arc.spline, spline.info, y.margin)
            
            spline.arc@components <- components
            splines.list[[spline.arc@id]] <- spline.arc
            if (if.plot.svg) {
                spline.arc.svg <- plot.arc(spline.arc)
            } else {
                spline.arc.svg <- ""
            }
            splines[i] <- spline.arc.svg
        }
        svg.splines <- paste(splines, collapse = "\n")
    } else {
        svg.splines <- ""
        splines.list <- arcs.user
        for (i in seq_len(length.out = length(arcs.user))) {
            #spline.arc <- arcs.user[[i]]
            if (if.plot.svg) {
                spline.arc.svg <- plot.arc(arcs.user[[i]])
            } else {
                spline.arc.svg <- ""
            }
            svg.splines <- paste(svg.splines, spline.arc.svg, sep = "\n")
        }
    }
    return(list(svg.splines = svg.splines, splines.list = splines.list))
}

#########################################################################################################
# get child attributes of arc splines
get.spline.compo <- function(arc.spline, spline.info, y.margin) {
    
    children <- xml2::xml_children(arc.spline)
    components <- list()
    for (i in seq_len(length(children))) {
        child <- children[i]
        child.attrs <- xml2::xml_attrs(child)[[1]]
        spline.i <- 0
        if (xml2::xml_name(child) == "source.arc") {
            # state variables has different shape in PD and ER, so we need to switch the
            # shape when needed for interaction arcs, they have shapes at both ends, so we
            # need to change it: use 'assignment' arc to just plot one shape; also for source
            # arc, the orientation should be reversed
            if (spline.info["class"] == "interaction") {
                child.attrs["class"] <- "assignment"
            }
            arc <- new(paste(child.attrs["class"], ".sbgn.arc", sep = ""), 
                       id = paste(spline.info["id"], "source.arc", sep = ":"), start.x = as.numeric(child.attrs["start.x"]), 
                       start.y = as.numeric(child.attrs["start.y"]) + y.margin, end.x = as.numeric(child.attrs["end.x"]), 
                       end.y = as.numeric(child.attrs["end.y"]) + y.margin, stroke.opacity = 0)
        } else if (xml2::xml_name(child) == "target.arc") {
            # state variables has different shape in PD and ER, so we need to switch the
            # shape when needed for interaction arcs, they have shapes at both ends, so we
            # need to change it: use 'assignment' arc to just plot one shape; also for source
            # arc, the orientation should be reversed
            if (spline.info["class"] == "interaction") {
                child.attrs["class"] <- "assignment"
            }
            arc <- new(paste(child.attrs["class"], ".sbgn.arc", sep = ""), 
                       id = paste(spline.info["id"], "target.arc", sep = ":"), start.x = as.numeric(child.attrs["start.x"]), 
                       start.y = as.numeric(child.attrs["start.y"]) + y.margin, end.x = as.numeric(child.attrs["end.x"]), 
                       end.y = as.numeric(child.attrs["end.y"]) + y.margin, stroke.opacity = 0)
        } else if (xml2::xml_name(child) == "spline") {
            spline.i <- spline.i + 1
            if.y <- grepl("\\.y", names(child.attrs))
            child.attrs[if.y] <- as.numeric(child.attrs[if.y]) + y.margin
            arc <- new("spline", id = paste(spline.info["id"], "spline", spline.i, 
                                            sep = ":"), spline.coords = child.attrs)
        }
        components <- c(components, arc)
    }
    return(components)
}

#########################################################################################################
# Use users input arcs list or get arcs information from sbgn file, and create objects of arc class. 
# generate svg code for plotting arcs
parse.arcs <- function(sbgn.xml, glyphs, if.plot.svg = TRUE, y.margin = 0, 
                       global.parameters.list, arcs.user = list()) {
    
    arcs.list <- list()
    if (length(arcs.user) == 0) {
        edge.paras <- list(line.stroke.color = "black", line.stroke.opacity = 1, 
                           line.width = 2, tip.stroke.opacity = 1, tip.stroke.color = "black", tip.stroke.width = 1, 
                           tip.fill.color = "black", tip.fill.opacity = 1)
        all.arcs <- xml2::xml_find_all(sbgn.xml, ".//arc")
        for (i in seq_len(length(all.arcs))) {
            arc <- all.arcs[i]
            arc.info <- xml2::xml_attrs(arc)[[1]]  # extract attributes of this glyph
            if (is.na(arc.info["id"])) {
                arc.info["id"] <- paste(arc.info["source"], arc.info["target"], sep = "->")
            }
            arc.class <- gsub(" ", "_", arc.info["class"])
            arc.segments <- get.arc.segments(arc, arcs.list, arc.class, y.margin, 
                                             arc.info, edge.paras, global.parameters.list, glyphs)
            svg.arc <- arc.segments$svg.arc
            arcs.list <- arc.segments$arcs.list
        }
    } else {
        arcs.list <- arcs.user
        svg.arc <- ""
        for (i in seq_len(length.out = length(arcs.list))) {
            svg.arc <- paste(svg.arc, plot.arc(arcs.list[[i]]), sep = "\n")
        }
    }
    return(list(svg.arc = svg.arc, arcs.list = arcs.list))
}

#########################################################################################################
# get child attributes of arcs
get.arc.segments <- function(arc, arcs.list, arc.class, y.margin, arc.info, edge.paras, 
                             global.parameters.list, glyphs) {
    
    svg.arc <- ""
    arc.line <- c(arc.info["source"], arc.info["target"], arc.info["id"], start.x = "", 
                  start.y = "", end.x = "", end.y = "")
    children <- xml2::xml_children(arc)
    for (i in seq_len(length(children))) {
        child <- children[i]
        coordinates <- xml2::xml_attrs(child)[[1]]
        coordinates["y"] <- as.numeric(coordinates["y"]) + y.margin
        
        if (xml2::xml_name(child) == "start") {
            # state variables has different shape in PD and ER, so we need to switch the
            # shape when needed
            arc.line["start.x"] <- coordinates["x"]
            arc.line["start.y"] <- coordinates["y"]
            
        } else if (xml2::xml_name(child) == "next") {
            arc.line["end.x"] <- coordinates["x"]
            arc.line["end.y"] <- coordinates["y"]
            if (is.na(arc.line["id"])) {
                arc.line["id"] <- paste("next", arc.info["source"], arc.info["target"], 
                                        sep = "->")
            }
            arc <- new("next.sbgn.arc", id = paste(arc.line["id"], arc.line["start.x"], sep = "_"), 
                       start.x = as.numeric(arc.line["start.x"]), start.y = as.numeric(arc.line["start.y"]), 
                       end.x = as.numeric(arc.line["end.x"]), end.y = as.numeric(arc.line["end.y"]))
            #arc@parameters.list <- global.parameters.list
            arc@parameters.list <- list()
            arc@arc.class <- "next"
            arc@edge <- edge.paras
            
            arcs.list[[arc@id]] <- arc
            svg.arc <- paste(svg.arc, plot.arc(arc), sep = "\n")
            arc.line["start.x"] <- coordinates["x"]  # after ploting this 'next', set the new_arc's start coordinates
            arc.line["start.y"] <- coordinates["y"]
            
        } else if (xml2::xml_name(child) == "end") {
            arc.line["end.x"] <- coordinates["x"]
            arc.line["end.y"] <- coordinates["y"]
            if (sum(as.numeric(arc.line[c("start.x", "start.y", "end.x", "end.y")])) == 0) {
                arc.line <- find.arc.coordinates(arc.line, glyphs)
            }
            if (is.na(arc.line["id"])) {
                arc.line["id"] <- paste(arc.info["source"], arc.info["target"], sep = "->")
            }
            arc <- new(paste(arc.class, ".sbgn.arc", sep = ""), id = paste(arc.line["id"], arc.line["start.x"], sep = "_"), 
                       start.x = as.numeric(arc.line["start.x"]), start.y = as.numeric(arc.line["start.y"]), 
                       end.x = as.numeric(arc.line["end.x"]), end.y = as.numeric(arc.line["end.y"]))
            arc@arc.class <- arc.class
            arc@source <- arc.info["source"]
            arc@target <- arc.info["target"]
            #arc@parameters.list <- global.parameters.list
            arc@parameters.list <- list()
            arc@edge <- edge.paras
            arcs.list[[arc@id]] <- arc
            svg.arc <- paste(svg.arc, plot.arc(arc), sep = "\n")
        }
    }
    return(list(arcs.list = arcs.list, svg.arc = svg.arc))
}

#########################################################################################################
find.arc.coordinates <- function(arc.line, glyphs) {
    
    arc.source <- arc.line["source"]
    arc.target <- arc.line["target"]
    
    node.source <- glyphs[[arc.source]]
    source.points <- list()
    
    if (is(node.source, "port")) {
        source.points[["port"]] <- c(x = node.source@x, y = node.source@y)
    } else {
        s.h <- node.source@h
        s.w <- node.source@w
        s.x <- node.source@x - s.w/2
        s.y <- node.source@y - s.h/2
        source.points[["n.s.1"]] <- c(x = s.x + s.w/2, y = s.y)
        source.points[["n.s.2"]] <- c(x = s.x + s.w/2, y = s.y + s.h)
        source.points[["n.s.3"]] <- c(x = s.x, y = s.y + s.h/2)
        source.points[["n.s.4"]] <- c(x = s.x + s.w, y = s.y + s.h/2)
    }
    
    node.target <- glyphs[[arc.target]]
    
    target.points <- list()
    if (is(node.target, "port")) {
        target.points[["port"]] <- c(x = node.target@x, y = node.target@y)
    } else {
        t.h <- node.target@h
        t.w <- node.target@w
        t.x <- node.target@x - t.w/2
        t.y <- node.target@y - t.h/2
        target.points[["n.t.1"]] <- c(x = t.x + t.w/2, y = t.y)
        target.points[["n.t.2"]] <- c(x = t.x + t.w/2, y = t.y + t.h)
        target.points[["n.t.3"]] <- c(x = t.x, y = t.y + t.h/2)
        target.points[["n.t.4"]] <- c(x = t.x + t.w, y = t.y + t.h/2)
    }
    min.dist <- Inf
    min.node.pairs <- data.frame()
    all.pairs <- expand.grid(source.points, target.points)
    for (i in seq_len(nrow(all.pairs))) {
        node.pairs <- all.pairs[i, ]
        s.xy <- node.pairs[1][[1]]
        t.xy <- node.pairs[2][[1]]
        distance <- (s.xy["x"] - t.xy["x"])^2 + (s.xy["y"] - t.xy["y"])^2
        if (distance < min.dist) {
            min.dist <- distance
            min.node.pairs <- node.pairs
        }
    }
    
    arc.line["start.x"] <- min.node.pairs[1][[1]]["x"]
    arc.line["start.y"] <- min.node.pairs[1][[1]]["y"]
    arc.line["end.x"] <- min.node.pairs[2][[1]]["x"]
    arc.line["end.y"] <- min.node.pairs[2][[1]]["y"]
    return(arc.line)
    
}

#########################################################################################################
# Use user input list of glyphs or generate glyph objects. Generate svg code for plotting glyph objects
parse.glyph <- function(sbgn.xml, user.data, if.plot.svg = TRUE, y.margin = 0, max.x, 
                        global.parameters.list, sbgn.id.attr, glyphs.user = list(), 
                        compartment.layer.info, if.plot.cardinality) {
    
    if.plot.annotation.nodes <- global.parameters.list$if.plot.annotation.nodes
    if.use.number.for.long.label <- global.parameters.list$if.use.number.for.long.label
    # parse glyphs and plot glyphs and ports
    map.language <- xml2::xml_attrs(xml2::xml_find_all(sbgn.xml, ".//map")[[1]])["language"]
    if (is.null(map.language)) {
        map.language <- ""
    }
    message("\nMap language is ", map.language, "\n")
    
    # glyph set list
    glyph.names <- c()
    # find plot parameters svg contents
    svg.ports <- ""
    svg.nodes <- ""
    svg.nodes.complex <- ""
    svg.nodes.compartment <- list()
    svg.cardinality <- ""  # the cadinality node are supposed to be in front of the arcs, so need to print it again at the end of the svg file
    node.set.list <- list(molecules.with.state.variables = list())
    if.has.chemical.nodes <- FALSE
    if.has.non.chemical.nodes <- FALSE
    glyph.coors <- list()
    user.defined.glyphs <- names(glyphs.user)
    
    shorter.label.mapping.list <- c("shorter label", "original label")
    long.words.count.list <- list()
    no.id.node.count.list <- list()
    all.glyphs.xml <- xml2::xml_find_all(sbgn.xml, ".//glyph")
    nodes.objs <- list()
    for (i in seq_along(all.glyphs.xml)) {
        glyph <- all.glyphs.xml[[i]]
        glyph.info <- xml2::xml_attrs(glyph)  # extract attributes of this glyph
        # Add compartment to free nodes
        if (is.na(glyph.info["compartmentRef"])) {
            glyph.info["compartmentRef"] <- "free.node"
        }
        
        glyph.class <- gsub(" ", "_", glyph.info["class"])
        if (glyph.class %in% c("simple_chemical_multimer", "simple_chemical")) {
            if.has.chemical.nodes <- TRUE
        }
        if (glyph.class %in% c("macromolecule", "nucleic_acid_feature")) {
            if.has.non.chemical.nodes <- TRUE
        }
        node <- new(paste(glyph.class, ".sbgn", sep = ""))
        node@glyph.class <- glyph.info["class"]
        
        # Some glyphs don't have ID, generate IDs for them
        parsed.ids.record <- generate.glyph.id(glyph.id = glyph.info["id"], glyph.class = glyph.info["class"], 
                                               glyph, node.set.list, no.id.node.count.list)
        node@id <- parsed.ids.record$glyph.id
        node.set.list <- parsed.ids.record$node.set.list
        no.id.node.count.list <- parsed.ids.record$no.id.node.count.list
        
        # use user glyph if found
        if (node@id %in% user.defined.glyphs) {  
            node <- glyphs.user[[node@id]]
            label.margin <- node@label.margin
            svg.port <- node@svg.port
            if (length(node@clone) > 0) {
                node.clone <- node@clone[[1]]
            } else {
                node.clone <- ""
            }
            if (is.null(node)) {
                message("No node found in user defined glyphs!\n")
                print(user.defined.glyphs)
                print(node@id)
                print(node)
            }
        } else {
            # create glyph object
            parse.result <- generate.node.obj(glyph, glyph.class, glyph.info, node, if.plot.svg, y.margin, 
                                              sbgn.id.attr, user.data, max.x, global.parameters.list, 
                                              if.use.number.for.long.label, if.plot.annotation.nodes, 
                                              map.language, long.words.count.list, shorter.label.mapping.list)
            node <- parse.result$node
            long.words.count.list <- parse.result$long.words.count.list
            shorter.label.mapping.list <- parse.result$shorter.label.mapping.list
            node.clone <- node@clone
        }
        glyph.names <- c(glyph.names, node@id)
        if (if.plot.svg) {
            if (glyph.class == "annotation" & !if.plot.annotation.nodes) {
            } else if (is(node, "complex.sbgn")) {
                svg.nodes.complex <- paste(svg.nodes.complex, plot.glyph(node), sep = "\n")
            } else if (is(node, "compartment.sbgn")) {
                # if this node is a clone marker
                svg.nodes.compartment <- c(svg.nodes.compartment, plot.glyph(node))
                names(svg.nodes.compartment)[length(svg.nodes.compartment)] <- node@id
            } else if (!is.character(node.clone)) {
                # if this node is a clone marker
                svg.nodes <- paste(svg.nodes, plot.glyph(node), plot.glyph(node.clone), 
                                   sep = "\n")
            } else if (is(node, "cardinality.sbgn") | is(node, "stoichiometry.sbgn")) {
                svg.cardinality <- paste(svg.cardinality, plot.glyph(node), sep = "\n")
            } else {
                svg.nodes <- paste(svg.nodes, plot.glyph(node), sep = "\n")
            }
            svg.ports <- paste(svg.ports, svg.port, sep = "\n")
        }
        
        ############################################################################### 
        if (glyph.class %in% c("complex", "compartment")) {
            label.margin <- node@label.margin
            glyph.coor <- c(node@x - node@w/2, node@y - node@h/2 - label.margin - 
                                6, node@x + node@w/2, node@y + node@h/2)  # find the node boundaries. We'll use them to find a place to put color panel
        } else {
            glyph.coor <- c(node@x - node@w/2, node@y - node@h/2, node@x + node@w/2, 
                            node@y + node@h/2)  # find the node boundaries. We'll use them to find a place to put color panel
        }
        glyph.coors[[node@id]] <- glyph.coor
        nodes.objs <- c(nodes.objs, node)
    }
    glyph.coors <- do.call(rbind, glyph.coors)
    colnames(glyph.coors) <- c("x", "y", "xw", "yh")
    names(nodes.objs) <- glyph.names
    # parse compartment
    if (length(compartment.layer.info) > 0) {
        if (!identical(compartment.layer.info, "original")) {
            svg.nodes.compartment <- svg.nodes.compartment[compartment.layer.info]
        }
    }
    svg.nodes.compartment <- paste(svg.nodes.compartment, collapse = "\n")
    # check cardinality
    if (!if.plot.cardinality) {
        svg.cardinality <- ""
    }
    out.list <- list(glyphs = nodes.objs, svg.ports = svg.ports, svg.nodes = svg.nodes, 
                     svg.nodes.complex = svg.nodes.complex, svg.nodes.compartment = svg.nodes.compartment, 
                     svg.cardinality = svg.cardinality, if.has.chemical.nodes = if.has.chemical.nodes, 
                     if.has.non.chemical.nodes = if.has.non.chemical.nodes, glyph.coors = glyph.coors, 
                     shorter.label.mapping.list = shorter.label.mapping.list)
    return(out.list)
}

#########################################################################################################
# generate IDs for glyphs for ones with no IDs
generate.glyph.id <- function(glyph.id, glyph.class, glyph, node.set.list, no.id.node.count.list) {
    
    if (is.na(glyph.id)) {
        if (glyph.class %in% c("unit of information", "state variable")) {
            glyph.parent <- xml2::xml_parent(glyph)
            parent.id <- xml2::xml_attr(glyph.parent, "id")
            if (is.null(node.set.list[["molecules.with.state.variables"]][[parent.id]])) {
                node.set.list[["molecules.with.state.variables"]][[parent.id]] <- 1
            } else {
                node.set.list[["molecules.with.state.variables"]][[parent.id]] <- node.set.list[["molecules.with.state.variables"]][[parent.id]] + 1
            }
            v.i <- node.set.list[["molecules.with.state.variables"]][[parent.id]]
            glyph.id <- paste(parent.id, v.i, sep = ".info.")
            
        } else if (!glyph.class %in% names(no.id.node.count.list)) {
            # some nodes don't have id(cardinality), we need to gennerate an id for them
            no.id.node.count.list[[glyph.class]] <- 0
        } else {
            no.id.node.count.list[[glyph.class]] <- no.id.node.count.list[[glyph.class]] + 1
        }
        index.id <- no.id.node.count.list[[glyph.class]]
        glyph.id <- paste(glyph.class, index.id, sep = ":")
    }
    return(list(glyph.id = glyph.id, node.set.list = node.set.list, 
                no.id.node.count.list = no.id.node.count.list))
}

#########################################################################################################
# used in generate.node.obj function from mapping.utilities.R
# generate children glyph objects
parse.glyph.children <- function(map.language, glyph, glyph.class, glyph.info, 
                                 node, if.plot.svg, y.margin) {
    
    if.complex.empty <- TRUE
    node.clone <- ""
    node.label <- ""
    svg.port <- ""
    glyph.port.info <- c()
    children <- xml2::xml_children(glyph)
    for (i in seq_len(length(children))) {
        child <- children[[i]]
        
        if (xml2::xml_name(child) == "state") {
            # state variables has different shape in PD and ER, so we need to switch the
            # shape when needed
            if (map.language == "entity relationship") {
                original.id <- node@id
                node <- new(paste(glyph.class, ".ER.sbgn", sep = ""))
                node@id <- original.id
                node@glyph.class <- glyph.info["class"]
            }
            glyph.label <- xml2::xml_attrs(child)  # find the label information for this glyph
            if (is.na(glyph.label["variable"])) {
                # sometimes there is no variable name, in this case we just show the value
                node@label <- glyph.label["value"]
            } else {
                node@label <- paste(glyph.label["value"], "@", glyph.label["variable"], sep = "")
            }
        } else if (xml2::xml_name(child) == "clone") {
            # if this parent node has a clone marker
            node.clone <- new("clone.sbgn")  # create a node for the marker
            if (glyph.class == "simple_chemical") {
                node.clone <- new("clone_simple_chemical.sbgn")
            }
            clone.children <- xml2::xml_children(child)
            if (length(clone.children) == 0) {
                # if the marker has no label
                node.clone@label <- ""
            } else {
                child.clone.label <- clone.children[1]  # if the marker has a label, find it
                child.clone.label.label <- xml2::xml_attrs(child.clone.label)  # find the label information for this glyph
                node.clone@label <- child.clone.label.label[[1]]
            }
        } else if (xml2::xml_name(child) == "label") {
            glyph.label <- xml2::xml_attrs(child)  # find the label information for this glyph
            node.label <- glyph.label["text"]
            node@label <- glyph.label["text"]
            glyph.info["label"] <- glyph.label
            
        } else if (xml2::xml_name(child) == "entity") {
            glyph.entity <- xml2::xml_attrs(child)  # find the label information for this glyph
            if (map.language == "activity flow" & glyph.class == "unit_of_information") {
                glyph.class <- gsub(" ", "_", glyph.entity)
                original.id <- node@id
                node <- new(paste(glyph.class, ".sbgn", sep = ""))
                node@id <- original.id
                node@glyph.class <- glyph.entity
                node@label <- node.label
                node@label_location <- "center"
            }
        } else if (xml2::xml_name(child) == "bbox") {
            glyph.box.information <- xml2::xml_attrs(child)  # find the box information for this glyph
            glyph.box.information["y"] <- as.numeric(glyph.box.information["y"]) + y.margin
            node@x <- as.numeric(glyph.box.information["x"])
            node@y <- as.numeric(glyph.box.information["y"])
            node@h <- as.numeric(glyph.box.information["h"])
            node@w <- as.numeric(glyph.box.information["w"])
            
            node@x <- node@x + node@w/2
            node@y <- node@y + node@h/2
        } else if (xml2::xml_name(child) == "port") {
            glyph.port.info <- xml2::xml_attrs(child)  # find the box information for this glyph
            glyph.port.info["y"] <- as.numeric(glyph.port.info["y"]) + y.margin
            # if port has coordinates, plot it
            if (as.numeric(glyph.port.info["x"]) + as.numeric(glyph.port.info["y"]) !=  0) {
                svg.port <- paste(svg.port, plot.arc.ports(glyph.port.info, node), sep = "\n")
            }
            node@svg.port <- svg.port
        } else if (xml2::xml_name(child) == "glyph" & 
                   xml2::xml_attr(child, "class") != "annotation") {
            if.complex.empty <- FALSE
        }
    }
    
    return(list(node = node, node.clone = node.clone, glyph.port.info = glyph.port.info, 
                svg.port = svg.port, node.label = node.label, glyph.info = glyph.info, 
                if.complex.empty = if.complex.empty))
}

#########################################################################################################
# plot arcs for ports of a glyph
plot.arc.ports <- function(glyph.port.info, node) {
    
    port.x <- as.numeric(glyph.port.info["x"])
    port.y <- as.numeric(glyph.port.info["y"])
    node.x <- node@x
    node.y <- node@y
    node.h <- node@h
    node.w <- node@w
    if (abs(port.y - node.y) < node.h/2 & port.x < node.x) {
        # sometimes there is small difference between the coordinates of port and it's
        # glyph. So we need to use abs(port.y- node.y)<node.h/2.
        x1 <- node.x - node.w/2
        y1 <- node.y
        x2 <- port.x
        y2 <- port.y
    } else if (abs(port.y - node.y) < node.h/2 & port.x > node.x) {
        x1 <- node.x + node.w/2
        y1 <- node.y
        x2 <- port.x
        y2 <- port.y
    } else if (port.y < node.y & abs(port.x - node.x) < node.w/2) {
        x1 <- node.x
        y1 <- node.y - node.h/2
        x2 <- port.x
        y2 <- port.y
    } else if (port.y > node.y & abs(port.x - node.x) < node.w/2) {
        x1 <- node.x
        y1 <- node.y + node.h/2
        x2 <- port.x
        y2 <- port.y
    } else {
        x1 <- port.x
        y1 <- port.y
        x2 <- node.x
        y2 <- node.y
    }
    LineCoordinates <- paste("x1=", x1, " y1=", y1, " x2=", x2, " y2=", y2, "", sep = "\"")
    svg.port <- plot.line(LineCoordinates, glyph.port.info["id"], stroke.color = "black")
    return(svg.port)
}

#########################################################################################################
#' Extract glyph information
#' 
#' This function will extract node information such as complex members, compartment members, node class, nodes with state variables etc. 
#' 
#' @param input.sbgn A character string. required. The pathway ID of pre-collected pathways or a path to a local SBGN-ML file. 
#' @param database A character string. If the SBGN-ML file is from one of these databases: 'MetaCyc' and 'pathwayCommons', this parameter should be set to the corresponding string. For these two databases, this function can output other ID types instead of the original IDs in the SBGN-ML files. Otherwise, the output IDs are the oritinal IDs in the 'id' attribure in the 'glyph' element.
#' @param output.gene.id.type A character string. The desired output gene ID type. It only works when the SBGN-ML file is from one of these databases: 'MetaCyc' and 'pathwayCommons'. Currently, only 'macromolecule' glyphs are supported. Please check \code{data('mapped.ids')} for the types accepted. If omitted, the output IDs are the oritinal IDs in the SBGN-ML file.
#' @param show.ids.for.multiHit Vector. When there are multiple output IDs mapped to a glyph, we only show the ones in this vector.
#' @param output.cpd.id.type A character string.  The desired output compound ID type. It only works when the SBGN-ML file is from one of these databases: 'MetaCyc' and 'pathwayCommons'. Currently, only 'simple chemical' glyphs are supported. Please check \code{data('mapped.ids')} for the types accepted. If omitted, the output IDs are the oritinal IDs in the SBGN-ML file.
#' @param species A character string. Only output IDs from this particular species. It only works when the SBGN-ML file is from one of these databases: 'MetaCyc' and 'pathwayCommons'. Please check data('supported.species') for supported species. If omitted, the function will output IDs from all species.
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files.
#' @param sbgn.dir A character string. Default: "./". The path to a folder that will hold created SBGN-ML files, if the input.sbgn are IDs of pre-collected pathways.
#' @return A list containing extracted glyph information.
#' @details  The following glyph information is extracted: complex members, compartment members,submap members, node class, nodes with state variables, class of state variables, edges with cardinality, nodes with ports, 'source and sink' nodes, process nodes.\cr When trying to output other ID types, sometimes multiple output IDs are mapped to one glyph. In this situation, the IDs are concatenated by '; ' to represent the glyph.
#' @examples 
#' data(mapped.ids)
#' data(sbgn.xmls)
#' data(pathways.info)
#' node.list <- sbgnNodes(input.sbgn = 'P00001',
#'                        output.gene.id.type = 'ENTREZID',
#'                        output.cpd.id.type = 'compound.name',
#'                        species = 'hsa')
#' @export

sbgnNodes <- function(input.sbgn, output.gene.id.type = NA, output.cpd.id.type = NA,
                      database = NA, species = NA, show.ids.for.multiHit = NULL,
                      SBGNview.data.folder = "./SBGNview.tmp.data", sbgn.dir = "./") {
    
    if (!file.exists(SBGNview.data.folder)) {
        dir.create(SBGNview.data.folder)
    }
    if(!is.na(output.cpd.id.type)) {
        if(output.cpd.id.type == "CompoundName") output.cpd.id.type = "compound.name"   
    }
    
    result.list <- list()
    input.sbgns <- unique(as.vector(as.character(input.sbgn)))
    all.pairs.id.mapping.list <- list()
    for (i in seq_len(length.out = length(input.sbgns))) {
        input.sbgn <- input.sbgns[i]
        input.sbgn <- gsub("\"", "", input.sbgn)
        input.sbgn.original <- input.sbgn
        if (is.null(input.sbgn)) {
            stop("Must provide either 'input.sbgn' !!!")
        } else {
            input.sbgn.full.path <- paste(sbgn.dir, input.sbgn, sep = "/")
            if (!file.exists(input.sbgn.full.path[i])) {
                message("\n", input.sbgn, " does not look like an existing local file.\n Using it as pathway ID in 'pathways.info'\n\n ")
                database <- pathways.info[pathways.info[, "pathway.id"] == input.sbgn, "database"]
                input.sbgn <- downloadSbgnFile(pathway.id = input.sbgn, download.folder = sbgn.dir)
                message("SBGN-ML files downloaded to: ", input.sbgn)
            } else {
                message("\n", input.sbgn, "looks like an existing local file. Using it directly.\n")
                if (input.sbgn %in% pathways.info[, "file.name"]) {
                  database <- pathways.info[pathways.info[, "file.name"] == input.sbgn, "database"]
                }
                input.sbgn <- input.sbgn.full.path
            }
        }
        
        if.other.id.types.available <- database %in% c("MetaCrop", "MetaCyc", "pathwayCommons")
        if (if.other.id.types.available) {
            all.id.mapping.result <- load.all.ids.mapping(database, all.pairs.id.mapping.list, 
                                                          species, output.gene.id.type, output.cpd.id.type, 
                                                          SBGNview.data.folder = SBGNview.data.folder)
            id.mapping.all.list <- all.id.mapping.result$id.mapping.all.list
            all.pairs.id.mapping.list <- all.id.mapping.result$all.pairs.id.mapping.list
            output.cpd.id.type.use <- all.id.mapping.result$output.cpd.id.type.use
            output.gene.id.type.use <- all.id.mapping.result$output.gene.id.type.use
            SBGN.file.cpd.id.type <- all.id.mapping.result$SBGN.file.cpd.id.type
            SBGN.file.gene.id.type <- all.id.mapping.result$SBGN.file.gene.id.type
        }
        
        message("reading SBGN-ML file for node set: ", input.sbgn)
        sbgn <- xml2::read_xml(input.sbgn)
        xml2::xml_attrs(sbgn) <- NULL  # Remove root node attribute. This is necessary Otherwise xml2 won't find the nodes when using xml_find_all.
        # node.set.list <- get.all.nodes.info(sbgn, if.other.id.types.available, output.cpd.id.type.use, 
        #                                     output.gene.id.type.use, SBGN.file.cpd.id.type, SBGN.file.gene.id.type, 
        #                                     id.mapping.all.list, show.ids.for.multiHit)
        # result.list[[input.sbgn.original]] <- node.set.list[["all.nodes"]]
        result.list[[input.sbgn.original]] <- get.all.nodes.info(sbgn, if.other.id.types.available, output.cpd.id.type.use, 
                                                                 output.gene.id.type.use, SBGN.file.cpd.id.type, SBGN.file.gene.id.type, 
                                                                 id.mapping.all.list, show.ids.for.multiHit) # returns node.set.matrix
    }
    return(result.list)
}

#########################################################################################################
# get node information (id, class, compartmentRef, etc...) from sbgn file
get.all.nodes.info <- function(sbgn, if.other.id.types.available, output.cpd.id.type.use, 
                               output.gene.id.type.use, SBGN.file.cpd.id.type, SBGN.file.gene.id.type, 
                               id.mapping.all.list, show.ids.for.multiHit) {
    
    #node.set.list <- list(all.nodes = matrix(ncol = 8, nrow = 0))
    node.set.matrix <- matrix(ncol = 8, nrow = 0)
    all.glyphs <- xml2::xml_find_all(sbgn, ".//glyph")
    for (i in seq_len(length(all.glyphs))) {
        glyph <- all.glyphs[[i]]
        nodes.info <- c(glyph.id = "", id = "", class = "", parent.complex = "", 
                        parent.submap = "", parent.compartment = "", parent.macromolecule = "", 
                        parent.node = "")
        nodes.info
        glyph.class <- xml2::xml_attr(glyph, "class")
        glyph.id <- xml2::xml_attr(glyph, "id")
        glyph.id.original <- glyph.id
        
        if (if.other.id.types.available & glyph.class %in% c("macromolecule", "simple chemical")) {
            # print(glyph.id)
            if (xml2::xml_attr(xml2::xml_parent(glyph), "class") %in% c("complex", 
                                                                        "submap") & !is.na(glyph.id)) {
                parent.id <- xml2::xml_attr(xml2::xml_parent(glyph), "id")
                glyph.id <- gsub(paste("_", parent.id, sep = ""), "", glyph.id)  # get rid of the complex name. That's how pathwayCommons name molecules within complex
            }
            glyph.id <- change.glyph.id(glyph.id, glyph.class, id.mapping.all.list, 
                                        output.gene.id.type.use, output.cpd.id.type.use, SBGN.file.cpd.id.type, 
                                        SBGN.file.gene.id.type, show.ids.for.multiHit = show.ids.for.multiHit)
        }
        nodes.info["id"] <- glyph.id
        nodes.info["class"] <- glyph.class
        
        glyph.compartment <- xml2::xml_attr(glyph, "compartmentRef")
        if (!is.na(glyph.compartment) & !is.na(glyph.id) & glyph.class != "unit of information" & 
            glyph.class != "state variable" & glyph.class != "source and sink") {
            # 'source and sink' will be assigned to their neighbor's compartment
            nodes.info["parent.compartment"] <- glyph.compartment
        }
        # generate nodes for layout, for the nodes not in any compartment. exclude
        # state.variable, nodes with no id, , still needs complex, because the network is
        # generated from edge list, so that step determins if a node is included in the
        # network, not this step
        if (is.na(glyph.compartment) & !is.na(glyph.id) & glyph.class != "unit of information" & 
            glyph.class != "state variable") {
            nodes.info["parent.compartment"] <- "free.node.compartment"
        }
        # test if the parent of this node is complex if the node's parent have a class
        # attribute
        if (!is.na(xml2::xml_attr(xml2::xml_parent(glyph), "class")) & !glyph.class %in% 
            c("annotation", "unit of information", "state variable")) {
            if (xml2::xml_attr(xml2::xml_parent(glyph), "class") %in% c("complex", 
                                                                        "submap") & !is.na(glyph.id)) {
                parent.id <- xml2::xml_attr(xml2::xml_parent(glyph), "id")
                if (xml2::xml_attr(xml2::xml_parent(glyph), "class") == "submap") {
                    nodes.info["parent.submap"] <- parent.id
                } else {
                    nodes.info["parent.complex"] <- parent.id
                }
            }
        }
        # generate state variable to molecule mapping
        if ((glyph.class == "unit of information" | glyph.class == "state variable")) {
            parent.id <- xml2::xml_attr(xml2::xml_parent(glyph), "id")
            nodes.info["parent.macromolecule"] <- parent.id
        }
        nodes.info[["glyph.id"]] <- glyph.id.original
        #node.set.list[["all.nodes"]] <- rbind(node.set.list[["all.nodes"]], nodes.info)
        node.set.matrix <- rbind(node.set.matrix, nodes.info)
    }
    #return(node.set.list)
    return(node.set.matrix)
}

#########################################################################################################
# generate data using the glyphs ids from sbgn file
simulate.user.data <- function(sbgn.file, n.samples = 3, ...) {
    
    all.nodes <- node.ids.from.sbgn(sbgn.file, ...)
    if (is.null(all.nodes)) {
        return("no.nodes")
    }
    user.data <- runif(length(all.nodes), -11, 11)
    user.data <- as.data.frame(matrix(runif(length(all.nodes) * n.samples, -1, 1), 
                                      ncol = n.samples))
    row.names(user.data) <- all.nodes
    user.data <- generate.user.data(user.data)
    return(user.data)
}

#########################################################################################################
# find unique glyphs in sbgn file
node.ids.from.sbgn <- function(sbgn.file, output.glyph.class = c("macromolecule", "simple chemical"), 
                               if.include.complex.member = FALSE) {
    
    message("reading SBGN-ML file for node ids: ", sbgn.file)
    sbgn <- xml2::read_xml(sbgn.file)
    xml2::xml_attrs(sbgn) <- NULL  # Remove root node attribute. This is necessary Otherwise xml2 won't find the nodes when using xml_find_all.
    # glyph.info = do.call(rbind,xml_attr(xml_find_all(sbgn,'.//glyph'),'id'))
    glyph.info <- xml2::xml_attr(xml2::xml_find_all(sbgn, ".//glyph"), "id")
    
    glyph.class <- xml2::xml_attr(xml2::xml_find_all(sbgn, ".//glyph"), "class")
    if (!identical(output.glyph.class, "all")) {
        glyph.info <- glyph.info[glyph.class %in% output.glyph.class]
    }
    if (!if.include.complex.member) {
        # this function need to be improved. Now we assume complex members have 'Complex'
        # in their IDs, which may only apply to pathwayCommons naming.
        glyph.info <- glyph.info[!grepl("Complex", glyph.info)]
    }
    glyph.info <- unique(glyph.info)
    glyph.info <- glyph.info[!is.na(glyph.info)]
    return(glyph.info)
}

#########################################################################################################
generate.user.data <- function(user.data) {
    
    for (j in seq_len(ncol(user.data))) {
        # if the variance is zero, scale will produce NaN and cause error
        print("simulate user data")
        print(j)
        print(head(user.data))
        print(head(user.data[, j]))
        vari <- var(user.data[, j])
        if (vari == 0) {
            user.data[1, j] <- mean(user.data[, j] + j)
        }
    }
    user.data <- scale(user.data)
    user.data <- apply(user.data, 2, function(x) {
        mean.x <- mean(x)
        max.x.user.data <- max(x)
        min.x <- min(x)
        vapply(x, function(i) {
            if (i < mean.x) {
                normalized <- (i - mean.x)/abs(mean.x - min.x)
                return(normalized)
            } else {
                normalized <- (i - mean.x)/abs(mean.x - max.x.user.data)
                return(normalized)
            }
        }, numeric(1) )
    })
    return(user.data)
}

#########################################################################################################
# find ports in sbgn file and return objects of class port
# used by renderSbgn function
xml.to.port.glyphs <- function(sbgn.xml, y.margin = 0) {
    
    ports <- list()
    ports.info <- do.call(rbind, xml2::xml_attrs(xml2::xml_find_all(sbgn.xml, ".//port")))
    if (!is.null(ports.info)) {
        ports <- apply(ports.info, 1, function(port) {
            port.object <- new("port", x = as.numeric(port["x"]), y = as.numeric(port["y"]) + 
                                   y.margin, id = port["id"])
            return(port.object)
        })
        names(ports) <- ports.info[, "id"]
    }
    return(ports)
}

#########################################################################################################
# get arc infromation from sbgn file
# used in SBGNview function
get.arcs.info <- function(sbgn) {
    
    # Retrieve edge information
    edge.spline.info.node <- xml2::xml_find_all(sbgn, ".//edge.spline.info")
    # in the original design, the routed edges are in the format of svg. In updated
    # version, the edges are recorded by their information in element
    # 'edge.spline.info'(location etc.)
    if (length(edge.spline.info.node) > 0) {
        arcs.info <- "parse splines"
    } else {
        arcs.info <- "straight"
    }
    return(arcs.info)
}

#########################################################################################################
# find compartments to check if any overlap
# used in SBGNview function
get.compartment.layer <- function(sbgn) {
    
    compartment.layer.info <- xml2::xml_find_all(sbgn, ".//compartment.layer.info")
    if (length(compartment.layer.info) > 0) {
        compartment.layer.info <- strsplit(xml2::xml_text(compartment.layer.info[1]), 
                                           "-compartment.layer.info-")[[1]]
    } else {
        compartment.layer.info <- "original"
    }
    return(compartment.layer.info)
}

#########################################################################################################
# used in renderSbgn function
# new version of function. removed max.xw when assigning y.margin
# incrementally increase/adjust y.margin (top margin) based on other parameters that affect overall svg height
find.max.xy <- function(sbgn.xml, arcs.info, color.panel.scale, global.parameters.list) {
    
    box.info <- do.call(rbind, xml2::xml_attrs(xml2::xml_find_all(sbgn.xml, ".//bbox")))
    x <- as.numeric(box.info[, "x"])
    w <- as.numeric(box.info[, "w"])
    max.xw <- max(x + w)
    
    y <- as.numeric(box.info[, "y"])
    h <- as.numeric(box.info[, "h"])
    max.yh <- max(y + h)
    
    min.y <- min(y)
    min.x <- min(x)
    
    y.margin <- max(70, max.yh/20 * color.panel.scale)
    
    add.scale.factor <- 14
    if(global.parameters.list$if.scale.compartment.font.size) {
        y.margin <- y.margin + (add.scale.factor * 100*global.parameters.list$node.width.adjust.factor.compartment)
    }
    if (global.parameters.list$if.scale.complex.font.size) {
        y.margin <- y.margin + (add.scale.factor * 100*global.parameters.list$node.width.adjust.factor.complex)
    }
    if (global.parameters.list$font.size.scale.compartment != 1.6) { # 1.6 default for font.size.scale.compartment
        y.margin <- y.margin + (add.scale.factor*global.parameters.list$font.size.scale.compartment)
    }
    if (global.parameters.list$font.size.scale.complex != 1.1) { # 1.1 default for font.size.scale.complex
        y.margin <- y.margin + (add.scale.factor*global.parameters.list$font.size.scale.complex)
    }
    if (global.parameters.list$font.size != 3) { # 3 default for font.size
        y.margin <- y.margin + (4*global.parameters.list$font.size)
    }
    if (global.parameters.list$pathway.name.font.size != 1) { # 1 default for pathway.name.font.size
        y.margin <- y.margin + (add.scale.factor*global.parameters.list$pathway.name.font.size)
    }
    
    return(list(max.xw = max.xw, max.yh = max.yh, min.y = min.y, min.x = min.x, y.margin = y.margin))
}

#########################################################################################################
## old version. assigns too much margin on the top
# find.max.xy <- function(sbgn.xml, arcs.info, color.panel.scale) {
#     box.info <- do.call(rbind, xml2::xml_attrs(xml2::xml_find_all(sbgn.xml, ".//bbox")))
#     x <- as.numeric(box.info[, "x"])
#     w <- as.numeric(box.info[, "w"])
#     max.xw <- max(x + w)
#     
#     y <- as.numeric(box.info[, "y"])
#     h <- as.numeric(box.info[, "h"])
#     max.yh <- max(y + h)
#     
#     min.y <- min(y)
#     min.x <- min(x)
#     
#     if (arcs.info == "straight") {
#         # if there is no spline edges, calculate margin to move both nodes and edges
#         # coordinates. This will give room for color legend
#         y.margin <- max(100, max(max.yh, max.xw)/22) * max(1, color.panel.scale)
#     } else {
#         # if there is routed edges' svg in the SBGN-ML file, we can't move the nodes
#         y.margin <- max(100, max(max.xw, max.yh)/15 * color.panel.scale)
#     }
#     return(list(max.xw = max.xw, max.yh = max.yh, min.y = min.y, min.x = min.x, y.margin = y.margin))
# }

#########################################################################################################
