
#########################################################################################################
# needs improvement. Currently, input in a set of nodes, 
# it will highlight all edges among the input nodes
highlight.edges <- function(node.set, arcs = NULL, node.set.id.type = NULL, arc.node.id.type = NULL, 
                            mol.type = NULL, arc.highlight.stroke.color = "red", 
                            arc.highlight.stroke.width = 4, arc.highlight.tip.size = 4) {
    
    if (node.set.id.type != arc.node.id.type) {
        new.ids <- changeIds(input.ids = node.set, input.type = node.set.id.type, 
            output.type = node.set.id.type)
        unique.n.ids <- c()
        for (i in seq_len(length.out = length(new.ids))) {
            n.ids <- new.ids[[i]]
            if (length(n.ids) == 0) {
                warning("\n", names(new.ids[i]), " can't be mapped to ", node.set.id.type, 
                  "!!!")
            } else if (length(n.ids) > 1) {
                warning("\n", names(new.ids[i]), " is mapped to multiple ", node.set.id.type, 
                  " IDs!!!\nUsing the first one!!")
                print(paste(n.ids, collapse = "\t"))
                unique.n.ids <- c(unique.n.ids, n.ids[1])
            }
        }
        new.ids <- unique.n.ids
        new.ids
    } else {
        new.ids <- node.set
    }
    
    for (i in seq_along(arcs)) {
        arc <- arcs[[i]]
        source.id <- try(arc@source)
        target.id <- arc@target
        if (all(c(source.id, target.id) %in% new.ids)) {
            arcs[[i]]@edge$line.width <- arc.highlight.stroke.width
            arcs[[i]]@edge$line.stroke.color <- arc.highlight.stroke.color
            arcs[[i]]@edge$tip.stroke.color <- arc.highlight.stroke.color
            arcs[[i]]@edge$tip.fill.color <- arc.highlight.stroke.color
            arcs[[i]]@parameters.list$edge.tip.size <- arc.highlight.tip.size
        }
    }
    return(arcs)
}

#########################################################################################################
#' Highlight input nodes
#' 
#' Change node properties such as border color and width to highlight a list of input nodes. This function should be used as the second argument to function. Run \code{help("+.SBGNview")} for more information. 
#' 
#' @param node.set A vector of character strings. Default: "all". Input molecule IDs whose nodes are to be highlighted. It can be any ID types supported by SBGNview. 
#' @param node.set.id.type A character string. Default: "compound.name". ID type of input nodes.
#' @param glyphs.id.type A character string. Default: "pathwayCommons". ID type of nodes in SBGN-ML file.
#' @param mol.type A character string. One of 'gene' or 'cpd' (default). 
#' @param stroke.color A character string. Default: "black". The border color of highlighted nodes.
#' @param stroke.width Integer. Default: 10. The border stroke width of highlighted nodes.
#' @param stroke.opacity Numeric. Default: 1. The border stroke opacity of highlighted nodes.
#' @param show.glyph.id  Logical. Default: F. If set to 'TRUE', glyph IDs are plotted instead of their labels.
#' @param select.glyph.class Character vector. Select glyphs by class. It should be one of the values of XML attribute 'class' of a 'glyph' element. For example 'macromolecule', "simple chemical"
#' @param label.x.shift Integer. Default: 0. Move labels horizontally.
#' @param label.y.shift  Integer. Default: 0. Move labels vertically. 
#' @param label.color  A character string. Default: "black". Change label color.
#' @param label.font.size  Integer. Adjust label font size. 
#' @param label.spliting.string  A character vector. When we split text into multiple lines, these characters will be used to split label (where a new line can be created). For example: label.spliting.string = c(' ','-',';','/','_'). 
#' @param labels  A vector of character strings. The labels to display on the nodes of input molecules. The names of this vector should be the vector of 'node.set'. The values are the new labels to display.
#' @return A SBGNview obj
#' @examples 
#' \dontrun{
#' obj.new <- SBGNview.obj + highlightNodes(node.set = c('tyrosine', '(+-)-epinephrine'),
#'                                          stroke.width = 4, 
#'                                          stroke.color = 'green') 
#' }
#' @export
 
highlightNodes <- function(node.set = "all", node.set.id.type = "compound.name", glyphs.id.type = "pathwayCommons", 
                           mol.type = "cpd", stroke.color = "black", stroke.width = 10, stroke.opacity = 1, 
                           show.glyph.id = FALSE, select.glyph.class = c(), label.x.shift = 0, label.y.shift = 0, 
                           label.color = "black", label.font.size = NULL, 
                           label.spliting.string = NULL, labels = NULL) {
    
    # changed mapping file names from CompondName to compound.name and kegg.ligand to kegg
    if(node.set.id.type == "CompoundName") node.set.id.type = "compound.name"
    if(node.set.id.type %in% c("kegg.ligand", "KEGG")) node.set.id.type <- "kegg"
    
    function(SBGNview.obj) {
        glyphs.arcs.list <- SBGNview.obj$data
        sbgns <- names(glyphs.arcs.list)
        for (s in seq_len(length.out = length(glyphs.arcs.list))) {
            # for each sbgn file
            glyphs <- glyphs.arcs.list[[s]]$glyphs.list
            glyphs <- highlight.nodes.each.sbgn(node.set = node.set, select.glyph.class = select.glyph.class, 
                glyphs = glyphs, pathway.id = sbgns[s], node.set.id.type = node.set.id.type, 
                glyphs.id.type = glyphs.id.type, mol.type = mol.type, stroke.width = stroke.width, 
                stroke.color = stroke.color, stroke.opacity = stroke.opacity, show.glyph.id = show.glyph.id, 
                label.x.shift = label.x.shift, label.y.shift = label.y.shift, label.color = label.color, 
                label.font.size = label.font.size, label.spliting.string = label.spliting.string, 
                labels = labels)
            glyphs.arcs.list[[s]]$glyphs.list <- glyphs
        }
        SBGNview.obj$data <- glyphs.arcs.list
        return(SBGNview.obj)
    }
}

#########################################################################################################
highlight.nodes.each.sbgn <- function(node.set = "all", select.glyph.class = c(), glyphs = NULL, pathway.id = NULL,
                                      node.set.id.type = NULL, glyphs.id.type = NULL, mol.type = NULL,
                                      stroke.width = 4, stroke.color = "red", stroke.opacity = 1, show.glyph.id = FALSE,
                                      label.x.shift = 0, label.y.shift = 0, label.color = "black", label.font.size = NULL,
                                      label.spliting.string = NULL, labels = NULL) {
    
    if (!identical(node.set, "all")) {
        # if we are not interested in all nodes, we need to get the node IDs for input
        # IDs if 1. input id type is not glyph id type, we need to change ids
        if (node.set.id.type != glyphs.id.type) {
            input.to.node.ids <- changeIds(input.ids = node.set, input.type = node.set.id.type, 
                output.type = glyphs.id.type, mol.type = mol.type)
            input.to.node.ids
            if (!is.null(labels)) {
                input.ids <- names(input.to.node.ids)
                input.ids
                input.node.id.to.user.labels <- matrix(ncol = 2, nrow = 0)
                for (i in seq_len(length.out = length(input.to.node.ids))) {
                  node.ids <- input.to.node.ids[[i]]
                  node.ids
                  input.id <- input.ids[i]
                  input.id
                  user.label <- labels[input.id]
                  user.label
                  input.node.id.to.user.labels <- rbind(input.node.id.to.user.labels, 
                    cbind(node.ids, user.label))
                }
                input.node.id.to.user.labels
                row.names(input.node.id.to.user.labels) <- input.node.id.to.user.labels[, 
                  "node.ids"]
            }
            input.node.ids <- unlist(input.to.node.ids)
            input.node.ids
        } else {
            # if input ids are node ids, we don't need to change them
            input.node.ids <- node.set
        }
    } else {
        # if we are interested in all nodes
        input.node.ids <- "all"
    }
    
    for (i in seq_len(length.out = length(glyphs))) {
        # i=49
        glyph <- glyphs[[i]]
        glyph.id <- glyph@id
        if (length(glyph.id) == 0) {
            glyph.id <- paste("no.id;x=", glyph@x, ",y=", glyph@y, sep = "")
        }
        glyph.class <- glyph@glyph.class
        if (length(glyph.class) == 0) {
            # some nodes don't have class
            glyph.class <- "no.class"
        }
        glyph.id <- gsub("_Complex.+", "", glyph.id, perl = TRUE)  # id-to-pathwayCommons mapping file was generated from pathwayCommons Biopax file, the ID is 'SmallMolecule_80910e7b27a1c71cb89e71508781c018'. But in the SBGN-ML file (P00001.sbgn), the ID is this 'SmallMolecule_80910e7b27a1c71cb89e71508781c018_Complex_33b98f178303ceef166819f17daaed43'. Complex ID is added, so we need to 
        if (glyph.id %in% input.node.ids | glyph.class %in% select.glyph.class | 
            (length(select.glyph.class) == 0 & "all" %in% node.set)) {
            # if this is the node we are interested in, or, we are interested in all nodes
            if (glyphs[[i]]@label == "SBGNhub Pathway Collection") {
            } else {
                glyphs[[i]]@shape$stroke.width <- stroke.width
                glyphs[[i]]@shape$stroke.color <- stroke.color
                glyphs[[i]]@shape$stroke.opacity <- stroke.opacity
                glyphs[[i]]@text$x.shift <- label.x.shift
                glyphs[[i]]@text$y.shift <- label.y.shift
                if (!is.null(label.font.size)) {
                  glyphs[[i]]@text$font.size <- label.font.size
                }
                glyphs[[i]]@text$color <- label.color
                if (!is.null(label.spliting.string)) {
                  glyphs[[i]]@parameters.list$label.spliting.string <- label.spliting.string
                }
                if (!is.null(labels)) {
                  glyphs[[i]]@label <- input.node.id.to.user.labels[glyph.id, "user.label"]
                }
                if (show.glyph.id) {
                  if (glyphs[[i]]@glyph.class %in% c("process", "uncertain process", 
                    "omitted process")) {
                    glyphs[[i]]@if.show.label <- TRUE
                  }
                  glyphs[[i]]@label <- glyph.id
                }
            }
        }
    }
    return(glyphs)
}

#########################################################################################################
#' Highlight shortest path between two nodes
#' 
#' Given two nodes, find the shortest path between them and highlight the path. Other molecules/nodes and edges involved in reactions in the path are also highlighted. This function should be used as the second argument to function. Run \code{help("+.SBGNview")} for more information.
#' 
#' @param from.node A character string. The molecule ID of source node.
#' @param to.node A character string. The molecule ID of target node.
#' @param directed.graph Logical. If treat the SBGN map as a directed graph. If treat it as directed graph, the direction of an arc is determined by the <arc> XML element in the input SBGN-ML file: from 'source' node to 'target' node.
#' @param node.set.id.type A character string. Default: "compound.name". ID type of input nodes.
#' @param glyphs.id.type A character string. Default: "pathwayCommons". ID type of nodes in SBGN-ML file.
#' @param mol.type A character string. One of 'gene' or 'cpd' (default). 
#' @param path.node.color A character string. Default: "black". Border color of nodes in the path.
#' @param path.node.stroke.width Integer. Default: 10. Adjust stroke width of path between nodes.
#' @param input.node.stroke.width Integer. Default: 10. Border stroke width of input nodes (affects both from.node and to.node)
#' @param from.node.color A character string. Default: "red". Border color of 'source'/from.node
#' @param to.node.color A character string. Default: "blue". Border color of 'target'/to.node
#' @param n.shortest.paths.plot Integer. Default: 1. There could be more than one shortest paths between two nodes. User can choose how many of them to plot.
#' @param shortest.paths.cols A vector of character string. Default: c("purple"). The colors of arcs in different shortest paths. The length of this vector (number of colors) should be the value of n.shortest.paths.plot. If one arc is in multiple shortest paths, its color will be set to the color that appears later in this vector.
#' @param path.stroke.width Integer. Default: 2. The width of line in edges in the shortest paths.
#' @param tip.size Integer. Default: 4. The size of arc tip in edges of the shortest paths.
#' @return A SBGNview obj
#' @examples 
#' \dontrun{
#' obj.new <- SBGNview.obj + highlightPath(from.node = c('tyrosine'), 
#'                                         to.node = c('dopamine'),
#'                                         from.node.color = 'red',
#'                                         to.node.color = 'blue')
#'}
#' @export

highlightPath <- function(from.node = NULL, to.node = NULL, directed.graph = TRUE, 
                          node.set.id.type = "compound.name", glyphs.id.type = "pathwayCommons", 
                          mol.type = "cpd", input.node.stroke.width = 10, from.node.color = "red", 
                          to.node.color = "blue", path.node.color = "black", path.node.stroke.width = 10, 
                          n.shortest.paths.plot = 1, shortest.paths.cols = c("purple"), 
                          path.stroke.width = 2, tip.size = 4) {
    
    # changed mapping file names from CompondName to compound.name and kegg.ligand to kegg
    if(node.set.id.type == "CompoundName") node.set.id.type = "compound.name"
    if(node.set.id.type %in% c("kegg.ligand", "KEGG")) node.set.id.type <- "kegg"
    
    function(SBGNview.obj) {
        glyphs.arcs.list <- SBGNview.obj$data
        sbgns <- names(glyphs.arcs.list)
        for (s in seq_len(length.out = length(glyphs.arcs.list))) {
            # for each sbgn file
            arcs <- glyphs.arcs.list[[s]]$arcs.list
            glyphs <- glyphs.arcs.list[[s]]$glyphs.list
            result.list <- highlight.path.each.sbgn(from.node = from.node, to.node = to.node, 
                glyphs = glyphs, arcs = arcs, pathway.id = sbgns[s], directed.graph = directed.graph, 
                node.set.id.type = node.set.id.type, glyphs.id.type = glyphs.id.type, 
                mol.type = mol.type, from.node.color = from.node.color, to.node.color = to.node.color, 
                input.node.stroke.width = input.node.stroke.width, path.node.color = path.node.color, 
                path.node.stroke.width = path.node.stroke.width, n.shortest.paths.plot = n.shortest.paths.plot, 
                path.stroke.width = path.stroke.width, shortest.paths.cols = shortest.paths.cols, 
                tip.size = tip.size)
            
            glyphs.arcs.list[[s]]$arcs.list <- result.list$arcs
            glyphs.arcs.list[[s]]$glyphs.list <- result.list$glyphs
        }
        SBGNview.obj$data <- glyphs.arcs.list
        return(SBGNview.obj)
    }
}

#########################################################################################################
# generate a new node set containing shortest paths connecting all nodes in the
# input nodeset
# 'arcs.to.graph' function removed and code from that function is added to this function.
shortest.path.from.arcs <- function(arcs, from.node, to.node, directed.graph = TRUE) {

    # code from 'arcs.to.graph(arcs, if.directed = directed.graph)' function which was removed
    edge.df <- lapply(arcs, function(arc) { return(c(source = arc@source, target = arc@target)) })
    edge.df <- as.data.frame(do.call(rbind, edge.df), stringsAsFactors = FALSE)
    grf <- igraph::graph_from_data_frame(edge.df, directed = directed.graph)
    
    from.node.ids <- which(igraph::V(grf)$name %in% from.node)
    from.node.ids <- unique(from.node.ids)
    from.node.ids
    to.node.ids <- which(igraph::V(grf)$name %in% to.node)
    to.node.ids <- unique(to.node.ids)
    to.node.ids
    node.names <- igraph::V(grf)$name
    
    node.paths <- list()
    for (i in seq_len(length.out = length(from.node.ids))) {
        from.node.name <- node.names[from.node.ids[i]]
        from.node.name
        for (j in seq_len(length.out = length(to.node.ids))) {
            to.node.name <- node.names[to.node.ids[j]]
            to.node.name
            pair.name <- paste(from.node.name, to.node.name, sep = "__to__")
            node.paths[[pair.name]] <- list()
            
            all.paths <- igraph::all_shortest_paths(grf, from = from.node.ids[i], 
                                                    to = to.node.ids[j])
            node.paths[[pair.name]][["from.node"]] <- from.node.name
            node.paths[[pair.name]][["to.node"]] <- to.node.name
            node.paths[[pair.name]][["path"]] <- all.paths
        }
    }
    names(node.paths)
    return(node.paths)
}

#########################################################################################################
highlight.path.each.sbgn <- function(from.node, to.node, glyphs = NULL, arcs = NULL, pathway.id = NULL, 
                                     directed.graph = TRUE, node.set.id.type = NULL, glyphs.id.type = NULL, 
                                     mol.type = NULL, from.node.color = "red", to.node.color = "blue", 
                                     input.node.stroke.width = 4, path.node.color = "black", path.node.stroke.width = 10, 
                                     n.shortest.paths.plot = 1, path.stroke.width = 3, 
                                     shortest.paths.cols = c("purple", "green"), tip.size = 4) {
    
    if (node.set.id.type != glyphs.id.type) {
        from.node <- changeIds(input.ids = from.node, input.type = node.set.id.type, 
            output.type = glyphs.id.type, mol.type = mol.type)
        from.node <- from.node[[1]]
        to.node <- changeIds(input.ids = to.node, input.type = node.set.id.type, 
            output.type = glyphs.id.type, mol.type = mol.type)
        to.node <- to.node[[1]]
        from.node
        to.node
        # new.ids =
        # c('SmallMolecule_96737c854fd379b17cb3b7715570b733','SmallMolecule_7753c3822ee83d806156d21648c931e6')
    } else {
    }
    all.pairs <- shortest.path.from.arcs(arcs = arcs, from.node = from.node, to.node = to.node, 
        directed.graph = directed.graph)
    names(all.pairs)
    if.printed <- 0  # sometimes one pair of input ids can map to multiple pairs of nodes in SBGN map, we want to plot only one pair that a shortest path can be found.
    for (i in seq_len(length.out = length(all.pairs))) {
        all.paths <- all.pairs[[i]]$path  # when only dealing with two nodes, just pick the paths for the first node is good
        if (length(all.paths$res) == 0) {
            warning("no shortest path found!")
            print(names(all.pairs)[i])
            (next)()
        } else {
            if (if.printed == 1) {
                (next)()
            }
            if.printed <- 1
            message("\n found shortest paths for", names(all.pairs)[i], "\n")
            for (p in seq_len(length.out = min(n.shortest.paths.plot, length(all.paths$res)))) {
                vs <- all.paths$res[[p]]
                path.nodes <- vs$name
                arcs <- highlight.edges(node.set = path.nodes, arcs = arcs, node.set.id.type = glyphs.id.type, 
                  arc.node.id.type = glyphs.id.type, arc.highlight.stroke.width = path.stroke.width, 
                  arc.highlight.stroke.color = shortest.paths.cols[p], arc.highlight.tip.size = tip.size)
            }
            glyphs <- highlight.nodes.each.sbgn(node.set = path.nodes, glyphs = glyphs, 
                node.set.id.type = glyphs.id.type, glyphs.id.type = glyphs.id.type, 
                stroke.color = path.node.color, stroke.width = path.node.stroke.width)
            glyphs <- highlight.nodes.each.sbgn(node.set = all.pairs[[i]]$to.node, 
                glyphs = glyphs, node.set.id.type = glyphs.id.type, glyphs.id.type = glyphs.id.type, 
                stroke.color = to.node.color, stroke.width = input.node.stroke.width)
            glyphs <- highlight.nodes.each.sbgn(node.set = all.pairs[[i]]$from.node, 
                glyphs = glyphs, node.set.id.type = glyphs.id.type, glyphs.id.type = glyphs.id.type, 
                stroke.color = from.node.color, stroke.width = input.node.stroke.width)
        }
    }
    return(list(glyphs = glyphs, arcs = arcs))
}

#########################################################################################################
#' Highlight arcs by arc class
#' 
#' This function can modify a SBGNview object's arc. Normally we use it as an argument to \code{+.SBGNview} and modify a SBGNview object. Run \code{help("+.SBGNview")} for more information.
#' 
#' @param class A character string. Default: "arc". The arc class to modify.
#' @param color A character string. Default: "black". The color of arcs with selected 'class'.
#' @param line.width Numeric. Default: 2. Line thickness.
#' @param tip.size Numeric. Default: 6. Tip size. 
#' @return A SBGNview object
#' @examples
#' \dontrun{
#' obj.new <- SBGNview.obj + highlightArcs(class = 'production', color = 'red') 
#' }
#' @export

highlightArcs <- function(class = "all", color = "black", line.width = 2, tip.size = 6) {
    
    function(SBGNview.obj) {
        glyphs.arcs.list <- SBGNview.obj$data
        sbgns <- names(glyphs.arcs.list)
        f <- 0
        for (s in seq_len(length.out = length(glyphs.arcs.list))) {
            # for each sbgn file
            arcs <- glyphs.arcs.list[[s]]$arcs.list
            message("\nSetting edge color for ", sbgns[s], "\n")
            a <- 0
            for (i in seq_len(length.out = length(arcs))) {
                if (arcs[[i]]@arc.class %in% class | class == "all") {
                  arcs[[i]]@edge$line.stroke.color <- color
                  arcs[[i]]@edge$tip.stroke.color <- color
                  arcs[[i]]@edge$tip.fill.color <- color
                  arcs[[i]]@edge$line.width <- line.width
                  arcs[[i]]@edge$tip.stroke.width <- line.width
                  arcs[[i]]@parameters.list$edge.tip.size <- tip.size
                }
            }
            glyphs.arcs.list[[s]]$arcs.list <- arcs
            glyphs.arcs.list[[s]]$render.sbgn.parameters.list$arcs.user <- arcs
        }
        SBGNview.obj$data <- glyphs.arcs.list
        return(SBGNview.obj)
    }
}

#########################################################################################################
