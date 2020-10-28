
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
 
#########################################################################################################
parse.input.sbgn <- function(input.sbgn, output.file, show.pathway.name, sbgn.dir, 
    sbgn.gene.id.type, sbgn.cpd.id.type, sbgn.id.attr, SBGNview.data.folder = "./SBGNview.tmp.data") {
    database <- NULL
    output.file.sbgn <- paste(output.file, "_", input.sbgn, sep = "")
    input.sbgn.full.path <- paste(sbgn.dir, input.sbgn, sep = "/")
    if (!file.exists(input.sbgn.full.path)) {
        # if the SBGN-ML file is not downloaded, in this case, input.sbgn is pathway.id
        message(input.sbgn, " does not look like a local file.\n Using it as pathway ID in 'pathways.info'\n\n ")
        is.input <- pathways.info[, "pathway.id"] == input.sbgn
        if (!any(is.input)) {
            stop("The input ID is not included in the pre-generated SBGN-ML file database! Please check 'data(\"pathways.info\")' for supported pathway IDs!\n\n")
        } else {
            database <- pathways.info[is.input, "database"]
            sub.database <- pathways.info[is.input, "sub.database"]
            
            pathway.name <- pathways.info[is.input, "pathway.name"]
            pathway.name <- gsub("[\\&\\/\\*\\?\\<\\>\\|\\:\"]", "_", pathway.name)
            pathway.name.on.graph <- paste(sub.database, input.sbgn, pathway.name, 
                sep = "::")
            if.file.in.collection <- ""
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
        pathway.name.on.graph <- input.sbgn
        message(input.sbgn, "looks like a local file. Using it directly.\n ")
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
    }
    if (database == "MetaCrop") {
        sbgn.id.attr <- "label"
    } else if (database == "MetaCyc") {
        sbgn.id.attr <- "id"
    } else if (database == "pathwayCommons") {
        sbgn.id.attr <- "id"
    }
    if.file.in.collection <- " Rendered by SBGNview"
    pathway.name.on.graph <- list(pathway.name.on.graph = pathway.name.on.graph, 
        if.file.in.collection = if.file.in.collection)
    
    return(list(input.sbgn.full.path = input.sbgn.full.path, output.file.sbgn = output.file.sbgn, 
        sbgn.gene.id.type = sbgn.gene.id.type, sbgn.cpd.id.type = sbgn.cpd.id.type, 
        pathway.name.on.graph = pathway.name.on.graph, if.file.in.collection = if.file.in.collection, 
        sbgn.id.attr = sbgn.id.attr, database = database))
}

#########################################################################################################
parse.omics.data <- function(gene.data, cpd.data, input.sbgn.full.path, database, 
    user.data.recorded, gene.id.type, cpd.id.type, id.mapping.gene, id.mapping.cpd, 
    node.sum, org, sbgn.gene.id.type, sbgn.cpd.id.type, simulate.data, SBGNview.data.folder = "./SBGNview.tmp.data") {
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
        gene.data.converted <- simulate.user.data(sbgn.file = input.sbgn.full.path, 
            n.samples = 3)
        gene.data.converted <- as.list(as.data.frame(t(gene.data.converted)))
        gene.data.converted[["max.gene"]] <- max(unlist(gene.data.converted), na.rm = TRUE)
        gene.data.converted[["min.gene"]] <- min(unlist(gene.data.converted), na.rm = TRUE)
        gene.data.converted[["max.cpd"]] <- max(unlist(gene.data.converted), na.rm = TRUE)
        gene.data.converted[["min.cpd"]] <- min(unlist(gene.data.converted), na.rm = TRUE)
        cpd.data.converted <- gene.data.converted
        
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
                    output.type = sbgn.gene.id.type, sum.method = node.sum, cpd.or.gene = "gene", 
                    org = org, id.mapping.table = id.mapping.gene, SBGNview.data.folder = SBGNview.data.folder)
                  if (identical(gene.data.converted, "no.id.mapped")) {
                    warning("no id mapped!")
                  }
                }
                gene.data.converted <- as.list(as.data.frame(t(gene.data.converted)))
                gene.data.converted[["max.gene"]] <- max(unlist(gene.data.converted), 
                  na.rm = TRUE)
                gene.data.converted[["min.gene"]] <- min(unlist(gene.data.converted), 
                  na.rm = TRUE)
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
                      output.type = sbgn.cpd.id.type, sum.method = node.sum, cpd.or.gene = "compound", 
                      id.mapping.table = id.mapping.cpd, SBGNview.data.folder = SBGNview.data.folder)
                    if (cpd.data.converted == "no.id.mapped") {
                      warning("no id mapped!")
                    }
                  }
                  cpd.data.converted <- as.list(as.data.frame(t(cpd.data.converted)))
                  cpd.data.converted[["max.cpd"]] <- max(unlist(cpd.data.converted), 
                    na.rm = TRUE)
                  cpd.data.converted[["min.cpd"]] <- min(unlist(cpd.data.converted), 
                    na.rm = TRUE)
                }
            }
            user.data.recorded[[database]] <- list()
            user.data.recorded[[database]][["gene.data"]] <- gene.data.converted
            user.data.recorded[[database]][["cpd.data"]] <- cpd.data.converted
        }
    }
    user.data <- c(gene.data.converted, cpd.data.converted)  # input data is a list, name is sbgn entity id, content is a vector of data values. Therefore the number of values can be different
    return(list(user.data = user.data, user.data.recorded = user.data.recorded))
}

#########################################################################################################
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
find.max.xy <- function(sbgn.xml, arcs.info, color.panel.scale) {
    box.info <- do.call(rbind, xml2::xml_attrs(xml2::xml_find_all(sbgn.xml, ".//bbox")))
    x <- as.numeric(box.info[, "x"])
    w <- as.numeric(box.info[, "w"])
    max.xw <- max(x + w)
    
    y <- as.numeric(box.info[, "y"])
    h <- as.numeric(box.info[, "h"])
    max.yh <- max(y + h)
    
    min.y <- min(y)
    min.x <- min(x)
    
    if (arcs.info == "straight") {
        # if there is no spline edges, calculate margin to move both nodes and edges
        # coordinates. This will give room for color legend
        y.margin <- max(100, max(max.yh, max.xw)/10) * max(1, color.panel.scale)
    } else {
        # if there is routed edges' svg in the SBGN-ML file, we can't move the nodes
        y.margin <- max(100, max(max.xw, max.yh)/5 * color.panel.scale)
    }
    return(list(max.xw = max.xw, max.yh = max.yh, min.y = min.y, min.x = min.x, y.margin = y.margin))
}

#########################################################################################################
download.mapping.file <- function(input.type, output.type, species = NULL, SBGNview.data.folder = "./SBGNview.tmp.data/") {
    options(warn = -1)
    # R CMD check will give different sort result if you didn't specify 'method': the
    # default sort method depends on the locale of your system. And R CMD check uses
    # a different locale than my interactive session.  The issue resided in this: R
    # interactively used:LC_COLLATE=en_US.UTF-8; R CMD check used: LC_COLLATE=C;
    # https://stackoverflow.com/questions/42272119/r-cmd-check-fails-devtoolstest-works-fine
    type.pair.name.1 <- paste(sort(c(input.type, output.type), method = "radix", 
        decreasing = FALSE), collapse = "_")
    type.pair.name.2 <- paste(sort(c(input.type, output.type), method = "radix", 
        decreasing = TRUE), collapse = "_")
    type.pair.name.1.org <- paste(species, type.pair.name.1, sep = "_")
    type.pair.name.2.org <- paste(species, type.pair.name.2, sep = "_")
    try.file.names <- c(type.pair.name.1.org, type.pair.name.2.org, type.pair.name.1, 
        type.pair.name.2)
    if.online.mapping.avai <- FALSE
    if.mapping.in.data.package <- FALSE
    for (i in seq_len(length(try.file.names))) {
        type.pair.name.try <- try.file.names[i]
        if.data.exists <- tryCatch(data(list = type.pair.name.try), warning = function(w) "no such data")
        if (if.data.exists %in% c(type.pair.name.try, "mapping.list", "mapping.table")) {
            mapping.file.name <- type.pair.name.try
            message("ID mapping ", type.pair.name.try, " is included in SBGNview.data package. Loading it.")
            if.mapping.in.data.package <- TRUE
            (break)()
        } else {
            print("ID mapping not in data package for this pair")
            print(type.pair.name.try)
            print(if.data.exists)
        }
    }
    if (!if.mapping.in.data.package) {
        # stop("Can't map this pair of IDs with SBGNview.data!", type.pair.name.try,  "\n")
        warning("Can't map this pair of IDs with SBGNview.data! Generating mapping with Pathview!", type.pair.name.try,  "\n")
        mapping.file.name = "online mapping not available"
    }
    options(warn = 1)
    return(mapping.file.name)
}

#########################################################################################################
#' Generate ID mapping table from input and output ID types. If provided a vector of input IDs (limit.to.ids argument), the function will output mapping table only containing the input IDs. Otherwise, the function will output all IDs of input and output types (restricted to a species if working on gene types and specified  the 'species' parameter).
#' @param input.type A character string. Gene or compound ID type
#' @param output.type A character string. Gene or compound ID type
#' @param species A character string. Three letter KEGG species code.
#' @param cpd.or.gene A character string. Either 'gene' or 'compound'
#' @param limit.to.ids Vector. Molecule IDs of 'input.type'.
#' @param SBGNview.data.folder A character string. 
#' @return A list containing the mapping table. 
#' @examples 
#'  data(mapped.ids)
#'  entrez.to.pathwayCommons = loadMappingTable(
#'                                 input.type = 'ENTREZID'
#'                                 ,output.type = 'pathwayCommons'
#'                                 ,species = 'hsa'
#'                                 ,cpd.or.gene = 'gene'
#'                              )
#'                              
#' @export

loadMappingTable <- function(output.type, input.type, species = NULL, cpd.or.gene, 
                             limit.to.ids = NULL, SBGNview.data.folder = "./SBGNview.tmp.data"){
    
    input.type = gsub("entrez","ENTREZID",input.type)
    output.type = gsub("entrez","ENTREZID",output.type)
    if(input.type == output.type){
        stop("Input type and output types are the same!")
    }
    species = gsub("Hs","hsa",species)
    type.pair.name = paste(sort(c(input.type,output.type),method = "radix",decreasing = TRUE),collapse="_") # R CMD check will give different sort result if we didn't specify "method":  the default sort method depends on the locale of your system. And R CMD check uses a different locale to R interactive session. The issue resided in this: R interactively used:LC_COLLATE=en_US.UTF-8; R CMD check used: LC_COLLATE=C; https://stackoverflow.com/questions/42272119/r-cmd-check-fails-devtoolstest-works-fine
    if(!file.exists(SBGNview.data.folder)){
        dir.create(SBGNview.data.folder)
    }
    
    mapping.file.name = download.mapping.file(input.type, output.type, species, SBGNview.data.folder = SBGNview.data.folder)
    
    if(file.exists(mapping.file.name) | tryCatch(data(list = mapping.file.name), 
                                                 warning = function(w) "no such data") %in% c(mapping.file.name, "mapping.list", "mapping.table" )
       ){
        message("Loading ID mapping file: ",mapping.file.name," \n")
        .GlobalEnv$mapping.list = NULL; .GlobalEnv$mapping.table = NULL
        if.from.file = FALSE
        if(file.exists(mapping.file.name)){
            if.from.file = TRUE
            var.name = load(mapping.file.name)
        }else{
            if(!exists(mapping.file.name)){
                data(list = mapping.file.name)
            }
        }
        if(!is.null(mapping.list)){ # if one of "output.type" or "input.type" is NOT "pathway.id", then the data's names is "mapping.list"
            id.map = mapping.list[[1]][[1]]
        }else if(!is.null(mapping.table)){ # if one of "output.type" or "input.type" is "pathway.id", then the data's names is "mapping.table"
            id.map = mapping.table
        }else if(if.from.file){ # if one of "output.type" or "input.type" is "pathway.id", then the data's names is "mapping.table"
            id.map = get(var.name)
        }else{
            # print("loading mapping file from data package")
            # print(mapping.file.name)
            # data(list = mapping.file.name)
            id.map = get(mapping.file.name)
        }
        message("Finished loading")
        
        if("species" %in% colnames(id.map) &  !is.null(id.map) ){
            if(any(id.map[,"species"] %in% species) ){
                id.map = id.map[id.map[,"species"] %in% species,c(input.type,output.type)]
            }
        }
    }else {
        if(cpd.or.gene == "gene"){
            if(any(c(input.type,output.type) %in% c("pathwayCommons","metacyc.SBGN","pathway.id")) & 
               ! any(c(input.type,output.type) %in% c("KO")) 
            ){
                ko.to.glyph.id = loadMappingTable(
                    input.type = output.type
                    ,output.type = "KO"
                    ,cpd.or.gene = "gene"
                    ,species = species
                    ,SBGNview.data.folder = SBGNview.data.folder
                )
                ko.to.glyph.id = ko.to.glyph.id[[1]][[1]]
                input.to.ko = loadMappingTable(
                    input.type = input.type
                    ,output.type = "KO"
                    ,cpd.or.gene = "gene"
                    ,limit.to.ids = limit.to.ids
                    ,species = species
                    ,SBGNview.data.folder = SBGNview.data.folder
                )
                input.to.ko = input.to.ko$gene[[1]]
                input.to.glyph.id = merge(input.to.ko,ko.to.glyph.id,all= FALSE)
                id.map = input.to.glyph.id[,c(input.type,output.type)]
            }else {
                message("\n\n ID mapping not pre-generated. Try to use Pathview!!!\n\n")
                id.map = geneannot.map.ko(in.ids = limit.to.ids
                                          ,in.type=input.type
                                          ,out.type=output.type
                                          ,species = species
                                          ,unique.map= FALSE
                                          ,SBGNview.data.folder = SBGNview.data.folder
                )
                if(is.vector(id.map)){
                    id.map = as.matrix(t(id.map))
                }
            }
        }else if(cpd.or.gene == "compound"){
            message("\n\n ID mapping not pre-generated. Try to use Pathview cpd!!!\n\n")
            if(is.null(limit.to.ids)){
                stop("Must provide input IDs when using pathview mapping!")
            }
            print(input.type)
            print(output.type)
            id.map =  pathview::cpdidmap(in.ids = limit.to.ids, in.type = input.type, out.type = output.type)
        }else{
            message("\n\n\nCouldn't fine ID mapping table between ",input.type," and ",output.type,"!!!\n")
            message("Tried online resource ",online.mapping.file,"\n")
            message("Please provide ID mapping table using \"id.mapping.table\"!!\n\n")
            stop()
        }
        id.map = id.map[!is.na(id.map[,2]),]
        if(is.vector(id.map)){
            id.map = as.matrix(t(id.map))
        }
        id.map = id.map[!is.na(id.map[,1]),]
        if(is.vector(id.map)){
            id.map = as.matrix(t(id.map))
        }
        id.map = unique(id.map)
    }
    
    if(!is.null(limit.to.ids)){
        # some IDs are quoted, need to remove the quote characters
        id.map[,input.type] = gsub("^\"","",id.map[,input.type])
        id.map[,input.type] = gsub("\"$","",id.map[,input.type])
        id.map = id.map[id.map[,input.type] %in% limit.to.ids,]
    }
    # add additional mapping using KO to glyph.id
    mapping.list = list()
    if(is.vector(id.map)){
        id.map = as.matrix(t(id.map))
    }
    id.map[,1] = as.character(id.map[,1])
    id.map[,2] = as.character(id.map[,2])
    id.map = as.matrix(id.map)
    mapping.list[[cpd.or.gene]]= list()
    mapping.list[[cpd.or.gene]][[type.pair.name]] = id.map
    message("Generated ID mapping list")
    
    return(mapping.list)
}

#########################################################################################################
geneannot.map.ko <- function(in.ids = NULL, in.type, out.type, species = "all", unique.map, 
    SBGNview.data.folder = "./SBGNview.tmp.data") {
    # pathview's geneannot.map can't map KO, so here included KO mapping
    if (!is.null(in.ids)) {
        in.ids <- gsub("\"", "", in.ids)
    }
    message("\nusing pathview for id mapping:", in.type, " to ", out.type, "\n\n")
    if (any(c(in.type, out.type) %in% "KO")) {
        filter.type <- in.type
        out.type <- setdiff(c(in.type, out.type), "KO")
        in.type <- "KO"
        mapping.list <- loadMappingTable(input.type = "KO", output.type = "ENTREZID", 
            cpd.or.gene = "gene", species = species, SBGNview.data.folder = SBGNview.data.folder)
        mapping.table <- mapping.list[[1]][[1]]
        if (out.type %in% c("ENTREZID", "ez", "entrezid")) {
            print("filter species")
            id.map <- mapping.table[mapping.table[, "species"] == species, c(in.type, 
                out.type)]
            if (!is.null(in.ids)) {
                mapping.table <- mapping.table[mapping.table[, filter.type] %in% 
                  in.ids, ]
            }
        } else {
            if (species == "mmu") {
                species <- "mouse"
            }
            output.to.eg <- pathview::eg2id(eg = mapping.table[, "ENTREZID"], category = out.type, 
                org = species, unique.map = unique.map)
            
            output.to.ko <- merge(output.to.eg, mapping.table, all.x = TRUE)
            id.map <- output.to.ko[, c("KO", out.type)]
            id.map <- id.map[!is.na(id.map[, out.type]), ]
            # filter to output only input IDs
            if (!is.null(in.ids)) {
                id.map <- id.map[id.map[, filter.type] %in% in.ids, ]
            }
        }
    } else {
        if (species == "mmu") {
            species <- "mouse"
        }
        if (is.null(in.ids)) {
            stop("Must provide input IDs when using pathview mapping!")
        }
        id.map <- geneannot.map.all(in.ids = in.ids, in.type = in.type, out.type = out.type, 
            org = species, unique.map = unique.map)
        
    }
    return(id.map)
}

#########################################################################################################
load.id.mapping.list.all <- function(SBGN.file.cpd.id.type = NULL, SBGN.file.gene.id.type = NULL, 
    output.gene.id.type = NULL, output.cpd.id.type = NULL, species = NULL, SBGNview.data.folder = "./SBGNview.tmp.data") {
    idmapping.all.list <- list()
    cpd.id.mapping.list <- list()
    gene.id.mapping.list <- list()
    if (!is.null(output.cpd.id.type)) {
        if (SBGN.file.cpd.id.type != output.cpd.id.type) {
            input.cpd.type <- SBGN.file.cpd.id.type
            cpd.id.mapping.list <- loadMappingTable(output.type = output.cpd.id.type, 
                input.type = input.cpd.type, cpd.or.gene = "compound", species = species, 
                SBGNview.data.folder = SBGNview.data.folder)
        }
    }
    
    if (!is.null(output.gene.id.type)) {
        if (SBGN.file.gene.id.type != output.gene.id.type) {
            input.gene.type <- SBGN.file.gene.id.type
            gene.id.mapping.list <- loadMappingTable(output.type = output.gene.id.type, 
                input.type = input.gene.type, cpd.or.gene = "gene", species = species, 
                SBGNview.data.folder = SBGNview.data.folder)
        }
    }
    id.mapping.all.list <- c(cpd.id.mapping.list, gene.id.mapping.list)
    return(id.mapping.all.list)
}

#########################################################################################################
change.id <- function(input.id, cpd.or.gene, input.type, output.type, id.mapping.all.list, 
                      show.ids.for.multiHit = NULL) {
    
    mapping.table <- id.mapping.all.list[[cpd.or.gene]][[1]]
    output.ids <- as.character(mapping.table[mapping.table[, input.type] == input.id, 
        output.type])
    if (any(is.na(output.ids))) {
        print(mapping.table[mapping.table[, input.type] == input.id, ])
        stop("IDs have na")
    }
    if (!is.null(show.ids.for.multiHit)) {
        if (any(tolower(output.ids) %in% tolower(show.ids.for.multiHit))) {
            output.ids <- output.ids[tolower(output.ids) %in% tolower(show.ids.for.multiHit)]
            output.ids <- unique(tolower(output.ids))
            output.ids <- c("mapped", output.ids)
        }
    }
    output.ids <- unique(output.ids)
    output.ids <- paste(output.ids, collapse = "; ")  # If there are multiple id mapped, we combine all of them
    
    return(output.ids)
}

#########################################################################################################
change.glyph.id <- function(glyph.id, glyph.class, id.mapping.all.list, output.gene.id.type = NA, 
    output.cpd.id.type = NA, SBGN.file.cpd.id.type, SBGN.file.gene.id.type, show.ids.for.multiHit = NULL) {
    cpd.or.gene <- glyph.class
    if (glyph.class == "macromolecule") {
        cpd.or.gene <- "gene"
        output.type <- output.gene.id.type
        input.type <- SBGN.file.gene.id.type
    } else if (glyph.class == "simple chemical") {
        cpd.or.gene <- "compound"
        output.type <- output.cpd.id.type
        input.type <- SBGN.file.cpd.id.type
    }
    if (input.type == output.type) {
        return(glyph.id)
    }
    glyph.id <- change.id(input.id = glyph.id, cpd.or.gene, input.type, output.type, 
        id.mapping.all.list, show.ids.for.multiHit = show.ids.for.multiHit)
}

#########################################################################################################
load.all.ids.mapping <- function(database, all.pairs.id.mapping.list, species, output.gene.id.type, 
    output.cpd.id.type, SBGNview.data.folder = "./SBGNview.tmp.data") {
    if (database == "MetaCrop") {
        sbgn.id.attr <- "label"
        SBGN.file.gene.id.type <- "ENZYME"
        SBGN.file.cpd.id.type <- "CompoundName"
    } else if (database == "MetaCyc") {
        sbgn.id.attr <- "id"
        SBGN.file.gene.id.type <- "metacyc.SBGN"
        SBGN.file.cpd.id.type <- "metacyc.SBGN"
    } else if (database == "pathwayCommons") {
        sbgn.id.attr <- "id"
        SBGN.file.gene.id.type <- "pathwayCommons"
        SBGN.file.cpd.id.type <- "pathwayCommons"
    }
    
    if (is.na(output.gene.id.type)) {
        output.gene.id.type.use <- SBGN.file.gene.id.type
    } else {
        output.gene.id.type.use <- output.gene.id.type
    }
    if (is.na(output.cpd.id.type)) {
        output.cpd.id.type.use <- SBGN.file.cpd.id.type
    } else {
        output.cpd.id.type.use <- output.cpd.id.type
    }
    mapped.id.names <- paste(c(SBGN.file.cpd.id.type, SBGN.file.gene.id.type, output.gene.id.type.use, 
        output.cpd.id.type.use), collapse = "_")
    if (mapped.id.names %in% names(all.pairs.id.mapping.list)) {
        id.mapping.all.list <- all.pairs.id.mapping.list[[mapped.id.names]]
    } else {
        id.mapping.all.list <- load.id.mapping.list.all(SBGN.file.cpd.id.type, SBGN.file.gene.id.type, 
            output.gene.id.type.use, output.cpd.id.type.use, species = species, SBGNview.data.folder = SBGNview.data.folder)
        all.pairs.id.mapping.list[[mapped.id.names]] <- id.mapping.all.list
    }
    return(list(all.pairs.id.mapping.list = all.pairs.id.mapping.list, id.mapping.all.list = id.mapping.all.list, 
        output.cpd.id.type.use = output.cpd.id.type.use, output.gene.id.type.use = output.gene.id.type.use, 
        SBGN.file.cpd.id.type = SBGN.file.cpd.id.type, SBGN.file.gene.id.type = SBGN.file.gene.id.type))
}

#########################################################################################################
get.all.nodes.info <- function(sbgn, if.other.id.types.available, output.cpd.id.type.use, 
    output.gene.id.type.use, SBGN.file.cpd.id.type, SBGN.file.gene.id.type, id.mapping.all.list, 
    show.ids.for.multiHit) {
    node.set.list <- list(all.nodes = matrix(ncol = 8, nrow = 0))
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
        node.set.list[["all.nodes"]] <- rbind(node.set.list[["all.nodes"]], nodes.info)
    }
    return(node.set.list)
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
#' @param SBGNview.data.folder A character string. 
#' @param sbgn.dir A character string. The path to a folder that will hold created SBGN-ML files, if the input.sbgn are IDs of pre-collected pathways.
#' @return A list, containing extracted glyph information.
#' @details  The following glyph information is extracted: complex members, compartment members,submap members, node class, nodes with state variables, class of state variables, edges with cardinality, nodes with ports, 'source and sink' nodes, process nodes.\cr When trying to output other ID types, sometimes multiple output IDs are mapped to one glyph. In this situation, the IDs are concatenated by '; ' to represent the glyph.
#' @examples 
#'  data(mapped.ids)
#'  data(sbgn.xmls)
#'  data(pathways.info)
#' node.list = sbgnNodes(
#'                        input.sbgn = 'P00001',
#'                        output.gene.id.type = 'ENTREZID',
#'                        output.cpd.id.type = 'CompoundName',
#'                        species = 'hsa'
#'                     )
#' @export

sbgnNodes <- function(input.sbgn
                      ,output.gene.id.type = NA
                      ,output.cpd.id.type = NA
                      ,database = NA
                      ,species = NA
                      ,show.ids.for.multiHit = NULL
                      ,SBGNview.data.folder = "./SBGNview.tmp.data"
                      ,sbgn.dir = "./SBGNview.tmp.data") {
    if (!file.exists(SBGNview.data.folder)) {
        dir.create(SBGNview.data.folder)
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
                database <- pathways.info[pathways.info[, "pathway.id"] == input.sbgn, 
                  "database"]
                input.sbgn <- downloadSbgnFile(pathway.id = input.sbgn, download.folder = sbgn.dir)
                message("SBGN-ML files downloaded to: ", input.sbgn.full.path)
            } else {
                message("\n", input.sbgn, "looks like an existing local file. Using it directly. \n\n ")
                if (input.sbgn %in% pathways.info[, "file.name"]) {
                  database <- pathways.info[pathways.info[, "file.name"] == input.sbgn, 
                    "database"]
                }
                input.sbgn <- input.sbgn.full.path
            }
        }
        
        if.other.id.types.available <- database %in% c("MetaCrop", "MetaCyc", "pathwayCommons")
        if (if.other.id.types.available) {
            all.id.mapping.result <- load.all.ids.mapping(database, all.pairs.id.mapping.list, 
            species, output.gene.id.type, output.cpd.id.type, SBGNview.data.folder = SBGNview.data.folder)
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
        node.set.list <- get.all.nodes.info(sbgn, if.other.id.types.available, output.cpd.id.type.use, 
            output.gene.id.type.use, SBGN.file.cpd.id.type, SBGN.file.gene.id.type, 
            id.mapping.all.list, show.ids.for.multiHit)
        result.list[[input.sbgn.original]] <- node.set.list[["all.nodes"]]
    }
    return(result.list)
}

#########################################################################################################
# parse ports
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
            arc <- new(paste(child.attrs["class"], ".sbgn.arc", sep = ""), id = paste(spline.info["id"], 
                "source.arc", sep = ":"), start.x = as.numeric(child.attrs["start.x"]), 
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
            arc <- new(paste(child.attrs["class"], ".sbgn.arc", sep = ""), id = paste(spline.info["id"], 
                "target.arc", sep = ":"), start.x = as.numeric(child.attrs["start.x"]), 
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
parse.splines <- function(sbgn.xml, glyphs, if.plot.svg = TRUE, y.margin = 0, global.parameters.list, 
    arcs.user = list()) {
    if (length(arcs.user) == 0) {
        svg.splines <- ""
        splines.list <- list()
        arc.splines <- xml2::xml_find_all(sbgn.xml, ".//edge.spline.info")
        arc.splines <- xml2::xml_find_all(arc.splines[[length(arc.splines)]], ".//arc.spline")  # if there are more than one version of routed edges, we use the latest version
        
        splines <- rep(NULL, times = length(arc.splines))
        for (i in seq_along(arc.splines)) {
            arc.spline <- arc.splines[[i]]
            spline.info <- xml2::xml_attrs(arc.spline)  # extract attributes of this glyph
            spline.arc <- new("spline.arc")
            spline.arc@id <- paste(spline.info["id"], spline.info["class"], sep = "==")
            # spline.arc@id = spline.info['id']
            spline.arc@source <- spline.info["source"]
            spline.arc@target <- spline.info["target"]
            spline.arc@arc.class <- spline.info["class"]
            spline.arc@global.parameters.list <- global.parameters.list
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
            spline.arc <- arcs.user[[i]]
            if (if.plot.svg) {
                spline.arc.svg <- plot.arc(spline.arc)
            } else {
                spline.arc.svg <- ""
            }
            svg.splines <- paste(svg.splines, spline.arc.svg, sep = "\n")
        }
    }
    return(list(svg.splines = svg.splines, splines.list = splines.list))
}

#########################################################################################################
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
            arc <- new("next.sbgn.arc", id = paste(arc.line["id"], arc.line["start.x"], 
                sep = "_"), start.x = as.numeric(arc.line["start.x"]), start.y = as.numeric(arc.line["start.y"]), 
                end.x = as.numeric(arc.line["end.x"]), end.y = as.numeric(arc.line["end.y"]))
            arc@global.parameters.list <- global.parameters.list
            arc@arc.class <- "next"
            arc@edge <- edge.paras
            
            arcs.list[[arc@id]] <- arc
            svg.arc <- paste(svg.arc, plot.arc(arc), sep = "\n")
            arc.line["start.x"] <- coordinates["x"]  # after ploting this 'next', set the new_arc's start coordinates
            arc.line["start.y"] <- coordinates["y"]
            
        } else if (xml2::xml_name(child) == "end") {
            arc.line["end.x"] <- coordinates["x"]
            arc.line["end.y"] <- coordinates["y"]
            if (sum(as.numeric(arc.line[c("start.x", "start.y", "end.x", "end.y")])) == 
                0) {
                arc.line <- find.arc.coordinates(arc.line, glyphs)
            }
            if (is.na(arc.line["id"])) {
                arc.line["id"] <- paste(arc.info["source"], arc.info["target"], sep = "->")
            }
            arc <- new(paste(arc.class, ".sbgn.arc", sep = ""), id = paste(arc.line["id"], 
                arc.line["start.x"], sep = "_"), start.x = as.numeric(arc.line["start.x"]), 
                start.y = as.numeric(arc.line["start.y"]), end.x = as.numeric(arc.line["end.x"]), 
                end.y = as.numeric(arc.line["end.y"]))
            arc@arc.class <- arc.class
            arc@source <- arc.info["source"]
            arc@target <- arc.info["target"]
            arc@global.parameters.list <- global.parameters.list
            arc@edge <- edge.paras
            arcs.list[[arc@id]] <- arc
            svg.arc <- paste(svg.arc, plot.arc(arc), sep = "\n")
        }
    }
    return(list(arcs.list = arcs.list, svg.arc = svg.arc))
    
}

#########################################################################################################
# parse arcs
parse.arcs <- function(sbgn.xml, glyphs, if.plot.svg = TRUE, y.margin = 0, global.parameters.list, 
    arcs.user = list()) {
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
generate.glyph.id <- function(glyph.id, glyph.class, glyph, node.set.list, no.id.node.count.list) {
    if (is.na(glyph.id)) {
        if (glyph.class %in% c("unit of information", "state variable")) {
            glyph.parent <- xml2::xml_parent(glyph)
            parent.id <- xml2::xml_attr(glyph.parent, "id")
            if (is.null(node.set.list[["molecules.with.state.variables"]][[parent.id]])) {
                node.set.list[["molecules.with.state.variables"]][[parent.id]] <- 1
            } else {
                node.set.list[["molecules.with.state.variables"]][[parent.id]] <- node.set.list[["molecules.with.state.variables"]][[parent.id]] + 
                  1
            }
            v.i <- node.set.list[["molecules.with.state.variables"]][[parent.id]]
            glyph.id <- paste(parent.id, v.i, sep = ".info.")
            
        } else if (!glyph.class %in% names(no.id.node.count.list)) {
            # some nodes don't have id(cardinality), we need to gennerate an id for them
            no.id.node.count.list[[glyph.class]] <- 0
        } else {
            no.id.node.count.list[[glyph.class]] <- no.id.node.count.list[[glyph.class]] + 
                1
        }
        index.id <- no.id.node.count.list[[glyph.class]]
        glyph.id <- paste(glyph.class, index.id, sep = ":")
    }
    return(list(glyph.id = glyph.id, node.set.list = node.set.list, no.id.node.count.list = no.id.node.count.list))
}

#########################################################################################################
add.omics.data.to.glyph <- function(glyph.info, glyph, node, sbgn.id.attr, user.data) {
    node.omics.data.id <- glyph.info[sbgn.id.attr]
    # remove complex name from omics.data.id for metacyc
    node.omics.data.id.without.complex <- gsub("_Complex.+:@:", ":@:", node.omics.data.id)
    node.omics.data.id.without.complex <- gsub("_Complex_.+", "", node.omics.data.id)
    
    if (!xml2::xml_attr(xml2::xml_parent(glyph), "class") %in% c("complex", "submap")) {
        # molecules within a complex sometimes have different ids, so we don't count them
        # when calculating the mapped nodes
    }
    ################################################################ 
    if (node.omics.data.id %in% names(user.data)) {
        node@user.data <- user.data[[node.omics.data.id]]
        if (!xml2::xml_attr(xml2::xml_parent(glyph), "class") %in% c("complex", "submap")) {
            # molecules within a complex sometimes have different ids, so we don't count them
            # when calculating the mapped nodes
            
        }
    } else if (node.omics.data.id.without.complex %in% names(user.data)) {
        node@user.data <- user.data[[node.omics.data.id.without.complex]]
        if (!xml2::xml_attr(xml2::xml_parent(glyph), "class") %in% c("complex", "submap")) {
            # molecules within a complex sometimes have different ids, so we don't count them
            # when calculating the mapped nodes
            
        }
    } else {
        node@user.data <- c("no.user.data")
    }
    if (length(node@user.data) == 1) {
        node@user.data <- as.matrix(t(c(node@user.data, node@user.data)))
    }
    return(list(node = node))
}

#########################################################################################################
parse.glyph <- function(sbgn.xml, user.data, if.plot.svg = TRUE, y.margin = 0, max.x, 
    global.parameters.list, sbgn.id.attr, glyphs.user = list(), compartment.layer.info, 
    if.plot.cardinality) {
    
    if.plot.annotation.nodes <- global.parameters.list$if.plot.annotation.nodes
    if.use.number.for.long.label <- global.parameters.list$if.use.number.for.long.label
    # parse glyphs and plot glyphs and ports
    map.language <- xml2::xml_attrs(xml2::xml_find_all(sbgn.xml, ".//map")[[1]])["language"]
    message("\n Map language is ", map.language, "\n")
    if (is.null(map.language)) {
        map.language <- ""
    }
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
        
        ##### use user glyph if found
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
            parse.result <- generate.node.obj(glyph, glyph.class, glyph.info, node, 
                if.plot.svg, y.margin, sbgn.id.attr, user.data, max.x, global.parameters.list, 
                if.use.number.for.long.label, if.plot.annotation.nodes, map.language, 
                long.words.count.list, shorter.label.mapping.list)
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
generate.node.obj <- function(glyph, glyph.class, glyph.info, node, if.plot.svg, 
    y.margin, sbgn.id.attr, user.data, max.x, global.parameters.list, if.use.number.for.long.label, 
    if.plot.annotation.nodes, map.language, long.words.count.list, shorter.label.mapping.list) {
    
    # parse children of the node, get information like label, coordinates, if complex
    # empty etc.
    parsed.children <- parse.glyph.children(map.language, glyph, glyph.class, glyph.info, 
        node, if.plot.svg, y.margin)
    node <- parsed.children$node
    node.clone <- parsed.children$node.clone
    # glyph.box.information = parsed.children$glyph.box.information
    glyph.port.info <- parsed.children$glyph.port.info
    svg.port <- parsed.children$svg.port
    node.label <- parsed.children$node.label
    glyph.info <- parsed.children$glyph.info
    if.complex.empty <- parsed.children$if.complex.empty
    
    node@compartment <- glyph.info["compartmentRef"]
    if (length(node@label) == 0) {
        # if the node has no label (character(o)), assign it to ''.
        node@label <- ""
    }
    
    # add omics data to node object, and record mapping result
    mapping.result <- add.omics.data.to.glyph(glyph.info, glyph, node, sbgn.id.attr, 
        user.data)
    node <- mapping.result$node
    
    if (glyph.class == "compartment") {
        node@max.x <- max.x
    }
    
    if (!is.na(glyph.info["orientation"])) {
        node@orientation <- glyph.info["orientation"]
    }
    
    result.list <- break.text.into.segments(label = node@label, w = node@w, glyph.class = glyph.class, 
        global.parameters.list = global.parameters.list, max.x = max.x, glyph = node)
    if.long.word <- result.list$if.long.word
    label.margin <- result.list$label.margin
    node@label.margin <- label.margin
    if (if.use.number.for.long.label & if.long.word) {
        if (!glyph.class %in% names(long.words.count.list)) {
            long.words.count.list[[glyph.class]] <- 1
        } else {
            long.words.count.list[[glyph.class]] <- long.words.count.list[[glyph.class]] + 
                1
        }
        index.long.words <- long.words.count.list[[glyph.class]]
        shorter.label <- paste(glyph.class, index.long.words, sep = "_")
        shorter.label.mapping.list <- rbind(shorter.label.mapping.list, c(shorter.label, 
            node@label))
        node@label <- shorter.label
    }
    
    node@global.parameters.list <- global.parameters.list
    # handle clone markers
    
    ##### set node specific parameters
    node@shape$stroke.width <- min(1, max.x/900)
    if (glyph.class == "annotation" & !if.plot.annotation.nodes) {
        node@stroke.opacity <- 0
    } else if (is(node, "complex.sbgn")) {
        node@if.complex.empty <- if.complex.empty
        node@shape$stroke.width <- min(3, max.x/300)
    } else if (is(node, "compartment.sbgn")) {
        # if this node is a clone marker
        node@shape$stroke.width <- min(6, max.x/150)
        # node@w = node@w + 6 node@h = node@h + 6 if this node is a clone marker
    } else if (!is.character(node.clone)) {
        node.clone@x <- node@x
        node.clone@y <- node@y
        node.clone@glyph.class <- glyph.class
        node.clone@global.parameters.list <- global.parameters.list
        node.clone@w <- node@w
        node.clone@h <- node@h
        node.clone@compartment <- glyph.info["compartmentRef"]
        node@clone <- list(node.clone)
    } else if (node@glyph.class %in% c("process", "uncertain process", "omitted process")) {
        node@if.show.label <- FALSE
    }
    return(list(node = node, long.words.count.list = long.words.count.list, shorter.label.mapping.list = shorter.label.mapping.list))
}

#########################################################################################################
parse.glyph.children <- function(map.language, glyph, glyph.class, glyph.info, node, 
    if.plot.svg, y.margin) {
    
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
                node@label <- paste(glyph.label["value"], "@", glyph.label["variable"], 
                  sep = "")
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
            glyph.box.information["y"] <- as.numeric(glyph.box.information["y"]) + 
                y.margin
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
            if (as.numeric(glyph.port.info["x"]) + as.numeric(glyph.port.info["y"]) != 
                0) {
                svg.port <- paste(svg.port, plot.arc.ports(glyph.port.info, node), 
                  sep = "\n")
            }
            node@svg.port <- svg.port
        } else if (xml2::xml_name(child) == "glyph" & xml2::xml_attr(child, "class") != 
            "annotation") {
            if.complex.empty <- FALSE
        }
    }
 
    return(list(node = node, node.clone = node.clone, glyph.port.info = glyph.port.info, 
        svg.port = svg.port, node.label = node.label, glyph.info = glyph.info, if.complex.empty = if.complex.empty))
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
        }, numeric(1))
    })
    return(user.data)
}

#########################################################################################################
node.ids.from.sbgn <- function(sbgn.file, output.glyph.class = c("macromolecule", 
    "simple chemical"), if.include.complex.member = FALSE) {
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
# the first column of id.map is input ids(need to convert FROM), second column is output ids(converting TO)
mol.sum.multiple.mapping = function(mol.data, id.map, sum.method, input.sbgn.file = "na"){   
    # function: change input data matrix with input/uniprot id to data matrix with pathwaycommons id
    # if there are multiple input ids mapped to a node, collapse their values using user specified function(sum.method)
    id.map[,1] = as.character(id.map[,1])
    id.map[,2] = as.character(id.map[,2])
    if(is.vector(mol.data)){
        mol.data = as.matrix(mol.data)
    }
    
    # mol.data = mol.data[!duplicated(tolower(row.names(mol.data))),]
    if(!is.matrix(mol.data)){
        mol.data = as.matrix(mol.data)
    }
    # row.names(mol.data) = tolower(row.names(mol.data))
    input.ids = row.names(mol.data)
    if(is.vector(id.map)){
        id.map = as.matrix(t(id.map))
    }
    # id.map[,1] = tolower(as.character(id.map[,1]))
    # id.map[,2] = as.character(id.map[,2])
    
    # if input sbgn file is provided, we just convert the ids that are in the file, to speed up
    if(input.sbgn.file != "na"){
        all.nodes = node.ids.from.sbgn(input.sbgn.file)
    }else{
        all.nodes = id.map[,2]
    }
    if.id.mapping.in.sbgn.file =id.map[,2] %in% all.nodes 
    # head(if.id.mapping.in.sbgn.file)
    if(!any(if.id.mapping.in.sbgn.file)){
        print("no sbgn ids in mapping file!")
        print(head(id.map))
        print(head(all.nodes))
        return("no.id.mapped")
    }
    id.map.in.target = id.map[if.id.mapping.in.sbgn.file,]  # if input sbgn file is provided, we just convert the ids that are in the file, to speed up
    # head(id.map.in.target)
    
    if(is.vector(id.map.in.target)){
        id.map.in.target = as.matrix(t(id.map.in.target))
    }
    # if.id.mapping.in.source = tolower(id.map.in.target[,1]) %in% tolower(input.ids)
    if.id.mapping.in.source = id.map.in.target[,1] %in% input.ids
    
    if(!any(if.id.mapping.in.source)){
        print("no source ids in mapping file!")
        return("no.id.mapped")
    }
    #############################################
    
    id.map.in.target.and.source = id.map.in.target[if.id.mapping.in.source,]
    # head(id.map.in.target.and.source)
    if(is.vector(id.map.in.target.and.source)){
        id.map.in.target.and.source = as.matrix(t(id.map.in.target.and.source))
    }
    
    class(id.map.in.target.and.source[,1])
    # head(mol.data)
    in.data.source.id = mol.data[id.map.in.target.and.source[,1],]  # the rows are in sequence of sbgn.id
    if(is.vector(in.data.source.id)){
        in.data.source.id = as.matrix(in.data.source.id)
    }
    message("Merging molecules with the same output ID")
    in.data.target.id = merge.molecule.data(in.data.source.id
                                            ,id.map.in.target.and.source
                                            ,sum.method
    )
    message("Merging finished")
    return(in.data.target.id)
}

#########################################################################################################
merge.molecule.data = function(in.data.source.id, id.map.in.target.and.source, sum.method){
    
    target.ids = id.map.in.target.and.source[,2]
    target.id.count = table(target.ids)
    
    if.target.id.single = target.id.count[target.ids] == 1
    if.target.id.multi = target.id.count[target.ids] > 1
    
    in.data.single.target.id = in.data.source.id[if.target.id.single,]
    if(is.vector(in.data.single.target.id)){
        in.data.single.target.id = as.matrix(in.data.single.target.id)
    }
    single.target.ids = target.ids[if.target.id.single]
    row.names(in.data.single.target.id) = single.target.ids
    
    # handel multiple output IDs to a single input ID
    if(any(if.target.id.multi)){
        in.data.multi.target.id = in.data.source.id[if.target.id.multi,]
        multi.target.ids = target.ids[if.target.id.multi]
        if(is.vector(in.data.multi.target.id)){
            in.data.multi.target.id = as.matrix(in.data.multi.target.id)
        }
        
        in.data.target.id = by(
            in.data.multi.target.id  # the rows are in sequence of sbgn.id
            ,as.factor(multi.target.ids)  # the rows are in sequence of sbgn.id
            ,function(data.same.id){
                # data.same.id = as.numeric(data.same.id)
                if(is.vector(data.same.id)){
                    return(data.same.id)
                }else{
                    sumed = apply(
                        data.same.id,2
                        ,FUN = sum.method
                    )
                }
                return(sumed)
            }
            ,simplify= TRUE
        )
        in.data.target.id = as.list(in.data.target.id)
        if(!is.list(in.data.target.id)){
            message("No data mapped! From ",colnames(id.map),"\n")
            return("no.id.mapped")
        }
        message("Combinding merging result")
        in.data.target.id = do.call(rbind,in.data.target.id)
        in.data.target.id = rbind(in.data.target.id,in.data.single.target.id)
    }else{
        in.data.target.id = in.data.single.target.id
    }
    
    return(in.data.target.id)
}

#########################################################################################################
# mol.sum.multiple.mapping <- function(mol.data, id.map, sum.method, input.sbgn.file = "na") {
#     # the first column of id.map is input ids(need to convert FROM), second column is
#     # output ids(converting TO) function: change input data matrix with input/uniprot
#     # id to data matrix with pathwaycommons id if there are multiple input ids mapped
#     # to a node, collapse their values using user specified function(sum.method)
#     if (is.vector(mol.data)) {
#         mol.data <- as.matrix(mol.data)
#     }
#     
#     # mol.data = mol.data[!duplicated(tolower(row.names(mol.data))),]
#     if (!is.matrix(mol.data)) {
#         mol.data <- as.matrix(mol.data)
#     }
#     # row.names(mol.data) = tolower(row.names(mol.data))
#     input.ids <- row.names(mol.data)
#     if (is.vector(id.map)) {
#         id.map <- as.matrix(t(id.map))
#     }
#     # id.map[,1] = tolower(as.character(id.map[,1])) id.map[,2] =
#     # as.character(id.map[,2])
#     
#     # if input sbgn file is provided, we just convert the ids that are in the file,
#     # to speed up
#     if (input.sbgn.file != "na") {
#         all.nodes <- node.ids.from.sbgn(input.sbgn.file)
#     } else {
#         all.nodes <- id.map[, 2]
#     }
#     if.id.mapping.in.sbgn.file <- id.map[, 2] %in% all.nodes
#     if (!any(if.id.mapping.in.sbgn.file)) {
#         print("no sbgn ids in mapping file!")
#         print(head(id.map))
#         print(head(all.nodes))
#         return("no.id.mapped")
#     }
#     id.map.in.target <- id.map[if.id.mapping.in.sbgn.file, ]  # if input sbgn file is provided, we just convert the ids that are in the file, to speed up
#     
#     if (is.vector(id.map.in.target)) {
#         id.map.in.target <- as.matrix(t(id.map.in.target))
#     }
#     # if.id.mapping.in.source = tolower(id.map.in.target[,1]) %in% tolower(input.ids)
#     if.id.mapping.in.source <- id.map.in.target[, 1] %in% input.ids
#     if (!any(if.id.mapping.in.source)) {
#         print("no source ids in mapping file!")
#         return("no.id.mapped")
#     }
#     id.map.in.target.and.source <- id.map.in.target[if.id.mapping.in.source, ]
#     if (is.vector(id.map.in.target.and.source)) {
#         id.map.in.target.and.source <- as.matrix(t(id.map.in.target.and.source))
#     }
#     in.data.source.id <- mol.data[id.map.in.target.and.source[, 1], ]  # the rows are in sequence of sbgn.id
#     if (is.vector(in.data.source.id)) {
#         in.data.source.id <- as.matrix(in.data.source.id)
#     }
#     message("Merging molecules with the same output ID")
#     in.data.target.id <- merge.molecule.data(in.data.source.id, id.map.in.target.and.source, 
#         sum.method)
#     message("Merging finished")
#     return(in.data.target.id)
# }
# 
# merge.molecule.data.each <- function(in.data.source.id, id.map.in.target.and.source, 
#     sum.method) {
#     target.ids <- id.map.in.target.and.source[, 2]
#     target.id.count <- table(target.ids)
#     
#     if.target.id.single <- target.id.count[target.ids] == 1
#     if.target.id.multi <- target.id.count[target.ids] > 1
#     
#     in.data.single.target.id <- in.data.source.id[if.target.id.single, ]
#     if (is.vector(in.data.single.target.id)) {
#         in.data.single.target.id <- as.matrix(in.data.single.target.id)
#     }
#     single.target.ids <- target.ids[if.target.id.single]
#     row.names(in.data.single.target.id) <- single.target.ids
#     
#     # handel multiple output IDs to a single input ID
#     if (any(if.target.id.multi)) {
#         in.data.multi.target.id <- in.data.source.id[if.target.id.multi, ]
#         multi.target.ids <- target.ids[if.target.id.multi]
#         if (is.vector(in.data.multi.target.id)) {
#             in.data.multi.target.id <- as.matrix(in.data.multi.target.id)
#         }
#         
#         in.data.target.id <- by(in.data.multi.target.id  # the rows are in sequence of sbgn.id
# , 
#             as.factor(multi.target.ids)  # the rows are in sequence of sbgn.id
# , 
#             function(data.same.id) {
#                 # data.same.id = as.numeric(data.same.id)
#                 if (is.vector(data.same.id)) {
#                   return(data.same.id)
#                 } else {
#                   sumed <- apply(data.same.id, 2, FUN = sum.method)
#                 }
#                 return(sumed)
#             }, simplify = TRUE)
#         in.data.target.id <- as.list(in.data.target.id)
#         if (!is.list(in.data.target.id)) {
#             message("No data mapped! From ", colnames(id.map.in.target.and.source), 
#                 "\n")
#             return("no.id.mapped")
#         }
#         message("Combinding merging result")
#         in.data.target.id <- do.call(rbind, in.data.target.id)
#         in.data.target.id <- rbind(in.data.target.id, in.data.single.target.id)
#     } else {
#         in.data.target.id <- in.data.single.target.id
#     }
#     
#     
#     return(in.data.target.id)
# }
# 
# 
# merge.molecule.data <- function(in.data.source.id, id.map.in.target.and.source, sum.method) {
#     in.data.target.id <- by(in.data.source.id  # the rows are in sequence of sbgn.id
# , 
#         as.factor(id.map.in.target.and.source[, 2])  # the rows are in sequence of sbgn.id
# , 
#         function(data.same.id) {
#             # data.same.id = as.numeric(data.same.id)
#             if (is.vector(data.same.id)) {
#                 return(data.same.id)
#             } else {
#                 apply(data.same.id, 2, FUN = sum.method)
#             }
#         }, simplify = TRUE)
#     in.data.target.id <- as.list(in.data.target.id)
#     if (!is.list(in.data.target.id)) {
#         message("No data mapped! From ", colnames(id.map.in.target.and.source), "\n")
#         return("no.id.mapped")
#     }
#     message("Combinding merging result")
#     in.data.target.id <- do.call(rbind, in.data.target.id)
#     
#     return(in.data.target.id)
# }

#########################################################################################################
#' Download pre-generated SBGN-ML file from GitHub
#' 
#' This function can generate a SBGN-ML file of our pre-collected SBGN-ML files 
#' 
#' @param pathway.id The ID of pathway. For accepted pathway IDs, please check \code{data('pathways.info')}. IDs are in column 'pathway.id' (pathways.info[,'pathway.id'])
#' @param download.folder The output folder to store created SBGN-ML files.
#' @return A vector of character strings. The path to the created SBGN-ML files.
#' @examples 
#' data('pathways.info')
#' data(sbgn.xmls)
#' input.sbgn = downloadSbgnFile(
#'                   pathway.id = pathways.info[1,'pathway.id'],
#'                   download.folder = './')
#' @export

downloadSbgnFile <- function(pathway.id, download.folder = ".") {
    if (!file.exists(download.folder)) {
        dir.create(download.folder)
    }
    if (any(pathway.id %in% c("AF_Reference_Card.sbgn", "PD_Reference_Card.sbgn", 
        "ER_Reference_Card.sbgn"))) {
        sbgn.file.names <- pathway.id
    } else {
        sbgn.file.names <- pathways.info[pathways.info[, "pathway.id"] %in% pathway.id, 
            "file.name"]
    }
    if (length(sbgn.file.names) == 0) {
        print("pathway.id")
        stop("Only pathway IDs in pathways.info[,\"pathway.id\"] are supported!!!\n ")
    }
    database.name <- pathways.info[pathways.info[, "pathway.id"] %in% pathway.id, 
        "database"]
    output.files <- c()
    for (i in seq_len(length.out = length(sbgn.file.names))) {
        sbgn.file.name <- gsub("\"", "", sbgn.file.names[i])
        output.file <- paste(download.folder, "/", sbgn.file.name, sep = "")
        if (!file.exists(output.file)) {
            if (sbgn.file.name %in% names(sbgn.xmls)) {
                write(sbgn.xmls[[sbgn.file.name]], file = output.file)
            } else {
                stop(sbgn.file.name, " is not included in SBGNview.data package!")
            }
            
        } else {
            message("SBGN file exists:")
            message(output.file)
        }
        output.files <- c(output.files, (as.character(output.file)))
    }
    output.files
    return(output.files)
}

#########################################################################################################
#' Change the data IDs of input omics data
#' 
#' This function changes the IDs of input omics data from one type to another. 
#' 
#' @param data.input.id A matrix. Input omics data. Rows are genes or compounds, columns are measurements. Row names are the original IDs that need to be transformed.
#' @param input.type A character string. The type of input IDs. Please check \code{data('mapped.ids')} for supported types.
#' @param output.type A character string. The type of output IDs. Please check \code{data('mapped.ids')} for supported types. 
#' @param sum.method  A character string. In some cases multiple input IDs are mapped to one output ID. In this situation ,we may need to derive only one value from them. This parameter is a function that can derive a single numeric value from a vector of numeric values (e.g. 'sum','max','min','mean'), including a User Defined Function (UDF).
#' @param org  A character string. The species source of omics data. 'changeDataId' uses pathview to map between some gene ID types. Please use '?geneannot.map' to check the detail. Pathview needs species information to do the job. This parameter is a two-letter abbreviation of organism name, or KEGG species code, or the common species name, used to determine the gene annotation package. For all potential values check: data(bods); bods. Default org='Hs', and can also be 'hsa' or 'human' (case insensitive). 
#' @param cpd.or.gene  A character string. Either 'compound' or 'gene' -- the type of input omics data. 
#' @param id.mapping.table A matrix.  Mapping table between input.type and output.type. This matrix should have two columns for input.type and output.type, respectively.  Column names should be the values of parameters 'input.type' and 'output.type'. See example section for an example. 
#' @param SBGNview.data.folder A character string. The path to a folder that will hold download ID mapping files and pathway information data files. 
#' @return A matrix, row names are changed to IDs of 'output.type'. Note the number of rows may be different from input matrix, because multiple input IDs could be collapsed to a single output ID. Also a single input ID could map to multiple output IDs.
#' @details  This function maps between various gene/compound ID types.
#' 
#'  1. Map other ID types to glyph IDs in SBGN-ML files of pathwayCommons database,and MetaCyc database:
#'     Use output.type = 'pathwayCommons' or output.type = 'metacyc.SBGN'. 
#'     Please check \code{data('mapped.ids')} for supported input ID types.
#'     
#'  2. Map between other ID types:
#'  
#'  2.1 ID types pairs can be mapped by pathview.
#'     Currently SBGNview uses pathview to do this mapping. Please check pathview functions 'geneannot.map' and 'cpdidmap' for more details.
#'     
#'  2.2 Other ID type pairs
#'     
#'     In this case, users need to provide id.mapping.table.
#' @examples 
#' # Change gene ID
#' data(mapped.ids)
#' library(pathview)
#' data('gse16873.d')
#' gene.data = gse16873.d[c('7157','1032'),]
#' mapping.table = data.frame(
#'                    ENTREZID = c('7157','1032')
#'                    ,SYMBOL = c('TP53','CDKN2D')
#'                    ,stringsAsFactors=FALSE
#' )
#' new.dt = changeDataId(
#'       data.input.id = gene.data,
#'       output.type='SYMBOL',
#'       input.type='ENTREZID',
#'       cpd.or.gene = 'gene',
#'       id.mapping.table = mapping.table
#'       )
#'       
#' 
#' @export    

changeDataId <- function(data.input.id, input.type, output.type, sum.method = "sum", 
    org = "hsa", cpd.or.gene, id.mapping.table = NULL, SBGNview.data.folder = "./SBGNview.tmp.data") {
    # function: change input data matrix with input/uniprot ID to data matrix with
    # pathwaycommons id if there are multiple input IDs mapped to a node, collapse
    # their values using user specified function(sum.method)
    input.type <- gsub("entrez", "ENTREZID", input.type)
    output.type <- gsub("entrez", "ENTREZID", output.type)
    if (is.null(id.mapping.table)) {
        # if user didn't provide mapping table, we try to download one.
        mapping.list <- loadMappingTable(output.type = output.type, input.type = input.type, 
            species = org, cpd.or.gene = cpd.or.gene, limit.to.ids = row.names(data.input.id), 
            SBGNview.data.folder = SBGNview.data.folder)
        id.map <- mapping.list[[1]][[1]]
        id.map <- as.matrix(id.map[, c(input.type, output.type)])
    } else {
        if (!all(c(input.type) %in% colnames(id.mapping.table))) {
            message(input.type, "must be in column names of id.mapping.table!\n")
            message("Column names of id.mapping.table are:\n")
            cat(colnames(id.mapping.table), "\n\n\n\n")
            stop()
        }
        if (!all(c(output.type) %in% colnames(id.mapping.table))) {
            message(output.type, "must be in column names of id.mapping.table!\n")
            message("Column names of id.mapping.table are:\n")
            cat(colnames(id.mapping.table), "\n\n\n\n")
            stop()
        }
        id.map <- id.mapping.table
    }
    message("Changing data IDs")
    in.data.target.id <- mol.sum.multiple.mapping(mol.data = data.input.id, id.map = id.map, 
        sum.method = sum.method)
    message("Finished changing data IDs")
    return(in.data.target.id)
}

#########################################################################################################
geneannot.map.all <- function(in.ids, in.type, out.type, org = "Hs", pkg.name = NULL, 
    unique.map = TRUE, na.rm = TRUE, keep.order = FALSE) {
    
    if (is.null(pkg.name)) {
        # pkg.name=paste('org', org, 'eg.db', sep='.')
        data(bods)
        ridx <- grep(tolower(paste0(org, "[.]")), tolower(bods[, 1]))
        if (length(ridx) == 0) {
            ridx <- grep(tolower(org), tolower(bods[, 2:3]))%%nrow(bods)
            if (length(ridx) == 0) 
                stop("Wrong org value!")
            if (any(ridx == 0)) 
                ridx[ridx == 0] <- nrow(bods)
        }
        pkg.name <- bods[ridx, 1]
    }
    all.mappings <- character(2)
    for (i in seq_len(length.out = nrow(bods))) {
        if (bods[i, 1] != pkg.name) {
            (next)()
        }
        pkg.name <- bods[i, 1]
        
        if (!pkg.name %in% rownames(installed.packages())) {
            (next)()
        }
        message("Using package:", pkg.name, "\n")
        db.obj <- eval(parse(text = paste0(pkg.name, "::", pkg.name)))
        id.types <- AnnotationDbi::columns(db.obj)  #columns(eval(as.name(pkg.name)))
        
        in.type <- toupper(in.type)
        out.type <- toupper(out.type)
        eii <- in.type == toupper("entrez") | in.type == toupper("eg")
        if (any(eii)) 
            in.type[eii] <- "entrez"
        eio <- out.type == toupper("entrez") | out.type == toupper("eg")
        if (any(eio)) 
            out.type[eio] <- "entrez"
        if (in.type == out.type) 
            stop("in.type and out.type are the same, no need to map!")
        
        nin <- length(in.type)
        if (nin != 1) 
            stop("in.type must be of length 1!")
        out.type <- out.type[!out.type %in% in.type]
        nout <- length(out.type)
        
        msg <- paste0("must from: ", paste(id.types, collapse = ", "), "!")
        if (!in.type %in% id.types) 
            stop("'in.type' ", msg)
        if (!all(out.type %in% id.types)) 
            stop("'out.type' ", msg)
        
        in.ids0 <- in.ids
        in.ids <- unique(as.character(in.ids))  #unique necessary for select()# if(unique.map)
        out.ids <- character(length(in.ids))
        res <- try(suppressWarnings(AnnotationDbi::select(db.obj, keys = in.ids, 
            keytype = in.type, columns = c(in.type, out.type))))
        all.mappings <- rbind(all.mappings, res)
    }
    res <- as.data.frame(all.mappings[-1, ])
    
    res <- res[, c(in.type, out.type)]
    
    na.idx <- is.na(res[, 2])
    if (sum(na.idx) > 0) {
        n.na <- length(unique(res[na.idx, 1]))
        if (n.na > 0) {
            print(paste("Note:", n.na, "of", length(in.ids), "unique input IDs unmapped."))
        }
        if (na.rm) 
            res <- res[!na.idx, ]
    }
    
    cns <- colnames(res)
    
    res <- as.matrix(res)
    rownames(res) <- NULL
    return(res)
}

#########################################################################################################
find.pathways.by.keywords <- function(keywords, keyword.type, keywords.logic, mol.name.match, 
    SBGNview.data.folder = "./SBGNview.tmp.data/") {
    if (is.null(keywords)) {
        pathways <- pathways.info[, c("pathway.id", "pathway.name", "sub.database", 
            "database")]
    } else if (keyword.type == "pathway.name") {
        keywords <- tolower(keywords)
        pathways.info[, "pathway.name"] <- tolower(pathways.info[, "pathway.name"])
        
        
        if.selected <- grepl(pattern = keywords[1], ignore.case = TRUE, pathways.info$pathway.name)
        if (length(keywords) > 1) {
            for (i in 2:length(keywords)) {
                if (keywords.logic == "and") {
                  if.selected <- if.selected & grepl(pattern = keywords[i], ignore.case = TRUE, 
                    pathways.info$pathway.name)
                } else if (keywords.logic == "or") {
                  if.selected <- if.selected | grepl(pattern = keywords[i], ignore.case = TRUE, 
                    pathways.info$pathway.name)
                } else {
                  stop("keywords.logic must be one of 'and' or 'or'!! ")
                }
            }
        }
        pathways <- pathways.info[if.selected, c("pathway.id", "pathway.name", "sub.database", 
            "database")]
        
        
    } else {
        # using IDs to search for pathway type.pair.name =
        # paste(keyword.type,'_pathway.id',sep='')
        data(mapped.ids)
        # if (keyword.type %in% mapped.ids$gene) {
        if (keyword.type %in% unique(mapped.ids$gene[,2])) {
            cpd.or.gene <- "gene"
        } else if (keyword.type %in% mapped.ids$cpd) {
            cpd.or.gene <- "compound"
        }
        mapping.list <- loadMappingTable(input.type = keyword.type, output.type = "pathway.id", 
            species = NULL, SBGNview.data.folder = SBGNview.data.folder, cpd.or.gene = cpd.or.gene)
        mapping.table <- mapping.list[[1]][[1]]
        if (keyword.type == "CompoundName") {
            if (keywords.logic == "and") {
                keywords.logic <- "or"
                warning("Searching pathways by Compound Names. Using keywords.logic='or'. Please don't set keywords.logic to 'and'!!!\n")
            }
            keywords <- tolower(keywords)
            mapping.table[, keyword.type] <- tolower(as.character(mapping.table[, 
                keyword.type]))
            if (mol.name.match == "exact.match") {
                if.selected <- tolower(mapping.table[, keyword.type]) %in% keywords
                pathways <- mapping.table[if.selected, ]
            } else if (mol.name.match == "presence.of.input-string.in.target-name") {
                if.selected <- grepl(pattern = paste(keywords, collapse = "|"), ignore.case = TRUE, 
                  mapping.table[, keyword.type])
                sum(if.selected)
                pathways <- mapping.table[if.selected, ]
            } else if (mol.name.match == "jaccard") {
                best.match <- match.names(input.names = tolower(keywords), output.names = tolower(mapping.table[, 
                  keyword.type]))
                pathways <- merge(best.match, mapping.table, by.x = "output.name", 
                  by.y = keyword.type)
                names(pathways)[1] <- c("CompoundName")
            }
        } else {
            pathways <- mapping.table[tolower(mapping.table[, keyword.type]) %in% 
                tolower(keywords), ]
        }
    }
    return(pathways = pathways)
}

#########################################################################################################
filter.pathways.by.org <- function(pathways, org, pathway.completeness, pathways.info.file.folder = "./SBGNview.tmp.data") {
    org <- tolower(org)
    data("pathway.species.pct_Mapped")
    if (org == "all") {
        org <- unique(pathway.species.pct_Mapped$species)
    }
    pathway.species.pct_Mapped[, "species"] <- as.character(pathway.species.pct_Mapped[, 
        "species"])
    pathway.species.pct_Mapped[, "pathway"] <- as.character(pathway.species.pct_Mapped[, 
        "pathway"])
    
    if (!is.null(pathway.completeness)) {
        message("Using user provided pathway completeness cutoff: the same cutoff for different pathways!\n")
        
        pathway.ids <- pathway.species.pct_Mapped[pathway.species.pct_Mapped$species %in% 
            org & pathway.species.pct_Mapped$pct.mapped.species.pathway >= pathway.completeness, 
            ]
    } else {
        message("Using pre-generated pathway-specific completeness cutoff: different cutoff for different pathways! Pathway 'exist' and 'not exist' accross all species defined by this cutoff have the largest ANOVA F statistic when comparing completeness between 'exist' and 'not exist' groups. \n")
        pathway.ids <- pathway.species.pct_Mapped[pathway.species.pct_Mapped$species %in% 
            org, ]
        data("pathway.completeness.cutoff.info")
        pathway.specific.cutoff <- pathway.completeness.cutoff.info$cutoff
        names(pathway.specific.cutoff) <- pathway.completeness.cutoff.info$pathway
        if.pass.cutoff <- pathway.ids$pct.mapped.species.pathway > pathway.specific.cutoff[pathway.ids$pathway]
        pathway.ids <- pathway.ids[if.pass.cutoff, ]
    }
    colnames(pathway.ids)[1] <- "pathway.id"
    pathways <- merge(pathways, pathway.ids, all = FALSE)
    return(pathways)
}

#########################################################################################################
#' Retrieve pathways by keywords
#' 
#' This function searches for pathways by input keywords.
#' 
#' @param keywords A character string or vector. The search is case in-sensitive.
#' @param keywords.logic A character string. Options are 'and' or 'or'. This will tell the function if the search require 'all' or 'any' of the keywords to be present. It only makes difference when keyword.type is 'pathway.name'.
#' @param keyword.type A character string. Either 'pathway.name' or one of the ID types in \code{data('mapped.ids')}
#' @param org  A character string. The KEGG species code.
#' @param SBGNview.data.folder A character string. 
#' @details If 'keyword.type' is 'pathway.name' (default), this function will search for the presence of any keyword in the pathway.name column of data(pathways.info). The search is case in-sensitive. If 'keyword.type' is one of the identifier types and 'keywords' are corresponding identifiers, this function will return pathways that include nodes mapped to input identifiers. 
#' @return A dataframe. Contains information of pathways found.
#' @examples 
#' data(pathways.info)
#' input.pathways <- findPathways('Adrenaline and noradrenaline biosynthesis')
#' @export

findPathways <- function(keywords = NULL, keywords.logic = "or", keyword.type = "pathway.name", 
    org = NULL,   SBGNview.data.folder = "./SBGNview.tmp.data/") {
    pathway.completeness = NULL
    mol.name.match = c("exact.match", "jaccard", 
        "presence.of.input-string.in.target-name")[1]
    if (!file.exists(SBGNview.data.folder)) {
        dir.create(SBGNview.data.folder)
    }
    
    keywords <- as.vector(keywords)
    
    pathways <- find.pathways.by.keywords(keywords, keyword.type, keywords.logic, 
        mol.name.match, SBGNview.data.folder = SBGNview.data.folder)
    
    if (!is.null(org)) {
        pathways <- filter.pathways.by.org(pathways, org, pathway.completeness, pathways.info.file.folder = "./SBGNview.tmp.data")
    }
    pathways.with.name <- pathways.info[, c("pathway.id", "pathway.name")]
    pathways <- merge(pathways, pathways.with.name, all.x = TRUE)
    pathways <- unique(pathways)
    
    return(pathways)
}

#########################################################################################################
# given two vectors of characters/words, this function find the best match
# between words in the two vectors, by spliting a word into sub-strings and use
# 'jaccard' to calculate two words' similarity (similarity between thier
# sub-string vectors). It for one name, the function will output the pair with
# largest similarity score.
match.names <- function(input.names, output.names, output.file = NULL) {
    input.names <- tolower(input.names)
    output.names <- tolower(output.names)
    input.names <- unique(input.names)
    output.names <- unique(output.names)
    
    all.pairs <- expand.grid(input.names, output.names)
    
    # calculate jaccard of all possible matches
    ot <- apply(all.pairs, 1, function(pair) {
        p1n <- pair[1]
        p2n <- pair[2]
        p1 <- strsplit(p1n, "[^a-zA-Z0-9]+", perl = TRUE)[[1]]
        p1 <- unique(p1)
        p1
        p2 <- strsplit(p2n, "[^a-zA-Z0-9]+", perl = TRUE)[[1]]
        p2 <- unique(p2)
        p2
        inter <- length(intersect(p1, p2))
        uni <- length(unique(union(p1, p2)))
        jc <- inter/uni
        return(c(p1n, p2n, jc))
    })
    ot <- t(ot)
    
    ot <- as.data.frame(ot, stringsAsFactors = FALSE)
    names(ot) <- c("input.name", "output.name", "jaccard")
    ot$jaccard <- as.numeric(ot$jaccard)
    
    # extract the pairs with max jaccard
    ot1 <- ot
    mapped.table <- by(ot1, as.factor(ot1$input.name), function(df) {
        if.max <- df$jaccard == max(df$jaccard)
        mx.df <- df[if.max, ]
    })
    mapped.table <- do.call(rbind, mapped.table)
    row.names(mapped.table) <- c()
    if (!is.null(output.file)) {
        write.table(mapped.table, output.file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
    return(mapped.table)
}

#########################################################################################################
#' Retrieve gene list or compound list from collected databases
#' 
#' @param database Character string. The database where gene list will be extracted. Acceptable values: 'MetaCyc', 'pathwayCommons', 'MetaCrop'. The value is case in-sensitive.
#' @param mol.list.ID.type Character string. The ID type of output gene list. One of the supported types in \code{data('mapped.ids')}
#' @param org Character string. The three letter species code used by KEGG. E.g. 'hsa','mmu'
#' @param cpd.or.gene Character string. One of 'gene' or 'compound'
#' @param output.pathway.name Logical. If set to 'TRUE', the names of returned list are in the format: 'pathway.id::pathway.name'. If set to 'FALSE', the format is 'pahtway.id'
#' @param combine.duplicated.set Logical.  Some pathways have the same geneset. If this parameter is set to 'TRUE', the output list will combine pathways that have the same gene set. The name in the list will be pathway names concatinated with '||'
#' @param truncate.name.length Integer. The pathway names will be truncated to at most that length. 
#' @param SBGNview.data.folder A character string.
#' @return A list. Each element is a genelist of a pathway.
#' @examples 
#' data(pathways.info)
#' mol.list <- getMolList(
#'                  database = 'pathwayCommons',
#'                  mol.list.ID.type = 'ENTREZID',
#'                  org = 'hsa'
#' )
#'   
#' @export

getMolList <- function(database = "pathwayCommons", mol.list.ID.type = "ENTREZID", 
    org = "hsa", cpd.or.gene = "gene", output.pathway.name = TRUE, combine.duplicated.set = TRUE, 
    truncate.name.length = 50, SBGNview.data.folder = "./SBGNview.tmp.data") {
    if (tolower(database) == "metacrop") {
        if (cpd.or.gene == "gene") {
            id.in.pathway <- "ENZYME"
        } else {
            id.in.pathway <- "CompoundName"
        }
        output.pathways <- subset(pathways.info, database == "MetaCrop", select = "pathway.id")
    } else if (tolower(database) %in% c("pathwaycommons", "metacyc")) {
        if (cpd.or.gene == "gene") {
            id.in.pathway <- "KO"
        } else {
            id.in.pathway <- "chebi"
        }
        output.pathways <- subset(pathways.info, database != "MetaCrop", select = "pathway.id")
    }
    # metacrop initial list is using enzyme
    mapping.list <- loadMappingTable(input.type = id.in.pathway, output.type = "pathway.id", 
        cpd.or.gene = cpd.or.gene, species = org, SBGNview.data.folder = SBGNview.data.folder)
    ref.to.pathway <- mapping.list[[1]][[1]]
    
    
    if (mol.list.ID.type == id.in.pathway) {
        out.id.to.pathway <- ref.to.pathway
        # change KO to output id
    } else {
        message(id.in.pathway, mol.list.ID.type, "\n\n")
        out.id.type.to.ref <- loadMappingTable(input.type = id.in.pathway, output.type = mol.list.ID.type, 
            cpd.or.gene = cpd.or.gene, limit.to.ids = ref.to.pathway[, id.in.pathway], 
            species = org, SBGNview.data.folder = SBGNview.data.folder)
        out.id.type.to.ref <- out.id.type.to.ref[[1]][[1]]
        # merge KO to pathway and KO to output id
        out.id.to.pathway <- merge(out.id.type.to.ref, ref.to.pathway, all.x = TRUE)
        out.id.to.pathway <- out.id.to.pathway[, c(mol.list.ID.type, "pathway.id")]
        out.id.to.pathway <- unique(out.id.to.pathway)
        out.id.to.pathway <- out.id.to.pathway[!is.na(out.id.to.pathway[, 2]), ]
    }
    out.id.to.pathway <- unique(out.id.to.pathway)
    out.id.to.pathway <- out.id.to.pathway[out.id.to.pathway[, "pathway.id"] %in% 
        output.pathways$pathway.id, ]
    out.id.to.pathway <- out.id.to.pathway[out.id.to.pathway[, mol.list.ID.type] != 
        "", ]
    out.id.to.pathway <- split(as.character(out.id.to.pathway[, mol.list.ID.type]), 
        out.id.to.pathway[, "pathway.id"])
    out.id.to.pathway <- out.id.to.pathway[!is.na(names(out.id.to.pathway))]
    
    if (output.pathway.name) {
        # merge pathway names
        pathway.id.to.name <- pathways.info[, c("pathway.id", "pathway.name")]
        row.names(pathway.id.to.name) <- pathway.id.to.name[, "pathway.id"]
        pathway.ids <- paste(pathway.id.to.name[, "pathway.id"], pathway.id.to.name[, 
            "pathway.name"], sep = "::")
        names(pathway.ids) <- pathway.id.to.name[, "pathway.id"]
        names(out.id.to.pathway) <- pathway.ids[names(out.id.to.pathway)]
        
    }
    
    if (combine.duplicated.set) {
        sets <- lapply(out.id.to.pathway, function(ids) {
            ids.str <- paste(sort(ids), collapse = "||")
        })
        sets <- unlist(sets)
        pathways.same.set <- tapply(names(sets), as.factor(sets), function(pathways) {
            pathways.joint <- paste(sort(pathways), collapse = "||")
        })
        out.id.to.pathway <- as.list(names(pathways.same.set))
        names(out.id.to.pathway) <- pathways.same.set
        out.id.to.pathway <- lapply(out.id.to.pathway, function(ids) {
            strsplit(ids, "\\|\\|")[[1]]
        })
    }
    names(out.id.to.pathway) <- substr(names(out.id.to.pathway), 1, truncate.name.length)
    
    sorted.names <- sort(names(out.id.to.pathway), method = "radix", decreasing = FALSE)
    out.id.to.pathway <- out.id.to.pathway[sorted.names]
    return(out.id.to.pathway)
}

#########################################################################################################
#' Change input IDs to another ID type
#' 
#' @param input.ids A vector of character strings. Input IDs that need to be converted.
#' @param input.type A character string. The type of input IDs. Supported ID types can be found in data(mapped.ids)
#' @param output.type A character string. The type of output IDs. Supported ID types can be found in data(mapped.ids)
#' @param cpd.or.gene A character string. One of 'gene' or 'compound'. 
#' @param limit.to.pathways A vector of character strings. A vector of pathways IDs. When 'output.type' is one of 'pathwayCommons' or 'metacyc.SBGN', one input ID (e.g. gene symbol) can map to multiple nodes in different pathways (SBGN-ML files). In this case, we can limit output to the specified pathways. When 'output.type' is NOT one of 'pathwayCommons' or 'metacyc.SBGN', this argument is ignored.
#' @param org Character string. Three letter KEGG species code.
#' @param SBGNview.data.folder A character string. The path to a folder that will hold download ID mapping files and pathway information data files. 
#' @return A list. Each element is the IDs in output.type that are mapped to one input ID.
#' 
#' @examples 
#'  data(mapped.ids)
#' mapping = changeIds(
#'   input.ids = c(100048912),
#'   input.type = 'ENTREZID',
#'   output.type = 'pathwayCommons',
#'   cpd.or.gene = 'gene',
#' )
#' 
#' @export

changeIds <- function(input.ids, input.type, output.type, cpd.or.gene, limit.to.pathways = NULL, 
    org = "hsa", SBGNview.data.folder = "./SBGNview.tmp.data") {
    if (!is.null(limit.to.pathways) & output.type %in% c("pathwayCommons", "metacyc.SBGN")) {
        ids.in.pathways <- sbgnNodes(limit.to.pathways, SBGNview.data.folder = SBGNview.data.folder)
        limit.to.output.ids <- unlist(lapply(ids.in.pathways, function(all.nodes) {
            # pathway.nodes$all.nodes[,'glyph.id']
            all.nodes[, "glyph.id"]
        }))
        limit.to.output.ids
    } else {
        limit.to.output.ids <- NULL
    }
    if (cpd.or.gene == "gene") {
        id.mapping.all.list <- load.id.mapping.list.all(SBGN.file.gene.id.type = output.type, 
            output.gene.id.type = input.type, species = org, SBGNview.data.folder = SBGNview.data.folder)
    } else if (cpd.or.gene == "compound") {
        id.mapping.all.list <- load.id.mapping.list.all(SBGN.file.cpd.id.type = output.type, 
            output.cpd.id.type = input.type, species = org, SBGNview.data.folder = SBGNview.data.folder)
    } else {
        stop("cpd or gene much be one of 'gene' or 'compound'!!")
    }
    new.ids <- sapply(input.ids, function(x) {
        mapped.ids <- change.id(input.id = x, cpd.or.gene = cpd.or.gene, input.type = input.type, 
            output.type = output.type, id.mapping.all.list = id.mapping.all.list)
        mapped.ids
        mapped.ids <- strsplit(mapped.ids, "; ")[[1]]  # one input id can map to multiple glyph ids(from the same pathway or different pathways), we select the ones that are in the glyphs (same pathway).
        if (!is.null(limit.to.output.ids)) {
            mapped.ids <- intersect(mapped.ids, limit.to.output.ids)
        }
        return(mapped.ids)
    })
    if (is.matrix(new.ids)) {
        new.ids <- as.list(as.data.frame(new.ids, stringsAsFactors = FALSE))
        new.ids
    }
    return(new.ids)
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
break.text.into.segments <- function(label, w, glyph.class, global.parameters.list, 
    max.x = 0, glyph) {
    if.long.word <- FALSE
    text.length.factor <- global.parameters.list$text.length.factor.macromolecule
    if (length(glyph@text$font.size) > 0) {
        font.size <- glyph@text$font.size
    } else {
        font.size <- global.parameters.list$font.size
        
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
            if (global.parameters.list$if.scale.compartment.font.size) {
                font.size <- max(glyph@shape$stroke.width * 3.5, font.size * global.parameters.list$node.width.adjust.factor.compartment * 
                  w)
            } else {
                font.size <- font.size * global.parameters.list$node.width.adjust.factor
                font.size <- font.size * global.parameters.list$font.size.scale.gene * 
                  global.parameters.list$font.size.scale.compartment
            }
            if (max.x > 0) {
                # font.size = font.size * max(1,max.x/2000)
            }
            text.length.factor <- global.parameters.list$text.length.factor.compartment
        } else if (glyph.class == "complex") {
            if (global.parameters.list$if.scale.complex.font.size) {
                font.size <- font.size * global.parameters.list$node.width.adjust.factor.complex * 
                  w
            } else {
                font.size <- font.size * global.parameters.list$node.width.adjust.factor
                font.size <- font.size * global.parameters.list$font.size.scale.gene * 
                  global.parameters.list$font.size.scale.complex
            }
            text.length.factor <- global.parameters.list$text.length.factor.complex
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
            font.size <- font.size * global.parameters.list$status.node.font.scale
        } else {
            font.size <- font.size * 2 * glyph@h * global.parameters.list$node.width.adjust.factor/70
            font.size <- font.size * global.parameters.list$font.size.scale.gene
            
            
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
            
            if (current.line.to.be.length > text.length.factor * w & (length(label.words) - 
                i) > 2 & (word.previous %in% global.parameters.list$label.spliting.string | 
                identical(global.parameters.list$label.spliting.string, c("any")))) {
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
    return(list(words.segments = words.segments, label.margin = label.margin, font.size = font.size, 
        if.long.word = if.long.word, nline = nline))
}

#########################################################################################################
