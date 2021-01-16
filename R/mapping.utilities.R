
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
# this function is a slightly modified version of mol.sum from pathview (1.30.1/1.31.1). This is a transitional
# copy and will be merge into pathview for the future. This function replace the old mol.sum.multiple.mapping
# + merge.molecule.data functions, which are very slow, and affect SBGNview() and many higher level functions.
mol.sum.multiple.mapping <-function(mol.data, id.map, gene.annotpkg="org.Hs.eg.db", sum.method=c("sum","mean", "median", "max", "max.abs", "random")[1]){
  if(is.character(mol.data)){
    gd.names=mol.data
    mol.data=rep(1, length(mol.data))
    names(mol.data)=gd.names
    ng=length(mol.data)
  } else if(!is.null(mol.data)){
    if(length(dim(mol.data))==2){
    gd.names=rownames(mol.data)
    ng=nrow(mol.data)
    } else if(is.numeric(mol.data) & is.null(dim(mol.data))){
    gd.names=names(mol.data)
    ng=length(mol.data)
    } else stop("wrong mol.data format!")
  } else stop("NULL mol.data!")

  if(is.character(id.map) & length(id.map)==1){
    id.map=id2eg(gd.names, category=id.map, pkg.name=gene.annotpkg)
  }
  
  sel.idx=id.map[,2]>"" & !is.na(id.map[,2])
  id.map=id.map[sel.idx,]
#  eff.idx=gd.names %in% id.map[,1]
  eff.idx1=id.map[,1] %in% gd.names
  id.map1=id.map[eff.idx1,]

#  map.idx=match(id.map[,1], gd.names[eff.idx])
#  mapped.ids=id.map[match(gd.names[eff.idx], id.map[,1]),2]
#key difference between mol.sum current and old versions is indexing
#mol.data using id.map1[,1] vs eff.idx, i.e. IDs vs indicators T/F
#  mol.data=cbind(cbind(mol.data)[id.map1[,1],])
  nmapped=sum(eff.idx1) #sum(eff.idx)
  if(nmapped<1) stop("no ID can be mapped!")
  else if(nmapped==1){
#    mapped.data=rbind(cbind(mol.data)[eff.idx,])
#    rownames(mapped.data)=mapped.ids[1]
    mapped.data=rbind(cbind(mol.data)[id.map1[1],])
    rownames(mapped.data)=id.map1[2]
  }
  else{
      mapped.ids=id.map1[,2]
      if(sum.method %in% c("sum","mean")){
      sum.method=eval(as.name(sum.method))
#      mapped.data=apply(cbind(cbind(mol.data)[eff.idx,]),2,function(x){
      mapped.data=apply(cbind(cbind(mol.data)[id.map1[,1],]),2,function(x){
        sum.res=tapply(x, mapped.ids, sum.method, na.rm=T)
        return(sum.res)
      })
#      if(length(unique(mapped.ids))==1){
#        if(length(mapped.data)>1){
#          mapped.data=rbind(mapped.data)
#          rownames(mapped.data)=mapped.ids[1]
#        }
#        else names(mapped.data)=mapped.ids[1]
#      }
    } else{
    sum.method=eval(as.name(sum.method))
#    mol.data=cbind(cbind(mol.data)[eff.idx,])
#    mol.data=cbind(cbind(mol.data)[id.map1[,1],])
    if(all(mol.data>=0) | all(mol.data<=0)){
      vars=apply(cbind(mol.data), 1, IQR)
    } else vars=apply(cbind(mol.data), 1, sum, na.rm=T)
    
#    sel.rn=tapply(1:sum(eff.idx), mapped.ids, function(x){
#    sel.rn=tapply(which(eff.idx), mapped.ids, function(x){
    sel.rn=tapply(1:length(mapped.ids), mapped.ids, function(x){
      if(length(x)==1) return(x)
      else return(x[which.min(abs(vars[x]-sum.method(vars[x], na.rm=T)))])
    })
    if(length(sel.rn)>1) mapped.data=cbind(mol.data[sel.rn,])
    else mapped.data=rbind(mol.data[sel.rn,])
    rownames(mapped.data)=names(sel.rn)
  }
}
  return(mapped.data)

}

#########################################################################################################
# the first column of id.map is input ids(need to convert FROM), second column is output ids(converting TO)
mol.sum.multiple.mapping.old <- function(mol.data, id.map, sum.method, input.sbgn.file = "na"){   
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
merge.molecule.data <- function(in.data.source.id, id.map.in.target.and.source, sum.method){
  
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
#' Change vector of input IDs to another ID type
#' 
#' This function changes a vector of input IDs to another ID type. It returns a list where each element is the input ID type is mapped to the output ID type. Please note that not every input ID might be mapped to the output type. To change the input IDs to another ID type for a data matrix, use the \code{\link{changeDataId}} function. 
#' 
#' @param input.ids A vector of character strings. Input IDs that need to be converted.
#' @param input.type A character string. The type of input IDs. Supported ID types can be found in data(mapped.ids)
#' @param output.type A character string. The type of output IDs. Supported ID types can be found in data(mapped.ids)
#' @param cpd.or.gene A character string. One of 'gene' or 'compound'. 
#' @param limit.to.pathways A vector of character strings. A vector of pathways IDs. When 'output.type' is one of 'pathwayCommons' or 'metacyc.SBGN', one input ID (e.g. gene symbol) can map to multiple nodes in different pathways (SBGN-ML files). In this case, we can limit output to the specified pathways. When 'output.type' is NOT one of 'pathwayCommons' or 'metacyc.SBGN', this argument is ignored.
#' @param org A character string. Default: "hsa". Three letter KEGG species code.
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files. 
#' @return A list. Each element is the IDs in output.type that are mapped to one input ID.
#' 
#' @examples 
#' data(mapped.ids)
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
    # id.mapping.all.list <- load.id.mapping.list.all(SBGN.file.gene.id.type = output.type, 
    #                                                 output.gene.id.type = input.type, 
    #                                                 species = org, SBGNview.data.folder = SBGNview.data.folder)
    #### use loadMappingTable
    id.mapping.all.list <- loadMappingTable(input.type = input.type, output.type = output.type,
                                            species = org, cpd.or.gene = "gene", limit.to.ids = input.ids,
                                            SBGNview.data.folder = SBGNview.data.folder)
    
  } else if (cpd.or.gene %in% c("cpd", "compound")) {
    # id.mapping.all.list <- load.id.mapping.list.all(SBGN.file.cpd.id.type = output.type, 
    #                                                 output.cpd.id.type = input.type, 
    #                                                 species = org, SBGNview.data.folder = SBGNview.data.folder)
    ### use loadMappingTable
    id.mapping.all.list <- loadMappingTable(input.type = input.type, output.type = output.type,
                                            cpd.or.gene = "compound", limit.to.ids = input.ids,
                                            SBGNview.data.folder = SBGNview.data.folder)
    
  } else {
    stop("cpd.or.gene must be one of 'gene' or 'compound'!!")
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
  
  message("\nChanged IDs from ", input.type, " to ", output.type)
  
  ###### checking how many IDs were mapped
  not.mapped.count <- 0
  mapped.count <- 0
  for(idx in seq_along(new.ids)){
    if(length(new.ids[[idx]]) == 0) { not.mapped.count <- not.mapped.count + 1 
    } else { mapped.count <- mapped.count + 1 }
  }
  if(not.mapped.count == length(new.ids)){
    message("None of the input IDs were mapped to specified ouput type")
  } else if (mapped.count == length(new.ids)) {
    message("All input IDs were mapped")
  } else {
    message("**NOTE**: ", mapped.count, " of ", length(new.ids), " in 'input.ids' were mapped to the output type ")
    message("Please check the ouput list for more information")
  }
  ######
  
  return(new.ids)
}

#########################################################################################################
#' Change the data IDs of input omics data matrix
#' 
#' This function changes the IDs of input omics data from one type to another. It returns a data matrix with row names changed to the specified output ID type. To change a vector of input IDs to another type, use the \code{\link{changeIds}} function.
#' 
#' @param data.input.id A matrix. Input omics data. Rows are genes or compounds, columns are measurements. Row names are the original IDs that need to be transformed.
#' @param input.type A character string. The type of input IDs. Please check \code{data('mapped.ids')} for supported types.
#' @param output.type A character string. The type of output IDs. Please check \code{data('mapped.ids')} for supported types. 
#' @param sum.method  A character string. Default: "sum". In some cases multiple input IDs are mapped to one output ID. In this situation ,we may need to derive only one value from them. This parameter is a function that can derive a single numeric value from a vector of numeric values (e.g. 'sum','max','min','mean'), including a User Defined Function (UDF).
#' @param org  A character string. Default: "hsa". The species source of omics data. 'changeDataId' uses pathview to map between some gene ID types. Please use '?geneannot.map' to check the detail. Pathview needs species information to do the job. This parameter is a two-letter abbreviation of organism name, or KEGG species code, or the common species name, used to determine the gene annotation package. For all potential values check: data(bods); bods. Default org='Hs', and can also be 'hsa' or 'human' (case insensitive). 
#' @param cpd.or.gene  A character string. Either 'compound' or 'gene' -- the type of input omics data. 
#' @param id.mapping.table A matrix.  Mapping table between input.type and output.type. This matrix should have two columns for input.type and output.type, respectively.  Column names should be the values of parameters 'input.type' and 'output.type'. See example section for an example. 
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files. 
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
                         org = "hsa", cpd.or.gene, id.mapping.table = NULL, 
                         SBGNview.data.folder = "./SBGNview.tmp.data") {
  # function: change input data matrix with input/uniprot ID to data matrix with
  # pathwaycommons id if there are multiple input IDs mapped to a node, collapse
  # their values using user specified function(sum.method)
  input.type <- gsub("entrez", "ENTREZID", input.type)
  output.type <- gsub("entrez", "ENTREZID", output.type)
  if (is.null(id.mapping.table)) {
    # if user didn't provide mapping table, we try to download one.
    # mapping.list <- loadMappingTable(output.type = output.type, input.type = input.type, 
    #                                  species = org, cpd.or.gene = cpd.or.gene, limit.to.ids = row.names(data.input.id), 
    #                                  SBGNview.data.folder = SBGNview.data.folder)
    
    # dim(obj) = NULL then vector (names), else the matrix or data.frame (row.names)
    input.ids <- NULL
    if(is.null(dim(data.input.id))){
      input.ids <- names(data.input.id)
    } else {
      input.ids <- row.names(data.input.id)
    }
    mapping.list <- loadMappingTable(output.type = output.type, input.type = input.type, 
                                     species = org, cpd.or.gene = cpd.or.gene, 
                                     limit.to.ids = input.ids, 
                                     SBGNview.data.folder = SBGNview.data.folder)
    
    id.map <- mapping.list[[1]][[1]]
    #id.map <- as.matrix(id.map[, c(input.type, output.type)])
  } else {
    if (!all(c(input.type) %in% colnames(id.mapping.table))) {
      message(input.type, " must be in column names of id.mapping.table!\n")
      message("Column names of id.mapping.table are:\n")
      cat(colnames(id.mapping.table), "\n\n\n\n")
      stop()
    }
    if (!all(c(output.type) %in% colnames(id.mapping.table))) {
      message(output.type, " must be in column names of id.mapping.table!\n")
      message("Column names of id.mapping.table are:\n")
      cat(colnames(id.mapping.table), "\n\n\n\n")
      stop()
    }
    id.map <- id.mapping.table
  }
  message("Changing data IDs")
  in.data.target.id <- mol.sum.multiple.mapping(mol.data = data.input.id, id.map = id.map, 
                                                sum.method = sum.method)
#  in.data.target.id <- pathview::mol.sum(mol.data = data.input.id, id.map = id.map,
#                                                sum.method = sum.method)
  message("Finished changing data IDs")
  return(in.data.target.id)

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
