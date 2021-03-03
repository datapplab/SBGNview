
#########################################################################################################
# searches for different id mapping file name combinations in local, SBGNview.data, SBGNhub locations
# check local directory for any mapping files, if not found then check SBGNView.data, then SBGNhub
# returns list containing the mapping file name and the location where its found
download.mapping.file <- function(input.type, output.type, species = NULL, 
                                  SBGNview.data.folder = "./SBGNview.tmp.data") { 
  
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
  try.file.names <- gsub("[[:space:]]", "", try.file.names) # remove any spaces
  
  # KO_ENTREZID.RData in SBGNhub may not contain the all the information for the speices in it. 
  # For case when input/output is KO/ENTREZID, try.file.names contains KO_ENTREZID,
  # so if species specific file not found, KO_ENTREZID will be downloaded and used instead 
  # of generating the species specific mapping between KO and ENTREZID. Modify 'try.files.names'
  # to search only for species specific mapping b/w KO and ENTREZID so if file not found, 
  # mapping is generated from scratch using KEGGREST.
  if(!is.null(species) & any(c(input.type, output.type) %in% c("KO", "ko")) &
     any(c(input.type, output.type) %in% c("ENTREZID"))) {
    try.file.names <- c(type.pair.name.1.org, type.pair.name.2.org)
  }
  
  location <- ""
  mapping.file.name <- ""
  
  for(try.file.name in try.file.names) { # for each try.file.name -> check local, data, hub if not exists()
    
    # try.file.name already exists in environment
    if(exists(try.file.name)) {
      location <- "exists"
      mapping.file.name <- try.file.name
      break
      
    } else { # check local, SBGNview.data, or hub
      
      file.name.rdata <- paste(try.file.name, ".RData", sep = "")
      # check local
      if(file.exists(paste(SBGNview.data.folder, file.name.rdata, sep = "/"))) {
        location <- "local"
        mapping.file.name <- file.name.rdata
        message("\nLocal mapping table file found: ", mapping.file.name)
        (break)()
      }
      
      # check SBGNview.data
      if.exists.SBGNview.data <- tryCatch(data(list = try.file.name), warning = function(w) "none")
      if(if.exists.SBGNview.data != "none") { # found in SBGNview.data
        location <- "SBGNview.data"
        mapping.file.name <- try.file.name
        message("\n", if.exists.SBGNview.data, " mapping table loaded from SBGNview.data")
        break
      }
      
      # check SBGNhub for .RData file
      if(location == "") {
        mapping.file.name <- search.sbgnhub.id.mapping(try.file.name = file.name.rdata, SBGNview.data.folder)
        if(!is.null(mapping.file.name)) { # found in SBGNhub
          location <- "SBGNhub"
          message(mapping.file.name, " downloaded from SBGNhub")
          message("saved to ", SBGNview.data.folder)
          break
        } 
      }
      
    } # end else block
  } # end for loop
  
  if(location == "" & is.null(mapping.file.name)) { # mapping table not found in local, SBGNview.data, SBGNhub
    location <- "pathview"
    mapping.file.name <- ""
  }
  
  mapping.file.info <- list(location = location, mapping.file.name = mapping.file.name)
  return(mapping.file.info)
}

#########################################################################################################
# used by download.mapping.file function
# search SBGNhub.id.mapping.tables, if file exits download from SBGNhub
search.sbgnhub.id.mapping <- function(try.file.name, SBGNview.data.folder) {
  
  mapping.file.name <- NULL
  if(!exists("SBGNhub.id.mapping.tables")) data(SBGNhub.id.mapping.tables, package="SBGNview")
  if(try.file.name %in% SBGNhub.id.mapping.tables) {
    message("\n", try.file.name, " found in SBGNhub")
    file.url <- paste("https://raw.githubusercontent.com/datapplab/SBGNhub/master/data/id.mapping.unique.pair.name/", 
                      try.file.name, sep = "")
    file.destination <- paste(SBGNview.data.folder, try.file.name, sep = "/")
    download.file(file.url, destfile = file.destination, method = "auto", mode = "wb")
    mapping.file.name <- try.file.name
  }
  return(mapping.file.name)
}

#########################################################################################################
#' Generate ID mapping table from input and output ID types 
#' 
#' This function generates the ID mapping table between input and output ID types. If provided a vector of input IDs (limit.to.ids argument), the function will output mapping table only containing the input IDs. Otherwise, the function will output all IDs of input and output types (restricted to a species if working on gene types and specified  the 'species' parameter).
#' This function will check if a mapping table exists in local folder ('SBGNview.data.folder' argument), SBGNview.data package, and SBGNhub. If no mapping table is found in these locations, a mapping table is generated using Pathview and saved to path specified by 'SBGNview.data.folder' argument.
#' 
#' @param input.type A character string. Gene or compound ID type
#' @param output.type A character string. Gene or compound ID type
#' @param species A character string. Three letter KEGG species code.
#' @param mol.type A character string. Either 'gene' or 'cpd'.
#' @param limit.to.ids Vector. Molecule IDs of 'input.type'.
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files.
#' @return A list containing the mapping table. 
#' @examples 
#'  data(mapped.ids)
#'  entrez.to.pathwayCommons <- loadMappingTable(input.type = 'ENTREZID',
#'                                               output.type = 'pathwayCommons',
#'                                               species = 'hsa',
#'                                               mol.type = 'gene')
#'                              
#' @export

loadMappingTable <- function(input.type = NULL, output.type = NULL, species = NULL, mol.type = NULL, 
                             limit.to.ids = NULL, SBGNview.data.folder = "./SBGNview.tmp.data") { 
  
  # input validation 
  if(is.null(input.type) | is.null(output.type)) {
    stop("Please make sure both input.type and output.type arguments are specified.") 
  }
  if(input.type %in% c("entrez", "eg", "entrezid", "ENTREZ")) input.type = "ENTREZID"
  if(output.type %in% c("entrez", "eg", "entrezid", "ENTREZ")) output.type = "ENTREZID"
  # changed mapping file names from CompondName to compound.name and kegg.ligand to kegg
  # we handle if input/output is CompoundName/kegg.ligand/KEGG
  if(input.type == "CompoundName") input.type <- "compound.name"
  if(input.type %in% c("kegg.ligand", "KEGG")) input.type <- "kegg"
  if(output.type == "CompoundName") output.type <- "compound.name"
  if(output.type %in% c("kegg.ligand", "KEGG")) output.type <- "kegg"
  
  if(input.type == output.type) {
    stop("Input type and output types are the same!")
  }
  species = gsub("Hs","hsa",species)
  
  if(is.null(mol.type) | !mol.type %in% c("gene", "cpd")) { 
    stop("'mol.type' argument must be 'gene' or 'cpd'")
  }
  
  if(!file.exists(SBGNview.data.folder)) { 
    dir.create(SBGNview.data.folder) 
  }
  
  # 'download.mapping.file checks' local, SBGNview.data, SBGNhub for mapping file
  # location values ("local", "SBGNview.data", "SBGNhub", "pathview")
  # local and SBGNhub run location, file names have .RData extension
  mapping.file.info <- download.mapping.file(input.type = input.type,
                                             output.type = output.type,
                                             species = species,
                                             SBGNview.data.folder = SBGNview.data.folder)
  
  # mapping table object already exists in environment
  if(mapping.file.info$location == "exists") { 
    message("\n", mapping.file.info$mapping.file.name, " mapping table exists in environment")
    mapping.list <- get(mapping.file.info$mapping.file.name)
    
  } else if(mapping.file.info$location == "local" || mapping.file.info$location == "SBGNhub") {
    
    # read local file or file downloaded from SBGNhub
    path.to.local.file <- paste(SBGNview.data.folder, mapping.file.info$mapping.file.name, sep = "/")
    var.name <- load(file = path.to.local.file)
    mapping.list <- get(var.name) # mapping.list returned
    message(mapping.file.info$mapping.file.name, " loaded from local")
    
  } else if (mapping.file.info$location == "SBGNview.data") {
    mapping.list <- get(mapping.file.info$mapping.file.name) # get data loaded from SBGNview.data
    
  } else {  # use pathview
    message("\nGenerating mapping using Pathview")
    if(mol.type == "gene") {
      # pathway.id used for heterogeneous purposes and mapping gene to pathway.
      # pathwayCommons and metacyc.SBGN are databases mentioned in pathways.stats
      # unique(pathways.info$macromolecule.ID.type) returns "pathwayCommons" "metacyc.SBGN"  "ENZYME"
      # ENZYME not in if condition below b/c it can be mapped directly
      if(any(c(input.type,output.type) %in% c("pathwayCommons","metacyc.SBGN","pathway.id")) &
         !any(c(input.type,output.type) %in% c("KO")) ) { # if input and output don't contain KO
        
        # load mapping table for mapping KO to glyph id
        ko.to.glyph.id <- loadMappingTable(input.type = output.type,
                                           output.type = "KO",
                                           mol.type = "gene",
                                           species = species,
                                           SBGNview.data.folder = SBGNview.data.folder)
        # load mapping table for mapping user input to KO id
        input.to.ko <- loadMappingTable(input.type = input.type,
                                        output.type = "KO",
                                        mol.type = "gene",
                                        limit.to.ids = limit.to.ids,
                                        species = species,
                                        SBGNview.data.folder = SBGNview.data.folder)
        input.to.glyph.id <- merge(input.to.ko, ko.to.glyph.id, all = FALSE)
        mapping.list <- input.to.glyph.id[,c(input.type,output.type)]
        
      } else {
        message("\nID mapping not pre-generated")
        # if(is.null(limit.to.ids)){
        #   stop("\nMust provide input IDs to 'limit.to.ids' argument when using pathview mapping!")
        # }
        mapping.list <- geneannot.map.ko(in.ids = limit.to.ids,
                                         in.type = input.type,
                                         out.type = output.type,
                                         species = species,
                                         unique.map= FALSE,
                                         SBGNview.data.folder = SBGNview.data.folder)
        if(is.vector(mapping.list)) {
          mapping.list <- as.matrix(t(mapping.list))
        }
      }
      
    } else if(mol.type == "cpd") {
      mapping.list <- generate.cpd.mapping.table(input.type = input.type, output.type = output.type,
                                                 limit.to.ids = limit.to.ids)
    } 
    
  } # end else use pathview
  
  # if mapping.list contains a species column, filter by value of 'species' argument
  if(!identical(species, character(0))) {
    if("species" %in% colnames(mapping.list) & !is.null(mapping.list)) { 
      
      message("Filtering mapping list by species: '", species, "'")
      if(!exists("korg")) data(korg, package="pathview")
      # species value is KEGG code based on 'species' argument documentation
      if(tolower(species) %in% korg[,3]) { 
        message("Species kegg code in 'korg' dataset")
        ridx <- match(species, korg[, 3]) %% nrow(korg)
        species.sci.name <- tolower(korg[ridx, 4])   # get scientific.name for input species
        
        if(any(mapping.list[,"species"] %in% species) ) { # check kegg code in species column
          mapping.list <- mapping.list[mapping.list[,"species"] %in% species, c(input.type,output.type)]
        } else if (any(tolower(mapping.list[,"species"]) %in% species.sci.name) ) {
          mapping.list <- mapping.list[tolower(mapping.list[,"species"]) %in% species.sci.name, c(input.type,output.type)]
          message("Filtered mapping list by '", species.sci.name, "'")
        } else {
          message("Was not able to filter mapping.list by ", species,  " or ", species.sci.name, ".")
        }
      } else {
        stop(species, " kegg code not in 'korg' dataset or incorrect kegg code.")
      }
      
    } # end if: filter 'species' column
  } # end if !is.null(species)
  
  if(!is.matrix(mapping.list)) mapping.list <- as.matrix(mapping.list)
  return(mapping.list)
}

#########################################################################################################
# uses generate.ko.mapping.list to generate mapping table from scratch using KEGGREST
# otherwise will use pathview
geneannot.map.ko <- function(in.ids = NULL, in.type, out.type, species = "hsa", unique.map, 
                             SBGNview.data.folder = "./SBGNview.tmp.data") {
  # pathview's geneannot.map can't map KO, so here included KO mapping
  if (!is.null(in.ids)) {
    in.ids <- gsub("\"", "", in.ids)
  }
  
  message("Mapping: ", in.type, ", ", out.type, "\n")
  
  if (any(c(in.type, out.type) %in% c("KO", "ko"))) {
    # if input/output KO/ENTREZID combo, set out.type to KO and in.type to ENTREZ
    # otherwise in.type and out.type order matter. if in.type = KO, generate.ko.mapping.list calls mapping.ko.to.arbitrary.id.type
    if(any(c(in.type, out.type) %in% c("ENTREZID"))) {
      in.type <- "ENTREZID"
      out.type <- "KO" 
    }
    
    message("Generating mapping list\n")
    id.map <- generate.ko.mapping.list(in.type = in.type, out.type = out.type, 
                                       species = species, in.ids = in.ids)
    message("\nGenerated mapping list using KEGGREST")
    message("Saving mapping list to current working directory")
    
    file.name <- paste(species, toupper(in.type), toupper(out.type), sep = "_")
    file.save.path <- paste(SBGNview.data.folder, paste(file.name, ".RData", sep = ""), sep = "/")
    assign(file.name, id.map)
    save(list = file.name, file = file.save.path)
    message("\nSaved generated mapping table at: ", file.save.path)
    
  } else { # input/output not KO
    if (species == "mmu") {
      species <- "mouse"
    }
    if (is.null(in.ids)) {
      stop("\nMust provide input IDs when using pathview mapping!")
    }
    id.map <- pathview::geneannot.map(in.ids = in.ids, in.type = in.type, out.type = out.type,
                                      org = species, unique.map = unique.map)
  }
  return(id.map)
}

#########################################################################################################
### generate mapping list on the fly for mapping between KO and other IDs in gene.idtype.list
# Pseudocode for this function:
# input = in.type, out.type ("ko"), species, in.ids
#   if species not in bods 
#       if in.type = entrez id
#           if can't be map to entrez -> stop
#           if "kegg.geneid" == "ncbi.geneid" for species 
#               no need of 2nd mapping table. Map from KEGG ID to KO using keggLink
#           else use two mapping tables 
#                entrez id to kegg id and kegg id to ko and merge the two lists for (entrez, ko) table
#       else 
#           if between KO and ncbiprotinid / uniprot, do mapping
#           else stop(can't map)
#   else (species in bods) 
#       if in.type == entrez
#           map to entrez similar to above
#       else if in.type = KO. map ko to type in gene.id.type.list for bods species
#       else use pathview::id2eg to map to entrezid
#           if "kegg.geneid" == "ncbi.geneid" for species 
#               map directly from KEGG ID to KO
#               merge with id2eg map to get (in.type, ko) table
#           else use 2nd mapping table
#               map from entrez id to keggid merge with id2eg map list
#               map from kegg id to ko and merge with first mereged list
generate.ko.mapping.list <- function(in.type, out.type, species, in.ids = NULL) {
  
  in.type <- tolower(in.type)
  out.type <- tolower(out.type)
  species <- tolower(species)
  
  # korg cols:: 3 = kegg.code; 6 = entrez.gnodes; 7 = kegg.geneid; 8 = ncbi.geneid
  if(!exists("korg")) data(korg, package="pathview")
  if(!exists("bods")) data(bods, package="pathview")
  
  if(!species %in% korg[, 3] & !species %in% bods[, 3]){ # stop if species is NOT in korg or bods
    stop("Incorrect KEGG species code!")
  } 
  
  if(!species %in% bods[, 3]) { # species not in bods
    # get species row index from korg
    ridx <- match(species, korg[, 3]) %% nrow(korg)
    message("'", species, "'", " species in korg dataset from Pathview")
    
    if(any(in.type %in% c("entrez", "eg", "entrezid"))){ # in.type == entrez
      
      if(is.na(korg[ridx, 8])){ # if ncbi.geneid == NA, can't map to entrezid
        stop("Can't map between ENTREZID and KO for ", species)
      }
      if(korg[ridx, 6] == "1"){   # "kegg.geneid" == "ncbi.geneid" (entrez.gnodes == 1)
        message("Mapping directly from ENTREZID to KO")
        mapping.list <- get.keggrest.data(tar = out.type,  src = species, 
                                          link.or.conv = "link", src.is.species = T)
        colnames(mapping.list) <- c("ENTREZID", "KO")
        
      } else { # kegg ID is NOT entrez id. entrez.gnodes == 0. need 2 mapping tables
        
        message("Mapping from ENTREZID to KEGG ID") # ncbi-geneid to kegg id
        conv.list <- get.keggrest.data(tar = species, src = "ncbi-geneid",
                                       link.or.conv = "conv", tar.is.species = T)
        
        message("Mapping from KEGG ID to KO") # kegg id to ko
        link.list <- get.keggrest.data(tar = "ko", src = species, 
                                       link.or.conv = "link", src.is.species = T)
        
        message("Merging mapping lists")
        merge.list <- merge(conv.list, link.list) # merge list. cols = keggid, ncbi-geneid, ko
        mapping.list <- merge.list[,2:3] # take only ncbi-geneid, ko
      }
    } else { # in.type is NOT entrez and species NOT in bods
      # map from in.type to kegg gene id =>  ko
      ## if input/output is ncbi-proteinid or uniprot for species in korg
      if(any(c(in.type, out.type) %in% c("ko")) & 
         any(c(in.type, out.type) %in% c("ncbi-proteinid", "uniprot"))) {
        in.type <- setdiff(c(in.type, out.type), "ko")
        out.type <- "ko"
        
        message("Mapping from ", toupper(in.type), " to KEGG ID") # in.type to kegg id
        conv.list <- get.keggrest.data(tar = species, src = in.type, 
                                       link.or.conv = "conv", tar.is.species = T)
        
        message("Mapping from KEGG ID to KO") # kegg id to ko
        link.list <- get.keggrest.data(tar = "ko", src = species, 
                                       link.or.conv = "link", src.is.species = T)
        
        message("Merging mapping lists")
        merge.list <- merge(conv.list, link.list) # merge mapping lists
        mapping.list <- merge.list[,2:3]
        
      } else {
        stop("This mapping cannot be done")
      }
    }
    
  } else { #### species in bods
    
    message("'", species, "'", " species in bods dataset from Pathview")
    # get species info from korg
    ridx <- match(species, korg[, 3]) %% nrow(korg)
    
    if(any(in.type %in% c("entrez", "eg", "entrezid"))) { # input type is entrez for bods species
      
      message("ID mapping from ENTREZID to KO for bods species")
      
      if(korg[ridx, 6] == "1"){   # "kegg.geneid" == "ncbi.geneid". entrez.gnodes == 1
        # no need of 2nd mapping table
        message("Mapping directly from ENTREZID to KO for bods species")
        mapping.list <- get.keggrest.data(tar = out.type, src = species,
                                          link.or.conv = "link", src.is.species = T)
        colnames(mapping.list) <- c("ENTREZID", "KO")
        
      } else { # kegg ID is NOT entrez id. entrez.gnodes == 0. need 2 mapping tables
        
        message("Mapping from ENTREZID to KEGG ID for bods species") # ncbi-geneid to kegg id
        conv.list <- get.keggrest.data(tar = species, src = "ncbi-geneid",
                                       link.or.conv = "conv", tar.is.species = T)
        
        message("Mapping from KEGG ID to KO for bods species") # kegg id to ko
        link.list <- get.keggrest.data(tar = "ko", src = species,
                                       link.or.conv = "link", src.is.species = T)
        
        message("Merging mapping lists")
        merge.list <- merge(conv.list, link.list) # merge list. cols = keggid, ncbi-geneid, ko
        mapping.list <- merge.list[,2:3] # take only ncbi-geneid, ko
        
      } # end if condition for input type is entrez for species in bods
      
    } else if (in.type == "ko") { ### if input is ko, and species in bods
      # map from KO to output type in gene.idtype.list for bods species
      # input.ko.ids required
      message("KO to ", toupper(out.type))
      mapping.list <- mapping.ko.to.arbitrary.id.type(input.ko.ids = in.ids, 
                                                      output.type = toupper(out.type),
                                                      species = species, 
                                                      ridx = ridx)
      
    } else { # input type not entrez. other input type
      # map from in.type to ENTREZID using pathview
      message("Mapping from ", toupper(in.type), " to ENTREZID for '", 
              species, "' using pathview::id2eg")
      
      if(is.null(in.ids)){
        stop("Need vector of input ids to use pathview::id2eg")
      }
      
      in.to.eg <- pathview::id2eg(ids = in.ids, category = in.type, org = species) # in.type, entrez
      
      if(korg[ridx, 6] == "1") { # check if KEGG id == ENTREZID in korg for bods species
        
        message("Mapping directly from ENTREZID to KO for bods species") 
        map.list <- get.keggrest.data(tar = out.type, src = species,  
                                      link.or.conv = "link", src.is.species = T) # kegg, ko
        colnames(map.list) <- c("ENTREZID", "KO")
        
        message("Merging mapping lists for bods species")
        merge.list <- merge(in.to.eg, map.list) # merge in.type, entrezid with keggid, ko
        mapping.list <- merge.list[,2:3] 
        
      } else { # KEGG id != ENTREZID. need 2nd mapping table
        message("Mapping from ENTREZID to KEGG ID for bods species") # ncbi-geneid to kegg id
        conv.list <- get.keggrest.data(tar = species, src = "ncbi-geneid",
                                       link.or.conv = "conv", tar.is.species = T)
        
        # in.type, entrezid merge with entrezid, keggid
        merge.list.1 <- merge(in.to.eg, conv.list)
        map.list.1 <- merge.list.1[,2:3]  # in.type, keggid
        
        message("Mapping from KEGG ID to KO for bods species") # kegg id to ko
        link.list <- get.keggrest.data(tar = "ko", src = species,
                                       link.or.conv = "link", src.is.species = T)
        
        # merge in.type, keggid with  keggid, ko
        message("Merging mapping lists for bods species")
        merge.list.2 <- merge(map.list.1, link.list) 
        mapping.list <- merge.list.2[,2:3] # in.type, ko
      }
      
    } # end else in.type != entrez
    
  } # end else species in bods
  
  return(mapping.list)
}

#########################################################################################################
# mapping from input type KO to output tpye (in gene.idtype.list) for species in bods dataset
# using keggrest and pathview::eg2id. Need vector of input of KO ids
mapping.ko.to.arbitrary.id.type <- function(input.ko.ids = NULL, output.type, species, ridx) {
  
  if(is.null(input.ko.ids)){
    stop("Need vector of input KO IDs")
  }
  
  if(!toupper(output.type) %in% gene.idtype.list) {
    stop("Output type not in 'gene.id.type.list'")
  }
  
  if(korg[ridx, 6] == "1") {   # default "kegg.geneid" is "ncbi.geneid". entrez.gnodes == 1
    
    message("Mapping KO to ENTREZID for bods species") 
    ko.to.entrez <- get.keggrest.data(tar = species, src = "ko",
                                      link.or.conv = "link", tar.is.species = T) # ko, keggid/entrez 
    colnames(ko.to.entrez) <- c("KO", "ENTREZID")
    
  } else { # kegg Id is not deault. so map KO to KEGGID, KEGGID to ncbiID, merge
    
    message("Mapping KO to KEGGID")
    ko.to.kegg <- get.keggrest.data(tar = species, src = "ko",
                                    link.or.conv = "link", tar.is.species = T)
    message("Mapping KEGGID to NCBI-GENEID")
    kegg.to.ncbi <- get.keggrest.data(tar = "ncbi-geneid", src = species,
                                      link.or.conv = "conv", src.is.species = T)
    message("Merging mapping lists")
    ko.to.ncbi <- merge(ko.to.kegg, kegg.to.ncbi)
    ko.to.entrez <- ko.to.ncbi[, 2:3]
    colnames(ko.to.entrez) <- c("KO", "ENTREZID")
  }
  
  # match input KO ids with entrez IDS in generated mapping table
  input.ko.ids <- unique(input.ko.ids)
  ko.to.entrez <- as.matrix(ko.to.entrez)
  filterd.ko.to.entrez <- subset(ko.to.entrez, ko.to.entrez[,1] %in% input.ko.ids)
  
  # use pathview::eg2id to map entrez to output type and merge to get KO, SYMBOL mapping list
  filterd.ko.to.entrez <- as.matrix(filterd.ko.to.entrez)
  entrez.to.output <- pathview::eg2id(eg = filterd.ko.to.entrez[,2], 
                                      category = toupper(output.type),
                                      org = species)
  mapping.list <- merge(filterd.ko.to.entrez, entrez.to.output)
  mapping.list <- mapping.list[, 2:3]
  
  return(mapping.list)
}

#########################################################################################################
# get data using keggLink and keggConv to get mapping data and format to a list for generate.ko.mapping.list
get.keggrest.data <- function(tar, src, link.or.conv, tar.is.species = F, src.is.species = F) {
  # tar and src should be lower case
  if(link.or.conv == "link"){ # keggLink
    kegg.data <- KEGGREST::keggLink(target = tar, source = src)
  } else { # keggConv
    kegg.data <- KEGGREST::keggConv(target = tar, source = src)
  }
  
  ids1 <- gsub("^.+:", "", names(kegg.data))
  ids2 <- gsub("^.+:", "", kegg.data)
  map.list <- cbind(ids1, ids2)
  rownames(map.list) <- NULL
  
  if(src == "ncbi-geneid") { src <- "ENTREZID" }
  if(src.is.species == TRUE) {src <- paste("KEGGID", src, sep = "-")}
  if(tar.is.species == TRUE) {tar <- paste("KEGGID", tar, sep = "-")}
  
  colnames(map.list) <- c(toupper(src), toupper(tar))
  
  return(map.list)
}

#########################################################################################################
# used in loadMappingTable function
# using chebi as bridge ID to map from IDs in data(cpd.simtypes) from pathview to 
# glyph IDs of pathwayCommons, pathway.id, or metacyc.SBGN
# if input type and output type in data(cpd.simtypes), use pathview::cpdidmap
# else map to kegg using pathview::cpd2kegg, kegg to chebi, chebi to glyph ID
generate.cpd.mapping.table <- function(input.type, output.type, limit.to.ids) {

  if(!exists("cpd.simtypes")) data("cpd.simtypes", package = "pathview")
  # input/output validation
  valid.in.out.types <- c(cpd.simtypes, "KEGG", "kegg", "pathwayCommons", "pathway.id", "metacyc.SBGN", "compound.name")
  if(!input.type %in% valid.in.out.types | !output.type %in% valid.in.out.types) {
    stop("'input.type' needs to be one of values in data(cpd.simtypes) or 'KEGG' if input.type is KEGG COMPOUND/DRUG accession.
    'output.type' can be one of values in data(cpd.simtypes) or one of 'pathwayCommons','pathway.id', 'metacyc.SBGN', 'compound.name'.")
  }
  if(is.null(limit.to.ids)) {
    stop("\nMust provide input IDs to 'limit.to.ids' argument when using pathview mapping!")
  }
  if(input.type %in% c("pathwayCommons", "pathway.id", "metacyc.SBGN", "KEGG", "kegg") & output.type %in% cpd.simtypes) {
    temp <- input.type
    input.type <- output.type
    output.type <- temp
  }
  if(input.type %in% c("pathwayCommons", "pathway.id", "metacyc.SBGN") & output.type %in% c("KEGG", "kegg")) {
    temp <- input.type
    input.type <- output.type
    output.type <- temp
  }
  # map between IDs in cpd.simtypes
  if(input.type %in% cpd.simtypes & output.type %in% c(cpd.simtypes, "KEGG", "kegg")) {
    message("Using pathview::cpdidmap to map ", input.type, " and ", output.type)
    mapping.list <- pathview::cpdidmap(in.ids = limit.to.ids, in.type = input.type,
                                       out.type = output.type)
  } else { # map to glyph ID type
    input.cpd.to.kegg <- limit.to.ids
    # if input.type not native KEGG, map to kegg ID first
    if(!input.type %in% c("KEGG COMPOUND accession", "KEGG DRUG accession", "KEGG", "kegg")) {
      message("Using pathview::cpd2kegg to map ", input.type, " to KEGG")
      input.cpd.to.kegg <- pathview::cpd2kegg(in.ids = limit.to.ids, in.type = input.type)
      colnames(input.cpd.to.kegg) <- c(input.type, "kegg")
    } else {
      # native kegg
      message("Input is a native KEGG compound ID type")
      input.cpd.to.kegg <- matrix(input.cpd.to.kegg)
      colnames(input.cpd.to.kegg) <- c("kegg")
      input.type <- "kegg"
    }
    # chebi_kegg.ligand doesn't contain DRUG and COMPOUND accession data
    # using chebi_kegg that contains both accessions data
    message("Mapping ", input.type, " to ", output.type)
    if(!exists("chebi_kegg")) data("chebi_kegg", package = "SBGNview.data")
    chebi.to.kegg <- as.matrix(get("chebi_kegg"))
    
    input.type.to.chebi <- merge(input.cpd.to.kegg, chebi.to.kegg)
    input.type.to.chebi <- input.type.to.chebi[ , c(input.type, "chebi")]
    
    chebi.to.ouput.type.map.table <- paste("chebi_", output.type, sep = "")
    if(!exists(chebi.to.ouput.type.map.table)) { 
      data(list = chebi.to.ouput.type.map.table, package = "SBGNview.data") 
    }
    chebi.to.output.type <- as.matrix(get(chebi.to.ouput.type.map.table))
    
    mapping.list <- merge(input.type.to.chebi, chebi.to.output.type)
    mapping.list <- mapping.list[, c(input.type, output.type)]
    
    if(length(mapping.list) == 0 | dim(mapping.list)[1] == 0) {
      stop("Mapping from ", input.type, " to ", output.type, " did not work. Please check input IDs type and try again.")
    }
  } # end else map to glyph ID
  message("Finished mapping from ", input.type, " to ", output.type)
  
  return(mapping.list)
}

#########################################################################################################
# find all mapping tables for gene and compound id types and combine into a list
load.id.mapping.list.all <- function(SBGN.file.cpd.id.type = NULL, SBGN.file.gene.id.type = NULL, 
                                     output.gene.id.type = NULL, output.cpd.id.type = NULL, species = NULL, SBGNview.data.folder = "./SBGNview.tmp.data") {
  idmapping.all.list <- list()
  cpd.id.mapping.list <- list()
  gene.id.mapping.list <- list()
  if (!is.null(output.cpd.id.type)) {
    if (SBGN.file.cpd.id.type != output.cpd.id.type) {
      input.cpd.type <- SBGN.file.cpd.id.type
      cpd.mapping.table <-  loadMappingTable(output.type = output.cpd.id.type, 
                                             input.type = input.cpd.type, mol.type = "cpd", species = species, 
                                             SBGNview.data.folder = SBGNview.data.folder)
      type.pair.name <- paste(sort(c(input.cpd.type, output.cpd.id.type), method = "radix",
                                   decreasing = TRUE), collapse="_")
      cpd.id.mapping.list[["cpd"]] <- list()
      cpd.id.mapping.list[["cpd"]][[type.pair.name]] <- cpd.mapping.table
    }
  }
  
  if (!is.null(output.gene.id.type)) {
    if (SBGN.file.gene.id.type != output.gene.id.type) {
      input.gene.type <- SBGN.file.gene.id.type
      gene.mapping.table <- loadMappingTable(output.type = output.gene.id.type, 
                                             input.type = input.gene.type, mol.type = "gene", species = species, 
                                             SBGNview.data.folder = SBGNview.data.folder)
      type.pair.name <- paste(sort(c(input.gene.type, output.gene.id.type), method = "radix",
                                   decreasing = TRUE), collapse="_")
      gene.id.mapping.list[["gene"]] <- list()
      gene.id.mapping.list[["gene"]][[type.pair.name]] <- gene.mapping.table
    }
  }
  id.mapping.all.list <- c(cpd.id.mapping.list, gene.id.mapping.list)
  return(id.mapping.all.list)
}

#########################################################################################################
# get all ID mapping tables for gene and compound ID types
load.all.ids.mapping <- function(database, all.pairs.id.mapping.list, species, output.gene.id.type, 
                                 output.cpd.id.type, SBGNview.data.folder = "./SBGNview.tmp.data") {
  
  if (database == "MetaCrop") {
    sbgn.id.attr <- "label"
    SBGN.file.gene.id.type <- "ENZYME"
    SBGN.file.cpd.id.type <- "compound.name"  # "CompoundName"
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
#' Retrieve gene list or compound list from collected databases
#' 
#' @param id.type A character string. Default: "ENTREZID". The ID type of output gene list. One of the supported types in \code{data('mapped.ids')}
#' @param mol.type A character string. Default: "gene". One of 'gene' or 'cpd'
#' @param species A character string. Default: "hsa". The three letter species code used by KEGG. E.g. 'hsa','mmu'
#' @param database A character string. Default: "pathwayCommons". The database where gene list will be extracted. Acceptable values: 'MetaCyc', 'pathwayCommons', 'MetaCrop'. The value is case in-sensitive.
#' @param output.pathway.name Logical. Default: T. If set to 'TRUE', the names of returned list are in the format: 'pathway.id::pathway.name'. If set to 'FALSE', the format is 'pahtway.id'
#' @param combine.duplicated.set Logical. Default: T. Some pathways have the same geneset. If this parameter is set to 'TRUE', the output list will combine pathways that have the same gene set. The name in the list will be pathway names concatinated with '||'
#' @param truncate.name.length Integer. Default: 50. The pathway names will be truncated to at most that length. 
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files.
#' @return A list. Each element is a genelist of a pathway.
#' @examples 
#' data(pathways.info)
#' mol.list <- sbgn.gsets(id.type = 'ENTREZID',
#'                        species = 'hsa',
#'                        database = 'pathwayCommons')
#'   
#' @export

sbgn.gsets <- function(id.type = "ENTREZID", mol.type = "gene", species = "hsa", database = "pathwayCommons", 
                       output.pathway.name = TRUE, combine.duplicated.set = TRUE, 
                       truncate.name.length = 50, SBGNview.data.folder = "./SBGNview.tmp.data") {
  
  if(id.type %in% c("ENTREZID", "eg", "entrez", "entrezid")) id.type <- "ENTREZID"
  if(id.type == "CompoundName") id.type <- "compound.name"
  if(id.type %in% c("kegg.ligand", "KEGG")) id.type <- "kegg"
  
  if (tolower(database) == "metacrop") {
    if (mol.type == "gene") {
      id.in.pathway <- "ENZYME"
    } else {
      id.in.pathway <- "compound.name" # "CompoundName"
    }
    output.pathways <- subset(pathways.info, database == "MetaCrop", select = "pathway.id")
  } else if (tolower(database) %in% c("pathwaycommons", "metacyc")) {
    if (mol.type == "gene") {
      id.in.pathway <- "KO"
    } else {
      id.in.pathway <- "chebi"
    }
    output.pathways <- subset(pathways.info, database != "MetaCrop", select = "pathway.id")
  } else { # if wrong value for database
    stop("'database' argument value is incorrect. Acceptable values are: 'MetaCyc', 'pathwayCommons', 'MetaCrop'")
  }
  # metacrop initial list is using enzyme
  ref.to.pathway <- loadMappingTable(input.type = id.in.pathway, output.type = "pathway.id", 
                                   mol.type = mol.type, species = species, SBGNview.data.folder = SBGNview.data.folder)

  if (id.type == id.in.pathway) {
    out.id.to.pathway <- ref.to.pathway
    # change KO to output id
  } else {
    message("Mapping ", id.in.pathway, " to ", id.type)
    
    out.id.type.to.ref <- loadMappingTable(input.type = id.in.pathway, output.type = id.type, 
                                           mol.type = mol.type, limit.to.ids = ref.to.pathway[, id.in.pathway], 
                                           species = species, SBGNview.data.folder = SBGNview.data.folder)
    
    out.id.type.to.ref <- out.id.type.to.ref 
    # merge KO to pathway and KO to output id
    out.id.to.pathway <- merge(out.id.type.to.ref, ref.to.pathway, all.x = TRUE)
    out.id.to.pathway <- out.id.to.pathway[, c(id.type, "pathway.id")]
    out.id.to.pathway <- unique(out.id.to.pathway)
    out.id.to.pathway <- out.id.to.pathway[!is.na(out.id.to.pathway[, 2]), ]
  }
  out.id.to.pathway <- unique(out.id.to.pathway)
  out.id.to.pathway <- out.id.to.pathway[out.id.to.pathway[, "pathway.id"] %in% 
                                           output.pathways$pathway.id, ]
  out.id.to.pathway <- out.id.to.pathway[out.id.to.pathway[, id.type] != 
                                           "", ]
  out.id.to.pathway <- split(as.character(out.id.to.pathway[, id.type]), 
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
#' Download pre-generated SBGN-ML file from GitHub
#' 
#' This function will download a SBGN-ML file from our pre-collected SBGN-ML files given the 'pathway.id' argument.
#' 
#' @param pathway.id A character string. The ID of pathway. For accepted pathway IDs, please check \code{data('pathways.info')}. IDs are in column 'pathway.id' (pathways.info[,'pathway.id'])
#' @param download.folder A character string. Default: "./". The output folder to store created SBGN-ML files.
#' @return A vector of character strings. The path to the created SBGN-ML files.
#' @examples 
#' data("pathways.info")
#' data("sbgn.xmls")
#' input.sbgn <- downloadSbgnFile(pathway.id = pathways.info[1,'pathway.id'],
#'                               download.folder = './')
#' @export

downloadSbgnFile <- function(pathway.id, download.folder = "./") {
  
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
# find pathways in pathways.info based on input keywords and keyword matching logic
find.pathways.by.keywords <- function(keywords, keyword.type, keywords.logic, mol.name.match, 
                                      SBGNview.data.folder = "./SBGNview.tmp.data/") {
  if (is.null(keywords)) {
    pathways <- pathways.info[, c("pathway.id", "pathway.name", "sub.database", "database")]
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
    pathways <- pathways.info[if.selected, c("pathway.id", "pathway.name", "sub.database", "database")]
    
  } else {
    # using IDs to search for pathway type.pair.name =
    # paste(keyword.type,'_pathway.id',sep='')
    if(!exists("mapped.ids")) data(mapped.ids)
    # if (keyword.type %in% mapped.ids$gene) {
    if (keyword.type %in% unique(mapped.ids$gene[,2])) {
      mol.type <- "gene"
    } else if (keyword.type %in% mapped.ids$cpd) {
      mol.type <- "cpd"
    }
    mapping.table <- loadMappingTable(input.type = keyword.type, output.type = "pathway.id", 
                                     species = NULL, SBGNview.data.folder = SBGNview.data.folder, mol.type = mol.type)
    
    if (keyword.type == "compound.name") { # keyword.type == "CompoundName"
      
      if (keywords.logic == "and") {
        keywords.logic <- "or"
        warning("Searching pathways by Compound Names. Using keywords.logic='or'. Please don't set keywords.logic to 'and'!!!\n")
      }
      keywords <- tolower(keywords)
      mapping.table[, keyword.type] <- tolower(as.character(mapping.table[, keyword.type]))
      if (mol.name.match == "exact.match") {
        if.selected <- tolower(mapping.table[, keyword.type]) %in% keywords
        pathways <- mapping.table[if.selected, ]
      } else if (mol.name.match == "presence.of.input-string.in.target-name") {
        if.selected <- grepl(pattern = paste(keywords, collapse = "|"), ignore.case = TRUE, 
                             mapping.table[, keyword.type])
        sum(if.selected)
        pathways <- mapping.table[if.selected, ]
      } else if (mol.name.match == "jaccard") {
        best.match <- match.names(input.names = tolower(keywords), 
                                  output.names = tolower(mapping.table[, keyword.type]))
        pathways <- merge(best.match, mapping.table, by.x = "output.name", 
                          by.y = keyword.type)
        names(pathways)[1] <- c("compound.name")
      }
    } else {
      pathways <- mapping.table[tolower(mapping.table[, keyword.type]) %in% tolower(keywords), ]
    }
  }
  return(pathways = pathways)
}

#########################################################################################################
# given two vectors of characters/words, this function find the best match
# between words in the two vectors, by splitting a word into sub-strings and use
# 'jaccard' to calculate two words' similarity (similarity between their
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
# find pathways based on organism and pathway completeness in a species
filter.pathways.by.org <- function(pathways, org, pathway.completeness, 
                                   pathways.info.file.folder = "./SBGNview.tmp.data") {
  
  org <- tolower(org)
  if(!exists("pathway.species.pct_Mapped")) data("pathway.species.pct_Mapped")
  if (org == "all") {
    org <- unique(pathway.species.pct_Mapped$species)
  }
  pathway.species.pct_Mapped[, "species"] <- as.character(pathway.species.pct_Mapped[, "species"])
  pathway.species.pct_Mapped[, "pathway"] <- as.character(pathway.species.pct_Mapped[,  "pathway"])
  
  if (!is.null(pathway.completeness)) {
    message("Using user provided pathway completeness cutoff: the same cutoff for different pathways!\n")
    
    pathway.ids <- pathway.species.pct_Mapped[pathway.species.pct_Mapped$species %in% 
                                                org & pathway.species.pct_Mapped$pct.mapped.species.pathway >= pathway.completeness, ]
  } else {
    message("Using pre-generated pathway-specific completeness cutoff: different cutoff for different pathways! Pathway 'exist' and 'not exist' accross all species defined by this cutoff have the largest ANOVA F statistic when comparing completeness between 'exist' and 'not exist' groups. \n")
    pathway.ids <- pathway.species.pct_Mapped[pathway.species.pct_Mapped$species %in% org, ]
    if(!exists("pathway.completeness.cutoff.info")) data("pathway.completeness.cutoff.info")
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
#' @param keywords.logic A character string. Options are 'and' or 'or' (default). This will tell the function if the search require 'all' or 'any' of the keywords to be present. It only makes difference when keyword.type is 'pathway.name'.
#' @param keyword.type A character string. Either 'pathway.name' (default) or one of the ID types in \code{data('mapped.ids')}
#' @param org  A character string. The KEGG species code.
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold downloaded ID mapping files and pathway information data files.
#' @details If 'keyword.type' is 'pathway.name' (default), this function will search for the presence of any keyword in the pathway.name column of data(pathways.info). The search is case in-sensitive. If 'keyword.type' is one of the identifier types and 'keywords' are corresponding identifiers, this function will return pathways that include nodes mapped to input identifiers. 
#' @return A dataframe. Contains information of pathways found.
#' @examples 
#' data(pathways.info)
#' input.pathways <- findPathways('Adrenaline and noradrenaline biosynthesis')
#' @export

findPathways <- function(keywords = NULL, keywords.logic = "or", keyword.type = "pathway.name", 
                         org = NULL, SBGNview.data.folder = "./SBGNview.tmp.data") {
  
  pathway.completeness = NULL
  mol.name.match = c("exact.match", "jaccard", 
                     "presence.of.input-string.in.target-name")[1]
  if (!file.exists(SBGNview.data.folder)) {
    dir.create(SBGNview.data.folder)
  }
  
  keywords <- as.vector(keywords)
  # changed CompoundName to compound.name and kegg.ligand to kegg in mapped.ids and mapping file names
  if(keyword.type == "CompoundName") keyword.type <- "compound.name"
  if(keyword.type %in% c("kegg.ligand", "KEGG")) keyword.type <- "kegg"
  
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
############# Below are previous versions of some of the above functions. 
############# These are kept for reference purposes.
#########################################################################################################
### original (old version) function replaced with new version. 
# download.mapping.file <- function(input.type, output.type, 
#                                   species = NULL, SBGNview.data.folder = "./SBGNview.tmp.data/") {
#   options(warn = -1)
#   # R CMD check will give different sort result if you didn't specify 'method': the
#   # default sort method depends on the locale of your system. And R CMD check uses
#   # a different locale than my interactive session.  The issue resided in this: R
#   # interactively used:LC_COLLATE=en_US.UTF-8; R CMD check used: LC_COLLATE=C;
#   # https://stackoverflow.com/questions/42272119/r-cmd-check-fails-devtoolstest-works-fine
#   type.pair.name.1 <- paste(sort(c(input.type, output.type), method = "radix", 
#                                  decreasing = FALSE), collapse = "_")
#   type.pair.name.2 <- paste(sort(c(input.type, output.type), method = "radix", 
#                                  decreasing = TRUE), collapse = "_")
#   type.pair.name.1.org <- paste(species, type.pair.name.1, sep = "_")
#   type.pair.name.2.org <- paste(species, type.pair.name.2, sep = "_")
#   try.file.names <- c(type.pair.name.1.org, type.pair.name.2.org, type.pair.name.1, 
#                       type.pair.name.2)
#   if.online.mapping.avai <- FALSE
#   if.mapping.in.data.package <- FALSE
#   for (i in seq_len(length(try.file.names))) {
#     type.pair.name.try <- try.file.names[i]
#     if.data.exists <- tryCatch(data(list = type.pair.name.try), warning = function(w) "no such data")
#     if (if.data.exists %in% c(type.pair.name.try, "mapping.list", "mapping.table")) {
#       mapping.file.name <- type.pair.name.try
#       message("ID mapping ", type.pair.name.try, " is included in SBGNview.data package. Loading it.")
#       if.mapping.in.data.package <- TRUE
#       (break)()
#     } else {
#       print("ID mapping not in data package for this pair")
#       print(type.pair.name.try)
#       print(if.data.exists)
#     }
#   }
#   if (!if.mapping.in.data.package) {
#     # stop("Can't map this pair of IDs with SBGNview.data!", type.pair.name.try,  "\n")
#     warning("Can't map this pair of IDs with SBGNview.data! Generating mapping with Pathview!", type.pair.name.try,  "\n")
#     mapping.file.name = "online mapping not available"
#   }
#   options(warn = 1)
#   return(mapping.file.name)
# }

#########################################################################################################
### original (old version) of the function
# loadMappingTable <- function(output.type, input.type, species = NULL, cpd.or.gene, 
#                              limit.to.ids = NULL, SBGNview.data.folder = "./SBGNview.tmp.data"){ 
#   
#   input.type = gsub("entrez","ENTREZID",input.type)
#   output.type = gsub("entrez","ENTREZID",output.type)
#   if(input.type == output.type){
#     stop("Input type and output types are the same!")
#   }
#   species = gsub("Hs","hsa",species)
#   type.pair.name = paste(sort(c(input.type,output.type),method = "radix",decreasing = TRUE),collapse="_") # R CMD check will give different sort result if we didn't specify "method":  the default sort method depends on the locale of your system. And R CMD check uses a different locale to R interactive session. The issue resided in this: R interactively used:LC_COLLATE=en_US.UTF-8; R CMD check used: LC_COLLATE=C; https://stackoverflow.com/questions/42272119/r-cmd-check-fails-devtoolstest-works-fine
#   if(!file.exists(SBGNview.data.folder)){
#     dir.create(SBGNview.data.folder)
#   }
#   
#   mapping.file.name <- download.mapping.file(input.type, output.type, species, SBGNview.data.folder = SBGNview.data.folder)
#   
#   if(file.exists(mapping.file.name) | tryCatch(data(list = mapping.file.name), 
#                                                warning = function(w) "no such data") %in% c(mapping.file.name, "mapping.list", "mapping.table" )
#   ){
#     message("Loading ID mapping file: ",mapping.file.name," \n")
#     .GlobalEnv$mapping.list = NULL; .GlobalEnv$mapping.table = NULL
#     if.from.file = FALSE
#     if(file.exists(mapping.file.name)){
#       if.from.file = TRUE
#       var.name = load(mapping.file.name)
#     }else{
#       if(!exists(mapping.file.name)){
#         data(list = mapping.file.name)
#       }
#     }
#     if(!is.null(mapping.list)){ # if one of "output.type" or "input.type" is NOT "pathway.id", then the data's names is "mapping.list"
#       id.map = mapping.list[[1]][[1]]
#     }else if(!is.null(mapping.table)){ # if one of "output.type" or "input.type" is "pathway.id", then the data's names is "mapping.table"
#       id.map = mapping.table
#     }else if(if.from.file){ # if one of "output.type" or "input.type" is "pathway.id", then the data's names is "mapping.table"
#       id.map = get(var.name)
#     }else{
#       # print("loading mapping file from data package")
#       # print(mapping.file.name)
#       # data(list = mapping.file.name)
#       id.map = get(mapping.file.name)
#     }
#     message("Finished loading")
#     
#     if("species" %in% colnames(id.map) &  !is.null(id.map) ){
#       if(any(id.map[,"species"] %in% species) ){
#         id.map = id.map[id.map[,"species"] %in% species,c(input.type,output.type)]
#       }
#     }
#   } else {
#     if(cpd.or.gene == "gene"){
#       # pathway.id used for heterogeneous purposes and mapping gene to pathway
#       # pathwayCommons and metacyc.SBGN are databases mentioned in pathways.stats
#       # unique(pathways.info$macromolecule.ID.type) returns "pathwayCommons" "metacyc.SBGN"  "ENZYME"  
#       # ENZYME not in if condition below b/c it can be mapped directly 
#       if(any(c(input.type,output.type) %in% c("pathwayCommons","metacyc.SBGN","pathway.id")) & 
#          !any(c(input.type,output.type) %in% c("KO")) 
#       ){ # if input and output don't contain KO
#         
#         # load mapping table for mapping KO to glyph id
#         ko.to.glyph.id = loadMappingTable(input.type = output.type,
#                                           output.type = "KO",
#                                           cpd.or.gene = "gene",
#                                           species = species,
#                                           SBGNview.data.folder = SBGNview.data.folder)
#         ko.to.glyph.id = ko.to.glyph.id[[1]][[1]]
#         
#         # load mapping table for mapping user input to KO id
#         input.to.ko = loadMappingTable(input.type = input.type,
#                                        output.type = "KO",
#                                        cpd.or.gene = "gene",
#                                        limit.to.ids = limit.to.ids,
#                                        species = species,
#                                        SBGNview.data.folder = SBGNview.data.folder)
#         input.to.ko = input.to.ko$gene[[1]]
#         input.to.glyph.id = merge(input.to.ko,ko.to.glyph.id,all= FALSE)
#         id.map = input.to.glyph.id[,c(input.type,output.type)]
#       } else {
#         message("\nID mapping not pre-generated. Using Pathview!!!\n")
#         id.map = geneannot.map.ko(in.ids = limit.to.ids,
#                                   in.type = input.type,
#                                   out.type = output.type,
#                                   species = species,
#                                   unique.map= FALSE,
#                                   SBGNview.data.folder = SBGNview.data.folder)
#         if(is.vector(id.map)){
#           id.map = as.matrix(t(id.map))
#         }
#       }
#     } else if(cpd.or.gene == "compound"){
#       message("\nID mapping not pre-generated. Using Pathview cpd!!!\n")
#       if(is.null(limit.to.ids)){
#         stop("Must provide input IDs when using pathview mapping!")
#       }
#       print(input.type)
#       print(output.type)
#       id.map =  pathview::cpdidmap(in.ids = limit.to.ids, in.type = input.type, 
#                                    out.type = output.type)
#     } else {
#       message("\nCouldn't fine ID mapping table between ", input.type, " and ", output.type, "!!!\n")
#       #message("Tried online resource ", online.mapping.file, "\n")
#       message("Please provide ID mapping table using \"id.mapping.table\"!!\n")
#       stop("ID mapping table not provided")
#     }
#     
#     id.map = id.map[!is.na(id.map[,2]),]
#     if(is.vector(id.map)){
#       id.map = as.matrix(t(id.map))
#     }
#     id.map = id.map[!is.na(id.map[,1]),]
#     if(is.vector(id.map)){
#       id.map = as.matrix(t(id.map))
#     }
#     id.map = unique(id.map)
#   }
#   
#   # if(!is.null(limit.to.ids)){
#   #     # some IDs are quoted, need to remove the quote characters
#   #     id.map[,input.type] = gsub("^\"","",id.map[,input.type])
#   #     id.map[,input.type] = gsub("\"$","",id.map[,input.type])
#   #     id.map = id.map[id.map[,input.type] %in% limit.to.ids,]
#   # }
#   
#   # add additional mapping using KO to glyph.id
#   mapping.list = list()
#   if(is.vector(id.map)){
#     id.map = as.matrix(t(id.map))
#   }
#   id.map[,1] = as.character(id.map[,1])
#   id.map[,2] = as.character(id.map[,2])
#   id.map = as.matrix(id.map)
#   mapping.list[[cpd.or.gene]]= list()
#   mapping.list[[cpd.or.gene]][[type.pair.name]] = id.map
#   message("Generated ID mapping list")
#   
#   return(mapping.list)
# }

#########################################################################################################
### old version of function
# geneannot.map.ko <- function(in.ids = NULL, in.type, out.type, species = "all", unique.map, 
#     SBGNview.data.folder = "./SBGNview.tmp.data") {
#     # pathview's geneannot.map can't map KO, so here included KO mapping
#     if (!is.null(in.ids)) {
#         in.ids <- gsub("\"", "", in.ids)
#     }
#     message("\nusing pathview for id mapping: ", in.type, " to ", out.type, "\n\n")
#     
#     if (any(c(in.type, out.type) %in% "KO")) {
#         filter.type <- in.type
#         out.type <- setdiff(c(in.type, out.type), "KO")
#         in.type <- "KO"
#         
#         # break endless loop - this loop doesn't allow the rest of the code to run
#         # to break this, we need a mapping list. Since one doesn't exit, we generate mapping on the fly
#         # use KEGG API. KEGGREST Bioconductor Pkg
#         mapping.list <- loadMappingTable(input.type = "KO", output.type = "ENTREZID", 
#             cpd.or.gene = "gene", species = species, SBGNview.data.folder = SBGNview.data.folder)
#         #
#         message("loaded mapping list")
#         
#         mapping.table <- mapping.list[[1]][[1]]
#         if (out.type %in% c("ENTREZID", "ez", "entrezid")) {
#             print("filter species")
#             id.map <- mapping.table[mapping.table[, "species"] == species, c(in.type, 
#                 out.type)]
#             if (!is.null(in.ids)) {
#                 mapping.table <- mapping.table[mapping.table[, filter.type] %in% 
#                   in.ids, ]
#             }
#         } else {
#             if (species == "mmu") {
#                 species <- "mouse"
#             }
#             output.to.eg <- pathview::eg2id(eg = mapping.table[, "ENTREZID"], category = out.type, 
#                 org = species, unique.map = unique.map)
#             
#             output.to.ko <- merge(output.to.eg, mapping.table, all.x = TRUE)
#             id.map <- output.to.ko[, c("KO", out.type)]
#             id.map <- id.map[!is.na(id.map[, out.type]), ]
#             # filter to output only input IDs
#             if (!is.null(in.ids)) {
#                 id.map <- id.map[id.map[, filter.type] %in% in.ids, ]
#             }
#         }
#     } else {
#         if (species == "mmu") {
#             species <- "mouse"
#         }
#         if (is.null(in.ids)) {
#             stop("Must provide input IDs when using pathview mapping!")
#         }
#         
#         # id.map <- geneannot.map.all(in.ids = in.ids, in.type = in.type, out.type = out.type, 
#         #     org = species, unique.map = unique.map)
#         print(out.type)
#         id.map <- pathview::geneannot.map(in.ids = in.ids, 
#                                           in.type = in.type, 
#                                           out.type = out.type,
#                                           org = species, 
#                                           unique.map = unique.map)
#     }
#     return(id.map)
# }

#########################################################################################################
### this function was previously used in the older version geneannot.map.ko
### the functionality provided by this function was replaced with 
### pathview::geneannot.map in the updated geneannot.map.ko function 
# geneannot.map.all <- function(in.ids, in.type, out.type, org = "Hs", pkg.name = NULL, 
#                               unique.map = TRUE, na.rm = TRUE, keep.order = FALSE) {
#   
#   if (is.null(pkg.name)) {
#     # pkg.name=paste('org', org, 'eg.db', sep='.')
#     message("data bods")
#     data(bods)
#     ridx <- grep(tolower(paste0(org, "[.]")), tolower(bods[, 1]))
#     if (length(ridx) == 0) {
#       ridx <- grep(tolower(org), tolower(bods[, 2:3]))%%nrow(bods)
#       if (length(ridx) == 0) 
#         stop("Wrong org value!")
#       if (any(ridx == 0)) 
#         ridx[ridx == 0] <- nrow(bods)
#     }
#     pkg.name <- bods[ridx, 1]
#   }
#   all.mappings <- character(2) 
#   message("all mappings: ", all.mappings, " items")
#   for (i in seq_len(length.out = nrow(bods))) {
#     if (bods[i, 1] != pkg.name) {
#       (next)()
#     }
#     pkg.name <- bods[i, 1]
#     
#     if (!pkg.name %in% rownames(installed.packages())) { 
#       #(next)() # if package not installed, evals to True, but doesn't install pkg.
#       # install pkg.name
#       message(pkg.name, " is not installed. Installing ", pkg.name)
#       BiocManager::install(pkg.name, update=FALSE)
#       
#     }
#     message("Using package: ", pkg.name, "\n")
#     db.obj <- eval(parse(text = paste0(pkg.name, "::", pkg.name)))
#     id.types <- AnnotationDbi::columns(db.obj)  #columns(eval(as.name(pkg.name)))
#     
#     in.type <- toupper(in.type)
#     out.type <- toupper(out.type)
#     eii <- in.type == toupper("entrez") | in.type == toupper("eg")
#     if (any(eii)) 
#       in.type[eii] <- "entrez"
#     eio <- out.type == toupper("entrez") | out.type == toupper("eg")
#     if (any(eio)) 
#       out.type[eio] <- "entrez"
#     if (in.type == out.type) 
#       stop("in.type and out.type are the same, no need to map!")
#     
#     nin <- length(in.type)
#     if (nin != 1) 
#       stop("in.type must be of length 1!")
#     out.type <- out.type[!out.type %in% in.type]
#     nout <- length(out.type)
#     
#     msg <- paste0("must from: ", paste(id.types, collapse = ", "), "!")
#     if (!in.type %in% id.types) 
#       stop("'in.type' ", msg)
#     if (!all(out.type %in% id.types)) 
#       stop("'out.type' ", msg)
#     
#     in.ids0 <- in.ids
#     in.ids <- unique(as.character(in.ids))  #unique necessary for select()# if(unique.map)
#     out.ids <- character(length(in.ids))
#     res <- try(suppressWarnings(AnnotationDbi::select(db.obj, keys = in.ids, 
#                                                       keytype = in.type, columns = c(in.type, out.type))))
#     all.mappings <- rbind(all.mappings, res)
#   }
#   res <- as.data.frame(all.mappings[-1, ])
#   
#   res <- res[, c(in.type, out.type)]
#   
#   na.idx <- is.na(res[, 2])
#   if (sum(na.idx) > 0) {
#     n.na <- length(unique(res[na.idx, 1]))
#     if (n.na > 0) {
#       print(paste("Note:", n.na, "of", length(in.ids), "unique input IDs unmapped."))
#     }
#     if (na.rm) 
#       res <- res[!na.idx, ]
#   }
#   
#   cns <- colnames(res)
#   
#   res <- as.matrix(res)
#   rownames(res) <- NULL
#   return(res)
# }

#########################################################################################################
