
#########################################################################################################
#' Overlay omics data on SBGN pathway diagram and output image files.
#' 
#' This is the main function to map, integrate and render omics data on pathway graphs. 
#' Two inputs are needed: 1. A pathway file in SBGN-ML format and 2. gene and/or compound omics data. 
#' The function generates image file of a pathway graph with the omics data mapped to the glyphs and rendered as pseudo-colors. 
#' If no gene and/or compound omics data is supplied to the function, the function will output the SVG image file (and other selected file formats) of the parsed input file . This is useful for viewing the pathway graph without overlaid omics data. 
#' This function is similar to Pathview except the pathways are rendered with SBGN notation.
#' In addition, users can control more graph properties including node/edge attributes. We collected SBGN-ML files from several pathway databases: Reactome, MetaCyc, MetaCrop, PANTHER and SMPDB. Given a vector of patway IDs, SBGNview can automatically download and use these SBGN-ML files. To map omics data to glyphs, user just needs to specify omics data ID types. When using user customized SBGN-ML files, users need to provide a mapping file from omics data's molecule IDs to SBGN-ML file's glyph IDs.\cr
#' 
#' 
#' @param gene.data A matrix, vector or SummarizedExperiment object. The same as 'gene.data' argument in package 'pathview', it is either a vector (single measurement) or a matrix-like data (multiple measurements). If the data is a vector, the entries should be numeric and names of entries should be gene IDs. Matrix data structure has genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a generic concept, including multiple types: gene, transcript or protein. Default gene.data=NULL.
#' @param cpd.data A matrix, vector or SummarizedExperiment object. The same as 'gene.data', excpet named with compound IDs. Default cpd.data=NULL. 
#' @param simulate.data Logical. SBGNview can simulate a dataset. If set to TRUE, SBGNview will simulate a gene data set and a compound dataset and user input 'gene.data' and 'cpd.data' are ignored.
#' @param input.sbgn  A character vector. Can be either names of local SBGN files or pathway IDs of our pre-collected pathways. For pre-collected pathway IDs, run 'data(pathways.info)'
#' @param sbgn.dir  A character string. Default: ".". The path to the folder that holds SBGN-ML files. If 'input.sbgn' is a vector of pathway IDs in data 'pathways.info', the SBGN-ML files will be downloaded into this folder. 
#' @param org A character string. Default: "hsa". The species of the gene omics data. It is used for species specific gene ID mapping. Currently only supports three letters KEGG code (e.g. hsa, mmu, ath).  For a complete list of KEGG codes, see this page:\cr \href{https://www.genome.jp/kegg/catalog/org_list.html}{KEGG Organisms: Complete Genomes} 
#' @param output.formats   A character vector. It specifies the formats of output image files. The vector should be a subset of c('pdf' , 'ps', 'png'). By default the function will always output a svg file. SBGNview uses rsvg to convert svg file to other formats. If other 'output.formats' is set but 'rsvg' package is not installed, an error will occur. See this page for how to \href{https://github.com/jeroen/rsvg}{install 'rsvg'}
#' @param output.file   A character string. Default: "./output.svg". Path to the output image files. Because we often work with multiple pathways, each pathway will have its own image files. Each string in 'input.sbgn' will be added to the end of 'output.file'. Depending on the image format specified by the 'output.formats' parameter, extentions will be added to the end (e.g. .pdf, .png etc.).
#' @param gene.id.type  A character string. The type of gene ID in 'gene.data'. This parameter is used for ID mapping. It should be one of the IDs in data 'mapped.ids'. For details, run: \code{data('mapped.ids')}
#' @param cpd.id.type  A character string. The type of compound ID in 'cpd.data'.  For details, run: \code{data('mapped.ids')}
#' @param sbgn.id.attr  A character string. This tells SBGNview where to find the ID of a glyph in SBGN-ML file for ID mapping. This ID is used to map omics data to the glyph. It is normally the name of an attribute in the 'glyph' element . For example : <glyph class='macromolecule' id='p53'> </glyph>. We can specify: sbgn.id.attr = 'id'; sbgn.gene.id.type = 'SYMBOL'. For \href{https://github.com/datapplab/SBGN-ML.files/tree/master/data/SBGN.with.stamp}{our pre-generated SBGN-ML files}, the ID attribute will be determined automatically thus can be omitted. 
#'          Accepted values: 
#'              1. Any attribute name in element 'glyph' For example : <glyph class='macromolecule' id='p53' protein='P04637'> </glyph>. We can specify: sbgn.id.attr = 'protein'; sbgn.gene.id.type = 'UNIPROT', then 'P04637' will be the glyph ID. 
#'              2. The string 'label', this will make SBGNview use the glyph label as glyph ID. For example: <glyph id='glyph14' class='simple chemical'> <label text='L-alanine'/>  </glyph>. We can specify: sbgn.id.attr = 'label'; sbgn.cpd.id.type = 'compound.name', then 'L-alanine' will be used as glyph ID.
#' @param sbgn.gene.id.type  A character string.  The ID type of "macromolecule" glyphs in SBGN-ML file (See parameter 'sbgn.id.attr' for more details). This parameter is used for ID mapping, i.e. either use our pre-generated mapping tables or find corresponding columns in user defined mapping tables in 'id.mapping.gene'.  For \href{https://github.com/datapplab/SBGN-ML.files/tree/master/data/SBGN.with.stamp}{our pre-generated SBGN-ML files}, this will be determined automatically according to the pathway IDs thus can be omitted. For user defined SBGN-ML file, this parameter should be one of the column names of the matrix 'id.mapping.gene'.
#' @param sbgn.cpd.id.type   A character string. Similar to 'sbgn.gene.id.type'. The corresponding glyphs are "simple chemicals"
#' @param id.mapping.gene A matrix.  Mapping table between gene.id.type and sbgn.gene.id.type. This table is needed if the ID pair of gene.id.type and sbgn.gene.id.type is NOT included in data 'mapped.ids' or not mappable by package 'pathview'. This matrix should have two columns for gene.id.type and sbgn.gene.id.type, respectively.  Column names should be the values of parameters 'sbgn.gene.id.type' and 'gene.id.type'.  See example section for an example.
#' @param id.mapping.cpd A matrix. See id.mapping.gene.
#' @param node.sum  A character string. Default: "sum". Sometimes multiple omics genes/compounds are mapped to one SBGN glyph. Therefore multiple values will be mapped to one measurement/slice on the glyph. In this situation, we may need to derive a single value for the slice on the glyph. This function can be any R function that takes a numeric vector as input and output a single numeric value (e.g. 'sum','max','min','mean'). It can also be a User Defined Function (UDF).
#' @param pathway.name A character string. Change/update pathway name displayed on the output graph. If 'input.sbgn' is a pathway ID in data(pathways.info), the pathway name and database associated with the pathway ID will be displayed on the output graph. If 'input.sbgn' is a SBGN-ML file not part of our pre-generated SBGN-ML files, nothing will be dispalyed for the pathway name unless set using this arugmnet. 
#' @param show.pathway.name Logical. Default: F. If set to TRUE and 'input.sbgn' are pre-collected pathway IDs, the pathway name will be added to the output file name.
#' @param SBGNview.data.folder A character string. Default: "./SBGNview.tmp.data". The path to a folder that will hold temp data files.            
#' @param ... Other parameters passed to function \code{\link{renderSbgn}}
#' @return  SBGNview object (S3 class object) which is a list containing three elements: 
#' 
#'      1. data: 
#       A list contains information for all input SBGN pathways. Each element is a list (e.g. 'P00001') holding information of a single SBGN pathway. The information include parsed information of glyphs, arcs and SBGNview parameters passed to SBGNview function
#'
#'      2. output.file:
#'      A string of the path to the output file. It is the string set by parameter 'output.file' in function SBGNview.     
#'       
#'      3. output.formats:
#'      A character vector specifying the formats of output image files. The vector should be a subset of c('pdf' , 'ps', 'png'). By default the function will always output a svg file.
#'      
#'      This S3 class of objects contains information generated by the SBGNview function. 
#'      The SBGNview object can be modified by using the "+" binary operator in conjunction with the \code{\link{highlightArcs}}, \code{\link{highlightNodes}}, and \code{\link{highlightPath}} functions. 
#'      This mechanism to allow layer-by-layer graph modification is similar to that of ggplot2, a widely used R package for data visualization. 
#'      Defining this class as a formal (S4) class will produce an S4 SBGNview object which doesn't work with the binary operator and doesn't allow for layer-by-layer graph modification. 
#'      Therefore, we decided to go with the S3 implementation of the SBGNview class since it works well with the binary operator making this functionality more intuitive and user friendly.     
#'                 
#' @details 
#'          1. About SBGNview()
#' 
#'          This function extracts glyph (node) and arc (edge) data from a SBGN-ML file and creates a SBGN graph from the extracted data (draws shapes etc. in SVG format). Then it maps omics data to the glyphs and renders data as colors. Currently it maps gene/protein omics data to 'macromolecule' glyphs and maps compound omics data to 'simple chemical' glyphs.
#'          
#'          2. About SBGN-ML files and curved arcs encoding
#'          
#'            SBGNview can parse SBGN-ML files with standard SBGN PD syntax. For arcs, SBGNview can handle both straight lines (standard SBGN PD syntax) and spline curves (a syntax add by us). Current SBGN-ML syntax supports straight lines. The coordinates of line start or end points are stored in element 'arc'. But straight lines often produce node-edge or edge-edge crossings. Therefore, we generated SBGN-ML files with pre-routed spline edges. 
#'            
#'            To store the routed splines, we added an XHTML element called 'edge.spline.info', which has children elements called 'arc.spline' . Each 'arc.spline' element has three types of children elements: 'source.arc', 'target.arc' and 'spline'. 'source.arc' and 'target.arc' will be rendered as straight line arcs controlled by attributes 'start.x','start.y',  'end.x', 'end.y' (line ends' coordinates) and 'class' (type of the straight line arc). These two arcs ensure the notation of the spline arc comply with its class. 'spline' will be rendered as splines connecting 'source.arc' and 'target.arc'.  Each 'spline' is defined by coordinates of four points: s (starting point of spline), c1 (the first control point of spline), c2 (the second control point of spline) and e (ending point of spline). In case of complicated routing, there could be multiple 'splines' in an 'arc.spline'. 
#'           
#'          The function first checks if the SBGN-ML file has spline arcs (XHTML element 'edge.spline.info')  and use it if found. When there are no spline arcs, it will use straight line arcs (XHTML element 'arc').  Please check out examples in  \href{ https://github.com/datapplab/SBGN-ML.files/tree/master/data/SBGN}{our SBGN-ML file collection.}
#'          
#'          3. About ID mapping
#'          
#'           SBGNview can automatically map several ID types to glyphs of pathwayCommons, MetaCyc and MetaCrop. For user defined SBGN-ML file, users need to provide information about how to map their omics data to glyphs. 
#'           
#'            3.1 How SBGNview finds glyph IDs in SBGN-ML file: Glyph IDs are recorded in attribute 'id' in XHTML element 'glyph'. But for ID mapping, user can use other attributes by changing parameter 'sbgn.id.attr'.   
#'           
#'            3.2 For \href{ https://github.com/datapplab/SBGN-ML.files/tree/master/data/SBGN}{our SBGN-ML file collection.}, SBGNview can do ID mapping automatically. It uses extracted mapping between 1) UNIPROT/uniref  and 'macromolecule' glyph IDs and 2) ChEBI and 'simple chemical' glyph IDs from biopax files in pathwayCommons and MetaCyc. For other ID types, we used pathview (gene/protein) and UniChem (compound) to map to UNIPROT and ChEBI, respectively, then map them to glyph IDs. For MetaCrop, we used pathview for ID mapping.
#'          
#'          
#'          4. Two common scenarios of using SBGNview
#'            
#'            4.1 Using \href{https://github.com/datapplab/SBGN-ML.files/tree/master/data/SBGN}{our pre-generated SBGN-ML files}.
#'            
#'             Supported pathways can be found using \code{data('pathways.info')}. This is a collection of SBGN-ML files for these databases: MetaCyc, MetaCrop and three databases collected by pathwayCommons (PANTHER Pathway, Reactome and SMPDB). For SBGN-ML each file, the glyph layout is based on fdp and optimized to eliminate glyph-glyph overlaps. The arcs are splines that are routed to eliminate arc-glyph crossings. 
#'             
#'          To use these data, SBGNview needs the following parameters: 
#'          
#'          -gene.id.type and/or cpd.id.type (at least one should be provided)
#'          
#'          
#'          SBGNview can map omics data to SBGN-ML glyphs automatically. Supported ID types can be found in \code{data('mapped.ids')}
#'          
#'          Input SBGN-ML files can be obtained by using function 'downloadSbgnFile'. '          
#'          
#'          
#'            4.2 Using SBGN-ML files from other sources.
#'            
#'          
#'              4.2.1 Input omics data have the SAME ID type as the glyph ID type in SBGN-ML file:
#'              
#'          In this scenario, SBGNview needs the following information to map omics data to SBGN-ML glyphs: 
#'          
#'          -ID type of input omics data (gene.id.type and/or cpd.id.type)
#'          
#'          -ID type of glyphs of input SBGN-ML file (sbgn.gene.id.type and/or sbgn.cpd.id.type).
#'          
#'          These ID types can be any characters, but gene.id.type must be the same as sbgn.gene.id.type, and cpd.id.type must be the same as sbgn.cpd.id.type.
#'          
#'          Users can use the function 'changeDataId' to change the omics IDs to the glyph IDs in SBGN-ML file. 
#'          
#'              4.2.2 Input omics data have DIFFERENT ID type as the glyph ID type in SBGN-ML file:
#'              
#'          In this scenario, SBGNview needs the following information to map omics data to SBGN-ML glyphs: 
#'          
#'          -ID type of input omics data (gene.id.type and/or cpd.id.type)
#'          
#'          -ID type of glyphs of input SBGN-ML file (sbgn.gene.id.type and/or sbgn.cpd.id.type).
#'          
#'          -A mapping table between input omics IDs and SBGN-ML glyph IDs (id.mapping.gene and/or id.mapping.cpd).
#'          
#'          For user's convinience, pathview can generate such tables for several ID types (functions 'geneannot.map' and 'cpdidmap'). But column names need to be changed to the values of 'gene.id.type' and 'sbgn.gene.id.type'.
#' 
#' @examples 
#' ### Use simulated data. Please see vignettes for more examples.
#' ### Run `browseVignettes(package = "SBGNview")` 
#' 
#' # load demo dataset, SBGN pathway data collection and info, which may take a few seconds
#' # we use a cancer microarray dataset 'gse16837.d' from package 'pathview'
#' data("gse16873.d", "pathways.info", "sbgn.xmls")
#' 
#' # search for pathways with user defined keywords
#' input.pathways <- findPathways("Adrenaline and noradrenaline biosynthesis")
#' 
#' # render SBGN pathway graph and output image files
#' SBGNview.obj <- SBGNview(gene.data = gse16873.d[,1:3],
#'                          gene.id.type = "entrez",
#'                          input.sbgn = input.pathways$pathway.id,
#'                          output.file = "quick.start",
#'                          output.formats = c("png", "pdf"),
#'                          min.gene.value = -1,
#'                          max.gene.value = 1) 
#'  SBGNview.obj           
#' 
#' @export

SBGNview <- function(gene.data = NULL, cpd.data = NULL, simulate.data = FALSE, input.sbgn = NULL, 
                     sbgn.dir = "./", output.file = "./output.svg", node.sum = "sum", gene.id.type = NA, 
                     cpd.id.type = NA, sbgn.id.attr = "id", sbgn.gene.id.type = NULL, sbgn.cpd.id.type = NA, 
                     id.mapping.gene = NULL, id.mapping.cpd = NULL, org = "hsa", output.formats = c("svg"), 
                     pathway.name = NULL, show.pathway.name = FALSE, 
                     SBGNview.data.folder = "./SBGNview.tmp.data", ...) {
    
    if (!dir.exists(sbgn.dir)) {
        warning("'sbgn.dir' folder does not exist! Creating: ", sbgn.dir, "\n")
        system(paste("mkdir -p", sbgn.dir))
    }
    # using changed mapping file names from CompoundName and compound.name and kegg.ligand to kegg
    if(!is.na(sbgn.cpd.id.type)) {
        if(sbgn.cpd.id.type == "CompoundName") sbgn.cpd.id.type = "compound.name"
        if(sbgn.cpd.id.type %in% c("kegg.ligand", "KEGG")) sbgn.cpd.id.type = "kegg"
    }
    
    # Parse all files in a folder when no input.sbgn is specified.
    if (is.null(input.sbgn)) {
        warning("'input.sbgn' is not provided, using all files in folder 'sbgn.dir':", 
                sbgn.dir, "\n Please make sure all files in ", sbgn.dir, " are SBGM-ML files.\n")
        input.sbgn <- list.files(sbgn.dir, full.names = FALSE)
        if (length(input.sbgn) == 0) {
            stop("Must provide 'input.sbgn' if 'sbgn.dir':", sbgn.dir, " is empty! ")
        }
    }
    input.sbgn <- unique(as.character(input.sbgn))

    user.data.recorded <- list()  # SBGNview can render multiple SBGN-ML files but the input omics data only need to be processed once. Here we record the converted data.
    SBGNview.obj.data <- list()
    
    # parse all SBGN-ML files
    for (i in seq_along(input.sbgn)) {
        input.sbgn.i <- input.sbgn[i]
        message("\n\nProcessing pathway: ", input.sbgn.i, "\n")
        # Find related information according to input sbgn
        parse.sbgn.list <- parse.input.sbgn(input.sbgn.i, output.file, show.pathway.name, 
                                            sbgn.dir, sbgn.gene.id.type, sbgn.cpd.id.type, 
                                            sbgn.id.attr, SBGNview.data.folder)
        input.sbgn.full.path <- parse.sbgn.list$input.sbgn.full.path
        output.file.sbgn <- parse.sbgn.list$output.file.sbgn
        sbgn.gene.id.type <- parse.sbgn.list$sbgn.gene.id.type
        sbgn.cpd.id.type <- parse.sbgn.list$sbgn.cpd.id.type
        pathway.name.on.graph <- parse.sbgn.list$pathway.name.on.graph
        if.file.in.collection <- parse.sbgn.list$if.file.in.collection
        sbgn.id.attr <- parse.sbgn.list$sbgn.id.attr
        database <- parse.sbgn.list$database
        
        if(!is.null(pathway.name)) {
            pathway.name.on.graph$pathway.name.on.graph <- paste(pathway.name, "user.named.pathway", sep = "::")
        }
        
        # Change omics data ID to SBGN-ML glyph ID. Input pathways may come from
        # different databases, so they may have different glyph ID types (e.g.
        # pathwayCommons, MetaCyc or EC number). But if two SBGN-ML files have the same
        # glyph type, we only need to do ID mapping once, record converted data in
        # 'user.data.recorded', then just reuse it next time we need the same glyph ID
        # type.
        parsed.data <- parse.omics.data(gene.data, cpd.data, input.sbgn.full.path, 
                                        database, user.data.recorded, gene.id.type, cpd.id.type, 
                                        id.mapping.gene, id.mapping.cpd, node.sum, org, sbgn.gene.id.type, 
                                        sbgn.cpd.id.type, simulate.data, SBGNview.data.folder = SBGNview.data.folder)
        user.data <- parsed.data$user.data  # A list holding both gene data and/or cpd data. The names of this list are glyph IDs (for 'macromolecule' and 'simple chemical') in 'input.sbgn'. Each element of 'user.data' is a vector containing omics data of the corresponding gene/macromolecule or compound/simple chemical. 
        user.data.recorded <- parsed.data$user.data.recorded  # Update user data list, recording data with different glyph ID types. The names of this list are glyph ID types. Each element is the list 'user.data'. It was updated if 'input.sbgn' has a new glyph ID type: a new 'user.data' with this new glyph ID type is added.
        
        # extract arc information and compartment layer information
        message("reading SBGN-ML file for arcs info: ", input.sbgn.full.path)
        sbgn.xml <- read_xml(input.sbgn.full.path)
        xml_attrs(sbgn.xml) <- NULL  # Remove root glyph attribute. This is necessary. Otherwise xml2 won't find the glyphs with xml_find_all() function.
        arcs.info <- get.arcs.info(sbgn.xml)
        # Retrieve compartment layer information. When compartments overlap, this will
        # ensure they are plotted in correct layers so nodes visible are always from
        # upper layer (i.e. no nodes are covered by other compartments)
        compartment.layer.info <- get.compartment.layer(sbgn.xml)
        
        message("Rendering SBGN graph")
        sbgn.result.list <- renderSbgn(input.sbgn = input.sbgn.full.path, output.file = output.file.sbgn, 
                                       arcs.info = arcs.info, compartment.layer.info = compartment.layer.info, 
                                       user.data = user.data, output.formats = output.formats, 
                                       sbgn.id.attr = sbgn.id.attr, pathway.name = pathway.name.on.graph, 
                                       if.plot.svg = FALSE, ...)
        # record all parameters. They might be used again when we later modify the 'SBGNview' object
        sbgn.result.list[["render.sbgn.parameters.list"]] <- list(input.sbgn = input.sbgn.full.path, 
                                                                  output.file = output.file.sbgn,
                                                                  arcs.info = arcs.info, 
                                                                  compartment.layer.info = compartment.layer.info, 
                                                                  user.data = user.data, sbgn.id.attr = sbgn.id.attr, 
                                                                  pathway.name = pathway.name.on.graph)
        SBGNview.obj.data[[input.sbgn.i]] <- sbgn.result.list
    }
    
    SBGNview.obj <- createSBGNviewObject(data = SBGNview.obj.data, output.file = output.file,
                                         output.formats = output.formats)
    message("SBGNview object generated")
    return(SBGNview.obj)
}

#########################################################################################################
### this function is for internal use only
# SBGNview S3 class
# 
# This S3 class of objects contains information generated by the SBGNview function.
# The SBGNview object can be modified by using the "+" binary operator in conjunction with the \code{\link{highlightArcs}}, \code{\link{highlightNodes}}, and \code{\link{highlightPath}} functions. 
# This mechanism to allow layer-by-layer graph modification is similar to that of ggplot2, a widely used R package for data visualization. 
# Defining this class as a formal (S4) class will produce an S4 SBGNview object which doesn't work with the binary operator and doesn't allow for layer-by-layer graph modification. 
# Therefore, we decided to go with the S3 implementation of the SBGNview class since it works well with the binary operator making this 
# functionality more intuitive and user friendly. 
# 
# @format A list containing three elements:
#           1. data: 
#              A list contains information for all input SBGN pathways. Each element is a list (e.g. 'P00001') holding information of a single SBGN pathway. The information include parsed information of glyphs, arcs and SBGNview parameters passed to SBGNview function
#              
#           2. output.file:
#              A string of the path to the output file. It is the string set by parameter 'output.file' in function SBGNview.
#              
#           3. output.formats:
#              A character vector specifying the formats of output image files. The vector should be a subset of c('pdf' , 'ps', 'png'). By default the function will always output a svg file.
# @examples 

createSBGNviewObject <- function(data, output.file, output.formats){
    structure(list(data = data, output.file = output.file, 
                   output.formats = output.formats), class = "SBGNview")
}

#########################################################################################################
