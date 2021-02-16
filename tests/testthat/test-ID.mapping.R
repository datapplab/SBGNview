library(SBGNview)
library(testthat)
# context("test-ID.mapping") # context() no longer recommended

test.id.mapping <- function(species, input.types, output.types, 
                           example.ids, mol.type) {
    n = 0
    for(i in seq_len(length(input.types))) {
       input.type = input.types[i]
       for(j in seq_len(length(output.types))) {
           output.type = output.types[j]
           for(k in seq_len(length(species))) {
               species.run = species[k]
               if(identical(input.type,output.type)) {
                  next()
               }
               n = n+1
               message("\n\nChecking: ",n," ",mol.type," ",input.type," ",output.type," ",species.run,"\n")
               id.map <- try(loadMappingTable(output.type = output.type,
                                             input.type = input.type,
                                             mol.type = mol.type,
                                             species = species.run,
                                             limit.to.ids = example.ids[[input.type]][[species.run]]
                                             )
                            , silent = TRUE)
               # expect_that(nrow(id.map) > 0, expect_true()) # old style of testing that's no longer encouraged.
               expect_true(nrow(id.map) > 0)
           }
       }
    }
  return(invisible(0))
}

example.ids <- list(ENTREZID = list(mmu = c(26395,14693), hsa = c(29085,10993)),
                    ENSEMBL = list(hsa = c("ENSG00000156006","ENSG00000244593"),
                                   mmu = c("ENSMUSG00000002997", "ENSMUSG00000004455")),
                    KO = list(hsa = c("K01847", "K00921"), mmu = c("K01847", "K00921")),
                    chebi = list(cpd = c("10036","10049","18420","48828")),
                    CompoundName = list(cpd = c("tyrosine", "(+-)-epinephrine",
                                                "1,3, 7-Trimethyluric acid",
                                                "(1S)-3,7,7-trimethylbicyclo[4.1.0]hept-3-ene")),
                    kegg.ligand = list(cpd = c("C00451","C00186","C11382","C06304"))
                    )

test.gene <- function() {
    # species = c("hsa","mmu")
    # input.types = c("ENTREZID","KO","ENSEMBL")
    # output.types = c("pathwayCommons","metacyc.SBGN","pathway.id")
    species = c("hsa")
    input.types = c("ENTREZID")
    output.types = c("pathwayCommons")
    mol.type = "gene"
    test.id.mapping(species = species,
                    input.types  = input.types,
                    output.types = output.types,
                    mol.type = mol.type,
                    example.ids = example.ids)
}

test.compound <- function() {
    species = c("cpd")
    input.types = c("chebi","CompoundName","kegg.ligand")
    output.types = c("pathwayCommons","metacyc.SBGN","pathway.id")
    mol.type = "cpd"
    test.id.mapping(species  = species,
                    input.types  = input.types,
                    output.types = output.types,
                    mol.type = mol.type,
                    example.ids = example.ids)
}


test.gene()
test.compound()
