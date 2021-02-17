library(testthat)

###################################################
test_that("changeDataId for compound", {
  
  cpd.sim.data <- sim.mol.data(mol.type = "cpd",
                               id.type = "KEGG COMPOUND accession",
                               nmol = 50000, 
                               nexp = 2)
  change.cpd.id <- changeDataId(data.input.id = cpd.sim.data,
                                input.type = "kegg.ligand",
                                output.type = "pathwayCommons",
                                mol.type = "cpd",
                                sum.method = "sum")
  
  expect_true(nrow(change.cpd.id) > 0)
  expect_true(ncol(cpd.sim.data) == ncol(change.cpd.id))
})

###################################################
test_that("downloadSbgn", {
  
  data("sbgn.xmls")
  files <- downloadSbgnFile(pathway.id = c("P00001", "P00002"))
  files <- gsub(".//", "", files)
  get.files <- list.files(path = ".", pattern = "*.sbgn")
  
  expect_identical(files, get.files)
})

###################################################
test_that("sbgn.gsets", {
  
  mol.list <- sbgn.gsets(database = "metacrop",
                         id.type = "ENZYME",
                         species = "ath",
                         output.pathway.name = FALSE,
                         truncate.name.length = 50)
  expect_gte(length(mol.list), 0)
})

###################################################
# test_that("changeIds", {
# 
#   gdata.bta <- sim.mol.data(mol.type = "gene", id.type = "ENSEMBLPROT",
#                             species = "bta", nmol = 2000)
#   ci.bta <- changeIds(input.ids = names(gdata.bta),
#                       input.type = "ENSEMBLPROT",
#                       output.type = "KO",
#                       mol.type = "gene",
#                       org = "bta")
# 
#   expect_true(class(ci.bta) == "list")
#   expect_true(length(ci.bta) == 2000)
# })
