
# SBGNview 
SBGNview is a tool set for visualizing omics data on SBGN pathway maps.  Given omics data and a SBGN-ML file with layout information, SBGNview can display omics data as colors on glyphs and output image files. SBGNview provides extensive options to control glyph and edge features  (e.g. color, line width etc.). To facilitate pathway based analysis, SBGNview also provides functions to extract glyph information and pairwise interactions from SBGN-ML files. SBGNview can map a large collection of gene, protein and compound ID typs to glyphs.  


# Introduction
Molecular pathways have been widely used in omics data analysis. We previously developed an R/BioConductor package called Pathview, which maps, integrates and visualizes omics data onto KEGG pathway graphs. Since its publication, Pathview has been widely used in numerous omics studies and analysis tools. Here we introduce the SBGNview package, which adopts Systems Biology Graphical Notation (SBGN)[@le2009systems] and greatly extends the Pathview project by supporting multiple major pathway databases besides KEGG.

Key features:

* Pathway diagram is drawn with SBGN notations 

* Supports major pathway databases and user defined pathways. 

* Extensive choices for graphical control. 

* Pathway related data extraction and analysis

# Installation

## Prerequisites
*SBGNview* depends or imports from the following R packages. Note these dependencies will be automatically installed when SBGNview is installed from BioConductor or GitHub

* xml2: parse SBGN-ML files
* rsvg: convert svg files to other formats (pdf, png, ps). librsvg2 is needed to install rsvg. See this page for more details: https://github.com/jeroen/rsvg
* igraph: find shortest paths
* httr: search [SBGNhub](https://github.com/datapplab/SBGNhub/tree/master/data/id.mapping.unique.pair.name) for mapping files
* KEGGREST: generate mapping tables from scratch when needed
* [pathview](https://bioconductor.org/packages/release/bioc/html/pathview.html): map between different ID types for gene and chemical compound
* [gage](https://bioconductor.org/packages/release/bioc/html/gage.html): R package for pathway enrichment analysis.
* [SBGNview.data](https://bioconductor.org/packages/release/data/experiment/html/SBGNview.data.html): demo gene expression datasets for SBGNview package
* [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html): main function accepts data as SummarizedExperiment objects
* AnnotationDbi: filter and select data for mapping between IDs

```{r setup, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)){
     install.packages("BiocManager")
}
BiocManager::install(c("xml2", "rsvg", "igraph", "httr", "KEGGREST", "pathview", "gage", "SBGNview.data", "SummarizedExperiment", "AnnotationDbi"))
```

Dependencies for Operating Systems:  
**Windows 10**: none

**Linux (Ubuntu)**: needs additional packages (libxml2-dev, libssl-dev, libcurl4-openssl-dev, librsvg2-dev) to be installed. Run the command below in a terminal to install the necessary packages. The same or similar packages can be found for other distributes of linux.
```{r depend, eval = FALSE}
sudo apt install libxml2-dev libssl-dev libcurl4-openssl-dev librsvg2-dev
```

## Install SBGNview
Install **SBGNview** through Bioconductor: 
```{r install, eval = FALSE}
BiocManager::install(c("SBGNview"))
```
Install **SBGNview** through GitHub:
```{r install.1, eval = FALSE}
install.packages("devtools")
devtools::install_github("datapplab/SBGNview")
```
Clone the Git repository:
```{r clone.git, eval = FALSE}
git clone https://github.com/datapplab/SBGNview.git
```


## Quick example
```{r, echo = TRUE, eval = TRUE, results = 'hide', message = FALSE, warning = FALSE}
# load demo dataset and pathway information of built-in collection of SBGN-ML files
library(SBGNview)
data("gse16873.d","pathways.info")
input.pathways <- findPathways("Adrenaline and noradrenaline biosynthesis")
SBGNview.obj <- SBGNview(
          gene.data = gse16873.d[,1:3], 
          gene.id.type = "entrez",
          input.sbgn = input.pathways$pathway.id,
          output.file = "quick.start", 
          output.formats =  c("png")
          ) 
print(SBGNview.obj)
```
Two image files (a .svg file by default and a .png file) will be created in the current working directory.
<img src="inst/app/www/quick.start_P00001.svg">   

# SVG Editing Tool
SBGNview provides extensive graphical control options for adjusting glyph/arc attributes such as line color/width and text size/color/wrapping/positioning. The SBGNview main function outputs SVG file format by default and other (png, pdf, ps) specified file formats. There might be cases where you would like to adjust specific graphical details more easily. In such cases, Inkscape can be used to modify the generated SVG file and save/export the SVG to other file formats provided by Inkscape. [Inkscape](https://inkscape.org/about/) is a free, open source, and easy to use vector graphics editor (like Illustrator) which is available on Linux, Windows and MacOS. Download Inkscape from [here](https://inkscape.org/release/inkscape-1.0.1/). 

# Additional information
SBGN website: https://sbgn.github.io/

For any questions, please contact Xiaoxi Dong(<dfdongxiaoxi@gmail.com>) or Weijun Luo(<Weijun.Luo@uncc.edu>)

# Citation
Nicolas Le Novère, Michael Hucka, Huaiyu Mi, Stuart Moodie, Falk Schreiber, Anatoly Sorokin, Emek Demir et al. The Systems Biology Graphical Notation Nature Biotechnology, 27(8):735-741, 2009

Luo, Weijun, and Cory Brouwer. "Pathview: an R/Bioconductor package for pathway-based data integration and visualization." Bioinformatics 29.14 (2013): 1830-1831.