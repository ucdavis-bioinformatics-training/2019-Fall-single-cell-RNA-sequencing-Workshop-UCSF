### Alignment lecture by [Dr. Gerald Quon](https://qlab.faculty.ucdavis.edu/)

Install Python 3.6.8 from here:

[https://www.python.org/downloads/release/python-368/](https://www.python.org/downloads/release/python-368/)

Set some options and make sure the packages cowplot, circlize, tensorflow, scAlign are installed (if not install it), and then load them and verify they all loaded correctly.

In the R console run the following commands
```r
if (!any(rownames(installed.packages()) == "cowplot")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("cowplot")
}
library(cowplot)

if (!any(rownames(installed.packages()) == "circlize")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("circlize")
}
library(circlize)

if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)

## if you have conda installed you may need to install tensorflow and initialize your environment fist.
## conda install tensorflow
## reticulate::use_condaenv("anaconda3",required = TRUE)
## reticulate::py_discover_config("tensorflow")

if (!any(rownames(installed.packages()) == "tensorflow")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("tensorflow")
}
library(tensorflow)

## installs CPU version of TensorFlow Python package, use this if you don't have a GPU.
install_tensorflow(version='1.15rc2')

#install scAlign
if (!any(rownames(installed.packages()) == "scAlign")){
  if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
  devtools::install_github(repo = 'quon-titative-biology/scAlign')
}
library(scAlign)

sessionInfo()
```

### Download the template Markdown workshop document PART1 and open it.

In the R console run the following command
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-Fall-single-cell-RNA-sequencing-Workshop-UCSF/master/scrnaseq_analysis/scRNA_Workshop-alignment.Rmd", "scRNA_Workshop-alignment.Rmd")
```


### Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "Single Cell RNAseq Alignment"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
---</code></pre>
