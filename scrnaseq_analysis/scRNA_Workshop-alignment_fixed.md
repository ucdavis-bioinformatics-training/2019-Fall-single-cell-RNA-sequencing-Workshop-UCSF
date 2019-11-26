---
title: "Single Cell RNAseq Part 7 - Alignment"
author: "Gerald Quon"
output:
    html_document:
      keep_md: TRUE
---



 Overview

In this section, we will learn how to use the method [scAlign](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1766-4) to take two separate datasets and "integrate" them, so that cells of the same type (across datasets) roughly fall into the same region of the scatterplots (instead of separating by dataset first). Integration is typically done in a few different scenarios, e.g., 1) if you collect data from across multiple conditions / days / batches / experimentalists / etc. and you want to remove these technical confounders, 2) if you are doing a case control study (as we are here) and you want to identify which cells match across condition, or 3) you have performed an experiment sequencing cells from a tissue (e.g. lung epithelium) and you want to label the cells by type, but you don't have marker genes available, however, you do have access to a database of annotated cells that you could map onto your dataset.

Here we will perform alignment as if we do not have any labels (case 3), but we will use the labels after alignment to check its accuracy. The following R markdown illustrates how to do integration with scAlign, and aligns two datasets pretty successfully.

This tutorial provides a guided alignment for two groups of cells from [Kowalczyk et al, 2015](https://www.ncbi.nlm.nih.gov/pubmed/26430063). In this experiment, single cell RNA (scRNA) sequencing profiles were generated from HSCs in young (2-3 mo) and old (>22 mo) C57BL/6 mice. Age related expression programs make a joint analysis of three isolated cell types long-term (LT), short-term (ST) and multipotent progenitors (MPPs). In this tutorial we demonstrate the unsupervised alignment strategy of scAlign described in Johansen et al, 2018 along with typical analysis utilizing the aligned dataset, and show how scAlign can identify and match cell types across age without using the labels as input.

 Alignment goals
The following is a walkthrough of a typical alignment problem for scAlign and has been designed to provide an overview of data preprocessing, alignment and finally analysis in our joint embedding space. Here, our primary goals include:

1. Learning a low-dimensional cell state space in which cells group by function and type, regardless of condition (age).
2. Computing a single cell paired differential expression map from paired cell projections.

First, we perform a typical scRNA preprocessing step using the Seurat package. Then, reduce to the top 3,000 highly variable genes from both datasets to improve convergence and reduce computation time.



```r
 User paths
working.dir = "." #where our data file, kowalcyzk_gene_counts.rda is located
results.dir = "." #where the output should be stored

 Load in data - either option
download.file('https://github.com/quon-titative-biology/examples/raw/master/scAlign_paired_alignment/kowalcyzk_gene_counts.rda','kowalcyzk_gene_counts.rda'); load('kowalcyzk_gene_counts.rda')
#load(url('https://github.com/quon-titative-biology/examples/raw/master/scAlign_paired_alignment/kowalcyzk_gene_counts.rda'))

 Extract age and cell type labels
cell_age = unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 1))
cell_type = gsub('HSC', '', unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 2)))

 Separate young and old data
young_data = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "young")]
old_data   = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "old")]

 Set up young mouse Seurat object
youngMouseSeuratObj <- CreateSeuratObject(counts = young_data, project = "MOUSE_AGE", min.cells = 0)
youngMouseSeuratObj <- NormalizeData(youngMouseSeuratObj)
youngMouseSeuratObj <- ScaleData(youngMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)
```

<div class='r_output'> Warning: The following arguments are not used: display.progress
</div>
<div class='r_output'> Suggested parameter: verbose instead of display.progress
</div>
<div class='r_output'> Centering and scaling data matrix
</div>
```r
 Set up old mouse Seurat object
oldMouseSeuratObj <- CreateSeuratObject(counts = old_data, project = "MOUSE_AGE", min.cells = 0)
oldMouseSeuratObj <- NormalizeData(oldMouseSeuratObj)
oldMouseSeuratObj <- ScaleData(oldMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)
```

<div class='r_output'> Warning: The following arguments are not used: display.progress
</div>
<div class='r_output'> Suggested parameter: verbose instead of display.progress

 Centering and scaling data matrix
</div>
```r
 Gene selection
youngMouseSeuratObj <- FindVariableFeatures(youngMouseSeuratObj, do.plot = F, nFeature=3000)
```

<div class='r_output'> Warning: The following arguments are not used: do.plot, nFeature
</div>
```r
oldMouseSeuratObj <- FindVariableFeatures(oldMouseSeuratObj, do.plot = F, nFeature=3000,)
```

<div class='r_output'> Warning: The following arguments are not used: do.plot, nFeature
</div>
```r
genes.use = Reduce(intersect, list(VariableFeatures(youngMouseSeuratObj),
                                   VariableFeatures(oldMouseSeuratObj),
                                   rownames(youngMouseSeuratObj),
                                   rownames(oldMouseSeuratObj)))

#[FUTURE] ComBat and Seurat (regress) batch correction results.
```

 Is alignment necessary?

Let's visualize the data from both young and old mice, to determine whether we need alignment/batch correction.


```r
 Combine our Seurat objects
hsc.combined <- merge(youngMouseSeuratObj, oldMouseSeuratObj, add.cell.ids = c("YOUNG", "OLD"), project = "HSC")
hsc.combined <- ScaleData(hsc.combined, do.scale=T, do.center=T, display.progress = T)
```

<div class='r_output'> Warning: The following arguments are not used: display.progress
</div>
<div class='r_output'> Suggested parameter: verbose instead of display.progress
</div>
<div class='r_output'> Centering and scaling data matrix
</div>
```r
VariableFeatures(hsc.combined) = genes.use

 Run PCA and TSNE
hsc.combined = RunPCA(hsc.combined, do.print=FALSE)
```

<div class='r_output'> PC_ 1
 Positive:  Mpl, Pdzk1ip1, Mmrn1, Ifitm1, Txnip, Car2, uc008aea.1,uc008aeb.1,uc008aec.2,uc012ajd.1, Shisa5, Cd63, uc008ewj.2
 	   uc012hdk.1, Zfp36, Sult1a1, Tbxas1, uc008xuk.1, Ifitm3, uc008ypk.1, Tgm2, H2-T23, Nupr1
 	   Clec1a, Procr, Uba7, Apoe, Ppp1r15a, Sat1, Tsc22d3, Cd74, Krt18, Pygm
 Negative:  Top2a, Pbk, Birc5, Mpo, Aurkb, Stmn1, uc009qeb.1, Cks1b, Cd48, Tk1
 	   Mcm5, Cdca8, Plac8, Kif11, Spc25, Ppia, Nusap1, H2afz, Kif22, Rrm2
 	   Tacc3, Prc1, Kif20a, Bub1b, Cdca3, Nkg7, Mcm7, Ckap2l, Hmgn2, Asf1b
 PC_ 2
 Positive:  Flt3, Cd34, Plac8, Wfdc17, Lat2, Ighv1-30, Gpr97, Serpinb1a, Emb, Tespa1
 	   Cd52, Tyrobp, Dntt, Il12a, Ccl3, Sell, Trf, Ncf1, Lsp1, H2-Ob
 	   uc008hln.1, Mpo, Bex6, Spns3, Cd48, uc007hxv.1, Haao, Cd33, Ubash3a, Ifi27l1
 Negative:  Nusap1, Prc1, Birc5, Cdca8, Fam64a, Ckap2l, Hmmr, Kif22, Ube2c, Kif20a
 	   Cdca3, Aurkb, Kif11, Bub1b, Mmrn1, Tpx2, Ccnb1, Ccnb2, Apoe, Pbk
 	   Shcbp1, Kif2c, Nupr1, Cenpa, Aurka, Spc25, Cdkn3, Cdc20, Spag5, Tubb5
 PC_ 3
 Positive:  Mcm2, Ldha, Mif, Cdca7, Mcm5, Mcm6, Tipin, Lig1, Mcm3, Uhrf1
 	   Slc25a5, Slbp, Mcm7, Pcna, Gins2, Mt1, Cd63, Fen1, Ppia, Myc
 	   Tubb5, Sdf2l1, Apoe, Hspa8, Itga2b, Gmnn, Rfc4, Tuba1b, Ppil1, Clu
 Negative:  Cdkn3, Ccnb2, Hmmr, Cenpa, Kif23, Ube2c, Kif2c, Kif20a, D2Ertd750e, Gimap6
 	   Ccnb1, Tpx2, uc007hxv.1, Fam64a, Mgst1, Cdc20, Ckap2l, Aurka, Cdca3, Nusap1
 	   Troap, Kif11, Cdc25b, Prc1, Ect2, D17H6S56E-5, Flt3, Cdca8, Psrc1, Kif22
 PC_ 4
 Positive:  Hba-a2, Hba-a1, Hbb-b2, Ebi3, Shisa5, Ly6a, Lig1, Tuba1b, Mcm3, Gins2
 	   Gbp6, Tk1, Uhrf1, Stat1, Mcm5, Tcf19, Mpl, Gbp9, Pcna, Ldhb
 	   Car2, Mcm7, Mcm6, Tipin, Rfc3, Tnfrsf18, Mamdc2, Gins1, Iigp1, Prim2
 Negative:  uc012hdk.1, Dhrs3, uc008ewj.2, uc008aea.1,uc008aeb.1,uc008aec.2,uc012ajd.1, Fos, Gpr97, uc008ypk.1, Prtn3, Zfp36, Lat2
 	   Mpo, Cenpa, Nupr1, uc008hln.1, Dntt, Sult1a1, Cd48, Sell, Nkg7, Btg2
 	   Atp8b4, Clec1a, Ccnb2, Mt1, Sat1, D2Ertd750e, Bex6, Ctsg, Ighv1-30, Plk2
 PC_ 5
 Positive:  Car2, Ctla2b, uc008xuk.1, Muc13, Acap1, Itga2b, Mfsd2b, Gata1, uc007mja.2,uc007mjb.2, Esam
 	   Rnase6, Pf4, Mc5r, Itm2a, Tuba8, Cpa3, Gsta4, Sdsl, Cd63, Fam109b
 	   uc007zwh.1, Arntl, Rab38, Cdc20, Zbtb48, Ppic, Rgs1, Stard4, Armcx1, Phlda2
 Negative:  Ighv1-30, Gbp4, Gbp6, Sult1a1, Gpx3, Dntt, Clec1a, Iigp1, Wfdc17, Clec12a
 	   Plk2, Emb, Stat1, Klrb1b, Tgtp2, Nupr1, Atp8b4, Gm4951, Tgtp1, Cd74
 	   Rorc, Clu, uc008muc.2, Cnn3, Gbp3, Flt3, Tpm4, Ehd3, Sdpr, Smtnl1
</div>
```r
hsc.combined = RunTSNE(hsc.combined, dims.use = 1:30, max_iter=2000)

#get cell type
celltypes = rownames(hsc.combined@meta.data)
celltypes[grep('MPP',celltypes)]='MPP'; celltypes[grep('ST_HSC',celltypes)]='ST_HSC'; celltypes[grep('LT_HSC',celltypes)]='LT_HSC'; celltypes=factor(celltypes);

 Plot tsne results
plot.me <- data.frame(x=hsc.combined@reductions$tsne@cell.embeddings[,1],
                      y=hsc.combined@reductions$tsne@cell.embeddings[,2],
                      labels=Idents(hsc.combined),
                      celltypes=celltypes,
                      stringsAsFactors=FALSE)


unaligned.plot <- ggplot(plot.me, aes(x=x, y=y, colour = labels)) +
                  geom_point(size=2) +
                  scale_colour_manual(values=c("blue", "red")) +
                  xlab('t-SNE 1') +
                  ylab('t-SNE 2') +
                  theme_bw() +
                  theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA),
                        axis.line = element_line(colour = 'black',size=1))
plot(unaligned.plot)
```

![](scRNA_Workshop-alignment_files/figure-html/dataprecheck-1.png)<!-- -->

```r
unaligned.plot <- ggplot(plot.me, aes(x=x, y=y, colour = celltypes)) +
                  geom_point(size=2) +
                  scale_colour_manual(values=c("blue", "red","green")) +
                  xlab('t-SNE 1') +
                  ylab('t-SNE 2') +
                  theme_bw() +
                  theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA),
                        axis.line = element_line(colour = 'black',size=1))
plot(unaligned.plot)
```

![](scRNA_Workshop-alignment_files/figure-html/dataprecheck-2.png)<!-- -->

Yes, there are definitely condition-specific clusters of cells in the tSNE, so we will go ahead with alignment.

 scAlign setup
The general design of `scAlign` makes it agnostic to the input RNA-seq data representation. Thus, the input data can either be
gene-level counts, transformations of those gene level counts or a preliminary step of dimensionality reduction such
as canonical correlates or principal component scores. Here we create the scAlign object from the previously defined
`Seurat` objects and perform both PCA and CCA on the unaligned data.


```r
 Create paired dataset SCE objects to pass into scAlignCreateObject
youngMouseSCE <- SingleCellExperiment(
    assays = list(counts = youngMouseSeuratObj@assays$RNA@counts[genes.use,],
                  logcounts  = youngMouseSeuratObj@assays$RNA@data[genes.use,],
                  scale.data = youngMouseSeuratObj@assays$RNA@scale.data[genes.use,])
)

oldMouseSCE <- SingleCellExperiment(
  assays = list(counts = oldMouseSeuratObj@assays$RNA@counts[genes.use,],
                logcounts  = oldMouseSeuratObj@assays$RNA@data[genes.use,],
                scale.data = oldMouseSeuratObj@assays$RNA@scale.data[genes.use,])
)
```
We now build the scAlign SCE object and compute PCs and/or CCs using Seurat for the assay defined by `data.use`. It is assumed that `data.use`, which is being used for the initial step of dimensionality reduction, is properly normalized and scaled.
Resulting combined matrices will always be ordered based on the sce.objects list order.



 Alignment of young and old HSCs
Now we align the young and old cpopulations for multiple input types which are specified by `encoder.data`. `scAlign` returns a
low-dimensional joint embedding space where the effect of age is removed allowing us to use the complete dataset for downstream analyses such as clustering or differential expression. For the gene level input we also run the decoder procedure which projects each cell into logcount space for both conditions to perform paired single cell differential expressional.


```r
 Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=500, steps.decoder=500, log.every=100, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                    encoder.data="scale.data",
                    decoder.data="logcounts",
                    supervised='none',
                    run.encoder=TRUE,
                    run.decoder=TRUE,
                    log.dir=file.path(results.dir, 'models','gene_input'),
                    device="CPU")
```

<div class='r_output'> [1] "============== Step 1/3: Encoder training ==============="
 [1] "Graph construction"
 [1] "Adding source walker loss"
 [1] "Adding target walker loss"
 [1] "Done random initialization"
 [1] "Step: 1    Loss: 16.1246"
 [1] "Step: 100    Loss: 11.7762"
 [1] "Step: 200    Loss: 10.7414"
 [1] "Step: 300    Loss: 10.3558"
 [1] "Step: 400    Loss: 10.2768"
 [1] "Step: 500    Loss: 10.0424"
 [1] "============== Alignment Complete =============="
 [1] "============== Step 2/3: YOUNG decoder training ==============="
 [1] "Graph construction"
 [1] "Done random initialization"
 [1] "Step: 1    Loss: 0.5452"
 [1] "Step: 100    Loss: 0.4272"
 [1] "Step: 200    Loss: 0.3511"
 [1] "Step: 300    Loss: 0.305"
 [1] "Step: 400    Loss: 0.2796"
 [1] "Step: 500    Loss: 0.2731"
 [1] "============== Step 3/3: OLD decoder training ==============="
 [1] "Graph construction"
 [1] "Done random initialization"
 [1] "Step: 1    Loss: 0.5484"
 [1] "Step: 100    Loss: 0.4364"
 [1] "Step: 200    Loss: 0.3506"
 [1] "Step: 300    Loss: 0.3035"
 [1] "Step: 400    Loss: 0.2752"
 [1] "Step: 500    Loss: 0.2679"
</div>
```r
 Additional run of scAlign with PCA, the early.stopping heuristic terminates the training procedure too early with PCs as input so it is disabled.
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=500, log.every=100, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE),
                    encoder.data="PCA",
                    supervised='none',
                    run.encoder=TRUE,
                    run.decoder=FALSE,
                    log.dir=file.path(results.dir, 'models','pca_input'),
                    device="CPU")
```

<div class='r_output'> [1] "============== Step 1/3: Encoder training ==============="
 [1] "Graph construction"
 [1] "Adding source walker loss"
 [1] "Adding target walker loss"
 [1] "Done random initialization"
 [1] "Step: 1    Loss: 51.2401"
 [1] "Step: 100    Loss: 44.7082"
 [1] "Step: 200    Loss: 37.4414"
 [1] "Step: 300    Loss: 37.5215"
 [1] "Step: 400    Loss: 37.8422"
 [1] "Step: 500    Loss: 32.8113"
 [1] "============== Alignment Complete =============="
</div>
```r
 Additional run of scAlign with CCA
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=500, log.every=100, norm=TRUE, batch.norm.layer=TRUE, early.stop=TRUE),
                    encoder.data="CCA",
                    supervised='none',
                    run.encoder=TRUE,
                    run.decoder=FALSE,
                    log.dir=file.path(results.dir, 'models','cca_input'),
                    device="CPU")
```

<div class='r_output'> [1] "============== Step 1/3: Encoder training ==============="
 [1] "Graph construction"
 [1] "Adding source walker loss"
 [1] "Adding target walker loss"
 [1] "Done random initialization"
 [1] "Step: 1    Loss: 47.7481"
 [1] "Step: 100    Loss: 40.5103"
 [1] "Step: 200    Loss: 36.0581"
 [1] "Step: 300    Loss: 34.1863"
 [1] "Step: 400    Loss: 32.3583"
 [1] "Step: 500    Loss: 33.182"
 [1] "============== Alignment Complete =============="
</div>
```r
 Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by input labels
set.seed(5678)

#try changing the line below to "ALIGNED-PCA" or "ALIGNED-CCA" to check how alignment using just PCs or CCAs worked
DATA_TO_PLOT = 'ALIGNED-GENE';

gene_plot_colorbytype = PlotTSNE(scAlignHSC, DATA_TO_PLOT, title="scAlign-Gene", perplexity=30)
legend = get_legend(PlotTSNE(scAlignHSC, DATA_TO_PLOT, title="scAlign-Gene", legend="right", max_iter=1))
combined_plot = grid.arrange(gene_plot_colorbytype, legend, nrow = 1, layout_matrix=cbind(1,1,1,1,1,2))
```

![](scRNA_Workshop-alignment_files/figure-html/runscAlign-1.png)<!-- -->

```r
set.seed(5678)

gene_plot_colorbyGroup = PlotTSNE(scAlignHSC, DATA_TO_PLOT, title="scAlign-Gene", cols=c("red","blue"), labels.use="group.by", perplexity=30)
legend = get_legend(PlotTSNE(scAlignHSC, DATA_TO_PLOT, title="scAlign-Gene", cols=c("red","blue"), labels.use="group.by", legend="right", max_iter=1))
combined_plot = grid.arrange(gene_plot_colorbyGroup, legend, nrow = 1, layout_matrix=cbind(1,1,1,1,1,2))
```

![](scRNA_Workshop-alignment_files/figure-html/runscAlign-2.png)<!-- -->


 Paired differential expression of young and old cells.
Since we have run the decoder for "ALIGNED-GENE" alignment we can also investigate the paired differential of cells with respect to the young and old conditions. The projections are saved as "YOUNG2OLD" and "OLD2YOUNG" in `scAlignHSC` reducededDims. "YOUNG2OLD" indicates the projection of young cells into the expression space of old cells. As a reminder the combined matrices are always in the order of the `sce.objects` list passed to `scAlignCreateObject`.


```r
library(ComplexHeatmap)
```

<div class='r_output'> Loading required package: grid
</div>
<div class='r_output'> ========================================
 ComplexHeatmap version 2.2.0
 Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
 Github page: https://github.com/jokergoo/ComplexHeatmap
 Documentation: http://jokergoo.github.io/ComplexHeatmap-reference

 If you use it in published research, please cite:
 Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
   genomic data. Bioinformatics 2016.
 ========================================
</div>
```r
library(circlize)
```

<div class='r_output'> ========================================
 circlize version 0.4.8
 CRAN page: https://cran.r-project.org/package=circlize
 Github page: https://github.com/jokergoo/circlize
 Documentation: http://jokergoo.github.io/circlize_book/book/

 If you use it in published research, please cite:
 Gu, Z. circlize implements and enhances circular visualization
   in R. Bioinformatics 2014.
 ========================================
</div>
```r
all_data_proj_2_old = reducedDim(scAlignHSC, "YOUNG2OLD")
all_data_proj_2_young = reducedDim(scAlignHSC, "OLD2YOUNG")
scDExpr = all_data_proj_2_old - all_data_proj_2_young
colnames(scDExpr) <- genes.use

 Identify most variant DE genes
gene_var <- apply(scDExpr,2,var)
high_var_genes <- sort(gene_var, decreasing=T, index.return=T)$ix[1:40]

 Subset to top variable genes
scDExpr = scDExpr[,high_var_genes]

 Cluster each cell type individually
hclust_LT <- as.dendrogram(hclust(dist(scDExpr[which(cell_type == "LT"),], method = "euclidean"), method="ward.D2"))
hclust_ST <- as.dendrogram(hclust(dist(scDExpr[which(cell_type == "ST"),], method = "euclidean"), method="ward.D2"))
hclust_MPP <- as.dendrogram(hclust(dist(scDExpr[which(cell_type == "MPP"),], method = "euclidean"), method="ward.D2"))

clustering <- merge(hclust_ST, hclust_LT)
clustering <- merge(clustering, hclust_MPP)

 Reorder based on cluster groups
reorder = c(which(cell_type == "LT"), which(cell_type == "ST"), which(cell_type == "MPP"))

 Colors for ComplexHeatmap annotation
cell_type_colors = c("red", "green", "blue")
names(cell_type_colors) = unique(cell_type)

ha = rowAnnotation(df = data.frame(type=cell_type[reorder]),
    						   col = list(type = cell_type_colors),
                   na_col = "white")

#png(file.path(results.dir,"scDExpr_heatmap.png"), width=28, height=12, units='in', res=150, bg="transparent")
Heatmap(as.matrix(scDExpr[reorder,]),
        col = colorRamp2(c(-0.75, 0, 0.75), c('blue', 'white', 'red')),
        row_dend_width = unit(1, "cm"),
        column_dend_height = unit(1, "cm"),
  	    cluster_columns = TRUE,
        cluster_rows = clustering,
        row_dend_reorder=FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        na_col = 'white',
        column_names_gp = gpar(fontsize = 10),
        split=3,
        gap = unit(0.3, "cm")) + ha
```

![](scRNA_Workshop-alignment_files/figure-html/plotPairedDE-1.png)<!-- -->

```r
#dev.off()
```

 Discussion points

When does alignment work well?
Statistical significance?
When does alignment make sense?
How do you know when alignment makes sense?

 More reading

For a larger list of alignment methods, as well as an evaluation of them, see our scAlign paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1766-4. As a side note, some of the benefits of scAlign compared to most other methods include:

* ability to 'project' data from aligned space back to UMI/count space (for post-analysis).
* ability to use partial labeling of cells to improve alignment - e.g. suppose only some cells can be labeled using high confidence markers.
