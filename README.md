# SupeRJump

## Overview

`SupeRJump` is an R package for cell-fate modeling of single cell RNA analysis through a supervised jump-diffusion strategy.

It is designed to help users:

- Perform robust fate predictions for lineages and perform batch correction on fates across samples.
- Help identify preferentially biased cells toward particular lineages.
- With hypothesis generation for underlying mechanisms causing cell type skewing.

The package is especially useful for single cell RNA lineage tracking with a strong apriori of the system biology allowing for supervised hypothesis testing.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start and Example Workflow](#quick-start-and-example-workflow)
- [Main Features](#main-features)
- [Vignettes](#vignettes)
- [License](#license)
- [Contact](#contact)

---

## Installation

You can install the development version of `SupeRJump` from GitHub:

```r
# install.packages("remotes")
remotes::install_github("namwob44/SupeRJump")
```

Then load the package:

```r
library(SupeRJump)
```

You can access and download the opensource seurat object we use in the vignette at [this googledrive link](https://drive.google.com/drive/folders/1xb-wKcAoJxsIBXGXyxoW0IRLb9YZbWkG?usp=sharing). The data we are working on today is from [one of our prior works](https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00397-0). This dataset explores mutation order for acute myeloid leukemia. We will already start with a Seurat object labeled with cell types.

---

## Quick Start and Example Workflow

Here is an example walkthrough of generating the core of SupeRJump:

```r
 library(SupeRJump)

# Example input
seurat_obj <- readRDS("./Data/seurat_object.rds")

# Getting Semi-supervised Signatures for Pseudotime
Idents(seurat_obj)<-"cluster_annotations"
DefaultAssay(seurat_obj)<-"RNA"
markers.differences<-FindAllMarkers(seurat_obj)

terminal_state_list<-c("pDC","DC_11bn","DC_11bp", "Neutrophil", "DC_Other", "Macrophage")
reduced_markers<-markers.differences%>%
  filter(cluster%in%terminal_state_list)

gene_list_from_findallmarkers<-list()
for(iter in 1:length(terminal_state_list)){
  gene_list_from_findallmarkers[[iter]]<-markers.differences%>%
  filter(cluster%in%terminal_state_list)%>%
    filter(p_val_adj<0.01)%>%
    filter(avg_log2FC>1)%>%
    filter(cluster==terminal_state_list[iter])%>%
    pull(gene)%>%
    unique
}

# Building Pseudotime Assay
seurat_obj<-GetPseudoOrdering(seurat_obj,gene_list_from_findallmarkers,name_of_score =terminal_state_list )
seurat_obj<-CombinePseudoOrdering(seurat_obj,terminal_state_list)
# Transform Data for Model Fitting
seurat_obj<-CosineTransformPCAData(seurat_obj)
Y<-GetSinkData(seurat_obj,state_grouping_column_name="cluster_annotations",pseudotime_column_name = "Combined_Ordering")

state_grouping_column_name="cluster_annotations"
pseudotime_column_name = "Combined_Ordering"
batch_correction_column_name="Sample_names"

# Model Fitting
Prob_reduced<-Reduced_Jump_Prob_fast(seurat_obj,
                                     Total_Eigen_Data = Y,
                                     start_point = 1,
                                     state_grouping_column_name = state_grouping_column_name,
                                     pseudotime_column_name = pseudotime_column_name,
                                     batch_correction_column_name = batch_correction_column_name)

# Getting Transition Probability Matrix (TPM)
seurat_obj<-GetCellToCellTPM_fast(seurat_obj,Y,Prob_reduced,state_grouping_column_name="cluster_annotations",pseudotime_column_name = "Combined_Ordering",batch_correction_column_name="Sample_names",n_cores = 4,cell_indices_to_use = NULL,order = 5,eps = 1e-300,normalize_rows = TRUE,verbose = TRUE)

# Getting Absorbing states
seurat_obj<-GetAutomaticAbsorbingCellAssignment(seurat_obj,markers_df = markers.differences,mode = "markers",
                                                 state_grouping_column_name = state_grouping_column_name,
                                                 terminal_state_list = terminal_state_list,
                                                 pseudotime_column_name = pseudotime_column_name,
                                                 batch_correction_column_name = batch_correction_column_name)


seurat_obj<-GetFateMatrixAndMetrics(seurat_obj,state_grouping_column_name="cluster_annotations",
                                     pseudotime_column_name="Combined_Ordering",
                                     batch_correction_column_name="Sample_names",
                                     absorbing_state_column_name="absorbing_states")


seurat_obj<-GetMembership(seurat_obj,state_grouping_column_name="cluster_annotations")

seurat_obj<-GetAllClassifiedPoisedCells(seurat_obj,state_grouping_column_name="cluster_annotations",lineages_to_compare=c("Neutrophil","pDC","DC_11bn","DC_11bp","Macrophage"))

```

---

## Main Features

`SupeRJump` provides tools for:

1. **Supervised pseudotime signatures**  
   We allow for different kernels to be used to represent the cell state, specifically to order them along a pseudotime. We allow supervised signatures through either custom gene sets curator by users, or through curated sets with msigdbr.

2. **Batch Correction of Fate Probabilities**  
   A drawback of fates is the sensitivity to determining the appropriate number of absorbing states and preserving imbalanced datasets. A solution for this is to perform batch correction on the fate calculations to account for effect size. We perform a modified mean matching inspired by ComBAT to achieve this.

3. **Visitation Probability**  
   Visitation probability allows us to determine which cells are traversed before being absorbed at all. This measure is helpful to identify cells likely navigation through intermediate states.

4. **Weighted Destination Time**
   Developed a heuristic called to inform whether cells of interest are uniquely adept at visiting certain populations over others. A unique aspect of weighted destination time is this measure allows for condition specific flux to clusters. Further we can use it to assess enrichment of conditions compared to a null distribution.

5. **Preferential lineage bias cell detection**
   First, we group cells by their cell type classification to determine preferential biases towards a cluster. We then use cells from this cluster to determine preferential biases for each lineage independently. To find preferentially biased cells we rely only on lineage values from aggregated fate. We build an empirical cumulative density function for each lineage of the users’ choice. We identify robust outlier scores and nominate candidate cells as preferentially biased toward lineages.
   
6. **Multi-variate linear modeling between fates and transcription factors**
   The output of the model produces scores and p-values that signal which transcription factors are active and inactive with significance for each lineage, which are then used for hypothesis generation.  

7. **Integration with existing R workflows**  
   This package relies on the `Seurat` framework to contain each component and additional assay.


---

The workflow is as follows:
1. The user starts with an annotated Seurat object.
2. The output is updated seurat objects.
3. The result can be summarized, plotted, exported, or reused using ggplot2, or scCustomize and Seurat.

---

## Vignettes

Please go to the vignette folder and download either the html or vignette.RMD for more thorough explaination of the flow and pipeline.

---

## License

This package is licensed under the [MIT License](LICENSE).

---

## Contact

Maintainer: Michael Bowman  
Email: michael.bowman@pennmedicine.upenn.edu
GitHub: [https://github.com/namwob44](https://github.com/namwob44)

For bugs, feature requests, or questions, please open an issue:

[https://github.com/namwob44/SupeRJump/issues](https://github.com/namwob44/SupeRJump/issues)
