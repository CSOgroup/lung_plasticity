# Plasticity of cancer and normal alveolar cells in lung adenocarcinoma

Repository containing the code for the analyses shown in the manuscript.

The clinical and molecular data for the various datasets was retreived from the sources cited in the manuscript and preprocessed with standard pipelines as outlined in the Methods section.

All the steps are performed following __lung_plasticity_pipeline.R__, which calls functions from __lung_plasticity_functions.R__ and __lung_plasticity_utils.R__

The functions to run Non-negative Matrix Factorization with patient cross-validation (cvNMF) to extract patient-independent transcriptional programs from malignant cells are stored in __cvNMF.R__ . The input is a standard count matrix (genes x cells) and an annotation file with cell IDs as row names and a 'Sample' column with patient names. The code is in R (v4) and requires the installation of R packages __Seurat__ (v4), __reshape2__, __RcppML__, __Matrix__, __corrplot__, __igraph__, __clustree__, __tidygraph__, __msigdbr__ and __fgsea__.

The demo below runs cvNMF on a subset of the dataset compendium used in the paper, which is available under __demo/data/__. 

```R
library(Seurat)
library(reshape2)
library(RcppML)
library(Matrix)
library(corrplot)
library(igraph)
library(clustree)
library(tidygraph)
library(msigdbr)
library(fgsea)

source("cvNMF.R")
DataDir = "demo/"

dataset = "qkwlxb"
load( file = paste0(DataDir,dataset,"_annot_subset.RData") ) # a random subset of 20 patients drawn from the harmonized and preprocessed metadata of qian, kim, wu, laughney, xing, bischoff (aka 'qkwlxb') datasets
load( file = paste0(DataDir,dataset,"_counts_subset.RData") ) # the corresponding counts matrix

# Setting parameters for demo
OutDir = paste0("demo/cvNMF/test_01/")
dir.create(OutDir,recursive=T,showWarnings = F)
kvec = c(5:30) # in the paper: c(5:100)
ntrials = 10 # same as in the paper
unbalanced_size_margin = 0.3 # in the paper: 0.1
cor_thresh = 0.75 # same as in the paper
recurrence_ratio = 0.2 # in the paper: 0.5

# Preprocessing
cvNMF_preprocessing( annot, counts, OutDir )

# cvNMF runs
dcat( "Running cvNMF" )
load(paste0(OutDir,"annot.RData"))
load(file=paste0(OutDir,"preprocessed_expression_mat.RData"))
for (chosen_rank in kvec)
{
	dcat( paste0("Current rank = ",chosen_rank),1 )
	this_OutDir = paste0(OutDir,"cvNMF_runs/","rank",chosen_rank,"/")
	dir.create(this_OutDir,showWarnings = F, recursive = T)

	# Computing 'ntrials' NMF models of ~evenly-partitioned (in terms of cell numbers) samples, and comparing latent factors
	cvNMF_partitioned_runs( annot, snp0, chosen_rank, ntrials, unbalanced_size_margin, cor_thresh, this_OutDir )

	# Computing 'ntrials' NMF models with all patients, comparing with results of partitioned runs, and extracting robust latent factors
	cvNMF_allpatients_run(snp0, chosen_rank, ntrials, cor_thresh, recurrence_ratio, this_OutDir)
}

# Aggregating results into co-occurrence matrix
ambient_genes = read.table(paste0(DataDir,"ambient_genes.txt"),header=T,stringsAsFactors=F,sep='\t')
exclude_genes = ambient_genes$ambient_genes
cvNMF_adjmat_construction(OutDir, kvec, exclude_genes)

# Clustering the co-occurrence matrix
res_range = seq(0.1,2,0.1)
cvNMF_adjmat_clustering(OutDir, res_range, minimum_n_lfs = 5)
clustering_solution = get_most_stable_clustering_solution(OutDir)

# Loading the list of transcriptional programs
load(file = paste0(OutDir,"louvain_clustering/",clustering_solution,"/tps.RData"))

# Gene ontology analysis
gsea_enrichment(OutDir, tps, exclude_genes)
```

The demo should yield 7 transcriptional programs with functional enrichments similar to a subset of those reported in the manuscript.

```
> sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /home/daniele/miniconda3/envs/r4/lib/libopenblasp-r0.3.17.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] fgsea_1.21.2       msigdbr_7.5.1      tidygraph_1.2.2    clustree_0.5.0    
 [5] ggraph_2.1.0       ggplot2_3.5.2      igraph_2.0.2       corrplot_0.92     
 [9] Matrix_1.5-3       RcppML_0.5.4       reshape2_1.4.4     SeuratObject_4.1.3
[13] Seurat_4.3.0      

loaded via a namespace (and not attached):
  [1] Rtsne_0.16             colorspace_2.0-3       RcppEigen_0.3.4.0.0   
  [4] deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4        
  [7] spatstat.data_3.0-0    leiden_0.4.3           listenv_0.8.0         
 [10] farver_2.1.1           graphlayouts_0.8.4     ggrepel_0.9.2         
 [13] fansi_1.0.3            codetools_0.2-18       splines_4.1.1         
 [16] knitr_1.41             polyclip_1.10-4        jsonlite_1.8.4        
 [19] ica_1.0-3              cluster_2.1.4          png_0.1-8             
 [22] uwot_0.1.14            ggforce_0.4.1          shiny_1.7.3           
 [25] sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.1.1        
 [28] httr_1.4.4             fastmap_1.1.0          lazyeval_0.2.2        
 [31] cli_3.6.4              later_1.3.0            tweenr_2.0.2          
 [34] htmltools_0.5.4        tools_4.1.1            gtable_0.3.1          
 [37] glue_1.6.2             RANN_2.6.1             dplyr_1.1.4           
 [40] fastmatch_1.1-3        Rcpp_1.0.8.3           scattermore_0.8       
 [43] vctrs_0.6.5            babelgene_22.9         spatstat.explore_3.0-5
 [46] nlme_3.1-160           progressr_0.12.0       lmtest_0.9-40         
 [49] spatstat.random_3.0-1  xfun_0.35              stringr_1.5.0         
 [52] globals_0.16.2         mime_0.12              miniUI_0.1.1.1        
 [55] lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3         
 [58] future_1.29.0          MASS_7.3-58.3          zoo_1.8-11            
 [61] scales_1.3.0           promises_1.2.0.1       spatstat.utils_3.1-3  
 [64] parallel_4.1.1         RColorBrewer_1.1-3     reticulate_1.26       
 [67] pbapply_1.6-0          gridExtra_2.3          stringi_1.7.8         
 [70] BiocParallel_1.28.3    rlang_1.1.6            pkgconfig_2.0.3       
 [73] matrixStats_0.63.0     lattice_0.20-45        ROCR_1.0-11           
 [76] purrr_1.0.2            tensor_1.5             patchwork_1.1.2       
 [79] htmlwidgets_1.6.2      cowplot_1.1.1          tidyselect_1.2.0      
 [82] parallelly_1.32.1      RcppAnnoy_0.0.20       plyr_1.8.8            
 [85] magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
 [88] DBI_1.1.3              pillar_1.9.0           withr_2.5.0           
 [91] fitdistrplus_1.1-8     survival_3.4-0         abind_1.4-5           
 [94] sp_1.5-1               tibble_3.2.1           future.apply_1.10.0   
 [97] KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_3.0-3   
[100] plotly_4.10.1          viridis_0.6.2          grid_4.1.1            
[103] data.table_1.14.6      digest_0.6.31          xtable_1.8-4          
[106] tidyr_1.2.1            httpuv_1.6.6           munsell_0.5.0         
[109] viridisLite_0.4.2     
> 
```
