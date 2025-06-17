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
