
source( "luad_plasticity_utils.R" )
source( "luad_plasticity_functions.R" )

### Define pipeline run name/ID and MainDir
MainDir = "processed/cvNMF_pipeline/run_05/"
dir.create(MainDir)

tps_discovery = FALSE
if (tps_discovery){

	### Preprocessing plots
	OutDir = paste0(MainDir,"cnv_score_plot/")
	dir.create(OutDir)
	cnv_score_plot( OutDir )

	### Step 1: cvNMF filtering and co-occurrence matrix construction
	OutDir = paste0(MainDir,"tps_discovery/")
	dir.create(OutDir)
	ambient_genes = extract_ambient_genes( min_cw = 0.2, min_cor = 0.3, include_ig = TRUE)
	cvNMF_AdjMat_Construction(OutDir, exclude_genes = ambient_genes)

	### Step 2: adjmat clustering
	adjmat_clustering(OutDir)

	### Step 3: tps splitting
	clustering_solution = get_most_stable_clustering_solution(OutDir)
	tps_splitting(OutDir, clustering_solution = clustering_solution)
	gsea_enrichment( OutDir, clustering_solution, ambient_genes )
	# c2 into pEMT-hypoxia, c3 into interferon (x2) + HLA/B2M, c7 into Stress_HSP,AP1,Secreted, c11: MHC-II separate, but also other alveolar TPs

	### Step 4: collating clusters into transcriptional programs
	load(file = paste0(OutDir,"louvain_clustering/",clustering_solution,"/wrComm_genes.RData"))
	tps = wrComm_genes
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c2/c2_c1.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c2_pEMT = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c2/c2_c2.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c2_hypoxia = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c3/c3_c2.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c3_MHCI = as.character(df$V1)
	tps$c3_interferon = tps$c3[!(tps$c3 %in% tps$c3_MHCI)]
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c7/c7_c1.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c7_stresshsp = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c7/c7_c2.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c7_stressap1 = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c7/c7_c3.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c7_stresssecreted = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c11/c11_c1.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c11_lepidic = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c11/c11_c2.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c11_AT2club = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c11/c11_c3.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c11_club = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c11/c11_c4.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c11_MHCII = as.character(df$V1)
	df = read.table(file=paste0(OutDir,"louvain_clustering/splitting_c11/c11_c5.txt"),sep="\t",quote='',stringsAsFactors = F, header = F)
	tps$c11_AT2sftp = as.character(df$V1)
	tps = tps[!(names(tps) %in% c( "c2","c3","c7","c11" ))] # removing the splitted ones
	tps = tps[sort(names(tps))]
	names(tps)
	names(tps) = c( "Cell_proliferation","Translation_initiation","AT2-Club-like","AT2-like","Club-like","lepidic-like","MHC-II","Unas_emp","Unas_ccnd1","Unas_manygenes","Unas_asna1","Unas_ahsa1","Unas_abcb8",
		"Unas_crebbp","Unas_dgcr2","Hypoxia","pEMT","Unas_bola1","EMT","Unas_camk","Unas_casc3","Unas_anap",
		"Unas_bri3","Unas_c2small","Interferon","MHC-I","Basal-like","OxPhos","RNA_processing","Stress_AP1","Stress_HSP","Stress_secreted","Metal","Unfolded_protein_response")
	save(tps,file = paste0(OutDir,"tps_unvalidated.RData"))

	### Step 5: validating tps
	tps_validation(OutDir,exclude_genes = ambient_genes)
	OutDir = paste0( MainDir,"tps_discovery/tps_enrichments/" )
	tps_enrichments_barplots( OutDir )
	tps_enrichments_heatmap( OutDir,order_tps ) # hardcoded selected representative enrichments
	OutDir = paste0( MainDir,"tps_discovery/tps_representation/" )
	dir.create(OutDir)
	tps_representation( OutDir )
}

cells_normalize = FALSE
if (cells_normalize){

	### Step 6: Normalization
	extended_atlas_normalization( MainDir,exclude_genes = ambient_genes )
	extended_normal_atlas_normalization( MainDir,exclude_genes = ambient_genes )
}

tps_cells_scoring = FALSE
if (tps_cells_scoring){
	
	### Step 7: Scoring of extended + normal epithelial
	OutDir = paste0(MainDir,"scoring/")
	dir.create(OutDir)
	tps_scoring_extended_normal( OutDir )
	load(file = paste0(MainDir,"scoring/tps_scores_extendedAtlasNormal.RData"))
	scores_original = scores
	prefix = "EANregrAT2"
	OutDir = paste0(MainDir,"scoring/")
	regress_scores_normal_epi( OutDir, scores_original, prefix=prefix, ct_estimate_coefficients=c("AT2") )
	
	### tps MECO
	OutDir = paste0(MainDir,"scoring/tps_MECO/")
	dir.create(OutDir)
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	tps_MECO( OutDir, scores, prefix = "EANregrAT2" )
}

tps_clinicalcharacteristics_mutations = FALSE
if (tps_clinicalcharacteristics_mutations){
	
	### tps vs clinical characteristics
	order_tps = c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia",
					"Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" )
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	OutDir = paste0(MainDir,"tps_across_ClinicalCharacteristics/")
	dir.create(OutDir)
	tps_across_ClinicalCharacteristics_singlecells(OutDir, scores, order_tps)
	tps_across_ClinicalCharacteristics_singlecells_patientlevel(OutDir, scores, order_tps)
	tps_across_ClinicalCharacteristics_bulk_datasets(OutDir, order_tps, scoring_method = 'singscore')
	tps_across_ClinicalCharacteristics_heatmap( OutDir,scores,order_tps )

	### tps across mutations 
	OutDir = paste0(MainDir,"tps_vs_mutations/")
	dir.create(OutDir)
	load(file = paste0("data/extendedAtlasNormal_clin.RData"))
	clin["bp023t","mut_TP53"] = "mut" # initially missclassified 
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scores = scores[scores$Epi_Cancer=="Cancer",]
	tps_vs_mutations_singlecells( OutDir, order_tps, clin, scores )
	load(file = paste0( OutDir,"InvasiveLuads_wt/wtdf_corrected.RData" ))
	tps_vs_mutations_bulk( OutDir, order_tps, pairs_to_validate = wtdf_corrected, scoring_method = 'singscore' )
}

cytotrace_analysis = FALSE
if (cytotrace_analysis){

	### tps vs cytotrace
	OutDir = paste0(MainDir,"tps_vs_cytotrace/")
	dir.create(OutDir)
	load( file = paste0(MainDir,"extendedAtlasNormal_counts.RData") )
	load( file = paste0(MainDir,"extendedAtlasNormal_annot.RData") )
	annot = annot[annot$Epi_Cancer %in% c( "Cancer","AT0","AT1","AT2" ),]
	datasets = list()
	for (dataset in unique(annot$Dataset)){ datasets[[ dataset ]] = counts[,annot[annot$Dataset==dataset,"CellID"]] }
	saves_dir = paste0(OutDir,"saves/")
	dir.create(saves_dir)
	source("iCytoTRACE_modified.R")
	iCytoTRACE_preprocess(datasets,saves_dir)
	system("python scanoramaCT_modified.py", intern = TRUE, wait = TRUE, input = saves_dir)
	iCytoTRACE_downstream(datasets,saves_dir)
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	tps_vs_cytotrace( OutDir, order_tps, scores )
}

cellstates_discovery = FALSE
if (cellstates_discovery){
	
	### Cell states
	load(paste0(MainDir,"tps_discovery/tps.RData"))
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	OutDir = paste0(MainDir,"CellStates/")
	dir.create(OutDir)

	# Silhouette analysis
	OutDir = paste0(MainDir,"CellStates/Silhouette_Analysis_CellStates/")
	dir.create(OutDir)
	Silhouette_Analysis_CellStates(OutDir, order_tps, range = c(2:10)) # n = 2 wins 

	OutDir = paste0(MainDir,"CellStates/")
	for (thisk1 in c(1:10)){
		dcat(thisk1)
		for (thisk2 in c(1:10)){
			pair = paste0(as.character(thisk1),as.character(thisk2))
			if (pair %in% c( "11" )){ next }
			dcat(thisk2)
			CellStates_inference_hierarchical(OutDir,tps, solution_step1 = 2, nsubstates=c(thisk1,thisk2))
		}
	}
	OutDir = paste0(MainDir,"CellStates/Silhouette_Analysis_CellStates/")
	Silhouette_Analysis_CellStates_Hierarchical(OutDir, order_tps, range1 = c(1:10), range2 = c(1:10)) # n = 1,2 wins

	# Plotting ComplexHeatmaps for cell state profiles
	OutDir = paste0(MainDir,"CellStates/formatted_heatmaps/")
	dir.create(OutDir)
	load(file=paste0(OutDir,"../cs2_TPprofile_dots_weighted.RData"))
	colorz_twoways = colorRampPalette(c("gold","white","forestgreen"))(100)
	pdf(paste0(OutDir,"cs2_TPprofile_dots_weighted.pdf"),4,1.5)
	col.lim = signif(c(-max(abs(c(max(weighted_profile),min(weighted_profile)))),0,+max(abs(c(max(weighted_profile),min(weighted_profile))))),2)
	plot=(ComplexHeatmap::Heatmap(reorderMatt(weighted_profile,columns_only=TRUE),rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(weighted_profile), "cm")/3.5,height=unit(nrow(weighted_profile), "cm")/3.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(col.lim,c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(at=col.lim,grid_height = unit(1, "mm"),title_position="topcenter",title="Mean score",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	print(draw(plot,heatmap_legend_side="bottom"))
	dev.off()
	load(file=paste0(OutDir,"../hierarchical_2_substates12_TPprofile_weighted.RData"))
	pdf(paste0(OutDir,"hierarchical_2_substates12_TPprofile_weighted.pdf"),4,2)
	col.lim = signif(c(-max(abs(c(max(weighted_profile),min(weighted_profile)))),0,+max(abs(c(max(weighted_profile),min(weighted_profile))))),2)
	plot=(ComplexHeatmap::Heatmap(reorderMatt(weighted_profile,columns_only=TRUE),rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(weighted_profile), "cm")/3.5,height=unit(nrow(weighted_profile), "cm")/3.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(col.lim,c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(at=col.lim,grid_height = unit(1, "mm"),title_position="topcenter",title="Mean score",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	print(draw(plot,heatmap_legend_side="bottom"))
	dev.off()
	save(weighted_profile,file=paste0(OutDir,"../hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_weighted.RData"))
	
	# Cell states vs cytotrace
	load(paste0( MainDir,"CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	OutDir = paste0(MainDir,"CellStates/cs_vs_cytotrace/")
	dir.create(OutDir)
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	cs_vs_cytotrace( OutDir, scores, prefix='hierarchical_2_substates12_', cs_map=cs_map )
	
	# Cell state representations
	OutDir = paste0(MainDir,"CellStates/PatientLevel_TriangularPlot/")
	dir.create(OutDir)
	CellStates_PatientLevel_TriangularPlot( OutDir, scores, prefix='hierarchical_2_substates12_', cs1_name="Alveolar", cs2_name="Proliferative", cs3_name="Hypoxic" )

	OutDir = paste0(MainDir,"CellStates/CellStates_2D_asymmetrical_TPprofile_weighted/")
	dir.create(OutDir)
	load(paste0( MainDir,"CellStates/cs2_TPprofile_dots_weighted.RData"))
	profile2 = weighted_profile
	load(paste0( MainDir,"CellStates/hierarchical_2_substates12_TPprofile_weighted.RData"))
	profile3 = weighted_profile
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	preprefix = ""
	CellStates_2D_asymmetrical(OutDir, scores, profile2, profile3, cs_map, order_tps)

	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	OutDir = paste0(MainDir,"CellStates/stackedbars_heatmaps/")
	dir.create(OutDir)
	CellStates_stackedbars_heatmaps( OutDir,scores,order_tps,cs_map )

	# Extracting cell state signatures
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	OutDir = paste0(MainDir,"CellStates/signatures/")
	dir.create( OutDir )
	CellStates_signatures( OutDir,scores )
}

normal_cells_analyses = FALSE
if (normal_cells_analyses){

	### Normal cells analyses
	OutDir = paste0(MainDir,"tps_atlases_NormalCells/")
	epithelial_cells_dotplot( OutDir )
	dir.create(OutDir)
	tps_atlases_NormalCells_preprocessing( OutDir )
	tps_atlases_NormalCells_scoring( OutDir )
	tps_atlases_NormalCells( OutDir,order_tps )
	NormalCells_smoking( OutDir, order_tps )
	OutDir = paste0(MainDir,"tps_atlases_NormalCells/DiffExpr/")
	dir.create(OutDir)
	differential_expression_NormalCells( OutDir )

	# Including preinvasive samples
	OutDir = paste0(MainDir,"tps_atlases_NormalCells_withPreinvasive/")
	dir.create(OutDir)
	tps_atlases_NormalCells_preprocessing( OutDir,withPreinvasive=T )
	tps_atlases_NormalCells_scoring( OutDir,withPreinvasive=T )
	tps_atlases_NormalCells_withPreinvasive( OutDir,order_tps )
	OutDir = paste0(MainDir,"tps_atlases_NormalCells_withPreinvasive/DiffExpr/")
	dir.create(OutDir)
	differential_expression_NormalCells( OutDir, withPreinvasive=TRUE )
}

bulk_deconvolutions = FALSE
if (bulk_deconvolutions){

	### Cell state enrichments on bulk
	OutDir = paste0(MainDir,"bulk_deconvolutions/")
	dir.create(OutDir)

	# extract tumor-specific marker genes from each 
	OutDir = paste0(MainDir,"bulk_deconvolutions/tumorspecific_enrichments/")
	dir.create( OutDir )
	extract_tumorspecific_genes( OutDir )

	# tumor-specific enrichments
	tumorspecific_enrichments( OutDir,universe="permissive",signature="restrictive",scoring_method='singscore' )
	tumorspecific_enrichments( OutDir,universe="permissive",signature="TwoStates",scoring_method='singscore' )
}

visium_downstream = FALSE
if (visium_downstream){
	
	### Downstream analyses on Visium samples
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	OutDir = paste0(MainDir,"visium_downstream/")
	dir.create(OutDir)

	# cell state and tps scoring
	visium_cellstates_scoring( OutDir,method="AmsCentered",canceronly_genes=FALSE )
	visium_tps_scoring( OutDir,method="AmsCentered",canceronly_genes=FALSE )

	# tumor-normal boundary analyses
	OutDir = paste0(MainDir,"visium_downstream/visium_boundary_analysis_tumor_spacing500/")
	dir.create(OutDir)
	visium_boundary_analysis_tumor( OutDir,spacing=500 )
	OutDir = paste0(MainDir,"visium_downstream/visium_boundary_analysis_tumor_spacing250/")
	dir.create(OutDir)
	visium_boundary_analysis_tumor( OutDir,spacing=250 )

	OutDir = paste0(MainDir,"visium_downstream/visium_boundary_analysis_normal_spacing500/")
	dir.create(OutDir)
	visium_boundary_analysis_normal( OutDir,spacing=500 )
	OutDir = paste0(MainDir,"visium_downstream/visium_boundary_analysis_normal_spacing250/")
	dir.create(OutDir)
	visium_boundary_analysis_normal( OutDir,spacing=250 )

	# tme analyses in Visium
	OutDir = paste0(MainDir,"visium_downstream/")
	visium_tme_correlation( OutDir )
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))	
	# visium_tme_PvsH( OutDir, subspa_all, cutoff = 25, use1hopNN = TRUE )
	# visium_tme_PvsH( OutDir, subspa_all, cutoff = 25, use1hopNN = FALSE )
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 10, use1hopNN = TRUE )
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 10, use1hopNN = FALSE )
}

run_xenium_analyses = FALSE
if (run_xenium_analyses){

	### Xenium analyses
	whichDataset = "xenium_18May2023"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/")
	dir.create(OutDir)
	xenium_analyses( OutDir, whichDataset )
	xenium_gridding( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/cancer_cells/")
	dir.create(OutDir)
	xenium_cancercells( OutDir, whichDataset )

	whichDataset = "xenium_05February2024"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/")
	dir.create(OutDir)
	xenium_analyses( OutDir, whichDataset )
	xenium_gridding( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	
	# OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/cancer_cells/")
	# dir.create(OutDir)
	# xenium_cancercells( OutDir, whichDataset ) # only 1 gene (!) for the hypoxic signature

	whichDataset = "xenium_nondiseasedlung"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_nondiseasedlung_vs_tumor/")
	dir.create(OutDir)
	xenium_analyses_nondiseasedlung_vs_tumor( OutDir, density_threshold = 0.002 )
}

tme_analyses = FALSE
if (tme_analyses){

	### Single cell TME analyses
	OutDir = paste0(MainDir,"tps_cs_vs_tme/")
	dir.create(OutDir)
	load(paste0( OutDir, "../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	prefix = "h12_"
	TME_vs_tps_cs(OutDir, scores, prefix, order_tps)
	OutDir = paste0(MainDir,"tps_cs_vs_tme/formatted_heatmaps/")
	dir.create(OutDir)
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	TME_formatted_heatmaps( OutDir,cs_map )

	# same but for two-cs
	OutDir = paste0(MainDir,"tps_cs_vs_tme/")
	load(paste0( OutDir, "../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	scores$cs = scores$cs_level1
	prefix = "twocs_"
	TME_vs_tps_cs(OutDir, scores, prefix, order_tps)
	OutDir = paste0(MainDir,"tps_cs_vs_tme/formatted_heatmaps/")
	dir.create(OutDir)
	cs_map = data.frame(row.names=c("cs1","cs2"),alias=c("Alveolar","Dedifferentiated" ),colorz = c("dodgerblue4","firebrick4"),stringsAsFactors=F)
	TME_formatted_heatmaps( OutDir,cs_map )

	# formatting UMAP for tme
	OutDir = paste0(MainDir,"tps_cs_vs_tme/formatted_umaps/")
	dir.create(OutDir)
	TME_formatted_umaps( OutDir )

	OutDir = paste0(MainDir,"tps_cs_vs_tme_NoMetastases/")
	dir.create(OutDir)
	load(paste0( OutDir, "../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	prefix = "h12_"
	TME_vs_tps_cs(OutDir, scores, prefix, order_tps, noMets = T)
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	OutDir = paste0(MainDir,"tps_cs_vs_tme_NoMetastases/formatted_heatmaps/")
	dir.create(OutDir)
	TME_formatted_heatmaps( OutDir,cs_map )
}

run_cin_analyses = FALSE
if (run_cin_analyses){

	### Tps vs chromosomal instability
	OutDir = paste0(MainDir,"cin_analyses/")
	dir.create(OutDir)
	cin_analyses( OutDir )
}

run_facs_data_analyses = FALSE 
if (run_facs_data_analyses){

	### FACS data analysis
	OutDir = paste0(MainDir, "facs_data_analyses/")
	dir.create(OutDir)
	facs_data_analyses( OutDir )
}

get_n_samples_patients_cells = FALSE
if (get_n_samples_patients_cells){

	### Get number of cells and samples/patients
	# including preinvasive
	CellTypeDir = paste0("data/cell_typing/non_malignant/")
	load( file = paste0(MainDir,"extendedAtlasS0normal_annot.RData") )
	N_cells_cancer = sum(annot$Epi_Cancer=="Cancer")
	annot_cancer = annot
	total = 0
	total_patients = 0
	total_samples = 0
	for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu","wangs0","zhus0")){
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_annot_annd.RData" ))
		if (dataset=="wang"){ 
			sampless = unique(annot$Sample)
			patients_wang = unique(annot$Patient)
		 }
		if (dataset=="wangs0"){ 
			annot = annot[!(annot$Sample %in% sampless), ] 
			# annot = annot[!(annot$SampleType %in% c( "AAH","AIS","MIA" )), ] 
			patients_wangs0 = unique(annot$Patient)
		}
		# annot = annot[annot$Patient %in% scores$Patient,]
		# dcat(length(unique(annot$Sample)))

		total_patients = total_patients + length(unique(annot$Patient))
		total_samples = total_samples + length(unique(annot$Sample))
		total = total + nrow(annot)
	}
	print(paste0( "N total cells = ",total+N_cells_cancer )) # 1124285
	print(paste0( "N cancer cells = ",N_cells_cancer ))	# 114984
	print(paste0( "N samples = ",total_samples )) # 243
	print(paste0( "N samples = ",total_patients-length(intersect(patients_wang,patients_wangs0)) )) # 170

	# not including preinvasive samples
	CellTypeDir = paste0("data/cell_typing/non_malignant/")
	load( paste0("processed/cvNMF_pipeline/run_05/scoring/EANregrAT2_tps_scores_DatasetRegressed.RData") )
	N_cells_cancer = sum(scores$Epi_Cancer=="Cancer")
	samples_cancer_primary = unique(scores$Sample[scores$SampleType=="Primary"])
	samples_cancer_met = unique(scores$Sample[scores$SampleType=="Metastasis"])
	samples_cancer_normal = unique(scores$Sample[scores$SampleType=="Normal"])
	total = 0
	total_patients = 0
	total_samples = 0
	samples_nonmalignant_primary = c()
	samples_nonmalignant_met = c()
	samples_nonmalignant_normal = c()
	for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu","wangs0","zhus0")){
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_annot_annd.RData" ))
		if (dataset=="wang"){ 
			sampless = unique(annot$Sample)
			patients_wang = unique(annot$Patient)
		 }
		if (dataset=="wangs0"){ 
			annot = annot[!(annot$Sample %in% sampless), ] 
			annot = annot[!(annot$SampleType %in% c( "AAH","AIS","MIA" )), ] 
			patients_wangs0 = unique(annot$Patient)
		}
		annot = annot[annot$Patient %in% scores$Patient,]
		dcat(length(unique(annot$Sample)))
		samples_nonmalignant_primary = c(samples_nonmalignant_primary,unique(annot$Sample[annot$SampleType=="Primary"]))
		samples_nonmalignant_met = c(samples_nonmalignant_met,unique(annot$Sample[annot$SampleType=="Metastasis"]))
		samples_nonmalignant_normal = c(samples_nonmalignant_normal,unique(annot$Sample[annot$SampleType=="Normal"]))
		total_patients = total_patients + length(unique(annot$Patient))
		total_samples = total_samples + length(unique(annot$Sample))
		total = total + nrow(annot)
	}
	print(paste0( "N total cells = ",total+N_cells_cancer )) # 925285
	print(paste0( "N cancer cells = ",N_cells_cancer ))	# 106149
	print(paste0( "N samples = ",total_samples )) # 210
	print(paste0( "N samples = ",total_patients )) # 153

	# mets: 22
	all(samples_cancer_met %in% samples_nonmalignant_met)
	length(samples_cancer_met)
	length(union(samples_cancer_met,samples_nonmalignant_met))

	# primary: 125
	all(samples_cancer_primary %in% samples_nonmalignant_primary)
	length(samples_cancer_primary)
	length(union(samples_cancer_primary,samples_nonmalignant_primary))

	# normal: 63
	all(samples_cancer_normal %in% samples_nonmalignant_normal)
	length(samples_cancer_normal)
	length(union(samples_cancer_normal,samples_nonmalignant_normal))
}

formatting_supplementary_tables = FALSE
if (formatting_supplementary_tables){

	### Formatting supplementary tables
	load(paste0(MainDir,"tps_discovery/tps.RData"))
	remap_tps = c("AT0-like","AT2-like","PSM - basal-like","Cell proliferation","Club-like","Migratory I","Hypoxia","Interferon","Lepidic-like","Metallothionine","MHC-I","MHC-II","OxPhos","Migratory II","RNA processing","Stress AP1","Stress HSP","Stress secreted","Translation initiation","Unassigned","UPR - ER stress")
	names(remap_tps) = names(tps)
	library(openxlsx)
	OutDir = paste0(MainDir,"supplementary_tables/")
	dir.create(OutDir)
	# tps
	load(paste0(MainDir,"tps_discovery/tps.RData"))
	names(tps) = remap_tps[names(tps)]
	df = dan.df(1:max(sapply(tps,length)),names(tps),data="")
	for (cn in colnames(df)){ df[1:length(tps[[cn]]),cn] = tps[[cn]] }
	write.xlsx(df,paste0(OutDir,"supplementary_table_tps.xlsx"))
	# tps enrichment atlases
	load(paste0(MainDir,"tps_discovery/tps.RData"))
	tp = "Hypoxia"
	load(paste0(MainDir,"tps_discovery/tps_enrichments/withAtlases/",tp,"_vs_PublishedMps_All.RData"))
	rownames(thisdf) = thisdf$mp
	df = dan.df(rownames(thisdf),names(tps),data=NA)
	for (cn in colnames(df)){
		load(paste0(MainDir,"tps_discovery/tps_enrichments/withAtlases/",cn,"_vs_PublishedMps_All.RData"))
		rownames(thisdf) = thisdf$mp
		df[rownames(thisdf),cn] = thisdf$qvals
	}
	colnames(df) = remap_tps[colnames(df)]
	df = df[!(substr(rownames(df),1,3) %in% c( "cao","cm_","coo","heF","mar" )),]
	rownames(df) = gsub( "kinker_","Kinker et al. ",rownames(df) )
	rownames(df) = gsub( "gavish_","Gavish et al. ",rownames(df) )
	rownames(df) = gsub( "barkley_","Barkley et al. ",rownames(df) )
	rownames(df) = gsub( "trav_","Travaglini et al. ",rownames(df) )
	rownames(df) = gsub( "hlca_","HLCA ",rownames(df) )
	rownames(df) = gsub( "han_","Han et al. ",rownames(df) )
	rownames(df) = gsub( "lepidic_classic","Tavernari et al. lepidic signature",rownames(df) )
	rownames(df) = gsub( "lepidic_augm","Tavernari et al. lepidic signature (augmented)",rownames(df) )
	rownames(df) = gsub( "solid_classic","Tavernari et al. solid signature",rownames(df) )
	rownames(df) = gsub( "solid_augm","Tavernari et al. solid signature (augmented)",rownames(df) )
	df[,"Gene sets"] = rownames(df)
	df = df[,c("Gene sets",colnames(df)[colnames(df)!="Gene sets"])]
	write.xlsx(df,paste0(OutDir,"supplementary_table_tpsEnrichmentAtlases.xlsx"))
	# tps enrichment msigdb
	load(paste0(MainDir,"tps_discovery/tps.RData"))
	df = dan.df(0,c( "pathway","pval","padj","overlap","size","Transcriptional program" ))
	for (tp in names(tps)){
		thisdf = dan.read(paste0(MainDir,"tps_discovery/tps_enrichments/",tp,"_MSigDB.txt"))
		if (nrow(thisdf)<2){ 
			dcat(tp) # in case I will add manually
			next 
		}
		thisdf[,"Transcriptional program"] = remap_tps[tp]
		df = rbind(df,thisdf)
	}
	thisdf = data.frame(pathway=c("GOBP_RESPIRATORY_GASEOUS_EXCHANGE_BY_RESPIRATORY_SYSTEM","GOBP_REGULATION_OF_PEPTIDASE_ACTIVITY"),padj=c(6.7143e-05,0.0070927),overlap=c(4,10),size=c(52,312),`Transcriptional program`=c("AT2-like","lepidic-like"),stringsAsFactors=F)
	colnames(thisdf) = gsub("\\."," ",colnames(thisdf))
	df = rbind(thisdf,df)
	df = df[,c("Transcriptional program",colnames(df)[colnames(df)!="Transcriptional program"])]
	write.xlsx(df,paste0(OutDir,"supplementary_table_tpsEnrichmentMsigdb.xlsx"))
	# visium characteristics: see table uploaded to Zenodo
	# dea normal cells + go
	df = dan.df(0,c("Cell type","Gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj"))
	df2 = dan.df(0,c( "Cell type","Genes upregulated in","pathway","padj","overlap","size" ))
	for (ct in c( "AT1","AT0","AT2" )){
		load(file=paste0(paste0(MainDir,"tps_atlases_NormalCells/DiffExpr/",ct,"_singlecellsAllGenes_pairing_patient/"),"mm.RData"))
		mm[,"Cell type"] = ct
		mm$Gene = rownames(mm)
		df = rbind(df,mm[,colnames(df)])
		enrUp = dan.read(paste0(MainDir,"tps_atlases_NormalCells/DiffExpr/",ct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
		enrUp = enrUp[substr(enrUp$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
		enrUp[,"Cell type"] = ct
		enrUp[,"Genes upregulated in"] = "Tumor"
		enrDown = dan.read(paste0(MainDir,"tps_atlases_NormalCells/DiffExpr/",ct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
		enrDown = enrDown[substr(enrDown$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
		enrDown[,"Cell type"] = ct
		enrDown[,"Genes upregulated in"] = "Normal"
		df2 = rbind(df2,enrUp[,colnames(df2)])
		df2 = rbind(df2,enrDown[,colnames(df2)])
	}
	write.xlsx(df,paste0(OutDir,"supplementary_table_deaNormalMast.xlsx"))
	write.xlsx(df2,paste0(OutDir,"supplementary_table_deaNormalGO.xlsx"))	
}










