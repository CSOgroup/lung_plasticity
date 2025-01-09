


source("utils/dan.functions.R")
source( "2.cvNMF_pipeline_functions.R" )

### Define pipeline run name/ID and MainDir
MainDir = "processed/cvNMF_pipeline/run_05/"
dir.create(MainDir)
order_tps = c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia",
				"Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" )
ambient_genes = extract_ambient_genes( min_cw = 0.2, min_cor = 0.3, include_ig = TRUE)

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

	# c2 into pEMT-hypoxia, c3 into interferon (x2) + HLA/B2M, c7 into Stress_HSP,AP1,Secreted, c11: MHC-II separate, but also other alveolar TPs, 

	# ### Step 3.1 (manual): gene ontology with MSigDB (H,GO:BP,C6,C8). I report splittings only if they show different enrichments
	# c1: Cell_proliferation, cell cycle G2M
	# c2: hallmark hypoxia, glycolisis, mTORC1
	## c2_c1: wound healing, oxidative stress, housekeeping, EMTIII Gavish
	## c2_c2: hallmark hypoxia, glycolisis, mTORC1
	# c3: interferon
	## c3_c1: interferon alpha, gamma, response to virus, defence response
	## c3_c2: MHC-I
	## c3_c3: interferon alpha, gamma, response to virus, defence response
	# c4: proliferating and proximal basal cell, MYC and E2F targets
	# c5: oxphos, respiration
	# c6: proximal ciliated cells, gobp RNA splicing/processing, some MYC targets
	# c7: stress
	## c7_c1: hallmark TNF signaling regulated by NFKB, stress response, heat response (HSP..)
	## c7_c2: hallmark TNF signaling regulated by NFKB, hallmark apoptosis
	## c7_c3: hallmark TNF signaling regulated by NFKB, response to cytokines/chemokines
	# c8: metal response
	# c9: unfolded protein response, protein maturation
	# c10: proliferating and proximal basal cell, translation initiation
	# c11: AT2, MHC-II
	## c11_c1: lepidic_augm, AT2
	## c11_c2: AT2 and club-like, lepidic_augm
	## c11_c3: secretory/club-like
	## c11_c4: MHC-II
	## c11_c5: purely AT2 (SFTPs)
	# c12: unclear/MALAT/EMP conserved
	# c13: 4 genes, nothing
	# c14: unclear, embryo development, metabolic process
	# c15: nothing
	# c16: nothing
	# c17: nothing
	# c18: nothing
	# c19: nothing
	# c20: nothing
	# c21: cell adhesion, hallmark EMT
	# c22: nothing
	# c23: nothing
	# c24: unclear/post translational modification
	# c25: nothing


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
	preprocessing_marjanovic_gene_sets( OutDir )
	tps_validation(OutDir,exclude_genes = ambient_genes)
	OutDir = paste0( MainDir,"tps_discovery/tps_enrichments/" )
	tps_enrichments_barplots( OutDir )
	tps_enrichments_heatmap( OutDir,order_tps ) # hardcoded selected representative enrichments
	OutDir = paste0( MainDir,"tps_discovery/tps_representation/" )
	dir.create(OutDir)
	tps_representation( OutDir )
	# comparing with Marjanovic et al. 

	# Check metastasis-specific programs
	OutDir = paste0( MainDir,"tps_discovery/metastasis_specific_programs/" )
	dir.create(OutDir)
	metastasis_specific_programs( OutDir )

}

cells_normalize = FALSE
if (cells_normalize){
	### Step 6: Normalization. Of extended atlas alone, then on extended + s0 + normal epithelial cells atlas, then of extended + normal epithelial (but no preinvasive)
	extended_atlas_normalization( MainDir,exclude_genes = ambient_genes )
	extended_normal_atlas_normalization( MainDir,exclude_genes = ambient_genes )

}

tps_cells_scoring = FALSE
if (tps_cells_scoring){
	### Step 7: Scoring of extended + normal epithelial. Then, regress out dataset and test if it worked (below, all tried approaches)
	OutDir = paste0(MainDir,"scoring/")
	dir.create(OutDir)
	# tps_scoring_extended_normal( OutDir )
	load(file = paste0(MainDir,"scoring/tps_scores_extendedAtlasNormal.RData"))
	scores_original = scores
	OutDir = paste0(MainDir,"scoring/test_dataset_patient_bias/")
	test_dataset_patient_bias( OutDir,scores_original,"EanOriginal" )
	prefix = "EANregrAllEpi"
	OutDir = paste0(MainDir,"scoring/")
	# regress_scores_normal_epi( OutDir, scores_original, prefix=prefix, ct_estimate_coefficients=c("AT1","AT2","AT0","preTB","club","ciliated","basal") )
	load(file = paste0(MainDir,"scoring/",prefix,"_tps_scores_DatasetRegressed.RData"))
	scores_regressed = scores
	OutDir = paste0(MainDir,"scoring/test_dataset_patient_bias/")
	test_dataset_patient_bias( OutDir,scores_regressed,paste0(prefix,"_Dataset") )
	compare_original_regressed( OutDir, prefix_original="EanOriginal", prefix_regressed="EANregrAllEpi_Dataset", prefix_output="Compare_EanOriginal_EANregrAllEpi" )
	prefix = "EANregrAT2"
	OutDir = paste0(MainDir,"scoring/")
	# regress_scores_normal_epi( OutDir, scores_original, prefix=prefix, ct_estimate_coefficients=c("AT2") )
	load(file = paste0(MainDir,"scoring/",prefix,"_tps_scores_DatasetRegressed.RData"))
	scores_regressed = scores
	OutDir = paste0(MainDir,"scoring/test_dataset_patient_bias/")
	test_dataset_patient_bias( OutDir,scores_regressed,paste0(prefix,"_Dataset") )
	compare_original_regressed( OutDir, prefix_original="EanOriginal", prefix_regressed="EANregrAT2_Dataset", prefix_output="Compare_EanOriginal_EANregrAT2" )
	### Step 8: tps MECO
	OutDir = paste0(MainDir,"scoring/tps_MECO/")
	dir.create(OutDir)
	load(file = paste0(MainDir,"scoring/tps_scores_extendedAtlasNormal.RData"))
	tps_MECO( OutDir, scores, prefix = "EanOriginal" )
	load(file = paste0(MainDir,"scoring/EANregrAllEpi_tps_scores_DatasetRegressed.RData"))
	tps_MECO( OutDir, scores, prefix = "EANregrAllEpi" )
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	tps_MECO( OutDir, scores, prefix = "EANregrAT2" )

	### scoring L2S
	OutDir = paste0(MainDir,"scoring/")
	load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelial.RData"))
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoS_SST/ClassicalSignature_29Oct_3PatientsSplitted/GeneSets/GeneSets_qval0.1_logfc1.RData")
	ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
	mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
	ribo.genes = ribo$Approved.symbol
	seu = subset(x = seu, features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))]) # 12677 x 101849
	db_rand = dan.barkley_MakeRand(seu,GeneSets, 3)
	scores = dan.Barkley_GeneToEnrichment_AmsCentered(seu, GeneSets, db_rand)
	regress_scores_normal_epi( OutDir, scores, prefix="L2S_classical", ct_estimate_coefficients=c("AT2"),tps=GeneSets )
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoS_SST/ClassicalSignature_29Oct_3PatientsSplitted/GeneSets/gsaug_10.RData")
	db_rand = dan.barkley_MakeRand(seu,gsaug, 3)
	scores = dan.Barkley_GeneToEnrichment_AmsCentered(seu, gsaug, db_rand)
	regress_scores_normal_epi( OutDir, scores, prefix="L2S_augmented", ct_estimate_coefficients=c("AT2"),tps=gsaug )

	method = "pearson"
	load(paste0(OutDir,"L2S_augmented_tps_scores_DatasetRegressed.RData"))
	l2s_scores = scores[scores$Epi_Cancer=="Cancer",]
	l2s_scores$subtract = l2s_scores$solid_up-l2s_scores$lepidic_up
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scores = scores[rownames(l2s_scores),]
	rdf = data.frame(row.names = order_tps,cs = order_tps, cor_with_lepidic = NA, cor_with_solid = NA, cor_with_subtract = NA)
	dir.create(paste0(OutDir,"tp_vs_l2s/"))
	ts = cbind( scores,l2s_scores[,c("lepidic_up","solid_up","subtract")] )
	for (n in order_tps )
	{
		dcat(n)
	   rdf[n,"cor_with_lepidic"] = cor(scores[,n],l2s_scores[rownames(scores),"lepidic_up"],method=method)
	   rdf[n,"cor_with_solid"] = cor(scores[,n],l2s_scores[rownames(scores),"solid_up"],method=method)
	   rdf[n,"cor_with_subtract"] = cor(scores[,n],l2s_scores[rownames(scores),"subtract"],method=method)
	   ts$x = ts$subtract
	   ts$y = ts[,n]
	   this_color = "gray33"
	   if (n %in% c("AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II" )) { this_color = "dodgerblue4" }
	   if (n %in% c("Cell_proliferation","Hypoxia","Basal-like","EMT","pEMT" )) { this_color = "firebrick" }
	   if (!(n %in% c( "AT2-like","AT2-Club-like","Cell_proliferation","Hypoxia","Basal-like" ))) { next }
	   fileName = paste0(OutDir,"tp_vs_l2s/augmented_tp_l2s_density_",n,"_loess.pdf")
	   plotTitle = paste0("Pearson R = ",signif(cor(scores[,n],l2s_scores[rownames(scores),"subtract"],method=method),2))
	   # dan.scatterplot( fileName, ts$x, ts$y, xlab = "lepidic-to-solid scores", ylab = paste0( n," scores" ), plotTitle = plotTitle, dotSize = 0.1, plotFitLine = T, FitLineMethod = "loess", FitLineColor = this_color, plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 2, fileHeight = 2 )
	   pdf(paste0(OutDir,"tp_vs_l2s/augmented_tp_l2s_density_",n,"_loess.pdf"),1.4,1.4)
	   p = ggplot(ts,aes(x=x,y=y)) + geom_point(aes(x, y),size = 0.1,color="gray") + geom_smooth(color=this_color) + xlab( "lepidic-to-solid scores" ) + ylab( paste0( n," scores" ) ) + ggtitle(plotTitle) + theme_classic(base_size=6)
	   p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	   print(p)
	   dev.off()
	}
	colorz_solid = colorRampPalette(c("dodgerblue4","gray77","firebrick3"))(100)
	save(rdf, file=paste0( OutDir,"tp_vs_l2s/tps_cor_with_LepidicToSolid_augmentedSignature.RData" ))
	for (metric in c("cor_with_lepidic", "cor_with_solid", "cor_with_subtract")){
	   fileName = paste0(OutDir,"tp_vs_l2s/",metric,"_augmented.pdf")
	   rdf = rdf[order(rdf[,metric]),]
	   pdf(fileName, 2.7,2)
	   p = ggplot(data=rdf, aes(x=factor(cs,levels = as.character(cs)), y=rdf[,metric],fill=rdf[,metric])) + geom_bar(stat="identity") + scale_fill_gradientn( colours=colorz_solid ) + ylab( paste0(method," R") ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	   p = p + ylim(c( -max(abs(rdf[,metric])),+max(abs(rdf[,metric])) ))
	   p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) + theme(legend.position="none")
	   print(p)
	   dev.off()
	}
	load(paste0(OutDir,"L2S_classical_tps_scores_DatasetRegressed.RData"))
	l2s_scores = scores[scores$Epi_Cancer=="Cancer",]
	l2s_scores$subtract = l2s_scores$solid_up-l2s_scores$lepidic_up
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scores = scores[rownames(l2s_scores),]
	rdf = data.frame(row.names = order_tps,cs = order_tps, cor_with_lepidic = NA, cor_with_solid = NA, cor_with_subtract = NA)
	ts = cbind( scores,l2s_scores[,c("lepidic_up","solid_up","subtract")] )
	for (n in order_tps )
	{
	   rdf[n,"cor_with_lepidic"] = cor(scores[,n],l2s_scores[rownames(scores),"lepidic_up"],method=method)
	   rdf[n,"cor_with_solid"] = cor(scores[,n],l2s_scores[rownames(scores),"solid_up"],method=method)
	   rdf[n,"cor_with_subtract"] = cor(scores[,n],l2s_scores[rownames(scores),"subtract"],method=method)
	   ts$x = ts$subtract
	   ts$y = ts[,n]
	   this_color = "gray33"
	   if (n %in% c("AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II" )) { this_color = "dodgerblue4" }
	   if (n %in% c("Cell_proliferation","Hypoxia","Basal-like","EMT","pEMT" )) { this_color = "firebrick" }
	   if (!(n %in% c( "AT2-like","AT2-Club-like","Cell_proliferation","Hypoxia","Basal-like" ))) { next }
	   # fileName = paste0(OutDir,"tp_vs_l2s/classical_tp_l2s_density_",n,"_loess.pdf")
	   plotTitle = paste0("Pearson R = ",signif(cor(scores[,n],l2s_scores[rownames(scores),"subtract"],method=method),2))
	   # dan.scatterplot( fileName, ts$x, ts$y, xlab = "lepidic-to-solid scores", ylab = paste0( n," scores" ), plotTitle = plotTitle, dotSize = 0.1, plotFitLine = T, FitLineMethod = "loess", FitLineColor = this_color, plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 2, fileHeight = 2 )
	   pdf(paste0(OutDir,"tp_vs_l2s/classical_tp_l2s_density_",n,"_loess.pdf"),1.4,1.4)
	   p = ggplot(ts,aes(x=x,y=y)) + geom_point(aes(x, y),size = 0.1,color="gray") + geom_smooth(color=this_color) + xlab( "lepidic-to-solid scores" ) + ylab( paste0( n," scores" ) ) + ggtitle(plotTitle) + theme_classic(base_size=6)
	   p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	   print(p)
	   dev.off()
	}
	save(rdf, file=paste0( OutDir,"tp_vs_l2s/tps_cor_with_LepidicToSolid_classicalSignature.RData" ))
	for (metric in c("cor_with_lepidic", "cor_with_solid", "cor_with_subtract")){
	   fileName = paste0(OutDir,"tp_vs_l2s/",metric,"_classical.pdf")
	   rdf = rdf[order(rdf[,metric]),]
	   pdf(fileName, 2.7,2)
	   p = ggplot(data=rdf, aes(x=factor(cs,levels = as.character(cs)), y=rdf[,metric],fill=rdf[,metric])) + geom_bar(stat="identity") + scale_fill_gradientn( colours=colorz_solid ) + ylab( paste0(method," R") ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	   p = p + ylim(c( -max(abs(rdf[,metric])),+max(abs(rdf[,metric])) ))
	   p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) + theme(legend.position="none")
	   print(p)
	   dev.off()
	}

	### scoring Marjanovic
	OutDir = paste0(MainDir,"scoring/")
	load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelial.RData"))
	load( "data/marjnmf_mp.RData" )
 	load( "data/marjcluComplete_mp.RData" )
	GeneSets = c(marjnmf_mp,marjcluComplete_mp)   
	ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
	mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
	ribo.genes = ribo$Approved.symbol
	seu = subset(x = seu, features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))]) # 12677 x 101849
	db_rand = dan.barkley_MakeRand(seu,GeneSets, 3)
	scores = dan.Barkley_GeneToEnrichment_AmsCentered(seu, GeneSets, db_rand)
	save(scores, file=paste0( OutDir,"marj_scores.RData" ))
	regress_scores_normal_epi( OutDir, scores, prefix="marj", ct_estimate_coefficients=c("AT2"),tps=GeneSets )
	OutDir = paste0(MainDir,"scoring/marj_cor/")
	dir.create(OutDir)
	load(paste0(OutDir,"../marj_tps_scores_DatasetRegressed.RData"))
	marj_scores = scores[scores$Epi_Cancer=="Cancer",]
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scores = scores[rownames(marj_scores),]
	rdf = data.frame(row.names = order_tps,cs = order_tps )
	for (n in order_tps )
	{
		for (mcn in names(GeneSets)) { rdf[n,paste0("cor_with_",mcn)] = cor(scores[,n],marj_scores[rownames(scores),mcn]) }
	}
	save(rdf, file=paste0( OutDir,"tps_cor_with_Marj.RData" ))
	for (metric in paste0("cor_with_",names(GeneSets))){
	   fileName = paste0(OutDir,metric,".pdf")
	   rdf = rdf[order(rdf[,metric]),]
	   pdf(fileName, 6,4)
	   p = ggplot(data=rdf, aes(x=factor(cs,levels = as.character(cs)), y=rdf[,metric])) + geom_bar(stat="identity") + ylab( "Pearson r" ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	   p = p + ylim(c( -max(abs(rdf[,metric])),+max(abs(rdf[,metric])) ))
	   print(p)
	   dev.off()
	}

	### Test association with Rb/Mt genes
	OutDir = paste0(MainDir,"tps_discovery/vs_Rb_Mt/")
	dir.create(OutDir)
	tps_association_Rb_Mt( OutDir )

	### Comparing variance explained
	OutDir = paste0(MainDir,"scoring/variance_explained_comparison/")
	dir.create(OutDir)
	variance_explained_comparison( OutDir,order_tps, variable = "Patient" )

}

tps_clinicalcharacteristics_mutations_progeny = FALSE
if (tps_clinicalcharacteristics_mutations_progeny){
	### Step 9: progeny analysis 
	order_tps = c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia",
					"Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" )
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	progeny_analysis( paste0(MainDir,"tps_discovery/"), scores, paste0(MainDir,"extendedAtlasNormal_seu_CancerEpithelial.RData"), order_tps )
	OutDir = paste0(MainDir,"tps_across_ClinicalCharacteristics/")
	dir.create(OutDir)
	tps_across_ClinicalCharacteristics_singlecells(OutDir, scores, order_tps)
	tps_across_ClinicalCharacteristics_singlecells_patientlevel(OutDir, scores, order_tps)
	tps_across_ClinicalCharacteristics_bulk_datasets(OutDir, order_tps, scoring_method = 'singscore')
	tps_across_ClinicalCharacteristics_heatmap( OutDir,scores,order_tps )
	### Step 11: tps across mutations 
	OutDir = paste0(MainDir,"tps_vs_mutations/")
	dir.create(OutDir)
	load(file = paste0("data/extendedAtlasNormal_clin.RData"))
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scores = scores[scores$Epi_Cancer=="Cancer",]
	# scores[,order_tps] = apply(scores[,order_tps],2,function(x) x-mean(x)) 
	tps_vs_mutations_singlecells( OutDir, order_tps, clin, scores )
	OutDir = paste0(MainDir,"tps_vs_mutations_correction/")
	dir.create(OutDir)
	clin["bp023t","mut_TP53"] = "mut"
	tps_vs_mutations_singlecells( OutDir, order_tps, clin, scores )
	tps_otherCellTypes_expression( OutDir, order_tps )
	load(file = paste0( OutDir,"InvasiveLuads_wt/wtdf_corrected.RData" ))
	tps_vs_mutations_bulk( OutDir, order_tps, pairs_to_validate = wtdf_corrected, scoring_method = 'singscore' )
	KrasG12_in_singlecells( OutDir, order_tps, scores )
	### tps across age/sex, in tumor
	OutDir = paste0(MainDir,"tps_vs_agesex_tumor/")
	dir.create(OutDir)
	clin_age_sex_annotation( OutDir )
	tps_vs_agesex_tumor( OutDir )

}

cytotrace_analysis = FALSE
if (cytotrace_analysis){
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
	### Step 12: cell states
	load(paste0(MainDir,"tps_discovery/tps.RData"))
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	OutDir = paste0(MainDir,"CellStates/")
	dir.create(OutDir)
	dir.create(paste0(OutDir,"CellStates_optimalK/"))
	CellStates_optimalK(paste0(OutDir,"CellStates_optimalK/"), scores,tps, datasets = NULL)
	for (thisk in c(5:10)){
		dcat(thisk)
		CellStates_inference(OutDir, scores, tps, nstates = thisk)
	}
	## silhouette analysis
	OutDir = paste0(MainDir,"CellStates/Silhouette_Analysis_CellStates/")
	dir.create(OutDir)
	Silhouette_Analysis_CellStates(OutDir, order_tps, range = c(2:10)) # n = 2 wins both in stability and in silhouette analyses

	OutDir = paste0(MainDir,"CellStates/")
	dir.create(paste0(OutDir,"CellStates_hierarchical_optimalK/"))
	CellStates_inference_hierarchical_optimalK(paste0(OutDir,"CellStates_hierarchical_optimalK/"),tps, solution_step1 = 2)
	CellStates_inference_hierarchical(OutDir, solution_step1 = 2, nsubstates=c(1,2))
	CellStates_inference_hierarchical(OutDir, solution_step1 = 2, nsubstates=c(2,1))
	CellStates_inference_hierarchical(OutDir, solution_step1 = 2, nsubstates=c(2,2))
	for (thisk1 in c(3:10)){
		dcat(thisk1)
		for (thisk2 in c(1:10)){
			pair = paste0(as.character(thisk1),as.character(thisk2))
			if (pair %in% c( "11","12","21","22" )){ next }
			dcat(thisk2)
			CellStates_inference_hierarchical(OutDir,tps, solution_step1 = 2, nsubstates=c(thisk1,thisk2))
		}
	}
	OutDir = paste0(MainDir,"CellStates/Silhouette_Analysis_CellStates/")
	Silhouette_Analysis_CellStates_Hierarchical(OutDir, order_tps, range1 = c(1:10), range2 = c(1:10))

	# plotting ComplexHeatmaps
	OutDir = paste0(MainDir,"CellStates/formatted_heatmaps/")
	dir.create(OutDir)
	load(file=paste0(OutDir,"../cs2_TPprofile_dots_weighted.RData"))
	colorz_twoways = colorRampPalette(c("gold","white","forestgreen"))(100)
	pdf(paste0(OutDir,"cs2_TPprofile_dots_weighted.pdf"),4,1.5)
	col.lim = signif(c(-max(abs(c(max(weighted_profile),min(weighted_profile)))),0,+max(abs(c(max(weighted_profile),min(weighted_profile))))),2)
	plot=(ComplexHeatmap::Heatmap(reorderMatt(weighted_profile,columns_only=TRUE),rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(weighted_profile), "cm")/3.5,height=unit(nrow(weighted_profile), "cm")/3.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(col.lim,c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(at=col.lim,grid_height = unit(1, "mm"),title_position="topcenter",title="Mean score",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	# plot=(ComplexHeatmap::Heatmap(reorderMatt(weighted_profile,columns_only=TRUE),rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(weighted_profile), "cm")/3.5,height=unit(nrow(weighted_profile), "cm")/3.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(col.lim,c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height = unit(1, "mm"),title_position="topcenter",title="Mean score",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	print(draw(plot,heatmap_legend_side="bottom"))
	dev.off()
	
	load(file=paste0(OutDir,"../hierarchical_2_substates12_TPprofile_weighted.RData"))
	pdf(paste0(OutDir,"hierarchical_2_substates12_TPprofile_weighted.pdf"),4,2)
	col.lim = signif(c(-max(abs(c(max(weighted_profile),min(weighted_profile)))),0,+max(abs(c(max(weighted_profile),min(weighted_profile))))),2)
	plot=(ComplexHeatmap::Heatmap(reorderMatt(weighted_profile,columns_only=TRUE),rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(weighted_profile), "cm")/3.5,height=unit(nrow(weighted_profile), "cm")/3.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(col.lim,c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(at=col.lim,grid_height = unit(1, "mm"),title_position="topcenter",title="Mean score",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	# plot=(ComplexHeatmap::Heatmap(reorderMatt(weighted_profile,columns_only=TRUE),rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(weighted_profile), "cm")/3.5,height=unit(nrow(weighted_profile), "cm")/3.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(col.lim,c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height = unit(1, "mm"),title_position="topcenter",title="Mean score",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	print(draw(plot,heatmap_legend_side="bottom"))
	dev.off()

	save(weighted_profile,file=paste0(OutDir,"../hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_weighted.RData"))
	

	OutDir = paste0(MainDir,"CellStates/Cluster_Quality/")
	dir.create(OutDir)
	CellStates_Quality( OutDir )

	load(paste0( MainDir,"CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	OutDir = paste0(MainDir,"CellStates/cs_vs_cytotrace/")
	dir.create(OutDir)
	cs_vs_cytotrace( OutDir, scores, prefix='hierarchical_2_substates12_', cs_map=cs_map )
	# phate construction
	library(phateR)
	OutDir = paste0(MainDir,"CellStates/PHATE_plots/")
	dir.create(OutDir)
	data_phate80 = CellStates_Phate( OutDir, scores, prefix = 'hierarchical_2_substates12_', order_tps = order_tps,default = FALSE, knn=5, decay=40, gamma=1, t=80 )
	save(data_phate80,file = paste0( OutDir,"hierarchical_2_substates12_phate_data_phate80.RData" ))
	
	data_phateAuto = CellStates_Phate( OutDir, scores, prefix = 'hierarchical_2_substates12_', order_tps = order_tps,default = FALSE, knn=5, decay=40, gamma=1, t='auto' )
	save(data_phateAuto,file = paste0( OutDir,"hierarchical_2_substates12_phate_data_phateAuto.RData" ))

	load(file = paste0( OutDir,"hierarchical_2_substates12_phate_data_phate80.RData" ))
	CellStates_Phate_colorby( OutDir, scores, data_phate80, prefix = "t80_hierarchical_2_substates12_",by = c( "cs_level1","cs","Dataset","SampleType","Stage_collapsed","cyto" ), order_tps = order_tps, cs_map = cs_map )

	load(file = paste0( OutDir,"hierarchical_2_substates12_phate_data_phateAuto.RData" ))
	CellStates_Phate_colorby( OutDir, scores, data_phateAuto, prefix = "tAuto_hierarchical_2_substates12_",by = c( "cs_level1","cs","Dataset","SampleType","Stage_collapsed","cyto" ), order_tps = order_tps, cs_map = cs_map )

	OutDir = paste0(MainDir,"CellStates/PHATE_tests/")
	dir.create(OutDir)
	CellStates_Phate_test(OutDir,scores,order_tps,cs_map)

	for (t in c(10,20,30,40,50,60,70)){
		dcat(t)
		data_phate = CellStates_Phate( OutDir, scores, prefix = 'hierarchical_2_substates12_', order_tps = order_tps,default = FALSE, knn=5, decay=40, gamma=1, t=t )
		save(data_phate,file = paste0( OutDir,"hierarchical_2_substates12_phate_data_phate",t,".RData" ))
		CellStates_Phate_colorby( OutDir, scores, data_phate, prefix = paste0("t",t,"_hierarchical_2_substates12_"),by = c( "cs_level1","cs","Dataset","SampleType","Stage_collapsed","cyto" ), order_tps = order_tps, cs_map = cs_map )
	}

	for (t in c(10,20,30,40,50,60,70)){
		dcat(t)
		data_phate = CellStates_Phate( OutDir, scores, prefix = 'hierarchical_2_substates12_', order_tps = order_tps,default = FALSE, knn=5, decay=40, gamma=1, t=t )
		save(data_phate,file = paste0( OutDir,"hierarchical_2_substates12_phate_data_phate",t,".RData" ))
		CellStates_Phate_colorby( OutDir, scores, data_phate, prefix = paste0("t",t,"_hierarchical_2_substates12_"),by = c( "cs_level1","cs","Dataset","SampleType","Stage_collapsed","cyto" ), order_tps = order_tps, cs_map = cs_map )
	}
	
	OutDir = paste0(MainDir,"CellStates/PatientLevel_TriangularPlot/")
	dir.create(OutDir)
	CellStates_PatientLevel_TriangularPlot( OutDir, scores, prefix='hierarchical_2_substates12_', cs1_name="Alveolar", cs2_name="Proliferative", cs3_name="Hypoxic" )

	OutDir = paste0(MainDir,"CellStates/CellStates_2D_asymmetrical_TPprofile_mean/")
	dir.create(OutDir)
	load(paste0( MainDir,"CellStates/cs2_TPprofile_dots.RData"))
	profile2 = cellstate_profile
	load(paste0( MainDir,"CellStates/hierarchical_2_substates12_TPprofile_dots.RData"))
	profile3 = cellstate_profile
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	preprefix = ""
	CellStates_2D_asymmetrical(OutDir, scores, profile2, profile3, cs_map, order_tps)

	OutDir = paste0(MainDir,"CellStates/CellStates_2D_asymmetrical_TPprofile_weighted/")
	dir.create(OutDir)
	load(paste0( MainDir,"CellStates/cs2_TPprofile_dots_weighted.RData"))
	profile2 = weighted_profile
	load(paste0( MainDir,"CellStates/hierarchical_2_substates12_TPprofile_weighted.RData"))
	profile3 = weighted_profile
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	preprefix = ""
	CellStates_2D_asymmetrical(OutDir, scores, profile2, profile3, cs_map, order_tps)

	OutDir = paste0(MainDir,"CellStates/CellStates_2D_asymmetrical_TPprofile_extremes/")
	dir.create(OutDir)
	load(paste0( MainDir,"CellStates/cs2_TPprofile_dots_extremes.RData"))
	profile2 = extremes_profile
	load(paste0( MainDir,"CellStates/hierarchical_2_substates12_TPprofile_extremes.RData"))
	profile3 = extremes_profile
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	preprefix = ""
	CellStates_2D_asymmetrical(OutDir, scores, profile2, profile3, cs_map, order_tps)

	# patient-driven heterogeneity
	OutDir = paste0(MainDir,"CellStates/patient_heterogeneity/")
	dir.create(OutDir)
	dataset_wise_patient_clustering( OutDir )
	load(file = paste0( MainDir,"CellStates/PHATE_plots/","hierarchical_2_substates12_phate_data_phate80.RData" ))
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	patient_cellstates_clustering( OutDir, scores )

	# archetype analysis (on positron)
	OutDir = paste0(MainDir,"CellStates/archetype_analysis/")
	dir.create(OutDir)
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	archetype_analysis( OutDir,scores,order_tps )

	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	OutDir = paste0(MainDir,"CellStates/stackedbars_heatmaps/")
	dir.create(OutDir)
	CellStates_stackedbars_heatmaps( OutDir,scores,order_tps,cs_map )

	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	OutDir = paste0(MainDir,"CellStates/signatures/")
	dir.create( OutDir )
	CellStates_signatures( OutDir,scores )

	### Trying hierarchical clustering
	OutDir = paste0(MainDir,"CellStates_hc/")
	dir.create(OutDir)
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	CellStates_hierarchical_clustering(OutDir, scores, order_tps)
	load(file = paste0(MainDir,"CellStates_hc/HieCluWard_tps_CellStates_10.RData"))
	OutDir = paste0(MainDir,"CellStates_hc/signatures/")
	dir.create( OutDir )
	scores$cs = scores$csn4
	scores$cs_level1 = scores$csn2
	cs_map = data.frame(row.names=c("cs1","cs2","cs3","cs4"), alias=c("Alveolar","Immunogenic","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","darkorange1","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4","firebrick4"),stringsAsFactors=F)
	CellStates_hc_signatures( OutDir,scores,cs_map )
	# vs Marjanovic
	OutDir = paste0(MainDir,"CellStates_hc/vs_Marjanovic/")
	dir.create(OutDir)
	CellStates_hc_vs_Marjanovic( OutDir )
}

monocle_analysis = FALSE
if (monocle_analysis){
	OutDir = paste0(MainDir,"monocle_analysis/")
	dir.create(OutDir)
}

ordering_analysis = FALSE
if (ordering_analysis){
	OutDir = paste0(MainDir,"ordering_analysis/")
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	dir.create( OutDir )
	load(paste0( MainDir,"CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	load(paste0( MainDir,"CellStates/hierarchical_2_substates12_TPprofile_dots.RData"))
	tps_ordering_analysis( OutDir, scores, cellstate_profile, topN = 'all', quantile_cutoff = 0.75 )
	tps_ordering_analysis( OutDir, scores, cellstate_profile, topN = 1, quantile_cutoff = 0.75 )
	tps_ordering_analysis( OutDir, scores, cellstate_profile, topN = 2, quantile_cutoff = 0.75 ) 
	tps_ordering_analysis( OutDir, scores, cellstate_profile, topN = 3, quantile_cutoff = 0.75 ) 
	cs_ordering_analysis( OutDir, scores, cs_map, order_tps )
	pca_ordering( OutDir )
	OutDir = paste0(MainDir,"ordering_analysis_hc/")
	dir.create(OutDir)
	load(file = paste0(MainDir,"CellStates_hc/HieCluWard_tps_CellStates_10.RData"))
	scores$cs = scores$csn4
	scores$cs_level1 = scores$csn2
	cs_map = data.frame(row.names=c("cs1","cs2","cs3","cs4"), alias=c("Alveolar","Immunogenic","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","darkorange1","firebrick3","sienna4"),stringsAsFactors=F)
	cs_ordering_analysis( OutDir, scores, cs_map, order_tps )
	pca_ordering( OutDir,scores,cs_map )
}

scenic_analyses = FALSE
if (scenic_analyses){
	# Regulon discovery was performed on qkwlxb. Then, regulon scoring and downstream analyses on extendedAtlasNormal
	OutDir = paste0( MainDir,"scenic_analyses/" )
	dir.create(OutDir)
	scenic_regulon_scoring( OutDir )
	load(file = paste0(MainDir,"scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scenic_regulons_vs_tps( OutDir,scores,order_tps )
	load(paste0( MainDir,"CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	scenic_regulons_vs_cs( OutDir,scores,cs_map )
	load(file = paste0( MainDir,"CellStates/PHATE_plots/hierarchical_2_substates12_phate_data_phate80.RData" ))
	scenic_regulons_on_phate( OutDir, data_phate80 )
	scenic_regulons_heatmap( OutDir )
	# Dataset-wise
	scenic_preprocess_dataset_wise( OutDir )
	# [ arboreto run on electron ]
	scenic_regulon_build_regulons_dataset_wise( OutDir )
}

normal_cells_analyses = FALSE
if (normal_cells_analyses){
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

	OutDir = paste0(MainDir,"tps_atlases_NormalCells_withPreinvasive/")
	dir.create(OutDir)
	tps_atlases_NormalCells_preprocessing( OutDir,withPreinvasive=T )
	tps_atlases_NormalCells_scoring( OutDir,withPreinvasive=T )
	tps_atlases_NormalCells_withPreinvasive( OutDir,order_tps )
	OutDir = paste0(MainDir,"tps_atlases_NormalCells_withPreinvasive/DiffExpr/")
	dir.create(OutDir)
	differential_expression_NormalCells( OutDir, withPreinvasive=TRUE )

	OutDir = paste0(MainDir,"tps_atlases_NormalCells_withHlca/")
	dir.create(OutDir)
	tps_atlases_NormalCells_preprocessing_scoring_withHlca( OutDir )
	tps_atlases_NormalCells_withHlca( OutDir,order_tps )

	OutDir = paste0(MainDir,"weak_AT2_signature/")
	dir.create(OutDir)
	weak_AT2_signature( OutDir )
	weak_AT2_analyses( OutDir )

	OutDir = paste0(MainDir,"NormalCells_trajectory_analysis/")
	dir.create(OutDir)
	NormalCells_trajectory_analysis( OutDir )

	### tps across age/sex, in normal
	OutDir = paste0(MainDir,"tps_vs_agesex_normal/")
	dir.create(OutDir)
	tps_vs_agesex_normal( OutDir )
}

bulk_deconvolutions = FALSE
if (bulk_deconvolutions){
	OutDir = paste0(MainDir,"bulk_deconvolutions/")
	dir.create(OutDir)
	# bayesprism_PrepareReference( OutDir, nstates = 2)
	# bayesprism_PrepareReference( OutDir, nstates = 3)
	### TCGA
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_htseq_counts/TCGA_LUAD_ge_rawcounts_GeneSymbols.RData")
	cn = colnames(ge)[substr(as.character(colnames(ge)),14,15) %in% c("01")]
	ge = ge[,cn]
	colnames(ge) = substr(cn,1,12)
	# bayesprism_deconvolution( OutDir,ge=ge,prefix="TCGA",nstates=2 ) # (on positron)
	# bayesprism_deconvolution( OutDir,ge=ge,prefix="TCGA",nstates=3 ) # (on positron)

	# for (dataset in c( "TCGA","ChenEAS","TRACERx" )){
	# 	load( paste0(MainDir,"bulk_deconvolutions/TCGA_InstaPrism.res_2states.RData" ) )
	# 	thetas_cs = InstaPrism.res@Post.ini.cs@theta
	# 	save(thetas_cs,file=paste0(MainDir,"bulk_deconvolutions/",dataset,"_InstaPrism.res_2states_thetascs.RData"))
	# 	load( paste0(MainDir,"bulk_deconvolutions/TCGA_InstaPrism.res_3states.RData" ) )
	# 	thetas_cs = InstaPrism.res@Post.ini.cs@theta
	# 	save(thetas_cs,file=paste0(MainDir,"bulk_deconvolutions/",dataset,"_InstaPrism.res_3states_thetascs.RData"))
	# }
	load(file=paste0(MainDir,"bulk_deconvolutions/TCGA_InstaPrism.res_3states_thetascs.RData"))
	bayesprism_downstream(OutDir=paste0(MainDir,"bulk_deconvolutions/TCGA/"), th=thetas_cs)

	### ChenEAS
	load("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_expcounts_GeneSymbols.RData")
	# bayesprism_deconvolution( OutDir,ge=ge,prefix="ChenEAS",nstates=2 ) # (on positron)
	# bayesprism_deconvolution( OutDir,ge=ge,prefix="ChenEAS",nstates=3 ) # (on positron)
	### TRACERx421
	Clin = readRDS(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/figurecode/data/20221109_TRACERx421_all_patient_df.rds"))
	Clin$Patient = Clin$cruk_id
	rownames(Clin) = Clin$cruk_id
	Clin = Clin[grepl("LUAD",Clin$histology_multi_full_genomically.confirmed),]
	Clin2 = readRDS(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/figurecode/data/20221109_TRACERx421_all_tumour_df.rds"))
	Clin2 = as.data.frame(Clin2,stringsAsFactors=F)
	Clin2$Patient = Clin2$cruk_id
	Clin2$Sample = Clin2$tumour_id_muttable_cruk
	rownames(Clin2) = Clin2$Sample
	Clin2 = Clin2[Clin2$Patient %in% Clin$Patient,]
	Clin = Clin[Clin$Patient %in% Clin2$Patient,]
	rownames(Clin) = Clin$Patient
	library(fst)
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_counts_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	ge = ge[,substr(colnames(ge),1,8) %in% rownames(Clin)]
	# bayesprism_deconvolution( OutDir,ge=ge,prefix="TRACERx",nstates=2 ) # (on positron)
	# bayesprism_deconvolution( OutDir,ge=ge,prefix="TRACERx",nstates=3 ) # (on positron)

	### enrichment on bulk (maybe tumor-specific bulk?)

	# step 1: extract tumor-specific marker genes from each 
	OutDir = paste0(MainDir,"bulk_deconvolutions/tumorspecific_enrichments/")
	dir.create( OutDir )
	# extract_tumorspecific_genes( OutDir )
	tumorspecific_enrichments( OutDir,universe="permissive",signature="restrictive",scoring_method='singscore' )
	# tumorspecific_enrichments( OutDir,universe="permissive",signature="permissive",scoring_method='singscore' )
	# tumorspecific_enrichments( OutDir,universe="restrictive",signature="restrictive",scoring_method='singscore' )
	# tumorspecific_enrichments( OutDir,universe="restrictive",signature="permissive",scoring_method='singscore' )
	tumorspecific_enrichments( OutDir,universe="permissive",signature="TwoStates",scoring_method='singscore' )
	# tumorspecific_enrichments( OutDir,universe="restrictive",signature="TwoStates",scoring_method='singscore' )
	# chosen: permissive-restrictive (the most reasonable)

	

}

visium_downstream = FALSE
if (visium_downstream){
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	# OutDir = paste0(MainDir,"bulk_deconvolutions/")
	# for (Sample in batch1_table$Sample){
	# 	print(Sample)
	# 	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/",Sample,"_subspa.RData"))
	# 	cells = rownames(subspa@meta.data[subspa@meta.data$to_drop=="Valid",])
	# 	ge = as.matrix(subspa[['Spatial']]@counts)[,cells]
	# 	bayesprism_deconvolution( OutDir,ge=ge,prefix=paste0("Visium_",Sample),nstates=2 ) # (on positron)
	# 	bayesprism_deconvolution( OutDir,ge=ge,prefix=paste0("Visium_",Sample),nstates=3 ) # (on positron)	
	# }
	OutDir = paste0(MainDir,"visium_downstream/")
	dir.create(OutDir)
	visium_cellstates_scoring( OutDir,method="AMS",canceronly_genes=FALSE )
	visium_cellstates_scoring( OutDir,method="AmsCentered",canceronly_genes=FALSE )
	visium_cellstates_scoring( OutDir,method="AmsCentered",canceronly_genes=TRUE )
	visium_tps_scoring( OutDir,method="AMS",canceronly_genes=FALSE )
	visium_tps_scoring( OutDir,method="AmsCentered",canceronly_genes=FALSE )

	# selected_boxplots
	visium_selected_boxplots( OutDir )

	for (Sample in batch1_table$Sample){
		print(Sample)
		load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/",Sample,"_subspa.RData"))

		cells = rownames(subspa@meta.data[subspa@meta.data$to_drop=="Valid",])
		ge = as.matrix(subspa[['Spatial']]@counts)[,cells]
		bayesprism_deconvolution( OutDir,ge=ge,prefix=paste0("Visium_",Sample),nstates=2 ) # (on positron)
		bayesprism_deconvolution( OutDir,ge=ge,prefix=paste0("Visium_",Sample),nstates=3 ) # (on positron)	
	}

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

	### with hc
	OutDir = paste0(MainDir,"visium_downstream_hc/")
	dir.create(OutDir)
	visium_cellstates_scoring_hc( OutDir,method="AmsCentered",canceronly_genes=FALSE )

	OutDir = paste0(MainDir,"visium_downstream/")
	visium_tme_correlation( OutDir )
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))	
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 25, use1hopNN = TRUE )
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 10, use1hopNN = TRUE )
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 25, use1hopNN = FALSE )
	# visium_tme_PvsH( OutDir, subspa_all, cutoff = 20, use1hopNN = FALSE )
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 10, use1hopNN = FALSE )
	# visium_tme_PvsH( OutDir, subspa_all, cutoff = 5, use1hopNN = FALSE )
	
}

visium_dezuani = FALSE
if (visium_dezuani){
	OutDir = paste0(MainDir,"visium_dezuani/")
	dir.create(OutDir)
	visium_dezuani_preprocessing( OutDir )
	visium_dezuani_deconvolution( OutDir )
	visium_dezuani_tps_cellstates_scoring( OutDir )
	visium_dezuani_tme_correlation( OutDir )
	load(file = paste0(OutDir,"preprocessing/","subspa_all.RData"))
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 25, use1hopNN = FALSE )
	visium_tme_PvsH( OutDir, subspa_all, cutoff = 10, use1hopNN = FALSE )
	
}

run_xenium_analyses = FALSE
if (run_xenium_analyses){
	whichDataset = "xenium_18May2023"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	# dir.create(OutDir)
	# xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/")
	# dir.create(OutDir)
	xenium_analyses( OutDir, whichDataset )
	xenium_gridding( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="CXCL14", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="CD24", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="MMP7", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="LGALS3BP", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/cancer_cells/")
	dir.create(OutDir)
	xenium_cancercells( OutDir, whichDataset )
	xenium_gridding_cancercells( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/spatialneighbourhood_tme/")
	dir.create(OutDir)
	xenium_spatialneighbourhood_tme( OutDir )

	whichDataset = "xenium_05February2024"
	# OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	# dir.create(OutDir)
	# xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/")
	# dir.create(OutDir)
	xenium_analyses( OutDir, whichDataset )
	xenium_gridding( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="CXCR4", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 ) # almost never expressed
	xenium_gridding_singlegene( OutDir, gene="GDF15", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 ) # almost never expressed
	
	# OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/cancer_cells/")
	# dir.create(OutDir)
	# xenium_cancercells( OutDir, whichDataset ) # only 1 gene (!) for the hypoxic signature

	# haga2023, https://kero.hgc.jp/Early_cancer.html
	whichDataset = "xenium_haga2023_TSU21"
	# OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	# dir.create(OutDir)
	# xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/")
	# dir.create(OutDir)
	xenium_analyses( OutDir, whichDataset )
	xenium_gridding( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="LGALS3BP", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="ISG20", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="OAS2", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="HLA-DPB1", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="CD74", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="MMP7", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	xenium_gridding_singlegene( OutDir, gene="CXCL14", whichDataset = whichDataset, grid_spacing = 200, density_threshold = 0.002 )

	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/cancer_cells/")
	dir.create(OutDir)
	xenium_cancercells( OutDir, whichDataset )
	xenium_gridding_cancercells( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	# whichDataset = "xenium_haga2023_TSU20"
	# OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	# dir.create(OutDir)
	# xenium_preprocessing( OutDir, whichDataset )

	whichDataset = "xenium_nondiseasedlung"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	OutDir = paste0(MainDir,"xenium_analyses_nondiseasedlung_vs_tumor/")
	dir.create(OutDir)
	xenium_analyses_nondiseasedlung_vs_tumor( OutDir, density_threshold = 0.002 )

	# takano2024
	whichDataset = "xenium_takano2024_luad2"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	whichDataset = "xenium_takano2024_luad3"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/") # not possible to discriminate normal alveolar and tumor cells
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	whichDataset = "xenium_takano2024_luad14"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	whichDataset = "xenium_takano2024_luad16"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )
	whichDataset = "xenium_takano2024_luad17"
	OutDir = paste0(MainDir,"xenium_preprocessing_",whichDataset,"/")
	dir.create(OutDir)
	xenium_preprocessing( OutDir, whichDataset )

	whichDataset = "xenium_takano2024_luad2"
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/")
	dir.create(OutDir)
	xenium_analyses( OutDir, whichDataset )
	xenium_gridding( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 )
	OutDir = paste0(MainDir,"xenium_analyses_",whichDataset,"/cancer_cells/")
	dir.create(OutDir)
	xenium_cancercells( OutDir, whichDataset )
	
}

run_cosmx_analyses = FALSE
if (run_cosmx_analyses){
	OutDir = paste0(MainDir,"cosmx_preprocessing/")
	dir.create(OutDir)
	# cosmx_preprocessing( OutDir )
	OutDir = paste0(MainDir,"cosmx_analyses/")
	dir.create(OutDir)
	cosmx_analyses( OutDir )
	cosmx_gridding( OutDir, grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="LGALS3BP", grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="CD74", grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="STAT1", grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="IFI27", grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="CXCL14", grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="CD24", grid_spacing = 200, density_threshold = 0.002 )
	cosmx_gridding_singlegene( OutDir, gene="MMP7", grid_spacing = 200, density_threshold = 0.002 )
	
	OutDir = paste0(MainDir,"cosmx_analyses/cancer_cells/")
	dir.create(OutDir)
	cosmx_cancercells( OutDir )
	cosmx_gridding_cancercells( OutDir, grid_spacing = 200, density_threshold = 0.002 )
}

run_logan_analyses = FALSE
if (run_logan_analyses){
	OutDir = paste0(MainDir,"logan_analyses/")
	dir.create(OutDir)
	logan_analyses( OutDir )
}

tme_analyses = FALSE
if (tme_analyses){
	OutDir = paste0(MainDir,"tps_cs_vs_tme/")
	dir.create(OutDir)
	load(paste0( OutDir, "../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	prefix = "h12_"
	TME_vs_tps_cs(OutDir, scores, prefix, order_tps)
	OutDir = paste0(MainDir,"tps_cs_vs_tme/formatted_heatmaps/")
	dir.create(OutDir)
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	TME_formatted_heatmaps( OutDir,cs_map )

	OutDir = paste0(MainDir,"tps_cs_vs_tme_13/")
	dir.create(OutDir)
	load(paste0( OutDir, "../CellStates/tps_CellStates_hierarchical_2_substates13.RData"))
	prefix = "h13_"
	TME_vs_tps_cs(OutDir, scores, prefix, order_tps)
	OutDir = paste0(MainDir,"tps_cs_vs_tme/formatted_heatmaps/")
	dir.create(OutDir)

	OutDir = paste0(MainDir,"tps_cs_vs_tme_CellStates_hc/")
	dir.create(OutDir)
	load(paste0( OutDir, "../CellStates_hc/HieCluWard_tps_CellStates_10.RData"))
	scores$cs = scores$csn4
	prefix = "hc4_"
	TME_vs_tps_cs(OutDir, scores, prefix, order_tps)

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
	OutDir = paste0(MainDir,"cin_analyses/")
	dir.create(OutDir)
	cin_analyses( OutDir )
	OutDir = paste0(MainDir,"cgas_sting_analyses/")
	dir.create(OutDir)
	cgas_sting_analyses( OutDir )

}

run_mellon_analyses = FALSE
if (run_mellon_analyses){
	OutDir = paste0(MainDir,"mellon_analyses/")
	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	dir.create(OutDir)
	mellon_analyses( OutDir,cs_map )
}

run_paga_analysis = FALSE
if (run_paga_analysis){
	OutDir = paste0(MainDir,"paga_analysis/")
	dir.create(OutDir)
}

run_tps_in_celllines = FALSE
if (run_tps_in_celllines){
	OutDir = paste0(MainDir,"tps_in_celllines/")
	dir.create(OutDir)
	tps_in_celllines( OutDir )	
}

run_facs_data_analyses = FALSE 
if (run_facs_data_analyses){
	OutDir = paste0(MainDir, "facs_data_analyses/")
	dir.create(OutDir)
	facs_data_analyses( OutDir )
}

run_hypometh_analyses = FALSE
if (run_hypometh_analyses){
	OutDir = paste0(MainDir,"hypometh_analyses/")
	dir.create(OutDir)
	hypometh_analyses( OutDir )
}

get_n_samples_patients_cells = FALSE
if (get_n_samples_patients_cells){

	#### including preinvasive
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

	#### not including preinvasive samples
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












