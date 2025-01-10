library(igraph)
library(tidygraph)
library(clustree)
library(ClusterR)
library(corrplot)
library(reshape2)
library(dendextend)
library(Seurat)
library(GSVA)
library(RColorBrewer)
library(colorRamp2)
library(ComplexHeatmap)
colorz_solid = colorRampPalette(c("dodgerblue4","gray77","firebrick3"))(100)
colorz_single = colorRampPalette(c("white","forestgreen"))(100)
colorz_twoways = colorRampPalette(c("gold","white","forestgreen"))(100)
lepidic_color = "dodgerblue4"
acinar_color = "orange"
papillary_color = "lightseagreen"
solid_color = "red"

cnv_score_plot = function( OutDir ){
	CellTypeDir = "data/cell_typing/"
	taall = dan.df(0,c( "Dataset","TN","SampleType","cnv_score","Epi_Cancer" ))
	for (dataset in c( "qian","kim","wu","laughney","xing","bischoff") )
	{
		dcat(dataset)
		bb = readRDS(paste0(CellTypeDir,"run.final.infercnv_obj_",dataset))
		cnv = bb@expr.data
		load( paste0("data/",dataset,"_annot_harmonized.RData") )
		ta = annot[colnames(cnv),]
		ta$cnv_score = colMeans(abs(cnv-1))
		load( paste0("data/",dataset,"_annot_EpiCancer_harmonized.RData") )
		annot = annot[annot$Epi_Cancer=="Cancer",]
		ta$Epi_Cancer = "Non-malignant"
		ta[intersect(rownames(ta),rownames(annot)),"Epi_Cancer"] = "Cancer"
		ta = ta[,c( "Dataset","TN","SampleType","cnv_score","Epi_Cancer" )]
		ta$Dataset = tolower(ta$Dataset)
		rownames(ta) = paste0(ta$Dataset,rownames(ta))
		taall = rbind(taall,ta)
	}
	for (dataset in c( "yang","zhu","salcher","he","wang","hu" ) )
	{
		dcat(dataset)
		bb = readRDS(paste0(CellTypeDir,"infercnv_",dataset,"/run.final.infercnv_obj"))
		cnv = bb@expr.data
		load( paste0("data/",dataset,"_annot_harmonized.RData") )
		ta = annot[colnames(cnv),]
		ta$cnv_score = colMeans(abs(cnv-1))
		load( paste0("data/",dataset,"_annot_EpiCancer_harmonized.RData") )
		annot = annot[annot$Epi_Cancer=="Cancer",]
		ta$Epi_Cancer = "Non-malignant"
		if (dataset %in% c( "he","salcher" )){
			rownames(ta) = paste0(substr(dataset,1,1),rownames(ta))
		}
		ta[intersect(rownames(ta),rownames(annot)),"Epi_Cancer"] = "Cancer"
		ta = ta[,c( "Dataset","TN","SampleType","cnv_score","Epi_Cancer" )]
		ta$Dataset = tolower(ta$Dataset)
		rownames(ta) = paste0(ta$Dataset,rownames(ta))
		taall = rbind(taall,ta)
	}
	x = taall$Dataset
	y = taall$cnv_score
	fill = taall$Epi_Cancer
	pdf(paste0(OutDir,"cnv_score_plot.pdf"),4,2.5)
	ylimLeft = 0
	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	plot = dan.boxplots.multipages( x, y, fill, xlab = "Dataset", ylab = "CNV score", filllab = "", signifTest = NULL, xColors = "black", ylimLeft=ylimLeft,ylimRight=ylimRight,fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL,legend_position="top" )
	print(plot)
	dev.off()
	dan.save(taall,file=paste0(OutDir,"cnv_score_plot.pdf"))
	load(file=paste0(OutDir,"cnv_score_plot.RData"))
	taall = object
	x = taall$Dataset
	y = taall$cnv_score
	fill = taall$Epi_Cancer
	pdf(paste0(OutDir,"cnv_score_plot_violins.pdf"),4,2.5)
	ylimLeft = 0
	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	plot=dan.violinplots.multipages( x, y, fill,xlab = "Dataset", ylab = "CNV score", filllab = "Cell type", signifTest = NULL, xColors = "black", ylimLeft=ylimLeft,ylimRight=ylimRight,fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL,legend_position="top" )
	# plot = dan.boxplots.multipages( x, y, fill, xlab = "Dataset", ylab = "CNV score", filllab = "Cell type", signifTest = NULL, xColors = "black", ylimLeft=ylimLeft,ylimRight=ylimRight,fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
}

cin_analyses = function( OutDir ){
	CellTypeDir = "data/cell_typing/"
	taall = dan.df(0,c( "Dataset","TN","SampleType","cnv_score","Epi_Cancer" ))
	load(paste0( MainDir, "CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	scores$cnv_scores = NA
	dtable(scores$Dataset)
	for (dataset in c( "qian","kim","wu","laughney","xing","bischoff") )
	{
		dcat(dataset)
		bb = readRDS(paste0(CellTypeDir,"run.final.infercnv_obj_",dataset))
		cnv = bb@expr.data
		load( paste0("data/",dataset,"_annot_harmonized.RData") )
		ta = annot[colnames(cnv),]
		ta$cnv_score = colMeans(abs(cnv-1))
		load( paste0("data/",dataset,"_annot_EpiCancer_harmonized.RData") )
		annot = annot[annot$Epi_Cancer=="Cancer",]
		ta$Epi_Cancer = "Non-malignant"
		ta[intersect(rownames(ta),rownames(annot)),"Epi_Cancer"] = "Cancer"
		ta = ta[,c( "Dataset","TN","SampleType","cnv_score","Epi_Cancer" )]
		ta$Dataset = tolower(ta$Dataset)
		taall = rbind(taall,ta)
	}
	for (dataset in c( "yang","zhu","salcher","he","wang","hu" ) )
	{
		dcat(dataset)
		bb = readRDS(paste0(CellTypeDir,"infercnv_",dataset,"/run.final.infercnv_obj"))
		cnv = bb@expr.data
		load( paste0("data/",dataset,"_annot_harmonized.RData") )
		ta = annot[colnames(cnv),]
		ta$cnv_score = colMeans(abs(cnv-1))
		load( paste0("data/",dataset,"_annot_EpiCancer_harmonized.RData") )
		annot = annot[annot$Epi_Cancer=="Cancer",]
		ta$Epi_Cancer = "Non-malignant"
		if (dataset %in% c( "he","salcher" )){
			rownames(ta) = paste0(substr(dataset,1,1),rownames(ta))
		}
		ta[intersect(rownames(ta),rownames(annot)),"Epi_Cancer"] = "Cancer"
		ta = ta[,c( "Dataset","TN","SampleType","cnv_score","Epi_Cancer" )]
		ta$Dataset = tolower(ta$Dataset)
		if (dataset=="wang") { rownames(ta) = paste0(ta$Dataset,rownames(ta)) }
		taall = rbind(taall,ta)
	}
	taall = taall[taall$Epi_Cancer=="Cancer",]
	scores$cnv_score = taall[rownames(scores),"cnv_score"]
	save(scores, file = paste0( OutDir,"scores_with_cnvscores.RData" ))
	load( file = paste0( OutDir,"scores_with_cnvscores.RData" ))

	rdf = data.frame(row.names = order_tps,tps = order_tps, cnv_score = as.numeric(cor(scores$cnv_score,scores[,order_tps],method="spearman")) )
	colorz_solid = colorRampPalette(c("dodgerblue4","gray77","firebrick3"))(100)
	fileName = paste0(OutDir,"cor_tps_cnvScore.pdf")
	rdf = rdf[order(rdf[,"cnv_score"]),]
	pdf(fileName, 3,2.2)
	p = ggplot(data=rdf, aes(x=factor(rdf$tps,levels = as.character(rdf$tps)), y=rdf[,"cnv_score"],fill=rdf[,"cnv_score"])) + geom_bar(stat="identity") + scale_fill_gradientn( colours=colorz_solid ) + ylab( paste0("Spearman R with CNV score") ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	p = p + ylim(c( -max(abs(rdf[,"cnv_score"])),+max(abs(rdf[,"cnv_score"])) ))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) + theme(legend.position="none")
	print(p)
	dev.off()

	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	scores$cs_new = NA
	for (rn in rownames(cs_map)){
		scores[scores$cs==rn,"cs_new"] = cs_map[rn,"alias"]
		scores[scores$cs==rn,"cs_colorz"] = cs_map[rn,"colorz"]
		scores[scores$cs==rn,"cs_colorz_level1"] = cs_map[rn,"colorz_level1"]
	}
	scores_save = scores

	pdf( paste0(OutDir,"CnvScore_across_cs_boxplots.pdf"),2.5,2 )
	x = factor(scores[,"cs_new"],levels=cs_map$alias)
	y = scores$cnv_score
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	plot=dan.boxplots.multipages( x, y, xlab = "", ylab = "CNV score", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = max(y), xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"CnvScore_across_cs_Violins.pdf"),2.5,2 )
	plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 0, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
	dan.densityPlot( paste0(OutDir,"CnvScore_across_cs_Density.pdf"), y, x, groupinglab = "", xlab = paste0("cyto score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = cs_map$colorz, fileWidth = 2.5, fileHeight = 1.5 )

	# agg = aggregate(cnv_score~cs_new,data=scores_save,FUN='mean')

	pdf( paste0(OutDir,"CnvScore_across_SampleType_boxplots.pdf"),2,2 )
	x = factor(scores[,"SampleType"],levels=c( "Primary","Metastasis" ))
	y = scores$cnv_score
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	plot=dan.boxplots.multipages( x, y, xlab = "", ylab = "CNV score", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = max(y), xColors = c( "orange","purple" ), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"CnvScore_across_SampleType_Violins.pdf"),2,2 )
	plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 0, xColors = c( "orange","purple" ), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
	dan.densityPlot( paste0(OutDir,"CnvScore_across_SampleType_Density.pdf"), y, x, groupinglab = "", xlab = paste0("cyto score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = c( "orange","purple" ), fileWidth = 2.5, fileHeight = 1.5 )

	# only primary
	scores = scores[scores$SampleType=="Primary",]
	rdf = data.frame(row.names = order_tps,tps = order_tps, cnv_score = as.numeric(cor(scores$cnv_score,scores[,order_tps],method="spearman")) )
	colorz_solid = colorRampPalette(c("dodgerblue4","gray77","firebrick3"))(100)
	fileName = paste0(OutDir,"onlyPrimary_cor_tps_cnvScore.pdf")
	rdf = rdf[order(rdf[,"cnv_score"]),]
	pdf(fileName, 3,2.2)
	p = ggplot(data=rdf, aes(x=factor(rdf$tps,levels = as.character(rdf$tps)), y=rdf[,"cnv_score"],fill=rdf[,"cnv_score"])) + geom_bar(stat="identity") + scale_fill_gradientn( colours=colorz_solid ) + ylab( paste0("Spearman R with CNV score") ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	p = p + ylim(c( -max(abs(rdf[,"cnv_score"])),+max(abs(rdf[,"cnv_score"])) ))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) + theme(legend.position="none")
	print(p)
	dev.off()

	pdf( paste0(OutDir,"onlyPrimary_CnvScore_across_cs_boxplots.pdf"),2.5,2 )
	x = factor(scores[,"cs_new"],levels=cs_map$alias)
	y = scores$cnv_score
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	plot=dan.boxplots.multipages( x, y, xlab = "", ylab = "CNV score", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = max(y), xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"onlyPrimary_CnvScore_across_cs_Violins.pdf"),2.5,2 )
	plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 0, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
	dan.densityPlot( paste0(OutDir,"onlyPrimary_CnvScore_across_cs_Density.pdf"), y, x, groupinglab = "", xlab = paste0("cyto score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = 0, xlimRight = 0.07, groupingColors = cs_map$colorz, fileWidth = 2.5, fileHeight = 1.5 )

	scores_save = scores
	scores = scores_save#[scores_save$cs_new %in% c( "Proliferative","Hypoxic" ),]
	patient_level_df = dan.df(unique(scores$Patient),c( "winner","loser","wt_pval","median_proliferative","median_hypoxic" ),NA)
	for (rn in rownames(patient_level_df)){
		dcat(rn)
		sc = scores[scores$Patient==rn,]
		patient_level_df[rn,"median_proliferative"] = median(sc[sc$cs_new=="Proliferative","cnv_score"])
		patient_level_df[rn,"median_hypoxic"] = median(sc[sc$cs_new=="Hypoxic","cnv_score"])
		patient_level_df[rn,"mean_proliferative"] = mean(sc[sc$cs_new=="Proliferative","cnv_score"])
		patient_level_df[rn,"mean_hypoxic"] = mean(sc[sc$cs_new=="Hypoxic","cnv_score"])
		patient_level_df[rn,"n_proliferative"] = length(sc[sc$cs_new=="Proliferative","cnv_score"])
		patient_level_df[rn,"n_hypoxic"] = length(sc[sc$cs_new=="Hypoxic","cnv_score"])
		for (tp in order_tps){
			patient_level_df[rn,tp] = as.numeric(cor(sc$cnv_score,sc[,tp],method="pearson"))
		}
		# if ((length(dtable(sc$cs_new))<2) | ( sum(dtable(sc$cs_new)>100)<2 )){ next }
		# agg = aggregate(cnv_score~cs_new,data=sc,FUN='median')
		# patient_level_df[rn,"winner"] = agg[agg$cnv_score==max(agg$cnv_score),"cs_new"]
		# patient_level_df[rn,"loser"] = agg[agg$cnv_score==min(agg$cnv_score),"cs_new"]
		# patient_level_df[rn,"wt_pval"] = wilcox.test(cnv_score~cs_new,data=sc)$p.value
		# # patient_level_df[rn,"wt_pval"] = kruskal.test(cnv_score~cs_new,data=sc)$p.value
		# patient_level_df[rn,"Ncells"] = nrow(sc)

	}
	pld = patient_level_df[,order_tps]
	dtable(colnames(pld)[apply(pld,1,which.max)])
	colMeans(patient_level_df[,order_tps])

	wm_proliferative = weighted.mean(patient_level_df$mean_proliferative,patient_level_df$n_proliferative)
	wm_hypoxic = weighted.mean(patient_level_df$mean_hypoxic,patient_level_df$n_hypoxic)
	patient_level_df = patient_level_df[order(patient_level_df$n_hypoxic),]

	mean(patient_level_df$mean_proliferative,na.rm=T)
	mean(patient_level_df$mean_hypoxic,na.rm=T)

	patient_level_df = patient_level_df[!is.na(patient_level_df$wt_pval),]
	patient_level_df = patient_level_df[order(patient_level_df$wt_pval),]

	### CIN signature scoring
	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelial.RData"))
	ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
	mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
	ribo.genes = ribo$Approved.symbol
	seu = subset(x = seu, features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))], cells = seu@meta.data[seu@meta.data$Epi_Cancer=="Cancer","CellID"]) # 12677 x 101849

	cin_signature = dan.read("data/cin_signature_Bakhoum2018.txt")
	cin_signature = cin_signature$cin_signature
	cin_signature = list(cin_signature=cin_signature)
	db_rand = dan.barkley_MakeRand(seu,cin_signature, 3)
	cin_scores = dan.Barkley_GeneToEnrichment_AmsCentered(seu, cin_signature, db_rand)
	colnames(cin_scores)[substr(colnames(cin_scores),1,5)=="enric"] = "cin_signature"
	save(cin_scores,file = paste0(OutDir,"cin_signature_scores.RData"))
	load(file = paste0(OutDir,"cin_signature_scores.RData"))
	scores = scores_save
	scores$cin_signature = cin_scores[rownames(scores),"cin_signature"]

	rdf = data.frame(row.names = order_tps,tps = order_tps, cin_signature = as.numeric(cor(scores$cin_signature,scores[,order_tps],method="spearman")) )
	colorz_solid = colorRampPalette(c("dodgerblue4","gray77","firebrick3"))(100)
	fileName = paste0(OutDir,"cor_tps_cin_signature.pdf")
	rdf = rdf[order(rdf[,"cin_signature"]),]
	pdf(fileName, 3,2.2)
	p = ggplot(data=rdf, aes(x=factor(rdf$tps,levels = as.character(rdf$tps)), y=rdf[,"cin_signature"],fill=rdf[,"cin_signature"])) + geom_bar(stat="identity") + scale_fill_gradientn( colours=colorz_solid ) + ylab( paste0("Spearman R with CIN signature") ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	p = p + ylim(c( -max(abs(rdf[,"cin_signature"])),+max(abs(rdf[,"cin_signature"])) ))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) + theme(legend.position="none")
	print(p)
	dev.off()

	cs_map = data.frame(row.names=c("cs1","cs2","cs3"),alias=c("Alveolar","Proliferative","Hypoxic" ),colorz = c("dodgerblue4","firebrick3","sienna4"), colorz_level1 = c("dodgerblue4","firebrick4","firebrick4"),stringsAsFactors=F)
	scores$cs_new = NA
	for (rn in rownames(cs_map)){
		scores[scores$cs==rn,"cs_new"] = cs_map[rn,"alias"]
		scores[scores$cs==rn,"cs_colorz"] = cs_map[rn,"colorz"]
		scores[scores$cs==rn,"cs_colorz_level1"] = cs_map[rn,"colorz_level1"]
	}

	pdf( paste0(OutDir,"cin_signature_across_cs_boxplots.pdf"),2.5,2 )
	x = factor(scores[,"cs_new"],levels=cs_map$alias)
	y = scores$cin_signature
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	plot=dan.boxplots.multipages( x, y, xlab = "", ylab = "CIN signature", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = max(y), xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	dan.densityPlot( paste0(OutDir,"cin_signature_across_cs_Density.pdf"), y, x, groupinglab = "", xlab = paste0("cyto score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = cs_map$colorz, fileWidth = 2.5, fileHeight = 1.5 )
}

extract_ambient_genes = function( min_cw = 0.2, min_cor = 0.3, include_ig = TRUE){
	load(file = paste0("processed/CellBender_AmbientGenesExtraction/cordf2_all.RData"))
	# in the meantime, I updated the epithelial cell definition. But here I am excluding them anyway, so it doesn't need to be updated.
	bb2 = cordf2_all[(!(cordf2_all$which_max_cor %in% c("n_level2_AT1","n_level2_AT2","n_level2_KAC","n_level2_basal","n_level2_ciliated","n_level2_club","n_level2_transitional_club_AT2")) ) & (((cordf2_all$CB_wm>=min_cw) %in% c(T) ) | ( is.na(cordf2_all$CB_wm) )),]
	bb2 = bb2[bb2$max_cor>min_cor,]
	ambient_genes = sort(rownames(bb2))
	if (include_ig){ ambient_genes = sort(unique(c(ambient_genes, rownames(cordf2_all)[substr(rownames(cordf2_all),1,3) %in% c( "IGH","IGL","IGK" ) ]) )) }
	return(ambient_genes)
}

cvNMF_AdjMat_Construction = function(OutDir, exclude_genes = NULL){
	trialPrefix = "SCTcenteredRmRBMT_"
	DataDir = "data/"
	dataset = "qkwlxb"
	load( file = paste0(DataDir,dataset,"_annot.RData"))
	annot = annot[(annot$Epi_Cancer=="Cancer") & (annot$TN=="Tumor"),]
	annot$CellID = rownames(annot)
	load(file = paste0("processed/scNMF_CV/",trialPrefix,"_rank5/snp0_genes.RData"))
	flat_list = list()
	flat_index = 1
	for (k in c(5:100))
	{
	   this_file = paste0("processed/scNMF_CV/",trialPrefix,"_rank",k,"/wrobust.RData")
	   if (!file.exists(this_file)) { next }
	   dcat(k)
	   load(this_file)
	   w = wrobust
	   a = dan.ExtractGenes2(w, snp0_genes, n = 100, method = "max_NMF_package", exclude_genes = exclude_genes)
	   for (nn in names(a))
	   {
	      if (any(is.na(a[[nn]]))) { a[[nn]] = NULL }
	   }
	   for (ii in names(a))
	   {
	      flat_list[[ flat_index ]] = a[[ii]]
	      flat_index = flat_index+1
	   }
	}
	# Constructing huge matrix
	flat_list = flat_list[sapply(flat_list,length)>1]
	allGenez = unique(as.character(unlist(flat_list)))
	dcat(length(allGenez))
	adjmat = matrix(0,nrow = length(allGenez), ncol = length(allGenez), dimnames = list(allGenez,allGenez))
	index = 1
	for (g1 in allGenez)
	{
	   dcat(paste0(index, ", out of ", length(allGenez)))
	   for (g2 in allGenez)
	   {
	      if (g1!=g2)
	      {
	         adjmat[g1,g2] = sum(sapply(flat_list, function(x) g1 %in% x) & sapply(flat_list, function(x) g2 %in% x))
	      }
	   }
	   index = index+1
	}
	save(adjmat, file = paste0(OutDir,"full_adjmat.RData"))
}

adjmat_clustering = function(OutDir, adjmat_filtering = NULL){
	load(file = paste0("processed/scNMF_CV/SCTcenteredRmRBMT__rank5/snp0_genes.RData"))
	this_OutDir = paste0(OutDir,"louvain_clustering/")
	dir.create(this_OutDir)
	load(file = paste0(OutDir,"full_adjmat.RData"))
	if (!is.null(adjmat_filtering)){
		if (adjmat_filtering=="at_least_in_five_lfs"){
			adjmat2 = adjmat
			adjmat2[adjmat2<5] = 0
			namez = sort(names(which(rowSums(adjmat2)>0)))
			adjmat = adjmat2[namez,namez]
		}
	}
	set.seed(123)
	gw = graph_from_adjacency_matrix(adjmat, mode = "undirected", weighted = T)
	comm_infomap = cluster_infomap(gw)
	comm_louvain = cluster_louvain(gw)
	modularity_df = data.frame(row.names = c("infomap","louvain1","louvain2" ), clustering = c("infomap","louvain1","louvain2" ), modularity = c(comm_infomap$modularity,comm_louvain$modularity[1],comm_louvain$modularity[2]) , stringsAsFactors = F)
	commdf = data.frame(gene = comm_infomap$names, infomap = comm_infomap$membership,stringsAsFactors = F)
	for (i in 1:nrow(comm_louvain$memberships)){
	   commdf[,paste0("louvain", i)] = comm_louvain$memberships[i,]
	}
	rownames(commdf) = commdf$gene

	for (res in seq(0.1,2,0.1)){
	   print(res)
	   nt = cluster_louvain(gw,resolution=res)
	   for (rnn in 1:nrow(nt$memberships)){ 
	      commdf[nt$names,paste0("ig", res,"_sol",rnn)] = as.numeric(nt$memberships[rnn,])
	      modularity_df[paste0("ig", res,"_sol",rnn),"clustering"] = paste0("ig", res,"_sol",rnn)
	      modularity_df[paste0("ig", res,"_sol",rnn),"modularity"] = max(nt$modularity[rnn])
	   }
	}	
	### clustree analysis
	clustree_df = commdf[,c( "infomap" , "louvain1","louvain2")]
	colnames(clustree_df) = paste0("algo",c(1:3) )
	pdf(paste0(this_OutDir,"algorithms_clustree.pdf"),12,12)
	print(clustree(clustree_df, prefix = "algo"))
	dev.off()
	this_commdf = commdf
	cnn = colnames(commdf)
	this_commdf = this_commdf[,cnn[substr(cnn,nchar(cnn)-3,nchar(cnn))=="sol1"]]
	colnames(this_commdf) = gsub( "_sol1","",colnames(this_commdf) )
	clustree_df = this_commdf
	pdf(paste0(this_OutDir,"louvain_nt_res_clustree_sol1.pdf"),12,12)
	print(clustree(clustree_df, prefix = "ig"))
	dev.off()
	pdf(paste0(this_OutDir,"louvain_nt_res_clustree_stability_sol1.pdf"),12,12)
	print(clustree(clustree_df, prefix = "ig", node_colour = "sc3_stability"))
	dev.off()
	# compute average stability
	stab = clustree(clustree_df, prefix = "ig", return = "graph")
	stab2 <-
	  stab %>%
	  activate(nodes) %>%
	  data.frame()
	stabdf = data.frame(row.names = paste0("ig", unique(stab2$ig)), ig = paste0("ig", unique(stab2$ig)), avg_stability = NA)
	for (ig in unique(stab2$ig))
	{
	   stabdf[paste0("ig",ig ),"avg_stability"] = weighted.mean(stab2[stab2$ig==ig,"sc3_stability"], stab2[stab2$ig==ig,"size"])
	}
	save(stabdf, file = paste0( this_OutDir,"stabdf_sol1.RData" ))

	this_commdf = commdf
	this_commdf = this_commdf[,cnn[substr(cnn,nchar(cnn)-3,nchar(cnn))=="sol2"]]
	colnames(this_commdf) = gsub( "_sol2","",colnames(this_commdf) )
	clustree_df = this_commdf
	pdf(paste0(this_OutDir,"louvain_nt_res_clustree_sol2.pdf"),12,12)
	print(clustree(clustree_df, prefix = "ig"))
	dev.off()
	pdf(paste0(this_OutDir,"louvain_nt_res_clustree_stability_sol2.pdf"),12,12)
	print(clustree(clustree_df, prefix = "ig", node_colour = "sc3_stability"))
	dev.off()
	# compute average stability
	stab = clustree(clustree_df, prefix = "ig", return = "graph")
	stab2 <-
	  stab %>%
	  activate(nodes) %>%
	  data.frame()
	stabdf = data.frame(row.names = paste0("ig", unique(stab2$ig)), ig = paste0("ig", unique(stab2$ig)), avg_stability = NA)
	for (ig in unique(stab2$ig))
	{
	   stabdf[paste0("ig",ig ),"avg_stability"] = weighted.mean(stab2[stab2$ig==ig,"sc3_stability"], stab2[stab2$ig==ig,"size"])
	}
	save(stabdf, file = paste0( this_OutDir,"stabdf_sol2.RData" ))

	pdf(paste0(this_OutDir,"modularity_barplot.pdf"),18,8)
	modularity_df$clustering = factor(modularity_df$clustering,levels = as.character(modularity_df$clustering))
	lbl = apply(commdf[,2:ncol(commdf)],2,max)
	modularity_df[names(lbl),"label"] = as.numeric(lbl)
	p1 = ggplot(data=modularity_df, aes(x=clustering, y=modularity)) + geom_bar(stat="identity")+ xlab("") + ylab("\n\nModularity") +
	   geom_text(aes(label=label), vjust=1.6, color="white", size=3.5)+theme_minimal() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
	print(p1)
	dev.off()

	save(commdf, file = paste0(this_OutDir,"commdf.RData"))
	load( file = paste0(this_OutDir,"commdf.RData"))
	for (which_clustering in colnames(commdf)[colnames(commdf)!="gene"]){
	   dir.create(paste0(this_OutDir,which_clustering))
	   wrComm_genes = list()
	   for (ncomm in unique(commdf[,which_clustering]) )
	   {
	      these_wrGenes = commdf[commdf[,which_clustering]==ncomm,"gene"]
	      wrComm_genes[[paste0("c",ncomm )]] = sort(these_wrGenes)
	      save(wrComm_genes,file = paste0(this_OutDir,which_clustering, "/wrComm_genes.RData"))
	   }
	}
	## Comparison with published MPs
	load( file = paste0(this_OutDir,"commdf.RData"))
	for (which_clustering in colnames(commdf)[colnames(commdf)!="gene"]){
	   print(which_clustering)
	   load(file = paste0(this_OutDir,which_clustering,"/wrComm_genes.RData"))
	   dan.compare_published_mp(wrComm_genes, OutFolder = paste0(this_OutDir,which_clustering,"/"), snp0_genes, printAllGenes = T)
	   for (cl in names(wrComm_genes))
	   {
	      write.table(data.frame(gene = wrComm_genes[[cl]]), file = paste0(this_OutDir,which_clustering,"/",cl,".txt"), sep = "\t",quote=F, row.names = F, col.names = F)
	   }
	}
}

get_most_stable_clustering_solution = function(OutDir){
	load(paste0(OutDir,"louvain_clustering/stabdf_sol1.RData"))
	max_sol1 = max(stabdf$avg_stability)
	whichmax_sol1 = stabdf[stabdf$avg_stability==max(stabdf$avg_stability),"ig"]
	load(paste0(OutDir,"louvain_clustering/stabdf_sol2.RData"))
	max_sol2 = max(stabdf$avg_stability)
	whichmax_sol2 = stabdf[stabdf$avg_stability==max(stabdf$avg_stability),"ig"]
	clustering_solution = ifelse(max_sol1>max_sol2,paste0( whichmax_sol1,"_sol1" ),paste0( whichmax_sol1,"_sol2" ))
	return(clustering_solution)
}

tps_splitting = function(OutDir, clustering_solution){
	load(file = paste0("processed/scNMF_CV/SCTcenteredRmRBMT__rank5/snp0_genes.RData"))
	load(file = paste0(OutDir,"louvain_clustering/",clustering_solution,"/wrComm_genes.RData"))
	inter_cs = wrComm_genes
	inter_cs = wrComm_genes[sapply(wrComm_genes,length)>10]
	load(file = paste0(OutDir,"full_adjmat.RData"))
	set.seed(123)
	for (cs in names(inter_cs))
	{
	   dcat(cs)
	   this_adjmat = adjmat[inter_cs[[cs]],inter_cs[[cs]]]
	   gw = graph_from_adjacency_matrix(this_adjmat, mode = "undirected", weighted = T)
	   comm_louvain = cluster_louvain(gw)
	   commdf = data.frame(gene = comm_louvain$names, clusters = comm_louvain$membership,stringsAsFactors = F)
	   rownames(commdf) = commdf$gene
	   wrComm_genes = list()
	   for (ncomm in unique(commdf[,"clusters"]) ){ wrComm_genes[[paste0(cs,"_c",ncomm )]] = sort(commdf[commdf[,"clusters"]==ncomm,"gene"]) }
	   dan.compare_published_mp(wrComm_genes, OutFolder = paste0(OutDir,"louvain_clustering/","splitting_",cs,"/"),snp0_genes, printAllGenes = T)
	   for (cl in names(wrComm_genes)){ write.table(data.frame(gene = wrComm_genes[[cl]]), file = paste0(OutDir,"louvain_clustering/","splitting_",cs,"/",cl,".txt"), sep = "\t",quote=F, row.names = F, col.names = F) }
	   pdf(paste0(paste0(OutDir,"louvain_clustering/","splitting_",cs,"/"),cs,"_adjmat_connectingGenes.pdf"),20,20)
	   heatmap.2(this_adjmat,tracecol = NA, col = colorz_solid)
	   dev.off()
	   dcat(sum(this_adjmat)/((nrow(this_adjmat)^2)-nrow(this_adjmat) ))
	}
}

gsea_enrichment = function(OutDir, clustering_solution, ambient_genes){
	library(msigdbr)
	library(fgsea)
	load(file = paste0("processed/scNMF_CV/SCTcenteredRmRBMT__rank5/snp0_genes.RData"))
	gene_universe = snp0_genes[!(snp0_genes %in% ambient_genes)]
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	msig_df_C6 = msigdbr::msigdbr(species = "Homo sapiens", category = "C6")
	msig_df_C8 = msigdbr::msigdbr(species = "Homo sapiens", category = "C8")
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP,msig_df_C6,msig_df_C8 ))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	## in clustering_solution
	all_files = list.files(paste0(OutDir,"louvain_clustering/",clustering_solution,"/" ),pattern="*txt")
	all_files = all_files[!grepl("MSigDB",all_files )]
	for (f in all_files){
		genez = dan.read(paste0(OutDir,"louvain_clustering/",clustering_solution,"/",f ),header=F)
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = genez$V1,
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0( OutDir,"louvain_clustering/",clustering_solution,"/",gsub(".txt","",f),"_MSigDB_complete.txt" ))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,2)
		dan.write(a,file=paste0( OutDir,"louvain_clustering/",clustering_solution,"/",gsub(".txt","",f),"_MSigDB.txt" ))
	}
	## in splitted solutions
	all_splitted_folders = list.files(paste0(OutDir,"louvain_clustering/" ),pattern="splitting_*",include.dirs = TRUE)
	for (splitted in all_splitted_folders){
		all_files = list.files(paste0(OutDir,"louvain_clustering/",splitted,"/" ),pattern="*txt")
		all_files = all_files[!grepl("MSigDB",all_files )]
		for (f in all_files){
			genez = dan.read(paste0(OutDir,"louvain_clustering/",splitted,"/",f ),header=F)
			fgRes = fgsea::fora(pathways = msig_list,
	                       genes = genez$V1,
	                       universe = gene_universe)
			fgRes = fgRes[fgRes$padj <= 0.01,]
			a = as.data.frame(fgRes)
			a = apply(a,2,as.character)
			dan.write(a,file=paste0( OutDir,"louvain_clustering/",splitted,"/",gsub(".txt","",f),"_MSigDB_complete.txt" ))
			a = as.data.frame(fgRes)
			a$overlapGenes = NULL
			a$pval = NULL
			a$padj = signif(a$padj,2)
			dan.write(a,file=paste0( OutDir,"louvain_clustering/",splitted,"/",gsub(".txt","",f),"_MSigDB.txt" ))
		}
	}
}

tps_validation = function(OutDir, exclude_genes){
	suffix = "_unfiltered"
	load(file = paste0(OutDir,"tps_unvalidated.RData"))
	load(paste0("processed/scNMF_intra/qkwlxbyzsch_SCTcenteredRmRBMT_intra/all_models.RData"))
	load( file = paste0("data/qkwlxbyzsh_annot.RData"))
	these_patients = unique(annot[annot$Dataset %in% c( "yang","zhu","salcher","he" ),"Patient"])
	all_models = all_models[these_patients]
	cs_list = tps
	selected_cs = names(cs_list)
	load(paste0("processed/scNMF_intra/qkwlxbyzsch_SCTcenteredRmRBMT_intra/snp0_genes.RData"))
	filtered_models = list()
	flat_models = list()
	for (p in names(all_models))
	{
	   dcat(p)
	   kmodels = all_models[[p]]
	   model = kmodels[["k20"]]
	   w = model@w
	   w = data.frame(w)
	   a = dan.ExtractGenes2(w, snp0_genes, n = 100, method = "max_NMF_package", exclude_genes = exclude_genes)
	   a = a[sapply(a, length)>1]
	   if (length(a)==0) { next }
	   filtered_models[[p]] = a
	   names(a) = paste0(p, "_", names(a))
	   flat_models = c(flat_models,a)
	}
	all(selected_cs %in% names(cs_list))
	cs_list = cs_list[selected_cs]
	intra_inter_mat = matrix(0, nrow = length(cs_list), ncol = length(flat_models), dimnames = list(names(cs_list),names(flat_models)))
	intra_inter_mat_hypergeometric = matrix(0, nrow = length(cs_list), ncol = length(flat_models), dimnames = list(names(cs_list),names(flat_models)))
	universe = length(snp0_genes)
	for (rn in rownames(intra_inter_mat))
	{
	   for (cn in colnames(intra_inter_mat))
	   {
	      intra_inter_mat[rn,cn] = length(intersect(cs_list[[rn]],flat_models[[cn]]))/length(union(cs_list[[rn]],flat_models[[cn]]))
	      mymp_genes = length(cs_list[[rn]])
	      overlap = length(intersect(cs_list[[rn]],flat_models[[cn]]))
	      othermp_genes = length(flat_models[[cn]])
	      intra_inter_mat_hypergeometric[rn,cn] = -log10(phyper(overlap-1, mymp_genes, universe-mymp_genes, othermp_genes, lower.tail = FALSE, log.p = FALSE))
	   }
	}
	intra_inter_mat = intra_inter_mat[,colSums(intra_inter_mat>0)>0] # from 605 to 560
	assigned = rownames(intra_inter_mat)[apply(intra_inter_mat,2,which.max)]
	names(assigned) = names(apply(intra_inter_mat,2,which.max))
	rowOrder = names(sort(table(assigned),decreasing=T))
	intra_inter_mat = intra_inter_mat[rowOrder,]
	intra_inter_mat2 = intra_inter_mat
	colnames(intra_inter_mat2) = paste0("c",1:ncol(intra_inter_mat) )
	index = 1
	for (cs in rownames(intra_inter_mat2)){
	   cs_assigned = names(assigned[assigned==cs])
	   if (length(cs_assigned)>1) { cs_assigned = names(sort(intra_inter_mat[cs,cs_assigned],decreasing=T)) }
	   for ( col in cs_assigned ){
	      intra_inter_mat2[,paste0("c",index)] = intra_inter_mat[,col]
	      colnames(intra_inter_mat2)[colnames(intra_inter_mat2)==paste0("c",index)] = col
	      index = index+1
	   }
	}
	pdf(paste0(OutDir,"intercs",suffix,"_Cs_vs_singlePrograms_jaccard0.pdf"),17,10)
	heatmap.2(intra_inter_mat2,tracecol = NA, dendrogram = 'none', Rowv = F, margins=c(10,10),Colv = F, col = colorRampPalette(c("white","aquamarine4"))(200))
	dev.off()
	flat_models = flat_models[colnames(intra_inter_mat2)]
	intra_mat = matrix(0, nrow = length(flat_models), ncol = length(flat_models), dimnames = list(names(flat_models),names(flat_models)))
	for (rn in rownames(intra_mat))
	{
	   for (cn in colnames(intra_mat))
	   {
	      intra_mat[rn,cn] = length(intersect(flat_models[[rn]],flat_models[[cn]]))/length(union(flat_models[[rn]],flat_models[[cn]]))
	   }
	}
	colz = dan.colors( length(rownames(intra_inter_mat2)) )
	names(colz) = rownames(intra_inter_mat2)
	ColSideColors = colnames(intra_mat)
	ColSideColors2 = ColSideColors
	for (c in ColSideColors ){ ColSideColors2[ColSideColors2==c] = colz[assigned[c]] }
	pdf(paste0(OutDir,"intercs",suffix,"_singlePrograms_vs_singlePrograms_jaccard0.pdf"),30,30)
	heatmap.2(intra_mat,tracecol = NA, dendrogram = 'none', Rowv = F, margins=c(10,10),Colv = F, col = colorRampPalette(c("white","aquamarine4"))(200),RowSideColors = ColSideColors2)
	dev.off()

	save(intra_inter_mat_hypergeometric,file=paste0(OutDir,"intra_inter_mat_hypergeometric",suffix,".RData"))
	load(file=paste0(OutDir,"intra_inter_mat_hypergeometric",suffix,".RData"))
	intra_inter_mat = intra_inter_mat_hypergeometric
	intra_inter_mat = intra_inter_mat[,colSums(intra_inter_mat>0)>0] # from 605 to 560
	assigned = rownames(intra_inter_mat)[apply(intra_inter_mat,2,which.max)]
	names(assigned) = names(apply(intra_inter_mat,2,which.max))
	rowOrder = names(sort(table(assigned),decreasing=T))
	intra_inter_mat = intra_inter_mat[rowOrder,]
	intra_inter_mat2 = intra_inter_mat
	colnames(intra_inter_mat2) = paste0("c",1:ncol(intra_inter_mat) )
	index = 1
	for (cs in rownames(intra_inter_mat2)){
	   cs_assigned = names(assigned[assigned==cs])
	   if (length(cs_assigned)>1) { cs_assigned = names(sort(intra_inter_mat[cs,cs_assigned],decreasing=T)) }
	   for ( col in cs_assigned ){
	      intra_inter_mat2[,paste0("c",index)] = intra_inter_mat[,col]
	      colnames(intra_inter_mat2)[colnames(intra_inter_mat2)==paste0("c",index)] = col
	      index = index+1
	   }
	}
	library(ComplexHeatmap)
	library(viridis)
	rev_magma_white = rev(magma(10))
	rev_magma_white[1] = "white"
	rev_magma_white = colorRampPalette(rev_magma_white)
	rev_magma_white = rev_magma_white(200)
	pdf(paste0(OutDir,"intercs",suffix,"_Cs_vs_singlePrograms_hypergeometric_ComplexHeatmap.pdf"),4,2.5)
	plot=ComplexHeatmap::Heatmap(intra_inter_mat2,rect_gp = gpar(col = "white", lwd = .1),cluster_rows=FALSE,column_title_side="bottom",column_title="Patient-wise single NMF latent factors (k=20)",column_title_gp = gpar(fontsize = 6),cluster_columns=FALSE, column_names_rot = 45,use_raster = F,raster_quality=10, col=rev_magma_white, heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="leftcenter",title="-log10(p-value)",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),show_column_names = FALSE,row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top")
	dev.off()

	# pdf(paste0(OutDir, "matt2_ComplexHeatmap.pdf"),ncol(matt2)/3,nrow(matt2)/3)
	# dev.off()

	pdf(paste0(OutDir,"intercs",suffix,"_Cs_vs_singlePrograms_hypergeometric.pdf"),17,10)
	heatmap.2(intra_inter_mat2,tracecol = NA, dendrogram = 'none', Rowv = F, margins=c(10,10),Colv = F, col = colorRampPalette(c("white","aquamarine4"))(200))
	dev.off()
	flat_models = flat_models[colnames(intra_inter_mat2)]
	intra_mat = matrix(0, nrow = length(flat_models), ncol = length(flat_models), dimnames = list(names(flat_models),names(flat_models)))
	for (rn in rownames(intra_mat))
	{
	   for (cn in colnames(intra_mat))
	   {
	      intra_mat[rn,cn] = length(intersect(flat_models[[rn]],flat_models[[cn]]))/length(union(flat_models[[rn]],flat_models[[cn]]))
	   }
	}
	ColSideColors = colnames(intra_mat)
	ColSideColors2 = ColSideColors
	for (c in ColSideColors ){ ColSideColors2[ColSideColors2==c] = colz[assigned[c]] }
	pdf(paste0(OutDir,"intercs",suffix,"_singlePrograms_vs_singlePrograms_hypergeometric.pdf"),30,30)
	heatmap.2(intra_mat,tracecol = NA, dendrogram = 'none', Rowv = F, margins=c(10,10),Colv = F, col = colorRampPalette(c("white","aquamarine4"))(200),RowSideColors = ColSideColors2)
	dev.off()


	intra_inter_mat_hypergeometric_bin = 1*(intra_inter_mat_hypergeometric>(-log10(0.05)))
	table(as.numeric(colSums(intra_inter_mat_hypergeometric_bin)))
	validated = sort(names(sort(rowSums(intra_inter_mat_hypergeometric_bin)[rowSums(intra_inter_mat_hypergeometric_bin)>5])))
	tps = tps[validated]
	save(tps, file=paste0(OutDir,"tps.RData"))
	tps_universe = snp0_genes
	save(tps_universe, file = paste0(OutDir,"tps_universe.RData")) # this is the universe that includes all datasets. Not the one for tps discovery (for enrichment analysis)

	dir.create(paste0(OutDir,"tps_enrichments/"))
	library(msigdbr)
	library(fgsea)
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	# msig_df_C6 = msigdbr::msigdbr(species = "Homo sapiens", category = "C6")
	# msig_df_C8 = msigdbr::msigdbr(species = "Homo sapiens", category = "C8")
	# msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP,msig_df_C6,msig_df_C8 ))
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP ))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	for (tp in names(tps)){
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = tps[[tp]],
                       universe = as.character(tps_universe))
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0( OutDir,"tps_enrichments/",tp,"_MSigDB_complete.txt" ))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,5)
		dan.write(a,file=paste0( OutDir,"tps_enrichments/",tp,"_MSigDB.txt" ))
	}
	dan.compare_published_mp(tps, OutFolder = paste0( OutDir,"tps_enrichments/withAtlases/"), tps_universe, printAllGenes = T)

	load(file = paste0("processed/scNMF_CV/SCTcenteredRmRBMT__rank5/snp0_genes.RData"))
	tps_universe = snp0_genes
	dir.create(paste0(OutDir,"tps_enrichments_discovery_universe/"))
	library(msigdbr)
	library(fgsea)
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	# msig_df_C6 = msigdbr::msigdbr(species = "Homo sapiens", category = "C6")
	# msig_df_C8 = msigdbr::msigdbr(species = "Homo sapiens", category = "C8")
	# msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP,msig_df_C6,msig_df_C8 ))
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP ))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	for (tp in names(tps)){
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = tps[[tp]],
                       universe = as.character(tps_universe))
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0( OutDir,"tps_enrichments_discovery_universe/",tp,"_MSigDB_complete.txt" ))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,5)
		dan.write(a,file=paste0( OutDir,"tps_enrichments_discovery_universe/",tp,"_MSigDB.txt" ))
	}
	dan.compare_published_mp(tps, OutFolder = paste0( OutDir,"tps_enrichments_discovery_universe/withAtlases/"), tps_universe, printAllGenes = T)
}

tps_enrichments_barplots = function( OutDir ){
	tps_map = data.frame(row.names=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		tps=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		colorz = c( "midnightblue","dodgerblue4","dodgerblue3","dodgerblue2","darkgreen","forestgreen","green3","firebrick4","firebrick3","sienna1","sienna3","sienna4","goldenrod1","goldenrod3","gray25","gray55","gray85","magenta4","magenta3","lightpink3","lightpink" ),rankz=c(1:21),stringsAsFactors=F)
	all_enrich = dan.df( 0,c( "tps","colorz","term","padj" ) )
	all_enrich_unfiltered = dan.df( 0,c( "tps","colorz","term","padj" ) )
	for (tp in rownames(tps_map)){
		load(paste0(OutDir,tp,"_vs_PublishedMps.RData"))
		if (nrow(thisdf)>0){
			thisdf = thisdf[order(thisdf$qvals),]
			all_enrich_unfiltered = rbind(all_enrich_unfiltered, data.frame(tps=tp,colorz=tps_map[tp,"colorz"],term=thisdf$mp,padj=thisdf$qvals,stringsAsFactors=F) )
			thisdf = thisdf[1:min(c(5,nrow(thisdf))),]
			all_enrich = rbind(all_enrich, data.frame(tps=tp,colorz=tps_map[tp,"colorz"],term=thisdf$mp,padj=thisdf$qvals,stringsAsFactors=F) )
		}
		thisdf = dan.read(paste0(OutDir,tp,"_MSigDB.txt"))
		if (nrow(thisdf)>0){
			thisdf = thisdf[order(thisdf$padj),]
			all_enrich_unfiltered = rbind(all_enrich_unfiltered, data.frame(tps=tp,colorz=tps_map[tp,"colorz"],term=thisdf$pathway,padj=thisdf$padj,stringsAsFactors=F) )
			thisdf = thisdf[1:min(c(5,nrow(thisdf))),]
			all_enrich = rbind(all_enrich, data.frame(tps=tp,colorz=tps_map[tp,"colorz"],term=thisdf$pathway,padj=thisdf$padj,stringsAsFactors=F) )
		}
	}
	all_enrich_unfiltered$minusLog10qval = -log10(all_enrich_unfiltered$padj)
	all_enrich_unfiltered = all_enrich_unfiltered[order(all_enrich_unfiltered$tps,all_enrich_unfiltered$padj),]
	save(all_enrich_unfiltered,file=paste0( OutDir,"all_enrichments_unfiltered.RData" ))

	all_enrich$minusLog10qval = -log10(all_enrich$padj)
	all_enrich = all_enrich[order(all_enrich$tps,all_enrich$padj),]
	all_enrich$xlabz = substr(as.character(all_enrich$term),1,20)
	all_enrich$term = paste0("tps",tps_map[all_enrich$tps,"rankz"],"_",as.character(all_enrich$term))
	all_enrich$term = factor(all_enrich$term,levels=as.character(all_enrich$term) )
	all_enrich$tps = factor(all_enrich$tps,levels=tps_map$tps)
	# all_enrich = all_enrich[order(all_enrich$tps,all_enrich$padj),]
	pdf( paste0( OutDir,"barplots_all_enrichments.pdf" ),30,30 )
	plot = ggplot(all_enrich, aes(x=term,y=minusLog10qval, fill=tps)) + geom_bar(stat='identity') + facet_wrap(~tps,scale="free_x") + scale_fill_manual(values=tps_map$colorz) + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=all_enrich$term,labels=all_enrich$xlabz) + ylab("-log10(q-value)" )
	print(plot)
	dev.off()
	dan.save(all_enrich,paste0( OutDir,"barplots_all_enrichments.pdf" ))	

	## Selected
	selected_tps = c( "AT2-like","AT2-Club-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Stress_AP1","OxPhos","Unfolded_protein_response" )
	all_enrich = dan.df( 0,c( "tps","colorz","term","padj" ) )
	for (tp in selected_tps){
		this_enrich = dan.df( 0,c( "tps","colorz","term","padj" ) )
		load(paste0(OutDir,tp,"_vs_PublishedMps.RData"))
		if (nrow(thisdf)>0){
			thisdf = thisdf[order(thisdf$qvals),]
			thisdf = thisdf[1:min(c(5,nrow(thisdf))),]
			this_enrich = rbind(this_enrich, data.frame(tps=tp,colorz=tps_map[tp,"colorz"],term=thisdf$mp,padj=thisdf$qvals,stringsAsFactors=F) )
		}
		thisdf = dan.read(paste0(OutDir,tp,"_MSigDB.txt"))
		if (nrow(thisdf)>0){
			thisdf = thisdf[order(thisdf$padj),]
			thisdf = thisdf[1:min(c(5,nrow(thisdf))),]
			this_enrich = rbind(this_enrich, data.frame(tps=tp,colorz=tps_map[tp,"colorz"],term=thisdf$pathway,padj=thisdf$padj,stringsAsFactors=F) )
		}
		this_enrich = this_enrich[order(this_enrich$padj),]
		this_enrich = this_enrich[1:3,]
		all_enrich = rbind(all_enrich,this_enrich)
	}
	all_enrich$minusLog10qval = -log10(all_enrich$padj)
	all_enrich = all_enrich[order(all_enrich$tps,all_enrich$padj),]
	all_enrich$xlabz = substr(as.character(all_enrich$term),1,20)
	all_enrich$term = paste0("tps",tps_map[all_enrich$tps,"rankz"],"_",as.character(all_enrich$term))
	all_enrich$term = factor(all_enrich$term,levels=as.character(all_enrich$term) )
	all_enrich$tps = factor(all_enrich$tps,levels=tps_map$tps[tps_map$tps %in% all_enrich$tps ] )
	pdf( paste0( OutDir,"barplots_selected_enrichments.pdf" ),16,5 )
	plot = ggplot(all_enrich, aes(x=term,y=minusLog10qval, fill=tps)) + geom_bar(stat='identity') + facet_wrap(~tps,scale="free_x",nrow=1) + scale_fill_manual(values=tps_map[levels(all_enrich$tps),"colorz"]) + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=all_enrich$term,labels=all_enrich$xlabz) + ylab("-log10(q-value)" )
	print(plot)
	dev.off()
	dan.save(all_enrich,paste0( OutDir,"barplots_selected_enrichments.pdf" ))
}

tps_enrichments_heatmap = function( OutDir,order_tps ){
	load(paste0(OutDir,"../tps.RData"))
	selected_mp = c( "gavish_MP31.Alveolar","lepidic_augm","trav_club","hlca_AT0","hlca_AT2","hlca_preTB","trav_alveolar_epithelial_type_2","trav_proliferating_basal","trav_proximal_basal","barkley_Cycle",
	"gavish_MP1..Cell.Cycle...G2.M","hlca_AT2_prolif","hlca_club","gavish_MP13.EMT.II","barkley_pEMT","barkley_Hypoxia","gavish_MP6.Hypoxia","barkley_Interferon","gavish_MP17.Interferon.MHC.II..I.","hlca_AT1",
	"gavish_MP39.Metal.response","barkley_Metal","gavish_MP18.Interferon.MHC.II..II.","gavish_MP21.Respiration","barkley_Oxphos","gavish_MP14.EMT.III","gavish_MP41.Unassigned.II","barkley_Stress","gavish_MP5.Stress","gavish_MP22.Secreted.I","gavish_MP11.Translation.initiation","gavish_MP9.Unfolded.protein.response","gavish_MP10.Protein.maturation" )
	selected_mp_aliases = c( "MP31 Alveolar (Gavish et al.)","lepidic signature (Tavernari et al.)","Club (Travaglini et al.)","AT0 (HLCA)","AT2 (HLCA)","preTB (HLCA)","Alveolar epithelial type 2 (Travaglini et al.)","Proliferating basal (Travaglini et al.)","Proximal basal (Travaglini et al.)","Cycle (Barkley et al.)",
	"MP1 Cell Cycle G2M (Gavish et al.)","AT2 proliferating (HLCA)","Club (HLCA)","MP13 EMT-II (Gavish et al.)","pEMT (Barkley et al.)","Hypoxia (Barkley et al.)","MP6 Hypoxia (Gavish et al.)","Interferon (Barkley et al.)","MP17 Interferon MHC-II (I) (Gavish et al.)","AT1 (HLCA)",
	"MP39 Metal response (Gavish et al.)","Metal (Barkley et al.)","MP18 Interferon MHC-II (II) (Gavish et al.)","MP21 Respiration (Gavish et al.)","Oxphos (Barkley et al.)","MP14 EMT-III (Gavish et al.)","MP41 Unassigned-II (Gavish et al.)","Stress (Barkley et al.)","MP5 Stress (Gavish et al.)","MP22 Secreted I (Gavish et al.)","MP11 Translation initiation (Gavish et al.)","MP9 Unfolded protein response (Gavish et al.)","MP10 Protein maturation (Gavish et al.)" )
	
	order_tps = c( "AT2-like","AT2-Club-like","lepidic-like","Club-like","MHC-II","MHC-I","Interferon","Stress_AP1","Stress_HSP","Stress_secreted", "OxPhos", "RNA_processing","Unas_emp","Unfolded_protein_response","Translation_initiation","Basal-like","Cell_proliferation","pEMT","EMT","Hypoxia","Metal")
	selected_mp = c( "gavish_MP31.Alveolar","lepidic_augm","trav_club","hlca_AT0","hlca_AT2","hlca_preTB","trav_alveolar_epithelial_type_2","trav_proliferating_basal","trav_proximal_basal","barkley_Cycle",
	"gavish_MP1..Cell.Cycle...G2.M","hlca_club","gavish_MP13.EMT.II","barkley_pEMT","barkley_Hypoxia","gavish_MP6.Hypoxia","barkley_Interferon","gavish_MP17.Interferon.MHC.II..I.","hlca_AT1","gavish_MP39.Metal.response","barkley_Metal",
	"gavish_MP18.Interferon.MHC.II..II.","gavish_MP21.Respiration","barkley_Oxphos","gavish_MP14.EMT.III","gavish_MP41.Unassigned.II","barkley_Stress","gavish_MP5.Stress","gavish_MP22.Secreted.I","gavish_MP11.Translation.initiation","gavish_MP9.Unfolded.protein.response","gavish_MP10.Protein.maturation" )
	selected_mp_aliases = c( "MP31 Alveolar (Gavish et al.)","lepidic signature (Tavernari et al.)","Club (Travaglini et al.)","AT0 (HLCA)","AT2 (HLCA)","preTB (HLCA)","Alveolar epithelial type 2 (Travaglini et al.)","Proliferating basal (Travaglini et al.)","Proximal basal (Travaglini et al.)","Cycle (Barkley et al.)",
	"MP1 Cell Cycle G2M (Gavish et al.)","Club (HLCA)","MP13 EMT-II (Gavish et al.)","pEMT (Barkley et al.)","Hypoxia (Barkley et al.)","MP6 Hypoxia (Gavish et al.)","Interferon (Barkley et al.)","MP17 Interferon MHC-II (I) (Gavish et al.)","AT1 (HLCA)",
	"MP39 Metal response (Gavish et al.)","Metal (Barkley et al.)","MP18 Interferon MHC-II (II) (Gavish et al.)","MP21 Respiration (Gavish et al.)","Oxphos (Barkley et al.)","MP14 EMT-III (Gavish et al.)","MP41 Unassigned-II (Gavish et al.)","Stress (Barkley et al.)","MP5 Stress (Gavish et al.)","MP22 Secreted I (Gavish et al.)","MP11 Translation initiation (Gavish et al.)","MP9 Unfolded protein response (Gavish et al.)","MP10 Protein maturation (Gavish et al.)" )
	
	matt = dan.df(selected_mp,names(tps))
	for (tp in names(tps)){
		load(paste0(OutDir,"withAtlases/",tp,"_vs_PublishedMps_All.RData"))
		rownames(thisdf) = thisdf$mp
		matt[selected_mp,tp] = thisdf[selected_mp,"minusLog10Pvalz"]
	}
	rownames(matt) = selected_mp_aliases
	# mattBin=(matt2>10)*1
	mattClipped = matt
	mattClipped[matt>15] = 15
	mattClipped = reorderMatt(t(mattClipped[,order_tps]),columns_only=T )
	
	# colnames(mattClipped) = selected_mp_aliases
	# Complete, binary, clipped
	library(viridis)
	library(ComplexHeatmap)
	colorz = (viridis(200))
	rev_magma_white = rev(magma(10))
	rev_magma_white[1] = "white"
	rev_magma_white = colorRampPalette(rev_magma_white)
	rev_magma_white = rev_magma_white(200)
	colorz = rev_magma_white
	# pdf(paste0(OutDir, "matt2_ComplexHeatmap.pdf"),ncol(matt2)/3,nrow(matt2)/3)
	# print(ComplexHeatmap::Heatmap(matt2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="-log10(p-value)",direction = "vertical"),column_names_side = "bottom",row_names_side = "left"))
	# dev.off()
	pdf(paste0(OutDir, "mattClipped_ComplexHeatmap.pdf"),5.5,3.5)	
	plot=ComplexHeatmap::Heatmap(mattClipped,rect_gp = gpar(col = "white", lwd = .1),cluster_rows=FALSE,column_title_side="bottom",column_title_gp = gpar(fontsize = 6),cluster_columns=FALSE, column_names_rot = 45,use_raster = F,raster_quality=10, col=colorz, heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="leftcenter",title="-log10(p-value) (clipped to 15)",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),show_column_names = T,row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top")
	dev.off()
	# pdf(paste0(OutDir, "mattBin_ComplexHeatmap.pdf"),ncol(matt2)/3,nrow(matt2)/3)
	# print(ComplexHeatmap::Heatmap(mattBin,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="-log10(p-value)",direction = "vertical"),column_names_side = "bottom",row_names_side = "left"))
	# dev.off()
}

tps_representation = function( OutDir ){
	library(viridis)
	library(gplots)
	tps_map = data.frame(row.names=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		tps=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		colorz = c( "midnightblue","dodgerblue4","dodgerblue3","dodgerblue2","darkgreen","forestgreen","green3","firebrick4","firebrick3","sienna1","sienna3","sienna4","goldenrod1","goldenrod3","gray25","gray55","gray85","magenta4","magenta3","lightpink3","lightpink" ),rankz=c(1:21),stringsAsFactors=F)
	load( paste0(OutDir,"../tps.RData" ) )
	load( paste0(OutDir,"../full_adjmat.RData" ) )
	inter_cs = tps[ c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia",
            "Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" )]
	adjmat = adjmat[as.character(unlist(inter_cs)),as.character(unlist(inter_cs))]
	save(adjmat,file = paste0(OutDir,"adjmat.RData"))
	tps_side_colors = c()
	for (tp in names(inter_cs)){
	   lenn = length(inter_cs[[tp]])
	   tps_side_colors = c(tps_side_colors, rep(tps_map[tp,"colorz"],times=lenn))
	}
	rev_magma_white = rev(magma(10))
	rev_magma_white[1] = "white"
	rev_magma_white = colorRampPalette(rev_magma_white)
	rev_magma_white = rev_magma_white(200)
	pdf(paste0(OutDir,"raw_mat_log.pdf"),14,11)
	heatmap.2(log2(adjmat+1),tracecol = NA, dendrogram = 'none', Rowv = F, margins=c(3,3),lhei=c(2,10), lwid=c(2,3.5),Colv = F, col = rev_magma_white,RowSideColors = tps_side_colors, key.title = "gene-gene network weight", key.xlab="", key.ylab="",labRow = FALSE, xlab = "Marker genes of NMF latent factors", ylab = "Marker genes of NMF latent factors", labCol = FALSE,density.info="none")
	legend("bottomleft",      
	    legend = names(inter_cs),
	    col = tps_map[names(inter_cs),"colorz"], 
	    pch=15,
	    cex=1.5
	    )
	dev.off()

	load( paste0(OutDir,"../tps_unvalidated.RData" ) )
	load( paste0(OutDir,"../full_adjmat.RData" ) )
	inter_cs = tps[ c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia",
            "Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing",names(tps)[substr(names(tps),1,4)=="Unas"] )]
	adjmat = adjmat[as.character(unlist(inter_cs)),as.character(unlist(inter_cs))]
	save(adjmat,file = paste0(OutDir,"adjmat_withUnvalidated.RData"))
	tps_side_colors = c()
	for (tp in names(inter_cs)){
	   lenn = length(inter_cs[[tp]])
	   if (tp %in% rownames(tps_map)){ 
	   	tps_side_colors = c(tps_side_colors, rep(tps_map[tp,"colorz"],times=lenn)) 
	   } else {
	   	tps_side_colors = c(tps_side_colors, rep("gray95",times=lenn)) 
	   }
	}
	rev_magma_white = rev(magma(10))
	rev_magma_white[1] = "white"
	rev_magma_white = colorRampPalette(rev_magma_white)
	rev_magma_white = rev_magma_white(200)
	pdf(paste0(OutDir,"raw_mat_log_withUnvalidated.pdf"),14,11)
	heatmap.2(log2(adjmat+1),tracecol = NA, dendrogram = 'none', Rowv = F, margins=c(3,3),lhei=c(2,10), lwid=c(2,3.5),Colv = F, col = rev_magma_white,RowSideColors = tps_side_colors, key.title = "gene-gene network weight", key.xlab="", key.ylab="",labRow = FALSE, xlab = "Marker genes of NMF latent factors", ylab = "Marker genes of NMF latent factors", labCol = FALSE,density.info="none")
	# legend("bottomleft",      
	#     legend = names(inter_cs),
	#     col = tps_map[names(inter_cs),"colorz"], 
	#     pch=15,
	#     cex=1.5
	#     )
	dev.off()
}

get_common_genez = function( datasets_vector, suffix = "_counts_EpiCancer_harmonized.RData" ){
	first = TRUE
	for (dataset in datasets_vector){
		dcat(dataset)
		load( paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/",dataset,"_counts_EpiCancer_harmonized.RData") )
		if (first){
			commonz = rownames(counts)
			first = FALSE
		} else {
			commonz = intersect(commonz, rownames(counts))
		}
		dcat( paste0("# genes in this dataset  = ",length(rownames(counts)) ) )
		dcat( paste0("# common genes remaining = ",length(commonz) ) )
	}
	return( commonz )
}

extended_atlas_normalization = function(MainDir, exclude_genes = NULL){

	OutDir = MainDir
	commonz = get_common_genez( datasets_vector = c( "qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu" ) )
	dcat(length(commonz))
	if (!is.null(exclude_genes)){ commonz = commonz[!(commonz %in% exclude_genes)] }
	dcat(length(commonz))
	# commonz = 12835. The biggest drop is caused by yang. If this dataset is not considered, commonz would be 13986
	first = T
	for (dataset in c( "qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu" ) ) 
	{
		dcat(dataset)
		load(file=paste0("data/",dataset,"_annot_EpiCancer_harmonized.RData") )
		annot$hlca_prediction = NULL
		load( paste0("data/",dataset,"_counts_EpiCancer_harmonized.RData") )
		if (dataset=="wang"){ 
			annot$CellType_theirpaper = NA
			annot$CellID = paste0( "wang",annot$CellID )
			rownames(annot) = annot$CellID
			colnames(counts) = paste0( "wang",colnames(counts))
		}
		counts = counts[commonz,annot$CellID]
		if (first){
			aal = annot
			coall = counts
			first = F
		} else {
			aal = rbind(aal,annot[,colnames(aal)])
			coall = cbind(coall, counts)
		}
	}
	counts = coall
	annot = aal
	### fixing Sample-Patient initials in "zhu","salcher","chen","he"
	annot$Dataset = tolower(annot$Dataset)
	for (dataset in c("zhu","salcher","he")){
		annot[annot$Dataset==dataset,"Sample"] = paste0( substr(dataset,1,1),annot[annot$Dataset==dataset,"Sample"] )
		annot[annot$Dataset==dataset,"Patient"] = paste0( substr(dataset,1,1),annot[annot$Dataset==dataset,"Patient"] )
	}
	for (dataset in c("wang","hu")){
		annot[annot$Dataset==dataset,"Sample"] = paste0( dataset,annot[annot$Dataset==dataset,"Sample"] )
		annot[annot$Dataset==dataset,"Patient"] = paste0( dataset,annot[annot$Dataset==dataset,"Patient"] )
	}
	
	### Let's keep only cancer cells here
	annot = annot[annot$Epi_Cancer=="Cancer",]
	annot[annot$SampleType=="Tumor","SampleType"] = "Primary" # fixing
	counts = counts[,rownames(annot)]
	### Adding additional clinical information
	load(paste0("data/extendedAtlas_clin.RData"))
	for (sample in clin[clin$Dataset=="wu","Sample"]){
		annot[annot$Sample==sample,"Stage"] = clin[clin$Sample==sample,"Stage"]
		annot[annot$Sample==sample,"Stage_collapsed"] = clin[clin$Sample==sample,"Stage_collapsed"]
	}
	annot[(annot$Stage_collapsed=="III/IV") %in% c(T),"Stage_collapsed"] = NA
	annot[(annot$Stage_collapsed=="III/IV") %in% c(T),"Stage"] = NA
	for (s in unique(annot$Sample)){
		dcat(s)
		for (col in c( "SampleSubType","Smoking","histologic_pattern","Grade","Treatment","mut_TP53","mut_KRAS","mut_KMT2C","mut_PIK3CA","mut_KMT2D","mut_PTEN","mut_KEAP1","mut_STK11","mut_EGFR") ){
			annot[annot$Sample==s,col] = clin[clin$Sample==s,col]
		}
	}
	save(annot, file = paste0(OutDir,"extendedAtlas_annot.RData") )
	save(counts, file = paste0(OutDir,"extendedAtlas_counts.RData") )
	seu = CreateSeuratObject(counts = counts, project = "extendedAtlas", meta.data = annot)
	seu = SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("percent.mt"), verbose = FALSE)
	seu = RunPCA(seu)
	save(seu, file = paste0(OutDir,"extendedAtlas_seu_Cancer.RData"))
}

extended_normal_atlas_normalization = function( MainDir,exclude_genes = ambient_genes ){
	OutDir = MainDir
	commonz = get_common_genez( datasets_vector = c( "zhus0","wangs0" ) )
	if (!is.null(exclude_genes)){ commonz = commonz[!(commonz %in% exclude_genes)] }
	
	load(file = paste0(OutDir,"extendedAtlas_counts.RData") )
	load(file = paste0(OutDir,"extendedAtlas_annot.RData") )
	commonz = intersect( commonz,rownames(counts) ) # commonz cancer
	aal = annot
	coall = counts[commonz,]

	## let's now add annot and counts from the s0 atlas
	for (dataset in c( "zhus0","wangs0" ) ) 
	{
		dcat(dataset)
		load(file=paste0("data/",dataset,"_annot_EpiCancer_harmonized.RData") )
		annot$hlca_prediction = NULL
		annot = annot[annot$Epi_Cancer=="Cancer",]
		annot$Dataset = tolower(annot$Dataset)
		if (dataset=="zhus0"){
			annot[annot$Dataset=="zhu","Sample"] = paste0( substr(dataset,1,1),annot[annot$Dataset=="zhu","Sample"] )
			annot[annot$Dataset=="zhu","Patient"] = paste0( substr(dataset,1,1),annot[annot$Dataset=="zhu","Patient"] )
		}
		load( paste0("data/",dataset,"_counts_EpiCancer_harmonized.RData") )
		if (dataset=="wangs0"){ 
			annot$CellType_theirpaper = NA
			annot$CellID = paste0( "wang",annot$CellID )
			rownames(annot) = annot$CellID
			colnames(counts) = paste0( "wang",colnames(counts))
			annot$Stage_collapsed = "0"
			annot[annot$Dataset=="wang","Sample"] = paste0( "wang",annot[annot$Dataset=="wang","Sample"] )
			annot[annot$Dataset=="wang","Patient"] = paste0( "wang",annot[annot$Dataset=="wang","Patient"] )
		}
		counts = counts[commonz,annot$CellID]
		for (cn in colnames(aal)[!(colnames(aal) %in% colnames(annot))] ){ annot[,cn] = NA }
		aal = rbind(aal,annot[,colnames(aal)])
		coall = cbind(coall, counts)
	}
	annot = aal
	counts = coall
	counts = counts[,rownames(annot)]
	### Adding additional clinical information
	load(paste0("data/extendedAtlasS0_clin.RData"))
	for (s in unique(annot$Sample)){
		dcat(s)
		for (col in c( "SampleSubType","Smoking","histologic_pattern","Grade","Treatment","mut_TP53","mut_KRAS","mut_KMT2C","mut_PIK3CA","mut_KMT2D","mut_PTEN","mut_KEAP1","mut_STK11","mut_EGFR") ){
			annot[annot$Sample==s,col] = clin[clin$Sample==s,col]
		}
	}
	dcat(dim(annot))
	dcat(dim(counts))
	dtable(annot$Dataset,annot$Stage_collapsed)
	dtable(annot$Dataset,annot$SampleType)
	dcat(length(unique(annot$Sample))) # 157
	dcat(length(unique(annot$Patient))) # 146
	### Adding normal cells
	CellTypeDir = paste0("data/cell_typing/non_malignant/")
	epi_cell_types = c( "AT1","AT2","AT0","preTB","club","ciliated","basal" )
	for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","zhus0","hu","wang","wangs0")){
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_seu.RData"))
		ta = seu@meta.data[seu@meta.data$annd_level_2 %in% epi_cell_types,]
		ta$Epi_Cancer = ta$annd_level_2
		for (cn in colnames(annot)[!(colnames(annot) %in% colnames(ta))] ){ ta[,cn] = NA }
		seu = subset(seu, cells = rownames(ta))
		tc = as.matrix(seu@assays$RNA@data)
		if (dataset %in% c( "wang" )){ ta$CellType_theirpaper = NA }
		if (dataset=="wangs0"){
			ta$CellType_theirpaper = NA
			ta$Stage_collapsed = "0"
			ta = ta[!(rownames(ta) %in% rownames(annot)),] # some of the normal samples from wangs0 were already included in wang. In this way I don't include the same cells twice.
			tc = tc[,rownames(ta)]
		}
		annot = rbind(annot,ta[,colnames(annot)])
		commonz = intersect(rownames(tc),rownames(counts))
		print(dtable(ta$annd_level_2))
		dcat(length(commonz))
		counts = cbind(counts[commonz,],tc[commonz,rownames(ta)])
	}
	save(annot, file = paste0(OutDir,"extendedAtlasS0normal_annot.RData") )
	save(counts, file = paste0(OutDir,"extendedAtlasS0normal_counts.RData") )

	load(paste0("data/extendedAtlasS0_clin.RData"))
	load( file = paste0(MainDir,"extendedAtlasS0normal_annot.RData") )
	load( file = paste0(MainDir,"extendedAtlasS0normal_counts.RData") )
	## Remove preinvasive
	clin = clin[!(clin$SampleType %in% c( "AAH","AIS","MIA" )),]
	annot = annot[!(annot$SampleType %in% c( "AAH","AIS","MIA" )),]
	counts = counts[,rownames(annot)]
	save(clin, file = paste0("data/extendedAtlasNormal_clin.RData"))
	save(annot, file = paste0(MainDir,"extendedAtlasNormal_annot.RData") )
	save(counts, file = paste0(MainDir,"extendedAtlasNormal_counts.RData") )

	seu = CreateSeuratObject(counts = counts, project = "extendedAtlasNormal", meta.data = annot)
	seu = SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("percent.mt"), verbose = FALSE)
	seu = RunPCA(seu)
	save(seu, file = paste0(MainDir,"extendedAtlasNormal_seu_CancerEpithelial.RData"))
	## With, instead, preinvasive
	load(paste0("data/extendedAtlasS0_clin.RData"))
	load( file = paste0(MainDir,"extendedAtlasS0normal_annot.RData") )
	load( file = paste0(MainDir,"extendedAtlasS0normal_counts.RData") )
	seu = CreateSeuratObject(counts = counts, project = "extendedAtlasNormalS0", meta.data = annot)
	seu = SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("percent.mt"), verbose = FALSE)
	seu = RunPCA(seu)
	save(seu, file = paste0(MainDir,"extendedAtlasNormal_seu_CancerEpithelialS0.RData"))
}

tps_scoring_extended_normal = function(OutDir){
	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelial.RData"))
	ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
	mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
	ribo.genes = ribo$Approved.symbol
	seu = subset(x = seu, features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))]) # 12677 x 101849
	# data = seu@assays$SCT@data
	db_rand = dan.barkley_MakeRand(seu,tps, 3)
	save(db_rand,file = paste0(OutDir,"db_rand_extendedAtlasNormal.RData"))
	scores = dan.Barkley_GeneToEnrichment_AmsCentered(seu, tps, db_rand)
	save(scores,file = paste0(OutDir,"tps_scores_extendedAtlasNormal.RData"))
}

regress_scores_normal_epi = function( OutDir,scores,prefix,ct_estimate_coefficients,tps=NULL ){
		scores_store = scores
		if (is.null(tps)){ load(file = paste0(OutDir,"../tps_discovery/tps.RData")) }
		scores_noncoo = scores[(scores$Epi_Cancer %in% ct_estimate_coefficients ),]
		scores_noncoo[,names(tps)] = apply(scores_noncoo[,names(tps)],2,function(x) x-mean(x))
		# scores_epi = scores[!(scores$Epi_Cancer %in% ct_estimate_coefficients ),]
		scc = scores
		for (tp in names(tps))
		{
		  dcat(tp)
		  scores_noncoo$tp = scores_noncoo[,tp]
		  lm = lm(tp~0+Dataset, data = scores_noncoo)
		  lmc = lm$coefficients
		  names(lmc) = gsub( "Dataset","",names(lmc) )
		  # lmc = lmc-mean(lmc)
		  print(round(sort(lmc,decreasing=T),3))
		  print(abs(max(lmc)-mean(lmc))/abs(mean(lmc)) )
		  print(abs(min(lmc)-mean(lmc))/abs(mean(lmc)) )
		  print(abs(max(lmc)-min(lmc)) )
		  partial_residuals = scc[,tp]
		  for (d in names(lmc)){
		  	partial_residuals = partial_residuals-(scc$Dataset==d)*lmc[d]
		  }
		  scc[,tp] = partial_residuals
		}
		scores = scc
		scores[,names(tps)] = apply(scores[,names(tps)],2,function(x) x-mean(x))
		save(scores, file = paste0( OutDir,prefix,"_tps_scores_DatasetRegressed.RData" ))
		# scores = scores_store
		# scores_noncoo = scores[(scores$Epi_Cancer %in% ct_estimate_coefficients ),]
		# # scores_epi = scores[!(scores$Epi_Cancer %in% ct_estimate_coefficients ),]
		# scc = scores
		# for (tp in names(tps))
		# {
		#   dcat(tp)
		#   scores_noncoo$tp = scores_noncoo[,tp]
		#   lm = lm(tp~0+Patient, data = scores_noncoo)
		#   lmc = lm$coefficients
		#   names(lmc) = gsub( "Patient","",names(lmc) )
		#   # lmc = lmc-mean(lmc)
		#   partial_residuals = scc[,tp]
		#   for (d in names(lmc) ){
		#   	partial_residuals = partial_residuals-(scc$Patient==d)*lmc[d]
		#   }
		#   scc[,tp] = partial_residuals
		# }
		# scores = scc
		# scores[,names(tps)] = apply(scores[,names(tps)],2,function(x) x-mean(x))
		# save(scores, file = paste0( OutDir,prefix,"_tps_scores_PatientRegressed.RData" ))
}

tps_MECO = function( OutDir, scores, prefix ){
	library(ComplexHeatmap)
	load(file = paste0(OutDir,"../../tps_discovery/tps.RData"))
	scores = scores[scores$Epi_Cancer=="Cancer",]
	scores_store = scores
	cordf = cor(scores[,names(tps)],method = "pearson")
	order_tps = c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","Stress_AP1","Stress_HSP","Stress_secreted","Unas_emp","RNA_processing","MHC-I","Interferon","OxPhos","Cell_proliferation","Basal-like","pEMT","EMT","Metal","Hypoxia","Translation_initiation","Unfolded_protein_response" )
	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_pearson.pdf"),7,7)
	corrplot(cordf[order_tps,order_tps],col=colorz_solid,tl.col="black",tl.srt = 45)#,tl.cex=1/.pt,cl.cex=1/.pt,pch.cex=3/.pt)
	dev.off()
	library(colorRamp2)
	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_pearson_ch.pdf"),5,5)
	plot=(ComplexHeatmap::Heatmap(cordf[order_tps,order_tps],rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_pearson_ch_dendrogram.pdf"),10,10)
	a=heatmap.2(cordf,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("dodgerblue4","white","firebrick3"))(200))
	dev.off()

	selected_cs = names(tps)
	patientStates = matrix(0,nrow = length(unique(scores$Patient)), ncol = length(selected_cs),dimnames = list(unique(scores$Patient),selected_cs))
	for (p in rownames(patientStates)){
		for (cs in selected_cs){
			patientStates[p,cs] = mean(scores[scores$Patient==p,cs])
		}
	}
	patient_level = matrix(nrow=length(selected_cs),ncol=length(selected_cs),dimnames = list(selected_cs,selected_cs))
	for (cs1 in names(tps)){
		for (cs2 in names(tps)){
			patient_level[cs1,cs2] = cor(patientStates[,cs1],patientStates[,cs2],method = "pearson")
		}  
	}
	# diag(patient_level) = 0
	pdf(paste0(OutDir,prefix,"_csMECO_Patients_pearson.pdf"),10,10)
	a=heatmap.2(patient_level,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("navy","white","red4"))(200))
	dev.off()

	patient_level = patient_level[order_tps,order_tps]
	cordf = cordf[order_tps,order_tps]
	pdf(paste0(OutDir,prefix,"_csMECO_Patients_pearson_ch.pdf"),5,5)
	plot=(ComplexHeatmap::Heatmap(patient_level[order_tps,order_tps],rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(patient_level), "cm")/2.5,height=unit(nrow(patient_level), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("midnightblue","white","chocolate3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, patient-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	pdf(paste0(OutDir,prefix,"_csMECO_Patients_pearson_ch_dendrogram.pdf"),10,10)
	a=heatmap.2(patient_level,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("dodgerblue4","white","firebrick3"))(200))
	dev.off()

	### single cell correlation averaged across patients (weighted for number of cells)
	cordf = dan.df(order_tps,order_tps)
	patient_cordf_list = list()
	patient_ncells = c()
	for (p in unique(scores$Patient)){
		ts = scores[scores$Patient==p,]
		patient_ncells = c(patient_ncells,nrow(ts))
		patient_cordf_list[[p]] = cor(ts[,names(tps)],method = "pearson")
	}
	names(patient_ncells) = unique(scores$Patient)
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			values = sapply(names(patient_cordf_list),FUN=function(x) patient_cordf_list[[x]][rn,cn])
			weights = patient_ncells[names(patient_cordf_list)]
			cordf[rn,cn] = weighted.mean(values,weights)
		}
	}
	cordf = data.matrix(cordf)
	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_PatientWiseAveraged_pearson_ch.pdf"),5,5)
	plot=(ComplexHeatmap::Heatmap(cordf[order_tps,order_tps],rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R across single cells,\npatient-wise and averaged",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_PatientWiseAveraged_pearson_ch_dendrogram.pdf"),10,10)
	a=heatmap.2(cordf,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("dodgerblue4","white","firebrick3"))(200))
	dev.off()


	these_scores = scores[(scores$Epi_Cancer=="Cancer") & ( scores$SampleType=="Primary" ),]
	cordf = cor(these_scores[,names(tps)],method = "pearson")
	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_pearson_PrimaryOnly.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()
	selected_cs = names(tps)
	patientStates = matrix(0,nrow = length(unique(these_scores$Patient)), ncol = length(selected_cs),dimnames = list(unique(these_scores$Patient),selected_cs))
	for (p in rownames(patientStates)){
		for (cs in selected_cs){
			patientStates[p,cs] = mean(these_scores[these_scores$Patient==p,cs])
		}
	}
	patient_level = matrix(nrow=length(selected_cs),ncol=length(selected_cs),dimnames = list(selected_cs,selected_cs))
	for (cs1 in names(tps)){
		for (cs2 in names(tps)){
			patient_level[cs1,cs2] = cor(patientStates[,cs1],patientStates[,cs2],method = "pearson")
		}  
	}
	diag(patient_level) = 0
	pdf(paste0(OutDir,prefix,"_csMECO_Patients_pearson_PrimaryOnly.pdf"),10,10)
	a=heatmap.2(patient_level,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("navy","white","red4"))(200))
	dev.off()

	these_scores = scores[(scores$Epi_Cancer=="Cancer") & ( scores$SampleType=="Metastasis" ),]
	cordf = cor(these_scores[,names(tps)],method = "pearson")
	pdf(paste0(OutDir,prefix,"_csMECO_SingleCells_pearson_MetastasisOnly.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()
	selected_cs = names(tps)
	patientStates = matrix(0,nrow = length(unique(these_scores$Patient)), ncol = length(selected_cs),dimnames = list(unique(these_scores$Patient),selected_cs))
	for (p in rownames(patientStates)){
		for (cs in selected_cs){
			patientStates[p,cs] = mean(these_scores[these_scores$Patient==p,cs])
		}
	}
	patient_level = matrix(nrow=length(selected_cs),ncol=length(selected_cs),dimnames = list(selected_cs,selected_cs))
	for (cs1 in names(tps)){
		for (cs2 in names(tps)){
			patient_level[cs1,cs2] = cor(patientStates[,cs1],patientStates[,cs2],method = "pearson")
		}  
	}
	diag(patient_level) = 0
	pdf(paste0(OutDir,prefix,"_csMECO_Patients_pearson_MetastasisOnly.pdf"),10,10)
	a=heatmap.2(patient_level,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("navy","white","red4"))(200))
	dev.off()
	dir.create(paste0(OutDir,prefix,"_dataset_wise/"))
	for (dataset in unique(scores$Dataset)){
		these_scores = scores[( scores$Dataset==dataset ),]
		cordf = cor(these_scores[,names(tps)],method = "pearson")
		pdf(paste0(OutDir,prefix,"_dataset_wise/csMECO_SingleCells_pearson_",dataset,".pdf"),7,7)
		corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
		dev.off()
		selected_cs = names(tps)
		patientStates = matrix(0,nrow = length(unique(these_scores$Patient)), ncol = length(selected_cs),dimnames = list(unique(these_scores$Patient),selected_cs))
		for (p in rownames(patientStates)){
			for (cs in selected_cs){
				patientStates[p,cs] = mean(these_scores[these_scores$Patient==p,cs])
			}
		}
		patient_level = matrix(nrow=length(selected_cs),ncol=length(selected_cs),dimnames = list(selected_cs,selected_cs))
		for (cs1 in names(tps)){
			for (cs2 in names(tps)){
				patient_level[cs1,cs2] = cor(patientStates[,cs1],patientStates[,cs2],method = "pearson")
			}  
		}
		diag(patient_level) = 0
		pdf(paste0(OutDir,prefix,"_dataset_wise/csMECO_Patients_pearson_",dataset,".pdf"),10,10)
		a=heatmap.2(patient_level,tracecol = NA, margins=c(10,10),col = colorRampPalette(c("navy","white","red4"))(200))
		dev.off()
	}
}

tps_across_ClinicalCharacteristics_singlecells = function(OutDir, scores, order_tps = NULL){

	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	if (!is.null(order_tps)){ tps = tps[order_tps] }

	# tps across stages in LUADs only
	scores_store = scores
	scores = scores[scores$Epi_Cancer=="Cancer",]
	scores = scores[!is.na(scores$Stage_collapsed),]
	xLevels = c( "Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Metastasis" )
	xColors = c( "gold","orange","tomato3","firebrick","brown" )
	scores$xlevel = NA
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$SampleType=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$SampleType=="Primary","Stage_collapsed"])
	pdf( paste0(OutDir,"singlecells_Tps_vs_SampleTypesStages_OnlyLUADs.pdf"),8,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"singlecells_Tps_vs_SampleTypesStages_OnlyLUADs_SplitDataset.pdf"),22,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   fill = scores$Dataset
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, fill, xlab = "", filllab="Dataset", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"singlecells_Tps_vs_SampleTypesStages_OnlyLUADs_ViolinPlots.pdf"),8,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.violinplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
	# dataset-wise
	dir.create(paste0(OutDir,"dataset_wise/"))
	for (dataset in unique(scores$Dataset)){
		tsc = scores[scores$Dataset==dataset,]
		pdf( paste0(OutDir,"dataset_wise/",dataset,"_singlecells_Tps_vs_SampleTypesStages_OnlyLUADs.pdf"),8,6 )
		for (tp in names(tps)){
		   x = tsc[,"xlevel"]
		   txColors = xColors[xLevels %in% unique(x)]
		   txLevels = xLevels[xLevels %in% unique(x)]
		   y = tsc[,tp]
		   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		   plot=dan.boxplots.multipages( factor(x,levels=txLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = txColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
		   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
		   print(plot)
		}
		dev.off()
	}
	
	scores = scores_store
	scores = scores[!is.na(scores$Stage_collapsed),]
	xLevels = c( "basal","AT1","ciliated","club","preTB", "AT2","AT0" ,"AAH","AIS","MIA","Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Metastasis" )
	xColors = c( "gray10","gray26","gray42","gray58","gray74","gray90","black","dodgerblue4","steelblue3","steelblue1","gold","orange","tomato3","firebrick","brown" )
	scores$xlevel = NA
	scores[scores$Epi_Cancer!="Cancer","xlevel"] = scores[scores$Epi_Cancer!="Cancer","Epi_Cancer"]
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$xlevel=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$xlevel=="Primary","Stage_collapsed"])
	xLevels_keep = xLevels %in% unique(scores$xlevel)
	xLevels = xLevels[xLevels_keep]
	xColors = xColors[xLevels_keep]
	pdf( paste0(OutDir,"singlecells_Tps_vs_CellSampleTypes_full.pdf"),10,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	pdf( paste0(OutDir,"singlecells_Tps_vs_CellSampleTypes_full_ViolinPlots.pdf"),10,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.violinplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	scores = scores_store
	scores = scores[!is.na(scores$Stage_collapsed),]
	xLevels = c( "AT1","AT0","AT2" ,"Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Mets, LN","Mets, distal" )
	xColors = c( "black","gray44","gray77","gold","orange","tomato3","firebrick","brown","sienna4" )
	scores$xlevel = NA
	scores[(scores$SampleSubType=="LN") %in% c(T),"SampleType"] = "Mets, LN"
	scores[(scores$SampleType=="Metastasis") %in% c(T),"SampleType"] = "Mets, distal"
	scores[scores$Epi_Cancer!="Cancer","xlevel"] = scores[scores$Epi_Cancer!="Cancer","Epi_Cancer"]
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$xlevel=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$xlevel=="Primary","Stage_collapsed"])
	xLevels_keep = xLevels %in% unique(scores$xlevel)
	xLevels = xLevels[xLevels_keep]
	xColors = xColors[xLevels_keep]
	scores = scores[scores$xlevel %in% xLevels,]
	pdf( paste0(OutDir,"singlecells_Tps_vs_CellSampleTypes_full_SplitMets.pdf"),4,3 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	pdf( paste0(OutDir,"singlecells_Tps_vs_CellSampleTypes_full_ViolinPlots_SplitMets.pdf"),4,3 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.violinplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	## patterns in single cells
	patients_patterns = unique(scores[(scores$histologic_pattern %in% c( "lepidic","papillary","acinar","micropapillary","solid" )) %in% c(T),"Patient"]) # excluding sarcomatoid (only 21 cells)
	ts = scores[scores$Patient %in% patients_patterns,]
	ts = ts[ts$Epi_Cancer %in% c( "AT2","club","preTB","AT0","Cancer" ),]
	ts$xlevel = NA
	ts[ts$Epi_Cancer!="Cancer","xlevel"] = ts[ts$Epi_Cancer!="Cancer","Epi_Cancer"]
	ts[ts$Epi_Cancer=="Cancer","xlevel"] = ts[ts$Epi_Cancer=="Cancer","histologic_pattern"]
	dtable(ts$xlevel)
	xLevels = c( "club","preTB","AT2","AT0", "lepidic" ,"papillary","acinar","solid","micropapillary")
	xColors = c( "gray58","gray74","gray90","black","dodgerblue4","lightseagreen","orange","red","brown" )
	pdf( paste0(OutDir,"singlecells_Tps_vs_HistologicPatterns_Bischoff.pdf"),10,6 )
	for (tp in names(tps)){
		dcat(tp)
	   x = ts[,"xlevel"]
	   y = ts[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	## smoking in single cells
	scores = scores_store
	ts = scores[(scores$Epi_Cancer=="Cancer") & (!is.na(scores$Smoking)),]
	ts$xlevel = ifelse(ts$Smoking=="Never","Never smoker","Ever smoker")
	dtable(ts$xlevel)
	xLevels = c("Never smoker","Ever smoker")
	xColors = c("forestgreen","gray22")
	pdf( paste0(OutDir,"singlecells_Tps_vs_Smoking.pdf"),6,6 )
	for (tp in names(tps)){
	   x = ts[,"xlevel"]
	   y = ts[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
	## smoking in single cells detailed
	scores = scores_store
	ts = scores[(scores$Epi_Cancer=="Cancer") & (!is.na(scores$Smoking)),]
	ts = ts[(ts$Smoking!="Ever"),]
	ts$xlevel = ts$Smoking
	xLevels = c("Never","Ex","Current")
	xColors = c("forestgreen","gray55","gray22")
	pdf( paste0(OutDir,"singlecells_Tps_vs_Smoking_detailed.pdf"),6,6 )
	for (tp in names(tps)){
	   x = ts[,"xlevel"]
	   y = ts[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	scores = scores_store
	scores = scores[!is.na(scores$Stage_collapsed),]
	xLevels = c( "AT2","Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Metastasis" )
	xColors = c( "gray58","gold","orange","tomato3","firebrick","brown" )
	scores$xlevel = NA
	scores[scores$Epi_Cancer!="Cancer","xlevel"] = scores[scores$Epi_Cancer!="Cancer","Epi_Cancer"]
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$xlevel=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$xlevel=="Primary","Stage_collapsed"])
	scores = scores[scores$xlevel %in% xLevels,]
	xLevels_keep = xLevels %in% unique(scores$xlevel)
	xLevels = xLevels[xLevels_keep]
	xColors = xColors[xLevels_keep]
	pdf( paste0(OutDir,"singlecells_Tps_vs_CellSampleTypes_full_ViolinPlots_withAT2.pdf"),8,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.violinplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
}

tps_across_ClinicalCharacteristics_singlecells_patientlevel = function(OutDir, scores, order_tps = NULL){

	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	if (!is.null(order_tps)){ tps = tps[order_tps] }
	scores[scores$SampleType=="Normal","Stage_collapsed"] = NA
	scores$PST_ec = paste0(scores$Patient,"_",scores$SampleType,"_",scores$Epi_Cancer)
	pst = dan.df(unique(scores$PST_ec),c( "Dataset","Patient","PST_ec","SampleType","SampleSubType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern", order_tps ))	
	for (this_pst in rownames(pst)){
		for (tp in order_tps){
			pst[this_pst,tp] = mean(scores[ scores$PST_ec==this_pst,tp ])
		}
		for (var in c( "Dataset","Patient","PST_ec","SampleType","SampleSubType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern" )){
			pst[this_pst,var] = unique(scores[ scores$PST_ec==this_pst,var ])[1]
		}
	}
	scores = pst
	# tps across stages in LUADs only
	scores_store = scores
	scores = scores[scores$Epi_Cancer=="Cancer",]
	scores = scores[(scores$SampleType=="Normal") | (!is.na(scores$Stage_collapsed)),]
	xLevels = c( "Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Metastasis" )
	xColors = c( "gold","orange","tomato3","firebrick","brown" )
	scores$xlevel = NA
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$SampleType=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$SampleType=="Primary","Stage_collapsed"])
	scores$xcolorz = dan.expand_colors(scores$xlevel,xLevels,xColors)
	pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_SampleTypesStages_OnlyLUADs.pdf"),8,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = scores$xcolorz, labelJitteredPoints = NULL, includeJitters = T )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_SampleTypesStages_OnlyLUADs_SplitDataset.pdf"),22,6 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   fill = scores$Dataset
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, fill, xlab = "", filllab="Dataset", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	# dataset-wise
	dir.create(paste0(OutDir,"dataset_wise/"))
	for (dataset in unique(scores$Dataset)){
		tsc = scores[scores$Dataset==dataset,]
		pdf( paste0(OutDir,"dataset_wise/",dataset,"_patientlevel_singlecells_Tps_vs_SampleTypesStages_OnlyLUADs.pdf"),8,6 )
		for (tp in names(tps)){
		   x = tsc[,"xlevel"]
		   txColors = xColors[xLevels %in% unique(x)]
		   txLevels = xLevels[xLevels %in% unique(x)]
		   y = tsc[,tp]
		   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		   plot=dan.boxplots.multipages( factor(x,levels=txLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = txColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T )
		   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
		   print(plot)
		}
		dev.off()
	}
	
	scores = scores_store
	scores = scores[(scores$SampleType=="Normal") | (!is.na(scores$Stage_collapsed)),]
	xLevels = c( "AT1","AT0","AT2" ,"Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Metastasis" )
	xColors = c( "black","gray44","gray77","gold","orange","tomato3","firebrick","brown" )
	scores$xlevel = NA
	scores[scores$Epi_Cancer!="Cancer","xlevel"] = scores[scores$Epi_Cancer!="Cancer","Epi_Cancer"]
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$xlevel=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$xlevel=="Primary","Stage_collapsed"])
	scores = scores[scores$xlevel %in% xLevels,]
	xLevels_keep = xLevels %in% unique(scores$xlevel)
	xLevels = xLevels[xLevels_keep]
	xColors = xColors[xLevels_keep]
	scores$xcolorz = dan.expand_colors(scores$xlevel,xLevels,xColors)
	pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_CellSampleTypes_full.pdf"),3,2 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = scores$xcolorz, labelJitteredPoints = NULL, includeJitters = T )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	scores = scores_store
	scores = scores[(scores$SampleType=="Normal") | (!is.na(scores$Stage_collapsed)),]
	xLevels = c( "AT1","AT0","AT2" ,"Primary, stage I","Primary, stage II","Primary, stage III","Primary, stage IV","Mets, LN","Mets, distal" )
	xColors = c( "black","gray44","gray77","gold","orange","tomato3","firebrick","brown","sienna4" )
	scores$xlevel = NA
	scores[(scores$SampleSubType=="LN") %in% c(T),"SampleType"] = "Mets, LN"
	scores[(scores$SampleType=="Metastasis") %in% c(T),"SampleType"] = "Mets, distal"
	scores[scores$Epi_Cancer!="Cancer","xlevel"] = scores[scores$Epi_Cancer!="Cancer","Epi_Cancer"]
	scores[scores$Epi_Cancer=="Cancer","xlevel"] = scores[scores$Epi_Cancer=="Cancer","SampleType"]
	scores[scores$xlevel=="Primary","xlevel"] = paste0("Primary, stage ",scores[scores$xlevel=="Primary","Stage_collapsed"])
	scores = scores[scores$xlevel %in% xLevels,]
	xLevels_keep = xLevels %in% unique(scores$xlevel)
	xLevels = xLevels[xLevels_keep]
	xColors = xColors[xLevels_keep]
	scores$xcolorz = dan.expand_colors(scores$xlevel,xLevels,xColors)
	pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_CellSampleTypes_full_SplitMets.pdf"),3,2 )
	for (tp in names(tps)){
	   x = scores[,"xlevel"]
	   y = scores[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = scores$xcolorz, labelJitteredPoints = NULL, includeJitters = T )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()

	# ## patterns in single cells
	# patients_patterns = unique(scores[(scores$histologic_pattern %in% c( "lepidic","papillary","acinar","micropapillary","solid" )) %in% c(T),"Patient"]) # excluding sarcomatoid (only 21 cells)
	# ts = scores[scores$Patient %in% patients_patterns,]
	# ts = ts[ts$Epi_Cancer %in% c( "AT2","club","preTB","AT0","Cancer" ),]
	# ts$xlevel = NA
	# ts[ts$Epi_Cancer!="Cancer","xlevel"] = ts[ts$Epi_Cancer!="Cancer","Epi_Cancer"]
	# ts[ts$Epi_Cancer=="Cancer","xlevel"] = ts[ts$Epi_Cancer=="Cancer","histologic_pattern"]
	# dtable(ts$xlevel)
	# xLevels = c( "club","preTB","AT2","AT0", "lepidic" ,"papillary","acinar","solid","micropapillary")
	# xColors = c( "gray58","gray74","gray90","black","dodgerblue4","lightseagreen","orange","red","brown" )
	# ts$xcolorz = dan.expand_colors(ts$xlevel,xLevels,xColors)
	# pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_HistologicPatterns_Bischoff.pdf"),10,6 )
	# for (tp in names(tps)){
	# 	dcat(tp)
	#    x = ts[,"xlevel"]
	#    y = ts[,tp]
	#    ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	#    ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	#    plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = ts$xcolorz, labelJitteredPoints = NULL, includeJitters = T )
	#    # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	#    print(plot)
	# }
	# dev.off()

	## smoking in single cells
	scores = scores_store
	ts = scores[(scores$Epi_Cancer=="Cancer") & (!is.na(scores$Smoking)),]
	ts$xlevel = ifelse(ts$Smoking=="Never","Never smoker","Ever smoker")
	dtable(ts$xlevel)
	xLevels = c("Never smoker","Ever smoker")
	xColors = c("forestgreen","gray22")
	ts$xcolorz = dan.expand_colors(ts$xlevel,xLevels,xColors)
	pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_Smoking.pdf"),2,2 )
	for (tp in names(tps)){
	   x = ts[,"xlevel"]
	   y = ts[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox", comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = ts$xcolorz, labelJitteredPoints = NULL, includeJitters = T )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
	## smoking in single cells
	scores = scores_store
	ts = scores[(scores$Epi_Cancer=="Cancer") & (!is.na(scores$Smoking)),]
	ts = ts[(ts$Smoking!="Ever"),]
	ts$xlevel = ts$Smoking
	xLevels = c("Never","Ex","Current")
	xColors = c("forestgreen","gray55","gray22")
	ts$xcolorz = dan.expand_colors(ts$xlevel,xLevels,xColors)
	pdf( paste0(OutDir,"patientlevel_singlecells_Tps_vs_Smoking_detailed.pdf"),2.1,1.7 )
	for (tp in names(tps)){
	   x = ts[,"xlevel"]
	   y = ts[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = ts$xcolorz, labelJitteredPoints = NULL, includeJitters = T )
	   # if (tp %in% c("Cell_proliferation","Metal","Basal-like","Translation_initiation","RNA_processing","Unfolded_protein_response") ){ plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = max(c(min(y),-2)), ylimRight = 1.2,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F ) }
	   print(plot)
	}
	dev.off()
}

tps_across_ClinicalCharacteristics_bulk_datasets = function(OutDir, order_tps = NULL, scoring_method='singscore'){
	library(singscore)
	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	if (!is.null(order_tps)){ tps = tps[order_tps] }

	## patterns in lumu
	library(singscore)
	load( file = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/ExprTables/lung_multiregion_log2TPMcombat_table_GeneSymbols.RData" )
	logtpm = ge[intersect(rownames(ge),tps_universe),]
	Clin = read.table(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Clinical_4assigned.txt"), sep = "\t", header = TRUE, quote = '', stringsAsFactors = FALSE)
	lepidic_color = "dodgerblue4"
	acinar_color = "orange"
	papillary_color = "lightseagreen"
	solid_color = "red"
	Clin$color = lepidic_color # c("#FF8000", "#6666FF", "#009900")
	Clin[Clin$Pattern=="acinar","color"] = acinar_color
	Clin[Clin$Pattern=="papillary","color"] = papillary_color
	Clin[Clin$Pattern=="solid","color"] = solid_color
	Clin = Clin[Clin$Region!="A",] # assigns it to parent scope
	rownames(Clin) = Clin$Sample
	rankData = rankGenes(logtpm)
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_CHUV.pdf"),3.5,2.5 )
	for (tp in names(tps)){
	   x = Clin[,"Pattern"]
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[tp]]),verbose=F,method='ssgsea'))) }
	   Clin[rownames(scoredf),tp] = scoredf$TotalScore
	   y = Clin[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=c("lepidic","papillary","acinar","solid")), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = c(lepidic_color,papillary_color,acinar_color,solid_color), jitterColors = Clin$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	save(Clin,file=paste0(OutDir,"Bulk_Tps_scores_CHUV.RData"))

	## patterns, stages, has_metastasis, age, smoking, survival in TCGA
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_TPM/TCGA_LUAD_ge_TPM_GeneSymbols.RData")
	cn = as.character(colnames(ge))
	cn = cn[substr(cn,14,15) %in% c("01","11")] # only primary
	ge = ge[intersect(rownames(ge),tps_universe),cn]
	cn_short = substr(cn,1,15)
	colnames(ge) = cn_short
	logtpm = log2(ge+1)
	logtpm_tcga = logtpm
	CommonDataDir = "/mnt/ed2/daniele/Common_Data/"
	Clinical_extended = as.data.frame(t(read.table(paste0(CommonDataDir,"Clinical_extended/LUAD.Clinical.txt"), quote = '', sep = "\t", header = FALSE, row.names = 1)), stringsAsFactors = FALSE)
	Clinical_extended[,"Patient"] <- toupper(as.character(Clinical_extended[,"bcr_patient_barcode"]))
	rownames(Clinical_extended) = Clinical_extended[,"Patient"]
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/TCGA_histology_formatted_permissive_withNormals.RData"))
	rownames(Clin) = Clin$Sample
	Clin_patterns = Clin
	Clin_patterns$Pattern = as.character(Clin_patterns$Pattern)
	Clin = dan.df(colnames(ge),c("Sample","Patient","Pattern","color","pathologic_stage","pathologic_m","pathologic_n","age_at_initial_pathologic_diagnosis","gender","tobacco_smoking_history","vital_status", "days_to_death", "days_to_last_followup"))
	Clin$Sample = rownames(Clin)
	Clin$Patient = substr(Clin$Sample,1,12)
	Clin[intersect(rownames(Clin),rownames(Clin_patterns)),"Pattern"] = Clin_patterns[intersect(rownames(Clin),rownames(Clin_patterns)),"Pattern"]
	Clin[intersect(rownames(Clin),rownames(Clin_patterns)),"color"] = Clin_patterns[intersect(rownames(Clin),rownames(Clin_patterns)),"color"]
	for (pat in rownames(Clinical_extended)){ 
		for (var in c("Patient","pathologic_stage","pathologic_m","pathologic_n","age_at_initial_pathologic_diagnosis","gender","tobacco_smoking_history","vital_status", "days_to_death", "days_to_last_followup"))
		Clin[Clin$Patient==pat,var] = Clinical_extended[pat,var]
	}
	Clin$stage_collapsed = NA
	Clin[Clin$pathologic_stage %in% c("stage i","stage ia", "stage ib"),"stage_collapsed"] = "I"
	Clin[Clin$pathologic_stage %in% c("stage iia", "stage iib"),"stage_collapsed"] = "II"
	Clin[Clin$pathologic_stage %in% c("stage iia", "stage iib"),"stage_collapsed"] = "II"
	Clin[Clin$pathologic_stage %in% c("stage iiia", "stage iiib"),"stage_collapsed"] = "III"
	Clin[Clin$pathologic_stage %in% c("stage iv"),"stage_collapsed"] = "IV"
	Clin = Clin[Clin$Sample %in% colnames(logtpm),]
	Clin$Sample = as.character(Clin$Sample)
	rownames(Clin) = Clin$Sample
	rankData = rankGenes(logtpm)
	for (tp in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[tp]]),verbose=F,method='ssgsea'))) }
	   Clin[rownames(scoredf),tp] = scoredf$TotalScore
	}
	save(Clin,file=paste0(OutDir,"Bulk_Tps_scores_TCGA.RData"))
	load(file=paste0(OutDir,"Bulk_Tps_scores_TCGA.RData"))
	# vs patterns
	tc = Clin[!is.na(Clin$Pattern),]
	tc$color
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_TCGA.pdf"),2.1,1.7 )
	for (tp in names(tps)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=c("normal","lepidic","papillary","acinar","solid")), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y), xColors = c("gray",lepidic_color,papillary_color,acinar_color,solid_color), jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_TCGA_byStage.pdf"),10,4 )
	# for (tp in names(tps)){
	#    x = tc[,"Pattern"]
	#    y = tc[,tp]
	#    plot=dan.boxplots.multipages( tc$stage_collapsed, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	#    print(plot)
	# }
	# dev.off()
	# vs stages
	tc = Clin
	tc[(tc$Pattern=="normal") %in% c(T),"stage_collapsed"] = "Normal"
	tc = tc[!is.na(tc$stage_collapsed),]
	xLevels = c("Normal","I","II","III","IV")
	xColors = c("gray","gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$stage_collapsed,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Stages_TCGA.pdf"),7,5 )
	for (tp in names(tps)){
	   x = tc[,"stage_collapsed"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "Stage", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs has_metastasis
	tc = Clin[((!is.na(Clin$pathologic_m)) & (!is.na(Clin$pathologic_n))) & (!((Clin$pathologic_m=="mx") & (Clin$pathologic_n=="nx"))),]
	tc = tc[substr(tc$Sample,14,15)=="01",]
	tc$has_metastasis = "Metastatic_LNonly"
	tc[(tc$pathologic_m %in% c("m0","mx")) & (tc$pathologic_n %in% c("n0","nx")),"has_metastasis"] = "Localized"
	tc[(tc$pathologic_m %in% c("m1","m1a","m1b")),"has_metastasis"] = "Metastatic_distal"
	xLevels = c("Localized","Metastatic_LNonly","Metastatic_distal")
	xColors = c("gold","tomato3","brown")
	tc$color = dan.expand_colors(tc$has_metastasis,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_HasMetastasis_TCGA.pdf"),2.5,2.5 )
	for (tp in names(tps)){
	   x = tc[,"has_metastasis"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal.test", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking
	tc = Clin[!is.na(Clin$tobacco_smoking_history),]
	tc = tc[substr(tc$Sample,14,15)=="01",]
	tc$smoking = ifelse(tc$tobacco_smoking_history==1,"Never smoker","Ever smoker")
	xLevels = c("Never smoker","Ever smoker")
	xColors = c("forestgreen","gray22")
	tc$color = dan.expand_colors(tc$smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_TCGA.pdf"),2.5,2.5 )
	for (tp in names(tps)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_TCGA_byStage.pdf"),5,2 )
	for (tp in names(tps)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( tc$stage_collapsed, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking detailed
	tc = Clin[!is.na(Clin$tobacco_smoking_history),]
	tc = tc[substr(tc$Sample,14,15)=="01",]
	tc = tc[tc$tobacco_smoking_history!=5,]
	tc$smoking = "Never"
	tc[tc$tobacco_smoking_history==2,"smoking"] = "Current"
	tc[tc$tobacco_smoking_history==3,"smoking"] = "Ex (>15y)"
	tc[tc$tobacco_smoking_history==4,"smoking"] = "Ex (<=15y)"
	xLevels = c("Never","Ex (>15y)", "Ex (<=15y)", "Current")
	xColors = c("forestgreen","gray77","gray50","gray22")
	tc$color = dan.expand_colors(tc$smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_SmokingDetailed_TCGA.pdf"),2.1,1.7 )
	for (tp in names(tps)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"Bulk_Tps_vs_SmokingDetailed_TCGA_byStage.pdf"),5,2 )
	for (tp in names(tps)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( tc$stage_collapsed, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs survival
	library(survival)
	this_OutDir = paste0(OutDir,"survival_tcga/")
	dir.create(this_OutDir)
	colorz = c("dodgerblue4","firebrick")
	legendz_half = c("Below mean","Above mean")
	legendz_extremes = c("Bottom 25%","Top 25%")
	this_Clin = Clin[substr(Clin$Sample,14,15)=="01",]
	rownames(this_Clin) = this_Clin$Patient
	a = as.numeric(as.character(this_Clin$days_to_death))
	b = as.numeric(as.character(this_Clin$days_to_last_followup))
	vital_status_num <- vector(mode="numeric", length=length(this_Clin$vital_status))
	times <- vector(mode="numeric", length=length(vital_status_num))
	for (v in 1:length(vital_status_num))
	{
	  if (this_Clin$vital_status[v]=="alive")
	  {
	     vital_status_num[v] <- 0
	     times[v] <- b[v]
	  }
	  else
	  {
	     vital_status_num[v] <- 1
	     times[v] <- a[v]
	  }
	}
	this_Clin$Times <- times
	this_Clin$vital_status_num <- vital_status_num
	this_Clin$stage_numeric = NA
	this_Clin[(this_Clin$stage_collapsed == "I") %in% c(T),"stage_numeric"] = 1
	this_Clin[(this_Clin$stage_collapsed == "II") %in% c(T),"stage_numeric"] = 2
	this_Clin[(this_Clin$stage_collapsed == "III") %in% c(T),"stage_numeric"] = 3
	this_Clin[(this_Clin$stage_collapsed == "IV") %in% c(T),"stage_numeric"] = 4
	for (tp in names(tps))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_SplittingHalf.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxUniVariate.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxMultiVariate.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_Extremes.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}

	## patterns, stages, age, smoking, survival in Chen
	load( file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_normalized_GeneSymbols.RData") )
	chen_indir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_onChen/"
	load(file = paste0(chen_indir,"ClinChen.RData"))
	logtpm = ge[intersect(rownames(ge),tps_universe),intersect(rownames(ClinChen),colnames(ge))]
	logtpm = logtpm[rowSums(is.na(logtpm))==0, ]
	logtpm_chen = logtpm
	Clin = ClinChen[intersect(rownames(ClinChen),colnames(ge)),]
	rankData = rankGenes(logtpm)
	for (tp in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   Clin[rownames(scoredf),tp] = scoredf$TotalScore
	}
	save(Clin,file=paste0(OutDir,"Bulk_Tps_scores_ChenEAS.RData"))
	# vs patterns
	tc = Clin[!is.na(Clin$Pattern),]
	tc = tc[tc$Pattern!="unknown",]
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_ChenEAS.pdf"),7,5 )
	for (tp in names(tps)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=c("lepidic","papillary","acinar","solid","micropapillary")), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = c(lepidic_color,papillary_color,acinar_color,solid_color,"purple"), jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_ChenEAS_byStage.pdf"),10,4 )
	for (tp in names(tps)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( tc$Stage, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs stages
	tc = Clin	
	tc = tc[!is.na(tc$Stage),]
	xLevels = c("I","II","III","IV")
	xColors = c("gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$Stage,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Stages_ChenEAS.pdf"),3.5,2.5 )
	for (tp in names(tps)){
	   x = tc[,"Stage"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "Stage", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking
	tc = Clin[!is.na(Clin$Smoker),]
	tc$smoking = ifelse(tc$Smoker=="No","Never smoker","Ever smoker")
	xLevels = c("Never smoker","Ever smoker")
	xColors = c("forestgreen","gray22")
	tc$color = dan.expand_colors(tc$smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_ChenEAS.pdf"),2.5,2.5 )
	for (tp in names(tps)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_ChenEAS_byStage.pdf"),5,2 )
	for (tp in names(tps)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( tc$Stage, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs survival
	this_OutDir = paste0(OutDir,"survival_chen/")
	dir.create(this_OutDir)
	colorz = c("dodgerblue4","firebrick")
	legendz_half = c("Below mean","Above mean")
	legendz_extremes = c("Bottom 25%","Top 25%")
	this_Clin = Clin
	this_Clin$Times = as.numeric(this_Clin$OS.Month)
	this_Clin$vital_status_num = NA
	this_Clin[this_Clin$OS.Status=="Dead","vital_status_num"] = 1
	this_Clin[this_Clin$OS.Status=="Alive","vital_status_num"] = 0
	for (tp in names(tps))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$Gender, age = as.numeric(this_Clin$Age), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_SplittingHalf.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxUniVariate.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxMultiVariate.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_Extremes.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}

	## patterns, stages, has_metastasis, age, smoking, survival in TRACERx
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
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	logtpm = log2(ge[intersect(rownames(ge),tps_universe),]+0.0001)
	logtpm = logtpm[,substr(colnames(logtpm),1,8) %in% rownames(Clin)]
	# logtpm = logtpm[,!grepl("_N",colnames(logtpm))]
	library(singscore)
	rankData = rankGenes(logtpm)
	clin = dan.df(colnames(logtpm),paste0(names(tps) ) )
	colnames(clin) = gsub( "\\.","-",colnames(clin) )
	for (n in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[n]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[n]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   clin[rownames(scoredf),paste0(n )] = scoredf$TotalScore
	}
	clin$Patient = substr(rownames(clin),1,8)
	clin$Sample = rownames(clin)
	for (rn in rownames(clin)){
		tpat = clin[rn,"Patient"]
		clin[rn,"Pattern"] = Clin[tpat,"LUAD_pred_subtype"]
		clin[rn,"Smoking"] = Clin[tpat,"smoking_status_merged"]
		clin[rn,"stage_collapsed"] = gsub("B","",gsub("A","",Clin[tpat,"pathologyTNM"]))
	}
	clin[grepl("_N",rownames(clin)),"stage_collapsed"] = "Normal"
	clin[grepl("_N",rownames(clin)),"Pattern"] = "normal"
	clin$Smoking = as.character(clin$Smoking)
	clin[clin$Smoking=="Never Smoked","Smoking"] = "Never"
	clin[clin$Smoking=="Ex-Smoker","Smoking"] = "Ex"
	clin[clin$Smoking=="Smoker","Smoking"] = "Current"

	save(clin,file=paste0(OutDir,"Bulk_Tps_scores_TRACERx.RData"))
	load(file=paste0(OutDir,"Bulk_Tps_scores_TRACERx.RData"))
	# vs patterns
	tc = clin[!is.na(clin$Pattern),]
	# xColors = c("gray",lepidic_color,papillary_color,acinar_color,"darkgoldenrod","brown",solid_color,"purple")
	# xLevels = c("normal","lepidic","papillary","acinar","invasive_mucinous","cribriform","solid","micropapillary")
	xColors = c("gray",lepidic_color,papillary_color,acinar_color,solid_color)
	xLevels = c("normal","lepidic","papillary","acinar","solid")
	tc = tc[tc$Pattern %in% xLevels,]
	tc$color = dan.expand_colors(tc$Pattern,xLevels,xColors)

	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_TRACERx.pdf"),2.1,1.7 )
	for (tp in names(tps)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_TRACERx_byStage.pdf"),6,2.5 )
	for (tp in names(tps)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( tc$stage_collapsed, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs stages
	tc = clin
	tc = tc[!is.na(tc$stage_collapsed),]
	xLevels = c("Normal","I","II","III")
	xColors = c("gray","gold","orange","tomato3")
	tc$color = dan.expand_colors(tc$stage_collapsed,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Stages_TRACERx.pdf"),3.5,2.5 )
	for (tp in names(tps)){
	   x = tc[,"stage_collapsed"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "Stage", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking
	tc = clin[clin$stage_collapsed!="Normal",]
	xLevels = c("Never","Ex","Current")
	xColors = c("forestgreen","gray55","gray22")
	tc$color = dan.expand_colors(tc$Smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_TRACERx.pdf"),2.1,1.7 )
	for (tp in names(tps)){
	   x = tc[,"Smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_TRACERx_byStage.pdf"),5,2 )
	for (tp in names(tps)){
	   x = tc[,"Smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( tc$stage_collapsed, y, fill=factor(x,levels=xLevels), xlab = "", ylab = paste0(tp," score"), filllab="", plotTitle = "", signifTest = "kruskal", fillColors = xColors, comparisons = NULL, labelycoo = 0.51, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()

	## primary vs mets in Met500, TCGA, Chen, Lee, Suda, TRACERx421. Including patient matching.
	# First, let's compile joint gene expression and clinical tables for all datasets
	ge = read.table(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/MET500/M.mx.log2.txt"))
	ge_save = ge[intersect(rownames(ge),tps_universe),]
	rownames(ge) = ge$V1
	colnames(ge) = as.character(t(ge[1,]))
	ge = ge[2:nrow(ge),2:ncol(ge)]
	luad_samples = c("MO_1108-capt-SI_6356-C1N0NACXX","MO_1111-capt-SI_5897-D1DPVACXX","MO_1194-capt-SI_6862-H77PRADXX",
	  "MO_1236-capt-SI_7208-C245WACXX","MO_1250-capt-SI_7353-D2GK3ACXX","MO_1253-capt-SI_8631-C5A4VACXX","MO_1284-capt-SI_7658-D2GGHACXX",
	  "MO_1294-capt-SI_7777-C3L7FACXX","MO_1325-capt-SI_8160-C4JLUACXX","MO_1433-capt-SI_9688-C5N0MANXX","MO_1558-capt-SI_12475-C7G91ANXX",
	  "TP_2042-capt-SI_7076-C2560ACXX") 
	# MO_1150 absent
	ge = ge[,luad_samples]
	colnames(ge) = substr(luad_samples,1,7)
	mapping = read.table(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/MET500/gencode.v23.annotation.gene.probemap"), header = T)
	rownames(mapping) = mapping$id
	common = intersect(rownames(mapping),rownames(ge))
	ge = ge[common,]
	ge$gene_symbol = as.character(mapping[rownames(ge),"gene"])
	ge = ge[!duplicated(ge$gene_symbol),]
	ge = ge[!(rownames(ge) %in% c("ENSGR0000124333.14","ENSGR0000169100.12","ENSGR0000214717.9","ENSG00000279483.1") ),]
	rownames(ge) = ge$gene_symbol
	ge$gene_symbol = NULL
	logtpm_met <- as.data.frame(sapply(ge, as.numeric))
	rownames(logtpm_met) = rownames(ge)
	peek(logtpm_met)
	logtpm_tcga_primary = logtpm_tcga[,substr(colnames(logtpm_tcga),14,15)=="01"]
	logtpm_tcga_normal = logtpm_tcga[,substr(colnames(logtpm_tcga),14,15)=="11"]
	clinall = data.frame(Sample=colnames(logtpm_tcga_normal), Patient=substr(colnames(logtpm_tcga_normal),1,12), Dataset="TCGA", SampleType="Normal", BiopsySite="lung", stringsAsFactors=F)
	clinall = rbind(clinall,data.frame(Sample=colnames(logtpm_tcga_primary), Patient=substr(colnames(logtpm_tcga_primary),1,12), Dataset="TCGA", SampleType="Primary", BiopsySite="lung",stringsAsFactors=F))
	clinall = rbind(clinall,data.frame(Sample=colnames(logtpm_chen), Patient=colnames(logtpm_chen), Dataset="ChenEAS", SampleType="Primary",BiopsySite="lung",stringsAsFactors=F))
	thisclin = data.frame(Sample=colnames(logtpm_met), Patient=colnames(logtpm_met), Dataset="MET500", SampleType="Metastasis",BiopsySite=NA,stringsAsFactors=F)
	rownames(thisclin) = colnames(logtpm_met)
	clin500 = dan.read(file="/mnt/ndata/daniele/lung_multiregion/Data/MET500/Clin_met500.txt")
	clin500 = clin500[clin500$Sample_id %in% luad_samples,]
	rownames(clin500) = clin500$sample_source
	thisclin[rownames(clin500),"BiopsySite"] = clin500$biopsy_tissue
	clinall = rbind(clinall,thisclin)
	# Reading TRACERx421 data. A lot of information: metastasis-seeding regions, mutations, patterns, ...
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
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	logtpm = log2(ge[intersect(rownames(ge),tps_universe),]+0.0001)
	logtpm = logtpm[,substr(colnames(logtpm),1,8) %in% rownames(Clin)]
	Clin3 = dan.read( "/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/sampleOverview.txt" )
	Clin3$region = gsub("\\.","-",Clin3$region )
	non_metastatic = dan.read( "/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/metsFigures/data/nonMetastaticPatients.txt",header=F )
	ClinTracerx = data.frame(Sample=colnames(logtpm), Patient=substr(colnames(logtpm),1,8), Dataset="TRACERx", SampleType=NA,BiopsySite=NA,stringsAsFactors=F)
	ClinTracerx[ ClinTracerx$Patient %in% non_metastatic$V1,"SampleType" ] = "Primary"
	ClinTracerx[ ClinTracerx$Sample %in% Clin3[Clin3$sampleType=="primary","region"],"SampleType" ] = "Primary"
	ClinTracerx[ ClinTracerx$Sample %in% Clin3[Clin3$sampleType=="metastasis","region"],"SampleType" ] = "Metastasis"
	ClinTracerx[ ClinTracerx$Patient %in% non_metastatic$V1,"BiopsySite" ] = "lung"
	ClinTracerx[ ClinTracerx$Sample %in% Clin3[Clin3$sampleType=="primary","region"],"BiopsySite" ] = "lung"
	ClinTracerx[ ClinTracerx$Sample %in% Clin3[Clin3$sampleTypeDetail=="LN","region"],"BiopsySite" ] = "lymph_node"
	ClinTracerx[ ClinTracerx$Sample %in% Clin3[Clin3$sampleTypeDetail %in% c("metachronousMet","synchronousMet"),"region"],"BiopsySite" ] = "non-LN"

	sum(is.na(ClinTracerx$SampleType))
	# additionally, every sample containing 'LN' is a metastasis
	ClinTracerx[grepl("LN",ClinTracerx$Sample),"SampleType"] = "Metastasis"
	ClinTracerx[grepl("LN",ClinTracerx$Sample),"BiopsySite"] = "lymph_node"
	# and _N is normal
	ClinTracerx[grepl("_N",ClinTracerx$Sample),"SampleType"] = "Normal"
	ClinTracerx[grepl("_N",ClinTracerx$Sample),"BiopsySite"] = "lung"
	ClinTracerx = ClinTracerx[!is.na(ClinTracerx$SampleType),]
	logtpm_tracerx = logtpm[,ClinTracerx$Sample]
	clinall = rbind(clinall,ClinTracerx)
	load("/mnt/ndata/daniele/lung_multiregion/Data/Lee2020/Lee_luad_logtpm.RData")
	logtpm_lee = logtpm
	load("/mnt/ndata/daniele/lung_multiregion/Data/Lee2020/Lee_luad_Clin.RData")
	Clin = Clin[Clin$SampleType!="NormalMet",]
	logtpm_lee = logtpm_lee[,Clin$Sample]
	Clin[is.na(Clin$Metastatic_locus),"Metastatic_locus"] = "lung"
	ClinLee = data.frame(Sample=Clin$Sample, Patient=Clin$Patient, Dataset="Lee", SampleType=Clin$SampleType,BiopsySite=tolower(Clin$Metastatic_locus),stringsAsFactors=F)
	clinall = rbind(clinall,ClinLee)
	ge = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/MetastasesDatasets/","../Suda2018/Suda_fpkm.txt"))
	ge = ge[!grepl(",",ge$gene), ]
	ge2 = aggregate(.~gene,data=ge,FUN="max")
	ge = ge2[!(substr(ge2$gene,1,1) %in% c(0,1)),]
	rownames(ge) = ge$gene
	ge$gene = NULL
	logtpm_suda = ge[intersect(rownames(ge),tps_universe),]
	ClinSuda = data.frame(Sample=colnames(logtpm_suda), Patient=substr(colnames(logtpm_suda),1,2), Dataset="Suda", SampleType=ifelse(substr(colnames(logtpm_suda),4,10)=="primary","Primary","Metastasis"),BiopsySite=NA,stringsAsFactors=F)
	ClinSuda[ClinSuda$Sample %in% c( "P1_primary","P1_metLung","P2_primary","P2_metLungRML","P2_metLungLUL" ),"BiopsySite"] = "lung"
	ClinSuda[ClinSuda$Sample %in% c( "P1_metLiver","P1_metPleuraVisceral","P1_metPleuraMediastinal" ),"BiopsySite"] = "non-LN"
	ClinSuda[ClinSuda$Sample %in% c( "P2_metLungLN" ),"BiopsySite"] = "lymph_node"
	clinall = rbind(clinall,ClinSuda)
	dim(clinall)
	dtable(clinall$SampleType)
	commonz = intersect(rownames(logtpm_tcga),rownames(logtpm_chen))
	logtpm_all = cbind(logtpm_tcga[commonz,],logtpm_chen[commonz,])
	commonz = intersect(rownames(logtpm_all),rownames(logtpm_met))
	logtpm_all = cbind(logtpm_all[commonz,],logtpm_met[commonz,])
	commonz = intersect(rownames(logtpm_all),rownames(logtpm_tracerx))
	logtpm_all = cbind(logtpm_all[commonz,],logtpm_tracerx[commonz,])
	commonz = intersect(rownames(logtpm_all),rownames(logtpm_lee))
	logtpm_all = cbind(logtpm_all[commonz,],logtpm_lee[commonz,])
	commonz = intersect(rownames(logtpm_all),rownames(logtpm_suda))
	logtpm_all = cbind(logtpm_all[commonz,],logtpm_suda[commonz,])
	all(colnames(logtpm_all) %in% clinall$Sample)
	logtpm_all = logtpm_all[,clinall$Sample]
	rownames(clinall) = clinall$Sample
	rankData = rankGenes(logtpm_all)
	for (tp in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   clinall[rownames(scoredf),tp] = scoredf$TotalScore
	}
	# ChenEAS has a very weird behaviour
	clinall[clinall$BiopsySite=="lymph_node","BiopsySite"] = "LN"
	clinall[!(clinall$BiopsySite %in% c("LN","lung")),"BiopsySite"] = "other"
	clinall = clinall[clinall$Dataset!="ChenEAS",]
	save(clinall,file=paste0(OutDir,"Bulk_Tps_scores_Metastases.RData"))
	save(logtpm_all,file=paste0(OutDir,"Bulk_logtpm_Metastases.RData"))
	logtpm_all = logtpm_all[,clinall$Sample]
	tc = clinall
	tc$SampleType = paste0(tc$SampleType, "_",tc$BiopsySite)
	xLevels = c("Normal_lung","Primary_lung","Metastasis_lung","Metastasis_LN","Metastasis_other")
	xColors = c("gray","gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$SampleType,xLevels,xColors)
	# pdf( paste0(OutDir,"Bulk_Tps_vs_SampleType_AllDatasets.pdf"),7,5 )
	# for (tp in names(tps)){
	#    x = tc[,"SampleType"]
	#    y = tc[,tp]
	#    plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "SampleType", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	#    print(plot)
	# }
	# dev.off()

	tc = clinall
	tc$SampleType = paste0(tc$SampleType, "_",tc$BiopsySite)
	tc[tc$SampleType=="Normal_lung","SampleType"] = "Normal lung"
	tc[tc$SampleType=="Primary_lung","SampleType"] = "Primary tumor"
	tc[tc$SampleType=="Metastasis_LN","SampleType"] = "Mets, LN"
	tc[tc$SampleType=="Metastasis_lung","SampleType"] = "Mets, lung"
	tc[tc$SampleType=="Metastasis_other","SampleType"] = "Mets, distal"
	xLevels = c( "Normal lung" ,"Primary tumor","Mets, LN","Mets, lung","Mets, distal" )
	xColors = c( "gray","gold","tomato3","firebrick","sienna4" )
	# xLevels = c("Normal_lung","Primary_lung","Metastasis_lung","Metastasis_LN","Metastasis_other")
	# xColors = c("gray","gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$SampleType,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_SampleType_AllDatasets.pdf"),3,2 )
	for (tp in names(tps)){
	   x = tc[,"SampleType"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "SampleType", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()


	
	## same but with matched patients in Lee, Suda, TRACERx421
	clinall = clinall[clinall$SampleType!="Normal",]
	patients_multiple_samples = rowSums(dtable(clinall$Patient,clinall$SampleType)>0)
	patients_multiple_samples = names(which(patients_multiple_samples==2))
	clinall = clinall[clinall$Patient %in% patients_multiple_samples,]
	tc = clinall
	xLevels = c("Primary","Metastasis")
	xColors = c("orange","firebrick")
	tc$color = dan.expand_colors(tc$SampleType,xLevels,xColors)
	this_OutDir = paste0(OutDir,"matched_primary_metastases/")
	dir.create(this_OutDir)
	for (tp in names(tps)){
		fileName = paste0(this_OutDir,"Bulk_",tp,"_vs_SampleType_Matched_PairedDotPlot.pdf")
		x = tc[,"SampleType"]
		y = tc[,tp]
		dan.pairedDotPlot( fileName, factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", labelLines = NULL, signifTest = "wilcox.test", comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = tc$color, labelPoints="", fileWidth = 5, fileHeight = 5, linesPairings = tc$Patient )
	}

	## same but with only seeding regions and matched mets in TRACERx421
	clinall = clinall[clinall$Dataset=="TRACERx",]
	load(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/metsFigures/data/seedingRegionInfo.rda"))
	seedingRegionInfo$PrimaryRegion = gsub("\\.","-",seedingRegionInfo$PrimaryRegion )
	seedingRegionInfo = seedingRegionInfo[seedingRegionInfo$PrimaryRegion %in% rownames(clinall),]
	seedingRegionInfo = seedingRegionInfo[seedingRegionInfo$Metastasizing,]
	mets_keep = unlist(strsplit(seedingRegionInfo$Metastasis,split=";"))
	all_regions_keep = c(seedingRegionInfo$PrimaryRegion, gsub( "\\.","-",mets_keep ) )
	clinall = clinall[rownames(clinall) %in% all_regions_keep,]
	tc = clinall
	xLevels = c("Primary","Metastasis")
	xColors = c("orange","firebrick")
	tc$color = dan.expand_colors(tc$SampleType,xLevels,xColors)
	this_OutDir = paste0(OutDir,"matched_primary_metastases/")
	dir.create(this_OutDir)
	for (tp in names(tps)){
		fileName = paste0(this_OutDir,"Bulk_",tp,"_vs_SampleType_MatchedTracerxMetastasizingRegionsOnly_PairedDotPlot.pdf")
		x = tc[,"SampleType"]
		y = tc[,tp]
		dan.pairedDotPlot( fileName, factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", labelLines = NULL, signifTest = "wilcox.test", comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = tc$color, labelPoints="", fileWidth = 2.5, fileHeight = 2.5, linesPairings = tc$Patient )
	}
}

tps_across_ClinicalCharacteristics_heatmap = function(OutDir,scores,order_tps){
	library(pheatmap)
	library(ComplexHeatmap)
	library(colorRamp2)
	scores_store = scores
	## try to include normal cells. Do it also for all bulk datasets. 

	### Patient-level
	scores = scores_store
	scores[scores$SampleType=="Normal","Stage_collapsed"] = NA
	scores = scores[scores$Epi_Cancer=="Cancer",]	
	scores = scores[!is.na(scores$Stage_collapsed),]
	scores$PST = paste0(scores$Patient,"_",scores$SampleType)
	pst = dan.df(unique(scores$PST),c( "Dataset","Patient","PST","SampleType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern", order_tps ))	
	for (this_pst in rownames(pst)){
		if (sum(scores$PST==this_pst)<5){ next }
		for (tp in order_tps){
			pst[this_pst,tp] = mean(scores[ scores$PST==this_pst,tp ])
		}
		for (var in c( "Dataset","Patient","PST","SampleType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern" )){
			pst[this_pst,var] = unique(scores[ scores$PST==this_pst,var ])[1]
		}
	}
	pst = pst[!is.na(pst[,order_tps[1]]),]
	pst$Smoker = ifelse(pst$Smoking=="Never","Never","Ever")
	save(pst,file=paste0(OutDir,"PatientLevel_table.RData"))
	aa = t(pst[,order_tps])
	# annot_cols = pst[,c("SampleType","Stage_collapsed")]
	# colnames(annot_cols) = c( "SampleType","Stage" )
	# annot_colorz = list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24"))
	# pdf(paste0(OutDir,"pheatmap_PatientLevel.pdf"),width=20,height=6)
	# pheatmap(aa, annotation_colors=annot_colorz, annotation_col=annot_cols, show_colnames=F )
	# dev.off()
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel.pdf"),16,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=TRUE, use_raster = FALSE, col=colorRamp2(c(-1.5,0,1.5),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()

	# aa = t(pst[,order_tps])
	# qq = quantile(aa,seq(0,1,0.01))
	# aa[aa<qq["1%"]] = qq["1%"]
	# aa[aa>qq["99%"]] = qq["99%"]
	# pdf(paste0(OutDir,"pheatmap_PatientLevel_99percentClipping.pdf"),width=20,height=6)
	# pheatmap(aa, annotation_colors=annot_colorz, annotation_col=annot_cols, show_colnames=F )
	# dev.off()
	pst$SampleType = factor(pst$SampleType,levels=c( "Primary","Metastasis" ))
	pst = pst[order(pst$Stage_collapsed,pst$SampleType),]
	aa = t(pst[,order_tps])
	# pdf(paste0(OutDir,"pheatmap_PatientLevel_sortedStage.pdf"),width=20,height=6)
	# pheatmap(aa, annotation_colors=annot_colorz, annotation_col=annot_cols, cluster_cols=FALSE,show_colnames=F )
	# dev.off()
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_sortedStage.pdf"),16,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=F, use_raster = FALSE, col=colorRamp2(c(-1.5,0,1.5),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
	# pst = pst[order(pst$Stage_collapsed),]
	# aa = t(pst[,order_tps])
	# qq = quantile(aa,seq(0,1,0.01))
	# aa[aa<qq["1%"]] = qq["1%"]
	# aa[aa>qq["99%"]] = qq["99%"]
	# pdf(paste0(OutDir,"pheatmap_PatientLevel_sortedStage_99percentClipping.pdf"),width=20,height=6)
	# pheatmap(aa, annotation_colors=annot_colorz, annotation_col=annot_cols, cluster_cols=FALSE,show_colnames=F )
	# dev.off()

	# aa = t(scale(pst[,order_tps]))
	# # qq = quantile(aa,seq(0,1,0.01))
	# # aa[aa<qq["1%"]] = qq["1%"]
	# # aa[aa>qq["99%"]] = qq["99%"]
	# colnames(aa) = NULL
	# column_ha = HeatmapAnnotation(SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	# pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_scaled.pdf"),16,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=T, use_raster = FALSE, col=colorRamp2(c(-2.5,0,2.5),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	# dev.off()
	
	# pst = pst[order(pst$Stage_collapsed),]
	# aa = t(scale(pst[,order_tps]))
	# qq = quantile(aa,seq(0,1,0.01))
	# aa[aa<qq["1%"]] = qq["1%"]
	# aa[aa>qq["99%"]] = qq["99%"]
	# pdf(paste0(OutDir,"pheatmap_PatientLevel_scaled_sortedStage_99percentClipping.pdf"),width=20,height=6)
	# pheatmap(aa, annotation_colors=annot_colorz, annotation_col=annot_cols, cluster_cols=FALSE,show_colnames=F )
	# dev.off()

	# including AT1/AT2/AT0
	scores = scores_store[scores_store$Epi_Cancer %in% c( "AT0","AT1","AT2","Cancer"),]
	scores[scores$SampleType=="Normal","Stage_collapsed"] = NA
	scores$PST = paste0(scores$Patient,"_",scores$SampleType,"_",scores$Epi_Cancer)
	pst = dan.df(unique(scores$PST),c( "Dataset","Patient","PST","SampleType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern", order_tps ))	
	for (this_pst in rownames(pst)){
		if (sum(scores$PST==this_pst)<5){ next }
		for (tp in order_tps){
			pst[this_pst,tp] = mean(scores[ scores$PST==this_pst,tp ])
		}
		for (var in c( "Dataset","Patient","PST","SampleType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern" )){
			pst[this_pst,var] = unique(scores[ scores$PST==this_pst,var ])[1]
		}
	}
	pst = pst[!is.na(pst[,order_tps[1]]),]
	pst$Smoker = ifelse(pst$Smoking=="Never","Never","Ever")
	# annot_cols = pst[,c("SampleType","Stage_collapsed","Epi_Cancer")]
	# colnames(annot_cols) = c( "SampleType","Stage","CellType" )
	# annot_colorz = list(SampleType=c(Normal="gray",Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24"),CellType=c( AT1="chocolate",AT2="dodgerblue4",AT0="purple",Cancer="firebrick" ))
	# pdf(paste0(OutDir,"pheatmap_PatientLevel_withNormal.pdf"),width=20,height=6)
	# pheatmap(t(pst[,order_tps]), annotation_colors=annot_colorz, annotation_col=annot_cols, show_colnames=F )
	# dev.off()
	pst[pst$Epi_Cancer!="Cancer","Stage_collapsed"] = '0'
	pst = pst[!((pst$Epi_Cancer=="Cancer") & (is.na(pst$Stage_collapsed)) ),]
	pst$SampleType = factor(pst$SampleType,levels=c( "Normal","Primary","Metastasis" ))
	pst$Epi_Cancer = factor(pst$Epi_Cancer,levels=c( "AT1","AT2","AT0","Cancer" ))
	aa = t(pst[,order_tps])
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=pst$SampleType, Stage=pst$Stage_collapsed, CellType=pst$Epi_Cancer, col=list(SampleType=c(Normal="gray",Primary="#73AB84",Metastasis="#323633"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24"),CellType=c( AT1="chocolate",AT2="dodgerblue4",AT0="purple",Cancer="firebrick" )) )
	pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_withNormal.pdf"),16,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=TRUE, use_raster = FALSE, col=colorRamp2(c(-1.5,0,1.5),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
	pst = pst[order(pst$Epi_Cancer,pst$Stage_collapsed,pst$SampleType),]
	tpst = pst[!((pst$Epi_Cancer=="Cancer") & (is.na(pst$Stage_collapsed)) ),]
	aa = t(tpst[,order_tps])
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=tpst$SampleType, Stage=tpst$Stage_collapsed, CellType=tpst$Epi_Cancer, col=list(SampleType=c(Normal="gray",Primary="#73AB84",Metastasis="#323633"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24"),CellType=c( AT1="chocolate",AT2="dodgerblue4",AT0="purple",Cancer="firebrick" )) )
	pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_withNormal_SortedCellType.pdf"),16,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE, use_raster = FALSE, col=colorRamp2(c(-1.5,0,1.5),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()

	# ### Single cell-level
	# ts = scores_store[!is.na(scores_store$Stage_collapsed),]
	# ts = ts[sample(rownames(ts),length(rownames(ts))),]
	# # ts = ts[sample(rownames(ts),10000),]
	# ts = ts[ts$Epi_Cancer=="Cancer",]
	# library(fastcluster)
	# # hc = fastcluster::hclust(dist(ts))
	# # ts$hc_order = hc$order
	# hc_cols = fastcluster::hclust(dist(t(ts[,order_tps])))
	# hc_cols$labels[hc_cols$order]
	# new_order = hc_cols$labels[hc_cols$order]
	# new_order =  c( "AT2-like","AT2-Club-like","MHC-II","Club-like","lepidic-like","Stress_AP1","Stress_HSP","Stress_secreted","Unas_emp","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Translation_initiation","Unfolded_protein_response","RNA_processing" )
	# fake_counts = t(ts[,new_order])
	# fake_seu = CreateSeuratObject(counts = fake_counts, project = "fake", meta.data = ts)
	# fake_seu@meta.data$SampleType = factor(fake_seu@meta.data$SampleType,levels=c( "Primary","Metastasis" ))
	# fake_seu@meta.data$Stage_collapsed = factor(fake_seu@meta.data$Stage_collapsed,levels=c( "I","II","III","IV" ))
	# fake_seu@assays$RNA@scale.data = as.matrix(fake_seu@assays$RNA@data)
	# library(grid)
	# library(Scillus)
	# pdf(paste0(OutDir, "pheatmap_SingleCellLevel_sortedStage.pdf"),20,nrow(fake_seu)/2)
	# Scillus::plot_heatmap(dataset = fake_seu, 
	#           markers = rownames(fake_seu),
	#           sort_var = c("Stage_collapsed"),
	#           anno_var = c("Stage_collapsed","SampleType"),
	#           hm_limit = c(-1,0,1),
	#           anno_colors = list(c("#ECD072","orangered","#92140C","#1E1E24","gray"),c("#73AB84","#323633")))
	# dev.off()
	colorz = colorRamp2(c(-1,-0.5,0,0.5,1),c("dodgerblue4","dodgerblue3","white","firebrick3","firebrick4"))
	ts = scores_store[!is.na(scores_store$Stage_collapsed),]
	ts = ts[ts$Epi_Cancer=="Cancer",]
	ts = ts[sample(rownames(ts),20000),]
	# ts = ts[sample(rownames(ts),nrow(ts)),]
	ts$SampleType = factor(ts$SampleType,levels=c( "Primary","Metastasis" ))	
	ts = ts[order(ts$SampleType,ts$Stage_collapsed),]
	new_order = c( "AT2-like","MHC-II","Stress_AP1","Club-like","AT2-Club-like","lepidic-like","Stress_HSP","Stress_secreted","Unas_emp","MHC-I","Basal-like","Cell_proliferation","Translation_initiation","OxPhos","Interferon","RNA_processing","Metal", "EMT","pEMT","Hypoxia","Unfolded_protein_response")
	aa = t(ts[,new_order])
	colnames(aa) = NULL
	# column_ha = HeatmapAnnotation(SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	# pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedStage_raster10_bis2.pdf"),14,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = F,  col=colorz,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=10)
	# dev.off()
	column_ha = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 6),simple_anno_size=unit(6, "pt"),annotation_legend_param=list(SampleType=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),Stage=list(nrow = 1,grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6))),
		SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),gp=gpar(fontsize = 6) )
	pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedStage_raster10_bis2.pdf"),5,3)
	plot=ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = T,  col=colorz,heatmap_legend_param=list(title="Score",grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)), top_annotation=column_ha, raster_quality=10,row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top",annotation_legend_side="top",merge_legend = TRUE)
	dev.off()
	# pdf(paste0(OutDir,"intercs",suffix,"_Cs_vs_singlePrograms_hypergeometric_ComplexHeatmap.pdf"),4,2.5)
	# plot=ComplexHeatmap::Heatmap(intra_inter_mat2,rect_gp = gpar(col = "white", lwd = .1),cluster_rows=FALSE,column_title_side="bottom",column_title="Patient-wise single NMF latent factors (k=20)",column_title_gp = gpar(fontsize = 6),cluster_columns=FALSE, column_names_rot = 45,use_raster = F,raster_quality=10, col=rev_magma_white, heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="leftcenter",title="-log10(p-value)",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),show_column_names = FALSE,row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6))
	# print(draw(plot,heatmap_legend_side="top"))
	# dev.off()

	colorz = colorRamp2(c(-1,-0.5,0,0.5,1),c("dodgerblue4","dodgerblue3","white","firebrick3","firebrick4"))
	ts = scores_store[!is.na(scores_store$Stage_collapsed),]
	ts = ts[ts$Epi_Cancer=="Cancer",]
	# ts = ts[sample(rownames(ts),20000),]
	# ts = ts[sample(rownames(ts),nrow(ts)),]
	ts$SampleType = factor(ts$SampleType,levels=c( "Primary","Metastasis" ))	
	ts = ts[order(ts$SampleType,ts$Stage_collapsed),]
	new_order = c( "AT2-like","MHC-II","Stress_AP1","Club-like","AT2-Club-like","lepidic-like","Stress_HSP","Stress_secreted","Unas_emp","MHC-I","Basal-like","Cell_proliferation","Translation_initiation","OxPhos","Interferon","RNA_processing","Metal", "EMT","pEMT","Hypoxia","Unfolded_protein_response")
	aa = t(ts[,new_order])
	colnames(aa) = NULL
	# column_ha = HeatmapAnnotation(SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	# pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedStage_raster10_bis2.pdf"),14,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = F,  col=colorz,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=10)
	# dev.off()
	colorz2 = colorRampPalette(c("dodgerblue4","dodgerblue3","white","firebrick3","firebrick4"))
	dataset_colors = dan.colors(length(unique(ts$Dataset)))
	patient_colors = dan.colors(length(unique(ts$Patient)))
	names(dataset_colors) = sample(unique(ts$Dataset),length(unique(ts$Dataset)))
	names(patient_colors) = sample(unique(ts$Patient),length(unique(ts$Patient)))
	column_ha = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 6),simple_anno_size=unit(6, "pt"),annotation_legend_param=list(SampleType=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 1),labels_gp = gpar(fontsize = 1)),Dataset=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 1),labels_gp = gpar(fontsize = 1)),Patient=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 1),labels_gp = gpar(fontsize = 1)),Stage=list(nrow = 1,grid_height=unit(1, "pt"),grid_width=unit(1, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 1),labels_gp = gpar(fontsize = 1))),
		SampleType = ts$SampleType, Stage = ts$Stage_collapsed, Dataset=ts$Dataset, Patient=ts$Patient, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Dataset=dataset_colors, Patient = patient_colors, Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),gp=gpar(fontsize = 6) )
	pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedStage_raster10_bis2_noSampling_DatasetAnnotation.pdf"),8,3)
	plot=ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = T,  col=colorz,heatmap_legend_param=list(title="Score",grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)), top_annotation=column_ha, raster_quality=10,row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top",annotation_legend_side="top",merge_legend = TRUE)
	dev.off()





	colorz = colorRamp2(c(-1,-0.5,0,0.5,1),c("dodgerblue4","dodgerblue3","white","firebrick3","firebrick4"))
	ts = scores_store[!is.na(scores_store$Stage_collapsed),]
	ts = ts[ts$Epi_Cancer=="Cancer",]
	# new_order = c( "AT2-like","MHC-II","Stress_AP1","Club-like","AT2-Club-like","lepidic-like","Stress_HSP","Stress_secreted","Unas_emp","MHC-I","Basal-like","Cell_proliferation","Translation_initiation","OxPhos","Interferon","RNA_processing","Metal", "EMT","pEMT","Hypoxia","Unfolded_protein_response")
	ts$SampleType = factor(ts$SampleType,levels=c( "Primary","Metastasis" ))
	agg = aggregate(.~Stage_collapsed,data=ts[,c("Stage_collapsed",order_tps )],FUN='mean')
	agg$Stage_collapsed = NULL
	new_order = colnames(agg)[order(as.numeric(agg[4,])-as.numeric(agg[1,]))]
	tp_I = new_order[1]
	tp_IV = new_order[length(new_order)]

	ts = ts[order(ts$SampleType,ts$Stage_collapsed,ts[,tp_IV]-ts[,tp_I]),]
	aa = t(ts[,new_order])
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedStage_raster10_bis.pdf"),14,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE,cluster_rows=F, use_raster = TRUE,  col=colorz,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=10)
	# rect_gp = gpar(col = "white", lwd = 1)
	# plot=(ComplexHeatmap::Heatmap(aa,cluster_rows=FALSE, cluster_columns=FALSE, use_raster = T, col=colorz, heatmap_legend_param=list(title="Score",title_position="leftcenter-rot",title="Spearman R, patient-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)),top_annotation=column_ha, raster_quality=10)
	# print(draw(plot,heatmap_legend_side="right"))
	dev.off()


	# including AT1/AT2/AT0
	ts = scores_store
	ts = ts[ts$Epi_Cancer %in% c( "AT0","AT1","AT2","Cancer"),]
	ts = ts[!((ts$Epi_Cancer=="Cancer") & (is.na(ts$Stage_collapsed)) ),]
	ts[ts$Epi_Cancer!="Cancer","Stage_collapsed"] = '0'
	ts = ts[sample(rownames(ts),length(rownames(ts))),]
	ts = ts[order(ts$Epi_Cancer,ts$Stage_collapsed,ts$SampleType),]
	ts$SampleType = factor(ts$SampleType,levels=c( "Normal","Primary","Metastasis" ))
	ts$Epi_Cancer = factor(ts$Epi_Cancer,levels=c( "AT1","AT2","AT0","Cancer" ))
	aa = t(ts[,order_tps])
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=ts$SampleType, Stage=ts$Stage_collapsed, CellType=ts$Epi_Cancer, col=list(SampleType=c(Normal="gray",Primary="#73AB84",Metastasis="#323633"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24"),CellType=c( AT1="chocolate",AT2="dodgerblue4",AT0="purple",Cancer="firebrick" )) )
	pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_withNormal_sortedStage_raster10.pdf"),16,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE, use_raster = TRUE,  col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=10)
	dev.off()

	### Bulk
	## TCGA
	# with normals
	load(paste0(OutDir,"Bulk_Tps_scores_TCGA.RData"))
	Clin$Smoker = ifelse(Clin$tobacco_smoking_history==1,"Never","Ever")
	Clin$SampleType = ifelse(substr(Clin$Sample,14,15)==11,"Normal","Primary")
	Clin$Stage_collapsed = Clin$stage_collapsed
	Clin[Clin$SampleType=="Normal","Stage_collapsed"] = "0"
	Clin = Clin[!is.na(Clin$Stage_collapsed),]
	aa = t(Clin[,order_tps])
	aa = t(scale(t(aa)))
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=Clin$SampleType, Stage=Clin$Stage_collapsed, col=list(SampleType=c(Normal="gray",Primary="#73AB84"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_TCGA_withNormal.pdf"),20,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=TRUE, use_raster = FALSE, col=colorRamp2(c(-2,0,2),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
	Clin = Clin[order(Clin$SampleType,Clin$Stage_collapsed),]
	aa = t(Clin[,order_tps])
	aa = t(scale(t(aa)))
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=Clin$SampleType, Stage=Clin$Stage_collapsed, col=list(SampleType=c(Normal="gray",Primary="#73AB84"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_TCGA_withNormal_sortedStage.pdf"),20,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE, use_raster = FALSE, col=colorRamp2(c(-2,0,2),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
	load(paste0(OutDir,"Bulk_Tps_scores_TCGA.RData"))
	Clin$Smoker = ifelse(Clin$tobacco_smoking_history==1,"Never","Ever")
	Clin$SampleType = ifelse(substr(Clin$Sample,14,15)==11,"Normal","Primary")
	Clin = Clin[Clin$SampleType!="Normal",]
	Clin$Stage_collapsed = Clin$stage_collapsed
	Clin[Clin$SampleType=="Normal","Stage_collapsed"] = "0"
	Clin = Clin[!is.na(Clin$Stage_collapsed),]
	aa = t(Clin[,order_tps])
	aa = t(scale(t(aa)))
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=Clin$SampleType, Stage=Clin$Stage_collapsed, col=list(SampleType=c(Normal="gray",Primary="#73AB84"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_TCGA.pdf"),20,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=TRUE, use_raster = FALSE, col=colorRamp2(c(-2,0,2),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
	Clin = Clin[order(Clin$SampleType,Clin$Stage_collapsed),]
	aa = t(Clin[,order_tps])
	aa = t(scale(t(aa)))
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(SampleType=Clin$SampleType, Stage=Clin$Stage_collapsed, col=list(SampleType=c(Normal="gray",Primary="#73AB84"), Stage=c('0'="gray",I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_TCGA_sortedStage.pdf"),20,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE, use_raster = FALSE, col=colorRamp2(c(-2,0,2),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
	# with patterns
	load(paste0(OutDir,"Bulk_Tps_scores_TCGA.RData"))
	Clin$Smoker = ifelse(Clin$tobacco_smoking_history==1,"Never","Ever")
	Clin$SampleType = ifelse(substr(Clin$Sample,14,15)==11,"Normal","Primary")
	Clin = Clin[Clin$SampleType!="Normal",]
	Clin$Stage_collapsed = Clin$stage_collapsed
	Clin[Clin$SampleType=="Normal","Stage_collapsed"] = "0"
	Clin = Clin[!is.na(Clin$Stage_collapsed),]
	Clin = Clin[!is.na(Clin$Pattern),]
	Clin$Pattern = factor(Clin$Pattern,levels=c( "lepidic","papillary","acinar","solid" ))
	Clin = Clin[order(Clin$Pattern),]
	aa = t(Clin[,order_tps])
	aa = t(scale(t(aa)))
	colnames(aa) = NULL
	column_ha = HeatmapAnnotation(Pattern=Clin$Pattern, Stage=Clin$Stage_collapsed, col=list(Pattern=c(lepidic=lepidic_color,papillary=papillary_color,acinar=acinar_color,solid=solid_color), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	pdf(paste0(OutDir, "ComplexHeatmap_TCGA_sortedPattern.pdf"),20,5)
	ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE, use_raster = FALSE, col=colorRamp2(c(-2,0,2),c("dodgerblue4","white","firebrick3")),heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	dev.off()
}

CellStates_inference = function(OutDir, scores, tps, nstates = 4, datasets = NULL, extreme_top_quantile = 0.1){
	####### Clustering with nstates classes
	this_OutDir = OutDir
	dir.create(this_OutDir)
	### GMM
	scores[is.na(scores$Stage_collapsed) & (scores$Dataset=="wu"),"Stage_collapsed"] = "III/IV"
	scores[is.na(scores$Stage_collapsed) & (scores$Dataset=="wang"),"Stage_collapsed"] = "unknown"
	scores = scores[(scores$Epi_Cancer=="Cancer"),]
	if (!is.null(datasets)){
		scores = scores[scores$Dataset %in% datasets,]
	}
	cl_scores = scores
	scorez = scores[,names(tps)]
	seed = 123 # doesn't change much if I use different seeds
	gmm = GMM(scorez, nstates, seed = seed)
	pr1 = predict(gmm, newdata = scorez)
	pr2 = predict_GMM(data = scorez,CENTROIDS=gmm$centroids, COVARIANCE=gmm$covariance_matrices, WEIGHTS=gmm$weights)
	pr2 = data.frame(pr2$cluster_proba)
	colnames(pr2) = paste0("cs",1:nstates)
	rownames(pr2) = rownames(scorez)
	pr = predict(gmm, newdata = scorez)
	cl_scores[rownames(scorez),paste0( "cs" )] = paste0( "cs",as.numeric(pr) )

	# remapping cluster names
	# aa = table(cl_scores$Stage_collapsed,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Stage_collapsed,cl_scores[,paste0( "cs" )]))
	patient_table = table(cl_scores$Patient,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Patient,cl_scores[,paste0( "cs" )]))
	s_table = dan.df(unique(cl_scores$Stage_collapsed),colnames(patient_table))
	for (rn in rownames(s_table)){
	     s_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$Stage_collapsed==rn,"Patient"]), ]))
	}
	s_table = as.matrix(s_table)
	diffs = as.numeric(s_table["IV",]-s_table["I",])
	remapdf = data.frame(row.names = colnames(s_table),old_ids = colnames(s_table), new_ids = paste0("cs",rank(diffs)))
	cl_scores[,paste0( "cs" )] = remapdf[cl_scores[,paste0( "cs" )],"new_ids"]
	scores = cl_scores
	colnames(pr2) = remapdf[colnames(pr2),"new_ids"]

	save(pr2,file = paste0(this_OutDir,"cs",nstates,"_cluster_probabilities.RData"))

	pdf(paste0(this_OutDir,"cs",nstates,"_nCells.pdf"),nstates+1,5)
	print(barplot(table(cl_scores[rownames(scorez),paste0( "cs" )]), main = paste0(nstates," states"), xlab = "State", ylab = "Number of cells" ))
	dev.off()

	cs_seed123 = cl_scores[,paste0( "cs" )]
	for (cl in sort(unique(cl_scores[,paste0( "cs" )]))){
	  tsc = melt(cl_scores[cl_scores[,paste0( "cs" )]==cl,c(names(tps))])
	  tsc$cluster = cl
	  if (cl=="cs1"){
	    tscall = tsc
	  } else {
	    tscall = rbind(tscall,tsc)
	  }
	}

	mean_tscall = aggregate(value~.,tscall,FUN="mean")
	matt = dcast(mean_tscall,variable~cluster)
	rownames(matt) = matt$variable
	matt$variable = NULL
	matt = t(as.matrix(matt*1))
	matt2 = reorderMatt(matt)
	col.lim = c(-max(abs(c(max(matt2),min(matt2)))),+max(abs(c(max(matt2),min(matt2)))))

	matt2 = reorderMatt(matt,columns_only=TRUE)

	cellstate_profile = matt2
	patient_table = table(cl_scores$Patient,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Patient,cl_scores[,paste0( "cs" )]))

	aa = table(cl_scores$Dataset,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Dataset,cl_scores[,paste0( "cs" )]))
	dataset_table = dan.df(unique(cl_scores$Dataset),colnames(patient_table))
	for (rn in rownames(dataset_table)){
	     dataset_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$Dataset==rn,"Patient"]), ]))
	}
	dataset_table = as.matrix(dataset_table[rownames(aa),])

	aa = table(cl_scores$SampleType,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$SampleType,cl_scores[,paste0( "cs" )]))
	st_table = dan.df(unique(cl_scores$SampleType),colnames(patient_table))
	for (rn in rownames(st_table)){
	     st_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$SampleType==rn,"Patient"]), ]))
	}
	st_table = as.matrix(st_table[rownames(aa),])

	aa = table(cl_scores$Stage_collapsed,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Stage_collapsed,cl_scores[,paste0( "cs" )]))
	s_table = dan.df(unique(cl_scores$Stage_collapsed),colnames(patient_table))
	for (rn in rownames(s_table)){
	     s_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$Stage_collapsed==rn,"Patient"]), ]))
	}
	s_table = as.matrix(s_table[rownames(aa),])

	pdf(paste0(this_OutDir,"cs",nstates,"_Dataset_dots.pdf"),5,5)
	corrplot(dataset_table,is.corr=F,col=colorRampPalette(c("white","purple"))(100),tl.srt = 45,tl.col='black',addCoef.col = 'black',cl.pos='n')
	dev.off()

	pdf(paste0(this_OutDir,"cs",nstates,"_SampleType_dots.pdf"),5,3)
	corrplot(st_table,is.corr=F,col=colorRampPalette(c("white","purple"))(100),tl.srt = 45,tl.col='black',addCoef.col = 'black',cl.pos='n')
	dev.off()

	s_table = s_table[c( "I","II","III","IV" ),]
	pdf(paste0(this_OutDir,"cs",nstates,"_Stage_collapsed_dots.pdf"),5,3)
	corrplot(s_table,is.corr=F,col=colorRampPalette(c("white","purple"))(100),tl.srt = 45,tl.col='black',addCoef.col = 'black',cl.pos='n')
	dev.off()

	save(scores,file = paste0(OutDir,"tps_CellStates_",nstates,".RData"))
	save(cellstate_profile,file=paste0(this_OutDir,"cs",nstates,"_TPprofile_dots.RData"))
	pdf(paste0(OutDir,"cs",nstates,"_TPprofile_dots.pdf"),8,3)
	corrplot(cellstate_profile,is.corr=F,col=colorz_twoways,tl.srt = 45,tl.col='black',col.lim = col.lim)
	dev.off()

	extremes_profile = cellstate_profile
	for (cs in rownames(extremes_profile)){
		if (ncol(pr2)>2){
			these_extremes = rowSums(pr2[,colnames(pr2)!=cs])
			these_extremes = these_extremes[scores[scores$cs==cs,"CellID" ]]
			this_pr2 = pr2[names(these_extremes),]
			n_top = round(length(these_extremes)*extreme_top_quantile)
			these_extremes = rownames(this_pr2[order(these_extremes),])[1:n_top]
		} else {
			these_extremes = pr2[,colnames(pr2)!=cs]
			names(these_extremes) = rownames(pr2)
			these_extremes = these_extremes[scores[scores$cs==cs,"CellID" ]]
			this_pr2 = pr2[names(these_extremes),]
			n_top = round(length(these_extremes)*extreme_top_quantile)
			these_extremes = rownames(this_pr2[order(these_extremes),])[1:n_top]
		}
		extremes_profile[cs,colnames(cellstate_profile)] = as.numeric(colMeans(scores[these_extremes,colnames(cellstate_profile)]))
	}
	col.lim = c(-max(abs(c(max(extremes_profile),min(extremes_profile)))),+max(abs(c(max(extremes_profile),min(extremes_profile)))))
	pdf(paste0(OutDir,"cs",nstates,"_TPprofile_dots_extremes.pdf"),8,3)
	corrplot(reorderMatt(extremes_profile,columns_only=TRUE),is.corr=F,col=colorz_twoways,tl.srt = 45,tl.col='black',col.lim = col.lim)
	dev.off()
	save(extremes_profile,file=paste0(this_OutDir,"cs",nstates,"_TPprofile_dots_extremes.RData"))

	weighted_profile = cellstate_profile
	for (cs in rownames(weighted_profile)){
		if (ncol(pr2)>2){
			weights = -log10(rowMeans(pr2[,colnames(pr2)!=cs]))
		} else {
			weights = -log10(pr2[,colnames(pr2)!=cs])
		}
		weights = weights[scores$cs==cs]
		weighted_profile[cs,colnames(cellstate_profile)] = as.numeric( apply(scores[scores$cs==cs,colnames(cellstate_profile)],2,weighted.mean,w=weights ) )
	}
	col.lim = c(-max(abs(c(max(weighted_profile),min(weighted_profile)))),+max(abs(c(max(weighted_profile),min(weighted_profile)))))
	pdf(paste0(OutDir,"cs",nstates,"_TPprofile_dots_weighted.pdf"),8,3)
	corrplot(reorderMatt(weighted_profile,columns_only=TRUE),is.corr=F,col=colorz_twoways,tl.srt = 45,tl.col='black',col.lim = col.lim)
	dev.off()
	save(weighted_profile,file=paste0(this_OutDir,"cs",nstates,"_TPprofile_dots_weighted.RData"))
}

reorderMatt = function( matt, columns_only=FALSE ){
    if (columns_only){
    	assigned = rownames(matt)[apply(matt,2,which.max)]
		names(assigned) = names(apply(matt,2,which.max))
		matt2 = matt
		colnames(matt2) = paste0("c",1:ncol(matt) )
		index = 1
		for (cs in rownames(matt2)){
		   cs_assigned = names(assigned[assigned==cs])
		   if (length(cs_assigned)>1) { cs_assigned = names(sort(matt[cs,cs_assigned],decreasing=T)) }
		   for ( col in cs_assigned ){
		      matt2[,paste0("c",index)] = matt[,col]
		      colnames(matt2)[colnames(matt2)==paste0("c",index)] = col
		      index = index+1
		   }
		}
    } else {
		assigned = rownames(matt)[apply(matt,2,which.max)]
		names(assigned) = names(apply(matt,2,which.max))
		rowOrder = names(sort(table(assigned),decreasing=T))
		matt = matt[rowOrder,]
		matt2 = matt
		colnames(matt2) = paste0("c",1:ncol(matt) )
		index = 1
		for (cs in rownames(matt2)){
		   cs_assigned = names(assigned[assigned==cs])
		   if (length(cs_assigned)>1) { cs_assigned = names(sort(matt[cs,cs_assigned],decreasing=T)) }
		   for ( col in cs_assigned ){
		      matt2[,paste0("c",index)] = matt[,col]
		      colnames(matt2)[colnames(matt2)==paste0("c",index)] = col
		      index = index+1
		   }
		}
	}
	return(matt2)
}

CellStates_inference_hierarchical = function( OutDir, tps, solution_step1, nsubstates, extreme_top_quantile = 0.1 ){
	this_OutDir = OutDir
	dir.create(this_OutDir)
	load(file = paste0(OutDir,"tps_CellStates_",solution_step1,".RData"))
	scores$cs_level1 = scores[,paste0("cs")]
	offset = 0
	for (i in 1:solution_step1){
		these_scores = scores[(scores[,paste0("cs_level1" )]==paste0("cs",i )),]
		scorez = these_scores[,names(tps)]
		seed = 123 # doesn't change much if I use different seeds
		gmm = GMM(scorez, nsubstates[i], seed = seed)
		pr = predict(gmm, newdata = scorez)
		pr2 = predict_GMM(data = scorez,CENTROIDS=gmm$centroids, COVARIANCE=gmm$covariance_matrices, WEIGHTS=gmm$weights)
		pr2 = data.frame(pr2$cluster_proba)
		colnames(pr2) = paste0("cs",(1:nsubstates[i])+offset )
		rownames(pr2) = rownames(scorez)
		save(pr2,file = paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_solutionStep1_",i,"_cluster_probabilities.RData"))
		scores[rownames(scorez),"cs"] = paste0("cs", as.numeric(pr)+offset)
		offset = offset + nsubstates[i]
	}

	# remapping cluster names
	cl_scores = scores
	patient_table = table(cl_scores$Patient,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Patient,cl_scores[,paste0( "cs" )]))
	s_table = dan.df(unique(cl_scores$Stage_collapsed),colnames(patient_table))
	for (rn in rownames(s_table)){
	     s_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$Stage_collapsed==rn,"Patient"]), ]))
	}
	s_table = as.matrix(s_table)
	diffs = as.numeric(s_table["IV",]-s_table["I",])
	remapdf = data.frame(row.names = colnames(s_table),old_ids = colnames(s_table), new_ids = paste0("cs",rank(diffs)))
	cl_scores[,paste0( "cs" )] = remapdf[cl_scores[,paste0( "cs" )],"new_ids"]
	scores = cl_scores
	for (i in 1:solution_step1){
		load(file = paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_solutionStep1_",i,"_cluster_probabilities.RData"))
		colnames(pr2) = remapdf[colnames(pr2),"new_ids"]
		save(pr2,file = paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_solutionStep1_",i,"_cluster_probabilities.RData"))
	}

	pdf(paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_nCells.pdf"),5,5)
	barplot(dtable(scores[,"cs"]), xlab = "State", ylab = "Number of cells" )
	dev.off()

	first = TRUE
	for (cl in sort(unique(scores$cs))){
	  tsc = melt(scores[scores$cs==cl,c(names(tps))])
	  tsc$cluster = cl
	  if (first){
	    tscall = tsc
	    first = FALSE
	  } else {
	    tscall = rbind(tscall,tsc)
	  }
	}

	mean_tscall = aggregate(value~.,tscall,FUN="mean")
	matt = dcast(mean_tscall,variable~cluster)
	rownames(matt) = matt$variable
	matt$variable = NULL
	matt = t(as.matrix(matt*1))
	col.lim = c(-max(abs(c(max(matt),min(matt)))),+max(abs(c(max(matt),min(matt)))))
	cellstate_profile = matt
	cl_scores = scores

	patient_table = dtable(cl_scores$Patient,cl_scores[,paste0( "cs" )])/rowSums(dtable(cl_scores$Patient,cl_scores[,paste0( "cs" )]))
	scores[is.na(scores$Stage_collapsed) & (scores$Dataset=="wu"),"Stage_collapsed"] = "III/IV"
	scores[is.na(scores$Stage_collapsed) & (scores$Dataset=="wang"),"Stage_collapsed"] = "unknown"
	aa = table(cl_scores$Stage_collapsed,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Stage_collapsed,cl_scores[,paste0( "cs" )]))
	s_table = dan.df(unique(cl_scores$Stage_collapsed),colnames(patient_table))
	for (rn in rownames(s_table)){
	     s_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$Stage_collapsed==rn,"Patient"]), ]))
	}
	s_table = as.matrix(s_table[rownames(aa),])

	save(cellstate_profile,file=paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_dots.RData"))

	pdf(paste0(OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_dots.pdf"),8,3)
	corrplot(reorderMatt(cellstate_profile,columns_only=TRUE),is.corr=F,col=colorz_twoways,tl.srt = 45,tl.col='black',col.lim = col.lim)
	dev.off()

	patient_table = table(cl_scores$Patient,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Patient,cl_scores[,paste0( "cs" )]))

	aa = table(cl_scores$Dataset,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$Dataset,cl_scores[,paste0( "cs" )]))
	dataset_table = dan.df(unique(cl_scores$Dataset),colnames(patient_table))
	for (rn in rownames(dataset_table)){
	     dataset_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$Dataset==rn,"Patient"]), ]))
	}
	dataset_table = as.matrix(dataset_table[rownames(aa),])

	aa = table(cl_scores$SampleType,cl_scores[,paste0( "cs" )])/rowSums(table(cl_scores$SampleType,cl_scores[,paste0( "cs" )]))
	st_table = dan.df(unique(cl_scores$SampleType),colnames(patient_table))
	for (rn in rownames(st_table)){
	     st_table[rn,] = as.numeric(colMeans(patient_table[ unique(cl_scores[cl_scores$SampleType==rn,"Patient"]), ]))
	}
	st_table = as.matrix(st_table[rownames(aa),])

	pdf(paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_Dataset_dots.pdf"),5,5)
	corrplot(dataset_table,is.corr=F,col=colorRampPalette(c("white","purple"))(100),tl.srt = 45,tl.col='black',addCoef.col = 'black',cl.pos='n')
	dev.off()

	pdf(paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_SampleType_dots.pdf"),5,3)
	corrplot(st_table,is.corr=F,col=colorRampPalette(c("white","purple"))(100),tl.srt = 45,tl.col='black',addCoef.col = 'black',cl.pos='n')
	dev.off()

	s_table = s_table[c( "I","II","III","IV" ),]
	pdf(paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_Stage_collapsed_dots.pdf"),4,3)
	corrplot(s_table,is.corr=F,col=colorRampPalette(c("white","purple"))(100),tl.srt = 45,tl.col='black',addCoef.col = 'black',cl.pos='n')
	dev.off()

	scores = cl_scores
	save(scores,file = paste0(OutDir,"tps_CellStates_hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),".RData"))

	load(file = paste0(this_OutDir,"cs",solution_step1,"_TPprofile_dots_extremes.RData")) # extremes_profile
	extremes_profile_step1 = extremes_profile
	extremes_profile = cellstate_profile
	load(file = paste0(this_OutDir,"cs",solution_step1,"_TPprofile_dots_weighted.RData")) # extremes_profile
	weighted_profile_step1 = weighted_profile
	weighted_profile = cellstate_profile
	for (i in 1:solution_step1){
		load(file = paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_solutionStep1_",i,"_cluster_probabilities.RData"))
		if (ncol(pr2)==1){
			if (colnames(pr2)=="cs1"){
				extremes_profile[colnames(pr2),] = extremes_profile_step1[colnames(pr2),colnames(extremes_profile)]
				weighted_profile[colnames(pr2),] = weighted_profile_step1[colnames(pr2),colnames(weighted_profile)]	
			} else {
				extremes_profile[colnames(pr2),] = extremes_profile_step1["cs2",colnames(extremes_profile)]
				weighted_profile[colnames(pr2),] = weighted_profile_step1["cs2",colnames(weighted_profile)]	
			}
		} else {
			for (cs in colnames(pr2)){
				if (ncol(pr2)>2){
					these_extremes = rowSums(pr2[,colnames(pr2)!=cs])
					these_extremes = these_extremes[scores[scores$cs==cs,"CellID" ]]
					this_pr2 = pr2[names(these_extremes),]
					n_top = round(length(these_extremes)*extreme_top_quantile)
					these_extremes = rownames(this_pr2[order(these_extremes),])[1:n_top]
					weights = -log10(rowMeans(pr2[,colnames(pr2)!=cs]))
					names(weights) = rownames(pr2)
				} else {
					these_extremes = pr2[,colnames(pr2)!=cs]
					names(these_extremes) = rownames(pr2)
					these_extremes = these_extremes[scores[scores$cs==cs,"CellID" ]]
					this_pr2 = pr2[names(these_extremes),]
					n_top = round(length(these_extremes)*extreme_top_quantile)
					these_extremes = rownames(this_pr2[order(these_extremes),])[1:n_top]
					weights = -log10(pr2[,colnames(pr2)!=cs])
					names(weights) = rownames(pr2)
				}
				weights = weights[scores[scores$cs==cs,"CellID"]]
				weighted_profile[cs,colnames(cellstate_profile)] = as.numeric( apply(scores[scores$cs==cs,colnames(cellstate_profile)],2,weighted.mean,w=weights ) )
				extremes_profile[cs,colnames(cellstate_profile)] = as.numeric(colMeans(scores[these_extremes,colnames(cellstate_profile)]))
			}
		}
	}
	col.lim = c(-max(abs(c(max(extremes_profile),min(extremes_profile)))),+max(abs(c(max(extremes_profile),min(extremes_profile)))))
	pdf(paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_extremes.pdf"),8,3)
	corrplot(reorderMatt(extremes_profile,columns_only=TRUE),is.corr=F,col=colorz_twoways,tl.srt = 45,tl.col='black',col.lim = col.lim)
	dev.off()
	save(extremes_profile,file=paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_extremes.RData"))
	col.lim = c(-max(abs(c(max(weighted_profile),min(weighted_profile)))),+max(abs(c(max(weighted_profile),min(weighted_profile)))))
	pdf(paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_weighted.pdf"),8,3)
	corrplot(reorderMatt(weighted_profile,columns_only=TRUE),is.corr=F,col=colorz_twoways,tl.srt = 45,tl.col='black',col.lim = col.lim)
	dev.off()
	save(weighted_profile,file=paste0(this_OutDir,"hierarchical_",solution_step1,"_substates",paste0(nsubstates,collapse=""),"_TPprofile_weighted.RData"))
}

euclidean_distance = function(v1,v2){
	sqrt(sum((v1 - v2)^2))
}

Silhouette_Analysis_CellStates = function(OutDir, order_tps, range){
	sildf = dan.df(as.character(range),c( "range","mean_silhouette" ))
	for (i in range){
		dcat(i)
		load(paste0(OutDir,"../tps_CellStates_",i,".RData"))
		load(paste0(OutDir,"../cs",i,"_TPprofile_dots_weighted.RData"))
		these_cs = sort(unique(scores$cs))
		simdf = dan.df(rownames(scores),these_cs)
		for (tcs in these_cs){
			dcat(tcs,1)
			simdf[,tcs] = apply(scores[,order_tps],1,function(x) euclidean_distance(as.numeric(x),as.numeric(weighted_profile[tcs,order_tps])))
		}
		a = c()
		b = c()
		for (tcs in these_cs){
			a = c(a,simdf[scores$cs==tcs,tcs])
			b1 = simdf[scores$cs==tcs,colnames(simdf)!=tcs]
			if (ncol(simdf)>2){ b1 = as.numeric(apply(b1,1,min)) }
			b = c(b,b1)
		}
		sil = (b-a)/apply(cbind(b,a),1,max)
		sil = mean(sil)
		sildf[as.character(i),"range"] = i
		sildf[as.character(i),"mean_silhouette"] = sil
	}	
	save(sildf,file=paste0( OutDir,"sildf.RData" ))
	dan.scatterplot( paste0( OutDir,"sildf_scatterplot.pdf" ), sildf$range, sildf$mean_silhouette, fill = NULL, xlab = "Number of clusters", ylab = "mean silhouette score", dotSize = 3, plotBisector = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 2.5, fileHeight = 2 )
}

Silhouette_Analysis_CellStates_Hierarchical = function(OutDir, order_tps, range1, range2){
	sildf = dan.df(0,c( "range1","range2","mean_silhouette" ))
	for (i in range1){
		dcat(i)
		for (j in range2){
			if ((i==1) & (j==1)){ next }
			dcat(j,1)
			load(paste0(OutDir,"../tps_CellStates_hierarchical_2_substates",i,j,".RData"))
			load(paste0(OutDir,"../hierarchical_2_substates",i,j,"_TPprofile_weighted.RData"))
			these_cs = sort(unique(scores$cs))
			simdf = dan.df(rownames(scores),these_cs)
			for (tcs in these_cs){
				dcat(tcs,1)
				simdf[,tcs] = apply(scores[,order_tps],1,function(x) euclidean_distance(as.numeric(x),as.numeric(weighted_profile[tcs,order_tps])))
			}
			a = c()
			b = c()
			for (tcs in these_cs){
				a = c(a,simdf[scores$cs==tcs,tcs])
				b1 = simdf[scores$cs==tcs,colnames(simdf)!=tcs]
				if (ncol(simdf)>2){ b1 = as.numeric(apply(b1,1,min)) }
				b = c(b,b1)
			}
			sil = (b-a)/apply(cbind(b,a),1,max)
			sil = mean(sil)
			this_sildf = data.frame(row.names=paste0("i",i,"j",j),range1=i,range2=j,mean_silhouette=sil)
			sildf = rbind(sildf,this_sildf)
		}
	}
	save(sildf,file=paste0( OutDir,"sildf_hierarchical.RData" ))
	sildf = sildf[order(-sildf$mean_silhouette),]	
	dan.scatterplot( paste0( OutDir,"sildf_hierarchical_scatterplot.pdf" ), sildf$range1, sildf$range2, fill = sildf$mean_silhouette, xlab = "Number of subclusters (h1)", ylab = "Number of subclusters (h2)", filllab="Mean silhouette", dotSize = 1, fillColors_continuous=colorz_solid, coord_fixed=T, fileWidth = 2.5, fileHeight = 2.5 )
}

CellStates_2D_asymmetrical_colorby = function(OutDir, scores, x, y, xlab, ylab, prefix, order_tps, cs_map){
	fixed_yl = min(y)
	fixed_yr = max(y)
	fixed_xl = min(x)
	fixed_xr = max(x)
	fixed_yl = -0.215
	fixed_yr = 0.28
	fixed_xl = -0.36
	fixed_xr = 0.3
	library(viridis)
	scores$x = x
	scores$y = y
	## by highest tp
	order_tpz = order_tps[ !( order_tps %in% "Unas_emp" ) ]
	scores$highest_tp = remapping_tps_names(order_tpz[as.numeric(apply( scores[,order_tpz],1,FUN=which.max  ))])
	fileName = paste0( OutDir,prefix,"HighestTp.pdf" )
	order_tps_colorz = c( "midnightblue","dodgerblue4","dodgerblue3","dodgerblue2","darkgreen","forestgreen","green3","firebrick4","firebrick3","sienna1","sienna3","sienna4","goldenrod1","goldenrod3","gray25","gray55","gray85","magenta4","magenta3","lightpink3" )
	dan.scatterplot( fileName, scores$x, scores$y, fill=factor(scores[,"highest_tp"],levels=order_tpz), xlab=xlab, ylab=ylab, fillColors = order_tps_colorz, filllab="", dotSize = 0.1,coord_fixed = T, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 5, fileHeight = 5 )
	x = scores$x
	y = scores$y
	fill=factor(scores[,"highest_tp"],levels=remapping_tps_names(order_tpz))
	fillColors = order_tps_colorz
	filllab=""
	pdf( fileName, width = 3.7, height = 2.5, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(size = 0.4, alpha = 0.7,stroke=0)
	p = p + scale_color_manual( values=fillColors,drop = FALSE )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p+ theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + guides(color=guide_legend(ncol =2,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	fileName = paste0( OutDir,prefix,"HighestTpNA.pdf" )
	pdf( fileName, width = 3.7, height = 2.5, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(shape=NA,size = 0, alpha = 0,stroke=0)
	p = p + scale_color_manual( values=fillColors,drop = FALSE )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p+ theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + guides(color=guide_legend(ncol =2,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()

	fileName = paste0( OutDir,prefix,"byStage.pdf" )
	stages = c( "I","II","III","IV" )
	stages_colorz = c( "#ECD072","orangered","#92140C","#1E1E24" )
	scorez = scores[scores$Stage_collapsed %in% c( "I","II","III","IV" ),]
	fill=factor(scorez[,"Stage_collapsed"],levels=c( "I","II","III","IV" ))
	x = scorez$x
	y = scorez$y
	fillColors = stages_colorz
	filllab="Stage"
	# dan.scatterplot( fileName, x=scores$x, y=scores$y, fill=fill, xlab=xlab, ylab=ylab, fillColors = fillColors, filllab=filllab, dotSize = 0.1,coord_fixed = T, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 5, fileHeight = 5 )
	pdf( fileName, width = 2.8, height = 2.8, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(size = 0.4, alpha = 0.7,stroke = 0)
	p = p + scale_color_manual( values=fillColors,drop = FALSE )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p + theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + guides(color=guide_legend(ncol =1,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()

	fileName = paste0( OutDir,prefix,"bySampleType.pdf" )
	scorez = scores
	sampletypes = c( "Primary","Metastasis" )
	st_colorz = c( "#73AB84","#323633" )
	x = scorez$x
	y = scorez$y
	fill=factor(scorez[,"SampleType"],levels=c( "Primary","Metastasis" ))
	fillColors = st_colorz
	filllab="SampleType"
	pdf( fileName, width = 2.8, height = 2.8, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(size = 0.4, alpha = 0.7,stroke = 0)
	p = p + scale_color_manual( values=fillColors,drop = FALSE )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p+ theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + guides(color=guide_legend(ncol =1,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()

	fileName = paste0( OutDir,prefix,"CellStates.pdf" )
	x = scores$x
	y = scores$y
	fill=factor(scores$cs_new,levels=cs_map$alias)
	fillColors = cs_map$colorz
	filllab="Cell state"
	pdf( fileName, width = 2.8, height = 2.8, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(size = 0.4, alpha = 0.7,stroke = 0)
	p = p + scale_color_manual( values=fillColors,drop = FALSE )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p+ theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + guides(color=guide_legend(ncol =1,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()

	load( file=paste0(OutDir,"../../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	scores$cyto = annot[rownames(scores),"cyto"]
	fileName = paste0( OutDir,prefix,"byCyto.pdf" )
	scorez = scores
	x = scorez$x
	y = scorez$y
	fill=scorez$cyto
	fillColors = cs_map$colorz
	filllab="CytoTRACE"
	pdf( fileName, width = 2.5, height = 2.5, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(size = 0.4, alpha = 0.7,stroke = 0)
	p = p + scale_colour_gradientn( colours=viridis(100) )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p+ theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	fileName = paste0( OutDir,prefix,"byCytoNA.pdf" )
	pdf( fileName, width = 2.5, height = 2.5, useDingbats = F )
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(shape=NA,size = 0, alpha = 0,stroke=0)
	p = p + scale_colour_gradientn( colours=viridis(100) )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p+ theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()

	

	#### Giovanni's 4 stages - primary only
	fileName = paste0( OutDir,prefix,"density_stageI.pdf" )
	scorez1 = scores[((scores$Stage_collapsed %in% c( "I" ))) & (scores$SampleType=="Primary"),]
	scorez2 = scores[(!(scores$Stage_collapsed %in% c( "I" ))) & (scores$SampleType=="Primary"),]
	pdf( fileName, width = 2, height = 2, useDingbats = F,pointsize=6 )
	par(lwd = 0.5)
	plot(scorez2$x,scorez2$y, pch=19, col="gray", cex=0.1, xlim=c(fixed_xl,fixed_xr), ylim=c(fixed_yl,fixed_yr), xlab="", ylab="",xaxt="n",yaxt="n") 
	densplot(scorez1$x, scorez1$y, points=T, pch = 19, cex = 0.1, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	axis(2, at = c(-0.2,0.0,0.2), las=2,lwd=0.5)
	axis(1, at = c(-0.3,0.0,0.3),lwd=0.5)
	dev.off()
	fileName = paste0( OutDir,prefix,"density_stageII.pdf" )
	scorez1 = scores[((scores$Stage_collapsed %in% c( "II" ))) & (scores$SampleType=="Primary"),]
	scorez2 = scores[(!(scores$Stage_collapsed %in% c( "II" ))) & (scores$SampleType=="Primary"),]
	pdf( fileName, width = 2, height = 2, useDingbats = F,pointsize=6 )
	par(lwd = 0.5)
	plot(scorez2$x,scorez2$y, pch=19, col="gray", cex=0.1, xlim=c(fixed_xl,fixed_xr), ylim=c(fixed_yl,fixed_yr), xlab="", ylab="",xaxt="n",yaxt="n") 
	densplot(scorez1$x, scorez1$y, points=T, pch = 19, cex = 0.1, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	axis(2, at = c(-0.2,0.0,0.2), las=2,lwd=0.5)
	axis(1, at = c(-0.3,0.0,0.3),lwd=0.5)
	dev.off()
	fileName = paste0( OutDir,prefix,"density_stageIII.pdf" )
	scorez1 = scores[((scores$Stage_collapsed %in% c( "III" ))) & (scores$SampleType=="Primary"),]
	scorez2 = scores[(!(scores$Stage_collapsed %in% c( "III" ))) & (scores$SampleType=="Primary"),]
	pdf( fileName, width = 2, height = 2, useDingbats = F,pointsize=6 )
	par(lwd = 0.5)
	plot(scorez2$x,scorez2$y, pch=19, col="gray", cex=0.1, xlim=c(fixed_xl,fixed_xr), ylim=c(fixed_yl,fixed_yr), xlab="", ylab="",xaxt="n",yaxt="n") 
	densplot(scorez1$x, scorez1$y, points=T, pch = 19, cex = 0.1, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	axis(2, at = c(-0.2,0.0,0.2), las=2,lwd=0.5)
	axis(1, at = c(-0.3,0.0,0.3),lwd=0.5)
	dev.off()
	fileName = paste0( OutDir,prefix,"density_stageIV.pdf" )
	scorez1 = scores[((scores$Stage_collapsed %in% c( "IV" ))) & (scores$SampleType=="Primary"),]
	scorez2 = scores[(!(scores$Stage_collapsed %in% c( "IV" ))) & (scores$SampleType=="Primary"),]
	pdf( fileName, width = 2, height = 2, useDingbats = F,pointsize=6 )
	par(lwd = 0.5)
	plot(scorez2$x,scorez2$y, pch=19, col="gray", cex=0.1, xlim=c(fixed_xl,fixed_xr), ylim=c(fixed_yl,fixed_yr), xlab="", ylab="",xaxt="n",yaxt="n") 
	densplot(scorez1$x, scorez1$y, points=T, pch = 19, cex = 0.1, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	axis(2, at = c(-0.2,0.0,0.2), las=2,lwd=0.5)
	axis(1, at = c(-0.3,0.0,0.3),lwd=0.5)
	dev.off()

	fileName = paste0( OutDir,prefix,"density_primary.pdf" )
	scorez1 = scores[(scores$SampleType %in% c( "Primary" )),]
	scorez2 = scores[!(scores$SampleType %in% c( "Primary" )),]
	pdf( fileName, width = 2, height = 2, useDingbats = F,pointsize=6 )
	par(lwd = 0.5)
	plot(scorez2$x,scorez2$y, pch=19, col="gray", cex=0.1, xlim=c(fixed_xl,fixed_xr), ylim=c(fixed_yl,fixed_yr), xlab="", ylab="",xaxt="n",yaxt="n") 
	densplot(scorez1$x, scorez1$y, points=T, pch = 19, cex = 0.1, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	axis(2, at = c(-0.2,0.0,0.2), las=2,lwd=0.5)
	axis(1, at = c(-0.3,0.0,0.3),lwd=0.5)
	dev.off()
	fileName = paste0( OutDir,prefix,"density_metastasis.pdf" )
	scorez1 = scores[(scores$SampleType %in% c( "Metastasis" )),]
	scorez2 = scores[!(scores$SampleType %in% c( "Metastasis" )),]
	pdf( fileName, width = 2, height = 2, useDingbats = F,pointsize=6 )
	par(lwd = 0.5)
	plot(scorez2$x,scorez2$y, pch=19, col="gray", cex=0.1, xlim=c(fixed_xl,fixed_xr), ylim=c(fixed_yl,fixed_yr), xlab="", ylab="",xaxt="n",yaxt="n") 
	densplot(scorez1$x, scorez1$y, points=T, pch = 19, cex = 0.1, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	axis(2, at = c(-0.2,0.0,0.2), las=2,lwd=0.5)
	axis(1, at = c(-0.3,0.0,0.3),lwd=0.5)
	dev.off()

	# dan.scatterplot( fileName, x=scores$x, y=scores$y, fill=fill, xlab=xlab, ylab=ylab, fillColors = fillColors, filllab=filllab, dotSize = 0.1,coord_fixed = T, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 5, fileHeight = 5 )
	
	p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( "" ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	p = p + geom_point(size = 0.4, alpha = 0.7,stroke = 0)
	p = p + scale_color_manual( values=fillColors,drop = FALSE )
	p = p + ylim(fixed_yl, fixed_yr)
	p = p + xlim(fixed_xl, fixed_xr)
	p = p + theme(text = element_text(size=14))
	p = p + coord_fixed()
	p = p + guides(color=guide_legend(ncol =1,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()


	


	## by stage
	fileName = paste0( OutDir,prefix,"byStage.pdf" )
	stages = c( "I","II","III","IV" )
	stages_colorz = c( "#ECD072","orangered","#92140C","#1E1E24" )
	scorez = scores[scores$Stage_collapsed %in% c( "I","II","III","IV" ),]
	dan.scatterplot( fileName, scorez$x, scorez$y, fill=factor(scorez[,"Stage_collapsed"],levels=c( "I","II","III","IV" )), xlab=xlab, ylab=ylab, fillColors = stages_colorz, filllab="", dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	for (this_st in stages){
		fileName = paste0( OutDir,prefix,"byStage_",this_st,".pdf" )
		scorez = scores[scores$Stage_collapsed %in% this_st,]
		dan.scatterplot( fileName, scorez$x, scorez$y, fill=factor(scorez[,"Stage_collapsed"],levels=c( "I","II","III","IV" )), xlab=xlab, ylab=ylab, fillColors = stages_colorz, filllab="", dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	}
	## by sample type
	fileName = paste0( OutDir,prefix,"bySampleType.pdf" )
	scorez = scores
	sampletypes = c( "Primary","Metastasis" )
	st_colorz = c( "#73AB84","#323633" )
	dan.scatterplot( fileName, scorez$x, scorez$y, fill=factor(scorez[,"SampleType"],levels=c( "Primary","Metastasis" )), xlab=xlab, ylab=ylab, fillColors = st_colorz, filllab="", dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	for (this_st in sampletypes){
		fileName = paste0( OutDir,prefix,"bySampleType_",this_st,".pdf" )
		scorez = scores[scores$SampleType %in% this_st,]
		dan.scatterplot( fileName, scorez$x, scorez$y, fill=factor(scorez[,"SampleType"],levels=c( "Primary","Metastasis" )), xlab=xlab, ylab=ylab, fillColors = st_colorz, filllab="", dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	}
	## by tp
	for (tp in order_tps){
		fileName = paste0( OutDir,prefix,"tp_",tp,".pdf" )
		scorez = scores[order(scores[,tp]),]
		dan.scatterplot( fileName, scorez$x, scorez$y, fill=scorez[,tp], xlab=xlab, ylab=ylab, filllab=tp, fillColors_continuous=colorz_solid, dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	}
	## by cs
	fileName = paste0( OutDir,prefix,"CellStates.pdf" )
	dan.scatterplot( fileName, scores$x, scores$y, fill=factor(scores$cs_new,levels=cs_map$alias), xlab=xlab, ylab=ylab,  filllab="Cell state", fillColors = cs_map$colorz, dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	fileName = paste0( OutDir,prefix,"CellStates_cs1.pdf" )
	dan.scatterplot( fileName, scores$x[scores$cs=="cs1"], scores$y[scores$cs=="cs1"], fill=factor(scores$cs_new[scores$cs=="cs1"],levels=cs_map$alias), xlab=xlab, ylab=ylab,  filllab="Cell state", fillColors = cs_map$colorz, dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	fileName = paste0( OutDir,prefix,"CellStates_cs2.pdf" )
	dan.scatterplot( fileName, scores$x[scores$cs=="cs2"], scores$y[scores$cs=="cs2"], fill=factor(scores$cs_new[scores$cs=="cs2"],levels=cs_map$alias), xlab=xlab, ylab=ylab,  filllab="Cell state", fillColors = cs_map$colorz, dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	fileName = paste0( OutDir,prefix,"CellStates_cs3.pdf" )
	dan.scatterplot( fileName, scores$x[scores$cs=="cs3"], scores$y[scores$cs=="cs3"], fill=factor(scores$cs_new[scores$cs=="cs3"],levels=cs_map$alias), xlab=xlab, ylab=ylab,  filllab="Cell state", fillColors = cs_map$colorz, dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	fileName = paste0( OutDir,prefix,"CellStates_shaded_cs1.pdf" )
	pdf( fileName, width = 10, height = 10 , useDingbats = F )
	thisy = scores[scores$x<0,"y"]
	thisx = scores[scores$x<0,"x"]
	fillColors_continuous = colorRampPalette(c(cs_map["cs1","colorz"],"gray77"))(100)
	p = ggplot(mapping = aes(y = thisy, x = thisx, color = thisx)) + geom_point(size = 0.1, alpha = 0.7) + scale_colour_gradientn( colours=fillColors_continuous )
	p = p + xlim(fixed_xl,fixed_xr) + ylim(fixed_yl,fixed_yr) + xlab( xlab ) + ylab( ylab ) + theme_classic() + theme(legend.position="none") + theme(text = element_text(size=14)) + coord_fixed()
	print(p)
	dev.off()
	fileName = paste0( OutDir,prefix,"CellStates_shaded_cs2.pdf" )
	pdf( fileName, width = 10, height = 10 , useDingbats = F )
	thisy = scores[(scores$x>=0) & (scores$y<0),"y"]
	thisx = scores[(scores$x>=0) & (scores$y<0),"x"]
	fillColors_continuous = colorRampPalette(c(cs_map["cs2","colorz"],"gray77"))(100)
	p = ggplot(mapping = aes(y = thisy, x = thisx, color = thisy*thisx)) + geom_point(size = 0.1, alpha = 0.7) + scale_colour_gradientn( colours=fillColors_continuous )
	p = p + xlim(fixed_xl,fixed_xr) + ylim(fixed_yl,fixed_yr) + xlab( xlab ) + ylab( ylab ) + theme_classic() + theme(legend.position="none") + theme(text = element_text(size=14)) + coord_fixed()
	print(p)
	dev.off()
	fileName = paste0( OutDir,prefix,"CellStates_shaded_cs3.pdf" )
	pdf( fileName, width = 10, height = 10 , useDingbats = F )
	thisy = scores[(scores$x>=0) & (scores$y>=0),"y"]
	thisx = scores[(scores$x>=0) & (scores$y>=0),"x"]
	fillColors_continuous = colorRampPalette(c(cs_map["cs3","colorz"],"gray77"))(100)
	p = ggplot(mapping = aes(y = thisy, x = thisx, color = -thisy*thisx)) + geom_point(size = 0.1, alpha = 0.7) + scale_colour_gradientn( colours=fillColors_continuous )
	p = p + xlim(fixed_xl,fixed_xr) + ylim(fixed_yl,fixed_yr) + xlab( xlab ) + ylab( ylab ) + theme_classic() + theme(legend.position="none") + theme(text = element_text(size=14)) + coord_fixed()
	print(p)
	dev.off()

	## by cyto
	load( file=paste0(OutDir,"../../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	scores$cyto = annot[rownames(scores),"cyto"]
	fileName = paste0( OutDir,prefix,"byCyto.pdf" )
	scorez = scores
	dan.scatterplot( fileName, scorez$x, scorez$y, fill=scorez$cyto, xlab=xlab, ylab=ylab, fillColors_continuous=viridis(100), filllab="CytoTRACE", dotSize = 0.1,coord_fixed = TRUE, xlimLeft=fixed_xl,xlimRight=fixed_xr,ylimLeft=fixed_yl,ylimRight=fixed_yr, fileWidth = 2.8, fileHeight = 2.8 )
	
	pdf(paste0( OutDir,prefix,"DensityPlot.pdf" ),2.8,2.8)
	densplot(scores$x, scores$y, pch = 19, cex = 0.25, xlab=xlab, ylab=ylab, palette = viridis_pal(option = 'B')(256), xlim=c(fixed_xl,fixed_xr),ylim=c(fixed_yl,fixed_yr))
	dev.off()
}

CellStates_stackedbars_heatmaps = function( OutDir,scores,order_tps,cs_map ){
	library(colorRamp2)
	library(ComplexHeatmap)
	scores_store = scores
	colorscalez = colorRamp2(c(-1,-0.25,0,0.25,1),c("dodgerblue4","dodgerblue3","white","red","firebrick3"))
	ts = scores_store[(scores_store$Stage_collapsed %in% c( "I","II","III","IV" )) %in% c(T),]
	ts = ts[ts$Epi_Cancer=="Cancer",]
	ts = ts[sample(rownames(ts),20000),]
	ts$SampleType = factor(ts$SampleType,levels=c( "Primary","Metastasis" ))
	ts$cs_level1 = ifelse(ts$cs_level1=="cs1","Alveolar","Dedifferentiated")
	ts$CellState = factor(ts$cs_level1,levels=c("Alveolar","Dedifferentiated"))
	ts = ts[order(ts$CellState,ts$Stage_collapsed,ts$SampleType),]
	new_order = c( "AT2-like","MHC-II","Stress_AP1","Club-like","AT2-Club-like","lepidic-like","Stress_HSP","Stress_secreted","Unas_emp","MHC-I","Basal-like","Cell_proliferation","Translation_initiation","OxPhos","Interferon","RNA_processing","Metal", "EMT","pEMT","Hypoxia","Unfolded_protein_response")
	aa = t(ts[,new_order])
	colnames(aa) = NULL
	# column_ha = HeatmapAnnotation(CellState=ts$CellState, SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(CellState=c('Alveolar'="dodgerblue4",'Dedifferentiated'="firebrick4" ),SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	# pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedCellStates_TwoStates_raster10.pdf"),16,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE,cluster_rows=FALSE, use_raster = TRUE,  col=colorscalez,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=10)
	# dev.off()
	column_ha = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 6),simple_anno_size=unit(6, "pt"),annotation_legend_param=list(CellState=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),SampleType=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),Stage=list(nrow = 1,grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6))),
		CellState=ts$CellState,SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(CellState=c('Alveolar'="dodgerblue4",'Dedifferentiated'="firebrick4" ),SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),gp=gpar(fontsize = 6) )
	pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedCellStates_TwoStates_raster10.pdf"),5,3)
	plot=ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = T,  col=colorscalez,heatmap_legend_param=list(title="Score",grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)), top_annotation=column_ha, raster_quality=10,row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top",annotation_legend_side="top",merge_legend = TRUE)
	dev.off()
	ts$CellState = "Alveolar"
	ts[ts$cs=="cs2","CellState"] = "Proliferative"
	ts[ts$cs=="cs3","CellState"] = "Hypoxic"
	ts$CellState = factor(ts$CellState,levels=c("Alveolar","Proliferative","Hypoxic"))
	ts = ts[order(ts$CellState,ts$Stage_collapsed,ts$SampleType),]
	# column_ha = HeatmapAnnotation(CellState=ts$CellState, SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(CellState=c('Alveolar'="dodgerblue4",'Proliferative'="firebrick3",'Hypoxic'="sienna4" ),SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")) )
	# pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedCellStates_ThreeStates_raster10.pdf"),16,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE,cluster_rows=FALSE, use_raster = TRUE,  col=colorscalez,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=10)
	# dev.off()
	column_ha = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 6),simple_anno_size=unit(6, "pt"),annotation_legend_param=list(CellState=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),SampleType=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),Stage=list(nrow = 1,grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6))),
		CellState=ts$CellState,SampleType = ts$SampleType, Stage = ts$Stage_collapsed, col=list(CellState=c('Alveolar'="dodgerblue4",'Proliferative'="firebrick3",'Hypoxic'="sienna4" ),SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),gp=gpar(fontsize = 6) )
	pdf(paste0(OutDir, "ComplexHeatmap_SingleCellLevel_sortedCellStates_ThreeStates_raster10.pdf"),8,3)
	plot=ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = T,  col=colorscalez,heatmap_legend_param=list(title="Score",grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)), top_annotation=column_ha, raster_quality=10,row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top",annotation_legend_side="top",merge_legend = TRUE)
	dev.off()

	# patient level with barplot annotation
	scores = scores_store
	scores = scores[(scores$Stage_collapsed %in% c( "I","II","III","IV" )) %in% c(T),]
	scores[scores$SampleType=="Normal","Stage_collapsed"] = NA
	scores = scores[scores$Epi_Cancer=="Cancer",]	
	scores = scores[!is.na(scores$Stage_collapsed),]
	scores$PST = paste0(scores$Patient,"_",scores$SampleType)
	pst = dan.df(unique(scores$PST),c( "Dataset","Patient","PST","SampleType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern", order_tps ))	
	for (this_pst in rownames(pst)){
		if (sum(scores$PST==this_pst)<5){ next }
		for (tp in order_tps){
			pst[this_pst,tp] = mean(scores[ scores$PST==this_pst,tp ])
		}
		for (var in c( "Dataset","Patient","PST","SampleType","Stage_collapsed","Epi_Cancer","Smoking","histologic_pattern" )){
			pst[this_pst,var] = unique(scores[ scores$PST==this_pst,var ])[1]
		}
	}
	pst = pst[!is.na(pst[,order_tps[1]]),]
	pst$Smoker = ifelse(pst$Smoking=="Never","Never","Ever")
	pst$ratio_cs1 = NA
	pst$ratio_cs2 = NA
	pst$ratio_cs3 = NA
	for (rn in rownames(pst)){
		tss = scores[scores$Patient==pst[rn,"Patient"],] # no duplicates here
		pst[rn,"ratio_cs1"] = sum(tss$cs=="cs1")/nrow(tss)
		pst[rn,"ratio_cs2"] = sum(tss$cs=="cs2")/nrow(tss)
		pst[rn,"ratio_cs3"] = sum(tss$cs=="cs3")/nrow(tss)
	}
	pst$ratio_cs23 = pst$ratio_cs2+pst$ratio_cs3
	# aa = t(pst[,new_order])
	# colnames(aa) = NULL
	# colorscalez = colorRamp2(c(-1.5,-0.5,0,0.5,1.5),c("dodgerblue4","dodgerblue3","white","red","firebrick3"))
	# column_ha = HeatmapAnnotation(SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),
	# 	CellState = anno_barplot(pst[,c( "ratio_cs1","ratio_cs23" )],gp=gpar(fill=c("dodgerblue4","firebrick4" )),bar_width=1,border=F ) )
	# pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_TwoStates.pdf"),16,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=TRUE,cluster_rows=FALSE, use_raster = FALSE, col=colorscalez,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	# dev.off()
	# column_ha = HeatmapAnnotation(SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),
	# 	CellState = anno_barplot(pst[,c( "ratio_cs1","ratio_cs2","ratio_cs3" )],gp=gpar(fill=c("dodgerblue4","firebrick3","sienna4" )),bar_width=1,border=F ) )
	# pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_ThreeStates.pdf"),16,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=TRUE,cluster_rows=FALSE, use_raster = FALSE, col=colorscalez,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	# dev.off()

	pst = pst[order(-pst$ratio_cs1),]
	aa = t(pst[,new_order])
	colnames(aa) = NULL
	colorscalez = colorRamp2(c(-1.5,-0.5,0,0.5,1.5),c("dodgerblue4","dodgerblue3","white","red","firebrick3"))
	column_ha = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 6),simple_anno_size=unit(6, "pt"),annotation_legend_param=list(SampleType=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),Stage=list(nrow = 1,grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6))),
		CellState= anno_barplot(pst[,c( "ratio_cs1","ratio_cs23" )],gp=gpar(fill=c("dodgerblue4","firebrick4" )),bar_width=1,border=F,height=unit(24, "pt") ),SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),gp=gpar(fontsize = 6) )
	pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_SortedCellStates_TwoStates.pdf"),8,3)
	plot=ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = F,  col=colorscalez,heatmap_legend_param=list(title="Score",grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)), top_annotation=column_ha, raster_quality=10,row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top",annotation_legend_side="top",merge_legend = TRUE)
	dev.off()

	pst = pst[order(-pst$ratio_cs1+pst$ratio_cs3),]
	aa = t(pst[,new_order])
	colnames(aa) = NULL
	# column_ha = HeatmapAnnotation(SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),
	# 	CellState = anno_barplot(pst[,c( "ratio_cs1","ratio_cs2","ratio_cs3" )],gp=gpar(fill=c("dodgerblue4","firebrick3","sienna4" )),bar_width=1,border=F ) )
	# pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_SortedCellStates_ThreeStates.pdf"),16,5)
	# ComplexHeatmap::Heatmap(aa,cluster_columns=FALSE,cluster_rows=FALSE, use_raster = FALSE, col=colorscalez,heatmap_legend_param=list(title="Score"), top_annotation=column_ha, raster_quality=1)
	# dev.off()

	column_ha = HeatmapAnnotation(annotation_name_gp= gpar(fontsize = 6),simple_anno_size=unit(6, "pt"),annotation_legend_param=list(SampleType=list(nrow = 1,title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),Stage=list(nrow = 1,grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6))),
		CellState= anno_barplot(pst[,c( "ratio_cs1","ratio_cs2","ratio_cs3" )],gp=gpar(fill=c("dodgerblue4","firebrick3","sienna4" )),bar_width=1,border=F,height=unit(24, "pt") ),SampleType = pst$SampleType, Stage = pst$Stage_collapsed, col=list(SampleType=c(Primary="#73AB84",Metastasis="#323633"), Stage=c(I="#ECD072",II="orangered",III="#92140C",IV="#1E1E24")),gp=gpar(fontsize = 6) )
	pdf(paste0(OutDir, "ComplexHeatmap_PatientLevel_SortedCellStates_ThreeStates.pdf"),8,3)
	plot=ComplexHeatmap::Heatmap(aa,cluster_columns=F,cluster_rows=F, use_raster = F,  col=colorscalez,heatmap_legend_param=list(title="Score",grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)), top_annotation=column_ha, raster_quality=10,row_names_gp = grid::gpar(fontsize = 6))
	draw(plot,heatmap_legend_side="top",annotation_legend_side="top",merge_legend = TRUE)
	dev.off()
}

densplot = function(x,y,points = FALSE, pch=19, cex=1, palette = viridis_pal()(256),
          xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), xlab='x', ylab='y',main='', log="", log_dens = FALSE){
  library(scales)
  # xlim2 = c( min(c(min(x),min(y))),max(c(max(x),max(y))) ) 
  df = data.frame(x,y)
  d = densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens = col2rgb(d)[1,] + 1L
  cols = palette
  df$logdens = log2(df$dens + 1)
  df$col = cols[df$dens]
  if(log_dens)
    df$col = cols[df$logdens]
  df=df[order(df$dens),]
  if(points)
    points(df$x,df$y, pch=pch, col=df$col, cex=cex)
  else
    plot(df$x,df$y, pch=pch, col=df$col, cex=cex, xlim=xlim, ylim=xlim, xlab=xlab, ylab=ylab, main=main, log=log) 
}

CellStates_2D_asymmetrical = function(OutDir, scores, profile2, profile3, cs_map, order_tps){
	if (!is.null(cs_map)){
		scores$cs_new = NA
		for (rn in rownames(cs_map)){
			scores[scores$cs==rn,"cs_new"] = cs_map[rn,"alias"]
			scores[scores$cs==rn,"cs_colorz"] = cs_map[rn,"colorz"]
			scores[scores$cs==rn,"cs_colorz_level1"] = cs_map[rn,"colorz_level1"]
		}
	}
	load(file = paste0(OutDir,"../../tps_discovery/tps.RData"))
	scorez = scores[,names(tps)]
	cordf2 = data.frame(cor(t(scorez[,names(tps)]),t(profile2[,names(tps)]),method="pearson" ),stringsAsFactors=F)
	cordf3 = data.frame(cor(t(scorez[,names(tps)]),t(profile3[,names(tps)]),method="pearson" ),stringsAsFactors=F)
	x = cordf2[,"cs2"]-cordf2[,"cs1"]
	y = cordf3[,"cs3"]-cordf3[,"cs2"]
	xlab = "- Cor(Alveolar), + Cor(Proliferative/Hypoxic)"
	ylab = "- Cor(Proliferative), + Cor(Hypoxic)"
	prefix = "cor_"
	CellStates_2D_asymmetrical_colorby(OutDir, scores, x, y, xlab, ylab, prefix, order_tps, cs_map)

	euclidean_similarity <- function(v1,v2){
		1/(1+sqrt(sum((v1 - v2)^2)))
	}

	prefix = "EuclSim_"
	load( file=paste0( OutDir,prefix,"_simdf.RData" ))
	simdf = dan.df(rownames(scorez),c( "cs2_1","cs2_2","cs3_2","cs3_3" ))
	simdf$cs2_1 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile2["cs1",names(tps)])))
	simdf$cs2_2 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile2["cs2",names(tps)])))
	simdf$cs3_2 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile3["cs2",names(tps)])))
	simdf$cs3_3 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile3["cs3",names(tps)])))
	x = simdf[,"cs2_2"]-simdf[,"cs2_1"]
	y = simdf[,"cs3_3"]-simdf[,"cs3_2"]
	xlab = "ES (DD - Alveolar-like)"
	ylab = "ES (DD2 - DD1)"
	prefix = "EuclSim_"
	CellStates_2D_asymmetrical_colorby(OutDir, scores, x, y, xlab, ylab, prefix, order_tps, cs_map)
	simdf$x = x
	simdf$y = y
	save(simdf, file=paste0( OutDir,prefix,"_simdf.RData" ))

	simdf = dan.df(rownames(scorez),c( "cs2_1","cs2_2","cs3_2","cs3_3" ))
	simdf$cs2_1 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile3["cs1",names(tps)])))
	simdf$cs2_2 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(colMeans(profile3[c("cs2","cs3"),names(tps)]))) )
	simdf$cs3_2 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile3["cs2",names(tps)])))
	simdf$cs3_3 = apply(scorez,1,function(x) euclidean_similarity(x,as.numeric(profile3["cs3",names(tps)])))
	x = simdf[,"cs2_2"]-simdf[,"cs2_1"]
	y = simdf[,"cs3_3"]-simdf[,"cs3_2"]
	xlab = "- EuclSim(Alveolar), + EuclSim(Proliferative/Hypoxic)"
	ylab = "- EuclSim(Proliferative), + EuclSim(Hypoxic)"
	prefix = "EuclSimProfile3only_"
	CellStates_2D_asymmetrical_colorby(OutDir, scores, x, y, xlab, ylab, prefix, order_tps, cs_map)
	simdf$x = x
	simdf$y = y
	save(simdf, file=paste0( OutDir,prefix,"_simdf.RData" ))
}

CellStates_TriangularPlot = function(OutDir, scores, cellstate_profile, cs1_name="cs1", cs2_name="cs2", cs3_name="cs3"){
	library(ggtern)
	load(file = paste0(OutDir,"../../tps_discovery/tps.RData"))
	scorez = scores[,names(tps)]
	cordf = data.frame(cor(t(scorez[,names(tps)]),t(cellstate_profile[,names(tps)]),method="spearman" ),stringsAsFactors=F)
	cordf = (cordf+1)/2
	cordf = cordf/rowSums(cordf)*100
	pdf( paste0(OutDir,prefix,"TriangularPlot_singlecells.pdf"),7,7 )
	print(ggtern(data=cordf, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=0.1) + theme_rgbw() + tern_limits(  ))
	dev.off()

	# finding top 100 for each state
	topp = as.matrix(scorez[,colnames(cellstate_profile)]) %*% t(cellstate_profile)
	topp_cs1 = rownames(topp[order(topp[,"cs1"],decreasing=TRUE),])[1:1000]
	topp_cs2 = rownames(topp[order(topp[,"cs2"],decreasing=TRUE),])[1:1000]
	topp_cs3 = rownames(topp[order(topp[,"cs3"],decreasing=TRUE),])[1:1000]

	for (var in c( "SampleType","Stage_collapsed","Dataset" )){
		dcat(var)
		print(dtable(scores[topp_cs1,var]))
		print(dtable(scores[topp_cs2,var]))
		print(dtable(scores[topp_cs3,var]))
	}

	commonz = intersect(topp_cs1,topp_cs2)
	commonz = c(commonz,intersect(topp_cs1,topp_cs3))
	commonz = c(commonz,intersect(topp_cs2,topp_cs3))
	topp_cs1 = topp_cs1[!( topp_cs1 %in% commonz )]
	topp_cs2 = topp_cs2[!( topp_cs2 %in% commonz )]
	topp_cs3 = topp_cs3[!( topp_cs3 %in% commonz )]

	extremes_profile = cellstate_profile
	extremes_profile["cs1",] = colMeans( scorez[topp_cs1,colnames(cellstate_profile)] )
	extremes_profile["cs2",] = colMeans( scorez[topp_cs2,colnames(cellstate_profile)] )
	extremes_profile["cs3",] = colMeans( scorez[topp_cs3,colnames(cellstate_profile)] )

	cordf = data.frame(cor(t(scorez[,names(tps)]),t(extremes_profile[,names(tps)]),method="spearman" ),stringsAsFactors=F)
	cordf = (cordf+1)/2
	cordf = cordf/rowSums(cordf)*100
	pdf( paste0(OutDir,prefix,"TriangularPlot_singlecells_extremesProfile.pdf"),7,7 )
	print(ggtern(data=cordf, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=0.1) + theme_rgbw() + tern_limits(  ))
	dev.off()
}

CellStates_PatientLevel_TriangularPlot = function( OutDir, scores, prefix, cs1_name="cs1", cs2_name="cs2", cs3_name="cs3" ){
	library(ggtern)
	# Patient-level aggregation
	load(paste0("data/extendedAtlas_clin.RData"))
	clin = clin[!duplicated(clin$Patient),]
	clin = clin[clin$Patient %in% unique(scores$Patient),]
	rownames(clin) = clin$Patient
	for (rn in rownames(clin)){
		clin[rn,"n_cs1"] = nrow(scores[(scores$Patient==rn) & (scores$cs=="cs1"),])
		clin[rn,"n_cs2"] = nrow(scores[(scores$Patient==rn) & (scores$cs=="cs2"),])
		clin[rn,"n_cs3"] = nrow(scores[(scores$Patient==rn) & (scores$cs=="cs3"),])
		clin[rn,"n_total"] = nrow(scores[(scores$Patient==rn),])
		clin[rn,"cs1"] = clin[rn,"n_cs1"]/clin[rn,"n_total"]*100
		clin[rn,"cs2"] = clin[rn,"n_cs2"]/clin[rn,"n_total"]*100
		clin[rn,"cs3"] = clin[rn,"n_cs3"]/clin[rn,"n_total"]*100
	}
	pdf( paste0(OutDir,prefix,"TriangularPlot.pdf"),4.5,4.5 )
	plot = ggtern(data=clin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2,color='gray33',alpha=0.8) + theme_nomask() +
		theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
		Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
	print(plot)
	dev.off()
	for (stage in c( "I","II","III","IV" )){
		tclin = clin[(clin$Stage_collapsed==stage) %in% c(T),]
		pdf( paste0(OutDir,prefix,"TriangularPlot_Stage",stage,".pdf"),5.6,5.6 )
		plot = ggtern(data=tclin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2.5,color='gray33',alpha=0.8) + theme_nomask() +
			theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
			Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
		print(plot)
		dev.off()
	}
	for (st in c( "Primary","Metastasis" )){
		tclin = clin[(clin$SampleType==st) %in% c(T),]
		pdf( paste0(OutDir,prefix,"TriangularPlot_st",st,".pdf"),5.6,5.6 )
		plot = ggtern(data=tclin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2.5,color='gray33',alpha=0.8) + theme_nomask() +
			theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
			Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
		print(plot)
		dev.off()
	}
	tclin = clin[(clin$Stage_collapsed %in% c( "I","II","III","IV" )) %in% c(T),]
	tclin$colorz = "#ECD072"
	tclin[tclin$Stage_collapsed=="II","colorz"] = "orangered"
	tclin[tclin$Stage_collapsed=="III","colorz"] = "#92140C"
	tclin[tclin$Stage_collapsed=="IV","colorz"] = "#1E1E24"
	pdf( paste0(OutDir,prefix,"TriangularPlot_byStage.pdf"),4.5,4.5 )
	plot = ggtern(data=tclin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2,alpha=0.8,color=tclin$colorz) + scale_fill_identity() + theme_nomask() +
		theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
		Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
	print(plot)
	dev.off()
	# st_colorz = c( "#73AB84","#323633" )
	tclin = clin
	tclin$colorz = "#73AB84"
	tclin[tclin$SampleType=="Metastasis","colorz"] = "#323633"
	tclin = tclin[order(tclin$SampleType,decreasing=T),]
	pdf( paste0(OutDir,prefix,"TriangularPlot_bySampleType.pdf"),4.5,4.5 )
	plot = ggtern(data=tclin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2,alpha=0.8,color=tclin$colorz) + scale_fill_identity() + theme_nomask() +
		theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
		Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
	print(plot)
	dev.off()

	tclin = clin[!is.na(clin$Smoking),]
	tclin[(tclin$Smoking!="Never"),"Smoking"] = "Ever"
	tclin$colorz = "gray66"
	tclin[tclin$Smoking=="Ever","colorz"] = "black"
	pdf( paste0(OutDir,prefix,"TriangularPlot_bySmoking.pdf"),5.6,5.6 )
	plot = ggtern(data=tclin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2.5,alpha=0.8,color=tclin$colorz) + scale_fill_identity() + theme_nomask() +
		theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
		Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
	print(plot)
	dev.off()

	for (gene in c( "EGFR","KRAS","TP53","STK11","KEAP1" )){
		tclin = clin[!is.na(clin[,paste0( "mut_",gene )]),]
		tclin$colorz = "gray66"
		tclin[tclin[,paste0( "mut_",gene )]=="mut","colorz"] = "black"
		pdf( paste0(OutDir,prefix,"TriangularPlot_byMut",gene,".pdf"),5.6,5.6 )
		plot = ggtern(data=tclin, aes(x=cs1,y=cs2,z=cs3)) + geom_point(size=2.5,alpha=0.8,color=tclin$colorz) + scale_fill_identity() + theme_nomask() +
			theme_custom(col.T="firebrick3",col.L="dodgerblue4",col.R="sienna4", tern.panel.background = "white") +
			Llab(cs1_name) + Tlab(cs2_name) + Rlab(cs3_name) + theme_showarrows() + Larrowlab("") + Tarrowlab("") + Rarrowlab("")
		print(plot)
		dev.off()
	}
}

CellStates_signatures = function( OutDir,scores,cs_map ){
	load( file = paste0(OutDir,"../../extendedAtlas_seu_Cancer.RData") )
	load( file = paste0(OutDir,"../../tps_discovery/tps_universe.RData") )
	all(rownames(seu@meta.data) %in% rownames(scores))
	seu@meta.data$cs = cs_map[scores[rownames(seu@meta.data),"cs"],"alias"]
	seu@meta.data$cs_level1 = scores[rownames(seu@meta.data),"cs_level1"]
	scores$cs_new = NA
	for (rn in rownames(cs_map)){
		scores[scores$cs==rn,"cs_new"] = cs_map[rn,"alias"]
		scores[scores$cs==rn,"cs_colorz"] = cs_map[rn,"colorz"]
		scores[scores$cs==rn,"cs_colorz_level1"] = cs_map[rn,"colorz_level1"]
	}
	Idents(seu) = seu@meta.data$cs_level1
	mm = FindAllMarkers(seu,features=intersect(tps_universe,rownames(seu)),only.pos=TRUE)
	mm = mm[ mm$p_val_adj<0.01,]
	save( mm,file=paste0(OutDir,"mm_Wilcox_cs_level1.RData") )
	cs_signatures = list()
	cs_signatures[[ "cs1" ]] = mm[mm$cluster=="cs1","gene"]
	cs_signatures[[ "cs2" ]] = mm[mm$cluster=="cs2","gene"]
	save(cs_signatures, file=paste0( OutDir,"wilcox_csLevel1_signatures_permissive.RData" ))
	# positive check
	load( file=paste0( OutDir,"wilcox_csLevel1_signatures_permissive.RData" ))
	seu@meta.data = seu@meta.data[,!(colnames(seu@meta.data) %in% c( "cs1","cs2" ))]
	seu = AddModuleScore(seu,cs_signatures)
	colnames(seu@meta.data)[substr(colnames(seu@meta.data),1,4)=="Clus"] = names(cs_signatures)
	pdf( paste0(OutDir,"wilcox_csLevel1_signatures_permissive_across_cs_boxplots.pdf"),7,5 )
	for (n in names(cs_signatures)){
		scores$sign = seu@meta.data[rownames(scores),n]
		x = factor(scores[,"cs_level1"],levels=c( "cs1","cs2" ))
		y = scores$sign
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(n," state score"), filllab = "", plotTitle = "", signifTest = NULL,comparisons = NULL, ylimLeft=ylimLeft,ylimRight=ylimRight,labelycoo = 1, xColors = cs_map$colorz_level1[c(1,2)], jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
	}
	dev.off()

	Idents(seu) = seu@meta.data$cs
	mm = FindAllMarkers(seu,features=intersect(tps_universe,rownames(seu)),only.pos=TRUE)
	mm = mm[ mm$p_val_adj<0.01,]
	save( mm,file=paste0(OutDir,"mm_Wilcox_cs_level2.RData") )
	mm_ph = FindMarkers(seu,features=intersect(tps_universe,rownames(seu)),ident.1="Proliferative",ident.2="Hypoxic")
	mm_ph = mm_ph[mm_ph$p_val_adj<0.01,]
	mm_ap = FindMarkers(seu,features=intersect(tps_universe,rownames(seu)),ident.1="Alveolar",ident.2="Proliferative")
	mm_ap = mm_ap[mm_ap$p_val_adj<0.01,]
	mm_ah = FindMarkers(seu,features=intersect(tps_universe,rownames(seu)),ident.1="Alveolar",ident.2="Hypoxic")
	mm_ah = mm_ah[mm_ah$p_val_adj<0.01,]
	cutoff_fc = 0
	mm = mm[!(mm$gene %in% mm$gene[duplicated(mm$gene)] ),]
	am = intersect(rownames(mm_ap[mm_ap$avg_log2FC>0,]),rownames(mm_ah[mm_ah$avg_log2FC>0,]))
	pm = intersect(rownames(mm_ap[mm_ap$avg_log2FC<0,]),rownames(mm_ph[mm_ph$avg_log2FC>0,]))
	hm = intersect(rownames(mm_ph[mm_ph$avg_log2FC<0,]),rownames(mm_ah[mm_ah$avg_log2FC<0,]))
	length(intersect(mm[mm$cluster=="Alveolar","gene"],am))/length(union(mm[mm$cluster=="Alveolar","gene"],am))
	length(intersect(mm[mm$cluster=="Proliferative","gene"],pm))/length(union(mm[mm$cluster=="Proliferative","gene"],pm))
	length(intersect(mm[mm$cluster=="Hypoxic","gene"],hm))/length(union(mm[mm$cluster=="Hypoxic","gene"],hm))

	cs_signatures = list()
	cs_signatures[[ "Alveolar" ]] = am
	cs_signatures[[ "Proliferative" ]] = pm
	cs_signatures[[ "Hypoxic" ]] = hm
	save(cs_signatures, file=paste0( OutDir,"wilcox_cs_signatures_restrictive.RData" ))
	cs_signatures = list()
	cs_signatures[[ "Alveolar" ]] = mm[mm$cluster=="Alveolar","gene"]
	cs_signatures[[ "Proliferative" ]] = mm[mm$cluster=="Proliferative","gene"]
	cs_signatures[[ "Hypoxic" ]] = mm[mm$cluster=="Hypoxic","gene"]
	save(cs_signatures, file=paste0( OutDir,"wilcox_cs_signatures_permissive.RData" ))

	# positive check
	load( file=paste0( OutDir,"wilcox_cs_signatures_permissive.RData" ))
	seu@meta.data = seu@meta.data[,!(colnames(seu@meta.data) %in% c( "Alveolar","Proliferative","Hypoxic" ))]
	seu = AddModuleScore(seu,cs_signatures)
	colnames(seu@meta.data)[substr(colnames(seu@meta.data),1,4)=="Clus"] = names(cs_signatures)
	pdf( paste0(OutDir,"wilcox_cs_signatures_permissive_across_cs_boxplots.pdf"),7,5 )
	for (n in names(cs_signatures)){
		scores$sign = seu@meta.data[rownames(scores),n]
		x = factor(scores[,"cs_new"],levels=cs_map$alias)
		y = scores$sign
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(n," state score"), filllab = "", plotTitle = "", signifTest = NULL,comparisons = NULL, ylimLeft=ylimLeft,ylimRight=ylimRight,labelycoo = 1, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
	}
	dev.off()
	load( file=paste0( OutDir,"wilcox_cs_signatures_restrictive.RData" ))
	seu@meta.data = seu@meta.data[,!(colnames(seu@meta.data) %in% c( "Alveolar","Proliferative","Hypoxic" ))]
	seu = AddModuleScore(seu,cs_signatures)
	colnames(seu@meta.data)[substr(colnames(seu@meta.data),1,4)=="Clus"] = names(cs_signatures)
	pdf( paste0(OutDir,"wilcox_cs_signatures_restrictive_across_cs_boxplots.pdf"),7,5 )
	for (n in names(cs_signatures)){
		scores$sign = seu@meta.data[rownames(scores),n]
		x = factor(scores[,"cs_new"],levels=cs_map$alias)
		y = scores$sign
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(n," state score"), filllab = "", plotTitle = "", signifTest = NULL,comparisons = NULL, ylimLeft=ylimLeft,ylimRight=ylimRight,labelycoo = 1, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
	}
	dev.off()

	load( file = paste0(OutDir,"../../extendedAtlas_counts.RData") )
	load( file = paste0(OutDir,"../../extendedAtlas_annot.RData") )
	annot$cs = cs_map[scores[rownames(annot),"cs"],"alias"]
	annot$cs_level1 = scores[rownames(annot),"cs_level1"]
	counts = counts[rownames(counts) %in% tps_universe,]
	tseu = CreateSeuratObject(counts = counts, project = "this", meta.data = annot)
	tseu = NormalizeData(tseu)
	tseu = FindVariableFeatures(tseu)
	tseu = ScaleData(tseu)
	tseu = RunPCA(tseu)
	Idents(tseu) = tseu@meta.data$cs_level1
	mm = FindAllMarkers(tseu,test.use='MAST',latent.vars='Dataset')
	mm = mm[ mm$p_val_adj<0.01,]
	save( mm,file=paste0(OutDir,"mm_Mast_CorrectingDataset_cs_level1.RData") )
	Idents(tseu) = tseu@meta.data$cs
	mm = FindAllMarkers(tseu,test.use='MAST',latent.vars='Dataset')
	mm = mm[ mm$p_val_adj<0.01,]
	save( mm,file=paste0(OutDir,"mm_Mast_CorrectingDataset_cs_level2.RData") )
	mm = mm[!(mm$gene %in% mm$gene[duplicated(mm$gene)] ),]
	mm = mm[ mm$avg_log2FC>0.25, ]
	dtable(mm$cluster)
	cs_signatures = list()
	cs_signatures[[ "Alveolar" ]] = mm[mm$cluster=="Alveolar","gene"]
	cs_signatures[[ "Proliferative" ]] = mm[mm$cluster=="Proliferative","gene"]
	cs_signatures[[ "Hypoxic" ]] = mm[mm$cluster=="Hypoxic","gene"]
	save(cs_signatures, file=paste0( OutDir,"MastDataset_cs_signatures_permissive.RData" ))
	seu@meta.data = seu@meta.data[,!(colnames(seu@meta.data) %in% c( "Alveolar","Proliferative","Hypoxic" ))]
	seu = AddModuleScore(seu,cs_signatures)
	colnames(seu@meta.data)[substr(colnames(seu@meta.data),1,4)=="Clus"] = names(cs_signatures)
	pdf( paste0(OutDir,"MastDataset_cs_signatures_permissive_across_cs_boxplots.pdf"),7,5 )
	for (n in names(cs_signatures)){
		scores$sign = seu@meta.data[rownames(scores),n]
		x = factor(scores[,"cs_new"],levels=cs_map$alias)
		y = scores$sign
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(n," state score"), filllab = "", plotTitle = "", signifTest = NULL,comparisons = NULL, ylimLeft=ylimLeft,ylimRight=ylimRight,labelycoo = 1, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
	}
	dev.off()
}

TME_vs_tps_cs = function(OutDir, scores, prefix, order_tps, noMets = F){
	library(compositions)
	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	tps = tps[order_tps]
	DataDir = "data/"
	load(paste0("data/extendedAtlas_clin.RData"))

	clin = clin[clin$SampleType!="Normal",]
	if (noMets) { clin = clin[clin$SampleType!="Metastasis",] }
	clin$PST = paste0(clin$Patient,"_",substr(clin$SampleType,1,1))
	scores$PST = paste0(scores$Patient,"_",substr(scores$SampleType,1,1))
	clin = clin[!duplicated(clin$PST),]
	rownames(clin) = clin$PST
	codf = dan.df(clin$PST,c(  "n_allcells", "n_cancercells", "n_tmecells", paste0("n_",sort(unique(scores$cs)) ) ))
	codf = cbind(clin,codf)

	for (pst in rownames(codf)){
	  ts = scores[scores$PST==pst,]
	  codf[ pst, "n_cancercells" ] = nrow(ts)
	  for (thiscs in sort(unique(scores$cs))){
	  	codf[ pst, paste0("n_", thiscs) ] = sum(ts$cs==thiscs)
	  }
	  for (tp in names(tps)){
	    codf[ pst, paste0( "tp_",tp ) ] = mean( ts[,tp] )
	  }
	}
	CellTypeDir = paste0("data/cell_typing/non_malignant/")
	ctdf = dan.df(0,c( "annd_level_1","annd_level_2","annd_level_3","stiched" ))
	for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu")){
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_annot_annd.RData"))
		annot$stiched = paste0(annot$annd_level_1,"_",annot$annd_level_2,"_",annot$annd_level_3)
		ctdf = rbind(ctdf,annot[,c("annd_level_1","annd_level_2","annd_level_3","stiched" )])
		ctdf = ctdf[!duplicated(ctdf$stiched),]	
	}
	ctdf = ctdf[!duplicated(ctdf$stiched),]
	ctdf = ctdf[!(ctdf$annd_level_1=="unclear"), ]
	ctdf = ctdf[order(ctdf$annd_level_1,ctdf$annd_level_2,ctdf$annd_level_3),]
	all_ct_1 = unique(ctdf$annd_level_1)
	all_ct_2 = unique(ctdf$annd_level_2)
	all_ct_3 = unique(ctdf$annd_level_3)

	# macro flavours
	macro_flavours = c( "IFN_TAMs","Inflam_TAMs","LA_TAMs","Angio_TAMs","Reg_TAMs","Prolif_TAMs","Alveolar_RTM_like","MT_RTM_like" )
	load("data/cell_typing/non_malignant/extended_atlas/integrating_macrophages/AllDatasets_integrated_MacroSubsets_md.RData")
	macrof = md[,c("Dataset","Patient","SampleType","annd_level_3",macro_flavours)]
	load("data/cell_typing/non_malignant/withIntegration/macrophages/AllDatasets_integrated_MacroSubsets_md.RData")
	macrof = rbind(macrof,md[,c("Dataset","Patient","SampleType","annd_level_3",macro_flavours)])
	macrof$Patient[macrof$Dataset %in% c( "hu","wang" )] = paste0( macrof$Dataset[macrof$Dataset %in% c( "hu","wang" )],macrof$Patient[macrof$Dataset %in% c( "hu","wang" )] )
	macrof$PST = paste0(macrof$Patient,"_",substr(macrof$SampleType,1,1))
	colnames(macrof)[colnames(macrof) %in% macro_flavours] = paste0("mf_",macro_flavours  )

	for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu")){
	  dcat(dataset)
	  load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_annot_annd.RData"))
	  annot = annot[annot$SampleType!="Normal",]
	  if (noMets) { clin = clin[clin$SampleType!="Metastasis",] }
	  annot[annot$SampleType=="Tumor","SampleType"] = "Primary"
	  annot$PST = paste0(annot$Patient,"_",substr(annot$SampleType,1,1))
	  for (pst in rownames(codf)[codf$Dataset==dataset] ){
	  	ma = macrof[macrof$PST==pst,]
	    an = annot[annot$PST==pst,]
	    an = an[!(an$annd_level_1=="unclear"),]
	    if (nrow(an)==0){
	      dcat( paste0(pst, " not found in tme") , 1)
	      next
	    }
	    if (nrow(ma)==0){
	      dcat( paste0(pst, " not found in macro macro flavours") , 1)
	      next
	    }
	    codf[ pst, "n_tmecells" ] = nrow(an)
	    for (ct in all_ct_1){
	      codf[ pst, paste0("n_level1_",ct) ] = sum(an$annd_level_1==ct)
	    }
	    for (ct in all_ct_2){
	      codf[ pst, paste0("n_level2_",ct) ] = sum(an$annd_level_2==ct)
	    }
	    for (ct in all_ct_3){
	      codf[ pst, paste0("n_level3_",ct) ] = sum(an$annd_level_3==ct)
	    }
	    for (thismf in paste0("mf_",macro_flavours  )){
	      codf[ pst, thismf ] = mean(ma[,thismf])
	    }
	  }
	}
	codf$n_allcells = codf$n_cancercells+codf$n_tmecells
	codf = codf[!is.na(rowSums(codf[,paste0( "tp_",order_tps )])),]
	sum(codf$n_allcells)
	sum(codf$n_cancercells)

	codf_save = codf

	#### 0), just N cancer cells vs N tme cells
	x = log10(codf$n_cancercells)
	y = log10(codf$n_tmecells)
	fileName = paste0(OutDir, prefix,"n_cancercells_vs_n_tmecells_fillDataset.pdf")
	plotTitle = paste0("Spearman r = ",signif(cor(x,y,method="spearman"),2),", p-val = ",signif(cor.test(x,y,method="spearman")$p.value,2))
	dan.scatterplot( fileName, x, y, fill = codf$Dataset, xlab = "n_cancercells, log10", ylab = "n_tmecells, log10", filllab = "Dataset", fillColors = NULL, plotTitle = plotTitle, dotLabels = NULL, dotSize = 3, plotBisector = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, coord_fixed = F, fileWidth = 7, fileHeight = 7 )
	x = log10(codf$n_cancercells)
	y = log10(codf$n_tmecells)
	fileName = paste0(OutDir, prefix,"n_cancercells_vs_n_tmecells_fillStage.pdf")
	plotTitle = paste0("Spearman r = ",signif(cor(x,y,method="spearman"),2),", p-val = ",signif(cor.test(x,y,method="spearman")$p.value,2))
	dan.scatterplot( fileName, x, y, fill = codf$Stage_collapsed, xlab = "n_cancercells, log10", ylab = "n_tmecells, log10", filllab = "Stage", fillColors = NULL, plotTitle = plotTitle, dotLabels = NULL, dotSize = 3, plotBisector = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, coord_fixed = F, fileWidth = 7, fileHeight = 7 )

	#### I), vs total
	radf = codf
	radf$n_allcells = NULL
	columnz = colnames(radf)[substr(colnames(radf),1,2)=="n_"]
	radf[,columnz] = sweep(radf[,columnz],1,codf$n_allcells,FUN="/")
	radf = radf[,substr(colnames(radf),nchar(colnames(radf))-7,nchar(colnames(radf)))!="_unclear" ]

	# ratio cancer cells vs n cancer cells
	x = log10(codf$n_cancercells)
	y = (radf$n_cancercells)
	fileName = paste0(OutDir, prefix,"n_cancercells_vs_r_cancercells_fillStage.pdf")
	plotTitle = paste0("Spearman r = ",signif(cor(x,y,method="spearman"),2),", p-val = ",signif(cor.test(x,y,method="spearman")$p.value,2))
	dan.scatterplot( fileName, x, y, fill = codf$Stage_collapsed, xlab = "n_cancercells, log10", ylab = "r_cancercells", filllab = "Stage", fillColors = NULL, plotTitle = plotTitle, dotLabels = NULL, dotSize = 3, plotBisector = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, coord_fixed = F, fileWidth = 7, fileHeight = 7 )

	this_OutDir = paste0(OutDir,"I_vsTotal/")
	dir.create(this_OutDir)
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level1.pdf"),ncol(cordf),nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_2
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs macro_flavours
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,3)=="mf_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_macro_flavours.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()

	# csRatios vs level_1
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level1.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# csRatios vs level_2
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# csRatios vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# csRatios vs macro_flavours
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,3)=="mf_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_macro_flavours.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()

	## with clr transformation, level1
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,2)=="n_"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_allcells,FUN="/")
	table(rowSums(clrdf[,c( "n_cancercells","n_tmecells" )]))
	table(rowSums(clrdf[,c( "n_cancercells", colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level1"]  )]))
	table(rowSums(clrdf[,c( "n_cancercells", colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level2"]  )]))

	table(rowSums(clrdf[,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level3"]  )]))
	# mean tps vs level_1
	for (rn in rownames(clrdf)){
	  clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level1"]  )] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level1"]  )]))
	}
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level1.pdf"),ncol(cordf),nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level1.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	## with clr transformation, level2
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,2)=="n_"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_allcells,FUN="/")
	for (rn in rownames(clrdf)){
	  clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level2"]  )] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level2"]  )]))
	}
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	## with clr transformation, level3
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,2)=="n_"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_allcells,FUN="/")
	for (rn in rownames(clrdf)){
	  clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level3"]  )] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))], colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level3"]  )]))
	}
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# vs macro_flavours
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,3)=="mf_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_macro_flavours.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()


	#### II), cs vs total cancer, tme vs total tme
	radf = codf
	radf$n_allcells = NULL
	columnz = colnames(radf)[substr(colnames(radf),1,4)=="n_cs"]
	radf[,columnz] = sweep(radf[,columnz],1,codf$n_cancercells,FUN="/")
	columnz = colnames(radf)[substr(colnames(radf),1,7)=="n_level"]
	radf[,columnz] = sweep(radf[,columnz],1,codf$n_tmecells,FUN="/")
	radf = radf[,substr(colnames(radf),nchar(colnames(radf))-7,nchar(colnames(radf)))!="_unclear" ]

	this_OutDir = paste0(OutDir,"II_vsTotalCancer_TotalTME/")
	dir.create(this_OutDir)
	# mean tps vs level_1
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level1.pdf"),ncol(cordf),nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_2
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()

	# mean tps vs level_1
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level1.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_2
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()

	## with clr transformation
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,4)=="n_cs"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_cancercells,FUN="/")
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,7)=="n_level"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_tmecells,FUN="/")
	table(rowSums(clrdf[,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))] )]))
	table(rowSums(clrdf[,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level1"]) ]))
	table(rowSums(clrdf[,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level2"]) ]))
	table(rowSums(clrdf[,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level3"]) ]))
	for (rn in rownames(clrdf)){
	  clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))])] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))]  )] ))
	  clrdf[rn,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level1"])] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level1"])] ))
	  clrdf[rn,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level2"])] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level2"])] ))
	  clrdf[rn,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level3"])] = clr(as.numeric(clrdf[rn,c( colnames(clrdf)[substr(colnames(clrdf),1,8)=="n_level3"])] ))
	}
	save(clrdf,file = paste0( this_OutDir,prefix,"clrdf.RData" ))
	# mean tps vs level_1
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level1.pdf"),ncol(cordf),nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_tps_vs_level1_table.RData"))
	# mean tps vs level_2
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_tps_vs_level2_table.RData"))
	# mean tps vs level_3
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_tps_vs_level3_table.RData"))
	# mean tps vs level_1
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level1.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level1_table.RData"))
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level1_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level1_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level1_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("CLR-transformed ", gsub("n_cs_","",rn)," proportion" ), ylab=paste0("CLR-transformed ", gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level1_scatterplots/ratios_",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cor(radf[,rn], radf[,cn],method="spearman",use="complete.obs"),2),", p-val = ",signif(cor.test(radf[,rn], radf[,cn],method="spearman",use="complete.obs")$p.value,2) )
	    dan.scatterplot( fileName, radf[,rn], radf[,cn], xlab=paste0(gsub("n_cs_","",rn)," proportion"), ylab=paste0(gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	  }
	}
	save(cordf_pvals,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level1_table_pvals.RData"))
	# mean tps vs level_2
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level2_table.RData"))
	# mean tps vs level_3
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level3_table.RData"))

	#### III), cs vs total cancer, tme vs total tme splitted by level1
	radf = codf
	radf$n_allcells = NULL
	columnz = colnames(radf)[substr(colnames(radf),1,4)=="n_cs"]
	radf[,columnz] = sweep(radf[,columnz],1,codf$n_cancercells,FUN="/")
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="epithelial","annd_level_2"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="endothelial","annd_level_2"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="immune","annd_level_2"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="stroma","annd_level_2"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="epithelial","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="endothelial","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="immune","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="stroma","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	radf = radf[,substr(colnames(radf),nchar(colnames(radf))-7,nchar(colnames(radf)))!="_unclear" ]
	radf = radf[,!(colnames(radf) %in% c( "n_level1_glial","n_level2_oligodendrocytes","n_level3_oligodendrocytes" ))]

	this_OutDir = paste0(OutDir,"III_vsTotalCancer_TotalTME_splittedbylevel1/")
	dir.create(this_OutDir)
	# mean tps vs level_2
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_2
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()

	## with clr transformation
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,4)=="n_cs"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_cancercells,FUN="/")
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="epithelial","annd_level_2"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="endothelial","annd_level_2"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="immune","annd_level_2"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level2_",unique(ctdf[ctdf$annd_level_1=="stroma","annd_level_2"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="epithelial","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="endothelial","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="immune","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_1=="stroma","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	clrdf = clrdf[,substr(colnames(clrdf),nchar(colnames(clrdf))-7,nchar(colnames(clrdf)))!="_unclear" ]
	clrdf = clrdf[,!(colnames(clrdf) %in% c( "n_level1_glial","n_level2_oligodendrocytes","n_level3_oligodendrocytes" ))]

	# mean tps vs level_2
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level2_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level2_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("mean ",gsub("tp_","",rn)," score"), ylab=paste0("CLR-transformed ", gsub("n_level2_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level2_scatterplots/ratios_",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cor(radf[,rn], radf[,cn],method="spearman",use="complete.obs"),2),", p-val = ",signif(cor.test(radf[,rn], radf[,cn],method="spearman",use="complete.obs")$p.value,2) )
	    dan.scatterplot( fileName, radf[,rn], radf[,cn], xlab=paste0("mean ",gsub("tp_","",rn)," score"), ylab=paste0("CLR-transformed ", gsub("n_level2_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	  }
	}
	pdf(paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,p.mat = cordf_pvals, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_tps_vs_level2_table.RData"))
	# mean tps vs level_3
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level3_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level3_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("mean ",gsub("tp_","",rn)," score"), ylab=paste0("CLR-transformed ", gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	  }
	}
	pdf(paste0(this_OutDir,prefix,"withPvals_clr_tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,p.mat = cordf_pvals, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_tps_vs_level3_table.RData"))
	# mean cs vs level_2
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level2_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level2_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level2_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("CLR-transformed ", gsub("n_cs_","",rn)," proportion" ), ylab=paste0("CLR-transformed ", gsub("n_level2_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level2_scatterplots/ratios_",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cor(radf[,rn], radf[,cn],method="spearman",use="complete.obs"),2),", p-val = ",signif(cor.test(radf[,rn], radf[,cn],method="spearman",use="complete.obs")$p.value,2) )
	    dan.scatterplot( fileName, radf[,rn], radf[,cn], xlab=paste0(gsub("n_cs_","",rn)," proportion"), ylab=paste0(gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )

	  }
	}
	pdf(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level2.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,p.mat = cordf_pvals, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level2_table.RData"))
	save(cordf_pvals,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level2_table_pvals.RData"))
	# mean cs vs level_3
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("CLR-transformed ", gsub("n_cs_","",rn)," proportion" ), ylab=paste0("CLR-transformed ", gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3_scatterplots/ratios_",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cor(radf[,rn], radf[,cn],method="spearman",use="complete.obs"),2),", p-val = ",signif(cor.test(radf[,rn], radf[,cn],method="spearman",use="complete.obs")$p.value,2) )
	    dan.scatterplot( fileName, radf[,rn], radf[,cn], xlab=paste0(gsub("n_cs_","",rn)," proportion"), ylab=paste0(gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	  }
	}
	pdf(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,p.mat = cordf_pvals, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level3_table.RData"))
	save(cordf_pvals,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level3_table_pvals.RData"))
	# macro_flavours
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_macro_flavours_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,3)=="mf_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_macro_flavours_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("CLR-transformed ", gsub("n_cs_","",rn)," proportion" ), ylab=paste0("mean ",cn," score" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_macro_flavours_scatterplots/ratios_",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cor(radf[,rn], radf[,cn],method="spearman",use="complete.obs"),2),", p-val = ",signif(cor.test(radf[,rn], radf[,cn],method="spearman",use="complete.obs")$p.value,2) )
	    dan.scatterplot( fileName, radf[,rn], radf[,cn], xlab=paste0(gsub("n_cs_","",rn)," proportion"), ylab=paste0("mean ",cn," score" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	  }
	}
	pdf(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_macro_flavours.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,p.mat = cordf_pvals, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_macro_flavours_table.RData"))
	save(cordf_pvals,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_macro_flavours_table_pvals.RData"))

	aa = clrdf[(!is.na(clrdf$mut_KEAP1)) | (!(is.na(clrdf$mut_STK11))),]
	fileName = paste0( OutDir,"STK11_vs_MHCII_boxplots.pdf" )
	dan.boxplots(fileName, x = aa$mut_STK11, y = aa[,"tp_MHC-II"])
	fileName = paste0( OutDir,"KEAP1_vs_MHCII_boxplots.pdf" )
	dan.boxplots(fileName, x = aa$mut_KEAP1, y = aa[,"tp_MHC-II"])
	aa$mut_KEAP1_STK11 = ifelse((aa$mut_KEAP1=="mut") | (aa$mut_STK11=="mut"),"KEAP1/STK11","Double wt" )
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_KEAP1_STK11",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"KEAP1orSTK11_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_KEAP1_STK11, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} # more T cells in double wt

	fileName = paste0( OutDir,"KEAP1orSTK11_vs_Tcells.pdf" )
	dan.boxplots( fileName, x = aa$mut_KEAP1_STK11, y = aa$n_level2_Tcells, xColors = c("black","black"), signifTest="wilcox.test", jitterColors = ifelse(aa$mut_KEAP1_STK11=="Double wt","gray","red") , xlab = "", ylab = "CLR", fileHeight = 4, fileWidth = 5 )

	aa = clrdf[(!is.na(clrdf$mut_EGFR)),]
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_EGFR",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"EGFR_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_EGFR, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} 

	aa = clrdf[(!is.na(clrdf$mut_KRAS)),]
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_KRAS",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"KRAS_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_KRAS, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} 

	aa = clrdf[(!is.na(clrdf$mut_TP53)),]
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_TP53",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"TP53_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_TP53, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} 

	# same but with ratios (and not CLR)
	clrdf = radf
	aa = clrdf[(!is.na(clrdf$mut_KEAP1)) | (!(is.na(clrdf$mut_STK11))),]
	fileName = paste0( OutDir,"ratios_STK11_vs_MHCII_boxplots.pdf" )
	dan.boxplots(fileName, x = aa$mut_STK11, y = aa[,"tp_MHC-II"])
	fileName = paste0( OutDir,"ratios_KEAP1_vs_MHCII_boxplots.pdf" )
	dan.boxplots(fileName, x = aa$mut_KEAP1, y = aa[,"tp_MHC-II"])
	aa$mut_KEAP1_STK11 = ifelse((aa$mut_KEAP1=="mut") | (aa$mut_STK11=="mut"),"KEAP1/STK11","Double wt" )
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_KEAP1_STK11",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"ratios_KEAP1orSTK11_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_KEAP1_STK11, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} # more T cells in double wt

	fileName = paste0( OutDir,"ratios_KEAP1orSTK11_vs_Tcells.pdf" )
	dan.boxplots( fileName, x = aa$mut_KEAP1_STK11, y = aa$n_level2_Tcells, xColors = c("black","black"), signifTest="wilcox.test", jitterColors = ifelse(aa$mut_KEAP1_STK11=="Double wt","gray","red") , xlab = "", ylab = "CLR", fileHeight = 4, fileWidth = 5 )

	aa = clrdf[(!is.na(clrdf$mut_EGFR)),]
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_EGFR",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"ratios_EGFR_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_EGFR, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} 

	aa = clrdf[(!is.na(clrdf$mut_KRAS)),]
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_KRAS",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"ratios_KRAS_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_KRAS, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} 

	aa = clrdf[(!is.na(clrdf$mut_TP53)),]
	library(reshape2)
	for (level in c(2,3)){
	  ab = aa[,c( "mut_TP53",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"ratios_TP53_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_TP53, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} 

	#### IV), cs vs total cancer, c(Bcells,Tcells,macrophages,monocytes,DC) subsets vs total c(Bcells,Tcells,macrophages,monocytes,DC)
	radf = codf
	radf$n_allcells = NULL
	columnz = colnames(radf)[substr(colnames(radf),1,4)=="n_cs"]
	radf[,columnz] = sweep(radf[,columnz],1,codf$n_cancercells,FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="Bcells","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="Tcells","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="macrophages","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="monocytes","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="DC","annd_level_3"]))
	radf[,columnz] = sweep(radf[,columnz],1,as.numeric(rowSums(radf[,columnz])),FUN="/")
	radf = radf[, c( c( colnames(radf)[grepl("n_cs",colnames(radf))] ) ,paste0( "tp_",names(tps) ) ,unique(paste0( "n_level3_",c( ctdf[ctdf$annd_level_2 %in% c( "Bcells","Tcells","macrophages","monocytes","DC" ),"annd_level_3"] )))) ]
	radf = radf[,substr(colnames(radf),nchar(colnames(radf))-7,nchar(colnames(radf)))!="_unclear" ]

	this_OutDir = paste0(OutDir,"IV_vsSubsets/")
	dir.create(this_OutDir)
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,3)=="tp_" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( radf[,substr(colnames(radf),1,4)=="n_cs" ],radf[,substr(colnames(radf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	### with clr
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,4)=="n_cs"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_cancercells,FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="Bcells","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="Tcells","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="macrophages","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="monocytes","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="DC","annd_level_3"]))
	clrdf[,columnz] = sweep(clrdf[,columnz],1,as.numeric(rowSums(clrdf[,columnz])),FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	clrdf = clrdf[, c( c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))] ) ,paste0( "tp_",names(tps) ) ,unique(paste0( "n_level3_",c( ctdf[ctdf$annd_level_2 %in% c( "Bcells","Tcells","macrophages","monocytes","DC" ),"annd_level_3"] )))) ]
	clrdf = clrdf[,substr(colnames(clrdf),nchar(colnames(clrdf))-7,nchar(colnames(clrdf)))!="_unclear" ]
	# mean tps vs level_3
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean cs vs level_3
	dir.create(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3_scatterplots/"))
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	cordf_pvals = cordf
	for (rn in rownames(cordf_pvals) ){
	  for (cn in colnames(cordf_pvals)){
	    cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn],method="spearman",use="complete.obs")$p.value
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3_scatterplots/",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cordf[rn,cn],2),", p-val = ",signif(cordf_pvals[rn,cn],2) )
	    dan.scatterplot( fileName, clrdf[,rn], clrdf[,cn], xlab=paste0("CLR-transformed ", gsub("n_cs_","",rn)," proportion" ), ylab=paste0("CLR-transformed ", gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	    fileName = paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3_scatterplots/ratios_",rn,"_vs_",cn,"_scatterplot.pdf")
	    plotTitle = paste0( "Spearman R = ",signif(cor(radf[,rn], radf[,cn],method="spearman",use="complete.obs"),2),", p-val = ",signif(cor.test(radf[,rn], radf[,cn],method="spearman",use="complete.obs")$p.value,2) )
	    dan.scatterplot( fileName, radf[,rn], radf[,cn], xlab=paste0(gsub("n_cs_","",rn)," proportion"), ylab=paste0(gsub("n_level3_","",cn)," proportion" ), plotTitle=plotTitle, dotSize = 3,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", fileWidth = 6, fileHeight = 6 )
	  }
	}
	pdf(paste0(this_OutDir,prefix,"withPvals_clr_csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,p.mat = cordf_pvals, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', tl.col="black")
	dev.off()
	save(cordf,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level3_table.RData"))
	save(cordf_pvals,file = paste0(this_OutDir,prefix,"clr_csRatios_vs_level3_table_pvals.RData"))



	clrdf$mut_KEAP1 = codf[rownames(clrdf),"mut_KEAP1"]
	clrdf$mut_STK11 = codf[rownames(clrdf),"mut_STK11"]
	aa = clrdf[codf[(!is.na(codf$mut_KEAP1)) | (!(is.na(codf$mut_STK11))),"PST"],]
	aa$mut_KEAP1_STK11 = ifelse((aa$mut_KEAP1=="mut") | (aa$mut_STK11=="mut"),"KEAP1 or STK11-mutant","Double wild-type" )
	library(reshape2)
	for (level in c(3)){
	  ab = aa[,c( "mut_KEAP1_STK11",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"III_KEAP1orSTK11_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_KEAP1_STK11, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} # more T cells in double wt

	#### V), cs vs total cancer, Tcells subsets vs CD4 and CD8 Tcells separately
	### with clr
	this_OutDir = paste0(OutDir,"V_vsCD4CD8/")
	dir.create(this_OutDir)
	clrdf = codf
	clrdf$n_allcells = NULL
	columnz = colnames(clrdf)[substr(colnames(clrdf),1,4)=="n_cs"]
	clrdf[,columnz] = sweep(clrdf[,columnz],1,codf$n_cancercells,FUN="/")
	table(rowSums(clrdf[,columnz ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
	columnz = paste0("n_level3_",unique(ctdf[ctdf$annd_level_2=="Tcells","annd_level_3"]))
	columnz_CD4 = columnz[substr(columnz,1,12)=="n_level3_CD4"]
	columnz_CD8 = columnz[substr(columnz,1,12)=="n_level3_CD8"]
	clrdf[,columnz_CD4] = sweep(clrdf[,columnz_CD4],1,as.numeric(rowSums(clrdf[,columnz_CD4])),FUN="/")
	table(rowSums(clrdf[,columnz_CD4 ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz_CD4] = clr(as.numeric(clrdf[rn,columnz_CD4])) }
	clrdf[,columnz_CD8] = sweep(clrdf[,columnz_CD8],1,as.numeric(rowSums(clrdf[,columnz_CD8])),FUN="/")
	table(rowSums(clrdf[,columnz_CD8 ]))
	for (rn in rownames(clrdf)){ clrdf[rn,columnz_CD8] = clr(as.numeric(clrdf[rn,columnz_CD8])) }

	clrdf = clrdf[, c( c( colnames(clrdf)[grepl("n_cs",colnames(clrdf))] ) ,paste0( "tp_",names(tps) ) ,columnz_CD4,columnz_CD8) ]
	clrdf = clrdf[,substr(colnames(clrdf),nchar(colnames(clrdf))-7,nchar(colnames(clrdf)))!="_unclear" ]
	# mean tps vs level_3
	cordf = cor( clrdf[,substr(colnames(clrdf),1,3)=="tp_" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_tps_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf)/2)
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()
	# mean tps vs level_3
	cordf = cor( clrdf[,substr(colnames(clrdf),1,4)=="n_cs" ],clrdf[,substr(colnames(clrdf),1,9)=="n_level3_" ], method="spearman",use="complete.obs" )
	pdf(paste0(this_OutDir,prefix,"clr_csRatios_vs_level3.pdf"),ncol(cordf)/2,nrow(cordf))
	corrplot(cordf,col=colorz_solid,tl.col="black")
	dev.off()

	clrdf$mut_KEAP1 = codf[rownames(clrdf),"mut_KEAP1"]
	clrdf$mut_STK11 = codf[rownames(clrdf),"mut_STK11"]
	aa = clrdf[codf[(!is.na(codf$mut_KEAP1)) | (!(is.na(codf$mut_STK11))),"PST"],]
	aa$mut_KEAP1_STK11 = ifelse((aa$mut_KEAP1=="mut") | (aa$mut_STK11=="mut"),"KEAP1 or STK11-mutant","Double wild-type" )
	library(reshape2)
	for (level in c(3)){
	  ab = aa[,c( "mut_KEAP1_STK11",colnames(aa)[substr(colnames(aa),1,8)==paste0("n_level",level )] )]
	  ac = melt(ab)
	  fileName = paste0( OutDir,"V_KEAP1orSTK11_vs_level",level,".pdf" )
	  dan.boxplots( fileName, fill = ac$mut_KEAP1_STK11, y = ac$value, x = ac$variable, xlab = "", ylab = "CLR", filllab="Status", fileWidth = ncol(ab)*2 )
	} # more T cells in double wt
}

cs_vs_cytotrace = function( OutDir, scores, prefix, cs_map ){
	load( file=paste0(OutDir,"../../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	scores$cyto = annot[rownames(scores),"cyto"]
	pdf( paste0(OutDir,"cyto_across_cs_level1_boxplots.pdf"),5,5 )
	x = factor(scores[,"cs_level1"],levels=c("cs1","cs2"))
	y = scores$cyto
	plot=dan.boxplots.multipages( x, y, xlab = "", ylab = "cyto score", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 1, xColors = unique(cs_map$colorz_level1), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"cyto_across_cs_level1_Violins.pdf"),5,5 )
	plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 0, xColors = unique(cs_map$colorz_level1), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
	x = factor(scores[,"cs_level1"],levels=c("cs1","cs2"))
	y = scores$cyto
	dan.densityPlot( paste0(OutDir,"cyto_across_cs_level1_Density.pdf"), y, x, groupinglab = "", xlab = paste0("cyto score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = unique(cs_map$colorz_level1), fileWidth = 2.2, fileHeight = 1.5 )
	
	scores$cs_new = NA
	for (rn in rownames(cs_map)){
		scores[scores$cs==rn,"cs_new"] = cs_map[rn,"alias"]
		scores[scores$cs==rn,"cs_colorz"] = cs_map[rn,"colorz"]
		scores[scores$cs==rn,"cs_colorz_level1"] = cs_map[rn,"colorz_level1"]
	}
	pdf( paste0(OutDir,"cyto_across_cs_boxplots.pdf"),7,5 )
	scores$cyto = annot[rownames(scores),"cyto"]
	x = factor(scores[,"cs_new"],levels=cs_map$alias)
	y = scores$cyto
	plot=dan.boxplots.multipages( x, y, xlab = "", ylab = "cyto score", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 1, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"cyto_across_cs_Violins.pdf"),7,5 )
	plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 0, xColors = cs_map$colorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
	dan.densityPlot( paste0(OutDir,"cyto_across_cs_Density.pdf"), y, x, groupinglab = "", xlab = paste0("cyto score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = cs_map$colorz, fileWidth = 2.5, fileHeight = 1.5 )

	for (csn in cs_map$alias){
		annot[scores[scores$cs_new==csn,"CellID"],"Epi_Cancer"] = csn
	}
	dtable(annot$Epi_Cancer)	
	x = factor(annot[,"Epi_Cancer"],levels=c( "AT1","AT2","AT0",cs_map$alias ))
	y = annot$cyto
	pdf( paste0(OutDir,"cyto_across_celltypes_cs_Violins.pdf"),9,5 )
	plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 0, xColors = c("chocolate","dodgerblue4","purple",cs_map$colorz), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
}

tps_vs_cytotrace = function( OutDir, order_tps, scores ){
	load(paste0(OutDir,"saves/final_results.RData"))
	load( file = paste0(OutDir,"../extendedAtlasNormal_annot.RData") )
	cyto = final_results$CytoTRACE
	annot = annot[names(cyto),]	
	annot$cyto = as.numeric(cyto)
	save(annot, file=paste0(OutDir,"annot_extendedAtlasNormal_Cyto.RData"))
	load( file=paste0(OutDir,"annot_extendedAtlasNormal_Cyto.RData"))
	scores = scores[scores$Epi_Cancer=="Cancer",] # only cancer
	scores$cyto = annot[rownames(scores),"cyto"]
	celltype_map = data.frame(row.names=c( "AT1","AT2","AT0","Cancer" ),
		colorz = c("chocolate","dodgerblue2","purple","firebrick") ,stringsAsFactors=F )
	xLevels = c( "Cancer","AT0","AT2","AT1" )
	xColors = celltype_map[xLevels,]
	pdf( paste0(OutDir,"cyto_across_EpiCancer.pdf"),8,6 )
	x = annot[,"Epi_Cancer"]
	y = annot$cyto
	plot=dan.violinplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = "cyto score", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"cyto_across_EpiCancer_BycMinimal.pdf"),14,6 )
	annot$cMinimal = NA
	annot[annot$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	annot[annot$SampleType %in% c( "Primary" ),"cMinimal"] = "Primary"
	annot[annot$SampleType %in% c( "Metastasis" ),"cMinimal"] = "Metastasis"
	x = annot[,"Epi_Cancer"]
	y = annot$cyto
	fill = factor(annot$cMinimal,levels=c( "Normal","Primary","Metastasis" ))
	fillColors = c( "gray44","firebrick", "sienna4" )
	plot=dan.boxplots.multipages( x, y, fill, xlab = "", ylab = "cyto score", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
	pdf( paste0(OutDir,"cyto_across_Datasets_BycMinimal_splittedByEpiCancer.pdf"),8,6 )
	for (ec in c( "Cancer","AT0","AT2","AT1" )){
		ta = annot[annot$Epi_Cancer==ec,]
		x = ta[,"Dataset"]
		y = ta$cyto
		fill = factor(ta$cMinimal,levels=c( "Normal","Primary","Metastasis" ))
		fillColors = c( "gray44","firebrick", "sienna4" )
		plot=dan.boxplots.multipages( x, y, fill, xlab = "", ylab = "cyto score", filllab = "", plotTitle = ec, signifTest = NULL, ylimLeft = 0, ylimRight = 1,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColors, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
	}
	dev.off()

	for (tp in order_tps){
		fileName = paste0(OutDir, "tps_cyto_vs_",tp,".pdf" )
		plotTitle = paste0( "Pearson R = ",signif(cor(scores[,tp],scores[,"cyto"]),2) )
		dan.scatterplot( fileName, scores[,tp], scores[,"cyto"], xlab=paste0(tp," score"), ylab="cyto score", plotTitle=plotTitle, dotSize = 0.1,coord_fixed = FALSE,plotFitLine=T,FitLineMethod="lm", plotFitLine_se=F, fileWidth = 6, fileHeight = 4 )
	}
	rdf = data.frame(row.names = order_tps,cs = order_tps, cor_with_cyto = NA)
	for (n in order_tps ){ rdf[n,"cor_with_cyto"] = cor(scores[,n],scores[,"cyto"]) }
	metric = "cor_with_cyto"
	fileName = paste0(OutDir,"tps_",metric,".pdf")
	rdf = rdf[order(rdf[,metric]),]
	pdf(fileName, 6,4)
	p = ggplot(data=rdf, aes(x=factor(cs,levels = as.character(cs)), y=rdf[,metric])) + geom_bar(stat="identity") + ylab( "Pearson r" ) + xlab( "" ) + theme_classic() + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	p = p + ylim(c( -max(abs(rdf[,metric])),+max(abs(rdf[,metric])) ))
	print(p)
	dev.off()
}

tps_vs_mutations_singlecells = function( OutDir, order_tps, clin, scores ){
	tps_map = data.frame(row.names=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		tps=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		colorz = c( "midnightblue","dodgerblue4","dodgerblue3","dodgerblue2","darkgreen","forestgreen","green3","firebrick4","firebrick3","sienna1","sienna3","sienna4","goldenrod1","goldenrod3","gray25","gray55","gray85","magenta4","magenta3","lightpink3","lightpink" ))

	variants = colnames(clin)[substr(colnames(clin),1,4)=="mut_"]
	for (var in variants){
		dcat(var)
		print(dtable(clin[,var]))
	}
	# unit of analysis: that on which mutations have been called. Which means, Patient X SampleType for all dataset except Wang, Sample x SampleType for wang.
	clin = clin[clin$SampleType!="Normal",]
	clin$pst = paste0(clin$Patient,"_",clin$SampleType)
	clin[clin$Dataset=="wang","pst"] = paste0(clin[clin$Dataset=="wang","Sample"],"_",clin[clin$Dataset=="wang","SampleType"])
	clin = clin[!duplicated(clin$pst),]
	rownames(clin) = clin$pst
	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	scores$pst = paste0(scores$Patient,"_",scores$SampleType)
	scores[scores$Dataset=="wang","pst"] = paste0(scores[scores$Dataset=="wang","Sample"],"_",scores[scores$Dataset=="wang","SampleType"])
	scores = scores[scores$Epi_Cancer=="Cancer",]
	for (tp in names(tps)){
	  dcat(tp)
	  for (pst in rownames(clin)){
	    dcat(pst)
	    clin[pst,paste0("tp_",tp )] = mean(scores[scores$pst==pst,tp])
	  }  
	}
	clin = clin[rowSums(is.na(clin[,paste0("tp_",names(tps))]) )==0,]
	save(clin,file=paste0(OutDir,"clin_pseudobulked.RData"))
	load(file=paste0(OutDir,"clin_pseudobulked.RData"))
	dan.write(clin,file=paste0(OutDir,"mutations_clin_pseudobulk_corrected.txt"),row.names=T)
	## all SampleTypes
	these_muts = c( "KRAS","TP53","EGFR", "KEAP1", "STK11" )
	pdf(paste0(OutDir,"tps_vs_muts_AllSampleTypes.pdf"),22,7)
	for (mut in these_muts){
	  tclin = clin[ (clin[,paste0("mut_",mut )] %in% c( "wt","mut" ) ) %in% c(T) ,c(paste0("mut_",mut ),paste0("tp_",names(tps)) ) ]
	  mtclin = melt(tclin)
	  plot=dan.boxplots.multipages(x=factor(mtclin$variable,levels=paste0( "tp_",order_tps )),y=mtclin$value,fill=mtclin[,paste0("mut_",mut )],signifTest="wilcox",xlab="",ylab="Score",filllab=paste0("mut_",mut ),includeJitters = T)
	  print(plot)
	}
	dev.off()
	### all tests
	# check normality
	for (tp in names(tps)){
		dcat(tp)
		print(shapiro.test(clin[,paste0("tp_",tp)])$p.value)
	} # None of them is normal.
	# wilcoxon test - correct p-values
	wtdf = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			wtdf[tp,mut] = wilcox.test( clin[clin[,paste0("mut_",mut )]=="wt",paste0("tp_",tp)],clin[clin[,paste0("mut_",mut )]=="mut",paste0("tp_",tp)] )$p.value
		}
	}
	wtdf_corrected = data.frame(matrix(p.adjust(unlist(wtdf),method="BH"),nrow=nrow(wtdf),ncol=ncol(wtdf),dimnames=list(rownames(wtdf),colnames(wtdf))))
	dir.create( paste0(OutDir,"AllSampleTypes_wt/") )
	save(wtdf,file = paste0(OutDir,"AllSampleTypes_wt/wtdf.RData"))
	save(wtdf_corrected,file = paste0(OutDir,"AllSampleTypes_wt/wtdf_corrected.RData"))
	mw = which(wtdf_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(wtdf_corrected)[mw[rn,"row"]]
		mut = colnames(wtdf_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"AllSampleTypes_wt/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Wilcoxon test\nBH-adjusted p-val = ",signif(wtdf_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 2, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL )
	}
	# two-ways anova
	twa = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			tc = clin[(!(is.na(clin[,paste0( "mut_",mut )]))) & (!(is.na(clin$Stage_collapsed))),]
			df = data.frame(marker_totest=tc[,paste0("tp_",tp)],stage=tc$Stage_collapsed,mut_status=tc[,paste0( "mut_",mut )], stringsAsFactors=F)
			lm1 = lm(marker_totest ~ stage + mut_status, data = df)
			twa[tp,mut] = anova(lm1)$"Pr(>F)"[2]
		}
	}
	twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
	dir.create( paste0(OutDir,"AllSampleTypes_TwoWaysAnova/") )
	save(twa,file = paste0(OutDir,"AllSampleTypes_TwoWaysAnova/twa.RData"))
	save(twa_corrected,file = paste0(OutDir,"AllSampleTypes_TwoWaysAnova/twa_corrected.RData"))
	mw = which(twa_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(twa_corrected)[mw[rn,"row"]]
		mut = colnames(twa_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"AllSampleTypes_TwoWaysAnova/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Two-ways anova, stage-corrected\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 2, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL )
	}

	## only invasive LUADs
	clin2 = clin[clin$SampleType %in% c( "Primary","Metastasis" ),]
	pdf(paste0(OutDir,"tps_vs_muts_InvasiveLuads.pdf"),22,7)
	for (mut in these_muts){
	  tclin = clin2[ (clin2[,paste0("mut_",mut )] %in% c( "wt","mut" ) ) %in% c(T) ,c(paste0("mut_",mut ),paste0("tp_",names(tps)) ) ]
	  mtclin = melt(tclin)
	  plot=dan.boxplots.multipages(x=factor(mtclin$variable,levels=paste0( "tp_",order_tps )),y=mtclin$value,fill=mtclin[,paste0("mut_",mut )],signifTest="wilcox",xlab="",ylab="Score",filllab=paste0("mut_",mut ),includeJitters = T)
	  print(plot)
	}
	dev.off()
	# wilcoxon test - correct p-values
	wtdf = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			xx = clin2[clin2[,paste0("mut_",mut )]=="wt",paste0("tp_",tp)]
			yy = clin2[clin2[,paste0("mut_",mut )]=="mut",paste0("tp_",tp)]
			wtdf[tp,mut] = wilcox.test( xx[!is.na(xx)] , yy[!is.na(yy)] )$p.value*ifelse(median(xx[!is.na(xx)])<median(yy[!is.na(yy)]),1,-1 )
		}
	}
	dan.write(wtdf,file=paste0(OutDir,"mutations_wilcoxon_signed_pvals_corrected.txt"),row.names=T)
	wtdf = abs(wtdf)
	five_colorz1 = c("#4357AD","#48A9A6","#CCC2B8","#CBA367","#C1666B")
	five_colorz2 = paste0("#",c("1446a0","db3069","f5d547","ebebd3","3c3c3b"))
	five_colorz2 = c("tomato","lavender","cornflowerblue","chartreuse4","gold")
	mm = melt(wtdf)
	mm$tp = rep(rownames(wtdf),ncol(wtdf))
	fileName = paste0( OutDir,"wtdf_summary_barplot.pdf" )
	x = factor(mm$tp,levels=rownames(wtdf))
	y = -log10(mm$value)
	fill = factor(mm$variable,levels=these_muts)
	dan.barplot( fileName, x, y, fill, xlab = "", ylab = "-log10(p-value), Wilcoxon test\nprogram scores in wildtype vs mutant", filllab = "", fillColors = five_colorz1, plotTitle = "", fileWidth = 4, fileHeight = 2.5, y_horizontalLine = 2, color_horizontalLine = "gray11" )

	wtdf_corrected = data.frame(matrix(p.adjust(unlist(wtdf),method="BH"),nrow=nrow(wtdf),ncol=ncol(wtdf),dimnames=list(rownames(wtdf),colnames(wtdf))))
	dir.create( paste0(OutDir,"InvasiveLuads_wt/") )
	save(wtdf,file = paste0(OutDir,"InvasiveLuads_wt/wtdf.RData"))
	save(wtdf_corrected,file = paste0(OutDir,"InvasiveLuads_wt/wtdf_corrected.RData"))
	mw = which(wtdf_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(wtdf_corrected)[mw[rn,"row"]]
		mut = colnames(wtdf_corrected)[mw[rn,"col"]]
		tc = clin2[!is.na(clin2[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray",tps_map[tp,"colorz"] )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"InvasiveLuads_wt/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Wilcoxon test\nBH-adjusted p-val = ",signif(wtdf_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 2, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL )
		dan.boxplots( paste0(OutDir,"InvasiveLuads_wt/",tp,"_vs_",mut,"_mutation_boxplot_nominalPvals.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), signifTest = "wilcox", labelycoo = max(y), xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 2, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL )
	}
	mw = which(wtdf_corrected<0.25,arr.ind=T)
	stage_colors = c( I="#C4DDEE",II="#8BBBDD",III="#4F98CC",IV="#0D76BC",Unknown="#B3B4B5" )
	dir.create(paste0(OutDir,"InvasiveLuads_ColorStage_SmokingShape/"))
	for (rn in 1:nrow(mw)){
		tp = rownames(wtdf_corrected)[mw[rn,"row"]]
		mut = colnames(wtdf_corrected)[mw[rn,"col"]]
		tc = clin2[!is.na(clin2[,paste0( "mut_",mut )]),]
		tc[is.na(tc$Stage_collapsed),"Stage_collapsed"] = "Unknown"
		tc[is.na(tc$Smoking),"Smoking"] = "Unknown"
		tc[!(tc$Smoking %in% c("Unknown","Never")),"Smoking"] = "Ever"
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		tc$colorz = dan.expand_colors(tc$Stage_collapsed,names(stage_colors),stage_colors)
		dan.boxplots( paste0(OutDir,"InvasiveLuads_ColorStage_SmokingShape/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, labelycoo = 1, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 0.5, jitterShape=tc$Smoking, fileWidth = 1.5, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL, legend_position = 'top' )
	}
	# two-ways anova, although not normal
	twa = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			tc = clin2[(!(is.na(clin2[,paste0( "mut_",mut )]))) & (!(is.na(clin2$Stage_collapsed))),]
			df = data.frame(marker_totest=tc[,paste0("tp_",tp)],stage=tc$Stage_collapsed,mut_status=tc[,paste0( "mut_",mut )], stringsAsFactors=F)
			lm1 = lm(marker_totest ~ stage + mut_status, data = df)
			twa[tp,mut] = anova(lm1)$"Pr(>F)"[2]
		}
	}
	twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
	dir.create( paste0(OutDir,"InvasiveLuads_TwoWaysAnova/") )
	save(twa,file = paste0(OutDir,"InvasiveLuads_TwoWaysAnova/twa.RData"))
	save(twa_corrected,file = paste0(OutDir,"InvasiveLuads_TwoWaysAnova/twa_corrected.RData"))
	mw = which(twa_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(twa_corrected)[mw[rn,"row"]]
		mut = colnames(twa_corrected)[mw[rn,"col"]]
		tc = clin2[!is.na(clin2[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"InvasiveLuads_TwoWaysAnova/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Two-ways anova, stage-corrected\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 2, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL )
	}

	# two-ways anova vs smoking
	clin2 = clin2[!is.na(clin2$Smoking),]
	clin2$smoking_collapsed = ifelse(clin2$Smoking=="Never","Never","Ever")
	twa = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			tc = clin2[(!(is.na(clin2[,paste0( "mut_",mut )]))) & (!(is.na(clin2$smoking_collapsed))),]
			df = data.frame(marker_totest=tc[,paste0("tp_",tp)],stage=tc$smoking_collapsed,mut_status=tc[,paste0( "mut_",mut )], stringsAsFactors=F)
			lm1 = lm(marker_totest ~ stage + mut_status, data = df)
			twa[tp,mut] = anova(lm1)$"Pr(>F)"[2]
		}
	}
	twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
	dir.create( paste0(OutDir,"InvasiveLuads_TwoWaysAnova_correctingSmoking/") )
	save(twa,file = paste0(OutDir,"InvasiveLuads_TwoWaysAnova_correctingSmoking/twa.RData"))
	save(twa_corrected,file = paste0(OutDir,"InvasiveLuads_TwoWaysAnova_correctingSmoking/twa_corrected.RData"))
	mw = which(twa_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(twa_corrected)[mw[rn,"row"]]
		mut = colnames(twa_corrected)[mw[rn,"col"]]
		tc = clin2[!is.na(clin2[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"InvasiveLuads_TwoWaysAnova/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Two-ways anova, smoking-corrected\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 2, fileHeight = 1.8, hlines_coo = NULL, hlines_labels = NULL )
	}

	clin2_never = clin2[clin2$smoking_collapsed=="Never",]
	wtdf = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			xx = clin2_never[clin2_never[,paste0("mut_",mut )]=="wt",paste0("tp_",tp)]
			yy = clin2_never[clin2_never[,paste0("mut_",mut )]=="mut",paste0("tp_",tp)]
			if (all(is.na(yy))){ next }
			wtdf[tp,mut] = wilcox.test( xx[!is.na(xx)] , yy[!is.na(yy)] )$p.value*ifelse(median(xx[!is.na(xx)])<median(yy[!is.na(yy)]),1,-1 )
		}
	}
}

tps_vs_mutations_bulk = function( OutDir, order_tps, pairs_to_validate, scoring_method='singscore' ){
	tps_map = data.frame(row.names=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		tps=c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","MHC-I","Interferon","Cell_proliferation","Basal-like","pEMT","EMT","Hypoxia","Metal","OxPhos","Stress_AP1","Stress_HSP","Stress_secreted","Translation_initiation","Unfolded_protein_response","RNA_processing","Unas_emp" ),
		colorz = c( "midnightblue","dodgerblue4","dodgerblue3","dodgerblue2","darkgreen","forestgreen","green3","firebrick4","firebrick3","sienna1","sienna3","sienna4","goldenrod1","goldenrod3","gray25","gray55","gray85","magenta4","magenta3","lightpink3","lightpink" ))
	load(file = paste0(OutDir,"../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	load(file = paste0("processed/scNMF_CV/SCTcenteredRmRBMT__rank5/snp0_genes.RData"))
	tps_universe = snp0_genes

	## TCGA
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_TPM/TCGA_LUAD_ge_TPM_GeneSymbols.RData")
	cn_full = colnames(ge)
	cn = as.character(colnames(ge))
	cn = cn[substr(cn,14,15) %in% c("01")] # only primary
	ge = ge[intersect(rownames(ge),tps_universe),cn]
	cn_short = substr(cn,1,12)
	colnames(ge) = cn_short
	logtpm = log2(ge+1)
	library(singscore)
	rankData = rankGenes(logtpm)
	clin = dan.df(colnames(logtpm),paste0("tp_",names(tps) ) )
	colnames(clin) = gsub( "\\.","-",colnames(clin) )
	for (n in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[n]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[n]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   clin[rownames(scoredf),paste0("tp_",n )] = scoredf$TotalScore
	}
	CommonDataDir = "/mnt/ed2/daniele/Common_Data/"
	load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
	maf_ss = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]

	maf = maf_ss[maf_ss$Hugo_Symbol=="KRAS",c("Patient","Tumor_Sample_Barcode","Variant_Classification","HGVSp_Short","Hugo_Symbol")]
	maf = maf[maf$Variant_Classification=="Missense_Mutation", ]
	clin[intersect(unique(maf_ss$Patient),rownames(clin)),"mut_KRAS"] = "wt"
	clin[intersect(unique(maf$Patient),rownames(clin)),"mut_KRAS"] = "mut"

	maf = maf_ss[maf_ss$Hugo_Symbol=="EGFR",c("Patient","Tumor_Sample_Barcode","Variant_Classification","HGVSc","Hugo_Symbol")]
	maf = maf[maf$Variant_Classification %in% c("Frame_Shift_Del","In_Frame_Del","In_Frame_Ins","Missense_Mutation"), ]
	clin[intersect(unique(maf_ss$Patient),rownames(clin)),"mut_EGFR"] = "wt"
	clin[intersect(unique(maf$Patient),rownames(clin)),"mut_EGFR"] = "mut"

	maf = maf_ss[maf_ss$Hugo_Symbol=="TP53",c("Patient","Tumor_Sample_Barcode","Variant_Classification","HGVSc","Hugo_Symbol")]
	maf = maf[maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site"), ]
	clin[intersect(unique(maf_ss$Patient),rownames(clin)),"mut_TP53"] = "wt"
	clin[intersect(unique(maf$Patient),rownames(clin)),"mut_TP53"] = "mut"

	maf = maf_ss[maf_ss$Hugo_Symbol=="KEAP1",c("Patient","Tumor_Sample_Barcode","Variant_Classification","HGVSc","Hugo_Symbol")]
	maf = maf[maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site"), ]
	clin[intersect(unique(maf_ss$Patient),rownames(clin)),"mut_KEAP1"] = "wt"
	clin[intersect(unique(maf$Patient),rownames(clin)),"mut_KEAP1"] = "mut"

	maf = maf_ss[maf_ss$Hugo_Symbol=="STK11",c("Patient","Tumor_Sample_Barcode","Variant_Classification","HGVSc","Hugo_Symbol")]
	maf = maf[maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site"), ]
	clin[intersect(unique(maf_ss$Patient),rownames(clin)),"mut_STK11"] = "wt"
	clin[intersect(unique(maf$Patient),rownames(clin)),"mut_STK11"] = "mut"

	clin = clin[!is.na(clin$mut_EGFR),]
	these_muts = c( "KRAS","TP53","EGFR", "KEAP1", "STK11" )
	for (mut in these_muts){
	  clin[is.na(clin[,paste0("mut_",mut )]), paste0("mut_",mut )] = "unknown"
	}
	pur = dan.read( paste0("/mnt/ndata/daniele/alfredo_egfr/Data/aran_purities.txt") )
	pur = pur[pur$Sample.ID %in% cn_full, ]
	pur = pur[substr(pur$Sample.ID,14,15)=="01", ]
	pur$Patient = substr(pur$Sample.ID,1,12)
	rownames(pur) = pur$Patient
	clin$Purity = NA
	clin[intersect(rownames(clin),rownames(pur)),"Purity"] = pur[intersect(rownames(clin),rownames(pur)),"CPE"]

	## Association with mutations
	pdf(paste0(OutDir,"TCGA_tps_vs_muts.pdf"),23,7)
	for (mut in these_muts){
	  tclin = clin[!(clin[,paste0("mut_",mut )]=="unknown"),c(paste0("mut_",mut ),paste0("tp_",names(tps)) ) ]
	  mtclin = melt(tclin)
	  plot=dan.boxplots.multipages(x=factor(mtclin$variable,levels=paste0( "tp_",order_tps )),y=mtclin$value,fill=mtclin[,paste0("mut_",mut )],signifTest="wilcox",xlab="",ylab="Score",labelycoo=max(mtclin$value+0.01),filllab=paste0("mut_",mut ),includeJitters = T)
	  print(plot)
	}
	dev.off()
	# wilcoxon test - correct p-values
	wtdf = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			wtdf[tp,mut] = wilcox.test( clin[clin[,paste0("mut_",mut )]=="wt",paste0("tp_",tp)],clin[clin[,paste0("mut_",mut )]=="mut",paste0("tp_",tp)] )$p.value
		}
	}
	wtdf[ pairs_to_validate[rownames(wtdf),colnames(wtdf)]>=0.1 ] = NA
	wtdf_corrected = data.frame(matrix(p.adjust(unlist(wtdf),method="BH"),nrow=nrow(wtdf),ncol=ncol(wtdf),dimnames=list(rownames(wtdf),colnames(wtdf))))
	dir.create( paste0(OutDir,"TCGA_wt/") )
	save(wtdf,file = paste0(OutDir,"TCGA_wt/wtdf.RData"))
	save(wtdf_corrected,file = paste0(OutDir,"TCGA_wt/wtdf_corrected.RData"))
	mw = which(wtdf_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(wtdf_corrected)[mw[rn,"row"]]
		mut = colnames(wtdf_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"TCGA_wt/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Wilcoxon test\nBH-adjusted p-val = ",signif(wtdf_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 3.5, fileWidth = 5, fileHeight = 5, hlines_coo = NULL, hlines_labels = NULL )
	}

	# two-ways anova, correcting for purity
	twa = dan.df(order_tps,these_muts)
	winner = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			tc = clin[(!(is.na(clin[,paste0( "mut_",mut )]))) & (!(is.na(clin$Purity))),]
			# df = data.frame(marker_totest=tc[,paste0("tp_",tp)],Purity=tc$Purity,mut_status=tc[,paste0( "mut_",mut )], stringsAsFactors=F)
			# lm1 = lm(marker_totest ~ Purity + mut_status, data = df)
			# winner[tp,mut] = ifelse(summary(lm1)$coefficients["mut_statuswt","Estimate"]>0,"wt","mut" )
			# twa[tp,mut] = anova(lm1)$"Pr(>F)"[2]

			# non-parametric using Quade method (see https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/Quade )
			tc$purityRank = rank(tc$Purity)
			tc$tpRank = rank(tc[,paste0("tp_",tp)])
			tc$am = tc[,paste0( "mut_",mut )]
			lm2 = lm(tpRank ~ purityRank, data = tc)
			raov = aov(lm2$residuals ~ am, tc)
			twa[tp,mut] = summary(raov)[[1]][["Pr(>F)"]][[1]]
			winner[tp,mut] = ifelse(median(lm2$residuals[rownames(tc[tc$am=="wt",])])>median(lm2$residuals[rownames(tc[tc$am=="mut",])]),"wt","mut" )
		}
	}
	twa[ pairs_to_validate[rownames(twa),colnames(twa)]>=0.1 ] = NA
	twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
	dir.create( paste0(OutDir,"TCGA_Quade/") )
	save(twa,file = paste0(OutDir,"TCGA_Quade/twa.RData"))
	save(twa_corrected,file = paste0(OutDir,"TCGA_Quade/twa_corrected.RData"))
	mw = which(twa_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(twa_corrected)[mw[rn,"row"]]
		mut = colnames(twa_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		tc = tc[order(-tc[,paste0( "tp_",tp )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","black" ) #c( "gray",tps_map[tp,"colorz"] )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"TCGA_Quade/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Quade ANCOVA, purity-corrected\nWinner: ",winner[tp,mut],", BH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 3.5, fileWidth = 5, fileHeight = 5, hlines_coo = NULL, hlines_labels = NULL )
		dan.scatterplot( paste0(OutDir,"TCGA_Quade/",tp,"_vs_",mut,"_mutation_scatterplot.pdf"), tc$Purity, y, fill = factor(x,levels=xLevels), xlab = "Tumor purity", ylab = tp, filllab = mut, fillColors = xColors, plotTitle = paste0("Pearson's R = ", signif(cor(tc$Purity,y),2),", p-value = ",signif(cor.test(tc$Purity,y)$p.value,2) ), dotLabels = NULL, dotSize = 2, plotBisector = NULL, plotFitLine = T, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 6, fileHeight = 5 )
		dan.barplot( paste0(OutDir,"TCGA_Quade/",tp,"_vs_",mut,"_mutation_barplot.pdf"), rank(y), y, fill = factor(x,levels=xLevels), sd = NULL, xlab = "", ylab = paste0(tp," score"), filllab = mut, fillColors = xColors, labelBarsBottom = NULL, textOnTop = NULL, plotTitle = paste0("Quade ANCOVA, purity-corrected\nWinner: ",winner[tp,mut],", BH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), fileWidth = 7, fileHeight = 5 )
		dan.densityPlot( paste0(OutDir,"TCGA_Quade/",tp,"_vs_",mut,"_mutation_DensityPlot.pdf"), y, factor(x,levels=xLevels), groupinglab = mut, xlab = paste0(tp," score"), ylab = "Density", show_medians = T, plotTitle = paste0("Purity-corrected Quade p = ",signif(twa[tp,mut],2)),xlimLeft = NULL, xlimRight = NULL, groupingColors = xColors, fileWidth = 1.9, fileHeight = 1 )
	}

	# twa but stage-wise
	a=load(file=paste0(OutDir,"../tps_across_ClinicalCharacteristics/Bulk_Tps_scores_TCGA.RData"))
	clin$stage_collapsed = Clin[rownames(clin),"stage_collapsed"]
	clin = clin[!is.na(clin$stage_collapsed),]
	for (stage in unique(clin$stage_collapsed) ){
		# two-ways anova, correcting for purity
		twa = dan.df(order_tps,these_muts)
		winner = dan.df(order_tps,these_muts)
		for (tp in order_tps){
			for (mut in these_muts){
				tc = clin[((!(is.na(clin[,paste0( "mut_",mut )]))) & (!(is.na(clin$Purity)))) & (clin$stage_collapsed==stage),]
				# non-parametric using Quade method (see https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/Quade )
				tc$purityRank = rank(tc$Purity)
				tc$tpRank = rank(tc[,paste0("tp_",tp)])
				tc$am = tc[,paste0( "mut_",mut )]
				lm2 = lm(tpRank ~ purityRank, data = tc)
				raov = aov(lm2$residuals ~ am, tc)
				twa[tp,mut] = summary(raov)[[1]][["Pr(>F)"]][[1]]
				winner[tp,mut] = ifelse(median(lm2$residuals[rownames(tc[tc$am=="wt",])])>median(lm2$residuals[rownames(tc[tc$am=="mut",])]),"wt","mut" )
			}
		}
		twa[ pairs_to_validate[rownames(twa),colnames(twa)]>=0.1 ] = NA
		twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
		save(twa,file = paste0(OutDir,"TCGA_Quade/twa_stage",stage,".RData"))
		save(twa_corrected,file = paste0(OutDir,"TCGA_Quade/twa_corrected_stage",stage,".RData"))
	}
	

	## ChenEAS
	load("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_normalized_GeneSymbols.RData")
	maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
	this_maf = maf[(maf$Variant_Classification %in% c("De_novo_Start_InFrame","De_novo_Start_OutOfFrame","frame_shift_del","frame_shift_ins",
	    "in_frame_del","in_frame_ins","nonsense_mutation","missense","nonstop_mutation","Splice_Site","Start_Codon_Del","Start_Codon_SNP","Stop_Codon_Del")),]
	cn = unique(maf$Tumor_Sample_Barcode)
	ge = ge[intersect(rownames(ge),tps_universe),intersect(cn,colnames(ge))]
	logtpm = ge[intersect(rownames(ge),tps_universe),]
	logtpm = logtpm[rowSums(is.na(logtpm))==0, ]
	library(singscore)
	rankData = rankGenes(logtpm)
	clin = dan.df(colnames(logtpm),paste0("tp_",names(tps) ) )
	colnames(clin) = gsub( "\\.","-",colnames(clin) )
	for (n in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[n]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[n]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   clin[rownames(scoredf),paste0("tp_",n )] = scoredf$TotalScore
	}

	clin$mut_KRAS = "wt"
	this_maf2 = this_maf[this_maf$Hugo_Symbol=="KRAS",] # all missense anyway
	mutants = unique(this_maf2[,"Tumor_Sample_Barcode"])
	clin[intersect(mutants,rownames(clin)),"mut_KRAS"] = "mut"

	clin$mut_EGFR = "wt"
	this_maf2 = this_maf[this_maf$Hugo_Symbol=="EGFR",] # all good
	mutants = unique(this_maf2[,"Tumor_Sample_Barcode"])
	clin[intersect(mutants,rownames(clin)),"mut_EGFR"] = "mut"

	clin$mut_TP53 = "wt"
	# this_maf2 = this_maf[(this_maf$Hugo_Symbol=="TP53") & (this_maf$Variant_Classification=="missense"),] 
	# mutants = unique(this_maf2[,"Tumor_Sample_Barcode"])
	# clin[intersect(mutants,rownames(clin)),"mut_TP53"] = "mut_missense"
	this_maf2 = this_maf[(this_maf$Hugo_Symbol=="TP53") & (this_maf$Variant_Classification!="missense"),] 
	mutants = unique(this_maf2[,"Tumor_Sample_Barcode"])
	clin[intersect(mutants,rownames(clin)),"mut_TP53"] = "mut"

	clin$mut_KEAP1 = "wt"
	this_maf2 = this_maf[ (this_maf$Variant_Classification!="missense") & (this_maf$Hugo_Symbol=="KEAP1"),] # all good
	mutants = unique(this_maf2[,"Tumor_Sample_Barcode"])
	clin[intersect(mutants,rownames(clin)),"mut_KEAP1"] = "mut"

	clin$mut_STK11 = "wt"
	this_maf2 = this_maf[ (this_maf$Variant_Classification!="missense") & (this_maf$Hugo_Symbol=="STK11"),] # all good
	mutants = unique(this_maf2[,"Tumor_Sample_Barcode"])
	clin[intersect(mutants,rownames(clin)),"mut_STK11"] = "mut"

	these_muts = c( "KRAS","TP53","EGFR", "KEAP1", "STK11" )
	for (mut in these_muts){
	  clin[is.na(clin[,paste0("mut_",mut )]), paste0("mut_",mut )] = "unknown"
	}
	## Association with mutations
	pdf(paste0(OutDir,"ChenEAS_tps_vs_muts.pdf"),23,7)
	for (mut in these_muts){
	  tclin = clin[!(clin[,paste0("mut_",mut )]=="unknown"),c(paste0("mut_",mut ),paste0("tp_",names(tps)) ) ]
	  mtclin = melt(tclin)
	  plot=dan.boxplots.multipages(x=factor(mtclin$variable,levels=paste0( "tp_",order_tps )),y=mtclin$value,fill=mtclin[,paste0("mut_",mut )],signifTest="wilcox",xlab="",ylab="Score",filllab=paste0("mut_",mut ),includeJitters = T)
	  print(plot)
	}
	dev.off()
	# wilcoxon test - correct p-values
	wtdf = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			wtdf[tp,mut] = wilcox.test( clin[clin[,paste0("mut_",mut )]=="wt",paste0("tp_",tp)],clin[clin[,paste0("mut_",mut )]=="mut",paste0("tp_",tp)] )$p.value
		}
	}
	wtdf[ pairs_to_validate[rownames(wtdf),colnames(wtdf)]>=0.1 ] = NA
	wtdf_corrected = data.frame(matrix(p.adjust(unlist(wtdf)[!is.na(unlist(wtdf))],method="BH"),nrow=nrow(wtdf),ncol=ncol(wtdf),dimnames=list(rownames(wtdf),colnames(wtdf))))
	dir.create( paste0(OutDir,"ChenEAS_wt/") )
	save(wtdf,file = paste0(OutDir,"ChenEAS_wt/wtdf.RData"))
	save(wtdf_corrected,file = paste0(OutDir,"ChenEAS_wt/wtdf_corrected.RData"))
	mw = which(wtdf_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(wtdf_corrected)[mw[rn,"row"]]
		mut = colnames(wtdf_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"ChenEAS_wt/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Wilcoxon test\nBH-adjusted p-val = ",signif(wtdf_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 3.5, fileWidth = 5, fileHeight = 5, hlines_coo = NULL, hlines_labels = NULL )
	}

	# adding purity
	pur = dan.read( paste0("/mnt/ndata/daniele/alfredo_egfr/Data/Chen_GIS031_clinical_data.tsv") )
	rownames(pur) = pur$Patient.ID
	clin$Purity = NA
	clin[intersect(rownames(clin),rownames(pur)),"Purity"] = pur[intersect(rownames(clin),rownames(pur)),"Purity"]
	# two-ways anova, correcting for purity
	twa = dan.df(order_tps,these_muts)
	winner = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			tc = clin[(!(is.na(clin[,paste0( "mut_",mut )]))) & (!(is.na(clin$Purity))),]
			# df = data.frame(marker_totest=tc[,paste0("tp_",tp)],Purity=tc$Purity,mut_status=tc[,paste0( "mut_",mut )], stringsAsFactors=F)
			# lm1 = lm(marker_totest ~ Purity + mut_status, data = df)
			# winner[tp,mut] = ifelse(summary(lm1)$coefficients["mut_statuswt","Estimate"]>0,"wt","mut" )
			# twa[tp,mut] = anova(lm1)$"Pr(>F)"[2]

			# non-parametric using Quade method (see https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/Quade )
			tc$purityRank = rank(tc$Purity)
			tc$tpRank = rank(tc[,paste0("tp_",tp)])
			tc$am = tc[,paste0( "mut_",mut )]
			lm2 = lm(tpRank ~ purityRank, data = tc)
			raov = aov(lm2$residuals ~ am, tc)
			twa[tp,mut] = summary(raov)[[1]][["Pr(>F)"]][[1]]
			winner[tp,mut] = ifelse(median(lm2$residuals[rownames(tc[tc$am=="wt",])])>median(lm2$residuals[rownames(tc[tc$am=="mut",])]),"wt","mut" )
		}
	}
	twa[ pairs_to_validate[rownames(twa),colnames(twa)]>=0.1 ] = NA
	twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
	dir.create( paste0(OutDir,"ChenEAS_Quade/") )
	save(twa,file = paste0(OutDir,"ChenEAS_Quade/twa.RData"))
	save(twa_corrected,file = paste0(OutDir,"ChenEAS_Quade/twa_corrected.RData"))
	mw = which(twa_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(twa_corrected)[mw[rn,"row"]]
		mut = colnames(twa_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray",tps_map[tp,"colorz"] )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"ChenEAS_Quade/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Quade ANCOVA, purity-corrected\nWinner: ",winner[tp,mut],"\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 3.5, fileWidth = 5, fileHeight = 5, hlines_coo = NULL, hlines_labels = NULL )
		dan.scatterplot( paste0(OutDir,"ChenEAS_Quade/",tp,"_vs_",mut,"_mutation_scatterplot.pdf"), tc$Purity, y, fill = factor(x,levels=xLevels), xlab = "Tumor purity", ylab = tp, filllab = "", fillColors = xColors, plotTitle = paste0("Pearson's R = ", signif(cor(tc$Purity,y),2),", p-value = ",signif(cor.test(tc$Purity,y)$p.value,2) ), dotLabels = NULL, dotSize = 2, plotBisector = NULL, plotFitLine = T, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 6, fileHeight = 5 )
		dan.barplot( paste0(OutDir,"ChenEAS_Quade/",tp,"_vs_",mut,"_mutation_barplot.pdf"), rank(y), y, fill = factor(x,levels=xLevels), sd = NULL, xlab = "", ylab = paste0(tp," score"), filllab = "", fillColors = xColors, labelBarsBottom = NULL, textOnTop = NULL, plotTitle = paste0("Quade ANCOVA, purity-corrected\nWinner: ",winner[tp,mut],"\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), fileWidth = 7, fileHeight = 5 )
	}

	## TRACERx
	# Reading TRACERx421 data. A lot of information: metastasis-seeding regions, mutations, patterns, ...
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
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	logtpm = log2(ge[intersect(rownames(ge),tps_universe),]+0.0001)
	logtpm = logtpm[,substr(colnames(logtpm),1,8) %in% rownames(Clin)]
	logtpm = logtpm[,!grepl("_N",colnames(logtpm))]
	library(singscore)
	rankData = rankGenes(logtpm)
	clin = dan.df(colnames(logtpm),paste0("tp_",names(tps) ) )
	colnames(clin) = gsub( "\\.","-",colnames(clin) )
	for (n in names(tps))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[n]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[n]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   clin[rownames(scoredf),paste0("tp_",n )] = scoredf$TotalScore
	}

	maf = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	mafr = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/figurecode/data/20221123_TRACERx421_mutation_table_region.fst"))
	mafr$RegionID = gsub(":","_",gsub("\\.","-",mafr$RegionID))
	mafr = mafr[mafr$RegionID %in% colnames(logtpm),]
	maf = maf[maf$mutation_id %in% mafr$mutation_id,c( "mutation_id","Hugo_Symbol","func","exonic.func","AAChange","DriverMut" )]
	maf_ss = merge(mafr,maf,by = "mutation_id",all.x=T)
	# maybe missing metastases?
	mafmet = load(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/metsFigures/data/patientMutTable.20220726.rda"))
	mafmet = mutTable[mutTable$patient_id %in% substr(colnames(logtpm),1,8),]
	# mafmet_regions_processed = c()
	# index = 1
	# for (rn in rownames(mafmet)){
	# 	aa = unlist(strsplit(mafmet[rn,"Is.present"],split=";"))
	# 	aa = gsub(":TRUE","",aa)
	# 	aa = gsub(":FALSE","",aa)
	# 	aa = gsub( "\\.","-",aa )
	# 	aa = paste0(mafmet[rn,"patient_id"],"_",aa)
	# 	mafmet_regions_processed = c(mafmet_regions_processed,aa)
	# 	index = index+1
	# 	if (index==1000) { dcat(index) }
	# 	if (index==10000) { dcat(index) }
	# }
	# mafmet_regions_processed = unique(mafmet_regions_processed)
	# save(mafmet_regions_processed,file = paste0( "/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/metsFigures/data/dan.mafmet_regions_processed.RData" ))
	load(file = paste0( "/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/metsFigures/data/dan.mafmet_regions_processed.RData" ))
	clin[intersect(rownames(clin),unique(mafmet_regions_processed)),"mut_KRAS"] = "wt" # for these regions, we are sure that they were WES/WGS. Then, we'll see which ones are mut
	clin[intersect(rownames(clin),unique(mafmet_regions_processed)),"mut_EGFR"] = "wt"
	clin[intersect(rownames(clin),unique(mafmet_regions_processed)),"mut_TP53"] = "wt"
	clin[intersect(rownames(clin),unique(mafmet_regions_processed)),"mut_STK11"] = "wt"
	clin[intersect(rownames(clin),unique(mafmet_regions_processed)),"mut_KEAP1"] = "wt"
	mafmet = mafmet[(mafmet$Hugo_Symbol %in% c( "KRAS","EGFR","TP53","STK11","KEAP1" )) & (!is.na(mafmet$exonic.func)),]

	maf = maf_ss[(maf_ss$Hugo_Symbol=="KRAS") & (!is.na(maf_ss$exonic.func)),c("RegionID","exonic.func","AAChange","DriverMut")]
	dtable(maf$exonic.func,maf$AAChange)
	maf = maf[maf$DriverMut, ] # all hotspots
	clin[intersect(unique(maf_ss$RegionID),rownames(clin)),"mut_KRAS"] = "wt"
	clin[intersect(unique(maf$RegionID),rownames(clin)),"mut_KRAS"] = "mut"
	tm = mafmet[(mafmet$Hugo_Symbol=="KRAS"),]
	dtable(tm$exonic.func,tm$AAChange)
	tm = tm[tm$DriverMut,]
	tm_regions = c()
	for (rn in rownames(tm)){
		aa = unlist(strsplit(tm[rn,"Is.present"],split=";"))
		aa = aa[grepl(":TRUE",aa)]
		aa = gsub(":TRUE","",aa)
		aa = gsub( "\\.","-",aa )
		aa = paste0(tm[rn,"patient_id"],"_",aa)
		tm_regions = c(tm_regions,aa)
	}
	clin[intersect( rownames(clin),unique(tm_regions) ),"mut_KRAS"] = "mut" 

	maf = maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & (!is.na(maf_ss$exonic.func)),c("RegionID","exonic.func","AAChange","DriverMut")]
	dtable(maf$exonic.func,maf$AAChange)
	maf = maf[maf$DriverMut, ] # all hotspots
	clin[intersect(unique(maf_ss$RegionID),rownames(clin)),"mut_EGFR"] = "wt"
	clin[intersect(unique(maf$RegionID),rownames(clin)),"mut_EGFR"] = "mut"
	tm = mafmet[(mafmet$Hugo_Symbol=="EGFR"),]
	dtable(tm$exonic.func,tm$AAChange)
	tm = tm[tm$DriverMut,]
	tm_regions = c()
	for (rn in rownames(tm)){
		aa = unlist(strsplit(tm[rn,"Is.present"],split=";"))
		aa = aa[grepl(":TRUE",aa)]
		aa = gsub(":TRUE","",aa)
		aa = gsub( "\\.","-",aa )
		aa = paste0(tm[rn,"patient_id"],"_",aa)
		tm_regions = c(tm_regions,aa)
	}
	clin[intersect( rownames(clin),unique(tm_regions) ),"mut_EGFR"] = "mut" 

	maf = maf_ss[(maf_ss$Hugo_Symbol=="TP53") & (!is.na(maf_ss$exonic.func)),c("RegionID","exonic.func","AAChange","DriverMut")]
	dtable(maf$exonic.func,maf$DriverMut)
	maf = maf[maf$DriverMut, ] # all hotspots
	clin[intersect(unique(maf_ss$RegionID),rownames(clin)),"mut_TP53"] = "wt"
	clin[intersect(unique(maf$RegionID),rownames(clin)),"mut_TP53"] = "mut"
	unassigned = rownames(clin[is.na(clin$mut_TP53),])
	tm = mafmet[(mafmet$Hugo_Symbol=="TP53"),]
	dtable(tm$exonic.func,tm$AAChange)
	tm = tm[tm$DriverMut,]
	tm_regions = c()
	for (rn in rownames(tm)){
		aa = unlist(strsplit(tm[rn,"Is.present"],split=";"))
		aa = aa[grepl(":TRUE",aa)]
		aa = gsub(":TRUE","",aa)
		aa = gsub( "\\.","-",aa )
		aa = paste0(tm[rn,"patient_id"],"_",aa)
		tm_regions = c(tm_regions,aa)
	}
	clin[intersect( rownames(clin),unique(tm_regions) ),"mut_TP53"] = "mut" 

	maf = maf_ss[(maf_ss$Hugo_Symbol=="KEAP1") & (!is.na(maf_ss$exonic.func)),c("RegionID","exonic.func","AAChange","DriverMut")]
	dtable(maf$exonic.func,maf$DriverMut)
	maf = maf[maf$DriverMut, ] # all hotspots
	clin[intersect(unique(maf_ss$RegionID),rownames(clin)),"mut_KEAP1"] = "wt"
	clin[intersect(unique(maf$RegionID),rownames(clin)),"mut_KEAP1"] = "mut"
	unassigned = rownames(clin[is.na(clin$mut_KEAP1),])
	tm = mafmet[(mafmet$Hugo_Symbol=="KEAP1"),]
	dtable(tm$exonic.func,tm$AAChange)
	tm = tm[tm$DriverMut,]
	tm_regions = c()
	for (rn in rownames(tm)){
		aa = unlist(strsplit(tm[rn,"Is.present"],split=";"))
		aa = aa[grepl(":TRUE",aa)]
		aa = gsub(":TRUE","",aa)
		aa = gsub( "\\.","-",aa )
		aa = paste0(tm[rn,"patient_id"],"_",aa)
		tm_regions = c(tm_regions,aa)
	}
	clin[intersect( rownames(clin),unique(tm_regions) ),"mut_KEAP1"] = "mut" 

	maf = maf_ss[(maf_ss$Hugo_Symbol=="STK11") & (!is.na(maf_ss$exonic.func)),c("RegionID","exonic.func","AAChange","DriverMut")]
	dtable(maf$exonic.func,maf$DriverMut)
	maf = maf[maf$DriverMut, ] # all hotspots
	clin[intersect(unique(maf_ss$RegionID),rownames(clin)),"mut_STK11"] = "wt"
	clin[intersect(unique(maf$RegionID),rownames(clin)),"mut_STK11"] = "mut"
	unassigned = rownames(clin[is.na(clin$mut_STK11),])
	tm = mafmet[(mafmet$Hugo_Symbol=="STK11"),]
	dtable(tm$exonic.func,tm$AAChange)
	tm = tm[tm$DriverMut,]
	tm_regions = c()
	for (rn in rownames(tm)){
		aa = unlist(strsplit(tm[rn,"Is.present"],split=";"))
		aa = aa[grepl(":TRUE",aa)]
		aa = gsub(":TRUE","",aa)
		aa = gsub( "\\.","-",aa )
		aa = paste0(tm[rn,"patient_id"],"_",aa)
		tm_regions = c(tm_regions,aa)
	}
	clin[intersect( rownames(clin),unique(tm_regions) ),"mut_STK11"] = "mut" 

	clin = clin[!is.na(clin$mut_EGFR),]
	these_muts = c( "KRAS","TP53","EGFR", "KEAP1", "STK11" )
	for (mut in these_muts){
	  clin[is.na(clin[,paste0("mut_",mut )]), paste0("mut_",mut )] = "unknown"
	}

	## Association with mutations
	pdf(paste0(OutDir,"TRACERx_tps_vs_muts.pdf"),23,7)
	for (mut in these_muts){
	  tclin = clin[!(clin[,paste0("mut_",mut )]=="unknown"),c(paste0("mut_",mut ),paste0("tp_",names(tps)) ) ]
	  mtclin = melt(tclin)
	  plot=dan.boxplots.multipages(x=factor(mtclin$variable,levels=paste0( "tp_",order_tps )),y=mtclin$value,fill=mtclin[,paste0("mut_",mut )],signifTest="wilcox",xlab="",ylab="Score",labelycoo=max(mtclin$value+0.01),filllab=paste0("mut_",mut ),includeJitters = T)
	  print(plot)
	}
	dev.off()
	# wilcoxon test - correct p-values
	wtdf = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			wtdf[tp,mut] = wilcox.test( clin[clin[,paste0("mut_",mut )]=="wt",paste0("tp_",tp)],clin[clin[,paste0("mut_",mut )]=="mut",paste0("tp_",tp)] )$p.value
		}
	}
	wtdf[ pairs_to_validate[rownames(wtdf),colnames(wtdf)]>=0.1 ] = NA
	wtdf_corrected = data.frame(matrix(p.adjust(unlist(wtdf)[!is.na(unlist(wtdf))],method="BH"),nrow=nrow(wtdf),ncol=ncol(wtdf),dimnames=list(rownames(wtdf),colnames(wtdf))))
	dir.create( paste0(OutDir,"TRACERx_wt/") )
	save(wtdf,file = paste0(OutDir,"TRACERx_wt/wtdf.RData"))
	save(wtdf_corrected,file = paste0(OutDir,"TRACERx_wt/wtdf_corrected.RData"))
	mw = which(wtdf_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(wtdf_corrected)[mw[rn,"row"]]
		mut = colnames(wtdf_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","firebrick" )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"TRACERx_wt/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Wilcoxon test\nBH-adjusted p-val = ",signif(wtdf_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 3.5, fileWidth = 5, fileHeight = 5, hlines_coo = NULL, hlines_labels = NULL )
	}

	# adding purity
	pur = dan.read(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/20221109_TRACERx421_manual_qc_sheet.tsv")) # contains loooots of metrics - can be used
	rownames(pur) = gsub("\\.","-",pur$Sample)
	clin$Purity = NA
	clin[intersect(rownames(clin),rownames(pur)),"Purity"] = pur[intersect(rownames(clin),rownames(pur)),"ASCAT.purity"]
	# two-ways anova, correcting for purity
	twa = dan.df(order_tps,these_muts)
	winner = dan.df(order_tps,these_muts)
	for (tp in order_tps){
		for (mut in these_muts){
			tc = clin[(!(is.na(clin[,paste0( "mut_",mut )]))) & (!(is.na(clin$Purity))),]
			# df = data.frame(marker_totest=tc[,paste0("tp_",tp)],Purity=tc$Purity,mut_status=tc[,paste0( "mut_",mut )], stringsAsFactors=F)
			# lm1 = lm(marker_totest ~ Purity + mut_status, data = df)
			# winner[tp,mut] = ifelse(summary(lm1)$coefficients["mut_statuswt","Estimate"]>0,"wt","mut" )
			# twa[tp,mut] = anova(lm1)$"Pr(>F)"[2]

			# non-parametric using Quade method (see https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/Quade )
			tc$purityRank = rank(tc$Purity)
			tc$tpRank = rank(tc[,paste0("tp_",tp)])
			tc$am = tc[,paste0( "mut_",mut )]
			lm2 = lm(tpRank ~ purityRank, data = tc)
			raov = aov(lm2$residuals ~ am, tc)
			twa[tp,mut] = summary(raov)[[1]][["Pr(>F)"]][[1]]
			winner[tp,mut] = ifelse(median(lm2$residuals[rownames(tc[tc$am=="wt",])])>median(lm2$residuals[rownames(tc[tc$am=="mut",])]),"wt","mut" )
		}
	}
	twa[ pairs_to_validate[rownames(twa),colnames(twa)]>=0.1 ] = NA
	twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
	dir.create( paste0(OutDir,"TRACERx_Quade/") )
	save(twa,file = paste0(OutDir,"TRACERx_Quade/twa.RData"))
	save(twa_corrected,file = paste0(OutDir,"TRACERx_Quade/twa_corrected.RData"))
	mw = which(twa_corrected<0.1,arr.ind=T)
	for (rn in 1:nrow(mw)){
		tp = rownames(twa_corrected)[mw[rn,"row"]]
		mut = colnames(twa_corrected)[mw[rn,"col"]]
		tc = clin[!is.na(clin[,paste0( "mut_",mut )]),]
		x = tc[,paste0( "mut_",mut )]
		y = tc[,paste0( "tp_",tp )]
		xLevels = c( "wt","mut" )
		xColors = c( "gray","black" )#c( "gray",tps_map[tp,"colorz"] )
		tc$colorz = dan.expand_colors(tc[,paste0( "mut_",mut )],xLevels,xColors)
		dan.boxplots( paste0(OutDir,"TRACERx_Quade/",tp,"_vs_",mut,"_mutation_boxplot.pdf"), factor(x,levels=xLevels), y, xlab=mut,ylab = paste0(tp," score"), plotTitle = paste0("Quade ANCOVA, purity-corrected\nWinner: ",winner[tp,mut],"\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), signifTest = NULL, labelycoo = 1, xColors = xColors, jitterColors = tc$colorz, labelJitteredPoints = NULL, jitterDotSize = 3.5, fileWidth = 5, fileHeight = 5, hlines_coo = NULL, hlines_labels = NULL )
		dan.scatterplot( paste0(OutDir,"TRACERx_Quade/",tp,"_vs_",mut,"_mutation_scatterplot.pdf"), tc$Purity, y, fill = factor(x,levels=xLevels), xlab = "Tumor purity", ylab = tp, filllab = "", fillColors = xColors, plotTitle = paste0("Pearson's R = ", signif(cor(tc$Purity,y),2),", p-value = ",signif(cor.test(tc$Purity,y)$p.value,2) ), dotLabels = NULL, dotSize = 2, plotBisector = NULL, plotFitLine = T, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 6, fileHeight = 5 )
		dan.barplot( paste0(OutDir,"TRACERx_Quade/",tp,"_vs_",mut,"_mutation_barplot.pdf"), rank(y), y, fill = factor(x,levels=xLevels), sd = NULL, xlab = "", ylab = paste0(tp," score"), filllab = "", fillColors = xColors, labelBarsBottom = NULL, textOnTop = NULL, plotTitle = paste0("Quade ANCOVA, purity-corrected\nWinner: ",winner[tp,mut],"\nBH-adjusted p-val = ",signif(twa_corrected[tp,mut],2)), fileWidth = 7, fileHeight = 5 )
		dan.densityPlot( paste0(OutDir,"TRACERx_Quade/",tp,"_vs_",mut,"_mutation_DensityPlot.pdf"), y, factor(x,levels=xLevels), groupinglab = mut, xlab = paste0(tp," score"), ylab = "Density", show_medians = T, plotTitle = paste0("Purity-corrected Quade p = ",signif(twa[tp,mut],2)),xlimLeft = NULL, xlimRight = NULL, groupingColors = xColors, fileWidth = 1.9, fileHeight = 1 )
	}

	# twa but stage-wise
	a=load(file=paste0(OutDir,"../tps_across_ClinicalCharacteristics/Bulk_Tps_scores_TRACERx.RData"))
	for (rn in rownames(clin)){
		tpat = substr(rn,1,8)
		clin[rn,"stage_collapsed"] = gsub("B","",gsub("A","",Clin[tpat,"pathologyTNM"]))
	}
	clin = clin[!is.na(clin$stage_collapsed),]
	for (stage in unique(clin$stage_collapsed) ){
		# two-ways anova, correcting for purity
		twa = dan.df(order_tps,these_muts)
		winner = dan.df(order_tps,these_muts)
		for (tp in order_tps){
			for (mut in these_muts){
				tc = clin[((!(is.na(clin[,paste0( "mut_",mut )]))) & (!(is.na(clin$Purity)))) & (clin$stage_collapsed==stage),]
				# non-parametric using Quade method (see https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/Quade )
				tc$purityRank = rank(tc$Purity)
				tc$tpRank = rank(tc[,paste0("tp_",tp)])
				tc$am = tc[,paste0( "mut_",mut )]
				lm2 = lm(tpRank ~ purityRank, data = tc)
				raov = aov(lm2$residuals ~ am, tc)
				twa[tp,mut] = summary(raov)[[1]][["Pr(>F)"]][[1]]
				winner[tp,mut] = ifelse(median(lm2$residuals[rownames(tc[tc$am=="wt",])])>median(lm2$residuals[rownames(tc[tc$am=="mut",])]),"wt","mut" )
			}
		}
		twa[ pairs_to_validate[rownames(twa),colnames(twa)]>=0.1 ] = NA
		twa_corrected = data.frame(matrix(p.adjust(unlist(twa),method="BH"),nrow=nrow(twa),ncol=ncol(twa),dimnames=list(rownames(twa),colnames(twa))))
		save(twa,file = paste0(OutDir,"TRACERx_Quade/twa_stage",stage,".RData"))
		save(twa_corrected,file = paste0(OutDir,"TRACERx_Quade/twa_corrected_stage",stage,".RData"))
	}
}

tps_atlases_NormalCells_preprocessing = function( OutDir, withPreinvasive=F ){
	
	if (!withPreinvasive){
		cMinimal_map = c( "gray44","firebrick" )
		names(cMinimal_map) = c( "Normal","Tumor" )
		cBroad_map = c( "gray44","firebrick" )
		names(cBroad_map) = c( "Normal","Tumor" )
		cDetailed_map = c( "gray44","orange","tomato3","firebrick","brown" )
		names(cDetailed_map) = c( "Normal","I","II","III","IV")
		# creating two unique annot tables with sample-level and pst-level information, n_<ct>, ... From that, we will extract at21. NAs where needed.
		DataDir = "data/"
		CellTypeDir = paste0("data/cell_typing/non_malignant/")
		load(file = paste0(OutDir,"../scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
		scores = scores[(scores$SampleType!="Metastasis"),]
		asl = dan.df(0, c( "Dataset","Sample","PST","Patient","SampleType","Stage_collapsed",paste0("n_",c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" )) ))
		
		tasl = scores[!duplicated(scores$Sample),c( "Dataset","Sample","Patient","SampleType","Stage_collapsed")]
		rownames(tasl) = tasl$Sample
		tasl$PST = paste0(tasl$Patient,tasl$SampleType)
		for (thisct in c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" )){
			for (s in rownames(tasl)){ tasl[s,paste0( "n_",thisct )] = nrow(scores[(scores$Sample==s) & (scores$Epi_Cancer==thisct),]) }
		}
		tasl = tasl[!(tasl$Sample %in% asl$Sample),]
		asl = rbind(asl,tasl[,colnames(asl)])

		asl[asl$SampleType %in% c( "AAH","AIS","MIA" ),"Stage_collapsed"] = 0
		asl[asl$SampleType=="Normal","Stage_collapsed"] = NA
		asl$cMinimal = "Normal"
		asl$cMinimal_colorz = "gray44"
		asl$cBroad = "Normal"
		asl$cBroad_colorz = "gray44"
		asl$cDetailed = "Normal"
		asl$cDetailed_colorz = "gray44"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cMinimal"] = "Preinvasive"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cMinimal_colorz"] = "orange"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cBroad"] = asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"SampleType"]
		asl[asl$cBroad %in% c( "AAH"),"cBroad_colorz"] = "dodgerblue4"
		asl[asl$cBroad %in% c( "AIS"),"cBroad_colorz"] = "steelblue4"
		asl[asl$cBroad %in% c( "MIA"),"cBroad_colorz"] = "gold"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cDetailed"] = asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"SampleType"]
		asl[asl$SampleType %in% c( "Primary"),"cMinimal"] = "Tumor"
		asl[asl$SampleType %in% c( "Primary"),"cMinimal_colorz"] = "firebrick"
		asl[asl$SampleType %in% c( "Primary"),"cBroad"] = "Tumor"
		asl[asl$cBroad %in% c( "Tumor"),"cBroad_colorz"] = "firebrick"
		asl[asl$SampleType %in% c( "Primary"),"cDetailed"] = asl[asl$SampleType %in% c( "Primary"),"Stage_collapsed"]
		asl[(asl$cDetailed=="AAH") %in% c(T),"cDetailed_colorz"] = "dodgerblue4"
		asl[(asl$cDetailed=="AIS") %in% c(T),"cDetailed_colorz"] = "steelblue4"
		asl[(asl$cDetailed=="MIA") %in% c(T),"cDetailed_colorz"] = "gold"
		asl[(asl$cDetailed=="I") %in% c(T),"cDetailed_colorz"] = "orange"
		asl[(asl$cDetailed=="II") %in% c(T),"cDetailed_colorz"] = "tomato3"
		asl[(asl$cDetailed=="III") %in% c(T),"cDetailed_colorz"] = "firebrick"
		asl[(asl$cDetailed=="IV") %in% c(T),"cDetailed_colorz"] = "brown"
		asl[is.na(asl$cDetailed) %in% c(T),"cDetailed_colorz"] = NA
		save(asl,file=paste0(OutDir,"asl.RData"))
		asl$Sample = NULL
		tasl = asl[!duplicated(asl$PST),c( "Dataset","Patient","SampleType","Stage_collapsed","PST","cMinimal","cBroad","cDetailed","cMinimal_colorz","cBroad_colorz","cDetailed_colorz" )]
		rownames(tasl) = tasl$PST
		asl = asl[,c("PST",paste0("n_",c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" )))]
		pst = aggregate(.~PST,data=asl,FUN='sum')
		rownames(pst) = pst$PST
		pst$PST = NULL
		pst = cbind(pst,tasl[rownames(pst),])
		pst$cMinimal = factor(pst$cMinimal, levels=c( "Normal","Tumor" ))
		pst$cBroad = factor(pst$cBroad, levels=c( "Normal","Tumor" ))
		pst$cDetailed = factor(pst$cDetailed, levels=c( "Normal","I","II","III","IV",NA ))
		save(pst,file=paste0(OutDir,"pst.RData"))
		load(file=paste0(OutDir,"asl.RData"))
		asl$cMinimal = factor(asl$cMinimal, levels=c( "Normal","Tumor" ))
		asl$cBroad = factor(asl$cBroad, levels=c( "Normal","Tumor" ))
		asl$cDetailed = factor(asl$cDetailed, levels=c( "Normal","I","II","III","IV",NA ))
		save(asl,file=paste0(OutDir,"asl.RData"))
		load(file=paste0(OutDir,"pst.RData"))
		pst = aggregate(.~cDetailed,pst[,c( "cDetailed",colnames(pst)[grepl("n_",colnames(pst) )] )],FUN='sum' )
		# % alveolar cells across stages
		pst$AT1 = 100*pst$n_AT1/rowSums(pst[,grepl("n_",colnames(pst) )])
		pst$AT2 = 100*pst$n_AT2/rowSums(pst[,grepl("n_",colnames(pst) )])
		pst$AT0 = 100*pst$n_AT0/rowSums(pst[,grepl("n_",colnames(pst) )])
		
		# barplot with jitters
		mpst = melt(pst[!is.na(pst$cDetailed),c( "cDetailed","AT1","AT0","AT2" )])
		xColorz = c( AT1="#594A42",AT0="#C2B59B",AT2="#39B54A" )
		mpst$colorz = dan.expand_colors(as.character(mpst$variable),names(xColorz),as.character(xColorz))
		mpst$variable = factor(mpst$variable,levels=names(xColorz))
		mpst$cDetailed = paste0(ifelse(mpst$cDetailed=="Normal","","Stage "),mpst$cDetailed)
		plot = ggplot(data = mpst, aes(x=cDetailed,y=value,fill=variable)) + geom_bar(color='black',lwd=0.25,stat="identity", position=position_dodge(.9)) +scale_fill_manual(values = xColorz) + theme_classic(base_size=6) + xlab( "" ) + ylab( "% of alveolar cells\n(out of all epithelial + cancer cells)" )
		# +geom_jitter( position = position_jitterdodge(jitter.width=0.2), color=mpst$colorz ) + 
		plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
		pdf(paste0(OutDir,"AlveolarCellsPercent_across_stages.pdf" ),4,1.6)
		print(plot)
		dev.off()

		# AT2/AT1
		load(file=paste0(OutDir,"pst.RData"))
		at21 = pst[pst$cMinimal!="Preinvasive",]
		at21$cMinimal = factor(at21$cMinimal, levels = c( "Normal","Tumor" ))
		at21f = at21[(at21$n_AT2>=5) & (at21$n_AT1>=5),]
		at21f$AT2_AT1_ratio = log10(at21f$n_AT2/at21f$n_AT1)
		fileName = paste0(OutDir,"Ratio_AT2_AT1_cMinimal.pdf")
		x = at21f$cMinimal
		y = at21f$AT2_AT1_ratio
		xColors = cMinimal_map[levels(x)]
		jitterColors = at21f$cMinimal_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio = 1.5", signifTest = "wilcox", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 0.5, fileWidth = 2, fileHeight = 1.5, hlines_coo = log10(1.5), hlines_labels = "" )
		at21 = pst
		at21f = at21[(at21$n_AT2>=5) & (at21$n_AT1>=5),]
		at21f = at21f[!is.na(at21f$cDetailed),]
		at21f$AT2_AT1_ratio = log10(at21f$n_AT2/at21f$n_AT1)
		fileName = paste0(OutDir,"Ratio_AT2_AT1_cDetailed.pdf")
		x = factor(at21f$cDetailed,levels=names(cDetailed_map)[names(cDetailed_map) %in% unique(at21f$cDetailed)])
		y = at21f$AT2_AT1_ratio
		xColors = cDetailed_map[levels(x)]
		jitterColors = at21f$cDetailed_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio of 1.5 is shown as dashed line", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 8, fileHeight = 6, hlines_coo = log10(1.5), hlines_labels = "" )

		# AT0 ratio vs other epithelial cells (non-malignant)
		at21 = pst[pst$cMinimal!="Preinvasive",]
		at21$cMinimal = factor(at21$cMinimal, levels = c( "Normal","Tumor" ))
		at21f = at21[as.numeric(rowSums(at21[,c("n_AT1","n_AT2","n_AT0","n_preTB","n_ciliated","n_club","n_basal")]))>0,]
		at21f$AT0_ratio = at21f$n_AT0/as.numeric(rowSums(at21f[,c("n_AT1","n_AT2","n_AT0","n_preTB","n_ciliated","n_club","n_basal")]))
		fileName = paste0(OutDir,"Ratio_AT0_AllEpi_cMinimal.pdf")
		x = at21f$cMinimal
		y = at21f$AT0_ratio
		xColors = cMinimal_map[levels(x)]
		jitterColors = at21f$cMinimal_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "# AT0s / # all epithelial cells", plotTitle = "", signifTest = "wilcox", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 5, fileHeight = 5)
		at21f = at21f[!is.na(at21f$cDetailed),]
		fileName = paste0(OutDir,"Ratio_AT0_AllEpi_cDetailed.pdf")
		x = factor(at21f$cDetailed,levels=names(cDetailed_map)[names(cDetailed_map) %in% unique(at21f$cDetailed)])
		y = at21f$AT0_ratio
		xColors = cDetailed_map[levels(x)]
		jitterColors = at21f$cDetailed_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "# AT0s / # all epithelial cells", plotTitle = "", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 8, fileHeight = 6)
	} else {
		cMinimal_map = c( "gray44","orange","firebrick" )
		names(cMinimal_map) = c( "Normal","Preinvasive","Tumor" )
		cBroad_map = c( "gray44","dodgerblue4","steelblue4","gold","firebrick" )
		names(cBroad_map) = c( "Normal","AAH", "AIS", "MIA","Tumor" )
		cDetailed_map = c( "gray44","dodgerblue4","steelblue4","gold","orange","tomato3","firebrick","brown" )
		names(cDetailed_map) = c( "Normal","AAH", "AIS", "MIA","I","II","III","IV")
		# creating two unique annot tables with sample-level and pst-level information, n_<ct>, ... From that, we will extract at21. NAs where needed.
		DataDir = "data/"
		CellTypeDir = paste0("data/cell_typing/non_malignant/")
		load(file = paste0(OutDir,"../extendedAtlasS0normal_annot.RData"))
		scores = annot[(annot$SampleType!="Metastasis"),]
		asl = dan.df(0, c( "Dataset","Sample","PST","Patient","SampleType","Stage_collapsed",paste0("n_",c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" )) ))
		
		tasl = scores[!duplicated(scores$Sample),c( "Dataset","Sample","Patient","SampleType","Stage_collapsed")]
		rownames(tasl) = tasl$Sample
		tasl$PST = paste0(tasl$Patient,tasl$SampleType)
		for (thisct in c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" )){
			for (s in rownames(tasl)){ tasl[s,paste0( "n_",thisct )] = nrow(scores[(scores$Sample==s) & (scores$Epi_Cancer==thisct),]) }
		}
		tasl = tasl[!(tasl$Sample %in% asl$Sample),]
		asl = rbind(asl,tasl[,colnames(asl)])
		asl[asl$SampleType %in% c( "AAH","AIS","MIA" ),"Stage_collapsed"] = 0
		asl[asl$SampleType=="Normal","Stage_collapsed"] = NA
		asl$cMinimal = "Normal"
		asl$cMinimal_colorz = "gray44"
		asl$cBroad = "Normal"
		asl$cBroad_colorz = "gray44"
		asl$cDetailed = "Normal"
		asl$cDetailed_colorz = "gray44"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cMinimal"] = "Preinvasive"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cMinimal_colorz"] = "orange"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cBroad"] = asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"SampleType"]
		asl[asl$cBroad %in% c( "AAH"),"cBroad_colorz"] = "dodgerblue4"
		asl[asl$cBroad %in% c( "AIS"),"cBroad_colorz"] = "steelblue4"
		asl[asl$cBroad %in% c( "MIA"),"cBroad_colorz"] = "gold"
		asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"cDetailed"] = asl[asl$SampleType %in% c( "AAH","AIS","MIA"),"SampleType"]
		asl[asl$SampleType %in% c( "Primary"),"cMinimal"] = "Tumor"
		asl[asl$SampleType %in% c( "Primary"),"cMinimal_colorz"] = "firebrick"
		asl[asl$SampleType %in% c( "Primary"),"cBroad"] = "Tumor"
		asl[asl$cBroad %in% c( "Tumor"),"cBroad_colorz"] = "firebrick"
		asl[asl$SampleType %in% c( "Primary"),"cDetailed"] = asl[asl$SampleType %in% c( "Primary"),"Stage_collapsed"]
		asl[(asl$cDetailed=="AAH") %in% c(T),"cDetailed_colorz"] = "dodgerblue4"
		asl[(asl$cDetailed=="AIS") %in% c(T),"cDetailed_colorz"] = "steelblue4"
		asl[(asl$cDetailed=="MIA") %in% c(T),"cDetailed_colorz"] = "gold"
		asl[(asl$cDetailed=="I") %in% c(T),"cDetailed_colorz"] = "orange"
		asl[(asl$cDetailed=="II") %in% c(T),"cDetailed_colorz"] = "tomato3"
		asl[(asl$cDetailed=="III") %in% c(T),"cDetailed_colorz"] = "firebrick"
		asl[(asl$cDetailed=="IV") %in% c(T),"cDetailed_colorz"] = "brown"
		asl[is.na(asl$cDetailed) %in% c(T),"cDetailed_colorz"] = NA
		save(asl,file=paste0(OutDir,"asl.RData"))
		asl$Sample = NULL
		tasl = asl[!duplicated(asl$PST),c( "Dataset","Patient","SampleType","Stage_collapsed","PST","cMinimal","cBroad","cDetailed","cMinimal_colorz","cBroad_colorz","cDetailed_colorz" )]
		rownames(tasl) = tasl$PST
		asl = asl[,c("PST",paste0("n_",c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" )))]
		pst = aggregate(.~PST,data=asl,FUN='sum')
		rownames(pst) = pst$PST
		pst$PST = NULL
		pst = cbind(pst,tasl[rownames(pst),])
		pst$cMinimal = factor(pst$cMinimal, levels=c( "Normal","Preinvasive","Tumor" ))
		pst$cBroad = factor(pst$cBroad, levels=c( "Normal","AAH", "AIS", "MIA","Tumor" ))
		pst$cDetailed = factor(pst$cDetailed, levels=c( "Normal","AAH", "AIS", "MIA","I","II","III","IV",NA ))
		save(pst,file=paste0(OutDir,"pst.RData"))
		load(file=paste0(OutDir,"asl.RData"))
		asl$cMinimal = factor(asl$cMinimal, levels=c( "Normal","Preinvasive","Tumor" ))
		asl$cBroad = factor(asl$cBroad, levels=c( "Normal","AAH", "AIS", "MIA","Tumor" ))
		asl$cDetailed = factor(asl$cDetailed, levels=c( "Normal","AAH", "AIS", "MIA","I","II","III","IV",NA ))
		save(asl,file=paste0(OutDir,"asl.RData"))

		# AT2/AT1
		load(file=paste0(OutDir,"pst.RData"))
		at21 = pst
		at21$cMinimal = factor(at21$cMinimal, levels=c( "Normal","Preinvasive","Tumor" ))
		at21f = at21[(at21$n_AT2>=5) & (at21$n_AT1>=5),]
		at21f$AT2_AT1_ratio = log10(at21f$n_AT2/at21f$n_AT1)
		fileName = paste0(OutDir,"Ratio_AT2_AT1_cMinimal.pdf")
		x = at21f$cMinimal
		y = at21f$AT2_AT1_ratio
		xColors = cMinimal_map[levels(x)]
		jitterColors = at21f$cMinimal_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio = 1.5", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 1, fileWidth = 2.2, fileHeight = 2, hlines_coo = log10(1.5), hlines_labels = "" )
		at21 = pst
		at21$cBroad = factor(at21$cBroad, levels=c( "Normal","AAH", "AIS", "MIA","Tumor" ))
		at21f = at21[(at21$n_AT2>=5) & (at21$n_AT1>=5),]
		at21f$AT2_AT1_ratio = log10(at21f$n_AT2/at21f$n_AT1)
		fileName = paste0(OutDir,"Ratio_AT2_AT1_cBroad.pdf")
		x = at21f$cBroad
		y = at21f$AT2_AT1_ratio
		xColors = cBroad_map[levels(x)]
		jitterColors = at21f$cBroad_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio = 1.5", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 7, fileHeight = 5, hlines_coo = log10(1.5), hlines_labels = "" )
		at21 = pst
		at21f = at21[(at21$n_AT2>=5) & (at21$n_AT1>=5),]
		at21f = at21f[!is.na(at21f$cDetailed),]
		at21f$AT2_AT1_ratio = log10(at21f$n_AT2/at21f$n_AT1)
		fileName = paste0(OutDir,"Ratio_AT2_AT1_cDetailed.pdf")
		x = factor(at21f$cDetailed,levels=names(cDetailed_map)[names(cDetailed_map) %in% unique(at21f$cDetailed)])
		y = at21f$AT2_AT1_ratio
		xColors = cDetailed_map[levels(x)]
		jitterColors = at21f$cDetailed_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio = 1.5", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 8, fileHeight = 6, hlines_coo = log10(1.5), hlines_labels = "" )

		# AT0 ratio vs other epithelial cells (non-malignant)
		at21 = pst
		at21$cMinimal = factor(at21$cMinimal, levels = c( "Normal","Preinvasive","Tumor" ))
		at21f = at21[as.numeric(rowSums(at21[,c("n_AT1","n_AT2","n_AT0","n_preTB","n_ciliated","n_club","n_basal")]))>0,]
		at21f$AT0_ratio = at21f$n_AT0/as.numeric(rowSums(at21f[,c("n_AT1","n_AT2","n_AT0","n_preTB","n_ciliated","n_club","n_basal")]))
		fileName = paste0(OutDir,"Ratio_AT0_AllEpi_cMinimal.pdf")
		x = at21f$cMinimal
		y = at21f$AT0_ratio
		xColors = cMinimal_map[levels(x)]
		jitterColors = at21f$cMinimal_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "# AT0s / # all epithelial cells", plotTitle = "", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 1, fileWidth = 2.2, fileHeight = 2)
		at21 = pst
		at21f = at21[as.numeric(rowSums(at21[,c("n_AT1","n_AT2","n_AT0","n_preTB","n_ciliated","n_club","n_basal")]))>0,]
		at21f$AT0_ratio = at21f$n_AT0/as.numeric(rowSums(at21f[,c("n_AT1","n_AT2","n_AT0","n_preTB","n_ciliated","n_club","n_basal")]))
		fileName = paste0(OutDir,"Ratio_AT0_AllEpi_cBroad.pdf")
		x = at21f$cBroad
		y = at21f$AT0_ratio
		xColors = cBroad_map[levels(x)]
		jitterColors = at21f$cBroad_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "# AT0s / # all epithelial cells", plotTitle = "", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 7, fileHeight = 5)
		at21f = at21f[!is.na(at21f$cDetailed),]
		fileName = paste0(OutDir,"Ratio_AT0_AllEpi_cDetailed.pdf")
		x = factor(at21f$cDetailed,levels=names(cDetailed_map)[names(cDetailed_map) %in% unique(at21f$cDetailed)])
		y = at21f$AT0_ratio
		xColors = cDetailed_map[levels(x)]
		jitterColors = at21f$cDetailed_colorz
		dan.boxplots( fileName, x = x, y = y, xlab = "", ylab = "# AT0s / # all epithelial cells", plotTitle = "", signifTest = "kruskal", labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, jitterDotSize = 3, fileWidth = 8, fileHeight = 6)
	}
}

tps_atlases_NormalCells_scoring = function( OutDir, withPreinvasive=F ){
	# gs construction
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/Travaglini_epithelial_markers.txt"))
	for (cn in colnames(aa)){ GeneSets[[paste0(cn)]] = aa[,cn][nchar(aa[,cn])>0] }
	names(GeneSets) = gsub("_CELL","",names(GeneSets))
	names(GeneSets)[names(GeneSets)=="ALVEOLAR_EPITHELIAL_TYPE_1"] = "AT1"
	names(GeneSets)[names(GeneSets)=="ALVEOLAR_EPITHELIAL_TYPE_2"] = "AT2"
	names(GeneSets)[names(GeneSets)=="SIGNALING_ALVEOLAR_EPITHELIAL_TYPE_2"] = "signaling_AT2"
	names(GeneSets) = paste0( "Trav_",names(GeneSets) )
	gs = GeneSets[sort(names(GeneSets))]
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/HLCA_epithelial_markers.txt"))
	for (cn in colnames(aa)){ GeneSets[[paste0(cn)]] = aa[,cn][nchar(aa[,cn])>0] }
	names(GeneSets) = paste0( "HLCA_",names(GeneSets) )
	GeneSets = GeneSets[sort(names(GeneSets))]
	gs = c(gs,GeneSets)
	cp = read.table(file=paste0("processed/for_Amaia_7Dec2023/","LocardPaulet_proliferationSignature.txt"),header=T)
	gs$LP_proliferation = cp$gene
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/Han_epithelial_markers.txt"))
	GeneSets = list()
	for (cn in colnames(aa)){ GeneSets[[paste0(cn)]] = aa[,cn][nchar(aa[,cn])>0] }
	names(GeneSets) = paste0("Han_",names(GeneSets) )
	gs = c(gs,GeneSets)
	hkc = dan.read("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/han_kac_deg_complete_table8.txt")
	hkc_top = hkc[order(hkc$avg_log2FC,decreasing=T),"gene"][1:100]
	gs$Han_KAC_signature = hkc_top
	save(gs,file = paste0( OutDir,"gs_TravHlcaLPHan.RData" ))
	load(file = paste0( OutDir,"gs_TravHlcaLPHan.RData" ))

	if (!withPreinvasive){
		load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelial.RData"))
		load(file = paste0(OutDir,"../extendedAtlasNormal_annot.RData"))
		ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
		mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
		ribo.genes = ribo$Approved.symbol
		seu = subset(x = seu, features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))]) # 12677 x 101849
		data = seu@assays$SCT@data
		rm(seu)
		gc()
		db_rand = dan.barkley_MakeRand_data(data,gs, 3)
		save(db_rand,file = paste0(OutDir,"db_rand_extendedAtlasNormal.RData"))
		load(file = paste0(OutDir,"db_rand_extendedAtlasNormal.RData"))
		scores_uncorrected = dan.Barkley_GeneToEnrichment_AmsCentered_data(data, annot, gs, db_rand)
		save(scores_uncorrected,file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
		load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
		regress_scores_normal_epi( OutDir, scores_uncorrected, prefix="gsAtlases", ct_estimate_coefficients=c("AT2"), tps=gs )
		load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelial.RData"))
		load(file = paste0(OutDir,"../extendedAtlasNormal_annot.RData"))
		ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
		mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
		ribo.genes = ribo$Approved.symbol
		seu = subset(x = seu, features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))], cells = annot[annot$Epi_Cancer!="Cancer","CellID"]) # 12677 x 101849
		data = seu@assays$SCT@data
		rm(seu)
		gc()
		db_rand = dan.barkley_MakeRand_data(data,gs, 3)
		save(db_rand,file = paste0(OutDir,"db_rand_onlyNormalEpi.RData"))
		scores_uncorrected = dan.Barkley_GeneToEnrichment_AmsCentered_data(data, annot, gs, db_rand)
		save(scores_uncorrected,file = paste0(OutDir,"gs_scores_onlyNormalEpi.RData"))
		load(file = paste0(OutDir,"gs_scores_onlyNormalEpi.RData"))
		regress_scores_normal_epi( OutDir, scores_uncorrected, prefix="gsAtlases_OnlyNormalEpi", ct_estimate_coefficients=c("AT2"), tps=gs )
	} else {
		load( file = paste0(OutDir,"../tps_discovery/tps.RData") )
		load( file = paste0(OutDir,"../extendedAtlasS0normal_annot.RData") )
		load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelialS0.RData"))
		ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
		mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
		ribo.genes = ribo$Approved.symbol
		data = seu@assays$SCT@data
		data = data[rownames(data)[!(rownames(data) %in% c(mito.genes,ribo.genes))],]
		rm(seu)
		gc()
		gs = gs[ c("HLCA_AT0","HLCA_AT1","HLCA_AT2","LP_proliferation","Han_AT2","Han_AT1","Han_KAC_signature") ]
		gs = c(gs,tps)
		db_rand = dan.barkley_MakeRand_data(data,gs,3)
		save(db_rand,file = paste0(OutDir,"db_rand_extendedAtlasNormalS0.RData"))
		scores_uncorrected = dan.Barkley_GeneToEnrichment_AmsCentered_data(data, annot, gs, db_rand)
		save(scores_uncorrected,file = paste0(OutDir,"gs_scores_extendedAtlasNormalS0.RData"))
		regress_scores_normal_epi( OutDir, scores_uncorrected, prefix="gsAtlases", ct_estimate_coefficients=c("AT2"), tps=gs )

		# load( file = paste0(OutDir,"../extendedAtlasS0normal_annot.RData") )
		# load( file = paste0(OutDir,"../extendedAtlasS0normal_counts.RData") )
		# annot = annot[(annot$Epi_Cancer!="Cancer") & (annot$SampleType!="Metastasis"),]
		# counts = counts[,rownames(annot)]
		# seu = CreateSeuratObject(counts = counts, project = "extendedAtlasNormalS0", meta.data = annot)
		# seu = SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("percent.mt"), verbose = FALSE)
		# seu = RunPCA(seu)
		# save(seu, file = paste0(OutDir,"extendedAtlasNormalS0_seu_AT120.RData"))
		# ribo = read.csv(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/","hgnc_ribosomial_proteins.csv"))
		# mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
		# ribo.genes = ribo$Approved.symbol
		# data = seu@assays$SCT@data
		# data = data[rownames(data)[!(rownames(data) %in% c(mito.genes,ribo.genes))],]
		# rm(seu)
		# gc()
		# gs = gs[ c("HLCA_AT0","HLCA_AT1","HLCA_AT2","LP_proliferation","Han_AT2","Han_AT1","Han_KAC_signature") ]
		# gs = c(gs,tps)
		# db_rand = dan.barkley_MakeRand_data(data,gs,3)
		# save(db_rand,file = paste0(OutDir,"db_rand_onlyNormalEpi.RData"))
		# scores_uncorrected = dan.Barkley_GeneToEnrichment_AmsCentered_data(data, annot, gs, db_rand)
		# save(scores_uncorrected,file = paste0(OutDir,"gs_scores_onlyNormalEpi.RData"))
		# load(file = paste0(OutDir,"gs_scores_onlyNormalEpi.RData"))
		# regress_scores_normal_epi( OutDir, scores_uncorrected, prefix="gsAtlases_OnlyNormalEpi", ct_estimate_coefficients=c("AT2"), tps=gs )
	}
}

tps_atlases_NormalCells = function( OutDir, order_tps ){
	
	celltype_map = data.frame(row.names=c( "AT1","AT2","AT0","preTB","ciliated","club","basal","Cancer" ),
		colorz = c("chocolate","dodgerblue2","purple","forestgreen","gold","orange","brown","firebrick") ,stringsAsFactors=F )
	cMinimal_map = c( "gray44","firebrick" )
	names(cMinimal_map) = c( "Normal","Tumor" )
	cBroad_map = c( "gray44","firebrick" )
	names(cBroad_map) = c( "Normal","Tumor" )
	cDetailed_map = c( "gray44","orange","tomato3","firebrick","brown" )
	names(cDetailed_map) = c( "Normal","I","II","III","IV")
	DataDir = "data/"
	CellTypeDir = paste0("data/cell_typing/non_malignant/")

	load(file=paste0(OutDir,"asl.RData"))
	load(file=paste0(OutDir,"pst.RData"))
	load(paste0(OutDir,"../tps_discovery/tps.RData"))
	load(paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	tps = tps[order_tps]
	load(file = paste0(OutDir,"../scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	
	# load(file = paste0(MainDir,"scoring/tps_scores_extendedAtlasNormal.RData"))
	# scores_cancer = aggregate(.~Sample,scores[scores$Epi_Cancer=="Cancer",c("AT2-like","Sample")],FUN='mean')
	# scores_at2 = aggregate(.~Sample,scores[scores$Epi_Cancer=="AT2",c("AT2-like","Sample")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Sample
	# rownames(scores_at2) = scores_at2$Sample
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, in_cancer=scores_cancer[commonz,"AT2-like"], in_at2=scores_at2[commonz,"AT2-like"])	
	# scores_cancer = scores[scores$Epi_Cancer=="Cancer",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Sample),]
	# rownames(scores_cancer) = scores_cancer$Sample
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc$in_cancer,sc$in_at2)


	### meansc tps:
	dir.create(paste0( OutDir,"meansc_tps/" ) )
	for (thisct in c( "AT1","AT2","AT0","preTB","ciliated","club","basal" )){
		dcat(thisct)
		load(file=paste0(OutDir,"pst.RData"))
		sco = scores[scores$Epi_Cancer==thisct,]
		if (nrow(sco)==0){ 
			dcat( "No rows" ) 
			next
		}
		for (rn in rownames(pst)){
			for (tp in names(tps)){
				pst[rn,tp] = mean(sco[(sco$Patient==pst[rn,"Patient"]) & (sco$SampleType==pst[rn,"SampleType"]),tp ])
			}
		}
		pst = pst[!is.na(pst[,names(tps)[1]]),]
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_Boxplots.pdf"),7,5 )
		for (tp in names(tps)){
			x = pst$cMinimal
			y = pst[,tp]
			xColors = cMinimal_map[levels(x)]
			jitterColors = pst$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_PairedDotPlot.pdf"),7,5 )
		for (tp in names(tps)){
			patients_multiple_samples = rowSums(dtable(pst$Patient,pst$cMinimal)>0)
			patients_multiple_samples = names(which(patients_multiple_samples>=2))
			tpst = pst[pst$Patient %in% patients_multiple_samples,]
			dftest = tpst[tpst$cMinimal %in% c( "Normal","Tumor" ),c( "Patient","cMinimal",tp )]
			patients_two = rowSums(dtable(dftest$Patient,dftest$cMinimal)>0)
			patients_two = names(which(patients_two==2))
			dftest = dftest[dftest$Patient %in% patients_two,]
			dftest = dftest[order(dftest$cMinimal,dftest$Patient),]
			wt = wilcox.test(x=dftest[dftest$cMinimal=="Normal",tp],y=dftest[dftest$cMinimal=="Tumor",tp],paired=TRUE)
			plotTitle = paste0( "Normal vs Tumor\nPaired Wilcoxon test, p-val = ",signif(wt$p.value,2) )
			x = tpst$cMinimal
			y = tpst[,tp]
			pairingz = tpst$Patient
			xColors = cMinimal_map[levels(x)]
			jitterColors = tpst$cMinimal_colorz
			plot=dan.pairedDotPlot.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = plotTitle, labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = jitterColors, labelPoints = "", linesPairings = pairingz )
			print(plot)
		}
		dev.off()
		pst2 = pst[pst$cMinimal!="Preinvasive",]
		volcano_df = dan.df( names(tps),c("patient_average_difference","paired_pval") )
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_OnlyNormalVsTumor_PairedDotPlot.pdf"),5,5 )
		for (tp in names(tps)){
			patients_multiple_samples = rowSums(dtable(pst2$Patient,pst2$cMinimal)>0)
			patients_multiple_samples = names(which(patients_multiple_samples>=2))
			tpst = pst2[pst2$Patient %in% patients_multiple_samples,]
			dftest = tpst[tpst$cMinimal %in% c( "Normal","Tumor" ),c( "Patient","cMinimal",tp )]
			patients_two = rowSums(dtable(dftest$Patient,dftest$cMinimal)>0)
			patients_two = names(which(patients_two==2))
			dftest = dftest[dftest$Patient %in% patients_two,]
			dftest = dftest[order(dftest$cMinimal,dftest$Patient),]
			volcano_df[tp,"patient_average_difference"] = mean( dftest[dftest$cMinimal=="Tumor",tp]-dftest[dftest$cMinimal=="Normal",tp] )
			wt = wilcox.test(x=dftest[dftest$cMinimal=="Normal",tp],y=dftest[dftest$cMinimal=="Tumor",tp],paired=TRUE)
			volcano_df[tp,"paired_pval"] = wt$p.value
			plotTitle = paste0( "Normal vs Tumor\nPaired Wilcoxon test, p-val = ",signif(wt$p.value,2) )
			x = tpst$cMinimal
			y = tpst[,tp]
			pairingz = tpst$Patient
			xColors = cMinimal_map[levels(x)]
			jitterColors = tpst$cMinimal_colorz
			plot=dan.pairedDotPlot.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = plotTitle, labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = jitterColors, labelPoints = "", linesPairings = pairingz )
			print(plot)
		}
		dev.off()
		## volcano paired 
		volcano_df$paired_qval = p.adjust(volcano_df$paired_pval,method="BH")
		fileName = paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_OnlyNormalVsTumor_VolcanoPaired.pdf")
		x = volcano_df$patient_average_difference
		y = -log10(volcano_df$paired_pval)
		volcano_df$fill = ifelse( (volcano_df$paired_pval>=0.05) | (abs(volcano_df$patient_average_difference)<0.1),"Not significant",ifelse( volcano_df$patient_average_difference<0,"Higher in normal","Higher in tumor"))
		volcano_df$fill = factor(volcano_df$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
		fillColorz = c( "gray44","firebrick","gray88" )
		volcano_df$repel_labelz = rownames(volcano_df)
		volcano_df[volcano_df$fill=="Not significant","repel_labelz"] = ""
		dan.scatterplot( fileName, x, y, fill = volcano_df$fill, xlab = "Average patient score difference", ylab = "-log10(p-value), paired Wilcoxon test", filllab = "", plotVline=c(-0.1,0.1),plotHline=0.05,plotTitle=paste0( "in ",thisct," cells" ) , fillColors = fillColorz, dotSize = 3, repel_labels=volcano_df$repel_labelz, coord_fixed = FALSE, fileWidth = 5, fileHeight = 4 )
		volcano_df = dan.df( names(tps),c("average_difference","pval") )
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_OnlyNormalVsTumor_Boxplots.pdf"),5,5 )
		for (tp in names(tps)){
			x = pst2$cMinimal
			y = pst2[,tp]
			wt = wilcox.test(x=y[x=="Normal"],y=y[x=="Tumor"])
			volcano_df[tp,"average_difference"] = mean(y[x=="Tumor"])-mean(y[x=="Normal"])
			volcano_df[tp,"pval"] = wt$p.value
			xColors = cMinimal_map[levels(x)[levels(x) %in% unique(as.character(x))]]
			jitterColors = pst2$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		## volcano unpaired 
		volcano_df$qval = p.adjust(volcano_df$pval,method='BH')
		fileName = paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_OnlyNormalVsTumor_VolcanoUnpaired.pdf")
		x = volcano_df$average_difference
		y = -log10(volcano_df$pval)
		volcano_df$fill = ifelse( (volcano_df$pval>=0.05) | (abs(volcano_df$average_difference)<0.1),"Not significant",ifelse( volcano_df$average_difference<0,"Higher in normal","Higher in tumor"))
		volcano_df$fill = factor(volcano_df$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
		fillColorz = c( "gray44","firebrick","gray88" )
		volcano_df$repel_labelz = rownames(volcano_df)
		volcano_df[volcano_df$fill=="Not significant","repel_labelz"] = ""
		dan.scatterplot( fileName, x, y, fill = volcano_df$fill, xlab = "Average patient-wise score difference", ylab = "-log10(p-value), Wilcoxon test", filllab = "", plotTitle=paste0( "Transcriptional programs in normal ",thisct," cells" ) , fillColors = fillColorz, dotSize = 3, repel_labels=volcano_df$repel_labelz, coord_fixed = FALSE, fileWidth = 7, fileHeight = 5 )
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cBroad_Boxplots.pdf"),8,5 )
		for (tp in names(tps)){
			x = pst$cBroad
			y = pst[,tp]
			xColors = cBroad_map[levels(x)]
			jitterColors = pst$cBroad_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cDetailed_Boxplots.pdf"),10,5 )
		for (tp in names(tps)){
			x = pst$cDetailed[!is.na(pst$cDetailed)]
			y = pst[!is.na(pst$cDetailed),tp]
			xColors = cDetailed_map[levels(x)]
			jitterColors = pst$cDetailed_colorz[!is.na(pst$cDetailed)]
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
	}

	### meansc cytotrace:
	dir.create(paste0( OutDir,"meansc_cytotrace/" ) )
	for (thisct in c( "AT1","AT2","AT0" )){
		dcat(thisct)
		load(file=paste0(OutDir,"pst.RData"))
		sco = scores[scores$Epi_Cancer==thisct,]
		load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
		annot = annot[rownames(sco),]
		if (nrow(sco)==0){ 
			dcat( "No rows" ) 
			next
		}
		for (rn in rownames(pst)){
			pst[rn,"cyto"] = mean(annot[(annot$Patient==pst[rn,"Patient"]) & (annot$SampleType==pst[rn,"SampleType"]),"cyto" ])
		}
		pst = pst[!is.na(pst[, "cyto" ]),]
		pdf( paste0(OutDir,"meansc_cytotrace/meansc_cyto_",thisct,"_cMinimal_Boxplots.pdf"),7,5 )
		for (tp in "cyto"){
			x = pst$cMinimal
			y = pst[,tp]
			xColors = cMinimal_map[levels(x)]
			jitterColors = pst$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_cytotrace/meansc_cyto_",thisct,"_cMinimal_PairedDotPlot.pdf"),7,5 )
		for (tp in "cyto"){
			patients_multiple_samples = rowSums(dtable(pst$Patient,pst$cMinimal)>0)
			patients_multiple_samples = names(which(patients_multiple_samples>=2))
			tpst = pst[pst$Patient %in% patients_multiple_samples,]
			dftest = tpst[tpst$cMinimal %in% c( "Normal","Tumor" ),c( "Patient","cMinimal",tp )]
			patients_two = rowSums(dtable(dftest$Patient,dftest$cMinimal)>0)
			patients_two = names(which(patients_two==2))
			dftest = dftest[dftest$Patient %in% patients_two,]
			dftest = dftest[order(dftest$cMinimal,dftest$Patient),]
			wt = wilcox.test(x=dftest[dftest$cMinimal=="Normal",tp],y=dftest[dftest$cMinimal=="Tumor",tp],paired=TRUE)
			plotTitle = paste0( "Normal vs Tumor\nPaired Wilcoxon test, p-val = ",signif(wt$p.value,2) )
			x = tpst$cMinimal
			y = tpst[,tp]
			pairingz = tpst$Patient
			xColors = cMinimal_map[levels(x)]
			jitterColors = tpst$cMinimal_colorz
			plot=dan.pairedDotPlot.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = plotTitle, labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = jitterColors, labelPoints = "", linesPairings = pairingz )
			print(plot)
		}
		dev.off()
		pst2 = pst[pst$cMinimal!="Preinvasive",]
		pdf( paste0(OutDir,"meansc_cytotrace/meansc_cyto_",thisct,"_cMinimal_OnlyNormalVsTumor_PairedDotPlot.pdf"),5,5 )
		for (tp in "cyto"){
			patients_multiple_samples = rowSums(dtable(pst2$Patient,pst2$cMinimal)>0)
			patients_multiple_samples = names(which(patients_multiple_samples>=2))
			tpst = pst2[pst2$Patient %in% patients_multiple_samples,]
			dftest = tpst[tpst$cMinimal %in% c( "Normal","Tumor" ),c( "Patient","cMinimal",tp )]
			patients_two = rowSums(dtable(dftest$Patient,dftest$cMinimal)>0)
			patients_two = names(which(patients_two==2))
			dftest = dftest[dftest$Patient %in% patients_two,]
			dftest = dftest[order(dftest$cMinimal,dftest$Patient),]
			wt = wilcox.test(x=dftest[dftest$cMinimal=="Normal",tp],y=dftest[dftest$cMinimal=="Tumor",tp],paired=TRUE)
			plotTitle = paste0( "Normal vs Tumor\nPaired Wilcoxon test, p-val = ",signif(wt$p.value,2) )
			x = tpst$cMinimal
			y = tpst[,tp]
			pairingz = tpst$Patient
			xColors = cMinimal_map[levels(x)]
			jitterColors = tpst$cMinimal_colorz
			plot=dan.pairedDotPlot.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = plotTitle, labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = jitterColors, labelPoints = "", linesPairings = pairingz )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_cytotrace/meansc_cyto_",thisct,"_cMinimal_OnlyNormalVsTumor_Boxplots.pdf"),5,5 )
		for (tp in "cyto"){
			x = pst2$cMinimal
			y = pst2[,tp]
			xColors = cMinimal_map[levels(x)[levels(x) %in% unique(as.character(x))]]
			jitterColors = pst2$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_cytotrace/meansc_cyto_",thisct,"_cBroad_Boxplots.pdf"),8,5 )
		for (tp in "cyto"){
			x = pst$cBroad
			y = pst[,tp]
			xColors = cBroad_map[levels(x)]
			jitterColors = pst$cBroad_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_cytotrace/meansc_cyto_",thisct,"_cDetailed_Boxplots.pdf"),10,5 )
		for (tp in "cyto"){
			x = pst$cDetailed[!is.na(pst$cDetailed)]
			y = pst[!is.na(pst$cDetailed),tp]
			xColors = cDetailed_map[levels(x)]
			jitterColors = pst$cDetailed_colorz[!is.na(pst$cDetailed)]
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
	}

	# # do scores correlate in AT2 and cancer cells?
	# load(paste0(OutDir,"../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	# states = scores

	# load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	# scores = scores_uncorrected
	# scores[rownames(states),"Epi_Cancer"] = states$cs
	# scores_cancer = aggregate(.~Sample,scores[scores$Epi_Cancer=="cs1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Sample")],FUN='mean')
	# scores_at2 = aggregate(.~Sample,scores[scores$Epi_Cancer=="AT2",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Sample")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Sample
	# rownames(scores_at2) = scores_at2$Sample
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, 
	# 	in_cancer_hlcaat2=scores_cancer[commonz,"HLCA_AT2"], in_at2_hlcaat2=scores_at2[commonz,"HLCA_AT2"],
	# 	in_cancer_hlcaat1=scores_cancer[commonz,"HLCA_AT1"], in_at2_hlcaat1=scores_at2[commonz,"HLCA_AT1"],
	# 	in_cancer_hlcaat0=scores_cancer[commonz,"HLCA_AT0"], in_at2_hlcaat0=scores_at2[commonz,"HLCA_AT0"])	
	# scores_cancer = scores[scores$Epi_Cancer=="cs1",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Sample),]
	# rownames(scores_cancer) = scores_cancer$Sample
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc[,"in_cancer_hlcaat2"],sc[,"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="I") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="I") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="II") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="II") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="III") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="III") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[,"in_cancer_hlcaat1"],sc[,"in_at2_hlcaat1"])
	# cor(sc[,"in_cancer_hlcaat0"],sc[,"in_at2_hlcaat0"])

	# load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	# scores = scores_uncorrected
	# scores[rownames(states),"Epi_Cancer"] = states$cs
	# scores_cancer = aggregate(.~Sample,scores[scores$Epi_Cancer=="cs1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Sample")],FUN='mean')
	# scores_at2 = aggregate(.~Sample,scores[scores$Epi_Cancer=="AT1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Sample")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Sample
	# rownames(scores_at2) = scores_at2$Sample
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, 
	# 	in_cancer_hlcaat2=scores_cancer[commonz,"HLCA_AT2"], in_at2_hlcaat2=scores_at2[commonz,"HLCA_AT2"],
	# 	in_cancer_hlcaat1=scores_cancer[commonz,"HLCA_AT1"], in_at2_hlcaat1=scores_at2[commonz,"HLCA_AT1"],
	# 	in_cancer_hlcaat0=scores_cancer[commonz,"HLCA_AT0"], in_at2_hlcaat0=scores_at2[commonz,"HLCA_AT0"])	
	# scores_cancer = scores[scores$Epi_Cancer=="cs1",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Sample),]
	# rownames(scores_cancer) = scores_cancer$Sample
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc[,"in_cancer_hlcaat1"],sc[,"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="I") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="I") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="II") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="II") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="III") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="III") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[,"in_cancer_hlcaat2"],sc[,"in_at2_hlcaat2"])
	# cor(sc[,"in_cancer_hlcaat0"],sc[,"in_at2_hlcaat0"])

	# load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	# scores = scores_uncorrected
	# scores[rownames(states),"Epi_Cancer"] = states$cs
	# scores_cancer = aggregate(.~Sample,scores[scores$Epi_Cancer=="cs1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Sample")],FUN='mean')
	# scores_at2 = aggregate(.~Sample,scores[scores$Epi_Cancer=="AT0",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Sample")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Sample
	# rownames(scores_at2) = scores_at2$Sample
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, 
	# 	in_cancer_hlcaat2=scores_cancer[commonz,"HLCA_AT2"], in_at2_hlcaat2=scores_at2[commonz,"HLCA_AT2"],
	# 	in_cancer_hlcaat1=scores_cancer[commonz,"HLCA_AT1"], in_at2_hlcaat1=scores_at2[commonz,"HLCA_AT1"],
	# 	in_cancer_hlcaat0=scores_cancer[commonz,"HLCA_AT0"], in_at2_hlcaat0=scores_at2[commonz,"HLCA_AT0"])	
	# scores_cancer = scores[scores$Epi_Cancer=="cs1",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Sample),]
	# rownames(scores_cancer) = scores_cancer$Sample
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc[,"in_cancer_hlcaat0"],sc[,"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="I") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="I") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="II") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="II") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="III") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="III") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[,"in_cancer_hlcaat2"],sc[,"in_at2_hlcaat2"])
	# cor(sc[,"in_cancer_hlcaat1"],sc[,"in_at2_hlcaat1"])

	# ## vs normal counterpart (expected correlations):
	# load(paste0(OutDir,"../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	# states = scores

	# load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	# scores = scores_uncorrected
	# scores[rownames(states),"Epi_Cancer"] = states$cs
	# scores_cancer = aggregate(.~Patient,scores[scores$Epi_Cancer=="cs1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Patient")],FUN='mean')
	# scores_at2 = aggregate(.~Patient,scores[(scores$Epi_Cancer=="AT2") & (scores$SampleType=="Normal"),c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Patient")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Patient
	# rownames(scores_at2) = scores_at2$Patient
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, 
	# 	in_cancer_hlcaat2=scores_cancer[commonz,"HLCA_AT2"], in_at2_hlcaat2=scores_at2[commonz,"HLCA_AT2"],
	# 	in_cancer_hlcaat1=scores_cancer[commonz,"HLCA_AT1"], in_at2_hlcaat1=scores_at2[commonz,"HLCA_AT1"],
	# 	in_cancer_hlcaat0=scores_cancer[commonz,"HLCA_AT0"], in_at2_hlcaat0=scores_at2[commonz,"HLCA_AT0"])	
	# scores_cancer = scores[scores$Epi_Cancer=="cs1",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Patient),]
	# rownames(scores_cancer) = scores_cancer$Patient
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc[,"in_cancer_hlcaat2"],sc[,"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="I") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="I") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="II") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="II") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="III") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="III") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_cancer_hlcaat2"],sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_at2_hlcaat2"])
	# cor(sc[,"in_cancer_hlcaat1"],sc[,"in_at2_hlcaat1"])
	# cor(sc[,"in_cancer_hlcaat0"],sc[,"in_at2_hlcaat0"])

	# load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	# scores = scores_uncorrected
	# scores[rownames(states),"Epi_Cancer"] = states$cs
	# scores_cancer = aggregate(.~Patient,scores[scores$Epi_Cancer=="cs1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Patient")],FUN='mean')
	# scores_at2 = aggregate(.~Patient,scores[(scores$Epi_Cancer=="AT1") & (scores$SampleType=="Normal"),c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Patient")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Patient
	# rownames(scores_at2) = scores_at2$Patient
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, 
	# 	in_cancer_hlcaat2=scores_cancer[commonz,"HLCA_AT2"], in_at2_hlcaat2=scores_at2[commonz,"HLCA_AT2"],
	# 	in_cancer_hlcaat1=scores_cancer[commonz,"HLCA_AT1"], in_at2_hlcaat1=scores_at2[commonz,"HLCA_AT1"],
	# 	in_cancer_hlcaat0=scores_cancer[commonz,"HLCA_AT0"], in_at2_hlcaat0=scores_at2[commonz,"HLCA_AT0"])	
	# scores_cancer = scores[scores$Epi_Cancer=="cs1",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Patient),]
	# rownames(scores_cancer) = scores_cancer$Patient
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc[,"in_cancer_hlcaat1"],sc[,"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="I") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="I") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="II") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="II") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="III") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="III") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_cancer_hlcaat1"],sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_at2_hlcaat1"])
	# cor(sc[,"in_cancer_hlcaat2"],sc[,"in_at2_hlcaat2"])
	# cor(sc[,"in_cancer_hlcaat0"],sc[,"in_at2_hlcaat0"])

	# load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	# scores = scores_uncorrected
	# scores[rownames(states),"Epi_Cancer"] = states$cs
	# scores_cancer = aggregate(.~Patient,scores[scores$Epi_Cancer=="cs1",c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Patient")],FUN='mean')
	# scores_at2 = aggregate(.~Patient,scores[(scores$Epi_Cancer=="AT0") & (scores$SampleType=="Normal"),c("HLCA_AT2","HLCA_AT1","HLCA_AT0","Patient")],FUN='mean')
	# rownames(scores_cancer) = scores_cancer$Patient
	# rownames(scores_at2) = scores_at2$Patient
	# commonz = intersect(rownames(scores_cancer),rownames(scores_at2))
	# sc = data.frame(row.names=commonz, samples=commonz, 
	# 	in_cancer_hlcaat2=scores_cancer[commonz,"HLCA_AT2"], in_at2_hlcaat2=scores_at2[commonz,"HLCA_AT2"],
	# 	in_cancer_hlcaat1=scores_cancer[commonz,"HLCA_AT1"], in_at2_hlcaat1=scores_at2[commonz,"HLCA_AT1"],
	# 	in_cancer_hlcaat0=scores_cancer[commonz,"HLCA_AT0"], in_at2_hlcaat0=scores_at2[commonz,"HLCA_AT0"])	
	# scores_cancer = scores[scores$Epi_Cancer=="cs1",]
	# scores_cancer = scores_cancer[!duplicated(scores_cancer$Patient),]
	# rownames(scores_cancer) = scores_cancer$Patient
	# sc$Stage_collapsed = scores_cancer[commonz,"Stage_collapsed"]
	# cor(sc[,"in_cancer_hlcaat0"],sc[,"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="I") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="I") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="II") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="II") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="III") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="III") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_cancer_hlcaat0"],sc[(sc$Stage_collapsed=="IV") %in% c(T),"in_at2_hlcaat0"])
	# cor(sc[,"in_cancer_hlcaat2"],sc[,"in_at2_hlcaat2"])
	# cor(sc[,"in_cancer_hlcaat1"],sc[,"in_at2_hlcaat1"])

	### meansc Trav HLCA LP:
	load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	scores = scores_uncorrected
	load(file = paste0( OutDir,"gs_TravHlcaLPHan.RData" ))
	gs = gs[!(grepl( "Trav_",names(gs) ))]
	dir.create(paste0( OutDir,"meansc_TravHlcaLPHan/" ) )
	for (thisct in c( "AT1","AT2","AT0","preTB","ciliated","club","basal" )){
		dcat(thisct)
		load(file=paste0(OutDir,"pst.RData"))
		sco = scores[scores$Epi_Cancer==thisct,]
		if (nrow(sco)==0){ 
			dcat( "No rows" ) 
			next
		}
		for (rn in rownames(pst)){
			for (tp in names(gs)){
				pst[rn,tp] = mean(sco[(sco$Patient==pst[rn,"Patient"]) & (sco$SampleType==pst[rn,"SampleType"]),tp ])
			}
		}
		pst = pst[!is.na(pst[,names(gs)[1]]),]
		pdf( paste0(OutDir,"meansc_TravHlcaLPHan/meansc_TravHlcaLPHan_",thisct,"_cMinimal_Boxplots.pdf"),7,5 )
		for (tp in names(gs)){
			x = pst$cMinimal
			y = pst[,tp]
			xColors = cMinimal_map[levels(x)]
			jitterColors = pst$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_TravHlcaLPHan/meansc_TravHlcaLPHan_",thisct,"_cMinimal_PairedDotPlot.pdf"),7,5 )
		for (tp in names(gs)){
			patients_multiple_samples = rowSums(dtable(pst$Patient,pst$cMinimal)>0)
			patients_multiple_samples = names(which(patients_multiple_samples>=2))
			tpst = pst[pst$Patient %in% patients_multiple_samples,]
			dftest = tpst[tpst$cMinimal %in% c( "Normal","Tumor" ),c( "Patient","cMinimal",tp )]
			patients_two = rowSums(dtable(dftest$Patient,dftest$cMinimal)>0)
			patients_two = names(which(patients_two==2))
			dftest = dftest[dftest$Patient %in% patients_two,]
			dftest = dftest[order(dftest$cMinimal,dftest$Patient),]
			wt = wilcox.test(x=dftest[dftest$cMinimal=="Normal",tp],y=dftest[dftest$cMinimal=="Tumor",tp],paired=TRUE)
			plotTitle = paste0( "Normal vs Tumor\nPaired Wilcoxon test, p-val = ",signif(wt$p.value,2) )
			x = tpst$cMinimal
			y = tpst[,tp]
			pairingz = tpst$Patient
			xColors = cMinimal_map[levels(x)]
			jitterColors = tpst$cMinimal_colorz
			plot=dan.pairedDotPlot.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = plotTitle, labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = jitterColors, labelPoints = "", linesPairings = pairingz )
			print(plot)
		}
		dev.off()
		pst2 = pst[pst$cMinimal!="Preinvasive",]
		pdf( paste0(OutDir,"meansc_TravHlcaLPHan/meansc_TravHlcaLPHan_",thisct,"_cMinimal_OnlyNormalVsTumor_PairedDotPlot.pdf"),5,5 )
		for (tp in names(gs)){
			patients_multiple_samples = rowSums(dtable(pst2$Patient,pst2$cMinimal)>0)
			patients_multiple_samples = names(which(patients_multiple_samples>=2))
			tpst = pst2[pst2$Patient %in% patients_multiple_samples,]
			dftest = tpst[tpst$cMinimal %in% c( "Normal","Tumor" ),c( "Patient","cMinimal",tp )]
			patients_two = rowSums(dtable(dftest$Patient,dftest$cMinimal)>0)
			patients_two = names(which(patients_two==2))
			dftest = dftest[dftest$Patient %in% patients_two,]
			dftest = dftest[order(dftest$cMinimal,dftest$Patient),]
			wt = wilcox.test(x=dftest[dftest$cMinimal=="Normal",tp],y=dftest[dftest$cMinimal=="Tumor",tp],paired=TRUE)
			plotTitle = paste0( "Normal vs Tumor\nPaired Wilcoxon test, p-val = ",signif(wt$p.value,2) )
			x = tpst$cMinimal
			y = tpst[,tp]
			pairingz = tpst$Patient
			xColors = cMinimal_map[levels(x)]
			jitterColors = tpst$cMinimal_colorz
			plot=dan.pairedDotPlot.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = plotTitle, labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = max(y), jitterColors = jitterColors, labelPoints = "", linesPairings = pairingz )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_TravHlcaLPHan/meansc_TravHlcaLPHan_",thisct,"_cMinimal_OnlyNormalVsTumor_Boxplots.pdf"),5,5 )
		for (tp in names(gs)){
			x = pst2$cMinimal
			y = pst2[,tp]
			xColors = cMinimal_map[levels(x)[levels(x) %in% unique(as.character(x))]]
			jitterColors = pst2$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_TravHlcaLPHan/meansc_TravHlcaLPHan_",thisct,"_cBroad_Boxplots.pdf"),8,5 )
		for (tp in names(gs)){
			x = pst$cBroad
			y = pst[,tp]
			xColors = cBroad_map[levels(x)]
			jitterColors = pst$cBroad_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_TravHlcaLPHan/meansc_TravHlcaLPHan_",thisct,"_cDetailed_Boxplots.pdf"),10,5 )
		for (tp in names(gs)){
			x = pst$cDetailed[!is.na(pst$cDetailed)]
			y = pst[!is.na(pst$cDetailed),tp]
			xColors = cDetailed_map[levels(x)]
			jitterColors = pst$cDetailed_colorz[!is.na(pst$cDetailed)]
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
	}

	### 2D plots with densities. see https://respiratory-research.biomedcentral.com/articles/10.1186/s12931-016-0358-z
	load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	scores = scores_uncorrected
	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	## AT1-AT2
	this_OutDir = paste0(OutDir,"DensityScatters_2D_AT1AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT12"
	these_cts = c( "AT1","AT2" )
	sco = scores[(scores$Epi_Cancer %in% these_cts) & (scores$SampleType!="Metastasis"),]
	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_AT1),ylimRight=max(sco$Han_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT1),ylimRight=max(sco$HLCA_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_KACsignature.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = sco$Han_KAC_signature, xlab = "Han_AT2", ylab = "Han_AT1", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_KACsignature.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = sco$Han_KAC_signature, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_CytoTRACE.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = annot$cyto, xlab = "Han_AT2", ylab = "Han_AT1", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_CytoTRACE.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = annot$cyto, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )

	## AT1-AT2-AT0
	this_OutDir = paste0(OutDir,"DensityScatters_2D_AT1AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	these_cts = c( "AT1","AT2","AT0" )
	sco = scores[(scores$Epi_Cancer %in% these_cts) & (scores$SampleType!="Metastasis"),]
	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_AT1),ylimRight=max(sco$Han_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT1),ylimRight=max(sco$HLCA_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_KACsignature.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = sco$Han_KAC_signature, xlab = "Han_AT2", ylab = "Han_AT1", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_KACsignature.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = sco$Han_KAC_signature, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_CytoTRACE.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = annot$cyto, xlab = "Han_AT2", ylab = "Han_AT1", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_CytoTRACE.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = annot$cyto, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )

	## AT1-AT2-AT0-Cancer
	this_OutDir = paste0(OutDir,"DensityScatters_withCancer_2D_AT1AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120Cancer"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	# downsample cancer and split by stage
	set.seed(42)
	nsampled = 10000
	ts = scores[((scores$Epi_Cancer %in% c( "Cancer" )) & (scores$SampleType!="Metastasis")),]
	ts = ts[sample(rownames(ts),nsampled),]
	sco = rbind(sco,ts)
	these_cts = c( "AT1","AT2","AT0","Cancer" )
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_AT1),ylimRight=max(sco$Han_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT1),ylimRight=max(sco$HLCA_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_KACsignature.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = sco$Han_KAC_signature, xlab = "Han_AT2", ylab = "Han_AT1", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_KACsignature.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = sco$Han_KAC_signature, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_CytoTRACE.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = annot$cyto, xlab = "Han_AT2", ylab = "Han_AT1", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_CytoTRACE.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = annot$cyto, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )

	## AT1-AT2-AT0 on AT2-KAC coordinates
	this_OutDir = paste0(OutDir,"DensityScatters_2D_KacAT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	these_cts = c( "AT1","AT2","AT0" )
	sco = scores[(scores$Epi_Cancer %in% these_cts) & (scores$SampleType!="Metastasis"),]
	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_KAC_signature, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==3))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_KAC_signature, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_KAC_signature),ylimRight=max(sco$Han_KAC_signature),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_CytoTRACE.pdf"), x=sco$Han_AT2, y=sco$Han_KAC_signature, fill = annot$cyto, xlab = "Han_AT2", ylab = "Han_KAC_signature", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_KAC_signature), ylimRight = max(sco$Han_KAC_signature), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	
	## AT1-AT2-AT0-Cancer on AT2-KAC coordinates
	this_OutDir = paste0(OutDir,"DensityScatters_withCancer_2D_KacAT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120Cancer"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	# downsample cancer and split by stage
	set.seed(42)
	nsampled = 10000
	ts = scores[((scores$Epi_Cancer %in% c( "Cancer" )) & (scores$SampleType!="Metastasis")),]
	ts = ts[sample(rownames(ts),nsampled),]
	sco = rbind(sco,ts)
	these_cts = c( "AT1","AT2","AT0","Cancer" )
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_KAC_signature, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_KAC_signature, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_KAC_signature),ylimRight=max(sco$Han_KAC_signature),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_CytoTRACE.pdf"), x=sco$Han_AT2, y=sco$Han_KAC_signature, fill = annot$cyto, xlab = "Han_AT2", ylab = "Han_KAC_signature", 
		xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_KAC_signature), ylimRight = max(sco$Han_KAC_signature), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )



	## AT1-AT2-AT0 on AT2-AT0 coordinates
	this_OutDir = paste0(OutDir,"DensityScatters_2D_AT0AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	these_cts = c( "AT1","AT2","AT0" )
	sco = scores[(scores$Epi_Cancer %in% these_cts) & (scores$SampleType!="Metastasis"),]
	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT0, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==3))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT0, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT0),ylimRight=max(sco$HLCA_AT0),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_CytoTRACE.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT0, fill = annot$cyto, xlab = "HLCA_AT2", ylab = "HLCA_AT0", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT0), ylimRight = max(sco$HLCA_AT0), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	
	## AT1-AT2-AT0-Cancer on AT2-AT0 coordinates
	this_OutDir = paste0(OutDir,"DensityScatters_withCancer_2D_AT0AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120Cancer"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	# downsample cancer and split by stage
	set.seed(42)
	nsampled = 10000
	ts = scores[((scores$Epi_Cancer %in% c( "Cancer" )) & (scores$SampleType!="Metastasis")),]
	ts = ts[sample(rownames(ts),nsampled),]
	sco = rbind(sco,ts)
	these_cts = c( "AT1","AT2","AT0","Cancer" )
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT0, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT0, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=factor(scop$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT0),ylimRight=max(sco$HLCA_AT0),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_CytoTRACE.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT0, fill = annot$cyto, xlab = "HLCA_AT2", ylab = "HLCA_AT0", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT0), ylimRight = max(sco$HLCA_AT0), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )




	## Then, 2D plots with tps (and tumor cells)
	load(file = paste0(OutDir,"../scoring/EANregrAT2_tps_scores_DatasetRegressed.RData"))
	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	scores[scores$SampleType %in% c( "AAH","AIS","MIA" ),"cMinimal"] = "Preinvasive"
	## AT2like-Interferon
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2likeInterferon/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	sco = sco[sco$cMinimal!="Preinvasive",]
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_like, y=sco$Interferon, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	## AT2_Club_like-Interferon
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2ClublikeInterferon/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	sco = sco[sco$cMinimal!="Preinvasive",]
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_Club_like, y=sco$Interferon, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	## AT2like-Unfolded_protein_response
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2likeUnfoldedProteinResponse/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	sco = sco[sco$cMinimal!="Preinvasive",]
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_like, y=sco$Unfolded_protein_response, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	## AT2like-MHC_II
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2likeMHCII/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	sco = sco[sco$cMinimal!="Preinvasive",]
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_like, y=sco$MHC_II, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=c( "Normal","Tumor" )), fillColors2=c( "gray44","firebrick" ), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)

	load(file = paste0(OutDir,"gs_scores_extendedAtlasNormal.RData"))
	scores = scores_uncorrected
	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	load(file=paste0(OutDir,"pst.RData"))
	for (rn in rownames(pst)){
		sc = scores[(scores$Patient==pst[rn,"Patient"]) & (scores$SampleType==pst[rn,"SampleType"]), ]
		pst[rn,"mean_AT1_in_AT1"] = mean(sc[sc$Epi_Cancer=="AT1","HLCA_AT1"])
		pst[rn,"median_AT1_in_AT1"] = median(sc[sc$Epi_Cancer=="AT1","HLCA_AT1"])
		pst[rn,"mean_AT2_in_AT2"] = mean(sc[sc$Epi_Cancer=="AT2","HLCA_AT2"])
		pst[rn,"median_AT2_in_AT2"] = median(sc[sc$Epi_Cancer=="AT2","HLCA_AT2"])
	}
	pst[pst$n_AT1<5,"mean_AT1_in_AT1"] = NA
	pst[pst$n_AT1<5,"median_AT1_in_AT1"] = NA
	pst[pst$n_AT2<5,"mean_AT2_in_AT2"] = NA
	pst[pst$n_AT2<5,"median_AT2_in_AT2"] = NA
	pst = pst[pst$Patient %in% unique(pst[pst$SampleType=="Normal","Patient"]),]
	tt = dtable(pst$SampleType,pst$Patient)
	tt = colnames(tt)[which(colSums(tt)>1)]
	pst = pst[pst$Patient %in% tt,]
	pstn = pst[pst$SampleType=="Normal",]
	rownames(pstn) = tt
	pstt = pst[pst$SampleType!="Normal",]
	rownames(pstt) = tt
	pt_at1 = data.frame( row.names=tt, diff_mean=pstt[tt,"mean_AT1_in_AT1"]-pstn[tt,"mean_AT1_in_AT1"], diff_median=pstt[tt,"median_AT1_in_AT1"]-pstn[tt,"median_AT1_in_AT1"] )
	pt_at1 = pt_at1[!is.na(rowSums(pt_at1)),]
	dcat( paste0( "mean HLCA_AT1 in AT1 decreased from normal to tumor in ",signif(100*sum(pt_at1$diff_mean<0)/nrow(pt_at1),2),"% of patients" ) )
	dcat( paste0( "median HLCA_AT1 in AT1 decreased from normal to tumor in ",signif(100*sum(pt_at1$diff_median<0)/nrow(pt_at1),2),"% of patients" ) )
	pt_at2 = data.frame( row.names=tt, diff_mean=pstt[tt,"mean_AT2_in_AT2"]-pstn[tt,"mean_AT2_in_AT2"], diff_median=pstt[tt,"median_AT2_in_AT2"]-pstn[tt,"median_AT2_in_AT2"] )
	pt_at2 = pt_at2[!is.na(rowSums(pt_at2)),]
	dcat( paste0( "mean HLCA_AT2 in AT2 decreased from normal to tumor in ",signif(100*sum(pt_at2$diff_mean<0)/nrow(pt_at2),2),"% of patients" ) )
	dcat( paste0( "median HLCA_AT2 in AT2 decreased from normal to tumor in ",signif(100*sum(pt_at2$diff_median<0)/nrow(pt_at2),2),"% of patients" ) )	
}

tps_atlases_NormalCells_withPreinvasive = function( OutDir, order_tps ){
	load(paste0(OutDir,"../CellStates/tps_CellStates_hierarchical_2_substates12.RData"))
	scores[scores$cs=="cs1","cs"] = "Alveolar"
	scores[scores$cs=="cs2","cs"] = "Proliferative"
	scores[scores$cs=="cs3","cs"] = "Hypoxic"
	cs_scores = scores
	celltype_map = data.frame(row.names=c( "AT1","AT2","AT0" ),
		colorz = c("chocolate","dodgerblue2","purple") ,stringsAsFactors=F )
	cMinimal_map = c( "gray44","orange","firebrick" )
	names(cMinimal_map) = c( "Normal","Preinvasive","Tumor" )
	cBroad_map = c( "gray44","dodgerblue2","steelblue4","gold","firebrick" )
	names(cBroad_map) = c( "Normal","AAH", "AIS", "MIA","Tumor" )
	cDetailed_map = c( "gray44","dodgerblue2","steelblue4","gold","orange","tomato3","firebrick","brown" )
	names(cDetailed_map) = c( "Normal","AAH", "AIS", "MIA","I","II","III","IV")
	DataDir = "data/"
	CellTypeDir = paste0("data/cell_typing/non_malignant/")
	load(file=paste0(OutDir,"asl.RData"))
	load(file=paste0(OutDir,"pst.RData"))
	load(paste0(OutDir,"../tps_discovery/tps.RData"))
	load(paste0(OutDir,"../tps_discovery/tps_universe.RData"))
	tps = tps[order_tps]
	load(file = paste0(OutDir,"gs_scores_extendedAtlasNormalS0.RData"))
	scores = scores_uncorrected

	### AT2-like density and tissue of origin
	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	scores[scores$SampleType %in% c( "AAH","AIS","MIA" ),"cMinimal"] = "Preinvasive"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco = sco[(sco$cMinimal %in% c( "Normal","Tumor" )),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Tumor" ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT0","AT2" ))
	plot = ggplot(sco, aes(x=`AT2-like`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="AT2-like scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_AT2likescores_splitcMinimal.pdf"),2.7,2 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_MHCIIscores_splitcMinimal.pdf"),2.7,2 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`AT2-Club-like`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="AT0-like scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_AT0likescores_splitcMinimal.pdf"),2.7,2 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`Han_KAC_signature`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_AT0likescores_splitcMinimal.pdf"),2.7,2 )
	print(plot)
	dev.off()

	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	scores[scores$SampleType %in% c( "AAH","AIS","MIA" ),"cMinimal"] = "Preinvasive"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco = sco[(sco$cMinimal %in% c( "Normal","Preinvasive","Tumor" )),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Preinvasive","Tumor" ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT0","AT2" ))
	plot = ggplot(sco, aes(x=`HLCA_AT2`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=c("gray44","orange","#1C75BC" )) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_hlcaAT2_withPreinvasive_splitcMinimal.pdf"),2.7,2 )
	print(plot)
	dev.off()

	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	scores[scores$SampleType %in% c( "AAH","AIS","MIA" ),"cMinimal"] = "Preinvasive"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco = sco[(sco$cMinimal %in% c( "Normal","Preinvasive","Tumor" )),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Preinvasive","Tumor" ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT0","AT2" ))
	plot = ggplot(sco, aes(x=`HLCA_AT1`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=c("gray44","orange","#1C75BC" )) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_hlcaAT1_withPreinvasive_splitcMinimal.pdf"),2.7,2 )
	print(plot)
	dev.off()

	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	scores[scores$SampleType %in% c( "AAH","AIS","MIA" ),"cMinimal"] = "Preinvasive"
	sco = scores[(scores$Epi_Cancer %in% c( "AT2","AT0","Cancer" )) & (scores$SampleType!="Metastasis"),]
	sco = sco[(sco$cMinimal %in% c( "Normal","Tumor" )),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Tumor" ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT0","AT2","Cancer" ))
	plot = ggplot(sco, aes(x=`Han_KAC_signature`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=c("gray44","#1C75BC" )) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_KAC_splitcMinimal.pdf"),3,2 )
	print(plot)
	dev.off()

	scores$cMinimal = NA
	scores[scores$Epi_Cancer %in% c( "Cancer" ),"cMinimal"] = "Tumor"
	scores[scores$Epi_Cancer %in% c( "AT2","AT0" ),"cMinimal"] = "Normal\nepithelial"
	sco = scores[((scores$Epi_Cancer %in% c( "AT2","AT0","Cancer" )) & (scores$SampleType %in% c("Normal","Primary") )) & (!is.na(scores$Stage_collapsed)),]
	sco = sco[(sco$cMinimal %in% c( "Normal\nepithelial","Tumor" )),]
	sco[sco$cMinimal=="Tumor","cMinimal"] = paste0( "Cancer,\nstage ",sco[sco$cMinimal=="Tumor","Stage_collapsed"] )
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal\nepithelial","Cancer,\nstage I","Cancer,\nstage II","Cancer,\nstage III","Cancer,\nstage IV" ))
	colorz = c( "gray44","#C4DDEE","#8BBBDD","#4F98CC","#0D76BC" )
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT0","AT2","Cancer" ))
	plot = ggplot(sco, aes(x=`Han_KAC_signature`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=colorz) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_KAC_splitStage.pdf"),3,2 )
	print(plot)
	dev.off()

	scores$cMinimal = NA
	scores[scores$Epi_Cancer %in% c( "Cancer" ),"cMinimal"] = "Tumor"
	scores[scores$Epi_Cancer %in% c( "AT2","AT0" ),"cMinimal"] = "Normal\nepithelial"
	sco = scores[((scores$Epi_Cancer %in% c( "AT2","AT0","Cancer" )) & (scores$SampleType %in% c("Normal","Primary") )),]
	sco = sco[(sco$cMinimal %in% c( "Normal\nepithelial","Tumor" )),]
	sco$CellState = NA
	sco[intersect(cs_scores[cs_scores$cs=="Alveolar","CellID"],rownames(sco)),"CellState"] = "Alveolar-like"
	sco[intersect(cs_scores[cs_scores$cs=="Proliferative","CellID"],rownames(sco)),"CellState"] = "Early-DD"
	sco[intersect(cs_scores[cs_scores$cs=="Hypoxic","CellID"],rownames(sco)),"CellState"] = "Advanced-DD"
	sco[sco$cMinimal=="Tumor","cMinimal"] = paste0( "Cancer,\n",sco[sco$cMinimal=="Tumor","CellState"] )
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal\nepithelial","Cancer,\nAlveolar-like","Cancer,\nEarly-DD","Cancer,\nAdvanced-DD" ))
	colorz = c( "gray44","#2B3990","#FBB040","#BE1E2D" )
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT0","AT2","Cancer" ))
	plot = ggplot(sco, aes(x=`Han_KAC_signature`,fill=cMinimal)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=colorz) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_KAC_splitStates.pdf"),3,2 )
	print(plot)
	dev.off()

	

	### meansc tps:
	dir.create(paste0( OutDir,"meansc_tps/" ) )
	for (thisct in c( "AT1","AT2","AT0" )){
		dcat(thisct)
		load(file=paste0(OutDir,"pst.RData"))
		sco = scores[scores$Epi_Cancer==thisct,]
		if (nrow(sco)==0){ 
			dcat( "No rows" ) 
			next
		}
		for (rn in rownames(pst)){
			for (tp in names(tps)){
				pst[rn,tp] = mean(sco[(sco$Patient==pst[rn,"Patient"]) & (sco$SampleType==pst[rn,"SampleType"]),tp ])
			}
		}
		pst = pst[!is.na(pst[,names(tps)[1]]),]
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cMinimal_Boxplots.pdf"),7,5 )
		for (tp in names(tps)){
			x = pst$cMinimal
			y = pst[,tp]
			xColors = cMinimal_map[levels(x)]
			jitterColors = pst$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cBroad_Boxplots.pdf"),7,5 )
		for (tp in names(tps)){
			x = pst$cBroad
			y = pst[,tp]
			xColors = cBroad_map[levels(x)]
			jitterColors = pst$cBroad_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_tps/meansc_tps_",thisct,"_cDetailed_Boxplots.pdf"),10,5 )
		for (tp in names(tps)){
			x = pst$cDetailed[!is.na(pst$cDetailed)]
			y = pst[!is.na(pst$cDetailed),tp]
			xColors = cDetailed_map[levels(x)]
			jitterColors = pst$cDetailed_colorz[!is.na(pst$cDetailed)]
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
	}

	### meansc Trav HLCA LP:
	load(file = paste0( OutDir,"gs_TravHlcaLPHan.RData" ))
	gs = gs[ c("HLCA_AT0","HLCA_AT1","HLCA_AT2","LP_proliferation","Han_AT2","Han_AT1","Han_KAC_signature") ]
	dir.create(paste0( OutDir,"meansc_HlcaLPHan/" ) )
	for (thisct in c( "AT1","AT2","AT0" )){
		dcat(thisct)
		load(file=paste0(OutDir,"pst.RData"))
		sco = scores[scores$Epi_Cancer==thisct,]
		if (nrow(sco)==0){ 
			dcat( "No rows" ) 
			next
		}
		for (rn in rownames(pst)){
			for (tp in names(gs)){
				pst[rn,tp] = mean(sco[(sco$Patient==pst[rn,"Patient"]) & (sco$SampleType==pst[rn,"SampleType"]),tp ])
			}
		}
		pst = pst[!is.na(pst[,names(gs)[1]]),]
		pdf( paste0(OutDir,"meansc_HlcaLPHan/meansc_HlcaLPHan",thisct,"_cMinimal_Boxplots.pdf"),7,5 )
		for (tp in names(gs)){
			x = pst$cMinimal
			y = pst[,tp]
			xColors = cMinimal_map[levels(x)]
			jitterColors = pst$cMinimal_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_HlcaLPHan/meansc_HlcaLPHan",thisct,"_cBroad_Boxplots.pdf"),8,5 )
		for (tp in names(gs)){
			x = pst$cBroad
			y = pst[,tp]
			xColors = cBroad_map[levels(x)]
			jitterColors = pst$cBroad_colorz
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
		pdf( paste0(OutDir,"meansc_HlcaLPHan/meansc_HlcaLPHan",thisct,"_cDetailed_Boxplots.pdf"),10,5 )
		for (tp in names(gs)){
			x = pst$cDetailed[!is.na(pst$cDetailed)]
			y = pst[!is.na(pst$cDetailed),tp]
			xColors = cDetailed_map[levels(x)]
			jitterColors = pst$cDetailed_colorz[!is.na(pst$cDetailed)]
			plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", ylimLeft = min(y), ylimRight = max(y),comparisons = NULL, labelycoo = max(y), xColors = xColors, jitterColors = jitterColors, labelJitteredPoints = NULL, includeJitters = T )
			print(plot)
		}
		dev.off()
	}

	### 2D plots with densities. see https://respiratory-research.biomedcentral.com/articles/10.1186/s12931-016-0358-z
	
	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	scores[scores$SampleType %in% c( "AAH","AIS","MIA" ),"cMinimal"] = "Preinvasive"

	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Preinvasive","Tumor" ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT2","AT0" ))
	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HlcaAT2scores.pdf"),2.5,1.5 )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT1,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HlcaAT1scores.pdf"),2.5,1.5 )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_Interferon_scores.pdf"),2.5,1.5 )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_MHCII_scores.pdf"),2.5,1.5 )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Han_KAC_signature,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HanKacSignature_scores.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=LP_proliferation,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="LP_proliferation scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_LP_proliferation_scores.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT0,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA_AT0 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HlcaAT0_scores.pdf"),7,4 )
	print(plot)
	dev.off()
	dan.save(scores,paste0(OutDir,"densities_singlecells_HlcaAT2scores.pdf"))

	load(paste0(OutDir,"densities_singlecells_HlcaAT2scores.RData"))
	sco = scores[(scores$Epi_Cancer %in% c( "AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco$SampleType[sco$SampleType=="Primary"] = "Tumor"
	sco$cBroad = factor(sco$SampleType,levels=c( names(cBroad_map) ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT2","AT0" ))
	ccolorz = c("gray44","gold","darkorange","chocolate3","#1C75BC" )
	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cBroad)) + geom_density(alpha=0.7,lwd=0.25) + scale_fill_manual(values=ccolorz) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_hlcaAT2scores_splitcBroad.pdf"),2.7,2 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cBroad)) + geom_density(alpha=0.7,lwd=0.25)+coord_cartesian(xlim=c(min(sco$Interferon), 1.2)) + scale_fill_manual(values=ccolorz) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.position="top", legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'),axis.line=element_line(size=0.25),axis.ticks = element_line(size = 0.25))
	pdf( paste0(OutDir,"densities_singlecells_interferonscores_splitcBroad.pdf"),2.7,2 )
	print(plot)
	dev.off()

	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Tumor" ))
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT2","AT0" ))
	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HlcaAT2scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT1,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HlcaAT1scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_Interferon_scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_MHCII_scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Han_KAC_signature,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_HanKacSignature_scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()

	### with cancer
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Preinvasive","Tumor" ))
	set.seed(123)
	nsampled = 20000
	ts = scores[((scores$Epi_Cancer %in% c( "Cancer" )) & (scores$SampleType=="Primary")),]
	ts = ts[sample(rownames(ts),nsampled),]
	ts[,"Epi_Cancer"] = cs_scores[rownames(ts),"cs"]
	sco = rbind(sco,ts)
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT2","AT0","Alveolar","Proliferative","Hypoxic" ))
	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_HlcaAT2scores.pdf"),2.5,3 )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT1,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_HlcaAT1scores.pdf"),2.5,3  )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_Interferon_scores.pdf"),2.5,3  )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_MHCII_scores.pdf"),2.5,3  )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Han_KAC_signature,fill=cMinimal)) + geom_density(alpha=0.7,size=0.2) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_HanKacSignature_scores.pdf"),2.5,3  )
	plot = plot + theme_classic(base_size=6)+ theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()

	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Tumor" ))
	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_HlcaAT2scores_noPreinvasive.pdf"),7,5.5 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT1,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_HlcaAT1scores_noPreinvasive.pdf"),7,5.5 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_Interferon_scores_noPreinvasive.pdf"),7,5.5 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_MHCII_scores_noPreinvasive.pdf"),7,5.5 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Han_KAC_signature,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_HanKacSignature_scores_noPreinvasive.pdf"),7,5.5 )
	print(plot)
	dev.off()

	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Preinvasive","Tumor" ))
	set.seed(123)
	nsampled = 10000
	ts = scores[((scores$Epi_Cancer %in% c( "Cancer" )) & (scores$SampleType=="Primary")),]
	ts = ts[sample(rownames(ts),nsampled),]
	ts[,"Epi_Cancer"] = cs_scores[rownames(ts),"cs"]
	sco = rbind(sco,ts)
	sco$Epi_Cancer = as.character(sco$Epi_Cancer)
	sco[sco$Epi_Cancer %in% c( "Alveolar","Proliferative","Hypoxic" ),"Epi_Cancer"] = "Cancer"
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=c( "AT1","AT2","AT0","Cancer" ))

	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_HlcaAT2scores.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT1,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_HlcaAT1scores.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_Interferon_scores.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_MHCII_scores.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Han_KAC_signature,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_HanKacSignature_scores.pdf"),7,4 )
	print(plot)
	dev.off()

	sco = sco[sco$cMinimal!="Preinvasive",]
	sco$cMinimal = factor(sco$cMinimal,levels=c( "Normal","Tumor" ))
	plot = ggplot(sco, aes(x=HLCA_AT2,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT2 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_HlcaAT2scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=HLCA_AT1,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="HLCA AT1 scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_HlcaAT1scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Interferon,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="Interferon scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_Interferon_scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=`MHC-II`,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="MHC-II scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_MHCII_scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()
	plot = ggplot(sco, aes(x=Han_KAC_signature,fill=cMinimal)) + geom_density(alpha=0.7) + scale_fill_manual(values=cMinimal_map[levels(sco$cMinimal)]) + facet_grid(Epi_Cancer ~ .) + labs(x="KAC signature scores",y="Density",fill="Sample type") + theme_classic()
	pdf( paste0(OutDir,"densities_singlecells_WithCancer_NoCs_HanKacSignature_scores_noPreinvasive.pdf"),7,4 )
	print(plot)
	dev.off()

	## AT1-AT2
	this_OutDir = paste0(OutDir,"DensityScatters_2D_AT1AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT12"
	these_cts = c( "AT1","AT2" )
	sco = scores[(scores$Epi_Cancer %in% these_cts) & (scores$SampleType!="Metastasis"),]
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=3, fileHeight=3)
	# dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_AT1, 
	# 	fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
	# 	fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
	# 	splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	# dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	# for (p in unique(sco$Patient)){
	# 	scop = sco[sco$Patient==p,]
	# 	scop$cMinimal = factor(scop$cMinimal,levels=intersect(names(cMinimal_map),unique(scop$cMinimal)) )
	# 	tt = dtable(scop$SampleType,scop$Epi_Cancer)
	# 	if ( any(scop$SampleType %in% c( "AAH","AIS","MIA" )) ){
	# 		dcat(p)
	# 		dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_AT1, 
	# 			fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=scop$cMinimal, fillColors2=as.character(cMinimal_map[levels(scop$cMinimal)]), filllab2="cMinimal",
	# 			xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_AT1),ylimRight=max(sco$Han_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
	# 	}
	# }
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		scop$cMinimal = factor(scop$cMinimal,levels=intersect(names(cMinimal_map),unique(scop$cMinimal)) )
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=scop$cMinimal, fillColors2=as.character(cMinimal_map[levels(scop$cMinimal)]), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT1),ylimRight=max(sco$HLCA_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	# dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_KACsignature.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = sco$Han_KAC_signature, xlab = "Han_AT2", ylab = "Han_AT1", 
	# 	xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	these_colorz = colorRampPalette(c("dodgerblue4","dodgerblue3","white","firebrick3","firebrick4"))(100)
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_KACsignature.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = sco$Han_KAC_signature, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "KAC_signature", fillColors_continuous = these_colorz, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )

	## AT1-AT2-AT0
	this_OutDir = paste0(OutDir,"DensityScatters_2D_AT1AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	these_cts = c( "AT1","AT2","AT0" )
	sco = scores[(scores$Epi_Cancer %in% these_cts) & (scores$SampleType!="Metastasis"),]
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.5, hline=0.5,vline=0.5,fileWidth=2.5, fileHeight=2.5)
	# dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han"), x=sco$Han_AT2, y=sco$Han_AT1, 
	# 	fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
	# 	fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
	# 	splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	# dir.create(paste0(this_OutDir,file_prefix,"_Han_PatientWise"))
	# for (p in unique(sco$Patient)){
	# 	scop = sco[sco$Patient==p,]
	# 	scop$cMinimal = factor(scop$cMinimal,levels=intersect(names(cMinimal_map),unique(scop$cMinimal)) )
	# 	tt = dtable(scop$SampleType,scop$Epi_Cancer)
	# 	if ( any(scop$SampleType %in% c( "AAH","AIS","MIA" )) ){
	# 		dcat(p)
	# 		dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Han_PatientWise/",p), x=scop$Han_AT2, y=scop$Han_AT1, 
	# 			fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=scop$cMinimal, fillColors2=as.character(cMinimal_map[levels(scop$cMinimal)]), filllab2="cMinimal",
	# 			xlimLeft=min(sco$Han_AT2),xlimRight=max(sco$Han_AT2),ylimLeft=min(sco$Han_AT1),ylimRight=max(sco$Han_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
	# 	}
	# }
	dir.create(paste0(this_OutDir,file_prefix,"_Hlca_PatientWise"))
	for (p in unique(sco$Patient)){ 
		scop = sco[sco$Patient==p,]
		scop$cMinimal = factor(scop$cMinimal,levels=intersect(names(cMinimal_map),unique(scop$cMinimal)) )
		tt = dtable(scop$SampleType,scop$Epi_Cancer)
		if ((all(tt>=5)) & ((nrow(tt)==2) & (ncol(tt)==2))){
			dcat(p)
			dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca_PatientWise/",p), x=scop$HLCA_AT2, y=scop$HLCA_AT1, 
				fill1=scop$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", fill2=scop$cMinimal, fillColors2=as.character(cMinimal_map[levels(scop$cMinimal)]), filllab2="cMinimal",
				xlimLeft=min(sco$HLCA_AT2),xlimRight=max(sco$HLCA_AT2),ylimLeft=min(sco$HLCA_AT1),ylimRight=max(sco$HLCA_AT1),splitByFill1=TRUE,splitByFill2=TRUE,dotSize=0.1,fileWidth=7,fileHeight=7)
		}
	}
	# dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_KACsignature.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = sco$Han_KAC_signature, xlab = "Han_AT2", ylab = "Han_AT1", 
	# 	xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_KACsignature.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = sco$Han_KAC_signature, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "KAC_signature", fillColors_continuous = these_colorz, dotSize = 0.01, coord_fixed = TRUE, fileWidth = 2.5, fileHeight = 2.5 )

	## AT1-AT2-AT0-Cancer
	celltype_map = data.frame(row.names=c( "AT1","AT2","AT0","Alveolar","Proliferative","Hypoxic","Cancer" ),
		colorz = c("chocolate","dodgerblue2","purple","dodgerblue4","firebrick3","sienna4","firebrick") ,stringsAsFactors=F )
	this_OutDir = paste0(OutDir,"DensityScatters_withCancer_2D_AT1AT2/")
	dir.create(this_OutDir)
	file_prefix = "AT120Cancer"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	# downsample cancer and split by stage
	set.seed(123)
	nsampled = 20000
	ts = scores[((scores$Epi_Cancer %in% c( "Cancer" )) & (scores$SampleType=="Primary")),]
	ts = ts[sample(rownames(ts),nsampled),]
	ts[,"Epi_Cancer"] = cs_scores[rownames(ts),"cs"]
	sco = rbind(sco,ts)
	these_cts = c( "AT1","AT2","AT0","Cancer" )
	these_cts = c( "AT1","AT2","AT0","Alveolar","Proliferative","Hypoxic" )
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.3, fileWidth=2.5, fileHeight=2.5)
	# dan.scatterplot( paste0(this_OutDir,file_prefix,"_Han_KACsignature.pdf"), x=sco$Han_AT2, y=sco$Han_AT1, fill = sco$Han_KAC_signature, xlab = "Han_AT2", ylab = "Han_AT1", 
		# xlimLeft = min(sco$Han_AT2), xlimRight = max(sco$Han_AT2), ylimLeft = min(sco$Han_AT1), ylimRight = max(sco$Han_AT1), filllab = "KAC_signature", fillColors_continuous = colorz_solid, dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_KACsignature.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = sco$Han_KAC_signature, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "KAC_signature", fillColors_continuous = these_colorz, dotSize = 0.01, coord_fixed = TRUE, fileWidth = 2.5, fileHeight = 2.5 )
	# color by cyto score
	library(viridis)
	load( file=paste0(OutDir,"../tps_vs_cytotrace/annot_extendedAtlasNormal_Cyto.RData"))
	annot = annot[rownames(sco),]
	dan.scatterplot( paste0(this_OutDir,file_prefix,"_Hlca_CytoTRACE.pdf"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, fill = annot$cyto, xlab = "HLCA_AT2", ylab = "HLCA_AT1", 
		xlimLeft = min(sco$HLCA_AT2), xlimRight = max(sco$HLCA_AT2), ylimLeft = min(sco$HLCA_AT1), ylimRight = max(sco$HLCA_AT1), filllab = "cyto score", fillColors_continuous = viridis(100), dotSize = 0.1, coord_fixed = TRUE, fileWidth = 7, fileHeight = 7 )
	sco$Epi_Cancer = as.character(sco$Epi_Cancer)	
	sco[sco$Epi_Cancer %in% c( "Alveolar","Proliferative","Hypoxic" ),"Epi_Cancer"] = "Cancer"
	these_cts = c( "AT1","AT2","AT0","Cancer" )
	sco$Epi_Cancer = factor(sco$Epi_Cancer,levels=these_cts)
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"NoCs_Hlca"), x=sco$HLCA_AT2, y=sco$HLCA_AT1, 
		fill1=sco$Epi_Cancer, fillColors1=celltype_map[these_cts,"colorz"], filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)


	## Then, 2D plots with tps (and tumor cells)
	## AT2like-Interferon
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2likeInterferon/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_like, y=sco$Interferon, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	## AT2_Club_like-Interferon
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2ClublikeInterferon/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	sco = sco[sco$cMinimal!="Preinvasive",]
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_Club_like, y=sco$Interferon, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	## AT2like-Unfolded_protein_response
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2likeUnfoldedProteinResponse/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_like, y=sco$Unfolded_protein_response, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)
	## AT2like-MHC_II
	this_OutDir = paste0(OutDir,"tps_DensityScatters_2D_AT2likeMHCII/")
	dir.create(this_OutDir)
	file_prefix = "AT120"
	sco = scores[(scores$Epi_Cancer %in% c( "AT1","AT2","AT0" )) & (scores$SampleType!="Metastasis"),]
	colnames(sco) = gsub("-","_",colnames(sco) )
	sco = sco[sco$cMinimal!="Preinvasive",]
	dan.ScatterDensity_2Dplot_2fills( fileNameRoot=paste0(this_OutDir,file_prefix,"_tpsNoPreinvasive"), x=sco$AT2_like, y=sco$MHC_II, 
		fill1=sco$Epi_Cancer, fillColors1=NULL, filllab1="Epi_Cancer", 
		fill2=factor(sco$cMinimal,levels=names(cMinimal_map)), fillColors2=as.character(cMinimal_map), filllab2="cMinimal",
		splitByFill1=TRUE,splitByFill2=TRUE, dotSize=0.1, fileWidth=7, fileHeight=7)

	load(file = paste0(OutDir,"gs_scores_extendedAtlasNormalS0.RData"))
	scores = scores_uncorrected
	nLimit = 5
	scores$cMinimal = NA
	scores[scores$SampleType %in% c( "Normal" ),"cMinimal"] = "Normal"
	scores[scores$SampleType %in% c( "Primary" ),"cMinimal"] = "Tumor"
	load(file=paste0(OutDir,"pst.RData"))
	for (rn in rownames(pst)){
		sc = scores[(scores$Patient==pst[rn,"Patient"]) & (scores$SampleType==pst[rn,"SampleType"]), ]
		pst[rn,"mean_AT1_in_AT1"] = mean(sc[sc$Epi_Cancer=="AT1","HLCA_AT1"])
		pst[rn,"median_AT1_in_AT1"] = median(sc[sc$Epi_Cancer=="AT1","HLCA_AT1"])
		pst[rn,"mean_AT2_in_AT2"] = mean(sc[sc$Epi_Cancer=="AT2","HLCA_AT2"])
		pst[rn,"median_AT2_in_AT2"] = median(sc[sc$Epi_Cancer=="AT2","HLCA_AT2"])
	}
	pst[pst$n_AT1<nLimit,"mean_AT1_in_AT1"] = NA
	pst[pst$n_AT1<nLimit,"median_AT1_in_AT1"] = NA
	pst[pst$n_AT2<nLimit,"mean_AT2_in_AT2"] = NA
	pst[pst$n_AT2<nLimit,"median_AT2_in_AT2"] = NA
	pst = pst[pst$Patient %in% unique(pst[pst$SampleType=="Normal","Patient"]),]
	tt = dtable(pst$SampleType,pst$Patient)
	tt = colnames(tt)[which(colSums(tt)>1)]
	pst = pst[pst$Patient %in% tt,]
	pstn = pst[pst$SampleType=="Normal",]
	rownames(pstn) = tt
	pstt = pst[pst$SampleType=="Primary",]
	rownames(pstt) = pstt$Patient
	pt_at1 = data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT1_in_AT1"]-pstn[rownames(pstt),"mean_AT1_in_AT1"], diff_median=pstt[rownames(pstt),"median_AT1_in_AT1"]-pstn[rownames(pstt),"median_AT1_in_AT1"] )
	pt_at2 = data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT2_in_AT2"]-pstn[rownames(pstt),"mean_AT2_in_AT2"], diff_median=pstt[rownames(pstt),"median_AT2_in_AT2"]-pstn[rownames(pstt),"median_AT2_in_AT2"] )
	pstt = pst[pst$SampleType=="AAH",]
	rownames(pstt) = pstt$Patient
	pt_at1 = rbind(pt_at1,data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT1_in_AT1"]-pstn[rownames(pstt),"mean_AT1_in_AT1"], diff_median=pstt[rownames(pstt),"median_AT1_in_AT1"]-pstn[rownames(pstt),"median_AT1_in_AT1"] ))
	pt_at2 = rbind(pt_at2,data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT2_in_AT2"]-pstn[rownames(pstt),"mean_AT2_in_AT2"], diff_median=pstt[rownames(pstt),"median_AT2_in_AT2"]-pstn[rownames(pstt),"median_AT2_in_AT2"] ))
	pstt = pst[pst$SampleType=="AIS",]
	rownames(pstt) = pstt$Patient
	pt_at1 = rbind(pt_at1,data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT1_in_AT1"]-pstn[rownames(pstt),"mean_AT1_in_AT1"], diff_median=pstt[rownames(pstt),"median_AT1_in_AT1"]-pstn[rownames(pstt),"median_AT1_in_AT1"] ))
	pt_at2 = rbind(pt_at2,data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT2_in_AT2"]-pstn[rownames(pstt),"mean_AT2_in_AT2"], diff_median=pstt[rownames(pstt),"median_AT2_in_AT2"]-pstn[rownames(pstt),"median_AT2_in_AT2"] ))
	pstt = pst[pst$SampleType=="MIA",]
	rownames(pstt) = pstt$Patient
	pt_at1 = rbind(pt_at1,data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT1_in_AT1"]-pstn[rownames(pstt),"mean_AT1_in_AT1"], diff_median=pstt[rownames(pstt),"median_AT1_in_AT1"]-pstn[rownames(pstt),"median_AT1_in_AT1"] ))
	pt_at2 = rbind(pt_at2,data.frame( row.names=pstt$Patient, diff_mean=pstt[rownames(pstt),"mean_AT2_in_AT2"]-pstn[rownames(pstt),"mean_AT2_in_AT2"], diff_median=pstt[rownames(pstt),"median_AT2_in_AT2"]-pstn[rownames(pstt),"median_AT2_in_AT2"] ))
	pt_at1 = pt_at1[!is.na(rowSums(pt_at1)),]
	dcat( paste0( "mean HLCA_AT1 in AT1 decreased from normal to tumor in ",signif(100*sum(pt_at1$diff_mean<0)/nrow(pt_at1),2),"% of patients" ) )
	dcat( paste0( "median HLCA_AT1 in AT1 decreased from normal to tumor in ",signif(100*sum(pt_at1$diff_median<0)/nrow(pt_at1),2),"% of patients" ) )
	pt_at2 = pt_at2[!is.na(rowSums(pt_at2)),]
	dcat( paste0( "mean HLCA_AT2 in AT2 decreased from normal to tumor in ",signif(100*sum(pt_at2$diff_mean<0)/nrow(pt_at2),2),"% of patients" ) )
	dcat( paste0( "median HLCA_AT2 in AT2 decreased from normal to tumor in ",signif(100*sum(pt_at2$diff_median<0)/nrow(pt_at2),2),"% of patients" ) )
	save(pt_at2,file=paste0( OutDir,"pt_at2.RData" ))
	save(pt_at1,file=paste0( OutDir,"pt_at1.RData" ))
}

epithelial_cells_dotplot = function( OutDir ){
	load(file = paste0(OutDir,"../extendedAtlasNormal_seu_CancerEpithelialS0.RData"))
	seu = subset(seu,cells = seu@meta.data[seu@meta.data$Epi_Cancer!="Cancer","CellID"])
	Idents(seu) = seu@meta.data$Epi_Cancer
	am = FindAllMarkers(seu,logfc.threshold=0.5,only.pos=T)
	am = am[am$p_val_adj<0.01,]
	am = am[order(-am$avg_log2FC),]
	tam = am[am$avg_log2FC>1,]
	# top3
	tg = c(  )
	for (ct in c( "AT1","AT2","AT0","preTB","club","ciliated","basal" )){
		tam = am[(am$cluster==ct) & ((am$pct.1-am$pct.2)>0.2 ),]
		tam = am[(am$cluster==ct),]
		tam = am[(am$cluster==ct) & ((am$pct.1)>0.5 ),]
		tg = c(tg,tam[1:4,"gene"])
	}
	tg = tg[!is.na(tg)]
	# AT0 from literature: CEACAM6, SCGB3A2, 
	# preTB from literature: 
	Idents(seu) = factor(seu@meta.data$Epi_Cancer,levels = rev(c( "AT1","AT2","AT0","preTB","club","ciliated","basal" ) ))
	library(viridis)
	pdf( paste0(OutDir,"dotplot_top4.pdf" ),9,4 )
	plot = DotPlot(seu, features = unique(tg), ) + RotatedAxis() + scale_colour_viridis(option="rocket",direction=-1)
	print(plot)
	dev.off()
}

differential_expression_NormalCells = function( OutDir, withPreinvasive = FALSE ){
	library(limma)
	library(edgeR)
	library(gplots)
	library(msigdbr)
	library(fgsea)
	if (!withPreinvasive){
		load( file = paste0("data/extendedAtlasNormal_clin.RData"))
		load( file = paste0(OutDir,"../../extendedAtlasNormal_annot.RData") )
		load( file = paste0(OutDir,"../../extendedAtlasNormal_counts.RData") ) # ambient_genes already removed	
	} else {
		load(paste0("data/extendedAtlasS0_clin.RData"))
		load(file = paste0(OutDir,"../../extendedAtlasS0normal_annot.RData"))
		load(file = paste0(OutDir,"../../extendedAtlasS0normal_counts.RData"))
	}
	
	gene_universe = rownames(counts)
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	msig_df_C6 = msigdbr::msigdbr(species = "Homo sapiens", category = "C6")
	msig_df_C8 = msigdbr::msigdbr(species = "Homo sapiens", category = "C8")
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP,msig_df_C6,msig_df_C8 ))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	for (thisct in unique(annot$Epi_Cancer[annot$Epi_Cancer!="Cancer"]) ){
		dcat(thisct)
		Ncutoff = 5
		aall = annot[(annot$Epi_Cancer==thisct) & (annot$SampleType!="Metastasis"),]
		aall$PST = paste0(aall$Patient,substr(aall$SampleType,1,3))
		coall = counts[,rownames(aall)]
		annot_pb = dan.df(unique(aall$PST),c( "Dataset","PST","Patient","TN","SampleType","Stage_collapsed","n_thisct" ))
		counts_pb = dan.df(rownames(coall),unique(aall$PST))
		for (s in unique(aall$PST)){
			# dcat(s)
			this_annot = aall[aall$PST==s,c( "Dataset","PST","Patient","TN","SampleType","Stage_collapsed" )]
			this_annot$n_thisct = nrow(this_annot)
			if (nrow(this_annot)>=Ncutoff){
				annot_pb[s,] = this_annot[1,]
				counts_pb[,s] = rowSums(coall[,rownames(this_annot)])
			}
		}
		annot_pb = annot_pb[!is.na(annot_pb$Dataset),]
		counts_pb = counts_pb[,rownames(annot_pb)]
		annot_pb$cMinimal = "Normal"
		annot_pb[annot_pb$SampleType %in% c( "Primary"),"cMinimal"] = "Tumor"
		dcat(dim(annot_pb))
		dcat(dim(counts_pb))
		dtable(annot_pb$cMinimal)
		save(annot_pb, file = paste0( OutDir,"annot_pb_",thisct,".RData" ))
		save(counts_pb, file = paste0( OutDir,"counts_pb_",thisct,".RData" ))
		condition_version = "cMinimal"
		ta = annot_pb
		ta = ta[order(ta[,condition_version],ta$Dataset),]
		ta$Patient = gsub("-","",ta$Patient)
		samples = rownames(ta)
		ge = counts_pb[,samples]
		groups = factor(ta[samples,condition_version])
		datasets = factor(ta[samples,"Dataset"])
		patients = factor(ta[samples,"Patient"])
		y = DGEList(counts = ge, group = groups, remove.zeros = TRUE , genes = rownames(ge))
		first = TRUE
		for (gg in unique(ta[,condition_version])){
			cols_y = rownames(ta[ta[,condition_version]==gg,] )
			isexpr = rowSums(cpm(y[,cols_y])>1) >= 3
			if (first){
				isexpr_all = isexpr
				first = FALSE
			} else {
				isexpr_all = isexpr_all | isexpr
			}
		}
		y = y[isexpr_all,]
		y = calcNormFactors(y)
		logcpm = edgeR::cpm(y, log=TRUE)
		design = model.matrix(~0+groups)
		colnames(design) = gsub("groups","",colnames(design))# c("acinar","lepidic","papillary","solid", "p2", "p3","p4", "p5", "p6", "p7", "p8", "p9", "p10")
		mc = voom(y, design, plot=F)
		contr.matrix <- makeContrasts(
		   NvsT = Tumor-Normal,
		   levels = colnames(design))
		fit = lmFit(mc,design)
		fit_pair = contrasts.fit(fit, contrasts = contr.matrix)
		fit_pair = eBayes(fit_pair)
		top = topTable(fit_pair, adjust="BH", num=Inf) # coef=2, , coef = c(1,2,3,4,5,6)
		topDecaying = top[(top$logFC<(-1) ) & (top$adj.P.Val<0.01),]
		topGrowing = top[(top$logFC>1) & (top$adj.P.Val<0.01),]
		this_OutDir = paste0(OutDir,thisct,"_pairing_none/")
		dir.create(this_OutDir)
		save(top,file=paste0(this_OutDir,condition_version,"_","topTable.RData"))
		dan.write(data.frame(gene=topGrowing$genes),file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing.txt"))
		dan.write(data.frame(gene=topDecaying$genes),file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying.txt"))
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = topGrowing$genes,
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing_MSigDB.txt"))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,2)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing_MSigDB_simplified.txt"))
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = topDecaying$genes,
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying_MSigDB.txt"))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,2)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying_MSigDB_simplified.txt"))
		## with patient pairing
		design = model.matrix(~0+groups+patients)
		colnames(design) = gsub("groups","",colnames(design))# c("acinar","lepidic","papillary","solid", "p2", "p3","p4", "p5", "p6", "p7", "p8", "p9", "p10")
		mc = voom(y, design, plot=F)
		contr.matrix <- makeContrasts(
		   NvsT = Tumor-Normal,
		   levels = colnames(design))
		fit = lmFit(mc,design)
		fit_pair = contrasts.fit(fit, contrasts = contr.matrix)
		fit_pair = eBayes(fit_pair)
		top = topTable(fit_pair, adjust="BH", num=Inf) # coef=2, , coef = c(1,2,3,4,5,6)
		topDecaying = top[(top$logFC<(-1) ) & (top$adj.P.Val<0.01),]
		topGrowing = top[(top$logFC>1) & (top$adj.P.Val<0.01),]
		this_OutDir = paste0(OutDir,thisct,"_pairing_patient/")
		dir.create(this_OutDir)
		save(top,file=paste0(this_OutDir,condition_version,"_","topTable.RData"))
		dan.write(data.frame(gene=topGrowing$genes),file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing.txt"))
		dan.write(data.frame(gene=topDecaying$genes),file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying.txt"))
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = topGrowing$genes,
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing_MSigDB.txt"))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing_MSigDB_simplified.txt"))
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = topDecaying$genes,
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying_MSigDB.txt"))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying_MSigDB_simplified.txt"))
	}
	for (thisct in c( "AT1","AT2","AT0","preTB","ciliated","club","basal" ) ){
		aall = annot[(annot$Epi_Cancer==thisct) & (annot$SampleType!="Metastasis"),]
		aall$PST = paste0(aall$Patient,substr(aall$SampleType,1,3))
		coall = counts[,rownames(aall)]
		## in single cells
		tseu = CreateSeuratObject(counts = coall, project = "this", meta.data = aall)
		tseu = NormalizeData(tseu)
		tseu = FindVariableFeatures(tseu)
		tseu = ScaleData(tseu)
		tseu = RunPCA(tseu)
		Idents(tseu) = tseu@meta.data$TN
		mm = FindMarkers(tseu,ident.1='Tumor',ident.2='Normal',logfc.threshold=0,test.use='MAST',latent.vars='Patient')
		this_OutDir = paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/")
		dir.create(this_OutDir)
		save(mm,file=paste0(this_OutDir,"mm.RData"))
		topGrowing = mm[(mm$avg_log2FC>(0.25) ) & (mm$p_val_adj<0.01),]
		topDecaying = mm[(mm$avg_log2FC<(-0.25)) & (mm$p_val_adj<0.01),]
		dan.write(data.frame(gene=rownames(topGrowing)),file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing.txt"))
		dan.write(data.frame(gene=rownames(topDecaying)),file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying.txt"))
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = rownames(topGrowing),
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing_MSigDB.txt"))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,2)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topGrowing_MSigDB_simplified.txt"))
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = rownames(topDecaying),
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying_MSigDB.txt"))
		a = as.data.frame(fgRes)
		a$overlapGenes = NULL
		a$pval = NULL
		a$padj = signif(a$padj,2)
		dan.write(a,file=paste0(this_OutDir,condition_version,"_","NvsT_topDecaying_MSigDB_simplified.txt"))
	}
	### intersecting pseudobulk and single-cell results
	# common genes
	for (thisct in unique(annot$Epi_Cancer[annot$Epi_Cancer!="Cancer"]) ){
		dcat(thisct)
		dcat("vs pairing none",1)
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_none/cMinimal_NvsT_topGrowing.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecells_pairing_patient/cMinimal_NvsT_topGrowing.txt"))
		iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		uz = union(b_g_genes$gene,sc_g_genes$gene)
		dcat(paste0( "growing sc (",length(sc_g_genes$gene),") vs bpn (",length(b_g_genes$gene),"): ",length(iz)," / ",length(uz) ))
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_none/cMinimal_NvsT_topDecaying.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecells_pairing_patient/cMinimal_NvsT_topDecaying.txt"))
		iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		uz = union(b_g_genes$gene,sc_g_genes$gene)
		dcat(paste0( "decaying sc (",length(sc_g_genes$gene),") vs bpn (",length(b_g_genes$gene),"): ",length(iz)," / ",length(uz) ))
		dcat("vs pairing dataset",1)
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_dataset/cMinimal_NvsT_topGrowing.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecells_pairing_patient/cMinimal_NvsT_topGrowing.txt"))
		iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		uz = union(b_g_genes$gene,sc_g_genes$gene)
		dcat(paste0( "growing sc (",length(sc_g_genes$gene),") vs bpd (",length(b_g_genes$gene),"): ",length(iz)," / ",length(uz) ))
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_dataset/cMinimal_NvsT_topDecaying.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecells_pairing_patient/cMinimal_NvsT_topDecaying.txt"))
		iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		uz = union(b_g_genes$gene,sc_g_genes$gene)
		dcat(paste0( "decaying sc (",length(sc_g_genes$gene),") vs bpd (",length(b_g_genes$gene),"): ",length(iz)," / ",length(uz) ))
		dcat("vs pairing patient",1)
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topGrowing.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecells_pairing_patient/cMinimal_NvsT_topGrowing.txt"))
		iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		uz = union(b_g_genes$gene,sc_g_genes$gene)
		dcat(paste0( "growing sc (",length(sc_g_genes$gene),") vs bpp (",length(b_g_genes$gene),"): ",length(iz)," / ",length(uz) ))
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topDecaying.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecells_pairing_patient/cMinimal_NvsT_topDecaying.txt"))
		iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		uz = union(b_g_genes$gene,sc_g_genes$gene)
		dcat(paste0( "decaying sc (",length(sc_g_genes$gene),") vs bpp (",length(b_g_genes$gene),"): ",length(iz)," / ",length(uz) ))
	}
	# common pathways
	for (thisct in c( "AT0","AT1","AT2" ) ){
		dcat(thisct)
		dcat("vs pairing patient",1)
		dcat("growing",2)
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
		# print(sc_g_genes$pathway[substr(sc_g_genes$pathway,1,4) %in% c( "TRAV","GOBP","HALL" )])
		# print(b_g_genes$pathway[substr(b_g_genes$pathway,1,4) %in% c( "TRAV","GOBP","HALL" )])
		iz = intersect(b_g_genes$pathway,sc_g_genes$pathway)
		iz = iz[substr(iz,1,4) %in% c( "TRAV","GOBP","HALL" )]
		b_g_genes = b_g_genes[b_g_genes$pathway %in% iz,]
		sc_g_genes = sc_g_genes[sc_g_genes$pathway %in% iz,]
		rownames(sc_g_genes) = sc_g_genes$pathway
		rownames(b_g_genes) = b_g_genes$pathway
		b_g_genes$sum_pvals = sc_g_genes[rownames(b_g_genes),"padj"]+b_g_genes[rownames(b_g_genes),"padj"]
		b_g_genes = b_g_genes[order(b_g_genes$sum_pvals),]
		print(b_g_genes$pathway)
		dcat("decaying",2)
		b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
		sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
		# print(sc_g_genes$pathway[substr(sc_g_genes$pathway,1,4) %in% c( "TRAV","GOBP","HALL" )])
		# print(b_g_genes$pathway[substr(b_g_genes$pathway,1,4) %in% c( "TRAV","GOBP","HALL" )])
		iz = intersect(b_g_genes$pathway,sc_g_genes$pathway)
		iz = iz[substr(iz,1,4) %in% c( "TRAV","GOBP","HALL" )]
		b_g_genes = b_g_genes[b_g_genes$pathway %in% iz,]
		sc_g_genes = sc_g_genes[sc_g_genes$pathway %in% iz,]
		rownames(sc_g_genes) = sc_g_genes$pathway
		rownames(b_g_genes) = b_g_genes$pathway
		b_g_genes$sum_pvals = sc_g_genes[rownames(b_g_genes),"padj"]+b_g_genes[rownames(b_g_genes),"padj"]
		b_g_genes = b_g_genes[order(b_g_genes$sum_pvals),]
		print(b_g_genes$pathway)
		# # saving common genes and pathways
		# b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topGrowing.txt"))
		# sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topGrowing.txt"))
		# iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		# dan.write(data.frame(gene=iz),file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_growing.txt"))
		# fgRes = fgsea::fora(pathways = msig_list,
        #                genes = iz,
        #                universe = gene_universe)
		# fgRes = fgRes[fgRes$padj <= 0.01,]
		# a = as.data.frame(fgRes)
		# a = apply(a,2,as.character)
		# dan.write(a,file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_growing_MSigDB.txt") )
		# a = as.data.frame(fgRes)
		# a$overlapGenes = NULL
		# a$pval = NULL
		# a$padj = signif(a$padj,2)
		# dan.write(a,file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_growing_MSigDB_simplified.txt") )
		# b_g_genes = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topDecaying.txt"))
		# sc_g_genes = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topDecaying.txt"))
		# iz = intersect(b_g_genes$gene,sc_g_genes$gene)
		# dan.write(data.frame(gene=iz),file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_decaying.txt"))
		# fgRes = fgsea::fora(pathways = msig_list,
        #                genes = iz,
        #                universe = gene_universe)
		# fgRes = fgRes[fgRes$padj <= 0.01,]
		# a = as.data.frame(fgRes)
		# a = apply(a,2,as.character)
		# dan.write(a,file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_decaying_MSigDB.txt") )
		# a = as.data.frame(fgRes)
		# a$overlapGenes = NULL
		# a$pval = NULL
		# a$padj = signif(a$padj,2)
		# dan.write(a,file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_decaying_MSigDB_simplified.txt") )
		a = dan.read(file=paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/common_bulkpp_scpp_growing_MSigDB_simplified.txt") )
		print(sc_g_genes$pathway[substr(sc_g_genes$pathway,1,4) %in% c( "TRAV","GOBP","HALL" )])
	}
	# https://bmccancer.biomedcentral.com/articles/10.1186/s12885-023-10523-z
	library(ggbreak)
	## MAST volcanos + selected enrichments
	thisct = "AT1"
	load(file=paste0(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/"),"mm.RData"))
	fileName = paste0(OutDir,"mast_",thisct,"_volcano.pdf")
	x = mm$avg_log2FC
	y = -log10(mm$p_val)
	mm$fill = ifelse( (mm$p_val_adj>=0.01) | (abs(mm$avg_log2FC)<=0.25),"Not significant",ifelse( mm$avg_log2FC<0,"Higher in normal","Higher in tumor"))
	mm$fill = factor(mm$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
	fillColorz = c( "gray44","firebrick","gray88" )
	mm$repel_labelz = rownames(mm)
	mm[!(rownames(mm) %in% c( "CAV1","AGER","NAPSA","KRT17","SFTPB" )),"repel_labelz"] = ""
	dan.scatterplot( fileName, x, y, fill = mm$fill, xlab = "log2(fold-change)", ylab = "-log10(MAST patient-corrected p-value)", filllab = "", fillColors = fillColorz, dotSize = 0.1, repel_labels=mm$repel_labelz, coord_fixed = FALSE, fileWidth = 2.5, fileHeight = 2.25, legend_position = "bottom" )
	# top enrichments
	enr = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(61,93,96,103,121)),]
	selected_enrich = enr[,c( "pathway","padj" )]	
	selected_enrich$minusLog10qval = -log10(selected_enrich$padj)
	selected_enrich$alias = c( "Travaglini AT2 cells","Hallmark interferon gamma response","GO:BP Antigen processing and presentation","GO:BP Immune response","Hallmark interferon alpha response" )
	selected_enrich$side = "Higher in tumor"
	enr = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(1,19,23,24,28)),]
	enr$minusLog10qval = -log10(enr$padj)
	enr$alias = c( "Travaglini AT1 cells","GO:BP Biological adhesion","GO:BP Cell morphogenesis","GO:BP Cell migration","GO:BP Regulation of cell adhesion")
	enr$side = "Higher in normal"
	selected_enrich = rbind(selected_enrich,enr[,colnames(selected_enrich)])
	selected_enrich = selected_enrich[order(selected_enrich$minusLog10qval),]
	ts = selected_enrich[selected_enrich$side=="Higher in normal",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"mast_",thisct,"_SelectedEnrichments_Normal.pdf"),2.6,nrow(ts)*0.2 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="gray44") + theme_classic() + ylab("-log10(p-value)" ) + scale_x_discrete(position = "top") + coord_flip() + scale_y_reverse() + xlab("")
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	# plot = plot + scale_y_break(c(8, 60),scales = 1)
	print(plot)
	dev.off()
	ts = selected_enrich[selected_enrich$side=="Higher in tumor",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"mast_",thisct,"_SelectedEnrichments_Tumor.pdf"),2.6,nrow(ts)*0.2 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="firebrick") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()


	thisct = "AT0"
	load(file=paste0(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/"),"mm.RData"))
	fileName = paste0(OutDir,"mast_",thisct,"_volcano.pdf")
	x = mm$avg_log2FC
	y = -log10(mm$p_val)
	mm$fill = ifelse( (mm$p_val_adj>=0.01) | (abs(mm$avg_log2FC)<=0.25),"Not significant",ifelse( mm$avg_log2FC<0,"Higher in normal","Higher in tumor"))
	mm$fill = factor(mm$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
	fillColorz = c( "gray44","firebrick","gray88" )
	mm$repel_labelz = rownames(mm)
	mm[!(rownames(mm) %in% c( "SLPI","CAV1","NAPSA","IFI27","SCGB3A2" )),"repel_labelz"] = ""
	dan.scatterplot( fileName, x, y, fill = mm$fill, xlab = "log2(fold-change)", ylab = "-log10(MAST patient-corrected p-value)", filllab = "", fillColors = fillColorz, dotSize = 0.1, repel_labels=mm$repel_labelz, coord_fixed = FALSE, fileWidth = 2.5, fileHeight = 2.25, legend_position = "bottom" )
	# top enrichments
	enr = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(37,44,61,131,137)),]
	selected_enrich = enr[,c( "pathway","padj" )]	
	selected_enrich$minusLog10qval = -log10(selected_enrich$padj)
	selected_enrich$alias = c( "GO:BP Peptide metabolic process","Travaglini AT2 cells","GO:BP Oxidative phosphorylation","GO:BP Antigen processing and presentation","Hallmark interferon gamma response" )
	selected_enrich$side = "Higher in tumor"
	enr = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(21,28,32,46,66)),]
	enr$minusLog10qval = -log10(enr$padj)
	enr$alias = c( "GO:BP Cell morphogenesis","Travaglini AT1 cells","GO:BP Cell migration","GO:BP Biological adhesion","GO:BP Tissue development")
	enr$side = "Higher in normal"
	selected_enrich = rbind(selected_enrich,enr[,colnames(selected_enrich)])
	selected_enrich = selected_enrich[order(selected_enrich$minusLog10qval),]
	ts = selected_enrich[selected_enrich$side=="Higher in normal",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"mast_",thisct,"_SelectedEnrichments_Normal.pdf"),2.6,nrow(ts)*0.2 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="gray44") + theme_classic() + ylab("-log10(p-value)" ) + scale_x_discrete(position = "top") + coord_flip() + scale_y_reverse() + xlab("")
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	# plot = plot + scale_y_break(c(8, 60),scales = 1)
	print(plot)
	dev.off()
	ts = selected_enrich[selected_enrich$side=="Higher in tumor",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"mast_",thisct,"_SelectedEnrichments_Tumor.pdf"),2.6,nrow(ts)*0.2 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="firebrick") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()

	thisct = "AT2"
	load(file=paste0(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/"),"mm.RData"))
	fileName = paste0(OutDir,"mast_",thisct,"_volcano.pdf")
	x = mm$avg_log2FC
	y = -log10(mm$p_val)
	mm$fill = ifelse( (mm$p_val_adj>=0.01) | (abs(mm$avg_log2FC)<=0.25),"Not significant",ifelse( mm$avg_log2FC<0,"Higher in normal","Higher in tumor"))
	mm$fill = factor(mm$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
	fillColorz = c( "gray44","firebrick","gray88" )
	mm$repel_labelz = rownames(mm)
	mm[!(rownames(mm) %in% c( "WFDC2","IFI27","CXCL14","IFI6","MSMO1","CTSE","FABP4" )),"repel_labelz"] = ""
	dan.scatterplot( fileName, x, y, fill = mm$fill, xlab = "log2(fold-change)", ylab = "-log10(MAST patient-corrected p-value)", filllab = "", fillColors = fillColorz, dotSize = 0.1, repel_labels=mm$repel_labelz, coord_fixed = FALSE, fileWidth = 2.5, fileHeight = 2.25, legend_position = "bottom" )
	# top enrichments
	enr = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(22,25,27,61,83)),]
	selected_enrich = enr[,c( "pathway","padj" )]	
	selected_enrich$minusLog10qval = -log10(selected_enrich$padj)
	selected_enrich$alias = c( "Hallmark interferon alpha response","Hallmark interferon gamma response","GO:BP Defense response","GO:BP Peptide metabolic process","GO:BP Antigen processing and presentation" )
	selected_enrich$side = "Higher in tumor"
	enr = dan.read(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(1,2,3,7,23)),]
	enr$minusLog10qval = -log10(enr$padj)
	enr$alias = c( "Travaglini AT2 cells","Hallmark cholesterol homeostasis","GO:BP Sterol metabolic process","GO:BP Lipid biosynthetic process","GO:BP Tissue development")
	enr$side = "Higher in normal"
	selected_enrich = rbind(selected_enrich,enr[,colnames(selected_enrich)])
	selected_enrich = selected_enrich[order(selected_enrich$minusLog10qval),]
	ts = selected_enrich[selected_enrich$side=="Higher in normal",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"mast_",thisct,"_SelectedEnrichments_Normal.pdf"),2.6,nrow(ts)*0.2 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="gray44") + theme_classic() + ylab("-log10(p-value)" ) + scale_x_discrete(position = "top") + coord_flip() + scale_y_reverse() + xlab("")
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	# plot = plot + scale_y_break(c(8, 60),scales = 1)
	print(plot)
	dev.off()
	ts = selected_enrich[selected_enrich$side=="Higher in tumor",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"mast_",thisct,"_SelectedEnrichments_Tumor.pdf"),2.6,nrow(ts)*0.2 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="firebrick") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()



	# all_enrich = dan.read(file=paste0(paste0(OutDir,thisct,"_singlecellsAllGenes_pairing_patient/"),"cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	# all_enrich = all_enrich[substr(all_enrich$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]

	# all_enrich$minusLog10qval = -log10(all_enrich$padj)
	# all_enrich = all_enrich[order(all_enrich$tps,all_enrich$padj),]
	# all_enrich$xlabz = substr(as.character(all_enrich$term),1,20)
	# all_enrich$term = paste0("tps",tps_map[all_enrich$tps,"rankz"],"_",as.character(all_enrich$term))
	# all_enrich$term = factor(all_enrich$term,levels=as.character(all_enrich$term) )
	# all_enrich$tps = factor(all_enrich$tps,levels=tps_map$tps[tps_map$tps %in% all_enrich$tps ] )
	# pdf( paste0( OutDir,"barplots_selected_enrichments.pdf" ),16,5 )
	# plot = ggplot(all_enrich, aes(x=term,y=minusLog10qval, fill=tps)) + geom_bar(stat='identity') + facet_wrap(~tps,scale="free_x",nrow=1) + scale_fill_manual(values=tps_map[levels(all_enrich$tps),"colorz"]) + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_x_discrete(breaks=all_enrich$term,labels=all_enrich$xlabz) + ylab("-log10(q-value)" )
	# print(plot)
	# dev.off()

	## pseudobulk paired patients volcanos
	thisct = "AT1"
	load(file=paste0(paste0(OutDir,thisct,"_pairing_patient/"),"cMinimal_topTable.RData"))
	fileName = paste0(OutDir,"PseudoBulkpp_",thisct,"_volcano.pdf")
	x = top$logFC
	y = -log10(top$P.Value)
	top$fill = ifelse( (top$adj.P.Val>=0.01) | (abs(top$logFC)<=1),"Not significant",ifelse( top$logFC<0,"Higher in normal","Higher in tumor"))
	top$fill = factor(top$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
	fillColorz = c( "gray44","firebrick","gray88" )
	top$repel_labelz = rownames(top)
	top[!(rownames(top) %in% c( "CAV1","AGER","CXCL14","MDK","NAPSA","ATOH8","LAMA3","UNC13D","RGCC","SCTR","ETV1" )),"repel_labelz"] = ""
	dan.scatterplot( fileName, x, y, fill = top$fill, xlab = "log2(fold-change)", ylab = "-log10(patient-corrected p-value)", filllab = "", fillColors = fillColorz, dotSize = 1, repel_labels=top$repel_labelz, coord_fixed = FALSE, fileWidth = 7, fileHeight = 7 )
	# top enrichments
	enr = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	selected_enrich = enr[,c( "pathway","pval" )]
	selected_enrich$minusLog10qval = -log10(selected_enrich$pval)
	selected_enrich$alias = c( "Travaglini AT2 cells","GO:BP Tube development","GO:BP Tube morphogenesis","Travaglini lung mucous cells" )
	selected_enrich$side = "Higher in tumor"
	enr = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(1,2,3,7,10,19)),]
	enr$minusLog10qval = -log10(enr$pval)
	enr$alias = c( "Travaglini AT1 cells","GO:BP Biological adhesion","GO:BP Regulation of cell adhesion","GO:BP Tissue development","GO:BP Tube morphogenesis","GO:BP Tube development" )
	enr$side = "Higher in normal"
	selected_enrich = rbind(selected_enrich,enr[,colnames(selected_enrich)])
	selected_enrich = selected_enrich[order(selected_enrich$minusLog10qval),]
	ts = selected_enrich[selected_enrich$side=="Higher in normal",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"PseudoBulkpp_",thisct,"_SelectedEnrichments_Normal.pdf"),5,nrow(ts)*0.5 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="gray44") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	print(plot)
	dev.off()
	ts = selected_enrich[selected_enrich$side=="Higher in tumor",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"PseudoBulkpp_",thisct,"_SelectedEnrichments_Tumor.pdf"),5,nrow(ts)*0.5 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="firebrick") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	print(plot)
	dev.off()

	thisct = "AT0"
	load(file=paste0(paste0(OutDir,thisct,"_pairing_patient/"),"cMinimal_topTable.RData"))
	fileName = paste0(OutDir,"PseudoBulkpp_",thisct,"_volcano.pdf")
	x = top$logFC
	y = -log10(top$P.Value)
	top$fill = ifelse( (top$adj.P.Val>=0.01) | (abs(top$logFC)<=1),"Not significant",ifelse( top$logFC<0,"Higher in normal","Higher in tumor"))
	top$fill = factor(top$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
	fillColorz = c( "gray44","firebrick","gray88" )
	top$repel_labelz = rownames(top)
	top[!(rownames(top) %in% c( "FABP4","EGR2","SLPI","IFI27","CAV1","COL1A1","SERINC2","CXCL14" )),"repel_labelz"] = ""
	dan.scatterplot( fileName, x, y, fill = top$fill, xlab = "log2(fold-change)", ylab = "-log10(patient-corrected p-value)", filllab = "", fillColors = fillColorz, dotSize = 1, repel_labels=top$repel_labelz, coord_fixed = FALSE, coord_flipped=FALSE, fileWidth = 7, fileHeight = 7 )
	# top enrichments
	enr = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	selected_enrich = enr[,c( "pathway","pval" )]
	selected_enrich$minusLog10qval = -log10(selected_enrich$pval)
	selected_enrich$alias = c( "Hallmark EMT","Travaglini lung goblet cells","Travaglini lung serous cells" )
	selected_enrich$side = "Higher in tumor"
	enr = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(14,34,35,44,46,56,64)),]
	enr$minusLog10qval = -log10(enr$pval)
	enr$alias = c( "Travaglini ciliated cells","GO:BP Defense response","GO:BP Immune response","Travaglini AT1 cells","GO:BP Inflammatory response","GO:BP Biological adhesion","GO:BP Response to cytokines" )
	enr$side = "Higher in normal"
	selected_enrich = rbind(selected_enrich,enr[,colnames(selected_enrich)])
	selected_enrich = selected_enrich[order(selected_enrich$minusLog10qval),]
	ts = selected_enrich[selected_enrich$side=="Higher in normal",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"PseudoBulkpp_",thisct,"_SelectedEnrichments_Normal.pdf"),5,nrow(ts)*0.5 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="gray44") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	print(plot)
	dev.off()
	ts = selected_enrich[selected_enrich$side=="Higher in tumor",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"PseudoBulkpp_",thisct,"_SelectedEnrichments_Tumor.pdf"),5,nrow(ts)*0.5 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="firebrick") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	print(plot)
	dev.off()

	thisct = "AT2"
	load(file=paste0(paste0(OutDir,thisct,"_pairing_patient/"),"cMinimal_topTable.RData"))
	fileName = paste0(OutDir,"PseudoBulkpp_",thisct,"_volcano.pdf")
	x = top$logFC
	y = -log10(top$P.Value)
	top$fill = ifelse( (top$adj.P.Val>=0.01) | (abs(top$logFC)<=1),"Not significant",ifelse( top$logFC<0,"Higher in normal","Higher in tumor"))
	top$fill = factor(top$fill,levels=c( "Higher in normal","Higher in tumor","Not significant" ))
	fillColorz = c( "gray44","firebrick","gray88" )
	top$repel_labelz = rownames(top)
	top[!(rownames(top) %in% c( "HMGCS1","CCL13","CXCL14","IFI44L","WFDC2","FABP4","CSF3","TTN","AGER" )),"repel_labelz"] = ""
	dan.scatterplot( fileName, x, y, fill = top$fill, xlab = "log2(fold-change)", ylab = "-log10(patient-corrected p-value)", filllab = "", fillColors = fillColorz, dotSize = 1, repel_labels=top$repel_labelz, coord_fixed = FALSE, coord_flipped=FALSE, fileWidth = 7, fileHeight = 6 )
	# top enrichments
	enr = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topGrowing_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	selected_enrich = enr[as.character(c(11,32,68,75,98,102,117)),c( "pathway","pval" )]
	selected_enrich$minusLog10qval = -log10(selected_enrich$pval)
	selected_enrich$alias = c( "GO:BP Inflammatory response","GO:BP Biological adhesion","GO:BP Immune response","GO:BP Response to chemokines","GO:BP Defense response","GO:BP Cell-cell adhesion","GO:BP Cell migration" )
	selected_enrich$side = "Higher in tumor"
	enr = dan.read(paste0(OutDir,thisct,"_pairing_patient/cMinimal_NvsT_topDecaying_MSigDB_simplified.txt"))
	enr = enr[substr(enr$pathway,1,4) %in% c( "TRAV","GOBP","HALL" ),]
	enr = enr[as.character(c(1,2,6,17,20)),]
	enr$minusLog10qval = -log10(enr$pval)
	enr$alias = c( "Travaglini AT2 cells","GO:BP Sterol metabolism","Hallmark cholesterol homeostasis","GO:BP Lipid metabolism","Hallmark MTORC1 signaling" )
	enr$side = "Higher in normal"
	selected_enrich = rbind(selected_enrich,enr[,colnames(selected_enrich)])
	selected_enrich = selected_enrich[order(selected_enrich$minusLog10qval),]
	ts = selected_enrich[selected_enrich$side=="Higher in normal",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"PseudoBulkpp_",thisct,"_SelectedEnrichments_Normal.pdf"),5,nrow(ts)*0.5 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="gray44") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	print(plot)
	dev.off()
	ts = selected_enrich[selected_enrich$side=="Higher in tumor",]
	ts$alias = factor(ts$alias,levels = as.character(ts$alias))
	pdf( paste0(OutDir,"PseudoBulkpp_",thisct,"_SelectedEnrichments_Tumor.pdf"),5,nrow(ts)*0.5 )
	plot = ggplot(ts, aes(x=alias,y=minusLog10qval)) + geom_bar(stat='identity',fill="firebrick") + theme_classic() + ylab("-log10(p-value)" ) + coord_flip() + xlab("")
	print(plot)
	dev.off()
}

extract_tumorspecific_genes = function( OutDir ){
	library(limma)
	library(edgeR)
	Ncutoff = 5
	load(file = paste0(OutDir,"../../tps_discovery/tps.RData"))
	load(file = paste0(OutDir,"../../tps_discovery/tps_universe.RData"))
	CellTypeDir = paste0("data/cell_typing/non_malignant/")
	load(file = paste0(OutDir,"../../extendedAtlas_annot.RData"))
	annot = annot[annot$SampleType!="Metastasis",]
	load(file = paste0(OutDir,"../../extendedAtlas_counts.RData"))
	annot$PCT = paste0(annot$Patient,"_",annot$Epi_Cancer)
	apb = dan.df(unique(annot$PCT),c( "PCT","CellType","Dataset","Patient","SampleType","ncells" ))
	cpb = dan.df(rownames(counts),unique(annot$PCT))
	dcat( "Cancer cells" )
	for (pct in unique(annot$PCT)){
		# dcat(pct,1)
		ta = annot[annot$PCT==pct,]
		if (nrow(ta)<Ncutoff) { next }
		apb[pct,"ncells"] = nrow(ta)
		apb[pct,c( "PCT","CellType","Dataset","Patient","SampleType" )] = as.character(ta[1,c( "PCT","Epi_Cancer","Dataset","Patient","SampleType" )])
		cpb[,pct] = rowSums(counts[,ta$CellID])
	}
	for (dataset in unique(annot$Dataset)){
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_seu.RData"))
		md2 = seu@meta.data
		md2 = md2[((md2$SampleType != "Metastasis") & (!(md2$annd_level_1 %in% c("glial","unclear")))) & ( md2$annd_level_2!="unclear" ),]
		matt2 = seu@assays$RNA@counts
		matt2 = matt2[intersect(rownames(matt2),tps_universe),rownames(md2)]
		md2$PCT = paste0(md2$Patient,"_",md2$annd_level_2)
		this_apb = dan.df(unique(md2$PCT),c( "PCT","CellType","Dataset","Patient","SampleType","ncells" ))
		this_cpb = dan.df(rownames(matt2),unique(md2$PCT))
		dcat( "Non-malignant cells",1 )
		for (pct in unique(md2$PCT)){
			ta = md2[md2$PCT==pct,]
			if (nrow(ta)<Ncutoff) { next }
			if (nrow(ta))
			this_apb[pct,"ncells"] = nrow(ta)
			this_apb[pct,c( "PCT","CellType","Dataset","Patient","SampleType" )] = as.character(ta[1,c( "PCT","annd_level_2","Dataset","Patient","SampleType" )])
			this_cpb[,pct] = rowSums(matt2[,ta$CellID])
		}
		apb = rbind(apb,this_apb[,colnames(apb)])
		common_genez = intersect(rownames(cpb),rownames(this_cpb))
		cpb = cbind(cpb[common_genez,],this_cpb[common_genez,])
	}
	apb = apb[!is.na(apb$ncells),]
	cpb = cpb[,rownames(apb)]
	apb$is_cancer = "Cancer"
	apb[apb$CellType!="Cancer","is_cancer"] = "Non_malignant"
	save(apb,file = paste0(OutDir,"apb.RData"))
	save(cpb,file = paste0(OutDir,"cpb.RData"))

	### Cancer vs non-cancer
	ta = apb
	ta$Patient = gsub("-","",ta$Patient)
	samples = rownames(ta)
	ge = cpb[,samples]
	groups = factor(ta[samples,"is_cancer"])
	datasets = factor(ta[samples,"Dataset"])
	patients = factor(ta[samples,"Patient"])
	y = DGEList(counts = ge, group = groups, remove.zeros = TRUE , genes = rownames(ge))
	first = TRUE
	for (gg in unique(ta[,"is_cancer"])){
		cols_y = rownames(ta[ta[,"is_cancer"]==gg,] )
		isexpr = rowSums(cpm(y[,cols_y])>1) >= 3
		if (first){
			isexpr_all = isexpr
			first = FALSE
		} else {
			isexpr_all = isexpr_all | isexpr
		}
	}
	y = y[isexpr_all,]
	y = calcNormFactors(y)
	logcpm = edgeR::cpm(y, log=TRUE)
	design = model.matrix(~0+groups)
	colnames(design) = gsub("groups","",colnames(design))# c("acinar","lepidic","papillary","solid", "p2", "p3","p4", "p5", "p6", "p7", "p8", "p9", "p10")
	mc = voom(y, design, plot=F)
	contr.matrix <- makeContrasts(
	   CvsN = Cancer-Non_malignant,
	   levels = colnames(design))
	fit = lmFit(mc,design)
	fit_pair = contrasts.fit(fit, contrasts = contr.matrix)
	fit_pair = eBayes(fit_pair)
	top = topTable(fit_pair, adjust="BH", num=Inf) # coef=2, , coef = c(1,2,3,4,5,6)
	save(top,file=paste0(OutDir,"is_cancer","_","topTable.RData"))
	cancer_markers = top[(top$logFC>1) & (top$adj.P.Val<0.01),"genes"]
	save(cancer_markers,file=paste0(OutDir,"is_cancer_markers_restrictive.RData"))
	cancer_markers = top[(top$logFC>0) & (top$adj.P.Val<0.01),"genes"]
	save(cancer_markers,file=paste0(OutDir,"is_cancer_markers_permissive.RData"))
	# now, correct for dataset
	design = model.matrix(~0+groups+datasets)
	colnames(design) = gsub("groups","",colnames(design))# c("acinar","lepidic","papillary","solid", "p2", "p3","p4", "p5", "p6", "p7", "p8", "p9", "p10")
	contr.matrix <- makeContrasts(
	   CvsN = Cancer-Non_malignant,
	   levels = colnames(design))
	fit = lmFit(mc,design)
	fit_pair = contrasts.fit(fit, contrasts = contr.matrix)
	fit_pair = eBayes(fit_pair)
	top = topTable(fit_pair, adjust="BH", num=Inf) # coef=2, , coef = c(1,2,3,4,5,6)
	save(top,file=paste0(OutDir,"is_cancer","_","topTable_DatasetCorrected.RData"))
	cancer_markers = top[(top$logFC>1) & (top$adj.P.Val<0.01),"genes"]
	save(cancer_markers,file=paste0(OutDir,"is_cancer_markers_restrictive_DatasetCorrected.RData"))
	cancer_markers = top[(top$logFC>0) & (top$adj.P.Val<0.01),"genes"]
	save(cancer_markers,file=paste0(OutDir,"is_cancer_markers_permissive_DatasetCorrected.RData"))
}

tumorspecific_enrichments = function( OutDir,universe="permissive",signature="TwoStates",scoring_method='singscore' ){

	library(singscore)
	OutDir = paste0( OutDir,"universe_",universe,"_signature_",signature,"/" )
	dir.create(OutDir)
	recalib = function(values,thMax,thMin,recalib_percent){
	   unit = recalib_percent/100
	   newv = ((values-(min(c(thMin,thMax))) )/abs(thMax-thMin))*(1/unit)
	   return(newv)
	}

	# defining universe and signatures
	if (universe=="restrictive"){ load( file=paste0(OutDir,"../is_cancer_markers_restrictive.RData") ) }
	if (universe=="permissive"){ load( file=paste0(OutDir,"../is_cancer_markers_permissive.RData") ) }
	gene_universe = cancer_markers
	if (signature=="restrictive"){ load( file=paste0(OutDir,"../../../CellStates/signatures/wilcox_cs_signatures_restrictive.RData" ) ) }
	if (signature=="permissive"){ load( file=paste0(OutDir,"../../../CellStates/signatures/wilcox_cs_signatures_permissive.RData" ) ) }
	if (signature=="TwoStates"){ load( file=paste0(OutDir,"../../../CellStates/signatures/wilcox_csLevel1_signatures_permissive.RData" ) ) }
	cs_signatures = sapply(cs_signatures, function(x) intersect(gene_universe,x))
	print(sapply(cs_signatures,length))
	cox_MultiVar_df = dan.df( c("TCGA","ChenEAS"),c( "dataset","Npatients","dataset_label",paste0("HR_",names(cs_signatures)), paste0("HR_lower95_",names(cs_signatures)), paste0("HR_upper95_",names(cs_signatures)),paste0("pval_",names(cs_signatures)) ) )
	cox_MultiVar_df_stageI = dan.df( c( "TCGA","ChenEAS" ),c( "dataset","Npatients","dataset_label",paste0("HR_",names(cs_signatures)), paste0("HR_lower95_",names(cs_signatures)), paste0("HR_upper95_",names(cs_signatures)),paste0("pval_",names(cs_signatures)) ) )
	
	load( paste0(OutDir,"../../../tps_discovery/tps.RData") )
	load( paste0(OutDir,"../../../tps_discovery/tps_universe.RData") )
	tps_universe = intersect(tps_universe,gene_universe)

	### restrictive, restrictive
	# Alveolar Proliferative       Hypoxic 
    #       172            47            31 
    ### restrictive, permissive
     # Alveolar Proliferative       Hypoxic 
     #      202            77            68
	### permissive, restrictive
	# Alveolar Proliferative       Hypoxic 
    #       230           133            50 
    ### permissive, permissive
    # Alveolar Proliferative       Hypoxic 
    #       277           293           139

    # two states, permissive
    # cs1 cs2 
	# 277 462 

	# two states, restrictive
    # cs1 cs2 
	# 202 161

	### cell states across bulk sample types
	load(file=paste0(OutDir,"../../../tps_across_ClinicalCharacteristics/Bulk_Tps_scores_Metastases.RData"))
	load(file=paste0(OutDir,"../../../tps_across_ClinicalCharacteristics/Bulk_logtpm_Metastases.RData"))
	logtpm_all = logtpm_all[,clinall$Sample]
	rankData = rankGenes(logtpm_all)
	for (tp in names(cs_signatures))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   clinall[rownames(scoredf),tp] = scoredf$TotalScore
	}
	tc = clinall
	tc$SampleType = paste0(tc$SampleType, "_",tc$BiopsySite)
	tc[tc$SampleType=="Normal_lung","SampleType"] = "Normal lung"
	tc[tc$SampleType=="Primary_lung","SampleType"] = "Primary tumor"
	tc[tc$SampleType=="Metastasis_LN","SampleType"] = "Mets, LN"
	tc[tc$SampleType=="Metastasis_lung","SampleType"] = "Mets, lung"
	tc[tc$SampleType=="Metastasis_other","SampleType"] = "Mets, distal"

	xLevels = c( "Normal lung" ,"Primary tumor","Mets, LN","Mets, lung","Mets, distal" )
	xColors = c( "gray","gold","tomato3","firebrick","sienna4" )

	# xLevels = c("Normal_lung","Primary_lung","Metastasis_lung","Metastasis_LN","Metastasis_other")
	# xColors = c("gray","gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$SampleType,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_SampleType_AllDatasets.pdf"),3,2 )
	for (tp in names(cs_signatures)){
	   x = tc[,"SampleType"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "SampleType", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()

	## patterns in lumu
	load( file = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/ExprTables/lung_multiregion_log2TPMcombat_table_GeneSymbols.RData" )
	logtpm = ge[intersect(rownames(ge),gene_universe),]
	Clin = read.table(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Clinical_4assigned.txt"), sep = "\t", header = TRUE, quote = '', stringsAsFactors = FALSE)
	lepidic_color = "dodgerblue4"
	acinar_color = "orange"
	papillary_color = "lightseagreen"
	solid_color = "red"
	Clin$color = lepidic_color # c("#FF8000", "#6666FF", "#009900")
	Clin[Clin$Pattern=="acinar","color"] = acinar_color
	Clin[Clin$Pattern=="papillary","color"] = papillary_color
	Clin[Clin$Pattern=="solid","color"] = solid_color
	Clin = Clin[Clin$Region!="A",] # assigns it to parent scope
	rownames(Clin) = Clin$Sample
	rankData = rankGenes(logtpm)
	pdf( paste0(OutDir,"Bulk_cs_vs_HistologicPatterns_CHUV.pdf"),7,5 )
	for (cs in names(cs_signatures)){
	   x = Clin[,"Pattern"]
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[cs]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[cs]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   Clin[rownames(scoredf),cs] = scoredf$TotalScore
	   y = Clin[,cs]
	   plot=dan.boxplots.multipages( factor(x,levels=c("lepidic","papillary","acinar","solid")), y, xlab = "", ylab = paste0(cs," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = c(lepidic_color,papillary_color,acinar_color,solid_color), jitterColors = Clin$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	for (cs in names(tps)){
	   scoredf = simpleScore(rankData, upSet = tps[[cs]])
	   Clin[rownames(scoredf),cs] = scoredf$TotalScore
	}
	save(Clin,file=paste0(OutDir,"Bulk_cs_scores_CHUV.RData"))
	cordf = cor(Clin[,names(tps)])
	cordf = cordf[rowSums(is.na(cordf))!=(nrow(cordf)-1),colSums(is.na(cordf))!=(nrow(cordf)-1)]
	library(corrplot)
	pdf(paste0(OutDir,"csMECO_CHUV_pearson.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()


	## patterns, stages, has_metastasis, age, smoking, survival in TCGA
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_TPM/TCGA_LUAD_ge_TPM_GeneSymbols.RData")
	cn = as.character(colnames(ge))
	cn = cn[substr(cn,14,15) %in% c("01")] # only primary
	ge = ge[intersect(rownames(ge),gene_universe),cn]
	cn_short = substr(cn,1,15)
	colnames(ge) = cn_short
	logtpm = log2(ge+1)
	logtpm_tcga = logtpm
	CommonDataDir = "/mnt/ed2/daniele/Common_Data/"
	Clinical_extended = as.data.frame(t(read.table(paste0(CommonDataDir,"Clinical_extended/LUAD.Clinical.txt"), quote = '', sep = "\t", header = FALSE, row.names = 1)), stringsAsFactors = FALSE)
	Clinical_extended[,"Patient"] <- toupper(as.character(Clinical_extended[,"bcr_patient_barcode"]))
	rownames(Clinical_extended) = Clinical_extended[,"Patient"]
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/TCGA_histology_formatted_permissive_withNormals.RData"))
	rownames(Clin) = Clin$Sample
	Clin_patterns = Clin
	Clin_patterns$Pattern = as.character(Clin_patterns$Pattern)
	Clin = dan.df(colnames(ge),c("Sample","Patient","Pattern","color","pathologic_stage","pathologic_m","pathologic_n","age_at_initial_pathologic_diagnosis","gender","tobacco_smoking_history","vital_status", "days_to_death", "days_to_last_followup"))
	Clin$Sample = rownames(Clin)
	Clin$Patient = substr(Clin$Sample,1,12)
	Clin[intersect(rownames(Clin),rownames(Clin_patterns)),"Pattern"] = Clin_patterns[intersect(rownames(Clin),rownames(Clin_patterns)),"Pattern"]
	Clin[intersect(rownames(Clin),rownames(Clin_patterns)),"color"] = Clin_patterns[intersect(rownames(Clin),rownames(Clin_patterns)),"color"]
	for (pat in rownames(Clinical_extended)){ 
		for (var in c("Patient","pathologic_stage","pathologic_m","pathologic_n","age_at_initial_pathologic_diagnosis","gender","tobacco_smoking_history","vital_status", "days_to_death", "days_to_last_followup"))
		Clin[Clin$Patient==pat,var] = Clinical_extended[pat,var]
	}
	Clin$stage_collapsed = NA
	Clin[Clin$pathologic_stage %in% c("stage i","stage ia", "stage ib"),"stage_collapsed"] = "I"
	Clin[Clin$pathologic_stage %in% c("stage iia", "stage iib"),"stage_collapsed"] = "II"
	Clin[Clin$pathologic_stage %in% c("stage iia", "stage iib"),"stage_collapsed"] = "II"
	Clin[Clin$pathologic_stage %in% c("stage iiia", "stage iiib"),"stage_collapsed"] = "III"
	Clin[Clin$pathologic_stage %in% c("stage iv"),"stage_collapsed"] = "IV"
	Clin = Clin[Clin$Sample %in% colnames(logtpm),]
	Clin$Sample = as.character(Clin$Sample)
	rownames(Clin) = Clin$Sample
	rankData = rankGenes(logtpm)
	for (cs in names(cs_signatures))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[cs]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[cs]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   Clin[rownames(scoredf),cs] = scoredf$TotalScore
	}
	for (cs in names(tps)){
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[cs]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[cs]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   Clin[rownames(scoredf),cs] = scoredf$TotalScore
	}
	save(Clin,file=paste0(OutDir,"Bulk_cs_scores_TCGA.RData"))
	cordf = cor(Clin[,names(tps)])
	cordf = cordf[rowSums(is.na(cordf))!=(nrow(cordf)-1),colSums(is.na(cordf))!=(nrow(cordf)-1)]
	library(corrplot)
	pdf(paste0(OutDir,"csMECO_TCGA_pearson.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()
	save(Clin,file=paste0(OutDir,"Bulk_cs_scores_TCGA.RData"))
	load(file=paste0(OutDir,"Bulk_cs_scores_TCGA.RData"))
	# vs patterns
	tc = Clin[!is.na(Clin$Pattern),]
	tc$color
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_TCGA.pdf"),7,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=c("lepidic","papillary","acinar","solid")), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = max(y), xColors = c(lepidic_color,papillary_color,acinar_color,solid_color), jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs stages
	tc = Clin
	tc[(tc$Pattern=="normal") %in% c(T),"stage_collapsed"] = "Normal"
	tc = tc[!is.na(tc$stage_collapsed),]
	xLevels = c("I","II","III","IV")
	xColors = c("gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$stage_collapsed,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Stages_TCGA.pdf"),7,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"stage_collapsed"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "Stage", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs has_metastasis
	tc = Clin[((!is.na(Clin$pathologic_m)) & (!is.na(Clin$pathologic_n))) & (!((Clin$pathologic_m=="mx") & (Clin$pathologic_n=="nx"))),]
	tc = tc[substr(tc$Sample,14,15)=="01",]
	tc$has_metastasis = "Metastatic_LNonly"
	tc[(tc$pathologic_m %in% c("m0","mx")) & (tc$pathologic_n %in% c("n0","nx")),"has_metastasis"] = "Localized"
	tc[(tc$pathologic_m %in% c("m1","m1a","m1b")),"has_metastasis"] = "Metastatic_distal"
	xLevels = c("Localized","Metastatic_LNonly","Metastatic_distal")
	xColors = c("gold","tomato3","brown")
	tc$color = dan.expand_colors(tc$has_metastasis,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_HasMetastasis_TCGA.pdf"),5,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"has_metastasis"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal.test", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking
	tc = Clin[!is.na(Clin$tobacco_smoking_history),]
	tc = tc[substr(tc$Sample,14,15)=="01",]
	tc$smoking = ifelse(tc$tobacco_smoking_history==1,"Never smoker","Ever smoker")
	xLevels = c("Never smoker","Ever smoker")
	xColors = c("forestgreen","gray22")
	tc$color = dan.expand_colors(tc$smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_TCGA.pdf"),5,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking detailed
	tc = Clin[!is.na(Clin$tobacco_smoking_history),]
	tc = tc[substr(tc$Sample,14,15)=="01",]
	tc = tc[tc$tobacco_smoking_history!=5,]
	tc$smoking = "Never"
	tc[tc$tobacco_smoking_history==2,"smoking"] = "Current"
	tc[tc$tobacco_smoking_history==3,"smoking"] = "Ex (>15y)"
	tc[tc$tobacco_smoking_history==4,"smoking"] = "Ex (<=15y)"
	xLevels = c("Never","Ex (>15y)", "Ex (<=15y)", "Current")
	xColors = c("forestgreen","gray77","gray50","gray22")
	tc$color = dan.expand_colors(tc$smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_SmokingDetailed_TCGA.pdf"),2.1,1.7 )
	for (tp in names(cs_signatures)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs survival
	library(survival)
	this_OutDir = paste0(OutDir,"survival_TCGA/")
	dir.create(this_OutDir)
	colorz = c("dodgerblue4","firebrick")
	legendz_half = c("Below mean","Above mean")
	legendz_extremes = c("Bottom 25%","Top 25%")
	this_Clin = Clin[substr(Clin$Sample,14,15)=="01",]
	rownames(this_Clin) = this_Clin$Patient
	a = as.numeric(as.character(this_Clin$days_to_death))
	b = as.numeric(as.character(this_Clin$days_to_last_followup))
	vital_status_num <- vector(mode="numeric", length=length(this_Clin$vital_status))
	times <- vector(mode="numeric", length=length(vital_status_num))
	for (v in 1:length(vital_status_num))
	{
	  if (this_Clin$vital_status[v]=="alive")
	  {
	     vital_status_num[v] <- 0
	     times[v] <- b[v]
	  }
	  else
	  {
	     vital_status_num[v] <- 1
	     times[v] <- a[v]
	  }
	}
	this_Clin$Times <- times
	this_Clin$vital_status_num <- vital_status_num
	this_Clin$stage_numeric = NA
	this_Clin[(this_Clin$stage_collapsed == "I") %in% c(T),"stage_numeric"] = 1
	this_Clin[(this_Clin$stage_collapsed == "II") %in% c(T),"stage_numeric"] = 2
	this_Clin[(this_Clin$stage_collapsed == "III") %in% c(T),"stage_numeric"] = 3
	this_Clin[(this_Clin$stage_collapsed == "IV") %in% c(T),"stage_numeric"] = 4
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_SplittingHalf.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxUniVariate.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   dataset="TCGA"
	   cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
		

	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxMultiVariate.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_Extremes.pdf")
	   pdf(fileName,4.5,4,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}
	this_Clin_save = this_Clin
	## vs survival, stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed == "I") %in% c(T),]
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_SplittingHalf_stageI.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxUniVariate_stageI.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   dataset="TCGA"
	   cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
	   cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxMultiVariate_stageI.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_Extremes_stageI.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}
	this_Clin = this_Clin_save[(this_Clin_save$stage_collapsed %in% c("I","II") ) %in% c(T),]
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_SplittingHalf_stage12.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxUniVariate_stage12.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxMultiVariate_stage12.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_Extremes_stage12.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}
	this_Clin = this_Clin_save[(this_Clin_save$stage_collapsed %in% c("III","IV") ) %in% c(T),]
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > median(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$gender, age = as.numeric(this_Clin$age_at_initial_pathologic_diagnosis), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_SplittingHalf_stage34.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxUniVariate_stage34.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_CoxMultiVariate_stage34.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "TCGA_survivalBy_",feature,"_Extremes_stage34.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}



	## patterns, stages, age, smoking, survival in Chen
	load( file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_normalized_GeneSymbols.RData") )
	chen_indir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_onChen/"
	load(file = paste0(chen_indir,"ClinChen.RData"))
	logtpm = ge[intersect(rownames(ge),gene_universe),intersect(rownames(ClinChen),colnames(ge))]
	logtpm = logtpm[rowSums(is.na(logtpm))==0, ]
	logtpm_chen = logtpm
	Clin = ClinChen[intersect(rownames(ClinChen),colnames(ge)),]
	rankData = rankGenes(logtpm)
	for (tp in names(cs_signatures))
	{
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   Clin[rownames(scoredf),tp] = scoredf$TotalScore
	}
	for (cs in names(tps)){
	   if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = tps[[cs]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=tps[[cs]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	   Clin[rownames(scoredf),cs] = scoredf$TotalScore
	}
	save(Clin,file=paste0(OutDir,"Bulk_Tps_scores_ChenEAS.RData"))
	# vs patterns
	tc = Clin[!is.na(Clin$Pattern),]
	tc = tc[tc$Pattern!="unknown",]
	pdf( paste0(OutDir,"Bulk_Tps_vs_HistologicPatterns_ChenEAS.pdf"),7,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"Pattern"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=c("lepidic","papillary","acinar","solid","micropapillary")), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = c(lepidic_color,papillary_color,acinar_color,solid_color,"purple"), jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs stages
	tc = Clin	
	tc = tc[!is.na(tc$Stage),]
	xLevels = c("I","II","III","IV")
	xColors = c("gold","orange","tomato3","firebrick")
	tc$color = dan.expand_colors(tc$Stage,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Stages_ChenEAS.pdf"),7,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"Stage"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "Stage", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, comparisons = NULL, labelycoo = 0, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs smoking
	tc = Clin[!is.na(Clin$Smoker),]
	tc$smoking = ifelse(tc$Smoker=="No","Never smoker","Ever smoker")
	xLevels = c("Never smoker","Ever smoker")
	xColors = c("forestgreen","gray22")
	tc$color = dan.expand_colors(tc$smoking,xLevels,xColors)
	pdf( paste0(OutDir,"Bulk_Tps_vs_Smoking_ChenEAS.pdf"),5,5 )
	for (tp in names(cs_signatures)){
	   x = tc[,"smoking"]
	   y = tc[,tp]
	   plot=dan.boxplots.multipages( factor(x,levels=xLevels), y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = "wilcox.test", comparisons = NULL, labelycoo = 0.51, xColors = xColors, jitterColors = tc$color, labelJitteredPoints = NULL, includeJitters = T )
	   print(plot)
	}
	dev.off()
	# vs survival
	this_OutDir = paste0(OutDir,"survival_ChenEAS/")
	dir.create(this_OutDir)
	colorz = c("dodgerblue4","firebrick")
	legendz_half = c("Below mean","Above mean")
	legendz_extremes = c("Bottom 25%","Top 25%")
	this_Clin = Clin
	this_Clin$Times = as.numeric(this_Clin$OS.Month)
	this_Clin$vital_status_num = NA
	this_Clin[this_Clin$OS.Status=="Dead","vital_status_num"] = 1
	this_Clin[this_Clin$OS.Status=="Alive","vital_status_num"] = 0
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$Gender, age = as.numeric(this_Clin$Age), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_SplittingHalf.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxUniVariate.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   dataset="ChenEAS"
	   cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	   cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxMultiVariate.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_Extremes.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}
	this_Clin_save = this_Clin
	## vs survival, stage I only
	this_Clin = this_Clin[(this_Clin$Stage == "I") %in% c(T),]
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$Gender, age = as.numeric(this_Clin$Age), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_SplittingHalf_StageI.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxUniVariate_StageI.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   dataset="ChenEAS"
	   cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
	   cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["Signature_value","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["Signature_value","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["Signature_value","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["Signature_value","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxMultiVariate_StageI.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_Extremes_StageI.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}
	this_Clin = this_Clin_save[(this_Clin_save$Stage %in% c("I","II") ) %in% c(T),]
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$Gender, age = as.numeric(this_Clin$Age), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_SplittingHalf_Stage12.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxUniVariate_Stage12.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxMultiVariate_Stage12.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_Extremes_Stage12.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}
	this_Clin = this_Clin_save[(this_Clin_save$Stage %in% c("III","IV") ) %in% c(T),]
	for (tp in names(cs_signatures))
	{
	   feature = tp
	   this_Clin$enrich = (this_Clin[,feature] > mean(this_Clin[,feature]))
	   Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Signature_value = this_Clin[,feature],
	      Signature = factor(this_Clin$enrich, levels = c(F,T) ), sex = this_Clin$Gender, age = as.numeric(this_Clin$Age), Stage = as.numeric(this_Clin$stage_numeric) )
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_SplittingHalf_Stage34.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs,mark.time=T, col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_half, lty=c(1,1), col=colorz)
	   dev.off()
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxUniVariate_Stage34.txt"))
	   coxr = coxph(as.formula(paste0("SurvObj ~ Signature_value + sex + age + Stage")), data = Surv_df) #  Age + Sex +
	   a = (summary(coxr))
	   capture.output(a, file = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_CoxMultiVariate_Stage34.txt"))
	   bottomThresh = quantile(this_Clin[,feature])["25%"]
	   topThresh = quantile(this_Clin[,feature])["75%"]
	   extreme_Clin = this_Clin[(this_Clin[,feature] >= topThresh) | (this_Clin[,feature] <= bottomThresh),]
	   extreme_Clin$enrich = (extreme_Clin[,feature] >= topThresh)
	   Surv_df = data.frame(Patient = extreme_Clin$Patient, vital_status = as.numeric(extreme_Clin$vital_status_num), Times = extreme_Clin$Times, Signature_value = extreme_Clin[,feature],
	      Signature = factor(extreme_Clin$enrich, levels = c(F,T) ))
	   Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	   km_gs = survfit(SurvObj~Signature, data = Surv_df)
	   km_gs_dif = survdiff(SurvObj~Signature, data = Surv_df, rho = 0)
	   p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	   fileName = paste0(this_OutDir, "ChenEAS_survivalBy_",feature,"_Extremes_Stage34.pdf")
	   pdf(fileName,6,6,useDingbats=F)
	   plot(km_gs, mark.time=T,col=colorz, main = paste0("Survival by ",feature,", - p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	   legend(x = "topright", legend = legendz_extremes, lty=c(1,1), col=colorz)
	   dev.off()
	}

	DataDir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternMicroEnv/ConsensusTME/All_datasets/"
	recalib_percent = 10

	# Beg
	dataset = "Beg"
	PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Beg2016/"
	load(file = paste0(PcDir,"Beg_LUAD_ge_microarray_GeneSymbols.RData"))
	load(file = paste0(PcDir, "Clin_Beg.RData"))
	logtpm = log2(ge[intersect(rownames(ge),gene_universe ),])
	load(paste0(DataDir,"ctme_",dataset,".RData"))
	rankData = rankGenes(logtpm)
	scores = data.frame( row.names = colnames(ctme), Patient = colnames(ctme), Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
	for (tp in names(cs_signatures))
	{
	  if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	  scores[,tp] = scoredf[rownames(scores),"TotalScore"]
	}
	this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender","Smoking","KRAS_status","EGFR_status","STK11_status","TP53_status")], by = "Patient")
	rownames(this_Clin) = this_Clin$Patient
	cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	# Harmonization
	this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
	this_Clin$gender = this_Clin$Gender
	this_Clin$Times = as.numeric(this_Clin$times)
	this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
	### Cox regressions
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}
	## Stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed==1) %in% c(T),]
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
		cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}

	### Micke
	dataset = "Micke"
	PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Micke2018/"
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Micke2018/", "Clin_Micke.RData")) # Clin
	load(file = paste0(PcDir,"Micke_LUAD_ge_FPKM_GeneSymbols.RData"))
	logtpm = log2(ge[intersect(rownames(ge),gene_universe ),]+1)
	load(paste0(DataDir,"ctme_",dataset,".RData"))
	rankData = rankGenes(logtpm)
	scores = data.frame( row.names = colnames(ctme), Patient = colnames(ctme), Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
	for (tp in names(cs_signatures))
	{
	  if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	  scores[,tp] = scoredf[rownames(scores),"TotalScore"]
	}
	this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender")], by = "Patient")
	rownames(this_Clin) = this_Clin$Patient
	cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	# Harmonization
	this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
	this_Clin$gender = this_Clin$Gender
	this_Clin$Times = as.numeric(this_Clin$times)
	this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
	### Cox regressions
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}
	## Stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed==1) %in% c(T),]
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
		cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}

	### Pintilie
	dataset = "Pintilie"
	PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Pintilie2013/"
	load(file = paste0(PcDir,"Pintilie_LUAD_ge_logcounts_GeneSymbols.RData"))
	load(file = paste0(PcDir, "Clin_Pintilie.RData"))
	logtpm = ge[intersect(rownames(ge),gene_universe ),]
	load(paste0(DataDir,"ctme_",dataset,".RData"))
	rankData = rankGenes(logtpm)
	scores = data.frame( row.names = colnames(ctme), Patient = colnames(ctme), Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
	for (tp in names(cs_signatures))
	{
	  if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	  scores[,tp] = scoredf[rownames(scores),"TotalScore"]
	}
	this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender","Smoking")], by = "Patient")
	rownames(this_Clin) = this_Clin$Patient
	cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	# Harmonization
	this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
	this_Clin$gender = this_Clin$Gender
	this_Clin$Times = as.numeric(this_Clin$times)
	this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
	### Cox regressions
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}
	## Stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed==1) %in% c(T),]
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
		cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}

	## Shedden
	dataset = "Shedden"
	PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Shedden2008/"
	load(file = paste0(PcDir,"Shedden_LUAD_ge_rawcounts_GeneSymbols.RData"))
	load(file = paste0(PcDir, "Clin_Shedden.RData"))
	logtpm = log2(ge[intersect(rownames(ge),gene_universe ),])
	load(paste0(DataDir,"ctme_",dataset,".RData"))
	rankData = rankGenes(logtpm)
	scores = data.frame( row.names = colnames(ctme), Patient = colnames(ctme), Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
	for (tp in names(cs_signatures))
	{
	  if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	  scores[,tp] = scoredf[rownames(scores),"TotalScore"]
	}
	this_Clin = merge(scores, Clin[,c("Patient","vital_status_num","times","Stage","Age","Gender", "Adjuvant_chemo", "Adjuvant_rt", "ever_smoked", "grade")], by = "Patient")
	rownames(this_Clin) = this_Clin$Patient
	cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	# Harmonization
	this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
	this_Clin$gender = this_Clin$Gender
	this_Clin$Times = as.numeric(this_Clin$times)
	this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
	### Cox regressions
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}
	## Stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed==1) %in% c(T),]
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
		cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}	

	### Yokota
	dataset = "Yokota"
	PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Yokota2012/"
	load(file = paste0(PcDir,"Yokota_LUAD_ge_rawcounts_GeneSymbols.RData"))
	load(file = paste0(PcDir, "Clin_Yokota.RData"))
	logtpm = log2(ge[intersect(rownames(ge),gene_universe ),])
	load(paste0(DataDir,"ctme_",dataset,".RData"))
	rankData = rankGenes(logtpm)
	scores = data.frame( row.names = colnames(ctme), Patient = colnames(ctme), Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
	for (tp in names(cs_signatures))
	{
	  if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	  scores[,tp] = scoredf[rownames(scores),"TotalScore"]
	}
	this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender")], by = "Patient")
	rownames(this_Clin) = this_Clin$Patient
	cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	# Harmonization
	this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
	this_Clin$gender = this_Clin$Gender
	this_Clin$Times = as.numeric(this_Clin$times)
	this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
	### Cox regressions
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}
	## Stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed==1) %in% c(T),]
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
		cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}	

	### Yokota
	dataset = "Yokota"
	PcDir = "/mnt/ndata/daniele/lung_multiregion/Data/Yokota2012/"
	load(file = paste0(PcDir,"Yokota_LUAD_ge_rawcounts_GeneSymbols.RData"))
	load(file = paste0(PcDir, "Clin_Yokota.RData"))
	logtpm = log2(ge[intersect(rownames(ge),gene_universe ),])
	load(paste0(DataDir,"ctme_",dataset,".RData"))
	rankData = rankGenes(logtpm)
	scores = data.frame( row.names = colnames(ctme), Patient = colnames(ctme), Dataset = paste0("Schabath et al. (N=",nrow(scoredf),")"))
	for (tp in names(cs_signatures))
	{
	  if (scoring_method=='singscore'){ scoredf = simpleScore(rankData, upSet = cs_signatures[[tp]]) }
	   if (scoring_method=='GSVA'){ scoredf = data.frame(t(gsva(as.matrix(logtpm), list(TotalScore=cs_signatures[[tp]]), mx.diff=TRUE, verbose=FALSE, parallel.sz=0, kcdf="Gaussian"))) }
	  scores[,tp] = scoredf[rownames(scores),"TotalScore"]
	}
	this_Clin = merge(scores, Clin[,c("Patient","vital_status_num", "times", "Stage", "Age","Gender")], by = "Patient")
	rownames(this_Clin) = this_Clin$Patient
	cox_MultiVar_df[dataset,"Npatients"] = nrow(this_Clin)
	cox_MultiVar_df[dataset,"dataset_label"] = scores[1,"Dataset"]
	for (cs in names(cs_signatures)){
		this_Clin[,cs] = recalib(values = this_Clin[,cs], thMax = 0.5, thMin = -0.5 , recalib_percent = 10)
	}
	# Harmonization
	this_Clin$age_at_initial_pathologic_diagnosis = as.numeric(this_Clin$Age)
	this_Clin$gender = this_Clin$Gender
	this_Clin$Times = as.numeric(this_Clin$times)
	this_Clin$stage_collapsed = as.numeric(this_Clin$Stage)
	### Cox regressions
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}
	## Stage I only
	this_Clin = this_Clin[(this_Clin$stage_collapsed==1) %in% c(T),]
	this_Clin$SurvObj = with(this_Clin, Surv(Times, as.numeric(as.character(vital_status_num))))
	for (feature in names(cs_signatures)){
		this_Clin$feature = this_Clin[,feature]
		coxr = coxph(as.formula(paste0("SurvObj ~ age_at_initial_pathologic_diagnosis+stage_collapsed+gender+feature")), data = this_Clin) #  Age + Sex +
		summary(coxr)
		a = (summary(coxr))
		cox_MultiVar_df_stageI[dataset,"Npatients"] = nrow(this_Clin)
		cox_MultiVar_df_stageI[dataset,paste0("HR_",feature)] = a$coefficients["feature","exp(coef)"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_lower95_",feature)] = a$conf.int["feature","lower .95"]
		cox_MultiVar_df_stageI[dataset,paste0("HR_upper95_",feature)] = a$conf.int["feature","upper .95"]
		cox_MultiVar_df_stageI[dataset,paste0("pval_",feature)] = a$coefficients["feature","Pr(>|z|)"]
		cox_MultiVar_df_stageI[dataset,"dataset_label"] = paste0( dataset," (N=", nrow(this_Clin),")")
	}

	cox_MultiVar_df = cox_MultiVar_df[cox_MultiVar_df$Npatients>=100,]
	library(forestplot)
	if (signature!="TwoStates"){
		## forest plot
		# clip = c(0.2,10)
		save(cox_MultiVar_df, file = paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllDatasets.RData"))
		load( file = paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllDatasets.RData"))
	   tabletext = cbind(
	     c("Dataset", substr(cox_MultiVar_df$dataset_label,1,nchar(cox_MultiVar_df$dataset_label)-8), "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(Alveolar)", signif(cox_MultiVar_df$HR_Alveolar,3), signif(weighted.mean(x = cox_MultiVar_df$HR_Alveolar, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Alveolar)", formatC(signif(cox_MultiVar_df$pval_Alveolar,2)), NA),
	     c("HR\n(Proliferative)", signif(cox_MultiVar_df$HR_Proliferative,3), signif(weighted.mean(x = cox_MultiVar_df$HR_Proliferative, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Proliferative)", formatC(signif(cox_MultiVar_df$pval_Proliferative,2)), NA),
	     c("HR\n(Hypoxic)", signif(cox_MultiVar_df$HR_Hypoxic,3), signif(weighted.mean(x = cox_MultiVar_df$HR_Hypoxic, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Hypoxic)", formatC(signif(cox_MultiVar_df$pval_Hypoxic,2)), NA),
	     rep("      ",length(c("Dataset", rownames(cox_MultiVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllDatasets.pdf"),6.5,3, useDingbats=F,onefile=F,pointsize=6)
	   fp=forestplot(tabletext, 
	            legend = c("Alveolar","Proliferative","Hypoxic"),
	            align = c("l","l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(2,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_Alveolar,NA),c(NA,cox_MultiVar_df$HR_Proliferative,NA),c(NA,cox_MultiVar_df$HR_Hypoxic,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_Alveolar,NA),c(NA,cox_MultiVar_df$HR_lower95_Proliferative,NA),c(NA,cox_MultiVar_df$HR_lower95_Hypoxic,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_Alveolar,NA),c(NA,cox_MultiVar_df$HR_upper95_Proliferative,NA),c(NA,cox_MultiVar_df$HR_upper95_Hypoxic,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("dodgerblue4","firebrick3","sienna4"),line=c("dodgerblue4","firebrick3","sienna4")))
	   print(fp)
	   dev.off()
	   cox_MultiVar_df_store = cox_MultiVar_df
		cox_MultiVar_df = cox_MultiVar_df_stageI	   
		cox_MultiVar_df = cox_MultiVar_df[cox_MultiVar_df$Npatients>=100,]
		save(cox_MultiVar_df, file = paste0(OutDir,"cox_MultiVar_df_stageIonly_ForestPlot_AllDatasets.RData"))
	   ## forest plot
		# clip = c(0.2,10)
	   tabletext = cbind(
	     c("Dataset", substr(cox_MultiVar_df$dataset_label,1,nchar(cox_MultiVar_df$dataset_label)-8), "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(Alveolar)", signif(cox_MultiVar_df$HR_Alveolar,3), signif(weighted.mean(x = cox_MultiVar_df$HR_Alveolar, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Alveolar)", formatC(signif(cox_MultiVar_df$pval_Alveolar,2)), NA),
	     c("HR\n(Proliferative)", signif(cox_MultiVar_df$HR_Proliferative,3), signif(weighted.mean(x = cox_MultiVar_df$HR_Proliferative, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Proliferative)", formatC(signif(cox_MultiVar_df$pval_Proliferative,2)), NA),
	     c("HR\n(Hypoxic)", signif(cox_MultiVar_df$HR_Hypoxic,3), signif(weighted.mean(x = cox_MultiVar_df$HR_Hypoxic, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Hypoxic)", formatC(signif(cox_MultiVar_df$pval_Hypoxic,2)), NA),
	     rep("      ",length(c("Dataset", rownames(cox_MultiVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_MultiVar_df_stageIonly_ForestPlot_AllDatasets.pdf"),6.5,3, useDingbats=F,onefile=F,pointsize=6)
	   fp=forestplot(tabletext, 
	            legend = c("Alveolar","Proliferative","Hypoxic"),
	            align = c("l","l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(2,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_Alveolar,NA),c(NA,cox_MultiVar_df$HR_Proliferative,NA),c(NA,cox_MultiVar_df$HR_Hypoxic,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_Alveolar,NA),c(NA,cox_MultiVar_df$HR_lower95_Proliferative,NA),c(NA,cox_MultiVar_df$HR_lower95_Hypoxic,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_Alveolar,NA),c(NA,cox_MultiVar_df$HR_upper95_Proliferative,NA),c(NA,cox_MultiVar_df$HR_upper95_Hypoxic,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("dodgerblue4","firebrick3","sienna4"),line=c("dodgerblue4","firebrick3","sienna4")))
	   print(fp)
	   dev.off()
	} else {
		## forest plot
		# clip = c(0.2,10)
		save(cox_MultiVar_df, file = paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllDatasets.RData"))
		load( file = paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllDatasets.RData"))
		cox_MultiVar_df["Micke","dataset_label"] = gsub("Micke","Mezheyeuski et al.",cox_MultiVar_df["Micke","dataset_label"])

	   tabletext = cbind(
	     c("Dataset", substr(cox_MultiVar_df$dataset_label,1,nchar(cox_MultiVar_df$dataset_label)-8), "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(Alveolar)", signif(cox_MultiVar_df$HR_cs1,3), signif(weighted.mean(x = cox_MultiVar_df$HR_cs1, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Alveolar)", formatC(signif(cox_MultiVar_df$pval_cs1,2)), NA),
	     c("HR\n(Dediff.)", signif(cox_MultiVar_df$HR_cs2,3), signif(weighted.mean(x = cox_MultiVar_df$HR_cs2, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Dediff.)", formatC(signif(cox_MultiVar_df$pval_cs2,2)), NA),
	     rep("      ",length(c("Dataset", rownames(cox_MultiVar_df), "Summary"))))

	   pdf(paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllDatasets.pdf"),5,3, useDingbats=F,onefile=F,pointsize=6)
	   fp=forestplot(tabletext, 
	            legend = c("Alveolar","Dedifferentiated"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(2,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_cs1,NA),c(NA,cox_MultiVar_df$HR_cs2,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_cs1,NA),c(NA,cox_MultiVar_df$HR_lower95_cs2,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_cs1,NA),c(NA,cox_MultiVar_df$HR_upper95_cs2,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("#0F4F8B","#C2B59B"),line=c("#0F4F8B","#C2B59B")))
	   print(fp)
	   dev.off()
	   cox_MultiVar_df_store = cox_MultiVar_df
		cox_MultiVar_df = cox_MultiVar_df_stageI	   
		cox_MultiVar_df = cox_MultiVar_df[cox_MultiVar_df$Npatients>=100,]
		save(cox_MultiVar_df, file = paste0(OutDir,"cox_MultiVar_df_stageIonly_ForestPlot_AllDatasets.RData"))
	   ## forest plot
		# clip = c(0.2,10)
	   tabletext = cbind(
	     c("Dataset", substr(cox_MultiVar_df$dataset_label,1,nchar(cox_MultiVar_df$dataset_label)-8), "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(Alveolar)", signif(cox_MultiVar_df$HR_cs1,3), signif(weighted.mean(x = cox_MultiVar_df$HR_cs1, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Alveolar)", formatC(signif(cox_MultiVar_df$pval_cs1,2)), NA),
	     c("HR\n(Dedifferentiated)", signif(cox_MultiVar_df$HR_cs2,3), signif(weighted.mean(x = cox_MultiVar_df$HR_cs2, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(Dedifferentiated)", formatC(signif(cox_MultiVar_df$pval_cs2,2)), NA),
	     rep("      ",length(c("Dataset", rownames(cox_MultiVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_MultiVar_df_stageIonly_ForestPlot_AllDatasets.pdf"),5.5,3, useDingbats=F,onefile=F,pointsize=6)
	   fp=forestplot(tabletext, 
	            legend = c("Alveolar","Dedifferentiated"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(2,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_cs1,NA),c(NA,cox_MultiVar_df$HR_cs2,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_cs1,NA),c(NA,cox_MultiVar_df$HR_lower95_cs2,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_cs1,NA),c(NA,cox_MultiVar_df$HR_upper95_cs2,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("dodgerblue4","firebrick4"),line=c("dodgerblue4","firebrick4")))
	   print(fp)
	   dev.off()
	}
}

visium_tps_scoring = function( OutDir,method="AMS",canceronly_genes=FALSE ){
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	this_OutDir = paste0(OutDir,"tps_",method,"_canceronlygenes",as.character(canceronly_genes),"/")
	dir.create(this_OutDir)
	load(paste0(OutDir,"../tps_discovery/tps.RData"))
	cs_signatures = tps
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))
	if (canceronly_genes){ 
		load( file=paste0(OutDir,"../bulk_deconvolutions/tumorspecific_enrichments/is_cancer_markers_permissive.RData") )
		subspa_all = subset(subspa_all,features=intersect(rownames(subspa_all),cancer_markers))
	}
	if (method=="AMS"){
		subspa_all = AddModuleScore(subspa_all,cs_signatures)
		colnames(subspa_all@meta.data)[substr(colnames(subspa_all@meta.data),1,4)=="Clus"] = names(cs_signatures)
		scores = subspa_all@meta.data
	}
	if (method=="AmsCentered"){
		dcat( "AmsCentered" )
		# db_rand = dan.barkley_MakeRand_data(subspa_all,cs_signatures, 3)
		# save(db_rand,file = paste0(this_OutDir,"db_rand_unintegrated.RData"))
		# scores = dan.Barkley_GeneToEnrichment_AmsCentered( subspa_all, cs_signatures, db_rand)
		# save(scores,file = paste0(this_OutDir,"tps_scores.RData"))
		# colz_to_remove = colnames(scores)[substr(colnames(scores),nchar(colnames(scores))-2,nchar(colnames(scores)))==".cv"]
		# colz_to_remove = c(colz_to_remove,gsub(".cv","",colz_to_remove))
		# scores = scores[,!(colnames(scores) %in% colz_to_remove)]
		# save(scores,file = paste0(this_OutDir,"tps_scores.RData"))
		load(file = paste0(this_OutDir,"tps_scores.RData"))
	}
	## CancerOnly
	DecDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Deconvolution/"
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024/",Sample,"_thetadf.RData")) 
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall = thetadf
	    first = FALSE
	  } else {
	    thetall = rbind(thetall,thetadf)
	  }
	}
	colnames(thetall) = paste0( "rctd2_",colnames(thetall) )
	scores = cbind(scores,thetall[rownames(scores),])
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024_fine_reference/",Sample,"_thetadf.RData"))  
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall_fine = thetadf
	    first = FALSE
	  } else {
	    thetall_fine = rbind(thetall_fine,thetadf)
	  }
	}
	colnames(thetall_fine) = paste0( "rctd3_",colnames(thetall_fine) )
	scores = cbind(scores,thetall_fine[rownames(scores),])

	thresh = 0.1
	scc = scores[scores$rctd2_Cancer>thresh,]
	cordf = cor(scc[,names(cs_signatures)],method = "pearson")
	library(corrplot)
	order_tps = c( "AT2-like","AT2-Club-like","Club-like","lepidic-like","MHC-II","Stress_AP1","Stress_HSP","Stress_secreted","Unas_emp","RNA_processing","MHC-I","Interferon","OxPhos","Cell_proliferation","Basal-like","pEMT","EMT","Metal","Hypoxia","Translation_initiation","Unfolded_protein_response" )
	library(colorRamp2)
	library(ComplexHeatmap)
	pdf(paste0(this_OutDir,"csMECO_SingleSpots_CancerRich",thresh,"_pearson.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()
	pdf(paste0(this_OutDir,"csMECO_SingleSpots_CancerRich",thresh,"_pearson_tps_ordersc.pdf"),5,5)
	plot=(ComplexHeatmap::Heatmap(cordf[order_tps,order_tps],rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	library(ppcor)
	cordf = dan.df(names(tps),names(tps),1)
	cordf_pvals = dan.df(names(tps),names(tps),1)
	for (rn in rownames(cordf)){
		for (cn in rownames(cordf)){
			if (rn==cn){ next }
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
			cordf_pvals[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$p.value
		}
	}
	pdf(paste0(this_OutDir,"csMECO_SingleSpots_CancerRich",thresh,"_pearson_tps_ppcor.pdf"),17,17)
	corrplot(as.matrix(cordf),order='AOE',col = colorRampPalette(c("dodgerblue4","white","firebrick3"))(100),tl.col="black",addCoef.col="black",)
	dev.off()
	pdf(paste0(this_OutDir,"csMECO_SingleSpots_CancerRich",thresh,"_pearson_tps_ppcor_ordersc.pdf"),5,5)
	plot=(ComplexHeatmap::Heatmap(cordf[order_tps,order_tps],rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	## Spatial feature plot for each transcriptional program
	library(RColorBrewer)
	this_subspa_all = subset(subspa_all,cells=rownames(scc))
	for (cs in names(cs_signatures)){
	  dcat(cs)
	  this_subspa_all@meta.data[rownames(scc),cs] = scc[,cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,20)
	  print(SpatialFeaturePlot(this_subspa_all, features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 4) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,20)
	  pp = SpatialFeaturePlot(this_subspa_all, features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 4)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off() 
	}
	pdf(file = paste0(this_OutDir,"CancerRich_","Images_SpatialPlot.pdf"),20,20)
	pp = SpatialFeaturePlot(this_subspa_all,features = cs, crop=FALSE, alpha=0) + patchwork::plot_layout(ncol = 4, nrow = 4)
	print(pp)
	dev.off()
	# Only solid samples
	solid_samples = c( "R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scc), subset = Sample %in% solid_samples)
	for (cs in names(cs_signatures)){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs] = scc[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),10,5)
	  print(SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),10,5)
	  pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_","Images_SpatialPlot.pdf"),10,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples),features = cs, crop=FALSE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"CancerRich_SolidSamples_csMECO.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()
	# Only lepidic and solid samples
	lepidic_solid_samples = c("R26097-1I-Z1","R26097-1E-Z2","R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scc), subset = Sample %in% lepidic_solid_samples)
	for (cs in names(cs_signatures)){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs] = scc[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,5)
	  print(SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,5)
	  pp = SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_","Images_SpatialPlot.pdf"),20,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",lepidic_solid_samples),features = cs, crop=FALSE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"CancerRich_LepidicSolidSamples_csMECO.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()

	# all spots
	## Spatial feature plot for each transcriptional program
	library(RColorBrewer)
	this_subspa_all = subset(subspa_all,cells=rownames(scores))
	for (cs in names(cs_signatures)){
	  dcat(cs)
	  this_subspa_all@meta.data[rownames(scores),cs] = scores[,cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"AllSpots_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,20)
	  print(SpatialFeaturePlot(this_subspa_all, features = cs, crop=TRUE) + patchwork::plot_layout(ncol = 4, nrow = 4) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"AllSpots_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,20)
	  pp = SpatialFeaturePlot(this_subspa_all, features = cs, crop=TRUE) + patchwork::plot_layout(ncol = 4, nrow = 4)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off() 
	}
	pdf(file = paste0(this_OutDir,"AllSpots_","Images_SpatialPlot.pdf"),20,20)
	pp = SpatialFeaturePlot(this_subspa_all,features = cs, crop=TRUE, alpha=0) + patchwork::plot_layout(ncol = 4, nrow = 4)
	print(pp)
	dev.off()
	# Only solid samples
	solid_samples = c( "R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scores), subset = Sample %in% solid_samples)
	for (cs in names(cs_signatures)){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs] = scores[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"AllSpots_SolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),10,5)
	  print(SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"AllSpots_SolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),10,5)
	  pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"AllSpots_SolidSamples_","Images_SpatialPlot.pdf"),10,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples),features = cs, crop=TRUE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"AllSpots_SolidSamples_csMECO.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()

	# Only lepidic and solid samples
	lepidic_solid_samples = c("R26097-1F-Z2","R26097-1E-Z2","R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scores), subset = Sample %in% lepidic_solid_samples)
	for (cs in names(cs_signatures)){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs] = scores[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"AllSpots_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,5)
	  print(SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"AllSpots_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,5)
	  pp = SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"AllSpots_LepidicSolidSamples_","Images_SpatialPlot.pdf"),20,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",lepidic_solid_samples),features = cs, crop=TRUE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"AllSpots_LepidicSolidSamples_csMECO.pdf"),7,7)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black")
	dev.off()
}

visium_cellstates_scoring = function( OutDir,method="AMS",canceronly_genes=FALSE,signature="ThreeStates" ){
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	this_OutDir = paste0(OutDir,method,"_canceronlygenes",as.character(canceronly_genes),"_",signature,"/")
	dir.create(this_OutDir)
	if (signature=="TwoStates"){ load( file=paste0(OutDir,"../CellStates/signatures/wilcox_csLevel1_signatures_permissive.RData" ) ) }
	if (signature=="ThreeStates"){ load( file=paste0(OutDir,"../CellStates/signatures/wilcox_cs_signatures_restrictive.RData" ) ) }
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))
	if (canceronly_genes){ 
		load( file=paste0(OutDir,"../bulk_deconvolutions/tumorspecific_enrichments/is_cancer_markers_permissive.RData") )
		subspa_all = subset(subspa_all,features=intersect(rownames(subspa_all),cancer_markers))
	}
	if (method=="AMS"){
		subspa_all = AddModuleScore(subspa_all,cs_signatures)
		colnames(subspa_all@meta.data)[substr(colnames(subspa_all@meta.data),1,4)=="Clus"] = names(cs_signatures)
		scores = subspa_all@meta.data
	}
	if (method=="AmsCentered"){
		dcat( "AmsCentered" )
		# db_rand = dan.barkley_MakeRand_data(subspa_all,cs_signatures, 3)
		# save(db_rand,file = paste0(this_OutDir,"db_rand_unintegrated.RData"))
		# scores = dan.Barkley_GeneToEnrichment_AmsCentered( subspa_all, cs_signatures, db_rand)
		# save(scores,file = paste0(this_OutDir,"tps_scores.RData"))
		# colz_to_remove = colnames(scores)[substr(colnames(scores),nchar(colnames(scores))-2,nchar(colnames(scores)))==".cv"]
		# colz_to_remove = c(colz_to_remove,gsub(".cv","",colz_to_remove))
		# scores = scores[,!(colnames(scores) %in% colz_to_remove)]
		# save(scores,file = paste0(this_OutDir,"tps_scores.RData"))
		load(file = paste0(this_OutDir,"tps_scores.RData"))
	}
	## CancerOnly
	DecDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Deconvolution/"
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024/",Sample,"_thetadf.RData")) 
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall = thetadf
	    first = FALSE
	  } else {
	    thetall = rbind(thetall,thetadf)
	  }
	}
	colnames(thetall) = paste0( "rctd2_",colnames(thetall) )
	save(thetall,file=paste0( MainDir,"visium_downstream/",method,"_canceronlygenes",as.character(canceronly_genes),"_",signature,"/thetall_level2.RData"))
	scores = cbind(scores,thetall[rownames(scores),])
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024_fine_reference/",Sample,"_thetadf.RData"))  
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall_fine = thetadf
	    first = FALSE
	  } else {
	    thetall_fine = rbind(thetall_fine,thetadf)
	  }
	}
	colnames(thetall_fine) = paste0( "rctd3_",colnames(thetall_fine) )
	save(thetall_fine,file=paste0( MainDir,"visium_downstream/",method,"_canceronlygenes",as.character(canceronly_genes),"_",signature,"/thetall_fine_level3.RData"))
	rowSums(thetall_fine)
	scores = cbind(scores,thetall_fine[rownames(scores),])
	scores$CellState = names(cs_signatures)[apply(scores[,names(cs_signatures)],1,which.max)]
	if (signature=="ThreeStates"){ 
		scores$CellState = factor(scores$CellState,levels=c( "Alveolar","Proliferative","Hypoxic" )) 
		state_colorz = c("Alveolar"="dodgerblue4","Proliferative"="firebrick3","Hypoxic"="sienna4")
	}
	if (signature=="TwoStates"){ 
		scores$CellState = factor(scores$CellState,levels=c( "cs1","cs2" )) 
		state_colorz = c("cs1"="dodgerblue4","cs2"="firebrick4")
	}

	thresh = 0.1
	scc = scores[scores$rctd2_Cancer>thresh,]
	cordf = cor(scc[,names(cs_signatures)],method = "pearson")
	library(corrplot)
	pdf(paste0(this_OutDir,"csMECO_SingleSpots_CancerRich",thresh,"_pearson.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()
	library(ppcor)
	cordf = dan.df(names(cs_signatures),names(cs_signatures),1)
	cordf_pvals = dan.df(names(cs_signatures),names(cs_signatures),1)
	for (rn in rownames(cordf)){
		for (cn in rownames(cordf)){
			if (rn==cn){ next }
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
			cordf_pvals[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$p.value
		}
	}
	pdf(paste0(this_OutDir,"csMECO_SingleSpots_CancerRich",thresh,"_pearson_cs_signatures_ppcor.pdf"),5,5)
	corrplot(as.matrix(cordf),order='AOE',col = colorRampPalette(c("dodgerblue4","white","firebrick3"))(100),tl.col="black",addCoef.col="black",)
	dev.off()

	## patient-wise
	all_patients = unique(substr( rownames(scc),1,6 ))
	for (patient in all_patients){
		dcat(patient)
		these_cells = rownames(scc)[substr(rownames(scc),1,6)==patient]
		cordf = cor(scc[these_cells,names(cs_signatures)],method = "pearson")
		pdf(paste0(this_OutDir,"csMECO_",patient,"_SingleSpots_CancerRich",thresh,"_pearson.pdf"),5,5)
		corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
		dev.off()
		library(ppcor)
		cordf = dan.df(names(cs_signatures),names(cs_signatures),1)
		cordf_pvals = dan.df(names(cs_signatures),names(cs_signatures),1)
		for (rn in rownames(cordf)){
			for (cn in rownames(cordf)){
				if (rn==cn){ next }
				cordf[rn,cn] = pcor.test(scc[these_cells,rn],scc[these_cells,cn],scc[these_cells,"rctd2_Cancer"],method='pearson')$estimate
				cordf_pvals[rn,cn] = pcor.test(scc[these_cells,rn],scc[these_cells,cn],scc[these_cells,"rctd2_Cancer"],method='pearson')$p.value
			}
		}
		pdf(paste0(this_OutDir,"csMECO_",patient,"_SingleSpots_CancerRich",thresh,"_pearson_cs_signatures_ppcor.pdf"),5,5)
		corrplot(as.matrix(cordf),order='AOE',col = colorRampPalette(c("dodgerblue4","white","firebrick3"))(100),tl.col="black",addCoef.col="black",)
		dev.off()
	}
	mm = melt(scc[,c( "Patient","Alveolar","Proliferative","Hypoxic" )])
	x = mm$Patient
	y = mm$value
	fill = factor(mm$variable,levels=c("Alveolar","Proliferative","Hypoxic"))
	fillColors = c("dodgerblue4","firebrick3","sienna4")
	dan.boxplots( paste0( this_OutDir,"PatientWise_scores.pdf" ), x, y, fill, xlab = "Patient", ylab = "Score", filllab = "Cell state", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColors, jitterColors = "black", labelJitteredPoints = NULL, jitterDotSize = 0.05, fileWidth = 4, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )

	## Spatial feature plot for each transcriptional program
	library(RColorBrewer)
	this_subspa_all = subset(subspa_all,cells=rownames(scc))
	for (cs in c(names(cs_signatures))){ #,"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[rownames(scc),cs] = scc[,cs]

	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  # pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,20)
	  # print(SpatialFeaturePlot(this_subspa_all, features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 4) )
	  # dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,20)
	  pp = SpatialFeaturePlot(this_subspa_all, features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 4)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()

	  pp = SpatialFeaturePlot(this_subspa_all, features = cs, crop=FALSE,pt.size.factor=1.2,stroke=0) + patchwork::plot_layout(ncol = 1, nrow = 16 )
	  pp = pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score))
	  pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_SameScale_plotOnly.pdf"),5,20)
	  print(pp & theme(legend.position="none") )
	  dev.off()
	}

	non_lepidic_solid_samples = unique(subspa_all$Sample)[!(unique(subspa_all$Sample) %in% c("R26097-1F-Z2","R26097-1E-Z2","R24586-3D-Z1","R24586-3D-Z2" ))]
	this_subspa_all = subset(subspa_all,cells=rownames(scc), subset = Sample %in% non_lepidic_solid_samples)
	for (cs in c(names(cs_signatures))){ #,"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs] = scc[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  # pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,20)
	  # print(SpatialFeaturePlot(this_subspa_all, features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 4) )
	  # dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,16)
	  pp = SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",non_lepidic_solid_samples),features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 3)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()

	  pp = SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",non_lepidic_solid_samples),features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 12, nrow = 1 )
	  pp = pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score))
	  pdf(file = paste0(this_OutDir,"CancerRich_",cs,"_SpatialFeaturePlot_SameScale_plotOnly.pdf"),60,10)
	  print(pp & theme(legend.position="none") )
	  dev.off()

	}

	pdf(file = paste0(this_OutDir,"CancerRich_","Images_SpatialPlot.pdf"),20,20)
	pp = SpatialFeaturePlot(this_subspa_all,features = cs, crop=FALSE, alpha=0) + patchwork::plot_layout(ncol = 4, nrow = 4)
	print(pp)
	dev.off()
	this_subspa_all@meta.data$CellState = scores[rownames(this_subspa_all@meta.data),"CellState"]
	pdf(file = paste0(this_OutDir,"CancerRich_","AssignedCellStates_SpatialPlot.pdf"),20,20)
	pp = SpatialDimPlot(this_subspa_all,group.by = "CellState", cols=state_colorz, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(ncol = 4, nrow = 4)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"CancerRich_csMECO.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()
	# Only solid samples
	solid_samples = c( "R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scc), subset = Sample %in% solid_samples)
	for (cs in c(names(cs_signatures),"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs] = scc[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),10,5)
	  print(SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),10,5)
	  pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_","Images_SpatialPlot.pdf"),10,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples),features = cs, crop=FALSE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	this_subspa_all@meta.data$CellState = scores[rownames(this_subspa_all@meta.data),"CellState"]
	pdf(file = paste0(this_OutDir,"CancerRich_SolidSamples_","AssignedCellStates_SpatialPlot.pdf"),10,5)
	pp = SpatialDimPlot(this_subspa_all,images=gsub("-",".",solid_samples),group.by = "CellState", cols=state_colorz, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"CancerRich_SolidSamples_csMECO.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()
	# Only lepidic and solid samples
	lepidic_solid_samples = c("R26097-1F-Z2","R26097-1E-Z2","R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scc), subset = Sample %in% lepidic_solid_samples)
	for (cs in c(names(cs_signatures),"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs] = scc[intersect(rownames(scc),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,5)
	  print(SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,5)
	  pp = SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_","Images_SpatialPlot.pdf"),20,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",lepidic_solid_samples),features = cs, crop=FALSE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	this_subspa_all@meta.data$CellState = scores[rownames(this_subspa_all@meta.data),"CellState"]
	pdf(file = paste0(this_OutDir,"CancerRich_LepidicSolidSamples_","AssignedCellStates_SpatialPlot.pdf"),20,5)
	pp = SpatialDimPlot(this_subspa_all,images=gsub("-",".",lepidic_solid_samples),group.by = "CellState", cols=state_colorz, crop=FALSE,pt.size.factor=1.2) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"CancerRich_LepidicSolidSamples_csMECO.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()

	# all spots
	## Spatial feature plot for each transcriptional program
	library(RColorBrewer)
	this_subspa_all = subset(subspa_all,cells=rownames(scores))
	for (cs in c(names(cs_signatures),"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[rownames(scores),cs] = scores[,cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"AllSpots_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,20)
	  print(SpatialFeaturePlot(this_subspa_all, features = cs, crop=TRUE) + patchwork::plot_layout(ncol = 4, nrow = 4) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"AllSpots_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,20)
	  pp = SpatialFeaturePlot(this_subspa_all, features = cs, crop=TRUE) + patchwork::plot_layout(ncol = 4, nrow = 4)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off() 
	}
	pdf(file = paste0(this_OutDir,"AllSpots_","Images_SpatialPlot.pdf"),20,20)
	pp = SpatialFeaturePlot(this_subspa_all,features = cs, crop=TRUE, alpha=0) + patchwork::plot_layout(ncol = 4, nrow = 4)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"AllSpots_csMECO.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()
	# Only solid samples
	solid_samples = c( "R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scores), subset = Sample %in% solid_samples)
	for (cs in c(names(cs_signatures),"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs] = scores[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"AllSpots_SolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),10,5)
	  print(SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"AllSpots_SolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),10,5)
	  pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples), features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"AllSpots_SolidSamples_","Images_SpatialPlot.pdf"),10,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",solid_samples),features = cs, crop=TRUE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"AllSpots_SolidSamples_csMECO.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()
	# Only lepidic and solid samples
	lepidic_solid_samples = c("R26097-1F-Z2","R26097-1E-Z2","R24586-3D-Z1","R24586-3D-Z2" )
	this_subspa_all = subset(subspa_all,cells=rownames(scores), subset = Sample %in% lepidic_solid_samples)
	for (cs in c(names(cs_signatures),"rctd3_alveolar_macro","rctd3_mono_derived_macro","rctd3_CD8.CM","rctd3_CD8.TEX") ){
	  dcat(cs)
	  this_subspa_all@meta.data[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs] = scores[intersect(rownames(scores),rownames(this_subspa_all@meta.data)),cs]
	  min_score = min(this_subspa_all@meta.data[,cs])
	  max_score = max(this_subspa_all@meta.data[,cs])
	  pdf(file = paste0(this_OutDir,"AllSpots_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_Autoscaled.pdf"),20,5)
	  print(SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1) )
	  dev.off()  
	  colorz = rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))
	  pdf(file = paste0(this_OutDir,"AllSpots_LepidicSolidSamples_",cs,"_SpatialFeaturePlot_SameScale.pdf"),20,5)
	  pp = SpatialFeaturePlot(this_subspa_all, images=gsub("-",".",lepidic_solid_samples),features = cs, crop=TRUE) + patchwork::plot_layout(nrow = 1)
	  print(pp & scale_fill_gradientn(colours=colorz,limits = c(min_score,max_score)) )
	  dev.off()
	}
	pdf(file = paste0(this_OutDir,"AllSpots_LepidicSolidSamples_","Images_SpatialPlot.pdf"),20,5)
	pp = SpatialFeaturePlot(this_subspa_all,images=gsub("-",".",lepidic_solid_samples),features = cs, crop=TRUE, alpha=0) + patchwork::plot_layout(nrow = 1)
	print(pp)
	dev.off()
	cordf = cor(this_subspa_all@meta.data[,names(cs_signatures)],method = "pearson")
	pdf(paste0(this_OutDir,"AllSpots_LepidicSolidSamples_csMECO.pdf"),5,5)
	corrplot(cordf,order='AOE',col=colorz_solid,tl.col="black",addCoef.col="black",)
	dev.off()
}

visium_selected_boxplots = function( OutDir,signature="ThreeStates" ){
	this_OutDir = paste0(OutDir,"AmsCentered_canceronlygenesFALSE_ThreeStates/")
	load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))	
	load(file = paste0(this_OutDir,"tps_scores.RData"))
	md = subspa_all@meta.data
	## CancerOnly
	DecDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Deconvolution/"
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024/",Sample,"_thetadf.RData")) 
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall = thetadf
	    first = FALSE
	  } else {
	    thetall = rbind(thetall,thetadf)
	  }
	}
	colnames(thetall) = paste0( "rctd2_",colnames(thetall) )
	scores = cbind(scores,thetall[rownames(scores),])
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024_fine_reference/",Sample,"_thetadf.RData"))  
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall_fine = thetadf
	    first = FALSE
	  } else {
	    thetall_fine = rbind(thetall_fine,thetadf)
	  }
	}
	colnames(thetall_fine) = paste0( "rctd3_",colnames(thetall_fine) )
	scores = cbind(scores,thetall_fine[rownames(scores),])
	scores$CellState = names(cs_signatures)[apply(scores[,names(cs_signatures)],1,which.max)]
	if (signature=="ThreeStates"){ 
		scores$CellState = factor(scores$CellState,levels=c( "Alveolar","Proliferative","Hypoxic" )) 
		state_colorz = c("Alveolar"="dodgerblue4","Proliferative"="firebrick3","Hypoxic"="sienna4")
	}
	if (signature=="TwoStates"){ 
		scores$CellState = factor(scores$CellState,levels=c( "cs1","cs2" )) 
		state_colorz = c("cs1"="dodgerblue4","cs2"="firebrick4")
	}
	thresh = 0.1
	scc = scores[scores$rctd2_Cancer>thresh,]
	lepidic_solid_samples = c("R26097-1F-Z2","R26097-1E-Z2","R24586-3D-Z1","R24586-3D-Z2" )
	scc = scc[substr(rownames(scc),1,12) %in% lepidic_solid_samples,]
	scc$Pattern = ifelse(substr(rownames(scc),1,12) %in% c("R26097-1F-Z2","R26097-1E-Z2"),"lepidic","solid" )
	# comparing cs scores
	tscc = scc[,c("Pattern","Alveolar","Proliferative","Hypoxic" )]
	tscc = melt(tscc)
	pdf( paste0(OutDir,"Visium_LvsS_cs.pdf"),2,1.5 )
	x = factor(tscc$variable,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	y = tscc$value
	fill = factor(tscc$Pattern,levels=c("lepidic","solid" ))
	plot=dan.boxplots.multipages( x, y, fill=fill,xlab = "", ylab = "Cell state score", filllab = "Pattern", plotTitle = "", signifTest = NULL, ylimLeft = min(y), ylimRight = 0.65,comparisons = NULL, labelycoo = 1, fillColors = c( "dodgerblue4","firebrick3" ), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()

	tscc = scc[,c("Pattern","rctd3_CD8.CM","rctd3_CD8.TEX","rctd3_alveolar_macro","rctd3_mono_derived_macro" )]
	colnames(tscc) = c("Pattern","CD8+ T cells\nCentral memory","CD8+ T cells\nExhausted","Macrophages\nAlveolar","Macrophages\nMonocyte-derived"  )
	tscc = melt(tscc)
	pdf( paste0(OutDir,"Visium_LvsS_tme.pdf"),2.5,1.8 )
	x = factor(tscc$variable)
	y = log10(tscc$value)
	fill = factor(tscc$Pattern,levels=c("lepidic","solid" ))
	plot=dan.boxplots.multipages( x, y, fill=fill,xlab = "", ylab = "Proportion, log10", filllab = "Pattern", plotTitle = "", signifTest = NULL, ylimLeft = -6, ylimRight = max(y),comparisons = NULL, labelycoo = 1, fillColors = c( "dodgerblue4","firebrick3" ), jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
	print(plot)
	dev.off()
}

visium_tme_correlation = function( OutDir ){
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	this_OutDir = paste0(OutDir,"tme_correlation/")
	dir.create(this_OutDir)
	# vs cs
	load(file = paste0(OutDir,"AmsCentered_canceronlygenesFALSE_ThreeStates/tps_scores.RData"))
	## CancerOnly
	DecDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Deconvolution/"
	first = TRUE
	for (Sample in batch1_table$Sample){
	  load(file = paste0(DecDir,"RCTD_Feb2024/",Sample,"_thetadf.RData")) 
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall = thetadf
	    first = FALSE
	  } else {
	    thetall = rbind(thetall,thetadf)
	  }
	}
	colnames(thetall) = paste0( "rctd2_",colnames(thetall) )
	scores = cbind(scores,thetall[rownames(scores),])
	first = TRUE
	for (Sample in batch1_table$Sample){
		rownames(tissue) = paste0(Sample,"_",tissue$barcode)
	  load(file = paste0(DecDir,"RCTD_Feb2024_fine_reference/",Sample,"_thetadf.RData"))  
	  rownames(thetadf) = paste0(Sample,"_",rownames(thetadf))
	  if (first){
	    thetall_fine = thetadf
	    first = FALSE
	  } else {
	    thetall_fine = rbind(thetall_fine,thetadf)
	  }
	}
	colnames(thetall_fine) = paste0( "rctd3_",colnames(thetall_fine) )
	scores = cbind(scores,thetall_fine[rownames(scores),])

	thresh = 0.1
	scc = scores[scores$rctd2_Cancer>thresh,]

	rctd2_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd2"]
	rctd3_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd3"]

	for (Sample in batch1_table$Sample){
		dcat(Sample)
		tissue = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/NucleiDetection_StarDist/NucleiPerSpot/",Sample,"_NucleiPerSpot.txt") )
		rownames(tissue) = paste0(Sample,"_",tissue$barcode)
		tissue = tissue[intersect(rownames(tissue),rownames(scores)),]
		dd = as.matrix(dist(tissue[,c( "um_x","um_y" )]))
		dist_limit = 136 # distance to 6-NN should be 100, but sometimes it is more. Distance to the next neighbour should be 172, but sometimes it is less. 136 is the midpoint.
		ddbin = dd[intersect(rownames(scc),rownames(dd)),]
		ddbin[ddbin<dist_limit] = 1		
		ddbin[ddbin>=dist_limit] = 0
		if ((min(rowSums(ddbin))<1) | (max(rowSums(ddbin))>7) ){ dcat( "Problems!!",1 ) }
		scc[rownames(ddbin),"nn_spots"] = rowSums(ddbin)
		for ( cc in c(rctd2_cols,rctd3_cols) ){
			scc[rownames(ddbin),paste0( "mean_",cc )] = sapply(rownames(ddbin), FUN=function(x) mean(scores[colnames(ddbin)[as.numeric(ddbin[x,])==1],cc]) )
		}
	}
	# check that 1st-hop-averaged proportions kind-of correlate with spot proportions
	ccdf = data.frame(row.names=c(rctd2_cols,rctd3_cols),cc=c(rctd2_cols,rctd3_cols),cor=NA)
	for ( cc in c(rctd2_cols) ){
		dcat(paste0( cc,", cor = ",signif(cor(scc[,cc],scc[,paste0("mean_",cc)])) ))
		ccdf[cc,"cor"] = cor(scc[,cc],scc[,paste0("mean_",cc)])
	}
	for ( cc in c(rctd3_cols) ){
		dcat(paste0( cc,", cor = ",signif(cor(scc[,cc],scc[,paste0("mean_",cc)])) ))
		ccdf[cc,"cor"] = cor(scc[,cc],scc[,paste0("mean_",cc)])
	} # the strongest correlating: AT2, fibroblasts, AT1, B, cancer, mo_mac. the weakest: T cells, NK cells (yet still > 0.39). Very reproducible between Valais and DeZuani

	# Adding tps scores
	names_cs_signatures = c( "Alveolar","Proliferative","Hypoxic" )
	load(paste0(OutDir,"../tps_discovery/tps.RData"))
	scores_save = scores
	load(file = paste0(OutDir,"tps_AmsCentered_canceronlygenesFALSE/tps_scores.RData"))
	scc = cbind(scc,scores[rownames(scc),names(tps)])
	save(scc,file=paste0( this_OutDir,"scc.RData" ))
	load(file=paste0( this_OutDir,"scc.RData" ))

	library(colorRamp2)
	library(ComplexHeatmap)
	
	### cor with spot-level proportions
	cordf = cor(scc[,names_cs_signatures],scc[,rctd2_cols], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	# cordf_pvals = dan.df(rownames(cordf),colnames(cordf),as_df=F)
	# for (rn in rownames(cordf_pvals)){
	# 	for (cn in colnames(cordf_pvals)){
	# 		cordf_pvals[rn,cn] = cor.test(scc[,rn],scc[,cn])$p.value
	# 	}
	# }
	pdf(paste0(this_OutDir,"cor_cellstates_level2_spotlevel.pdf"),3,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = cor(scc[,names(tps)],scc[,rctd2_cols], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level2_spotlevel.pdf"),5,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = cor(scc[,names_cs_signatures],scc[,rctd3_cols], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level3_spotlevel.pdf"),3,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = cor(scc[,names(tps)],scc[,rctd3_cols], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level3_spotlevel.pdf"),6,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	### cor with nei-averaged proportions
	cordf = cor(scc[,names_cs_signatures],scc[,paste0("mean_",rctd2_cols)], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level2_neiAveraged.pdf"),3,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = cor(scc[,names(tps)],scc[,paste0("mean_",rctd2_cols)], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level2_neiAveraged.pdf"),5,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = cor(scc[,names_cs_signatures],scc[,paste0("mean_",rctd3_cols)], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level3_neiAveraged.pdf"),3,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = cor(scc[,names(tps)],scc[,paste0("mean_",rctd3_cols)], method = 'pearson')
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level3_neiAveraged.pdf"),6,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	### ppcor with spot-level proportions
	library(ppcor)
	cordf = dan.df(names_cs_signatures,rctd2_cols)
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level2_spotlevel_ppcor.pdf"),3,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = dan.df(names(tps),rctd2_cols)
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level2_spotlevel_ppcor.pdf"),5,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = dan.df(names_cs_signatures,rctd3_cols)
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level3_spotlevel_ppcor.pdf"),3,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = dan.df(names(tps),rctd3_cols)
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level3_spotlevel_ppcor.pdf"),6,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	### ppcor with nei-averaged proportions
	library(ppcor)
	cordf = dan.df(names_cs_signatures,paste0("mean_",rctd2_cols))
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level2_neiAveraged_ppcor.pdf"),3,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = dan.df(names(tps),paste0("mean_",rctd2_cols))
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level2_neiAveraged_ppcor.pdf"),5,7)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = dan.df(names_cs_signatures,paste0("mean_",rctd3_cols))
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_cellstates_level3_neiAveraged_ppcor.pdf"),3,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	cordf = dan.df(names(tps),paste0("mean_",rctd3_cols))
	for (rn in rownames(cordf)){
		for (cn in colnames(cordf)){
			cordf[rn,cn] = pcor.test(scc[,rn],scc[,cn],scc$rctd2_Cancer,method='pearson')$estimate
		}
	}
	cordf = t(cordf)
	# cordf = t(reorderMatt(cordf,columns_only=T))
	pdf(paste0(this_OutDir,"cor_tps_level3_neiAveraged_ppcor.pdf"),6,8)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	names_cs_signatures = c( "Alveolar","Proliferative","Hypoxic" )
	load(file=paste0( this_OutDir,"scc.RData" ))
	rctd2_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd2"]
	rctd3_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd3"]
	this_OutDir2 = paste0(this_OutDir, "same_as_sc/")
	dir.create(this_OutDir2)	
	library(compositions)
	# transforming proportions to match single-cell-like
	load(file = paste0("data/cell_typing/non_malignant/final_seus_annots/kim_annot_annd.RData"))
	all(rctd2_cols[rctd2_cols!="rctd2_Cancer"] %in% paste0( "rctd2_",unique(annot$annd_level_2)))
	all(rctd3_cols[rctd3_cols!="rctd3_Cancer"] %in% paste0( "rctd3_",unique(annot$annd_level_3)))
	load(file=paste0( this_OutDir,"scc.RData" ))
	# vs level 1
	ta = annot
	scc$mean_rctd1_epithelial = rowSums( scc[,paste0( "mean_",c( "rctd2_AT0","rctd2_AT1","rctd2_AT2","rctd2_basal","rctd2_club","rctd2_ciliated","rctd2_preTB" ))] )
	scc$mean_rctd1_Cancer = scc$mean_rctd2_Cancer
	scc$mean_rctd1_immune = rowSums( scc[,paste0( "mean_",c( "rctd2_Bcells","rctd2_DC","rctd2_macrophages","rctd2_mast","rctd2_monocytes","rctd2_neutrophils","rctd2_NKcells","rctd2_Tcells" ))] )
	scc$mean_rctd1_endothelial = rowSums( scc[,paste0( "mean_",c( "rctd2_EC_lymphatic","rctd2_EC_blood" ))] )
	scc$mean_rctd1_stroma = rowSums( scc[,paste0( "mean_",c( "rctd2_fibroblasts","rctd2_mesothelial","rctd2_smooth_muscle" ))] )
	columnz = paste0( "mean_rctd1_", c( "epithelial","Cancer","immune","endothelial","stroma" ))
	scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")

	cordf = t(cor(scc[,names_cs_signatures],scc[,columnz], method = 'spearman'))
	pdf(paste0(this_OutDir2,"cor_cellstates_level1_neiAveraged.pdf"),5,3)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	
	columnz = paste0( "mean_rctd1_", c( "epithelial","immune","endothelial","stroma" ))
	scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")
	cordf = cor(scc[,names_cs_signatures],scc[,columnz], method = 'spearman')
	cordf_pvals = dan.df(rownames(cordf),colnames(cordf),as_df = F)
	for (rn in rownames(cordf_pvals)){
		for (cn in colnames(cordf_pvals)){
			cordf_pvals[rn,cn] = cor.test(scc[,rn],scc[,cn], method = 'spearman')$p.value
		}
	}
	cordf_pvals = pvals_to_signLevels( cordf_pvals )
	pdf(paste0(this_OutDir2,"cor_cellstates_level1_neiAveraged_noCancer.pdf"),4,2)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_side = "bottom",row_names_side = "left",column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),,cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))}))
	draw(plot,heatmap_legend_side="right")
	dev.off()

	clrdf = scc
	for (rn in rownames(clrdf)){
	  clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz]))
	}
	cordf = cor(clrdf[,names_cs_signatures],clrdf[,columnz], method = 'spearman')
	cordf = cordf[,c( "mean_rctd1_epithelial","mean_rctd1_endothelial","mean_rctd1_immune","mean_rctd1_stroma")]
	cordf_pvals = dan.df(rownames(cordf),colnames(cordf),as_df = F)
	for (rn in rownames(cordf_pvals)){
		for (cn in colnames(cordf_pvals)){
			cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn], method = 'spearman')$p.value
		}
	}
	cordf_pvals = pvals_to_signLevels( cordf_pvals )
	pdf(paste0(this_OutDir2,"corClr_cellstates_level1_neiAveraged_noCancer.pdf"),4,2)
	plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_side = "bottom",row_names_side = "left",column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="CLR-transformed cell\nproportions, Spearman R",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),,cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))}))
	draw(plot,heatmap_legend_side="right")
	dev.off()
	for (ctlevel1 in c( "epithelial","immune","endothelial","stroma" )){
		dcat(ctlevel1)
		ta = annot[annot$annd_level_1==ctlevel1,]
		columnz = intersect(paste0( "mean_rctd2_", unique(ta$annd_level_2)),paste0("mean_",rctd2_cols) )
		scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")
		cordf = t(cor(scc[,names_cs_signatures],scc[,columnz], method = 'spearman'))
		pdf(paste0(this_OutDir2,"cor_cellstates_level2_",ctlevel1,"_neiAveraged.pdf"),5,3)
		plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
		draw(plot,heatmap_legend_side="right")
		dev.off()
		clrdf = scc
		for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
		cordf = cor(clrdf[,names_cs_signatures],clrdf[,columnz], method = 'spearman')
		cordf = cordf[,order( -as.numeric(cordf["Alveolar",])+as.numeric(cordf["Hypoxic",]) )]
		cordf_pvals = dan.df(rownames(cordf),colnames(cordf),as_df = F)
		for (rn in rownames(cordf_pvals)){
			for (cn in colnames(cordf_pvals)){
				cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn], method = 'spearman')$p.value
			}
		}
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		pdf(paste0(this_OutDir2,"corClr_cellstates_level2_",ctlevel1,"_neiAveraged_noCancer.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_side = "bottom",row_names_side = "left",column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="CLR-transformed cell\nproportions, Spearman R",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),,cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))}))
		draw(plot,heatmap_legend_side="right")
		dev.off()
	}
	for (ctlevel2 in c( "DC","macrophages","monocytes","Tcells","Bcells" )){
		dcat(ctlevel2)
		ta = annot[annot$annd_level_2==ctlevel2,]
		columnz = intersect(paste0( "rctd3_", unique(ta$annd_level_3)),rctd3_cols)
		scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")
		cordf = t(cor(scc[,names_cs_signatures],scc[,columnz], method = 'spearman'))
		pdf(paste0(this_OutDir2,"cor_cellstates_level3_",ctlevel2,"_neiAveraged.pdf"),5,3)
		plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2.5,height=unit(nrow(cordf), "cm")/2.5,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="Pearson R, single cell-level",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "top",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6)))
		draw(plot,heatmap_legend_side="right")
		dev.off()
		clrdf = scc
		for (rn in rownames(clrdf)){ clrdf[rn,columnz] = clr(as.numeric(clrdf[rn,columnz])) }
		cordf = cor(clrdf[,names_cs_signatures],clrdf[,columnz], method = 'spearman')
		cordf = cordf[,order( -as.numeric(cordf["Alveolar",])+as.numeric(cordf["Hypoxic",]) )]
		cordf_pvals = dan.df(rownames(cordf),colnames(cordf),as_df = F)
		for (rn in rownames(cordf_pvals)){
			for (cn in colnames(cordf_pvals)){
				cordf_pvals[rn,cn] = cor.test(clrdf[,rn],clrdf[,cn], method = 'spearman')$p.value
			}
		}
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		pdf(paste0(this_OutDir2,"corClr_cellstates_level3_",ctlevel2,"_neiAveraged_noCancer.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_side = "bottom",row_names_side = "left",column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(title_position="leftcenter-rot",title="CLR-transformed cell\nproportions, Spearman R",direction = "vertical",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),,cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))}))
		draw(plot,heatmap_legend_side="right")
		dev.off()
	}
	
	# this_OutDir2 = paste0(this_OutDir, "extreme_spots/")
	# dir.create(this_OutDir2)
	# # extracting spots that are high (>percentile cutoff) for a cell state, and low for both the other two
	# cutoff = "25"
	# alow = quantile(scc$Alveolar)[paste0(cutoff,"%")]
	# ahigh = -quantile(-scc$Alveolar)[paste0(cutoff,"%")]
	# plow = quantile(scc$Proliferative)[paste0(cutoff,"%")]
	# phigh = -quantile(-scc$Proliferative)[paste0(cutoff,"%")]
	# hlow = quantile(scc$Hypoxic)[paste0(cutoff,"%")]
	# hhigh = -quantile(-scc$Hypoxic)[paste0(cutoff,"%")]
	# spotsA = rownames(scc[((scc$Alveolar>ahigh) & (scc$Proliferative<plow)) & (scc$Hypoxic<hlow),])
	# spotsP = rownames(scc[((scc$Proliferative>phigh) & (scc$Alveolar<alow)) & (scc$Hypoxic<hlow),])
	# spotsH = rownames(scc[((scc$Hypoxic>hhigh) & (scc$Proliferative<plow)) & (scc$Alveolar<alow),])
	# scce = scc[c(spotsA,spotsP,spotsH),]
	# scce[spotsA,"extreme_for"] = "Alveolar"
	# scce[spotsP,"extreme_for"] = "Proliferative"
	# scce[spotsH,"extreme_for"] = "Hypoxic"
	# dtable(scce$extreme_for,scce$Patient)
	# aggregate(rctd2_Cancer~extreme_for, data=scce, FUN='mean')

	# ## level 1
	# columnz = paste0( "mean_rctd1_", c( "epithelial","immune","endothelial","stroma" ))
	# scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# sccem = melt(scce[,c( "extreme_for",columnz )])
	# fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level1.pdf" )
	# x = sccem$variable
	# y = sccem$value
	# fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	# dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )

	# for (ctlevel1 in c( "epithelial","immune","endothelial","stroma" )){
	# 	dcat(ctlevel1)
	# 	ta = annot[annot$annd_level_1==ctlevel1,]
	# 	columnz = intersect(paste0( "mean_rctd2_", unique(ta$annd_level_2)),paste0("mean_",rctd2_cols) )
	# 	scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 	sccem = melt(scce[,c( "extreme_for",columnz )])
	# 	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level2_",ctlevel1,".pdf" )
	# 	x = sccem$variable
	# 	y = sccem$value
	# 	fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 	fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 	colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	# 	dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
	# }
	# for (ctlevel2 in c( "DC","macrophages","monocytes","Tcells","Bcells" )){
	# 	dcat(ctlevel2)
	# 	ta = annot[annot$annd_level_2==ctlevel2,]
	# 	columnz = intersect(paste0( "rctd3_", unique(ta$annd_level_3)),rctd3_cols)
	# 	scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 	sccem = melt(scce[,c( "extreme_for",columnz )])
	# 	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level3_",ctlevel2,".pdf" )
	# 	x = sccem$variable
	# 	y = sccem$value
	# 	fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 	fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 	colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	# 	dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
	# }

	# this_OutDir2 = paste0(this_OutDir, "extreme_spots_PatientWise_PvsH/")
	# dir.create(this_OutDir2)
	# load(file=paste0( this_OutDir,"scc.RData" ))
	# names_cs_signatures = c( "Alveolar","Proliferative","Hypoxic" )
	# rctd2_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd2"]
	# rctd3_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd3"]
	# load(file = paste0("data/cell_typing/non_malignant/final_seus_annots/kim_annot_annd.RData"))
	# all(rctd2_cols[rctd2_cols!="rctd2_Cancer"] %in% paste0( "rctd2_",unique(annot$annd_level_2)))
	# all(rctd3_cols[rctd3_cols!="rctd3_Cancer"] %in% paste0( "rctd3_",unique(annot$annd_level_3)))
	# # vs level 1
	# ta = annot
	# scc$mean_rctd1_epithelial = rowSums( scc[,paste0( "mean_",c( "rctd2_AT0","rctd2_AT1","rctd2_AT2","rctd2_basal","rctd2_club","rctd2_ciliated","rctd2_preTB" ))] )
	# scc$mean_rctd1_Cancer = scc$mean_rctd2_Cancer
	# scc$mean_rctd1_immune = rowSums( scc[,paste0( "mean_",c( "rctd2_Bcells","rctd2_DC","rctd2_macrophages","rctd2_mast","rctd2_monocytes","rctd2_neutrophils","rctd2_NKcells","rctd2_Tcells" ))] )
	# scc$mean_rctd1_endothelial = rowSums( scc[,paste0( "mean_",c( "rctd2_EC_lymphatic","rctd2_EC_blood" ))] )
	# scc$mean_rctd1_stroma = rowSums( scc[,paste0( "mean_",c( "rctd2_fibroblasts","rctd2_mesothelial","rctd2_smooth_muscle" ))] )
	# columnz = paste0( "mean_rctd1_", c( "epithelial","Cancer","immune","endothelial","stroma" ))
	# scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")

	# load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))	
	# sccSave = scc
	# for (patient in unique(scc$Patient)){
	# 	dcat(patient)
	# 	scc = sccSave[sccSave$Patient==patient,]
	# 	# extracting spots that are high (>percentile cutoff) for a cell state, and low for both the other two
	# 	cutoff = "25"
	# 	alow = quantile(scc$Alveolar)[paste0(cutoff,"%")]
	# 	spotsNotA = rownames(scc[(scc$Alveolar<alow),])
	# 	scce = scc[spotsNotA,]
	# 	spotsH = quantile(scce$Proliferative-scce$Hypoxic)
	# 	spotsH = rownames(scce[(scce$Proliferative-scce$Hypoxic)<spotsH["25%"],])
	# 	spotsP = quantile(scce$Proliferative-scce$Hypoxic)
	# 	spotsP = rownames(scce[(scce$Proliferative-scce$Hypoxic)>spotsP["75%"],])
	# 	scce = scc[c(spotsP,spotsH),]
	# 	scce[spotsP,"extreme_for"] = "Proliferative"
	# 	scce[spotsH,"extreme_for"] = "Hypoxic"
	# 	tsa = subset(subspa_all,cells=rownames(scce))
	# 	tsa@meta.data$extreme_for = scce[colnames(tsa),"extreme_for"]
	# 	these_images = names(tsa@images)[grepl(patient,names(tsa@images))]
	# 	pdf(paste0( this_OutDir2,"cutoff",cutoff,"_Patient",patient,"_extreme_for_DimPlot.pdf" ),7*length(these_images),7)
	# 	plot = SpatialDimPlot(tsa, group.by="extreme_for",images=these_images,ncol=length(these_images))
	# 	print(plot)
	# 	dev.off()

	# 	## level 1
	# 	columnz = paste0( "mean_rctd1_", c( "epithelial","immune","endothelial","stroma" ))
	# 	scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 	sccem = melt(scce[,c( "extreme_for",columnz )])
	# 	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level1_Patient",patient,".pdf" )
	# 	x = sccem$variable
	# 	y = sccem$value
	# 	fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 	fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 	colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	# 	dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = "wilcox", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )

	# 	for (ctlevel1 in c( "epithelial","immune","endothelial","stroma" )){
	# 		dcat(ctlevel1)
	# 		ta = annot[annot$annd_level_1==ctlevel1,]
	# 		columnz = intersect(paste0( "mean_rctd2_", unique(ta$annd_level_2)),paste0("mean_",rctd2_cols) )
	# 		scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 		sccem = melt(scce[,c( "extreme_for",columnz )])
	# 		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level2_",ctlevel1,"_Patient",patient,".pdf" )
	# 		x = sccem$variable
	# 		y = sccem$value
	# 		fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 		fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 		colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	# 		dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = "wilcox", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
	# 	}
	# 	for (ctlevel2 in c( "DC","macrophages","monocytes","Tcells","Bcells" )){
	# 		dcat(ctlevel2)
	# 		ta = annot[annot$annd_level_2==ctlevel2,]
	# 		columnz = intersect(paste0( "rctd3_", unique(ta$annd_level_3)),rctd3_cols)
	# 		scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 		sccem = melt(scce[,c( "extreme_for",columnz )])
	# 		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level3_",ctlevel2,"_Patient",patient,".pdf" )
	# 		x = sccem$variable
	# 		y = sccem$value
	# 		fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 		fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 		colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	# 		dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = "wilcox", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
	# 	}
	# }

	# this_OutDir2 = paste0(this_OutDir, "extreme_spots_PatientWise/")
	# dir.create(this_OutDir2)
	# load(file=paste0( this_OutDir,"scc.RData" ))
	# names_cs_signatures = c( "Alveolar","Proliferative","Hypoxic" )
	# rctd2_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd2"]
	# rctd3_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd3"]
	# load(file = paste0("data/cell_typing/non_malignant/final_seus_annots/kim_annot_annd.RData"))
	# all(rctd2_cols[rctd2_cols!="rctd2_Cancer"] %in% paste0( "rctd2_",unique(annot$annd_level_2)))
	# all(rctd3_cols[rctd3_cols!="rctd3_Cancer"] %in% paste0( "rctd3_",unique(annot$annd_level_3)))
	# # vs level 1
	# ta = annot
	# scc$mean_rctd1_epithelial = rowSums( scc[,paste0( "mean_",c( "rctd2_AT0","rctd2_AT1","rctd2_AT2","rctd2_basal","rctd2_club","rctd2_ciliated","rctd2_preTB" ))] )
	# scc$mean_rctd1_Cancer = scc$mean_rctd2_Cancer
	# scc$mean_rctd1_immune = rowSums( scc[,paste0( "mean_",c( "rctd2_Bcells","rctd2_DC","rctd2_macrophages","rctd2_mast","rctd2_monocytes","rctd2_neutrophils","rctd2_NKcells","rctd2_Tcells" ))] )
	# scc$mean_rctd1_endothelial = rowSums( scc[,paste0( "mean_",c( "rctd2_EC_lymphatic","rctd2_EC_blood" ))] )
	# scc$mean_rctd1_stroma = rowSums( scc[,paste0( "mean_",c( "rctd2_fibroblasts","rctd2_mesothelial","rctd2_smooth_muscle" ))] )
	# columnz = paste0( "mean_rctd1_", c( "epithelial","Cancer","immune","endothelial","stroma" ))
	# scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")

	# sccSave = scc
	# for (patient in unique(scc$Patient)){
	# 	dcat(patient)
	# 	scc = sccSave[sccSave$Patient==patient,]
	# 	# extracting spots that are high (>percentile cutoff) for a cell state, and low for both the other two
	# 	cutoff = "25"
	# 	alow = quantile(scc$Alveolar)[paste0(cutoff,"%")]
	# 	ahigh = -quantile(-scc$Alveolar)[paste0(cutoff,"%")]
	# 	plow = quantile(scc$Proliferative)[paste0(cutoff,"%")]
	# 	phigh = -quantile(-scc$Proliferative)[paste0(cutoff,"%")]
	# 	hlow = quantile(scc$Hypoxic)[paste0(cutoff,"%")]
	# 	hhigh = -quantile(-scc$Hypoxic)[paste0(cutoff,"%")]
	# 	spotsA = rownames(scc[((scc$Alveolar>ahigh) & (scc$Proliferative<plow)) & (scc$Hypoxic<hlow),])
	# 	spotsP = rownames(scc[((scc$Proliferative>phigh) & (scc$Alveolar<alow)) & (scc$Hypoxic<hlow),])
	# 	spotsH = rownames(scc[((scc$Hypoxic>hhigh) & (scc$Proliferative<plow)) & (scc$Alveolar<alow),])
	# 	scce = scc[c(spotsA,spotsP,spotsH),]
	# 	scce[spotsA,"extreme_for"] = "Alveolar"
	# 	scce[spotsP,"extreme_for"] = "Proliferative"
	# 	scce[spotsH,"extreme_for"] = "Hypoxic"
	# 	dtable(scce$extreme_for,scce$Patient)
	# 	aggregate(rctd2_Cancer~extreme_for, data=scce, FUN='mean')

	# 	## level 1
	# 	columnz = paste0( "mean_rctd1_", c( "epithelial","immune","endothelial","stroma" ))
	# 	scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 	fileName = paste0( this_OutDir2,"multipages_cutoff",cutoff,"_level1_Patient",patient,".pdf" )
	# 	x = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 	fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 	colorz = dan.expand_colors(scce$extreme_for,levels(fill),fillColorz)
	# 	pdf(fileName,3,2)
	# 	for (xx in columnz){
	# 		y = scce[,xx]
	# 		plot = dan.boxplots.multipages( x, y, xlab = "", ylab = "Ratios in state-rich spots", plotTitle = xx, signifTest = "kruskal", labelycoo = max(y), xColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
	# 		print(plot)
	# 	}
	# 	dev.off()
	# 	for (ctlevel1 in c( "epithelial","immune","endothelial","stroma" )){
	# 		dcat(ctlevel1)
	# 		ta = annot[annot$annd_level_1==ctlevel1,]
	# 		columnz = intersect(paste0( "mean_rctd2_", unique(ta$annd_level_2)),paste0("mean_",rctd2_cols) )
	# 		scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 		fileName = paste0( this_OutDir2,"multipages_cutoff",cutoff,"_level2_",ctlevel1,"_Patient",patient,".pdf" )
	# 		x = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 		fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 		colorz = dan.expand_colors(scce$extreme_for,levels(fill),fillColorz)
	# 		pdf(fileName,3,2)
	# 		for (xx in columnz){
	# 			y = scce[,xx]
	# 			plot = dan.boxplots.multipages( x, y, xlab = "", ylab = "Ratios in state-rich spots", plotTitle = xx, signifTest = "kruskal", labelycoo = max(y), xColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
	# 			print(plot)
	# 		}
	# 		dev.off()
	# 	}
	# 	for (ctlevel2 in c( "DC","macrophages","monocytes","Tcells","Bcells" )){
	# 		dcat(ctlevel2)
	# 		ta = annot[annot$annd_level_2==ctlevel2,]
	# 		columnz = intersect(paste0( "rctd3_", unique(ta$annd_level_3)),rctd3_cols)
	# 		scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	# 		fileName = paste0( this_OutDir2,"multipages_cutoff",cutoff,"_level3_",ctlevel2,"_Patient",patient,".pdf" )
	# 		x = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	# 		fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	# 		colorz = dan.expand_colors(scce$extreme_for,levels(fill),fillColorz)
	# 		pdf(fileName,3,2)
	# 		for (xx in columnz){
	# 			y = scce[,xx]
	# 			plot = dan.boxplots.multipages( x, y, xlab = "", ylab = "Ratios in state-rich spots", plotTitle = xx, signifTest = "kruskal", labelycoo = max(y), xColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
	# 			print(plot)
	# 		}
	# 		dev.off()
	# 	}
	# }
}

weighted_harmonic_mean = function( pvec, wvec ){
	return( sum(wvec)/sum(wvec/pvec) )
}

visium_tme_PvsH = function( OutDir, subspa_all, cutoff = 25, use1hopNN = FALSE ){
	library(compositions)
	library(metap)
	this_OutDir = paste0(OutDir,"tme_correlation/")
	dir.create(this_OutDir)
	this_OutDir2 = paste0(this_OutDir, "extreme_spots_PvsH/")
	dir.create(this_OutDir2)
	load(file=paste0( this_OutDir,"scc.RData" ))
	names_cs_signatures = c( "Alveolar","Proliferative","Hypoxic" )
	rctd2_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd2"]
	rctd3_cols = colnames(scc)[substr(colnames(scc),1,5)=="rctd3"]
	load(file = paste0("data/cell_typing/non_malignant/final_seus_annots/kim_annot_annd.RData"))
	all(rctd2_cols[rctd2_cols!="rctd2_Cancer"] %in% paste0( "rctd2_",unique(annot$annd_level_2)))
	all(rctd3_cols[rctd3_cols!="rctd3_Cancer"] %in% paste0( "rctd3_",unique(annot$annd_level_3)))
	if (use1hopNN){
		this_prefix = "mean_"
	} else {
		this_prefix = ""
	}
	# vs level 1
	ta = annot
	scc$mean_rctd1_epithelial = rowSums( scc[,paste0( this_prefix,c( "rctd2_AT0","rctd2_AT1","rctd2_AT2","rctd2_basal","rctd2_club","rctd2_ciliated","rctd2_preTB" ))] )
	scc$mean_rctd1_Cancer = scc[,paste0(this_prefix,"rctd2_Cancer")]
	scc$mean_rctd1_immune = rowSums( scc[,paste0( this_prefix,c( "rctd2_Bcells","rctd2_DC","rctd2_macrophages","rctd2_mast","rctd2_monocytes","rctd2_neutrophils","rctd2_NKcells","rctd2_Tcells" ))] )
	scc$mean_rctd1_endothelial = rowSums( scc[,paste0( this_prefix,c( "rctd2_EC_lymphatic","rctd2_EC_blood" ))] )
	scc$mean_rctd1_stroma = rowSums( scc[,paste0( this_prefix,c( "rctd2_fibroblasts","rctd2_mesothelial","rctd2_smooth_muscle" ))] )
	columnz = paste0( "mean_rctd1_", c( "epithelial","Cancer","immune","endothelial","stroma" ))
	scc[,columnz] = sweep(scc[,columnz],1,as.numeric(rowSums(scc[,columnz])),FUN="/")

	# extracting spots that are high (>percentile cutoff) for a cell state, and low for both the other two
	this_seq = seq(0,1,cutoff/100)
	alow = quantile(scc$Alveolar,this_seq)[paste0(cutoff,"%")]
	spotsNotA = rownames(scc[(scc$Alveolar<alow),])
	scce = scc[spotsNotA,]
	spotsH = quantile(scce$Proliferative-scce$Hypoxic,this_seq)
	spotsH = rownames(scce[(scce$Proliferative-scce$Hypoxic)<spotsH[paste0(cutoff,"%")],])
	spotsP = quantile(scce$Proliferative-scce$Hypoxic,this_seq)
	spotsP = rownames(scce[(scce$Proliferative-scce$Hypoxic)>spotsP[paste0(100-cutoff,"%")],])
	scce = scc[c(spotsP,spotsH),]
	scce[spotsP,"extreme_for"] = "Proliferative"
	scce[spotsH,"extreme_for"] = "Hypoxic"
	dtable(scce$extreme_for,scce$Patient)
	aggregate(rctd2_Cancer~extreme_for, data=scce, FUN='mean')
	tsa = subset(subspa_all,cells=rownames(scce))
	tsa@meta.data$extreme_for = scce[colnames(tsa),"extreme_for"]
	pdf(paste0( this_OutDir2,"cutoff",cutoff,"_allPatients_extreme_for_DimPlot.pdf" ),20,20)
	plot = SpatialDimPlot(tsa, group.by="extreme_for",ncol=4,crop=F)
	print(plot)
	dev.off()
	scceSave = scce

	## level 1
	columnz = paste0( "mean_rctd1_", c( "epithelial","immune","endothelial","stroma" ))
	scce = scceSave
	scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
	
	sccem = melt(scce[,c( "extreme_for",columnz )])
	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level1_use1hopNN",as.character(use1hopNN),".pdf" )
	x = sccem$variable
	y = sccem$value
	fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = "wilcox.test", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y), xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level1_multipages_use1hopNN",as.character(use1hopNN),".pdf" )
	dan.save(scce,fileName)
	pdf(fileName,3,2)
	for (xx in columnz){
		x = scce$Patient
		y = scce[,xx]
		fill = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
		plot = dan.boxplots.multipages( x, y, fill, xlab = "", ylab = "Ratios in state-rich spots",filllab="", plotTitle = xx, signifTest = NULL, labelycoo = max(y), fillColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
		print(plot)
	}
	dev.off()
	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level1_PairedDot_use1hopNN",as.character(use1hopNN),".pdf" )
	pdf(fileName,2,2)
	for (xx in columnz){
		x = scce$Patient
		y = scce[,xx]
		state = as.character(scce$extreme_for)
		pvecG = c()
		pvecL = c()
		directions = c()
		wvec = c()
		pairdf = dan.df(0,c("Patient","State","Value","Colors"))
		for (pat in unique(x)){
			if ( length(dtable(state[x==pat]))<2 ) { next }
			pvecG = c(pvecG,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='greater')$p.value)
			pvecL = c(pvecL,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='less')$p.value)
			directions = c(directions,ifelse((median(y[(x==pat) & (state=="Proliferative")])>median(y[(x==pat) & (state=="Hypoxic")]) ),1,-1))
			wvec = c(wvec,sum(x==pat))
			pairdf = rbind(pairdf,data.frame(Patient=pat,State=c( "Proliferative","Hypoxic" ),Value=c(median(y[(x==pat) & (state=="Proliferative")]),median(y[(x==pat) & (state=="Hypoxic")])),Colors=c( "#FBB040","#BE1E2D" ),stringsAsFactors=F))
		}
		pvalG = as.numeric(metap::sumz( pvecG, wvec )$p)
		pvalL = as.numeric(metap::sumz( pvecL, wvec )$p)
		pval_harmonic = signif(min(c(pvalG,pvalL),na.rm=T)*2,2)

		plot = dan.pairedDotPlot.multipages( factor(pairdf$State,levels=c( "Proliferative","Hypoxic" )), pairdf$Value, xlab="", ylab="Median patient-wise ratios", plotTitle = paste0(substr(xx,nchar(this_prefix)+12,nchar(xx)),", combined p = ",pval_harmonic), labelLines = pairdf$Patient, signifTest = NULL, comparisons = NULL, ylimLeft = min(pairdf$Value)-0.01, ylimRight = max(pairdf$Value)+0.01, labelycoo = 1, jitterColors = pairdf$Colors, labelPoints = "", linesPairings = pairdf$Patient )
		print(plot)
	}
	dev.off()

	for (ctlevel1 in c( "epithelial","immune","endothelial","stroma" )){
		dcat(ctlevel1)
		ta = annot[annot$annd_level_1==ctlevel1,]
		columnz = intersect(paste0( this_prefix,"rctd2_", unique(ta$annd_level_2)),paste0(this_prefix,rctd2_cols) )
		scce = scceSave
		scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
		sccem = melt(scce[,c( "extreme_for",columnz )])
		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level2_",ctlevel1,"_use1hopNN",as.character(use1hopNN),".pdf" )
		x = sccem$variable
		y = sccem$value
		fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
		fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
		colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
		dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = "wilcox.test", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y), xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level2_",ctlevel1,"_multipages_use1hopNN",as.character(use1hopNN),".pdf" )
		dan.save(scce,fileName)
		pdf(fileName,3,2)
		for (xx in columnz){
			x = scce$Patient
			y = scce[,xx]
			fill = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
			plot = dan.boxplots.multipages( x, y, fill, xlab = "", ylab = "Ratios in state-rich spots",filllab="", plotTitle = xx, signifTest = NULL, labelycoo = max(y), fillColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
			print(plot)
		}
		dev.off()
		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level2_",ctlevel1,"_PairedDot_use1hopNN",as.character(use1hopNN),".pdf" )
		pdf(fileName,2,2)
		for (xx in columnz){
			x = scce$Patient
			y = scce[,xx]
			state = as.character(scce$extreme_for)
			pvecG = c()
			pvecL = c()
			directions = c()
			wvec = c()
			pairdf = dan.df(0,c("Patient","State","Value","Colors"))
			for (pat in unique(x)){
				if ( length(dtable(state[x==pat]))<2 ) { next }
				pvecG = c(pvecG,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='greater')$p.value)
				pvecL = c(pvecL,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='less')$p.value)
				directions = c(directions,ifelse((median(y[(x==pat) & (state=="Proliferative")])>median(y[(x==pat) & (state=="Hypoxic")]) ),1,-1))
				wvec = c(wvec,sum(x==pat))
				pairdf = rbind(pairdf,data.frame(Patient=pat,State=c( "Proliferative","Hypoxic" ),Value=c(median(y[(x==pat) & (state=="Proliferative")]),median(y[(x==pat) & (state=="Hypoxic")])),Colors=c( "#FBB040","#BE1E2D" ),stringsAsFactors=F))
			}
			pvalG = signif(as.numeric(metap::sumz( pvecG, wvec )$p),2)
			pvalL = signif(as.numeric(metap::sumz( pvecL, wvec )$p),2)
			pval_harmonic = min(c(pvalG,pvalL),na.rm=T)*2
			plot = dan.pairedDotPlot.multipages( factor(pairdf$State,levels=c( "Proliferative","Hypoxic" )), pairdf$Value, xlab="", ylab="Median patient-wise ratios", plotTitle = paste0(substr(xx,nchar(this_prefix)+7,nchar(xx)),", combined p = ",pval_harmonic), labelLines = pairdf$Patient, signifTest = NULL, comparisons = NULL, ylimLeft = min(pairdf$Value)-0.01, ylimRight = max(pairdf$Value)+0.01, labelycoo = 1, jitterColors = pairdf$Colors, labelPoints = "", linesPairings = pairdf$Patient )
			print(plot)
		}
		dev.off()
	}
	for (ctlevel2 in c( "DC","macrophages","monocytes","Tcells","Bcells" )){
		dcat(ctlevel2)
		ta = annot[annot$annd_level_2==ctlevel2,]
		columnz = intersect(paste0( this_prefix,"rctd3_", unique(ta$annd_level_3)),paste0(this_prefix,rctd3_cols) )
		scce = scceSave
		scce[,columnz] = sweep(scce[,columnz],1,as.numeric(rowSums(scce[,columnz])),FUN="/")
		sccem = melt(scce[,c( "extreme_for",columnz )])
		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level3_",ctlevel2,"_use1hopNN",as.character(use1hopNN),".pdf" )
		x = sccem$variable
		y = sccem$value
		fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
		fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
		colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
		dan.boxplots( fileName, x, y, fill, xlab = "", ylab = "Ratios in state-rich spots", filllab = "", plotTitle = "", signifTest = "wilcox.test", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y), xColors = "black", fillColors = fillColorz, jitterColors = colorz, labelJitteredPoints = NULL, jitterDotSize = 0.01, fileWidth = 5, fileHeight = 2, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level3_",ctlevel2,"_multipages_use1hopNN",as.character(use1hopNN),".pdf" )
		dan.save(scce,fileName)
		pdf(fileName,3,2)
		for (xx in columnz){
			x = scce$Patient
			y = scce[,xx]
			fill = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
			plot = dan.boxplots.multipages( x, y, fill, xlab = "", ylab = "Ratios in state-rich spots",filllab="", plotTitle = xx, signifTest = NULL, labelycoo = max(y), fillColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
			print(plot)
		}
		dev.off()
		fileName = paste0( this_OutDir2,"cutoff",cutoff,"_level3_",ctlevel2,"_PairedDot_use1hopNN",as.character(use1hopNN),".pdf" )
		pdf(fileName,2,2)
		for (xx in columnz){
			x = scce$Patient
			y = scce[,xx]
			state = as.character(scce$extreme_for)
			pvecG = c()
			pvecL = c()
			directions = c()
			wvec = c()
			pairdf = dan.df(0,c("Patient","State","Value","Colors"))
			for (pat in unique(x)){
			if ( length(dtable(state[x==pat]))<2 ) { next }
				pvecG = c(pvecG,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='greater')$p.value)
				pvecL = c(pvecL,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='less')$p.value)
				directions = c(directions,ifelse((median(y[(x==pat) & (state=="Proliferative")])>median(y[(x==pat) & (state=="Hypoxic")]) ),1,-1))
				wvec = c(wvec,sum(x==pat))
				pairdf = rbind(pairdf,data.frame(Patient=pat,State=c( "Proliferative","Hypoxic" ),Value=c(median(y[(x==pat) & (state=="Proliferative")]),median(y[(x==pat) & (state=="Hypoxic")])),Colors=c( "#FBB040","#BE1E2D" ),stringsAsFactors=F))
			}
			pvalG = signif(as.numeric(metap::sumz( pvecG, wvec )$p),2)
			pvalL = signif(as.numeric(metap::sumz( pvecL, wvec )$p),2)
			pval_harmonic = min(c(pvalG,pvalL),na.rm=T)*2
			plot = dan.pairedDotPlot.multipages( factor(pairdf$State,levels=c( "Proliferative","Hypoxic" )), pairdf$Value, xlab="", ylab="Median patient-wise ratios", plotTitle = paste0(substr(xx,nchar(this_prefix)+7,nchar(xx)),", combined p = ",pval_harmonic), labelLines = pairdf$Patient, signifTest = NULL, comparisons = NULL, ylimLeft = min(pairdf$Value)-0.01, ylimRight = max(pairdf$Value)+0.01, labelycoo = 1, jitterColors = pairdf$Colors, labelPoints = "", linesPairings = pairdf$Patient )
			print(plot)
		}
		dev.off()
	}

	gene_list = c( "IGHD","IGHE","IGHM","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGLV3-1","IGLV6-57","IGKV4-1","IGKV1-12","IGLC7","IGLL5","CD38","MZB1","DERL3","JSRP1","TNFRSF17","SLAMF7" )
	ctlevel2 = "Bcells"
	ta = annot[annot$annd_level_2==ctlevel2,]
	scce = scceSave
	a = subspa_all@assays$SCT@data[,rownames(scce)]
	a = a[rowSums(a>0)>0,]
	gene_list = intersect(rownames(a),gene_list)
	scce = cbind(scce,t(as.matrix(a[gene_list,rownames(scce)])))
	sccem = melt(scce[,c( "extreme_for",gene_list )])
	x = sccem$variable
	y = sccem$value
	fill = factor(sccem$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
	fillColorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
	colorz = dan.expand_colors(sccem$extreme_for,levels(fill),fillColorz)
	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_PlasmaCells.pdf" )
	dan.save(scce,fileName)
	pdf(fileName,3,2)
	for (xx in gene_list){
		x = scce$Patient
		y = scce[,xx]
		fill = factor(scce$extreme_for,levels=c( "Alveolar","Proliferative","Hypoxic" ))
		plot = dan.boxplots.multipages( x, y, fill, xlab = "", ylab = "Expression in spots",filllab="", plotTitle = xx, signifTest = NULL, labelycoo = max(y), fillColors = fillColorz, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
		print(plot)
	}
	dev.off()
	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_PlasmaCells_PairedDotPlot.pdf" )
	pdf(fileName,2,2)
	for (xx in gene_list){
		x = scce$Patient
		y = scce[,xx]
		state = as.character(scce$extreme_for)
		pvecG = c()
		pvecL = c()
		directions = c()
		wvec = c()
		pairdf = dan.df(0,c("Patient","State","Value","Colors"))
		for (pat in unique(x)){
		if ( length(dtable(state[x==pat]))<2 ) { next }
			pvecG = c(pvecG,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='greater')$p.value)
			pvecL = c(pvecL,wilcox.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='less')$p.value)
			directions = c(directions,ifelse((median(y[(x==pat) & (state=="Proliferative")])>median(y[(x==pat) & (state=="Hypoxic")]) ),1,-1))
			wvec = c(wvec,sum(x==pat))
			pairdf = rbind(pairdf,data.frame(Patient=pat,State=c( "Proliferative","Hypoxic" ),Value=c(median(y[(x==pat) & (state=="Proliferative")]),median(y[(x==pat) & (state=="Hypoxic")])),Colors=c( "#FBB040","#BE1E2D" ),stringsAsFactors=F))
		}
		pvalG = signif(as.numeric(metap::sumz( pvecG, wvec )$p),2)
		pvalL = signif(as.numeric(metap::sumz( pvecL, wvec )$p),2)
		pval_harmonic = min(c(pvalG,pvalL),na.rm=T)*2
		plot = dan.pairedDotPlot.multipages( factor(pairdf$State,levels=c( "Proliferative","Hypoxic" )), pairdf$Value, xlab="", ylab="Median patient-wise expression", plotTitle = paste0(xx,", combined p = ",pval_harmonic), labelLines = pairdf$Patient, signifTest = NULL, comparisons = NULL, ylimLeft = min(pairdf$Value)-0.01, ylimRight = max(pairdf$Value)+0.01, labelycoo = 1, jitterColors = pairdf$Colors, labelPoints = "", linesPairings = pairdf$Patient )
		print(plot)
	}
	dev.off()
	fileName = paste0( this_OutDir2,"cutoff",cutoff,"_PlasmaCells_PairedDotPlot_means.pdf" )
	pdf(fileName,2,2)
	for (xx in gene_list){
		x = scce$Patient
		y = scce[,xx]
		state = as.character(scce$extreme_for)
		pvecG = c()
		pvecL = c()
		directions = c()
		wvec = c()
		pairdf = dan.df(0,c("Patient","State","Value","Colors"))
		for (pat in unique(x)){
		if ( length(dtable(state[x==pat]))<2 ) { next }
			pvecG = c(pvecG,t.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='greater')$p.value)
			pvecL = c(pvecL,t.test(y[(x==pat) & (state=="Proliferative")],y[(x==pat) & (state=="Hypoxic")],alternative='less')$p.value)
			directions = c(directions,ifelse((mean(y[(x==pat) & (state=="Proliferative")])>mean(y[(x==pat) & (state=="Hypoxic")]) ),1,-1))
			wvec = c(wvec,sum(x==pat))
			pairdf = rbind(pairdf,data.frame(Patient=pat,State=c( "Proliferative","Hypoxic" ),Value=c(mean(y[(x==pat) & (state=="Proliferative")]),mean(y[(x==pat) & (state=="Hypoxic")])),Colors=c( "#FBB040","#BE1E2D" ),stringsAsFactors=F))
		}
		wvec = wvec[!is.na(pvecG)]
		pvecG = pvecG[!is.na(pvecG)]
		pvecL = pvecL[!is.na(pvecL)]
		pvalG = signif(as.numeric(metap::sumz( pvecG, wvec )$p),2)
		pvalL = signif(as.numeric(metap::sumz( pvecL, wvec )$p),2)
		pval_harmonic = min(c(pvalG,pvalL),na.rm=T)*2
		plot = dan.pairedDotPlot.multipages( factor(pairdf$State,levels=c( "Proliferative","Hypoxic" )), pairdf$Value, xlab="", ylab="Median patient-wise expression", plotTitle = paste0(xx,", combined p = ",pval_harmonic), labelLines = pairdf$Patient, signifTest = NULL, comparisons = NULL, ylimLeft = min(pairdf$Value)-0.01, ylimRight = max(pairdf$Value)+0.01, labelycoo = 1, jitterColors = pairdf$Colors, labelPoints = "", linesPairings = pairdf$Patient )
		print(plot)
	}
	dev.off()
}

findThinBoundary = function(qannot, thresh = 120, whichSide = "tumor"){
  d = dist(qannot[,c("um_x","um_y" ) ], method = "euclidean")
  d1 = as.matrix(d)
  d1 = d1[qannot[qannot$tumor.normal.boundary=="Boundary","Barcode"],qannot[qannot$tumor.normal.boundary==whichSide,"Barcode"]]
  d1 = d1[rowSums(d1<thresh)>0,]
  qannot$TNB = NA
  qannot[rownames(d1),"TNB"] = "Boundary"
  cat("\n", "Number of boundary spots:",length(rownames(d1)),"\n" )
  return(qannot)
}

findSector = function( qannot, minDistFromBoundary_um = 0, maxDistFromBoundary_um = 500, name = "sector_0_0.5mm", column = "TNB", whichSide = "tumor" ){
  d = dist(qannot[,c("um_x","um_y" ) ], method = "euclidean")
  d1 = as.matrix(d)
  d1 = d1[qannot[(qannot[,column]=="Boundary") %in% c(T),"Barcode"],qannot[qannot$tumor.normal.boundary==whichSide,"Barcode"]]
  if (length(intersect(colnames(d1),qannot[is.na(qannot[,column]),"Barcode"] ))>1) {
    d1 = d1[,intersect(colnames(d1),qannot[is.na(qannot[,column]),"Barcode"] )]
    d1 = d1[,(colSums(d1>minDistFromBoundary_um & d1<maxDistFromBoundary_um)>0)]
    qannot[colnames(d1),"TNB"] = name
    cat("\n", "Number of sector spots:",length(colnames(d1)),"\n" )
  }
  return(qannot)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

visium_boundary_analysis_tumor = function( OutDir,spacing=500 ){
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	SrDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/"
	DecDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Deconvolution/"
	this_OutDir = paste0(OutDir,"../AmsCentered_canceronlygenesFALSE_ThreeStates/")
	load(file = paste0(this_OutDir,"tps_scores.RData"))
	cs_scores = scores
	load(file = paste0(OutDir,"../tps_AmsCentered_canceronlygenesFALSE/","tps_scores.RData"))
	tps_scores = scores
	load(paste0(OutDir,"../../tps_discovery/tps.RData"))
	names_cs = c( "Alveolar","Proliferative","Hypoxic" )
	scores = cbind(tps_scores,cs_scores[rownames(tps_scores),names_cs])
	first = T
	minn_score = -Inf
	maxx_score = Inf
	for (Sample in batch1_table$Sample){
		### Plotting boundaries and sectors
		filez = paste0(SrDir,Sample,"/outs/tumor-normal-boundary.csv")
		if (!file.exists(filez)) { next }
		dcat(Sample)
		load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/",Sample,"_subspa_norm.RData"))
		qannot = read.csv(file = filez)
		qannot$Barcode = paste0(Sample,"_",qannot$Barcode)
		rownames(qannot) = qannot$Barcode
		qannot = qannot[qannot$tumor.normal.boundary %in% c("Boundary","tumor"),]
		qannot = qannot[intersect(rownames(qannot),colnames(subspa)),]
		tissue = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/NucleiDetection_StarDist/NucleiPerSpot/",Sample,"_NucleiPerSpot.txt") )
		rownames(tissue) = paste0(Sample,"_",tissue$barcode)
		qannot = qannot[intersect(rownames(qannot),rownames(tissue)),]
		tissue = tissue[rownames(qannot),]
		qannot$um_x = tissue$um_x
		qannot$um_y = tissue$um_y
		# next functions to be called in this order (not great but ok)
		subspa@meta.data$TNB = NA
		qannot = findThinBoundary(qannot)
		if (spacing==500){
			qannot = findSector(qannot,0,500,"0-0.5mm") 
			qannot = findSector(qannot,500,1000,"0.5-1mm")
			qannot = findSector(qannot,1000,1500,"1-1.5mm")
			qannot = findSector(qannot,1500,2000,"1.5-2mm")
			qannot = findSector(qannot,2000,2500,"2-2.5mm")
			tnb_levelz = c( "Boundary","0-0.5mm","0.5-1mm","1-1.5mm","1.5-2mm","2-2.5mm")
		}
		if (spacing==250){
			qannot = findSector(qannot,0,250,"sector_0_0.25mm") 
			qannot = findSector(qannot,250,500,"sector_0.25_0.5mm")
			qannot = findSector(qannot,500,750,"sector_0.5_0.75mm")
			qannot = findSector(qannot,750,1000,"sector_0.75_1mm")
			qannot = findSector(qannot,1000,1250,"sector_1_1.25mm")
			tnb_levelz = c( "Boundary","sector_0_0.25mm","sector_0.25_0.5mm","sector_0.5_0.75mm","sector_0.75_1mm","sector_1_1.25mm")
		}
		qannot = qannot[!is.na(qannot$TNB),]
		subspa = subset(subspa,cells=rownames(qannot))
		subspa@meta.data[rownames(qannot),"TNB"] = qannot[,"TNB"]
		subspa@meta.data$TNB = factor(subspa@meta.data$TNB,levels=tnb_levelz)
		tnb_colorz = c( "black","firebrick4","tomato3","orange","gray66","gray33" )
		names(tnb_colorz) = tnb_levelz
		pdf(file = paste0(OutDir,"SpatialPlot_TNB_",Sample,"_spacing",spacing,".pdf"),7,7)
		print(SpatialPlot(subspa, group.by = "TNB",cols=tnb_colorz,stroke=0,crop=F,pt.size.factor=1.3))
		dev.off()
		### Transcriptional programs and cell states across sectors
		df = qannot[!is.na(qannot$TNB),]
		df = df[rownames(df) %in% rownames(scores),]
		df = cbind(df[,c("Barcode","TNB")],scores[rownames(df),c(names_cs) ])
		df$Barcode = NULL
		mee = melt(df)
		mee$variable = gsub( "Alveolar","Alveolar-like",mee$variable )
		mee$variable = gsub( "Proliferative","Early-DD",mee$variable )
		mee$variable = gsub( "Hypoxic","Advanced-DD",mee$variable )
		mee$TNB = factor(mee$TNB,levels=tnb_levelz)
		mee$variable = factor(mee$variable,levels = c( "Alveolar-like","Early-DD","Advanced-DD" ))
		minn_score = min(c(minn_score,min(mee$value)))
		maxx_score = max(c(maxx_score,max(mee$value)))
		fileName = paste0(OutDir,Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		# Overall
		pdf(fileName,2,0.74)
		plot = ggplot(mee, aes(x=variable, y=value, fill=TNB)) + geom_boxplot(lwd=0.1,outlier.shape = NA) + scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("State score") + theme_classic() + theme(text = element_text(size=16))#, axis.text.x = element_text(angle = 45, hjust = 1))
		plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt')) + labs(x=NULL) +theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),legend.position="none",axis.line=element_line(size=0.1),axis.ticks = element_line(size = 0.1))
		print(plot)
		dev.off()
		mee$Sample = Sample
		dan.save(mee,fileName)
		### Thetas across sectors
		df = qannot[!is.na(qannot$TNB),]
		load(file = paste0(DecDir,"RCTD_Feb2024/",Sample,"_thetadf.RData"))
		rownames(thetadf) = paste0( Sample,"_",rownames(thetadf) )
		df = df[rownames(df) %in% rownames(thetadf),]
		df = cbind(df[,c("Barcode","TNB")],thetadf[rownames(df),])
		df$Barcode = NULL
		tee = melt(df)
		fileName = paste0(OutDir,"thetas_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		# Overall
		pdf(fileName,24,6)
		plot = ggplot(tee, aes(x=variable, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
		tee$Sample = Sample
		dan.save(tee,fileName)
		### Thetas across sectors
		df = qannot[!is.na(qannot$TNB),]
		load(file = paste0(DecDir,"RCTD_Feb2024_fine_reference/",Sample,"_thetadf.RData"))
		rownames(thetadf) = paste0( Sample,"_",rownames(thetadf) )
		df = df[rownames(df) %in% rownames(thetadf),]
		df = cbind(df[,c("Barcode","TNB")],thetadf[rownames(df),])
		df$Barcode = NULL
		tee3 = melt(df)
		fileName = paste0(OutDir,"thetas_",Sample,"_BoundaryAnalysis_level3_spacing",spacing,".pdf")
		# Overall
		pdf(fileName,28,6)
		plot = ggplot(tee3, aes(x=variable, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
		tee3$Sample = Sample
		dan.save(tee3,fileName)

		if (first){
		    meeAll = mee
		    teeAll = tee
		    teeAll3 = tee3
		    first = F
		  } else {
		    meeAll = rbind(meeAll,mee)
		    teeAll = rbind(teeAll,tee)
		    teeAll3 = rbind(teeAll3,tee3)
		  }
	}
	meeAll$Sample = substr(meeAll$Sample,8,nchar(meeAll$Sample))
	meeAll$TNB = factor(meeAll$TNB,levels=tnb_levelz)
	for (cs in c(names_cs)){
		fileName = paste0(OutDir,"AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = meeAll[meeAll$variable==cs,]
		pdf(fileName,8,2)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(size=0.3,outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab(paste0(cs," score")) + theme_classic(base_size=6) + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top")+ guides(fill = guide_legend(nrow = 1))
		plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		print(plot)
		dev.off()
	}
	for (cs in unique(teeAll$variable)){
		fileName = paste0(OutDir,"thetas_level2_AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = teeAll[teeAll$variable==cs,]
		pdf(fileName,20,6)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
	}
	for (cs in unique(teeAll3$variable)){
		fileName = paste0(OutDir,"thetas_level3_AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = teeAll3[teeAll3$variable==cs,]
		pdf(fileName,24,6)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
	}
	# selected
	ms = meeAll
	meeAll = meeAll[ meeAll$Sample %in% c("3D-Z1","3D-Z2","1E-Z2","1F-Z2"), ]
	meeAll$Sample = factor(meeAll$Sample,levels = c("1F-Z2","1E-Z2","3D-Z1","3D-Z2"))
	for (cs in c(names_cs)){
		fileName = paste0(OutDir,"selected_AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = meeAll[meeAll$variable==cs,]
		pdf(fileName,4,1.5)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(size=0.3,outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab(paste0(cs, " score")) + theme_classic(base_size=6) + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		print(plot)
		dev.off()
	}
}

visium_boundary_analysis_normal = function( OutDir,spacing=500 ){
	batch1_table = read.csv("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/batch1_visium_table.csv")
	SrDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/spaceranger/batch1/"
	DecDir = "/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Deconvolution/"
	# load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/","subspa_all.RData"))
	load(file = paste0( OutDir,"../../tps_atlases_NormalCells/gs_TravHlcaLPHan.RData" ))
	gs = gs[c( "HLCA_AT1","HLCA_AT2","Han_KAC_signature" )]
	# db_rand = dan.barkley_MakeRand_data(subspa_all,gs, 3)
	# save(db_rand,file = paste0(OutDir,"../visium_normalSignatures_db_rand_unintegrated.RData"))
	# scores = dan.Barkley_GeneToEnrichment_AmsCentered( subspa_all, gs, db_rand)
	# save(scores,file = paste0(OutDir,"../visium_normalSignatures_scores.RData"))
	load(file = paste0(OutDir,"../visium_normalSignatures_scores.RData"))
	first = T
	for (Sample in batch1_table$Sample){
		### Plotting boundaries and sectors
		if (Sample %in% c( "R19134-1C-Z1","R19134-1D-Z3" )){ next } # no good normal lung tissue
		filez = paste0(SrDir,Sample,"/outs/tumor-normal-boundary.csv")
		normal_spotz = paste0(SrDir,Sample,"/outs/normal_spotz.csv")
		if (!file.exists(filez)) { next }
		dcat(Sample)
		load(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/Visium_exploration/",Sample,"_subspa_norm.RData"))
		qannot = read.csv(file = filez)
		if (file.exists(normal_spotz)) { 
			normal_spotz = read.csv(file = normal_spotz)
			qannot = qannot[qannot$Barcode %in% normal_spotz[normal_spotz$normal_spotz=="normal_spotz","Barcode"],]
		}
		qannot$Barcode = paste0(Sample,"_",qannot$Barcode)
		rownames(qannot) = qannot$Barcode
		qannot = qannot[qannot$tumor.normal.boundary %in% c("Boundary",""),]
		qannot = qannot[intersect(rownames(qannot),colnames(subspa)),]
		tissue = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/spatial/Processed/NucleiDetection_StarDist/NucleiPerSpot/",Sample,"_NucleiPerSpot.txt") )
		rownames(tissue) = paste0(Sample,"_",tissue$barcode)
		qannot = qannot[intersect(rownames(qannot),rownames(tissue)),]
		tissue = tissue[rownames(qannot),]
		qannot$um_x = tissue$um_x
		qannot$um_y = tissue$um_y
		# next functions to be called in this order (not great but ok)
		subspa@meta.data$TNB = NA
		qannot = findThinBoundary(qannot,whichSide = "")
		if (spacing==500){
			qannot = findSector(qannot,0,500,"sector_0_0.5mm",whichSide = "") 
			qannot = findSector(qannot,500,1000,"sector_0.5_1mm",whichSide = "") 
			qannot = findSector(qannot,1000,1500,"sector_1_1.5mm",whichSide = "") 
			qannot = findSector(qannot,1500,2000,"sector_1.5_2mm",whichSide = "") 
			qannot = findSector(qannot,2000,2500,"sector_2_2.5mm",whichSide = "") 
			tnb_levelz = c( "Boundary","sector_0_0.5mm","sector_0.5_1mm","sector_1_1.5mm","sector_1.5_2mm","sector_2_2.5mm")
		}
		if (spacing==250){
			qannot = findSector(qannot,0,250,"sector_0_0.25mm",whichSide = "") 
			qannot = findSector(qannot,250,500,"sector_0.25_0.5mm",whichSide = "") 
			qannot = findSector(qannot,500,750,"sector_0.5_0.75mm",whichSide = "") 
			qannot = findSector(qannot,750,1000,"sector_0.75_1mm",whichSide = "") 
			qannot = findSector(qannot,1000,1250,"sector_1_1.25mm",whichSide = "") 
			tnb_levelz = c( "Boundary","sector_0_0.25mm","sector_0.25_0.5mm","sector_0.5_0.75mm","sector_0.75_1mm","sector_1_1.25mm")
		}
		qannot = qannot[!is.na(qannot$TNB),]
		subspa = subset(subspa,cells=rownames(qannot))
		subspa@meta.data[rownames(qannot),"TNB"] = qannot[,"TNB"]
		subspa@meta.data$TNB = factor(subspa@meta.data$TNB,levels=tnb_levelz)
		tnb_colorz = c( "black","firebrick4","tomato3","orange","gray66","gray33" )
		names(tnb_colorz) = tnb_levelz
		pdf(file = paste0(OutDir,"SpatialPlot_TNB_",Sample,"_spacing",spacing,".pdf"),7,7)
		print(SpatialPlot(subspa, group.by = "TNB",cols=tnb_colorz,crop=F,pt.size.factor=1.3))
		dev.off()
		### Transcriptional programs and cell states across sectors
		df = qannot[!is.na(qannot$TNB),]
		df = df[rownames(df) %in% rownames(scores),]
		df = cbind(df[,c("Barcode","TNB")],scores[rownames(df),names(gs) ])
		df$Barcode = NULL
		mee = melt(df)
		fileName = paste0(OutDir,Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		# Overall
		pdf(fileName,12,6)
		plot = ggplot(mee, aes(x=variable, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
		mee$Sample = Sample
		dan.save(mee,fileName)
		### Thetas across sectors
		df = qannot[!is.na(qannot$TNB),]
		load(file = paste0(DecDir,"RCTD_Feb2024/",Sample,"_thetadf.RData"))
		rownames(thetadf) = paste0( Sample,"_",rownames(thetadf) )
		df = df[rownames(df) %in% rownames(thetadf),]
		df = cbind(df[,c("Barcode","TNB")],thetadf[rownames(df),])
		df$Barcode = NULL
		tee = melt(df)
		fileName = paste0(OutDir,"thetas_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		# Overall
		pdf(fileName,24,6)
		plot = ggplot(tee, aes(x=variable, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
		tee$Sample = Sample
		dan.save(tee,fileName)

		df$ratio_AT2_AT1 = log10(df$AT2/df$AT1)
		fileName = paste0(OutDir,"ratio_AT2_AT1_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		x = df$TNB
		y = df$ratio_AT2_AT1
		jitterColorz = dan.expand_colors(df$TNB,sort(unique(df$TNB)),tnb_colorz[sort(unique(df$TNB))])
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		dan.boxplots( fileName, x, y,xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio = 1.5",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(df$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 6, fileHeight = 5, hlines_coo = log10(1.5) )

		df$ratio_AT0 = (df$AT0/(df$AT0+df$AT1+df$AT2+df$basal+df$ciliated+df$club+df$preTB) )
		fileName = paste0(OutDir,"ratio_AT0_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		x = df$TNB
		y = df$ratio_AT0
		jitterColorz = dan.expand_colors(df$TNB,sort(unique(df$TNB)),tnb_colorz[sort(unique(df$TNB))])
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		dan.boxplots( fileName, x, y,xlab = "", ylab = "# AT0s / # all epithelial cells",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = 0, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(df$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 6, fileHeight = 5)

		tee = melt(df)
		tee$Sample = Sample
		dan.save(tee,paste0(OutDir,"thetas_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf"))

		df = df[rownames(df) %in% rownames(scores),]
		df = cbind(df,scores[rownames(df),names(gs) ])
		df$HLCA_AT2_normalized = df$HLCA_AT2/df$AT2
		fileName = paste0(OutDir,"HLCA_AT2_normalized_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		x = df$TNB
		y = df$HLCA_AT2_normalized
		jitterColorz = dan.expand_colors(df$TNB,sort(unique(df$TNB)),tnb_colorz[sort(unique(df$TNB))])
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		dan.boxplots( fileName, x, y,xlab = "", ylab = "HLCA_AT2 normalized by AT2 proportion",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(df$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 6, fileHeight = 5 )

		df = df[rownames(df) %in% rownames(scores),]
		df = cbind(df,scores[rownames(df),names(gs) ])
		df$HLCA_AT1_normalized = df$HLCA_AT1/df$AT1
		fileName = paste0(OutDir,"HLCA_AT1_normalized_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		x = df$TNB
		y = df$HLCA_AT1_normalized
		jitterColorz = dan.expand_colors(df$TNB,sort(unique(df$TNB)),tnb_colorz[sort(unique(df$TNB))])
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		dan.boxplots( fileName, x, y,xlab = "", ylab = "HLCA_AT1 normalized by AT1 proportion",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(df$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 6, fileHeight = 5 )

		df = df[rownames(df) %in% rownames(scores),]
		df = cbind(df,scores[rownames(df),names(gs) ])
		fileName = paste0(OutDir,"Han_KAC_signature_",Sample,"_BoundaryAnalysis_spacing",spacing,".pdf")
		x = df$TNB
		y = df$Han_KAC_signature
		jitterColorz = dan.expand_colors(df$TNB,sort(unique(df$TNB)),tnb_colorz[sort(unique(df$TNB))])
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		dan.boxplots( fileName, x, y,xlab = "", ylab = "Han_KAC_signature",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(df$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 6, fileHeight = 5 )

		nee = melt(df[,c( "TNB","HLCA_AT1_normalized","HLCA_AT2_normalized","Han_KAC_signature" )])
		nee$Sample = Sample


		### Thetas across sectors
		df = qannot[!is.na(qannot$TNB),]
		load(file = paste0(DecDir,"RCTD_Feb2024_fine_reference/",Sample,"_thetadf.RData"))
		rownames(thetadf) = paste0( Sample,"_",rownames(thetadf) )
		df = df[rownames(df) %in% rownames(thetadf),]
		df = cbind(df[,c("Barcode","TNB")],thetadf[rownames(df),])
		df$Barcode = NULL
		tee3 = melt(df)
		fileName = paste0(OutDir,"thetas_",Sample,"_BoundaryAnalysis_level3_spacing",spacing,".pdf")
		# Overall
		pdf(fileName,28,6)
		plot = ggplot(tee3, aes(x=variable, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
		tee3$Sample = Sample
		dan.save(tee3,fileName)

		if (first){
		    meeAll = mee
		    teeAll = tee
		    teeAll3 = tee3
		    neeAll = nee
		    first = F
		  } else {
		    meeAll = rbind(meeAll,mee)
		    teeAll = rbind(teeAll,tee)
		    teeAll3 = rbind(teeAll3,tee3)
		    neeAll = rbind(neeAll,nee)
		  }
	}
	for (cs in names(gs)){
		fileName = paste0(OutDir,"AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = meeAll[meeAll$variable==cs,]
		pdf(fileName,22,6)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
	}
	for (cs in unique(teeAll$variable)){
		fileName = paste0(OutDir,"thetas_level2_AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = teeAll[teeAll$variable==cs,]
		pdf(fileName,20,6)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
	}
	for (cs in unique(teeAll3$variable)){
		fileName = paste0(OutDir,"thetas_level3_AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = teeAll3[teeAll3$variable==cs,]
		pdf(fileName,24,6)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
	}
	for (cs in unique(neeAll$variable)){
		fileName = paste0(OutDir,"norm_AcrossSamples_",cs,"_BoundaryAnalysis_spacing",spacing,".pdf")
		this_mee = neeAll[neeAll$variable==cs,]
		pdf(fileName,20,6)
		plot = ggplot(this_mee, aes(x=Sample, y=value, fill=TNB)) + geom_boxplot(outlier.shape = NA)+ scale_fill_manual(values = tnb_colorz) + xlab("") + ylab("Score") + theme_classic() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
		print(plot)
		dev.off()
	}
	
	thisnee = neeAll[neeAll$variable=="HLCA_AT1_normalized",]
	fileName = paste0(OutDir,"norm_AllSamples_AT1_BoundaryAnalysis_spacing",spacing,".pdf")
	x = thisnee$TNB
	y = thisnee$value
	jitterColorz = dan.expand_colors(thisnee$TNB,sort(unique(thisnee$TNB)),tnb_colorz[sort(unique(thisnee$TNB))])
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	dan.boxplots( fileName, x, y,xlab = "", ylab = "HLCA_AT1 normalized by AT1 proportion",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(thisnee$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 7, fileHeight = 5 )

	thisnee = neeAll[neeAll$variable=="HLCA_AT2_normalized",]
	fileName = paste0(OutDir,"norm_AllSamples_AT2_BoundaryAnalysis_spacing",spacing,".pdf")
	x = thisnee$TNB
	y = thisnee$value
	jitterColorz = dan.expand_colors(thisnee$TNB,sort(unique(thisnee$TNB)),tnb_colorz[sort(unique(thisnee$TNB))])
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	dan.boxplots( fileName, x, y,xlab = "", ylab = "HLCA_AT2 normalized by AT2 proportion",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(thisnee$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 7, fileHeight = 5 )

	thisnee = neeAll[neeAll$variable=="Han_KAC_signature",]
	fileName = paste0(OutDir,"norm_AllSamples_KAC_BoundaryAnalysis_spacing",spacing,".pdf")
	x = thisnee$TNB
	y = thisnee$value
	jitterColorz = dan.expand_colors(thisnee$TNB,sort(unique(thisnee$TNB)),tnb_colorz[sort(unique(thisnee$TNB))])
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	dan.boxplots( fileName, x, y,xlab = "", ylab = "KAC signature",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(thisnee$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 7, fileHeight = 5 )

	thisnee = teeAll[teeAll$variable=="ratio_AT2_AT1",]
	fileName = paste0(OutDir,"norm_AllSamples_ratio_AT2_AT1_BoundaryAnalysis_spacing",spacing,".pdf")
	x = thisnee$TNB
	y = thisnee$value
	jitterColorz = dan.expand_colors(thisnee$TNB,sort(unique(thisnee$TNB)),tnb_colorz[sort(unique(thisnee$TNB))])
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	dan.boxplots( fileName, x, y,xlab = "", ylab = "log10(AT2/AT1)", plotTitle = "Healthy AT2/AT1 ratio = 1.5",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = ylimLeft, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(thisnee$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 7, fileHeight = 5,hlines_coo = log10(1.5) )

	thisnee = teeAll[teeAll$variable=="ratio_AT0",]
	fileName = paste0(OutDir,"norm_AllSamples_ratio_AT0_BoundaryAnalysis_spacing",spacing,".pdf")
	x = thisnee$TNB
	y = thisnee$value
	jitterColorz = dan.expand_colors(thisnee$TNB,sort(unique(thisnee$TNB)),tnb_colorz[sort(unique(thisnee$TNB))])
	ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
   	ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	dan.boxplots( fileName, x, y,xlab = "", ylab = "# AT0s / # all epithelial cells",signifTest = "kruskal",labelycoo = ylimRight, ylimLeft = 0, ylimRight = ylimRight,xColors = tnb_colorz[sort(unique(thisnee$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 7, fileHeight = 5)

	for (pat in unique(substr(teeAll$Sample,1,6))){
		thist = teeAll[substr(teeAll$Sample,1,6)==pat,]
		tdf2 = dan.df(unique(thist$TNB),c("TNB","AT2props","AT1props"))
		for (rn in rownames(tdf2)){
			tdf2[rn,"AT2props"] = sum(thist[(thist$TNB==rn) & (thist$variable=="AT2"),"value"])/(sum(thist[(thist$TNB==rn) & (thist$variable=="AT2"),"value"])+sum(thist[(thist$TNB==rn) & (thist$variable=="AT1"),"value"]))
			tdf2[rn,"AT1props"] = sum(thist[(thist$TNB==rn) & (thist$variable=="AT1"),"value"])/(sum(thist[(thist$TNB==rn) & (thist$variable=="AT2"),"value"])+sum(thist[(thist$TNB==rn) & (thist$variable=="AT1"),"value"]))
			tdf2[rn,"TNB"] = rn
			
		}
		tdf2 = melt(tdf2[,c("TNB","AT2props","AT1props" )])
		tdf2$variable = factor(tdf2$variable,levels=c( "AT1props","AT2props" ))
		tdf2$TNB = factor(tdf2$TNB, levels = tnb_levelz[tnb_levelz %in% unique(tdf2$TNB)])
		px = ggplot(data=tdf2, aes(x=TNB, y=value, fill=variable)) + geom_bar(stat="identity", colour="black", linewidth = 0.1) + labs(x = "Distance from tumor-normal boundary", y = "Proportion" ) + ggtitle("") + scale_fill_manual(name = "Alveolar\ncell type",values=c( "chocolate","dodgerblue2" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))
		pdf( paste0(OutDir,"ratio_AT2_AT1_vsDistance_tileSizeMicrons",spacing,"_BarplotProportions_",pat,".pdf"),5,3 )
		print(px)
		dev.off()
	}

	thisnee = neeAll[neeAll$variable=="Han_KAC_signature",]
	for (pat in unique(substr(thisnee$Sample,1,6))){
		fileName = paste0(OutDir,"norm_PatientWise_KAC_BoundaryAnalysis_spacing",spacing,"_",pat,".pdf")
		thist = thisnee[substr(thisnee$Sample,1,6)==pat,]
		x = factor(thist$TNB,levels = tnb_levelz[tnb_levelz %in% unique(tdf2$TNB)])
		y = thist$value
		jitterColorz = dan.expand_colors(thist$TNB,sort(unique(thisnee$TNB)),tnb_colorz[sort(unique(thisnee$TNB))])
		ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
   		ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
		dan.boxplots( fileName, x, y,xlab = "", ylab = "KAC signature",signifTest = "kruskal",labelycoo = ylimRight, xColors = tnb_colorz[sort(unique(thisnee$TNB))], jitterColors = jitterColorz, jitterDotSize = 1, fileWidth = 6, fileHeight = 5 ) # ylimLeft = ylimLeft, ylimRight = ylimRight,
	}
}

pvals_to_signLevels = function( cordf_pvals ){
	mat = cordf_pvals
	mat[cordf_pvals>0.01] = ""
	mat[cordf_pvals<0.01] = "*"
	mat[cordf_pvals<0.001] = "**"
	mat[cordf_pvals<0.0001] = "***"
	return( mat )
}

TME_formatted_heatmaps = function( OutDir,cs_map ){
	## quite hardcoded, careful
	library(colorRamp2)
	library(ComplexHeatmap)
	colorz = colorRamp2(c(-1,-0.5,0,0.5,1),c("dodgerblue4","dodgerblue3","white","firebrick3","firebrick4"))

	if (nrow(cs_map)==3){
		## Level 1
		load( paste0(OutDir,"../II_vsTotalCancer_TotalTME/h12_clr_csRatios_vs_level1_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../II_vsTotalCancer_TotalTME/h12_clr_csRatios_vs_level1_table.RData") )
		cordf = cordf[,c( "n_level1_epithelial","n_level1_endothelial","n_level1_immune","n_level1_stroma")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Epithelial","Endothelial","Immune","Stromal" )
		cordf = cordf[, c( "Epithelial","Endothelial","Immune","Stromal" ) ]
		cordf_pvals = cordf_pvals[,c( "n_level1_epithelial","n_level1_endothelial","n_level1_immune","n_level1_stroma")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Epithelial","Endothelial","Immune","Stromal" )
		cordf_pvals = cordf_pvals[, c( "Epithelial","Endothelial","Immune","Stromal" ) ]
		pdf(paste0(OutDir,"h12_ComplexHeatmap_level1.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 2, immune
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table.RData") )
		cordf = cordf[,c( "n_level2_mast","n_level2_DC","n_level2_Bcells","n_level2_NKcells","n_level2_monocytes","n_level2_Tcells","n_level2_macrophages","n_level2_neutrophils" )]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Mast cells","DC","B cells","NK cells","Monocytes","T cells","Macrophages","Neutrophils" )
		cordf = cordf[, c( "DC","T cells","Mast cells", "B cells", "Neutrophils", "NK cells","Macrophages","Monocytes" ) ]
		cordf_pvals = cordf_pvals[,c( "n_level2_mast","n_level2_DC","n_level2_Bcells","n_level2_NKcells","n_level2_monocytes","n_level2_Tcells","n_level2_macrophages","n_level2_neutrophils" )]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Mast cells","DC","B cells","NK cells","Monocytes","T cells","Macrophages","Neutrophils" )
		cordf_pvals = cordf_pvals[, c( "DC","T cells","Mast cells", "B cells", "Neutrophils", "NK cells","Macrophages","Monocytes" ) ]
		pdf(paste0(OutDir,"h12_ComplexHeatmap_level2_immune.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 2, endothelial
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table.RData") )
		cordf = cordf[,c( "n_level2_EC_blood","n_level2_EC_lymphatic")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "EC blood", "EC lymphatic" )
		cordf = cordf[, c( "EC blood", "EC lymphatic" ) ]
		cordf_pvals = cordf_pvals[,c( "n_level2_EC_blood","n_level2_EC_lymphatic")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "EC blood", "EC lymphatic" )
		cordf_pvals = cordf_pvals[, c( "EC blood", "EC lymphatic" ) ]
		pdf(paste0(OutDir,"h12_ComplexHeatmap_level2_endothelial.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 2, epithelial
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table.RData") )
		cordf = cordf[,c( "n_level2_AT0","n_level2_AT1","n_level2_AT2","n_level2_basal","n_level2_ciliated","n_level2_club","n_level2_preTB")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "AT0","AT1","AT2","Basal","Ciliated","Club","PreTB" )
		cordf = cordf[, c( "AT0","PreTB","AT2","Ciliated","Club","Basal","AT1" ) ]
		cordf_pvals = cordf_pvals[,c( "n_level2_AT0","n_level2_AT1","n_level2_AT2","n_level2_basal","n_level2_ciliated","n_level2_club","n_level2_preTB")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "AT0","AT1","AT2","Basal","Ciliated","Club","PreTB" )
		cordf_pvals = cordf_pvals[, c( "AT0","PreTB","AT2","Ciliated","Club","Basal","AT1" ) ]
		pdf(paste0(OutDir,"h12_ComplexHeatmap_level2_epithelial.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 2, stromal
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_level2_table.RData") )
		cordf = cordf[,c( "n_level2_fibroblasts","n_level2_mesothelial","n_level2_smooth_muscle")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Fibroblasts","Mesothelial","Smooth muscle" )
		cordf = cordf[, c( "Fibroblasts","Mesothelial","Smooth muscle" ) ]
		cordf_pvals = cordf_pvals[,c( "n_level2_fibroblasts","n_level2_mesothelial","n_level2_smooth_muscle")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Fibroblasts","Mesothelial","Smooth muscle" )
		cordf_pvals = cordf_pvals[, c( "Fibroblasts","Mesothelial","Smooth muscle" ) ]
		pdf(paste0(OutDir,"h12_ComplexHeatmap_level2_stromal.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 3, B cells
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_GC_Bcells","n_level3_memory_Bcells","n_level3_naive_Bcells","n_level3_plasma_Bcells")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Germinal center","Memory","Naive","Plasma" )
		cordf = cordf[,c( "Memory","Plasma","Naive","Germinal center" )]
		cordf_pvals = cordf_pvals[,c( "n_level3_GC_Bcells","n_level3_memory_Bcells","n_level3_naive_Bcells","n_level3_plasma_Bcells")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Germinal center","Memory","Naive","Plasma" )
		cordf_pvals = cordf_pvals[, c( "Memory","Plasma","Naive","Germinal center" )]
		pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_Bcells.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 3, macrophages
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_alveolar_macro","n_level3_mono_derived_macro")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Alveolar","Monocyte-derived")
		cordf_pvals = cordf_pvals[,c( "n_level3_alveolar_macro","n_level3_mono_derived_macro")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Alveolar","Monocyte-derived")
		pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_Macrophages.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 3, monocytes
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_non_classical_mono","n_level3_classical_mono")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Non-classical","Classical")
		cordf_pvals = cordf_pvals[,c( "n_level3_non_classical_mono","n_level3_classical_mono")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Non-classical","Classical")
		pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_Monocytes.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 3, CD4 T cells
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_CD4.CTL_EOMES","n_level3_CD4.CTL_Exh","n_level3_CD4.CTL_GNLY","n_level3_CD4.NaiveLike","n_level3_CD4.Tfh","n_level3_CD4.Th17","n_level3_CD4.Treg")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "CTL EOMES+","CTL exhausted","CTL GNLY+","Naive-like","Tfh","Th17","Treg")
		cordf = cordf[,c( "Naive-like","CTL EOMES+","CTL GNLY+","Tfh","Treg","Th17","CTL exhausted")]
		cordf_pvals = cordf_pvals[,c( "n_level3_CD4.CTL_EOMES","n_level3_CD4.CTL_Exh","n_level3_CD4.CTL_GNLY","n_level3_CD4.NaiveLike","n_level3_CD4.Tfh","n_level3_CD4.Th17","n_level3_CD4.Treg")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "CTL EOMES+","CTL exhausted","CTL GNLY+","Naive-like","Tfh","Th17","Treg")
		cordf_pvals = cordf_pvals[, c( "Naive-like","CTL EOMES+","CTL GNLY+","Tfh","Treg","Th17","CTL exhausted")]
		pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_CD4T.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 3, CD8 T cells
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_CD8.CM","n_level3_CD8.EM","n_level3_CD8.MAIT","n_level3_CD8.NaiveLike","n_level3_CD8.TEMRA","n_level3_CD8.TEX","n_level3_CD8.TPEX")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Central memory","Effector memory","MAIT","Naive-like","TEMRA","Exhausted","Precursor-exhausted")
		cordf = cordf[, c( "Central memory","Effector memory","TEMRA","Naive-like","Precursor-exhausted","MAIT","Exhausted") ]
		cordf_pvals = cordf_pvals[,c( "n_level3_CD8.CM","n_level3_CD8.EM","n_level3_CD8.MAIT","n_level3_CD8.NaiveLike","n_level3_CD8.TEMRA","n_level3_CD8.TEX","n_level3_CD8.TPEX")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "Central memory","Effector memory","MAIT","Naive-like","TEMRA","Exhausted","Precursor-exhausted")
		cordf_pvals = cordf_pvals[, c( "Central memory","Effector memory","TEMRA","Naive-like","Precursor-exhausted","MAIT","Exhausted") ]
		pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_CD8T.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		## Level 3, DC
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table_pvals.RData") )
		cordf_pvals = pvals_to_signLevels( cordf_pvals )
		load( paste0(OutDir,"../IV_vsSubsets/h12_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_DC1","n_level3_DC2","n_level3_DC3","n_level3_pDC")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "DC1","DC2","DC3","pDC")
		cordf = cordf[, c( "DC2","DC1","DC3","pDC") ]
		cordf_pvals = cordf_pvals[,c( "n_level3_DC1","n_level3_DC2","n_level3_DC3","n_level3_pDC")]
		rownames(cordf_pvals) = cs_map[gsub("n_","",rownames(cordf_pvals)) ,"alias"]
		colnames(cordf_pvals) = c( "DC1","DC2","DC3","pDC")
		cordf_pvals = cordf_pvals[, c( "DC2","DC1","DC3","pDC") ]
		pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_DC.pdf"),4,2)
		plot=(ComplexHeatmap::Heatmap(cordf,show_heatmap_legend=F,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		draw(plot,heatmap_legend_side="right")
		dev.off()

		# ## macrophage flavours
		# load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/h12_clr_csRatios_vs_macro_flavours_table.RData") )
		# rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		# colnames(cordf) = substr(colnames(cordf),4,nchar(colnames(cordf)))
		# cordf = cordf[,c( "Alveolar_RTM_like","Reg_TAMs","Prolif_TAMs", "LA_TAMs","MT_RTM_like","IFN_TAMs","Inflam_TAMs","Angio_TAMs" )]
		# pdf(paste0(OutDir, "h12_ComplexHeatmap_level3_macro_flavours.pdf"),4,2)
		# plot=(ComplexHeatmap::Heatmap(cordf,rect_gp = gpar(col = "white", lwd = 1),width=unit(ncol(cordf), "cm")/2,height=unit(nrow(cordf), "cm")/2,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorRamp2(c(-1,0,1),c("dodgerblue4","white","firebrick3")), heatmap_legend_param=list(grid_height=unit(6, "pt"),grid_width=unit(6, "pt"),title_position="topcenter",title="CLR-transformed cell\nproportions, Spearman R",direction = "horizontal",title_gp=gpar(fontsize = 6),labels_gp = gpar(fontsize = 6)),column_names_side = "bottom",row_names_side = "left",column_names_gp = grid::gpar(fontsize = 6),row_names_gp = grid::gpar(fontsize = 6),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(cordf_pvals[i,j], x, y, gp = gpar(fontsize = 6))} ))
		# draw(plot,heatmap_legend_side="right")
		# dev.off()
	}
	if (nrow(cs_map)==2){
		## Level 2
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/twocs_clr_csRatios_vs_level2_table.RData") )
		cordf = cordf[,c( "n_level2_mast","n_level2_DC","n_level2_Bcells","n_level2_NKcells","n_level2_monocytes","n_level2_Tcells","n_level2_macrophages","n_level2_neutrophils" )]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Mast cells","DC","B cells","NK cells","Monocytes","T cells","Macrophages","Neutrophils" )
		cordf = cordf[, c( "DC","T cells","Mast cells", "B cells", "Neutrophils", "NK cells","Macrophages","Monocytes" ) ]
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level2.pdf"),8,2)
		print(ComplexHeatmap::Heatmap(cordf,cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="CLR-transformed cell\nproportion correlation",direction = "horizontal"),column_names_side = "top",row_names_side = "left"))
		dev.off()

		## Level 3, B cells
		load( paste0(OutDir,"../IV_vsSubsets/twocs_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_GC_Bcells","n_level3_memory_Bcells","n_level3_naive_Bcells","n_level3_plasma_Bcells")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Germinal center","Memory","Naive","Plasma" )
		cordf = cordf[,c( "Memory","Plasma","Naive","Germinal center" )]
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level3_Bcells.pdf"),ncol(cordf)*1.5,nrow(cordf)*1.4)
		print(ComplexHeatmap::Heatmap(cordf,column_title="B cells",cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="CLR-transformed cell\nproportion correlation",direction = "horizontal"),column_names_side = "bottom",row_names_side = "left"))
		dev.off()

		## Level 3, macrophages
		load( paste0(OutDir,"../IV_vsSubsets/twocs_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_alveolar_macro","n_level3_mono_derived_macro")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Alveolar","Monocyte-derived")
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level3_Macrophages.pdf"),ncol(cordf)*2,nrow(cordf)*1.3)
		print(ComplexHeatmap::Heatmap(cordf,column_title="Macrophages",cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="CLR-transformed cell\nproportion correlation",direction = "horizontal"),column_names_side = "bottom",row_names_side = "left"))
		dev.off()

		## Level 3, monocytes
		load( paste0(OutDir,"../IV_vsSubsets/twocs_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_non_classical_mono","n_level3_classical_mono")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Non-classical","Classical")
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level3_Monocytes.pdf"),ncol(cordf)*2,nrow(cordf)*1.2)
		print(ComplexHeatmap::Heatmap(cordf,column_title="Monocytes",cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="CLR-transformed cell\nproportion correlation",direction = "horizontal"),column_names_side = "bottom",row_names_side = "left"))
		dev.off()

		## Level 3, CD4 T cells
		load( paste0(OutDir,"../IV_vsSubsets/twocs_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_CD4.CTL_EOMES","n_level3_CD4.CTL_Exh","n_level3_CD4.CTL_GNLY","n_level3_CD4.NaiveLike","n_level3_CD4.Tfh","n_level3_CD4.Th17","n_level3_CD4.Treg")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "CTL EOMES+","CTL exhausted","CTL GNLY+","Naive-like","Tfh","Th17","Treg")
		cordf = cordf[,c( "Naive-like","CTL EOMES+","CTL GNLY+","Tfh","Treg","Th17","CTL exhausted")]
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level3_CD4T.pdf"),ncol(cordf),nrow(cordf)*1.2)
		print(ComplexHeatmap::Heatmap(cordf,column_title="CD4+ T cells",cluster_rows=FALSE, cluster_columns=F, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="CLR-transformed cell\nproportion correlation",direction = "horizontal"),column_names_side = "bottom",row_names_side = "left"))
		dev.off()

		## Level 3, CD8 T cells
		load( paste0(OutDir,"../IV_vsSubsets/twocs_clr_csRatios_vs_level3_table.RData") )
		cordf = cordf[,c( "n_level3_CD8.CM","n_level3_CD8.EM","n_level3_CD8.MAIT","n_level3_CD8.NaiveLike","n_level3_CD8.TEMRA","n_level3_CD8.TEX","n_level3_CD8.TPEX")]
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = c( "Central memory","Effector memory","MAIT","Naive-like","TEMRA","Exhausted","Precursor-exhausted")
		cordf = cordf[, c( "Central memory","Effector memory","TEMRA","Naive-like","Precursor-exhausted","MAIT","Exhausted") ]
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level3_CD8T.pdf"),ncol(cordf),nrow(cordf)*1.3)
		print(ComplexHeatmap::Heatmap(cordf,column_title="CD8+ T cells",cluster_rows=FALSE, cluster_columns=F, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="CLR-transformed cell\nproportion correlation",direction = "horizontal"),column_names_side = "bottom",row_names_side = "left"))
		dev.off()

		## macrophage flavours
		load( paste0(OutDir,"../III_vsTotalCancer_TotalTME_splittedbylevel1/twocs_clr_csRatios_vs_macro_flavours_table.RData") )
		rownames(cordf) = cs_map[gsub("n_","",rownames(cordf)) ,"alias"]
		colnames(cordf) = substr(colnames(cordf),4,nchar(colnames(cordf)))
		cordf = cordf[,c( "Alveolar_RTM_like","Reg_TAMs","Prolif_TAMs", "LA_TAMs","MT_RTM_like","IFN_TAMs","Inflam_TAMs","Angio_TAMs" )]
		pdf(paste0(OutDir, "twocs_ComplexHeatmap_level3_macro_flavours.pdf"),ncol(cordf),nrow(cordf)*1.3)
		print(ComplexHeatmap::Heatmap(cordf,column_title="Macrophage subsets",cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE, col=colorz, heatmap_legend_param=list(title="Spearman correlation",direction = "horizontal"),column_names_side = "bottom",row_names_side = "left"))
		dev.off()
	}
}

xenium_preprocessing = function( OutDir, whichDataset ){ # coordinates should be in um
	xdir = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/",whichDataset,"/")
	if (whichDataset %in% c("xenium_18May2023")){ xe = LoadXenium(xdir, fov = "fov") }
	if (whichDataset %in% c("xenium_05February2024","xenium_haga2023_TSU21","xenium_haga2023_TSU20","xenium_nondiseasedlung") ){ 
		xe = ReadXenium(data.dir = xdir, type = c("centroids", "segmentations"))
		xe %>% names()
		xe %>% str()
		assay = "Xenium"
		segmentations.data = list("centroids" = CreateCentroids(xe$centroids),"segmentation" = CreateSegmentation(xe$segmentations))
		coords = CreateFOV(coords = segmentations.data,type = c("segmentation", "centroids"),molecules = xe$microns,assay = assay)
		xeobj = CreateSeuratObject(counts = xe$matrix[["Gene Expression"]], assay = assay)
		xeobj[["BlankCodeword"]] = CreateAssayObject(counts = xe$matrix[["Unassigned Codeword"]])
		xeobj[["ControlCodeword"]] = CreateAssayObject(counts = xe$matrix[["Negative Control Codeword"]])
		xeobj[["ControlProbe"]] = CreateAssayObject(counts = xe$matrix[["Negative Control Probe"]])
		xeobj[['fov']] = coords
		xe = xeobj
	}
	if (whichDataset %in% c("xenium_takano2024_luad2","xenium_takano2024_luad3","xenium_takano2024_luad14","xenium_takano2024_luad16","xenium_takano2024_luad17") ){ 
		nn = gsub("xenium_takano2024_luad","",whichDataset)
		xdir = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/xenium_takano2024/Xenium_LUAD_No",nn,"/")
		xe = ReadXenium(data.dir = xdir, type = c("centroids", "segmentations"))
		xe %>% names()
		xe %>% str()
		assay = "Xenium"
		segmentations.data = list("centroids" = CreateCentroids(xe$centroids),"segmentation" = CreateSegmentation(xe$segmentations))
		coords = CreateFOV(coords = segmentations.data,type = c("segmentation", "centroids"),molecules = xe$microns,assay = assay)
		xeobj = CreateSeuratObject(counts = xe$matrix[["Gene Expression"]], assay = assay)
		xeobj[["BlankCodeword"]] = CreateAssayObject(counts = xe$matrix[["Unassigned Codeword"]])
		xeobj[["ControlCodeword"]] = CreateAssayObject(counts = xe$matrix[["Negative Control Codeword"]])
		xeobj[["ControlProbe"]] = CreateAssayObject(counts = xe$matrix[["Negative Control Probe"]])
		xeobj[['fov']] = coords
		xe = xeobj
	}
	md = xe@meta.data
	xe = subset(xe, subset = nCount_Xenium > 10)
	pdf(paste0(OutDir,"qc.pdf" ),7,7)
	plot=VlnPlot(xe, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
	print(plot)
	dev.off()
	xe <- SCTransform(xe, assay = "Xenium")
	xe <- RunPCA(xe, npcs = 30, features = rownames(xe))
	xe <- RunUMAP(xe, dims = 1:30)
	xe <- FindNeighbors(xe, reduction = "pca", dims = 1:30)
	xe <- FindClusters(xe,res=1.5)
	pdf(paste0(OutDir,"dimplot_clusters.pdf" ),12,8)
	print(DimPlot(xe,label=T) + coord_fixed())
	dev.off()
	pdf(paste0(OutDir,"featureplot_nCounts_nFeatures.pdf" ),16,8)
	print(FeaturePlot(xe, features = c( "nCount_Xenium","nFeature_Xenium" ),max.cutoff=c('q99','q99' ) ) + coord_fixed())
	dev.off()

	am = FindAllMarkers( xe, only.pos=TRUE )
	am = am[am$p_val_adj<0.01,]
	save(am,file=paste0(OutDir,"AllMarkers.RData" ))

	OutPrefix = paste0( OutDir,"celltyping_level1_" )
	pred_full = dan.CellTypeAnalysis_Full_xenium(xe, OutPrefix)
	load(file=paste0(OutDir,"AllMarkers.RData" ))
	am = am[(am$p_val_adj<0.01),]
	am$cluster = as.character(am$cluster)
	dir.create(paste0(OutPrefix,"enrichments/"))
	cluster_markers = list()
	for (cl in unique(am$cluster)){
		cluster_markers[[paste0("c",cl )]] = am[am$cluster==cl,"gene"]
	}
	dan.compare_published_celltypes(cluster_markers, paste0(OutPrefix,"enrichments/"), rownames(xe), printAllGenes = T, printTopIntersections = F)	
	# pdf(paste0(OutDir,"imagefeatureplot_cancer.pdf" ),12,12)
	# ImageFeaturePlot(xe, features = c("MKI67", "KRT8", "NKX2-1"), size = 0.4, cols = c("white", "red"))
	# dev.off()
	pdf(paste0(OutDir,"featureplot_cancer.pdf" ),20,8)
	FeaturePlot(xe, features = c("MKI67", "TOP2A","EPCAM","KRT8", "NKX2-1"), cols = c("white", "red"))
	dev.off()
	pdf(paste0(OutDir,"imagedimplot_clusters.pdf" ),12,12)
	ImageDimPlot(xe, size = 0.7,,border.size=NA )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()
	md = xe@meta.data
	# Cancer: positive for proliferation (TOP2A,MKI67) and epithelial (EPCAM). The other markers inspected on proteinatlas.org
	if (whichDataset=="xenium_nondiseasedlung"){ lev1_list = list( immune=c( 3,11,13:17,19,21,24,29,30,32,34,35,37 ), epithelial=c( 2,4,10,25,28,33,36,38:41 ), endothelial=c( 0,1,5,18,22,23 ), stroma=c( 6,7,8,9,12,20,26,27,31,42 ) ) }
	if (whichDataset=="xenium_05February2024"){ lev1_list = list( immune=c(8,11:13,15,17:20,22,23,29), epithelial=c(27,30), endothelial=c(14,21,28), stroma=c(4,10,24,25,31), Cancer=c(0:3,5,6,7,9,16,26,32,33)  ) } # 19 not super clear 
	if (whichDataset=="xenium_18May2023"){ lev1_list = list( immune=c(1,2,3,8,14,16,17,18,21,22,28,31,32,34,37), epithelial=c(4,6,19,23,24,29,35), endothelial=c(7,20,25,26,36), stroma=c(5,9,10,11,12,15,27,33), Cancer=c(0,13,30,38,39)  ) } # 
	if (whichDataset=="xenium_haga2023_TSU21"){ lev1_list = list( immune=c( 5,6,8:12,17,21,22,23,24,25,28,33,34 ), epithelial=c( 0,2,3,7,13,16,18,19,26,30,31,35,36,37,38 ), endothelial=c( 4,15,27,29,32 ), stroma=c( 1,14,20 )  ) } # 
	if (whichDataset=="xenium_haga2023_TSU20"){ lev1_list = list( immune=c( 4,10,16,22,24,26,30 ), epithelial=c( 0,2,3,5,6,7,11,12,13,15,17,20,31 ), endothelial=c( 8,9,23,27,28,29 ), stroma=c( 1,14,18,19,21,25 )  ) } # 
	if (whichDataset=="xenium_takano2024_luad2"){ lev1_list = list( immune=c( 10,15,17,24,25,31,35,40,46 ), epithelial=c( 2,5,6,8,9,11,16,28,32,36,38,39,41,42,43,44,47,48,49,50 ), endothelial=c( 3,21,29,30,45 ), stroma=c( 4,14 ), Cancer=c( 0,1,7,13,18,19,20,22,23,26,27,33,37 ), unclear = c( 12,34,51,52 )  ) } # 
	if (whichDataset=="xenium_takano2024_luad3"){ lev1_list = list( immune=c( 3,4,5,7,8,10,11,15,22,25,28,29,30 ), epithelial=c( 0,12,13,14,16,17,31 ), endothelial=c( 9,18,32 ), stroma=c( 1,2,6,19,20,21,23,24,26,27 )  ) } # 
	if (whichDataset=="xenium_takano2024_luad14"){ lev1_list = list( immune=c( 3,4,9,13,16,27,28,31,32,36 ), epithelial=c( 0,1,2,10,11,12,14,17,18,23,24,29,30,33,34 ), endothelial=c( 8,20,21,26,35 ), stroma=c( 6,7,15,19,22,37,38 ), unclear=c( 5,25 )  ) } # 
	if (whichDataset=="xenium_takano2024_luad16"){ lev1_list = list( immune=c( 3,5,9,10,12,17,18,20,25,26,27,31,32 ), epithelial=c( 0,2,6,16,22,24,28,30,33,35,36 ), endothelial=c( 14,15,21 ), stroma=c( 1,7,8,11,13,19,29 ),unclear=c( 4,23,34 )  ) } # 
	if (whichDataset=="xenium_takano2024_luad17"){ lev1_list = list( immune=c( 0,1,2,4,9,11,15,17,18,30,31,33 ), epithelial=c( 3,5,12,14,19,25,27 ), endothelial=c( 13,16,28,29 ), stroma=c( 6,7,8,20,21:24,26,32 ) ) } # 
	md$seurat_clusters = as.numeric(as.character(md$seurat_clusters))
	md$annd_level_1 = NA
	if (any(duplicated(as.numeric(unlist(lev1_list))))) { dcat( "Cluster assigned twice" ) }
	if ( !(identical(sort(as.numeric(unlist(lev1_list))),sort(unique(md$seurat_clusters)))) ) { dcat( "Some clusters are unassigned, or viceversa" ) }
	for (cl in sort(unique(md$seurat_clusters))){
		assigned = names(lev1_list)[sapply(lev1_list,  FUN=function(x) cl %in% x)]
    	if (length(assigned)!=1) { dcat( "Not unique assignment? weird" ) }
    	cellz = rownames(md[md$seurat_clusters==cl,])
    	md[cellz,"annd_level_1"] = assigned
	}
	dtable(md$annd_level_1)
	xe@meta.data = md

    pdf(paste0(OutDir,"celltyping_level1_assigned_imagedimplot.pdf"),12,8)
	print(ImageDimPlot(xe, group.by = "annd_level_1",border.size=NA) + coord_fixed()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
    dev.off()

	mdAll = xe@meta.data
	mdAll$annd_level_2 = NA

    ## epithelial
    md = xe@meta.data
    subxe = subset(xe,cells=rownames(md[md$annd_level_1=="epithelial",]))
    pdf(paste0(OutDir,"test_celltyping_level1_epithelial_imagedimplot.pdf"),12,8)
	print(ImageDimPlot(subxe,border.size=NA) + coord_fixed()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
    dev.off()
    OutPrefix = paste0( OutDir,"celltyping_level1epithelial_level2_" )
    subxe = SCTransform(subxe, assay = 'Xenium') # see https://github.com/satijalab/seurat/issues/1679 , but not regressing for cell cycle
	subxe = RunPCA(subxe)
	subxe = RunUMAP(subxe, dims = 1:30)
	subxe = FindNeighbors(subxe)
	subxe = FindClusters(subxe,res=1.5)
	pred_full = dan.CellTypeAnalysis_Full_xenium(subxe, OutPrefix, ct1="epithelial")
	md = subxe@meta.data
	save(pred_full,md, file = paste0(OutPrefix,"predfull_md.RData") )
    pdf(paste0(OutDir,"dimplot_epithelial_clusters.pdf" ),12,8)
	print(DimPlot(subxe,label=T) + coord_fixed())
	dev.off()
	pdf(paste0(OutDir,"imagedimplot_epithelial_clusters.pdf" ),12,8)
	print(ImageDimPlot(subxe,border.size=NA) + coord_fixed()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
	dev.off()
	xe@meta.data$epi_clusters = NA
	xe@meta.data[rownames(md),"epi_clusters"] = as.character(md$seurat_clusters)
	pdf(paste0(OutDir,"dimplot_epithelial_clusters_onXe.pdf" ),12,8)
	print(DimPlot(xe,group.by='epi_clusters',label=T) + coord_fixed())
	dev.off()
	am = FindAllMarkers( subxe, only.pos=TRUE )
	save(am,file=paste0(OutPrefix,"AllMarkers.RData"))
	load(file=paste0(OutPrefix,"AllMarkers.RData" ))
	am = am[(am$p_val_adj<0.01),]
	am$cluster = as.character(am$cluster)
	dir.create(paste0(OutPrefix,"enrichments/"))
	cluster_markers = list()
	for (cl in unique(am$cluster)){
		cluster_markers[[paste0("c",cl )]] = am[am$cluster==cl,"gene"]
	}
	dan.compare_published_mp(cluster_markers, paste0(OutPrefix,"enrichments/"), rownames(xe), printAllGenes = T, printTopIntersections = F)	
	if (whichDataset=="xenium_nondiseasedlung"){ lev2_list = list( AT1=c( 3,6,9,10,12,14,15,25,27 ), AT2=c( 0,2,4,5,7,13,16,18,20,21,26 ), club=c( 19,23 ), basal=c( 11,24,28 ), ciliated=c( 17 ),unclear=c( 1,8,22 ) ) } # 1,8
	if (whichDataset=="xenium_18May2023"){ lev2_list = list( AT1=c( 4,17,20), AT2=c(0,1,2,5,9,10,16,19 ), club=c(22), basal=c( 12), ciliated=c(13,15,21 ),unclear=c( 3,6,7,8,11,14,18 ) ) } # 19 not super clear 
	if (whichDataset=="xenium_05February2024"){ lev2_list = list( AT1=c( 5,10 ), AT2=c( 1 ), club=c( 7 ), basal=c( 6 ), ciliated=c( 2,4,8,12,13 ), unclear=c(0,3,9,11) ) } 
	if (whichDataset=="xenium_haga2023_TSU21"){ lev2_list = list( AT1=c( 0,1,7,22 ), AT2=c( 8,13,15,16,28 ), club=c( 2,10,27 ), basal=c( 33 ), ciliated=c( 30,32 ), Cancer=c( 3,4,5,6,9,12,14,19,21 ), unclear=c(11,17,18,20,23,24,25,26,29,31) ) } # cancer cells have HMGA1 / KRT8  markers
	if (whichDataset=="xenium_takano2024_luad2"){ lev2_list = list( AT2=c( 1,4,9,16,19,20,21,22 ), club=c( 8,12,14,23,28 ), basal=c( 13 ), unclear=c( 0,2,3,5,6,7,10,11,15,17,18,24,25,26,27,29 ) ) } # cancer cells have HMGA1 / KRT8  markers
	# if (whichDataset=="xenium_takano2024_luad3"){ lev2_list = list( AT1=c( 7,14,21 ), AT2=c( 0,8, ), club=c( 1,11,18,22 ), basal=c(  ), ciliated=c(  ), Cancer=c(  ), unclear=c( 6 ) ) } # impossible to discriminate normal / cancer cells
	# if (whichDataset=="xenium_takano2024_luad14"){ lev2_list = list(AT1=c(),AT2=c(15,18,20,28,29),club=c(1,3,4,5,13,16,21,22),basal=c(23),ciliated=c(14,24),Cancer=c(),unclear=c(0,30) ) } # impossible to discriminate normal / cancer cells
	# if (whichDataset=="xenium_takano2024_luad16"){ lev2_list = list( AT1=c(  ), AT2=c(  ), club=c(  ), basal=c(  ), ciliated=c(  ), Cancer=c(  ), unclear=c(  ) ) } 
	# if (whichDataset=="xenium_takano2024_luad17"){ lev2_list = list( AT1=c(  ), AT2=c(  ), club=c(  ), basal=c(  ), ciliated=c(  ), Cancer=c(  ), unclear=c(  ) ) } 
	load(file = paste0(OutPrefix,"predfull_md.RData") )
	md$seurat_clusters = as.numeric(as.character(md$seurat_clusters))
	md$annd_level_2 = NA
	if (any(duplicated(as.numeric(unlist(lev2_list))))) { dcat( "Cluster assigned twice" ) }
	if ( !(identical(sort(as.numeric(unlist(lev2_list))),sort(unique(md$seurat_clusters)))) ) { dcat( "Some clusters are unassigned, or viceversa" ) }
	for (cl in sort(unique(md$seurat_clusters))){
		assigned = names(lev2_list)[sapply(lev2_list,  FUN=function(x) cl %in% x)]
    	if (length(assigned)!=1) { dcat( "Not unique assignment? weird" ) }
    	cellz = rownames(md[md$seurat_clusters==cl,])
    	md[cellz,"annd_level_2"] = assigned
    	if ((whichDataset=="xenium_haga2023_TSU21") & (assigned=="Cancer")){ md[cellz,"annd_level_1"] = assigned }
	}
	dtable(md$annd_level_2)
	mdAll[rownames(md),"annd_level_2"] = md$annd_level_2

	## immune
	md = xe@meta.data
    subxe = subset(xe,cells=rownames(md[md$annd_level_1=="immune",]))
    OutPrefix = paste0( OutDir,"celltyping_level1immune_level2_" )
    subxe = SCTransform(subxe, assay = 'Xenium') # see https://github.com/satijalab/seurat/issues/1679 , but not regressing for cell cycle
	subxe = RunPCA(subxe)
	subxe = RunUMAP(subxe, dims = 1:30)
	subxe = FindNeighbors(subxe)
	subxe = FindClusters(subxe,res=1.5)
	pred_full = dan.CellTypeAnalysis_Full_xenium(subxe, OutPrefix, ct1="immune")
	md = subxe@meta.data
	save(pred_full,md, file = paste0(OutPrefix,"predfull_md.RData") )
    pdf(paste0(OutDir,"dimplot_immune_clusters.pdf" ),12,8)
	print(DimPlot(subxe,label=T) + coord_fixed())
	dev.off()
	am = FindAllMarkers( subxe, only.pos=TRUE )
	save(am,file=paste0(OutPrefix,"AllMarkers.RData" ))
	load(file=paste0(OutPrefix,"AllMarkers.RData" ))
	am = am[(am$p_val_adj<0.01),]
	am$cluster = as.character(am$cluster)
	dir.create(paste0(OutPrefix,"enrichments/"))
	cluster_markers = list()
	for (cl in unique(am$cluster)){
		cluster_markers[[paste0("c",cl )]] = am[am$cluster==cl,"gene"]
	}
	dan.compare_published_celltypes(cluster_markers, paste0(OutPrefix,"enrichments/"), rownames(xe), printAllGenes = T, printTopIntersections = F)	
	if (whichDataset=="xenium_18May2023"){ lev2_list = list( macrophages=c( 0,4,11,22,26,33,36 ), Tcells=c( 3,5,10,15,18,20,24,34,35 ), NKcells=c( 6 ), DC=c( 2,7,9,12,16,17,19,27,30 ), Bcells=c( 28 ), mast=c( 25,29 ), monocytes=c( 1 ), neutrophils = c( 14,32 ), unclear = c( 8,13,21,23,31 )) } # 19 not super clear 
	if (whichDataset=="xenium_05February2024"){ lev2_list = list( macrophages=c( 10,13 ), Tcells=c( 0,2,5,6,7,9,22 ), NKcells=c( 12,19,26 ), DC=c( 1,4,8,14,16,17,20,24,27 ), Bcells=c( 3,21 ), mast=c( 11,15 ), neutrophils = c( 18 ), unclear = c( 23,25 )) } # 19 not super clear 
	if (whichDataset=="xenium_haga2023_TSU21"){ lev2_list = list( non_epi = unique(as.numeric(as.character(md$seurat_clusters))) ) }
	if (whichDataset=="xenium_nondiseasedlung"){ lev2_list = list( non_epi = unique(as.numeric(as.character(md$seurat_clusters))) ) }
	load(file = paste0(OutPrefix,"predfull_md.RData") )
	md$seurat_clusters = as.numeric(as.character(md$seurat_clusters))
	md$annd_level_2 = NA
	if (any(duplicated(as.numeric(unlist(lev2_list))))) { dcat( "Cluster assigned twice" ) }
	if ( !(identical(sort(as.numeric(unlist(lev2_list))),sort(unique(md$seurat_clusters)))) ) { dcat( "Some clusters are unassigned, or viceversa" ) }
	for (cl in sort(unique(md$seurat_clusters))){
		assigned = names(lev2_list)[sapply(lev2_list,  FUN=function(x) cl %in% x)]
    	if (length(assigned)!=1) { dcat( "Not unique assignment? weird" ) }
    	cellz = rownames(md[md$seurat_clusters==cl,])
    	md[cellz,"annd_level_2"] = assigned
	}
	dtable(md$annd_level_2)
	mdAll[rownames(md),"annd_level_2"] = md$annd_level_2

	## endothelial
	md = xe@meta.data
    subxe = subset(xe,cells=rownames(md[md$annd_level_1=="endothelial",]))
    OutPrefix = paste0( OutDir,"celltyping_level1endothelial_level2_" )
    subxe = SCTransform(subxe, assay = 'Xenium') # see https://github.com/satijalab/seurat/issues/1679 , but not regressing for cell cycle
	subxe = RunPCA(subxe)
	subxe = RunUMAP(subxe, dims = 1:30)
	subxe = FindNeighbors(subxe)
	subxe = FindClusters(subxe,res=1.5)
	pred_full = dan.CellTypeAnalysis_Full_xenium(subxe, OutPrefix, ct1="endothelial")
	md = subxe@meta.data
	save(pred_full,md, file = paste0(OutPrefix,"predfull_md.RData") )
    pdf(paste0(OutDir,"dimplot_endothelial_clusters.pdf" ),12,8)
	print(DimPlot(subxe,label=T) + coord_fixed())
	dev.off()
	am = FindAllMarkers( subxe, only.pos=TRUE )
	save(am,file=paste0(OutPrefix,"AllMarkers.RData" ))
	load(file=paste0(OutPrefix,"AllMarkers.RData" ))
	am = am[(am$p_val_adj<0.01),]
	am$cluster = as.character(am$cluster)
	dir.create(paste0(OutPrefix,"enrichments/"))
	cluster_markers = list()
	for (cl in unique(am$cluster)){
		cluster_markers[[paste0("c",cl )]] = am[am$cluster==cl,"gene"]
	}
	dan.compare_published_celltypes(cluster_markers, paste0(OutPrefix,"enrichments/"), rownames(xe), printAllGenes = T, printTopIntersections = F)	
	if (whichDataset=="xenium_18May2023"){ lev2_list = list( EC_lymphatic=c( 12,22 ), EC_blood=c( 0,1,4,5,7:11,13:17,19:21,23:25 ), pericytes = c( 2,3,6,18 ) ) } # 19 not super clear 
	if (whichDataset=="xenium_05February2024"){ lev2_list = list( EC_lymphatic=c( 11 ), EC_blood=c( 0,1,3:8,10,12:15 ), pericytes = c( 2,9,16,17 ) ) } # 19 not super clear 
	if (whichDataset=="xenium_haga2023_TSU21"){ lev2_list = list( non_epi = unique(as.numeric(as.character(md$seurat_clusters))) ) }
	load(file = paste0(OutPrefix,"predfull_md.RData") )
	md$seurat_clusters = as.numeric(as.character(md$seurat_clusters))
	md$annd_level_2 = NA
	if (any(duplicated(as.numeric(unlist(lev2_list))))) { dcat( "Cluster assigned twice" ) }
	if ( !(identical(sort(as.numeric(unlist(lev2_list))),sort(unique(md$seurat_clusters)))) ) { dcat( "Some clusters are unassigned, or viceversa" ) }
	for (cl in sort(unique(md$seurat_clusters))){
		assigned = names(lev2_list)[sapply(lev2_list,  FUN=function(x) cl %in% x)]
    	if (length(assigned)!=1) { dcat( "Not unique assignment? weird" ) }
    	cellz = rownames(md[md$seurat_clusters==cl,])
    	md[cellz,"annd_level_2"] = assigned
	}
	dtable(md$annd_level_2)
	mdAll[rownames(md),"annd_level_2"] = md$annd_level_2

	## stroma
	md = xe@meta.data
    subxe = subset(xe,cells=rownames(md[md$annd_level_1=="stroma",]))
    OutPrefix = paste0( OutDir,"celltyping_level1stroma_level2_" )
    subxe = SCTransform(subxe, assay = 'Xenium') # see https://github.com/satijalab/seurat/issues/1679 , but not regressing for cell cycle
	subxe = RunPCA(subxe)
	subxe = RunUMAP(subxe, dims = 1:30)
	subxe = FindNeighbors(subxe)
	subxe = FindClusters(subxe,res=1.5)
	pred_full = dan.CellTypeAnalysis_Full_xenium(subxe, OutPrefix, ct1="stroma")
	md = subxe@meta.data
	save(pred_full,md, file = paste0(OutPrefix,"predfull_md.RData") )
    pdf(paste0(OutDir,"dimplot_stroma_clusters.pdf" ),12,8)
	print(DimPlot(subxe,label=T) + coord_fixed())
	dev.off()
	am = FindAllMarkers( subxe, only.pos=TRUE )
	save(am,file=paste0(OutPrefix,"AllMarkers.RData" ))
	load(file=paste0(OutPrefix,"AllMarkers.RData" ))
	am = am[(am$p_val_adj<0.01),]
	am$cluster = as.character(am$cluster)
	dir.create(paste0(OutPrefix,"enrichments/"))
	cluster_markers = list()
	for (cl in unique(am$cluster)){
		cluster_markers[[paste0("c",cl )]] = am[am$cluster==cl,"gene"]
	}
	dan.compare_published_celltypes(cluster_markers, paste0(OutPrefix,"enrichments/"), rownames(xe), printAllGenes = T, printTopIntersections = F)	
	if (whichDataset=="xenium_18May2023"){ lev2_list = list( fibroblasts=c( 1:6,10,14:17,19,20,22,23,25,26,28 ), mesothelial=c( 7,8,21,24,32 ), smooth_muscle=c( 0,9,11,12,13,18,27,30 ), unclear = c(29,31) ) } # 19 not super clear 
	if (whichDataset=="xenium_05February2024"){ lev2_list = list( fibroblasts=c( 0:2,5,6,8:11,13,14,17:20 ), mesothelial=c( 15,16 ), smooth_muscle=c( 3,21 ), unclear = c( 4,7,12 ) ) }
	if (whichDataset=="xenium_haga2023_TSU21"){ lev2_list = list( non_epi = unique(as.numeric(as.character(md$seurat_clusters))) )}
	load(file = paste0(OutPrefix,"predfull_md.RData") )
	md$seurat_clusters = as.numeric(as.character(md$seurat_clusters))
	md$annd_level_2 = NA
	if (any(duplicated(as.numeric(unlist(lev2_list))))) { dcat( "Cluster assigned twice" ) }
	if ( !(identical(sort(as.numeric(unlist(lev2_list))),sort(unique(md$seurat_clusters)))) ) { dcat( "Some clusters are unassigned, or viceversa" ) }
	for (cl in sort(unique(md$seurat_clusters))){
		assigned = names(lev2_list)[sapply(lev2_list,  FUN=function(x) cl %in% x)]
    	if (length(assigned)!=1) { dcat( "Not unique assignment? weird" ) }
    	cellz = rownames(md[md$seurat_clusters==cl,])
    	md[cellz,"annd_level_2"] = assigned
	}
	dtable(md$annd_level_2)
	mdAll[rownames(md),"annd_level_2"] = md$annd_level_2
	mdAll[mdAll$annd_level_1=="Cancer","annd_level_2"] = "Cancer"
	mdAll[is.na(mdAll$annd_level_2),"annd_level_2"] = "others"
	dtable(mdAll$annd_level_2)
	xe@meta.data = mdAll[rownames(xe@meta.data),]
	pdf(paste0(OutDir,"full_imageplot_celltypes_level1.pdf" ),12,12)
	ImageDimPlot(xe,group.by="annd_level_1",border.size=NA)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()
	pdf(paste0(OutDir,"full_dimplot_celltypes_level1.pdf" ),12,8)
	print(DimPlot(xe,group.by="annd_level_1",label=T) + coord_fixed())
	dev.off()
	pdf(paste0(OutDir,"full_imageplot_celltypes_level2.pdf" ),12,12)
	ImageDimPlot(xe,group.by="annd_level_2",border.size=NA)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()
	pdf(paste0(OutDir,"full_dimplot_celltypes_level2.pdf" ),12,8)
	print(DimPlot(xe,group.by="annd_level_2",label=T) + coord_fixed())
	dev.off()

	
	ct_colorz = c( Cancer="firebrick",epithelial="tomato",stroma="gray",immune="dodgerblue3",endothelial="goldenrod1" )
	xe@meta.data$annd_level_1 = factor(xe@meta.data$annd_level_1,levels=names(ct_colorz))
    pdf(paste0(OutDir,"celltyping_level1_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
    plot=(DimPlot(xe,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_1", label = T,raster=T,cols=ct_colorz) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot)
    dev.off()
    pdf(paste0(OutDir,"celltyping_level1_assigned_unrastered.pdf"),2.5,2.5,pointsize=6)
    plot=(dan.DimPlot(xe,pt.size=0.1,label.size=6/.pt, group.by = "annd_level_1", label = F,raster=F,cols=ct_colorz) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()
    pdf(paste0(OutDir,"celltyping_level1_assigned_ImagePlot.pdf"),4,2.5,pointsize=6)
    coords = GetTissueCoordinates(xe, which = "centroids")
	plot=ImageDimPlot(xe,group.by="annd_level_1",axes=T,size=0.08,border.size=NA,cols=ct_colorz)+scale_y_continuous(breaks=c( 0,2500,5000,7500,10000 ))+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()

    thisxe = subset(xe, cells = rownames(xe@meta.data[xe@meta.data$annd_level_2!="unclear",]) )
    ct_colorz = c( EC_blood="gold",EC_lymphatic="goldenrod1",pericytes="gold3",Cancer="firebrick",
							AT1="darkorange1",AT2="orangered",basal="sienna4",ciliated="tan3",club="sandybrown",
							Bcells="forestgreen",DC="deeppink2",macrophages="purple",mast="purple4",monocytes="mediumpurple1",neutrophils="deeppink4",NKcells="midnightblue",Tcells="dodgerblue2",
							fibroblasts="gray22",mesothelial="gray55",smooth_muscle="gray88")
	thisxe@meta.data$annd_level_2 = factor(thisxe@meta.data$annd_level_2,levels=names(ct_colorz))
	pdf(paste0(OutDir,"celltyping_level2_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
	plot=(DimPlot(thisxe,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_2",label = T,raster=T,cols=ct_colorz) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) + guides(color=guide_legend(ncol =2,override.aes = list(size=2)))
    print(plot)
    dev.off()
    pdf(paste0(OutDir,"celltyping_level2_assigned_unrastered.pdf"),2.5,2.5,pointsize=6)
    plot=(dan.DimPlot(thisxe,pt.size=0.1,label.size=6/.pt, group.by = "annd_level_2", label = F,raster=F,cols=ct_colorz) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()
    pdf(paste0(OutDir,"celltyping_level2_assigned_ImagePlot.pdf"),4,2.5,pointsize=6)
    coords = GetTissueCoordinates(thisxe, which = "centroids")
	plot=ImageDimPlot(thisxe,group.by="annd_level_2",axes=T,size=0.08,border.size=NA,cols=ct_colorz)+scale_y_continuous(breaks=c( 0,2500,5000,7500,10000 ))+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()

    # alveolar and cancer cells
    thisxe = subset(xe, cells = rownames(xe@meta.data[xe@meta.data$annd_level_2 %in% c( "Cancer","AT1","AT2" ),]) )
    ct_colorz = c( Cancer="#1C75BC",AT1="#594A42", AT2="#39B54A")
	thisxe@meta.data$annd_level_2 = factor(thisxe@meta.data$annd_level_2,levels=names(ct_colorz))
	pdf(paste0(OutDir,"celltyping_CancerAlveolar_rastered_withLegend.pdf"),6,6,pointsize=6)
	plot=(DimPlot(thisxe,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_2",label = T,raster=T,cols=ct_colorz) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) + guides(color=guide_legend(ncol =2,override.aes = list(size=2)))
    print(plot)
    dev.off()
    pdf(paste0(OutDir,"celltyping_CancerAlveolar_unrastered.pdf"),2.5,2.5,pointsize=6)
    plot=(dan.DimPlot(thisxe,pt.size=0.1,label.size=6/.pt, group.by = "annd_level_2", label = F,raster=F,cols=ct_colorz) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()
    pdf(paste0(OutDir,"celltyping_CancerAlveolar_ImagePlot.pdf"),4,2.5,pointsize=6)
    coords = GetTissueCoordinates(thisxe, which = "centroids")
	plot=ImageDimPlot(thisxe,group.by="annd_level_2",axes=T,size=0.08,border.size=NA,cols=ct_colorz)+scale_y_continuous(breaks=c( 0,2500,5000,7500,10000 ))+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.1),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()

    pdf(paste0(OutDir,"celltyping_CancerAlveolar_ImagePlot_whiteBg.pdf"),3,2.2,pointsize=6)
    coords = GetTissueCoordinates(thisxe, which = "centroids")
	plot=ImageDimPlot(thisxe,group.by="annd_level_2",axes=T,size=0.15,border.size=NA,cols=ct_colorz,dark.background=F)+theme_classic()+scale_y_continuous(breaks=c( 0,2500,5000,7500,10000 ))+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.1),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()

    

    pdf(paste0(OutDir,"celltyping_level2_assigned.pdf"),2.1,2.1,pointsize=6)
	plot=(DimPlot(xe,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_2", label = T,raster=T) + coord_fixed())
	plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot+NoLegend())
    dev.off()
    # thisxe = subset(xe, cells = rownames(mdAll[mdAll$annd_level_1 %in% c( "endothelial","epithelial" ),]) )
    

	mdAll = xe@meta.data
	xe_EpiCancer = subset(xe, cells = rownames(mdAll[mdAll$annd_level_1 %in% c( "Cancer","epithelial" ),]) )
	pdf(paste0(OutDir,"fEpiCancer_imageplot_celltypes_level1.pdf" ),12,12)
	ImageDimPlot(xe_EpiCancer,group.by="annd_level_1",border.size=NA)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()
	pdf(paste0(OutDir,"fEpiCancer_imageplot_celltypes_level2.pdf" ),12,12)
	ImageDimPlot(xe_EpiCancer,group.by="annd_level_2",border.size=NA)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()

	xe_AT12Cancer = subset(xe, cells = rownames(mdAll[mdAll$annd_level_2 %in% c( "Cancer","AT1","AT2" ),]) )
	pdf(paste0(OutDir,"fAT12Cancer_imageplot_celltypes_level1.pdf" ),12,12)
	ImageDimPlot(xe_AT12Cancer,group.by="annd_level_1",border.size=NA)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()

	pdf(paste0(OutDir,"fAT12Cancer_imageplot_celltypes_level2.pdf" ),12,12)
	ImageDimPlot(xe_AT12Cancer,group.by="annd_level_2",border.size=NA)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()

	save(xe_AT12Cancer,file=paste0(OutDir,"xe_AT12Cancer.RData" ))
	save(xe,file=paste0(OutDir,"xe.RData" ))
	load(file=paste0(OutDir,"xe.RData" ))
	mdAll = xe@meta.data
	save(mdAll,file=paste0(OutDir,"mdAll.RData"))
	dtable(xe_AT12Cancer@annd_level_2)
	dtable(xe_AT12Cancer@annd_level_1)
}

xenium_analyses = function( OutDir, whichDataset ){
	load(file=paste0(OutDir,"../xenium_preprocessing_",whichDataset,"/xe_AT12Cancer.RData" ))
	pdf(paste0(OutDir,"fAT12Cancer_imageplot_celltypes_level1.pdf" ),7,7)
	coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	plot=ImageDimPlot(xe_AT12Cancer,group.by="annd_level_1",axes=T,border.size=NA,size=0.4)+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot)
	dev.off()
	ct_colorz = c( "Cancer"="#1C75BC","AT1"="#594A42","AT2"="#39B54A" )
	# ct_colorz = c("Cancer"="firebrick","AT1"="chocolate","AT2"="dodgerblue2")
	# pdf(paste0(OutDir,"fAT12Cancer_imageplot_celltypes_level2.pdf" ),7,7)
	# plot=ImageDimPlot(xe_AT12Cancer,group.by="annd_level_2",axes=T,cols=ct_colorz,border.size=NA,size=0.4)+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
	# print(plot)
	# dev.off()
	pdf(paste0(OutDir,"fAT12Cancer_imageplot_celltypes_level2.pdf"),3.5,2.5,pointsize=6)
    coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	plot=ImageDimPlot(xe_AT12Cancer,group.by="annd_level_2",axes=T,size=0.15,border.size=NA,cols=ct_colorz,dark.background=F)+labs(y="x (µm)",x="y (µm)",fill="Cell type")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot)
    dev.off()

	pdf(paste0(OutDir,"fAT12Cancer_DimPlot_celltypes_level1.pdf" ),12,12)
	plot=DimPlot(xe_AT12Cancer,group.by="annd_level_1")
	print(plot)
	dev.off()
	pdf(paste0(OutDir,"fAT12Cancer_DimPlot_celltypes_level2.pdf" ),12,12)
	plot=DimPlot(xe_AT12Cancer,group.by="annd_level_2",cols=ct_colorz)
	print(plot)
	dev.off()

	md = xe_AT12Cancer@meta.data
	coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	rownames(coords) = coords$cell
	md$x = coords[rownames(md),"x"]
	md$y = coords[rownames(md),"y"]

	library(rdist)
	md$nearest_cancer_cell = NA
	md$nearest_dist = NA
	# mdtest = md[sample(rownames(md),1000),]
	md_cancer = md[md$annd_level_2=="Cancer",c( "x","y" )]
	md_at2 = md[md$annd_level_2=="AT2",c( "x","y" )]
	dd_CancerAt2 = cdist(md_cancer,md_at2)
	# ### for takano2024
	# dd_CancerAt2_1 = cdist(md_cancer[1:70000,],md_at2)
	# dd_CancerAt2_2 = cdist(md_cancer[70001:nrow(md_cancer),],md_at2)
	# dd_CancerAt2 = rbind(dd_CancerAt2_1,dd_CancerAt2_2)

	rownames(dd_CancerAt2) = rownames(md_cancer)
	colnames(dd_CancerAt2) = rownames(md_at2)
	min_dist = apply(dd_CancerAt2,2,min)
	min_cancer_cell = rownames(dd_CancerAt2)[apply(dd_CancerAt2,2,which.min)]
	md[colnames(dd_CancerAt2),"nearest_cancer_cell"] = min_cancer_cell
	md[colnames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)

	md_cancer = md[md$annd_level_2=="Cancer",c( "x","y" )]
	md_at1 = md[md$annd_level_2=="AT1",c( "x","y" )]
	dd_CancerAt1 = cdist(md_cancer,md_at1)
	rownames(dd_CancerAt1) = rownames(md_cancer)
	colnames(dd_CancerAt1) = rownames(md_at1)
	min_dist = apply(dd_CancerAt1,2,min)
	min_cancer_cell = rownames(dd_CancerAt1)[apply(dd_CancerAt1,2,which.min)]
	md[colnames(dd_CancerAt1),"nearest_cancer_cell"] = min_cancer_cell
	md[colnames(dd_CancerAt1),"nearest_dist"] = as.numeric(min_dist)

	head(md)
	sum(is.na(md$nearest_cancer_cell))
	sum(is.na(md$nearest_dist))

	save(md, file = paste0( OutDir,"md_save.RData" ))
	load( file = paste0( OutDir,"md_save.RData" ))


	xe_AT12Cancer@meta.data = md
	hlca_mp = list()
	aa = dan.read(file = "data/cell_typing_resources/HLCA_epithelial_markers.txt")
	for (cn in colnames(aa)){
		hlca_mp[[paste0("hlca_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
	}
	names(hlca_mp)
	hlca_mp = hlca_mp[c( "hlca_AT1","hlca_AT2" )]
	load(paste0(OutDir,"../weak_AT2_signature/markers_correctingCellType_correctingPatient.RData" ))
	weakAT2_up = rownames(mm[(mm$avg_log2FC>(+0.25)) & (mm$p_val_adj<0.01),])
	weakAT2_down = rownames(mm[(mm$avg_log2FC<(-0.25)) & (mm$p_val_adj<0.01),])
	hkc = dan.read("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/han_kac_deg_complete_table8.txt")
	hkc_top = hkc[order(hkc$avg_log2FC,decreasing=T),"gene"][1:100]
	gs = c(hlca_mp,list(weakAT2_up=weakAT2_up,weakAT2_down=weakAT2_down,KAC_signature=hkc_top))
	save(gs,file=paste0(OutDir,"gs.RData"))
	# xe_AT12Cancer = AddModuleScore(xe_AT12Cancer,gs,pool = NULL,nbin = 10,ctrl = 20,k = FALSE,assay = 'SCT',name = "ams")
	# colnames(xe_AT12Cancer@meta.data)[substr(colnames(xe_AT12Cancer@meta.data),1,3)=="ams"] = names(gs)
	dcat( "AmsCentered" )
	db_rand = dan.barkley_MakeRand_data(xe_AT12Cancer,gs, 3)
	scores = dan.Barkley_GeneToEnrichment_AmsCentered( xe_AT12Cancer, gs, db_rand)
	xe_AT12Cancer@meta.data = scores
	md = xe_AT12Cancer@meta.data

	pdf( paste0(OutDir,"alveolar_signatures_byCellType.pdf"),6,5 )
	for (tp in names(gs)){
	   x = md$annd_level_2
	   y = md[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   print(plot)
	}
	dev.off()

	md_at1 = md[md$annd_level_2=="AT1",]
	pdf( paste0(OutDir,"AT1_alveolar_signatures_vsDistance.pdf"),6,5 )
	for (tp in names(gs)){
	   x = md_at1$nearest_dist
	   y = md_at1[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer cell", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	md_at2 = md[md$annd_level_2=="AT2",]
	pdf( paste0(OutDir,"AT2_alveolar_signatures_vsDistance.pdf"),6,5 )
	for (tp in names(gs)){
	   x = md_at2$nearest_dist
	   y = md_at2[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer cell", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()
	save(gs,file=paste0(OutDir,"gs.RData"))
	save(md,file = paste0( OutDir,"md_allCells_allValues.RData" ))

	grid_spacing = 100
	library(viridis)
	md_at2$dist_discrete = ">=1mm"
	md_at2[md_at2$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_at2[(md_at2$nearest_dist>=lowerbound) & (md_at2$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(OutDir,"AT2_alveolar_signatures_vsDistance_Boxplots_SameGridSpacing.pdf"),8,5 )
	for (tp in intersect(names(gs),colnames(md_at2)) ){
	   x = md_at2$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer cell", ylab = paste0("mean ",tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	md_at1$dist_discrete = ">=1mm"
	md_at1[md_at1$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_at1[(md_at1$nearest_dist>=lowerbound) & (md_at1$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(OutDir,"AT1_alveolar_signatures_vsDistance_Boxplots_SameGridSpacing.pdf"),8,5 )
	for (tp in intersect(names(gs),colnames(md_at1)) ){
	   x = md_at1$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer cell", ylab = paste0("mean ",tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	load(file=paste0(OutDir,"../xenium_preprocessing_",whichDataset,"/xe_AT12Cancer.RData" ))
	rn = rownames(xe_AT12Cancer)
	save(rn,file=paste0(OutDir,"rownames.RData"))
	load(paste0(OutDir,"gs.RData"))
	gs_hlca_available = list()
	for (nn in c( "hlca_AT1","hlca_AT2" )){ gs_hlca_available[[nn]] = intersect(gs[[nn]],rownames(xe_AT12Cancer)) }
	save(gs_hlca_available,file=paste0(OutDir,"gs_hlca_available.RData"))
	load(paste0(OutDir,"../tps_discovery/tps.RData"))
	tps_available = list()
	for (nn in c( names(tps) )){ tps_available[[nn]] = intersect(tps[[nn]],rownames(xe_AT12Cancer)) }
	save(tps_available,file=paste0(OutDir,"tps_available.RData"))
	load( file=paste0(OutDir,"../CellStates/signatures/wilcox_cs_signatures_restrictive.RData" ) )
	cs_available = list()
	for (nn in c( names(cs_signatures) )){ cs_available[[nn]] = intersect(cs_signatures[[nn]],rownames(xe_AT12Cancer)) }
	save(cs_available,file=paste0(OutDir,"cs_available.RData"))
	dir.create( paste0(OutDir,"cs_available_featurePlots/" ))
	coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	cells_zoomin = coords[((coords$x<9000) & (coords$x>500)) & ((coords$y<8000)),'cell' ]
	subxe = subset(xe_AT12Cancer,cells=intersect(rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="Cancer",]),cells_zoomin) )
	for (cs in names(cs_available)){
		dcat(cs)
		for (g in cs_available[[cs]]){
			dcat(g,1)
			pdf(paste0(OutDir,"cs_available_featurePlots/",cs,"_",g,".pdf"),3,2.2,pointsize=6)
		    coords = GetTissueCoordinates(subxe, which = "centroids")
			plot=ImageFeaturePlot(subxe,features=g,axes=T,size=0.07,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill=g)+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		    print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		    dev.off()
		}	
	}
	dir.create( paste0(OutDir,"cs_available_featurePlots_smoothed/" ))
	for (cs in names(cs_available)){
		dcat(cs)
		for (g in cs_available[[cs]]){
			dcat(g,1)
			subxe@meta.data[,g] = as.numeric(subxe@assays$SCT@data[g,])
			subxe = nn_smoothing(subxe, feature = g, radius = 50)
			pdf(paste0(OutDir,"cs_available_featurePlots_smoothed/",cs,"_",g,".pdf"),2.7,2.2,pointsize=6)
		    coords = GetTissueCoordinates(subxe, which = "centroids")
			plot=ImageFeaturePlot(subxe,features=paste0("smoothed_",g),axes=T,size=0.07,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		    print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		    dev.off()
		}	
	}
	dir.create( paste0(OutDir,"tps_available_featurePlots_smoothed/" ))
	for (cs in names(tps_available)){
		dcat(cs)
		for (g in tps_available[[cs]]){
			dcat(g,1)
			subxe@meta.data[,g] = as.numeric(subxe@assays$SCT@data[g,])
			subxe = nn_smoothing(subxe, feature = g, radius = 50)
			pdf(paste0(OutDir,"tps_available_featurePlots_smoothed/",cs,"_",g,".pdf"),3,2.2,pointsize=6)
		    coords = GetTissueCoordinates(subxe, which = "centroids")
			plot=ImageFeaturePlot(subxe,features=paste0("smoothed_",g),axes=T,size=0.07,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		    print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		    dev.off()
		}	
	}
}

xenium_analyses_nondiseasedlung_vs_tumor = function( OutDir,density_threshold = 0.002 ){

	load(file=paste0(OutDir,"../xenium_preprocessing_xenium_nondiseasedlung/xe_AT12Cancer.RData" ))
	nondiseased = xe_AT12Cancer
	nondiseased@meta.data$condition = "nondiseased"
	load(file=paste0(OutDir,"../xenium_preprocessing_xenium_18May2023/xe_AT12Cancer.RData" ))
	diseased = xe_AT12Cancer
	diseased@meta.data$condition = "diseased"

	nondiseased_data = as.matrix(nondiseased@assays$SCT@data)
	colnames(nondiseased_data) = paste0( "n",colnames(nondiseased_data) )
	diseased_data = as.matrix(diseased@assays$SCT@data)
	colnames(diseased_data) = paste0( "d",colnames(diseased_data) )
	data = cbind(nondiseased_data,diseased_data[rownames(nondiseased_data),])

	nondiseased_md = nondiseased@meta.data
	rownames(nondiseased_md) = paste0( "n",rownames(nondiseased_md) )
	diseased_md = diseased@meta.data
	rownames(diseased_md) = paste0( "d",rownames(diseased_md) )
	annot = rbind(nondiseased_md[,c( "annd_level_1","annd_level_2","condition" )],diseased_md[,c( "annd_level_1","annd_level_2","condition" )])

	# hlca_mp = list()
	# aa = dan.read(file = "data/cell_typing_resources/HLCA_epithelial_markers.txt")
	# for (cn in colnames(aa)){
	# 	hlca_mp[[paste0("hlca_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
	# }
	# names(hlca_mp)
	# hlca_mp = hlca_mp[c( "hlca_AT1","hlca_AT2" )]
	# load(paste0(OutDir,"../weak_AT2_signature/markers_correctingCellType_correctingPatient.RData" ))
	# weakAT2_up = rownames(mm[(mm$avg_log2FC>(+0.25)) & (mm$p_val_adj<0.01),])
	# weakAT2_down = rownames(mm[(mm$avg_log2FC<(-0.25)) & (mm$p_val_adj<0.01),])
	# hkc = dan.read("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/han_kac_deg_complete_table8.txt")
	# hkc_top = hkc[order(hkc$avg_log2FC,decreasing=T),"gene"][1:100]
	# gs = c(hlca_mp,list(weakAT2_up=weakAT2_up,weakAT2_down=weakAT2_down,KAC_signature=hkc_top))
	# save(gs,file=paste0(OutDir,"gs.RData"))

	# db_rand = dan.barkley_MakeRand_data(data,gs, 3)
	# scores = dan.Barkley_GeneToEnrichment_AmsCentered_data( data, annot, gs, db_rand)
	# save(scores,file=paste0(OutDir,"scores.RData"))

	load(file=paste0(OutDir,"gs.RData"))
	load(file=paste0(OutDir,"scores.RData"))

	md = scores
	pdf( paste0(OutDir,"alveolar_signatures_byCellType.pdf"),3,2.5 )
	for (tp in names(gs)){
	   x = md$annd_level_2
	   y = md[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0(tp," score"), plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = 0, jitterColors = "black", labelJitteredPoints = NULL, includeJitters = F )
	   print(plot)
	}
	dev.off()

	mean(md[(md$condition=="nondiseased") & (md$annd_level_2=="AT1"), "hlca_AT1" ])
	mean(md[(md$condition=="nondiseased") & (md$annd_level_2=="AT2"), "hlca_AT2" ])

	coords = GetTissueCoordinates(diseased, which = "centroids")
	rownames(coords) = paste0( "d",coords$cell)
	diseased_md$x = coords[rownames(diseased_md),"x"]
	diseased_md$y = coords[rownames(diseased_md),"y"]

	coords = GetTissueCoordinates(nondiseased, which = "centroids")
	rownames(coords) = paste0( "n",coords$cell)
	nondiseased_md$x = coords[rownames(nondiseased_md),"x"]
	nondiseased_md$y = coords[rownames(nondiseased_md),"y"]
	
	### diseased
	grid_spacing = 200
	newscores = scores
	diseased_md = cbind(diseased_md,newscores[rownames(diseased_md),names(gs)])
	md = diseased_md
	xmax = ceiling(max(md$x)/grid_spacing)
  	ymax = ceiling(max(md$y)/grid_spacing)
	md$xg = factor(ceiling(md$x/grid_spacing), levels = c(1:xmax))
  	md$yg = factor(ceiling(md$y/grid_spacing), levels = c(1:ymax))
  	md$grid_coo = paste0( "xg",md$xg,"_yg",md$yg )
	library(rdist)
	md$nearest_cancer_cell = NA
	md$nearest_dist = NA
	md_cancer = md[md$annd_level_2=="Cancer",c( "x","y" )]
	md_at2 = md[md$annd_level_2=="AT2",c( "x","y" )]
	dd_CancerAt2 = cdist(md_cancer,md_at2)
	rownames(dd_CancerAt2) = rownames(md_cancer)
	colnames(dd_CancerAt2) = rownames(md_at2)
	min_dist = apply(dd_CancerAt2,2,min)
	min_cancer_cell = rownames(dd_CancerAt2)[apply(dd_CancerAt2,2,which.min)]
	md[colnames(dd_CancerAt2),"nearest_cancer_cell"] = min_cancer_cell
	md[colnames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)
	md_cancer = md[md$annd_level_2=="Cancer",c( "x","y" )]
	md_at1 = md[md$annd_level_2=="AT1",c( "x","y" )]
	dd_CancerAt1 = cdist(md_cancer,md_at1)
	rownames(dd_CancerAt1) = rownames(md_cancer)
	colnames(dd_CancerAt1) = rownames(md_at1)
	min_dist = apply(dd_CancerAt1,2,min)
	min_cancer_cell = rownames(dd_CancerAt1)[apply(dd_CancerAt1,2,which.min)]
	md[colnames(dd_CancerAt1),"nearest_cancer_cell"] = min_cancer_cell
	md[colnames(dd_CancerAt1),"nearest_dist"] = as.numeric(min_dist)

	celltype_map = data.frame(row.names=c( "AT1","AT2","Cancer" ),colorz = c("chocolate","dodgerblue2","firebrick") ,stringsAsFactors=F )
	this_OutDir = OutDir
	# AT1
  	md_at1 = md[md$annd_level_2=="AT1",]
  	md_at1$xg = factor(ceiling(md_at1$x/grid_spacing), levels = c(1:xmax))
  	md_at1$yg = factor(ceiling(md_at1$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at1 = (table(md_at1$xg, md_at1$yg))
  	md_at1$grid_coo = paste0( "xg",md_at1$xg,"_yg",md_at1$yg )
  	md_at1_grid = aggregate( .~grid_coo,data=md_at1[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_at1_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at1_grid$grid_coo)))
  	md_at1_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at1_grid$grid_coo)))
  	rownames(md_at1_grid) = md_at1_grid$grid_coo
  	md_at1_grid$n_AT1 = table(md_at1$grid_coo)[rownames(md_at1_grid)]
  	# AT2
  	md_at2 = md[md$annd_level_2=="AT2",]
  	md_at2$xg = factor(ceiling(md_at2$x/grid_spacing), levels = c(1:xmax))
  	md_at2$yg = factor(ceiling(md_at2$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at2 = (table(md_at2$xg, md_at2$yg))
  	md_at2$grid_coo = paste0( "xg",md_at2$xg,"_yg",md_at2$yg )
  	md_at2_grid = aggregate( .~grid_coo,data=md_at2[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_at2_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at2_grid$grid_coo)))
  	md_at2_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at2_grid$grid_coo)))
  	rownames(md_at2_grid) = md_at2_grid$grid_coo
  	md_at2_grid$n_AT2 = table(md_at2$grid_coo)[rownames(md_at2_grid)]
  	# Cancer
  	md_cancer = md[md$annd_level_2=="Cancer",]
  	md_cancer$xg = factor(ceiling(md_cancer$x/grid_spacing), levels = c(1:xmax))
  	md_cancer$yg = factor(ceiling(md_cancer$y/grid_spacing), levels = c(1:ymax))
  	meaMat_cancer = (table(md_cancer$xg, md_cancer$yg))
  	md_cancer$grid_coo = paste0( "xg",md_cancer$xg,"_yg",md_cancer$yg )
  	md_cancer_grid = aggregate( .~grid_coo,data=md_cancer[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_cancer_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_cancer_grid$grid_coo)))
  	md_cancer_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_cancer_grid$grid_coo)))
  	rownames(md_cancer_grid) = md_cancer_grid$grid_coo
  	md_cancer_grid$n_Cancer = table(md_cancer$grid_coo)[rownames(md_cancer_grid)]

	md_at1_grid$n_AT2 = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"]
	md_at1_grid$n_Cancer = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_Cancer"] = md_cancer_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_Cancer"]
	md_at2_grid$n_AT1 = 0
	md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"]
	md_at2_grid$n_Cancer = 0
	md_at2_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_Cancer"] = md_cancer_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_Cancer"]
	md_cancer_grid$n_AT1 = 0
	md_cancer_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_AT1"]
	md_cancer_grid$n_AT2 = 0
	md_cancer_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_AT2"]
	md_cancer_grid$ratio_Cancer = md_cancer_grid$n_Cancer/(md_cancer_grid$n_Cancer+md_cancer_grid$n_AT1+md_cancer_grid$n_AT2)
	md_cancer_grid$spatial_density_Cancer = md_cancer_grid$n_Cancer/(grid_spacing*grid_spacing)
	md_cancer_grid$cancer_rich = 0
	md_cancer_grid[(md_cancer_grid$spatial_density_Cancer>density_threshold),"cancer_rich"] = 1
	# md_cancer_rich = md_cancer_grid[(md_cancer_grid$ratio_Cancer>cancer_rich_ratio) & (md_cancer_grid$spatial_density_Cancer>density_threshold),c( "xg","yg" )]
	md_cancer_rich = md_cancer_grid[(md_cancer_grid$spatial_density_Cancer>density_threshold),c( "xg","yg" )]
	library(rdist)
	dd_CancerAt2 = cdist(md_cancer_rich,md_at2_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt2) = rownames(md_cancer_rich)
	colnames(dd_CancerAt2) = rownames(md_at2_grid)
	min_dist = apply(dd_CancerAt2,2,min)
	md_at2_grid[colnames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)*grid_spacing
	dd_CancerAt1 = cdist(md_cancer_rich,md_at1_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt1) = rownames(md_cancer_rich)
	colnames(dd_CancerAt1) = rownames(md_at1_grid)
	min_dist = apply(dd_CancerAt1,2,min)
	md_at1_grid[colnames(dd_CancerAt1),"nearest_dist"] = as.numeric(min_dist)*grid_spacing

	fileName = paste0(this_OutDir,"LUAD_Heatmap_tileSizeMicrons",grid_spacing,"_inAT1_hlcaAT1scores.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_at1_grid),"value"] = md_at1_grid$hlca_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="hlca_AT1",colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),na.value='black')+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("HLCA_AT1 scores in AT1 cells, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	xe_AT12Cancer = diseased
	thiss = md
	rownames(thiss) = substr(rownames(md),2,nchar(rownames(md)))
	xe_AT12Cancer@meta.data = thiss[rownames(xe_AT12Cancer@meta.data),]

	subxe = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="AT1",]))
	subxe = nn_smoothing(subxe, feature = "hlca_AT1", radius = 50)
	md = subxe@meta.data
	dan.save(md,paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf"))

	pdf(paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf"),3,2.2,pointsize=6)
    coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT1",axes=T,size=0.3,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
    dev.off()
	# pdf(paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf" ),7,7)
	# coords = GetTissueCoordinates(subxe, which = "centroids")
	# plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT1",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	# print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	# dev.off()

	pdf(paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50_whiteBg.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="hlca_AT1",axes=T,size=0.7,border.size=NA,dark.background=F)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+theme_classic()+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()

	fileName = paste0(this_OutDir,"LUAD_Heatmap_tileSizeMicrons",grid_spacing,"_inAT2_hlcaAT2scores.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_at2_grid),"value"] = md_at2_grid$hlca_AT2
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="hlca_AT2",colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),na.value='black')+
	  # scale_fill_gradientn(name="hlca_AT2",colours=colorRampPalette(c("white",as.character(celltype_map["AT2","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("HLCA_AT2 scores in AT2 cells, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()
	subxe = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="AT2",]))
	subxe = nn_smoothing(subxe, feature = "hlca_AT2", radius = 50)
	md = subxe@meta.data
	dan.save(md,paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50.pdf"))
	pdf(paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50.pdf"),3,2.2,pointsize=6)
    coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT2",axes=T,size=0.25,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
    dev.off()
	pdf(paste0(this_OutDir,"LUAD_ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50_whiteBg.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT2",axes=T,size=0.7,border.size=NA,dark.background=F)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+theme_classic()+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()

	md_at2_grid$dist_discrete = ">1mm"
	md_at2_grid[md_at2_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_at2_grid[(md_at2_grid$nearest_dist>250) & (md_at2_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_at2_grid[(md_at2_grid$nearest_dist>500) & (md_at2_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_at2_grid[(md_at2_grid$nearest_dist>750) & (md_at2_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"
	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"LUAD_AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),2.5,2.1 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	LUAD_md_at2_grid = md_at2_grid

	md_at1_grid$dist_discrete = ">1mm"
	md_at1_grid[md_at1_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_at1_grid[(md_at1_grid$nearest_dist>250) & (md_at1_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_at1_grid[(md_at1_grid$nearest_dist>500) & (md_at1_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_at1_grid[(md_at1_grid$nearest_dist>750) & (md_at1_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"
	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"LUAD_AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),2.5,2.1 )
	for (tp in intersect(names(gs),colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	LUAD_md_at1_grid = md_at1_grid

	tdf2 = aggregate(.~dist_discrete,data=md_at2_grid[,c( "dist_discrete","n_AT1","n_AT2" )],FUN='sum')
	tdf2$AT2 = tdf2$n_AT2/(tdf2$n_AT2+tdf2$n_AT1)
	tdf2$AT1 = tdf2$n_AT1/(tdf2$n_AT2+tdf2$n_AT1)
	tdf2 = melt(tdf2[,c("dist_discrete","AT2","AT1" )])
	tdf2$variable = factor(tdf2$variable,levels=c( "AT1","AT2" ))
	tdf2$dist_discrete = factor(tdf2$dist_discrete, levels = levelz[levelz %in% unique(tdf2$dist_discrete)])
	px = ggplot(data=tdf2, aes(x=dist_discrete, y=value, fill=variable)) + geom_bar(stat="identity", colour="black", linewidth = 0.1) + labs(x = "Distance from nearest cancer-rich tile", y = "Proportion" ) + ggtitle("") + scale_fill_manual(name = "Alveolar\ncell type",values=c( "chocolate","dodgerblue2" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))
	pdf( paste0(this_OutDir,"LUAD_ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,"_BarplotProportions.pdf"),5,3 )
	print(px)
	dev.off()

	tdf2 = md_at2_grid[((md_at2_grid$n_AT1)+(md_at2_grid$n_AT2))>0,]
	tdf2$ratio_AT2_AT1 = tdf2$n_AT2/(tdf2$n_AT1+tdf2$n_AT2)
	fileName = paste0(this_OutDir,"LUAD_Heatmap_tileSizeMicrons",grid_spacing,"_AT2AT1_proportions_.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(tdf2),"value"] = tdf2$ratio_AT2_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="AT2 /\n(AT2+AT1)",colours=colorRampPalette(c("chocolate","white","dodgerblue2" ))(n = 100),na.value='black')+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT2 proportion, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	### non-diseased
	grid_spacing = 200
	newscores = scores
	nondiseased_md = cbind(nondiseased_md,newscores[rownames(nondiseased_md),names(gs)])
	md = nondiseased_md
	xmax = ceiling(max(md$x)/grid_spacing)
  	ymax = ceiling(max(md$y)/grid_spacing)
	md$xg = factor(ceiling(md$x/grid_spacing), levels = c(1:xmax))
  	md$yg = factor(ceiling(md$y/grid_spacing), levels = c(1:ymax))
  	md$grid_coo = paste0( "xg",md$xg,"_yg",md$yg )

	celltype_map = data.frame(row.names=c( "AT1","AT2","Cancer" ),colorz = c("chocolate","dodgerblue2","firebrick") ,stringsAsFactors=F )
	this_OutDir = OutDir
	# AT1
  	md_at1 = md[md$annd_level_2=="AT1",]
  	md_at1$xg = factor(ceiling(md_at1$x/grid_spacing), levels = c(1:xmax))
  	md_at1$yg = factor(ceiling(md_at1$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at1 = (table(md_at1$xg, md_at1$yg))
  	md_at1$grid_coo = paste0( "xg",md_at1$xg,"_yg",md_at1$yg )
  	md_at1_grid = aggregate( .~grid_coo,data=md_at1[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_at1_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at1_grid$grid_coo)))
  	md_at1_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at1_grid$grid_coo)))
  	rownames(md_at1_grid) = md_at1_grid$grid_coo
  	md_at1_grid$n_AT1 = table(md_at1$grid_coo)[rownames(md_at1_grid)]
  	# AT2
  	md_at2 = md[md$annd_level_2=="AT2",]
  	md_at2$xg = factor(ceiling(md_at2$x/grid_spacing), levels = c(1:xmax))
  	md_at2$yg = factor(ceiling(md_at2$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at2 = (table(md_at2$xg, md_at2$yg))
  	md_at2$grid_coo = paste0( "xg",md_at2$xg,"_yg",md_at2$yg )
  	md_at2_grid = aggregate( .~grid_coo,data=md_at2[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_at2_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at2_grid$grid_coo)))
  	md_at2_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at2_grid$grid_coo)))
  	rownames(md_at2_grid) = md_at2_grid$grid_coo
  	md_at2_grid$n_AT2 = table(md_at2$grid_coo)[rownames(md_at2_grid)]

	md_at1_grid$n_AT2 = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"]
	md_at2_grid$n_AT1 = 0
	md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"]

	fileName = paste0(this_OutDir,"nondiseased_Heatmap_tileSizeMicrons",grid_spacing,"_inAT1_hlcaAT1scores.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_at1_grid),"value"] = md_at1_grid$hlca_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="hlca_AT1",colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),na.value='black')+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("HLCA_AT1 scores in AT1 cells, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	xe_AT12Cancer = nondiseased
	thiss = md
	rownames(thiss) = substr(rownames(md),2,nchar(rownames(md)))
	xe_AT12Cancer@meta.data = thiss[rownames(xe_AT12Cancer@meta.data),]

	subxe = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="AT1",]))
	subxe = nn_smoothing(subxe, feature = "hlca_AT1", radius = 50)
	pdf(paste0(this_OutDir,"nondiseased_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT1",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
	pdf(paste0(this_OutDir,"nondiseased_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50_whiteBG.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT1",axes=T,size=0.7,border.size=NA, dark.background=F)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+theme_classic()+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
	# pdf(paste0(OutDir,"LUAD_ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf" ),2.1,2.5)
	# 	coords = GetTissueCoordinates(subxe, which = "centroids")
	# 	plot=ImageFeaturePlot(subxe,features=paste0("smoothed_hlca_AT1"),axes=T,size=0.05,border.size=NA,max.cutoff='q99',min.cutoff='q1')+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
	# 	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
	# 	plot = plot + theme(legend.position = "bottom", legend.direction = "horizontal")
	# 	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	# 	dev.off()
	fileName = paste0(this_OutDir,"nondiseased_Heatmap_tileSizeMicrons",grid_spacing,"_inAT2_hlcaAT2scores.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_at2_grid),"value"] = md_at2_grid$hlca_AT2
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="hlca_AT2",colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),na.value='black')+
	  # scale_fill_gradientn(name="hlca_AT2",colours=colorRampPalette(c("white",as.character(celltype_map["AT2","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("HLCA_AT2 scores in AT2 cells, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()
	subxe = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="AT2",]))
	subxe = nn_smoothing(subxe, feature = "hlca_AT2", radius = 50)
	pdf(paste0(this_OutDir,"nondiseased_ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT2",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
	pdf(paste0(this_OutDir,"nondiseased_ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50_whiteBG.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT2",axes=T,size=0.7,border.size=NA,dark.background=F)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+theme_classic()+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()

	tdf2 = md_at2_grid[((md_at2_grid$n_AT1)+(md_at2_grid$n_AT2))>0,]
	tdf2$ratio_AT2_AT1 = tdf2$n_AT2/(tdf2$n_AT1+tdf2$n_AT2)
	fileName = paste0(this_OutDir,"nondiseased_Heatmap_tileSizeMicrons",grid_spacing,"_AT2AT1_proportions_.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(tdf2),"value"] = tdf2$ratio_AT2_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="AT2 /\n(AT2+AT1)",colours=colorRampPalette(c("chocolate","white","dodgerblue2" ))(n = 100),na.value='black')+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT2 proportion, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	### nondiseased vs LUAD
	md_at2_grid$dist_discrete = "healthy lung"
	commonz = intersect(colnames(LUAD_md_at2_grid),colnames(md_at2_grid))
	md_at2_grid = rbind(LUAD_md_at2_grid[,commonz],md_at2_grid[,commonz])
	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm","healthy lung" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33","gray11" )
	pdf( paste0(this_OutDir,"compare_AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),2.5,2.1 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	dan.save(md_at2_grid,paste0(this_OutDir,"compare_AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"))

	md_at1_grid$dist_discrete = "healthy lung"
	commonz = intersect(colnames(LUAD_md_at1_grid),colnames(md_at1_grid))
	md_at1_grid = rbind(LUAD_md_at1_grid[,commonz],md_at1_grid[,commonz])
	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm","healthy lung" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33","gray11" )
	pdf( paste0(this_OutDir,"compare_AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),2.5,2.1 )
	for (tp in intersect(names(gs),colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	dan.save(md_at1_grid,paste0(this_OutDir,"compare_AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"))

	tdf2 = aggregate(.~dist_discrete,data=md_at2_grid[,c( "dist_discrete","n_AT1","n_AT2" )],FUN='sum')
	tdf2$AT2 = tdf2$n_AT2/(tdf2$n_AT2+tdf2$n_AT1)
	tdf2$AT1 = tdf2$n_AT1/(tdf2$n_AT2+tdf2$n_AT1)
	tdf2 = melt(tdf2[,c("dist_discrete","AT2","AT1" )])
	tdf2$variable = factor(tdf2$variable,levels=c( "AT1","AT2" ))
	tdf2$dist_discrete = factor(tdf2$dist_discrete, levels = levelz[levelz %in% unique(tdf2$dist_discrete)])
	px = ggplot(data=tdf2, aes(x=dist_discrete, y=value, fill=variable)) + geom_bar(stat="identity", colour="black", linewidth = 0.1) + labs(x = "Distance from nearest cancer-rich tile", y = "Proportion" ) + ggtitle("") + scale_fill_manual(name = "Alveolar\ncell type",values=c( "chocolate","dodgerblue2" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))
	pdf( paste0(this_OutDir,"compare_ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,"_BarplotProportions.pdf"),5,3 )
	print(px)
	dev.off()
}

nn_smoothing = function( subxe, feature, radius = 50 ){
	library(dbscan)
	md = subxe@meta.data
	coords = GetTissueCoordinates(subxe, which = "centroids")
	coords$cell = NULL
	nn = frNN(x=coords,eps=radius)
	nnid = nn$id
	md[,paste0("smoothed","_",feature)] = sapply( 1:nrow(md), function(x) mean(md[c(x,nnid[[x]]),feature]) )
	subxe@meta.data = md
	return( subxe )
}

xenium_gridding = function( OutDir, whichDataset, grid_spacing = 200, cancer_rich_ratio = 0.75, density_threshold = 0.002 ){
	library(RColorBrewer)
	load(file=paste0(OutDir,"gs.RData"))
	celltype_map = data.frame(row.names=c( "AT1","AT2","Cancer" ),colorz = c("chocolate","dodgerblue2","firebrick") ,stringsAsFactors=F )
	this_OutDir = paste0(OutDir,"gridding_",grid_spacing,"/")
	dir.create(this_OutDir)
	load(file=paste0(OutDir,"../xenium_preprocessing_",whichDataset,"/xe_AT12Cancer.RData" ))
	load(file = paste0( OutDir,"md_allCells_allValues.RData" ))
	xmax = ceiling(max(md$x)/grid_spacing)
  	ymax = ceiling(max(md$y)/grid_spacing)
  	# AT1
  	md_at1 = md[md$annd_level_2=="AT1",]
  	md_at1$xg = factor(ceiling(md_at1$x/grid_spacing), levels = c(1:xmax))
  	md_at1$yg = factor(ceiling(md_at1$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at1 = (table(md_at1$xg, md_at1$yg))
  	md_at1$grid_coo = paste0( "xg",md_at1$xg,"_yg",md_at1$yg )
  	md_at1_grid = aggregate( .~grid_coo,data=md_at1[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_at1_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at1_grid$grid_coo)))
  	md_at1_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at1_grid$grid_coo)))
  	rownames(md_at1_grid) = md_at1_grid$grid_coo
  	md_at1_grid$n_AT1 = table(md_at1$grid_coo)[rownames(md_at1_grid)]
  	# AT2
  	md_at2 = md[md$annd_level_2=="AT2",]
  	md_at2$xg = factor(ceiling(md_at2$x/grid_spacing), levels = c(1:xmax))
  	md_at2$yg = factor(ceiling(md_at2$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at2 = (table(md_at2$xg, md_at2$yg))
  	md_at2$grid_coo = paste0( "xg",md_at2$xg,"_yg",md_at2$yg )
  	md_at2_grid = aggregate( .~grid_coo,data=md_at2[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_at2_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at2_grid$grid_coo)))
  	md_at2_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at2_grid$grid_coo)))
  	rownames(md_at2_grid) = md_at2_grid$grid_coo
  	md_at2_grid$n_AT2 = table(md_at2$grid_coo)[rownames(md_at2_grid)]
  	# Cancer
  	md_cancer = md[md$annd_level_2=="Cancer",]
  	md_cancer$xg = factor(ceiling(md_cancer$x/grid_spacing), levels = c(1:xmax))
  	md_cancer$yg = factor(ceiling(md_cancer$y/grid_spacing), levels = c(1:ymax))
  	meaMat_cancer = (table(md_cancer$xg, md_cancer$yg))
  	md_cancer$grid_coo = paste0( "xg",md_cancer$xg,"_yg",md_cancer$yg )
  	md_cancer_grid = aggregate( .~grid_coo,data=md_cancer[,c( "grid_coo",names(gs) )],FUN='mean' )
  	md_cancer_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_cancer_grid$grid_coo)))
  	md_cancer_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_cancer_grid$grid_coo)))
  	rownames(md_cancer_grid) = md_cancer_grid$grid_coo
  	md_cancer_grid$n_Cancer = table(md_cancer$grid_coo)[rownames(md_cancer_grid)]

	md_at1_grid$n_AT2 = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"]
	md_at1_grid$n_Cancer = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_Cancer"] = md_cancer_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_Cancer"]
	md_at2_grid$n_AT1 = 0
	md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"]
	md_at2_grid$n_Cancer = 0
	md_at2_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_Cancer"] = md_cancer_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_Cancer"]
	md_cancer_grid$n_AT1 = 0
	md_cancer_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_AT1"]
	md_cancer_grid$n_AT2 = 0
	md_cancer_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_AT2"]
	md_cancer_grid$ratio_Cancer = md_cancer_grid$n_Cancer/(md_cancer_grid$n_Cancer+md_cancer_grid$n_AT1+md_cancer_grid$n_AT2)
	md_cancer_grid$spatial_density_Cancer = md_cancer_grid$n_Cancer/(grid_spacing*grid_spacing)

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_cancer_spatialdensity.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_cancer)
	dfmelt$value = 0
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_cancer_grid),"value"] = md_cancer_grid$spatial_density_Cancer
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="Density",colours=colorRampPalette(c("white",as.character(celltype_map["Cancer","colorz"]) ))(n = 100))+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("Density (N/µm²) of cancer cells, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	md_cancer_grid$cancer_rich = 0
	md_cancer_grid[(md_cancer_grid$spatial_density_Cancer>density_threshold),"cancer_rich"] = 1
	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_is_cancer_rich_onlyDensity.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_cancer)
	dfmelt$value = 0
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_cancer_grid),"value"] = md_cancer_grid$cancer_rich
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="Rich",colours=colorRampPalette(c("white",as.character(celltype_map["Cancer","colorz"]) ))(n = 100))+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("Is cancer-rich, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	fileName = paste0(this_OutDir,"densityPlot_cancer_spatialdensity.pdf")
	pdf(fileName,5,3)
	plot = ggplot(md_cancer_grid, aes(x=spatial_density_Cancer)) + geom_density() + geom_vline(xintercept=density_threshold, linetype="dashed", color = "red") + labs(x="Cancer cell density (N/(µm²))",y="Density") + theme_classic()
	print(plot)
	dev.off()
	dan.save(md_cancer_grid,fileName)

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_AT1_counts.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at1)
	dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="n_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT1 cells, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()
	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_AT2_counts.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="n_AT2",colours=colorRampPalette(c("white",as.character(celltype_map["AT2","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT2 cells, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()
	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_Cancer_counts.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_cancer)
	dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="n_Cancer",colours=colorRampPalette(c("white",as.character(celltype_map["Cancer","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("Cancer cells, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+ coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_AT2AT1_ratios.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	meaMat_proportions = log10(meaMat_at2/meaMat_at1)
	meaMat_proportions[(meaMat_at1==0)] = max(meaMat_proportions[meaMat_proportions!=Inf],na.rm=T)
	meaMat_proportions[(meaMat_at2==0)] = min(meaMat_proportions[meaMat_proportions!=(-Inf)],na.rm=T)
	meaMat_proportions[(meaMat_at2==0) & (meaMat_at2==0)] = NA
	dfmelt = melt(meaMat_proportions)
	dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="log10(AT2/AT1)",colours=colorRampPalette(c("white","sienna4" ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT2/AT1 ratios, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	# md_cancer_rich = md_cancer_grid[(md_cancer_grid$ratio_Cancer>cancer_rich_ratio) & (md_cancer_grid$spatial_density_Cancer>density_threshold),c( "xg","yg" )]
	md_cancer_rich = md_cancer_grid[(md_cancer_grid$spatial_density_Cancer>density_threshold),c( "xg","yg" )]

	library(rdist)
	dd_CancerAt2 = cdist(md_cancer_rich,md_at2_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt2) = rownames(md_cancer_rich)
	colnames(dd_CancerAt2) = rownames(md_at2_grid)
	min_dist = apply(dd_CancerAt2,2,min)
	md_at2_grid[colnames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)*grid_spacing
	dd_CancerAt1 = cdist(md_cancer_rich,md_at1_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt1) = rownames(md_cancer_rich)
	colnames(dd_CancerAt1) = rownames(md_at1_grid)
	min_dist = apply(dd_CancerAt1,2,min)
	md_at1_grid[colnames(dd_CancerAt1),"nearest_dist"] = as.numeric(min_dist)*grid_spacing

	pdf( paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,".pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at1_grid))){
	   x = md_at1_grid$nearest_dist
	   y = md_at1_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,".pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$nearest_dist
	   y = md_at2_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	tdf = md_at2_grid[md_at2_grid$n_AT1>0,]
	tdf$ratio_AT2_AT1 = log10(tdf$n_AT2/tdf$n_AT1)
	pdf( paste0(this_OutDir,"ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,".pdf"),6,5 )
	for (tp in c( "ratio_AT2_AT1" ) ){
	   x = tdf$nearest_dist
	   y = tdf[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = "log10(n_AT2 / n_AT1)", plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	pdf( paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_loess.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at1_grid))){
	   x = md_at1_grid$nearest_dist
	   y = md_at1_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "loess", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_loess.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$nearest_dist
	   y = md_at2_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "loess", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	tdf = md_at2_grid[md_at2_grid$n_AT1>0,]
	tdf$ratio_AT2_AT1 = log10(tdf$n_AT2/tdf$n_AT1)
	pdf( paste0(this_OutDir,"ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,"_loess.pdf"),6,5 )
	for (tp in c( "ratio_AT2_AT1" ) ){
	   x = tdf$nearest_dist
	   y = tdf[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = "log10(n_AT2 / n_AT1)", plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "loess", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_AT2AT1_ratios_filtered.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(tdf),"value"] = tdf$ratio_AT2_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="log10(AT2/AT1)",colours=colorRampPalette(c("white","sienna4" ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT2/AT1 ratios, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_inAT1_hlcaAT1scores.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_at1_grid),"value"] = md_at1_grid$hlca_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="hlca_AT1",colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),na.value='black')+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("HLCA_AT1 scores in AT1 cells, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	xe_AT12Cancer@meta.data = md[rownames(xe_AT12Cancer@meta.data),]

	subxe = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="AT1",]))
	pdf(paste0(this_OutDir,"ImageFeaturePlot_inAT1_hlcaAT1scores.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="hlca_AT1",axes=T,border.size=NA,size=0.7)+labs(y="x (µm)",x="y (µm)")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
	subxe = nn_smoothing(subxe, feature = "hlca_AT1", radius = 50)
	pdf(paste0(this_OutDir,"ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT1",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
	pdf(paste0(OutDir,"ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed50.pdf" ),2.1,2.5)
		coords = GetTissueCoordinates(subxe, which = "centroids")
		plot=ImageFeaturePlot(subxe,features=paste0("smoothed_hlca_AT1"),axes=T,size=0.2,border.size=NA,max.cutoff='q99',min.cutoff='q1')+labs(y="x (µm)",x="y (µm)",fill=cs)+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
		plot = plot + theme(legend.position = "bottom", legend.direction = "horizontal")
		print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		dev.off()
	subxe = nn_smoothing(subxe, feature = "hlca_AT1", radius = 100)
	pdf(paste0(this_OutDir,"ImageFeaturePlot_inAT1_hlcaAT1scores_smoothed100.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT1",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT1")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_inAT2_hlcaAT2scores.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(md_at2_grid),"value"] = md_at2_grid$hlca_AT2
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="hlca_AT2",colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)),na.value='black')+
	  # scale_fill_gradientn(name="hlca_AT2",colours=colorRampPalette(c("white",as.character(celltype_map["AT2","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("HLCA_AT2 scores in AT2 cells, grid spacing = ",grid_spacing,"µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()

	subxe = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="AT2",]))
	pdf(paste0(this_OutDir,"ImageFeaturePlot_inAT2_hlcaAT2scores.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="hlca_AT2",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
	subxe = nn_smoothing(subxe, feature = "hlca_AT2", radius = 50)
	pdf(paste0(this_OutDir,"ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT2",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()
		pdf(paste0(OutDir,"ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed50.pdf" ),2.1,2.5)
		coords = GetTissueCoordinates(subxe, which = "centroids")
		plot=ImageFeaturePlot(subxe,features=paste0("smoothed_hlca_AT2"),axes=T,size=0.15,border.size=NA,max.cutoff='q99',min.cutoff='q1')+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
		plot = plot + theme(legend.position = "bottom", legend.direction = "horizontal")
		print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		dev.off()
	subxe = nn_smoothing(subxe, feature = "hlca_AT2", radius = 100)
	pdf(paste0(this_OutDir,"ImageFeaturePlot_inAT2_hlcaAT2scores_smoothed100.pdf" ),7,7)
	coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageFeaturePlot(subxe,features="smoothed_hlca_AT2",axes=T,size=0.7,border.size=NA)+labs(y="x (µm)",x="y (µm)",fill="hlca_AT2")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
	dev.off()

	md_at2_grid$dist_discrete = ">1mm"
	md_at2_grid[md_at2_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_at2_grid[(md_at2_grid$nearest_dist>250) & (md_at2_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_at2_grid[(md_at2_grid$nearest_dist>500) & (md_at2_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_at2_grid[(md_at2_grid$nearest_dist>750) & (md_at2_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"
	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	pdf( paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_ViolinPlots.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   plot=dan.violinplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo=max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	dan.save(md_at2_grid,paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"))
	# density plots
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   dan.densityPlot( paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_DensityPlot_",tp,".pdf"), y, x, groupinglab = "", xlab = paste0(tp," score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = these_colorz, fileWidth = 6, fileHeight = 3 )
	}

	md_at1_grid$dist_discrete = ">1mm"
	md_at1_grid[md_at1_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_at1_grid[(md_at1_grid$nearest_dist>250) & (md_at1_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_at1_grid[(md_at1_grid$nearest_dist>500) & (md_at1_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_at1_grid[(md_at1_grid$nearest_dist>750) & (md_at1_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"

	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	pdf( paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_ViolinPlots.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   plot=dan.violinplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal",labelycoo=max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	dan.save(md_at1_grid,paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"))
	# density plots
	for (tp in intersect(names(gs),colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   dan.densityPlot( paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_DensityPlot_",tp,".pdf"), y, x, groupinglab = "", xlab = paste0(tp," score"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = these_colorz, fileWidth = 6, fileHeight = 3 )
	}

	tdf$dist_discrete = ">1mm"
	tdf[tdf$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	tdf[(tdf$nearest_dist>250) & (tdf$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	tdf[(tdf$nearest_dist>500) & (tdf$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	tdf[(tdf$nearest_dist>750) & (tdf$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"

	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),6,5 )
	for (tp in c( "ratio_AT2_AT1" ) ){
	   x = tdf$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = tdf[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	tdf2 = aggregate(.~dist_discrete,data=md_at2_grid[,c( "dist_discrete","n_AT1","n_AT2" )],FUN='sum')
	tdf2$AT2 = tdf2$n_AT2/(tdf2$n_AT2+tdf2$n_AT1)
	tdf2$AT1 = tdf2$n_AT1/(tdf2$n_AT2+tdf2$n_AT1)
	tdf2 = melt(tdf2[,c("dist_discrete","AT2","AT1" )])
	tdf2$variable = factor(tdf2$variable,levels=c( "AT1","AT2" ))
	tdf2$dist_discrete = factor(tdf2$dist_discrete, levels = levelz[levelz %in% unique(tdf2$dist_discrete)])
	px = ggplot(data=tdf2, aes(x=dist_discrete, y=value, fill=variable)) + geom_bar(stat="identity", colour="black", linewidth = 0.1) + labs(x = "Distance from nearest cancer-rich tile", y = "Proportion" ) + ggtitle("") + scale_fill_manual(name = "Alveolar\ncell type",values=c( "chocolate","dodgerblue2" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))
	pdf( paste0(this_OutDir,"ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,"_BarplotProportions.pdf"),5,3 )
	print(px)
	dev.off()

	tdf2 = md_at2_grid[((md_at2_grid$n_AT1)+(md_at2_grid$n_AT2))>0,]
	tdf2$ratio_AT2_AT1 = tdf2$n_AT2/(tdf2$n_AT1+tdf2$n_AT2)
	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_AT2AT1_proportions_.pdf")
	pdf(fileName, 6, 6, useDingbats = F )
	dfmelt = melt(meaMat_at2)
	dfmelt$value = NA
	rownames(dfmelt) = paste0( "xg",dfmelt$Var1,"_yg",dfmelt$Var2 )
	dfmelt[rownames(tdf2),"value"] = tdf2$ratio_AT2_AT1
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="AT2 /\n(AT2+AT1)",colours=colorRampPalette(c("chocolate","white","dodgerblue2" ))(n = 100),na.value='black')+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("AT2 proportion, grid spacing = ",grid_spacing," µm"), x="x (grid)",y="y (grid)")+
	  coord_fixed() + scale_y_discrete(limits=rev)
	print(p)
	dev.off()



	library(viridis)

	md_at2_grid$dist_discrete = ">=1mm"
	md_at2_grid[md_at2_grid$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_at2_grid[(md_at2_grid$nearest_dist>=lowerbound) & (md_at2_grid$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(this_OutDir,"AT2_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots_SameGridSpacing.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   jitterColorz = dan.expand_colors(x,levels(x),these_colorz)
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = T,jitterColors = jitterColorz)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	md_at1_grid$dist_discrete = ">=1mm"
	md_at1_grid[md_at1_grid$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_at1_grid[(md_at1_grid$nearest_dist>=lowerbound) & (md_at1_grid$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(this_OutDir,"AT1_alveolar_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots_SameGridSpacing.pdf"),6,5 )
	for (tp in intersect(names(gs),colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	tdf$dist_discrete = ">=1mm"
	tdf[tdf$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		tdf[(tdf$nearest_dist>=lowerbound) & (tdf$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(this_OutDir,"ratio_AT2_AT1_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots_SameGridSpacing.pdf"),6,5 )
	for (tp in c( "ratio_AT2_AT1" ) ){
	   x = tdf$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = tdf[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
}

xenium_gridding_singlegene = function( OutDir, gene, whichDataset, grid_spacing = 200, cancer_rich_ratio = 0.75, density_threshold = 0.002 ){
	library(RColorBrewer)
	celltype_map = data.frame(row.names=c( "AT1","AT2","Cancer" ),colorz = c("chocolate","dodgerblue2","firebrick") ,stringsAsFactors=F )	
	this_OutDir = paste0(OutDir,"singlegene_gridding_",gene,"_",grid_spacing,"/")
	dir.create(this_OutDir)
	load(file=paste0(OutDir,"../xenium_preprocessing_",whichDataset,"/xe_AT12Cancer.RData" ))
	load(file = paste0( OutDir,"md_allCells_allValues.RData" ))
	xmax = ceiling(max(md$x)/grid_spacing)
  	ymax = ceiling(max(md$y)/grid_spacing)
  	md[,gene] = as.numeric(xe_AT12Cancer@assays$SCT@data[gene,rownames(md)])
  	# AT1
  	md_at1 = md[md$annd_level_2=="AT1",]
  	md_at1$xg = factor(ceiling(md_at1$x/grid_spacing), levels = c(1:xmax))
  	md_at1$yg = factor(ceiling(md_at1$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at1 = t(table(md_at1$xg, md_at1$yg))
  	md_at1$grid_coo = paste0( "xg",md_at1$xg,"_yg",md_at1$yg )
  	md_at1_grid = aggregate( .~grid_coo,data=md_at1[,c( "grid_coo",gene )],FUN='mean' )
  	md_at1_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at1_grid$grid_coo)))
  	md_at1_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at1_grid$grid_coo)))
  	rownames(md_at1_grid) = md_at1_grid$grid_coo
  	md_at1_grid$n_AT1 = table(md_at1$grid_coo)[rownames(md_at1_grid)]
  	# AT2
  	md_at2 = md[md$annd_level_2=="AT2",]
  	md_at2$xg = factor(ceiling(md_at2$x/grid_spacing), levels = c(1:xmax))
  	md_at2$yg = factor(ceiling(md_at2$y/grid_spacing), levels = c(1:ymax))
  	meaMat_at2 = t(table(md_at2$xg, md_at2$yg))
  	md_at2$grid_coo = paste0( "xg",md_at2$xg,"_yg",md_at2$yg )
  	md_at2_grid = aggregate( .~grid_coo,data=md_at2[,c( "grid_coo",gene )],FUN='mean' )
  	md_at2_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_at2_grid$grid_coo)))
  	md_at2_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_at2_grid$grid_coo)))
  	rownames(md_at2_grid) = md_at2_grid$grid_coo
  	md_at2_grid$n_AT2 = table(md_at2$grid_coo)[rownames(md_at2_grid)]
  	# Cancer
  	md_cancer = md[md$annd_level_2=="Cancer",]
  	md_cancer$xg = factor(ceiling(md_cancer$x/grid_spacing), levels = c(1:xmax))
  	md_cancer$yg = factor(ceiling(md_cancer$y/grid_spacing), levels = c(1:ymax))
  	meaMat_cancer = t(table(md_cancer$xg, md_cancer$yg))
  	md_cancer$grid_coo = paste0( "xg",md_cancer$xg,"_yg",md_cancer$yg )
  	md_cancer_grid = aggregate( .~grid_coo,data=md_cancer[,c( "grid_coo",gene )],FUN='mean' )
  	md_cancer_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_cancer_grid$grid_coo)))
  	md_cancer_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_cancer_grid$grid_coo)))
  	rownames(md_cancer_grid) = md_cancer_grid$grid_coo
  	md_cancer_grid$n_Cancer = table(md_cancer$grid_coo)[rownames(md_cancer_grid)]

	md_at1_grid$n_AT2 = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT2"]
	md_at1_grid$n_Cancer = 0
	md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_Cancer"] = md_cancer_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_Cancer"]
	md_at2_grid$n_AT1 = 0
	md_at2_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_at2_grid)),"n_AT1"]
	md_at2_grid$n_Cancer = 0
	md_at2_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_Cancer"] = md_cancer_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_Cancer"]
	md_cancer_grid$n_AT1 = 0
	md_cancer_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_AT1"] = md_at1_grid[intersect(rownames(md_at1_grid),rownames(md_cancer_grid)),"n_AT1"]
	md_cancer_grid$n_AT2 = 0
	md_cancer_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_AT2"] = md_at2_grid[intersect(rownames(md_at2_grid),rownames(md_cancer_grid)),"n_AT2"]
	md_cancer_grid$ratio_Cancer = md_cancer_grid$n_Cancer/(md_cancer_grid$n_Cancer+md_cancer_grid$n_AT1+md_cancer_grid$n_AT2)
	md_cancer_grid$spatial_density_Cancer = md_cancer_grid$n_Cancer/(grid_spacing*grid_spacing)

	# md_cancer_rich = md_cancer_grid[(md_cancer_grid$ratio_Cancer>cancer_rich_ratio) & (md_cancer_grid$spatial_density_Cancer>density_threshold),c( "xg","yg" )]
	md_cancer_rich = md_cancer_grid[(md_cancer_grid$spatial_density_Cancer>density_threshold),c( "xg","yg" )]

	library(rdist)
	dd_CancerAt2 = cdist(md_cancer_rich,md_at2_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt2) = rownames(md_cancer_rich)
	colnames(dd_CancerAt2) = rownames(md_at2_grid)
	min_dist = apply(dd_CancerAt2,2,min)
	md_at2_grid[colnames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)*grid_spacing
	dd_CancerAt1 = cdist(md_cancer_rich,md_at1_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt1) = rownames(md_cancer_rich)
	colnames(dd_CancerAt1) = rownames(md_at1_grid)
	min_dist = apply(dd_CancerAt1,2,min)
	md_at1_grid[colnames(dd_CancerAt1),"nearest_dist"] = as.numeric(min_dist)*grid_spacing

	pdf( paste0(this_OutDir,"AT1_expression",gene,"_vsDistance_tileSizeMicrons",grid_spacing,".pdf"),6,5 )
	for (tp in intersect(gene,colnames(md_at1_grid))){
	   x = md_at1_grid$nearest_dist
	   y = md_at1_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()
	pdf( paste0(this_OutDir,"AT2_expression",gene,"_vsDistance_tileSizeMicrons",grid_spacing,".pdf"),6,5 )
	for (tp in intersect(gene,colnames(md_at2_grid)) ){
	   x = md_at2_grid$nearest_dist
	   y = md_at2_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman"),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	md_at2_grid$dist_discrete = ">1mm"
	md_at2_grid[md_at2_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_at2_grid[(md_at2_grid$nearest_dist>250) & (md_at2_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_at2_grid[(md_at2_grid$nearest_dist>500) & (md_at2_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_at2_grid[(md_at2_grid$nearest_dist>750) & (md_at2_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"

	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"AT2_expression",gene,"_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),6,5 )
	for (tp in intersect(gene,colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	md_at1_grid$dist_discrete = ">1mm"
	md_at1_grid[md_at1_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_at1_grid[(md_at1_grid$nearest_dist>250) & (md_at1_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_at1_grid[(md_at1_grid$nearest_dist>500) & (md_at1_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_at1_grid[(md_at1_grid$nearest_dist>750) & (md_at1_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"

	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"AT1_expression",gene,"_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),6,5 )
	for (tp in intersect(gene,colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	library(viridis)

	md_at2_grid$dist_discrete = ">=1mm"
	md_at2_grid[md_at2_grid$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_at2_grid[(md_at2_grid$nearest_dist>=lowerbound) & (md_at2_grid$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(this_OutDir,"AT2_expression",gene,"_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots_SameGridSpacing.pdf"),6,5 )
	for (tp in intersect(gene,colnames(md_at2_grid)) ){
	   x = md_at2_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at2_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	md_at1_grid$dist_discrete = ">=1mm"
	md_at1_grid[md_at1_grid$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_at1_grid[(md_at1_grid$nearest_dist>=lowerbound) & (md_at1_grid$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(this_OutDir,"AT1_expression",gene,"_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots_SameGridSpacing.pdf"),6,5 )
	for (tp in intersect(gene,colnames(md_at1_grid)) ){
	   x = md_at1_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md_at1_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
}

xenium_cancercells = function( OutDir,whichDataset ){
	load(file=paste0(OutDir,"../../xenium_preprocessing_",whichDataset,"/xe_AT12Cancer.RData" ))
	# scoring of cell state signatures
	md = xe_AT12Cancer@meta.data
	coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	rownames(coords) = coords$cell
	md$x = coords[rownames(md),"x"]
	md$y = coords[rownames(md),"y"]
	library(rdist)
	md$nearest_dist = NA
	# mdtest = md[sample(rownames(md),1000),]
	md_cancer = md[md$annd_level_2=="Cancer",c( "x","y" )]
	md_at2 = md[md$annd_level_2 %in% c("AT2","AT1"),c( "x","y" )]
	dd_CancerAt2 = cdist(md_cancer,md_at2)
	rownames(dd_CancerAt2) = rownames(md_cancer)
	colnames(dd_CancerAt2) = rownames(md_at2)
	min_dist = apply(dd_CancerAt2,1,min)
	md[rownames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)
	mdraw = md
	xe_AT12Cancer@meta.data = mdraw
	load( file=paste0(OutDir,"../../CellStates/signatures/wilcox_cs_signatures_restrictive.RData" ) )
	xe_AT12Cancer = subset(xe_AT12Cancer,cells=rownames(xe_AT12Cancer@meta.data[xe_AT12Cancer@meta.data$annd_level_2=="Cancer",]) )
	# xe_AT12Cancer = AddModuleScore(xe_AT12Cancer,cs_signatures,pool = NULL,nbin = 10,ctrl = 20,k = FALSE,assay = 'SCT',name = "ams")
	# colnames(xe_AT12Cancer@meta.data)[substr(colnames(xe_AT12Cancer@meta.data),1,3)=="ams"] = names(cs_signatures)
	dcat( "AmsCentered" )
	db_rand = dan.barkley_MakeRand_data(xe_AT12Cancer,cs_signatures, 3)
	scores = dan.Barkley_GeneToEnrichment_AmsCentered( xe_AT12Cancer, cs_signatures, db_rand)
	xe_AT12Cancer@meta.data = scores
	md = xe_AT12Cancer@meta.data
	md$CellState = names(cs_signatures)[apply(md[,names(cs_signatures)],1,which.max)]
	# load(file = paste0( OutDir,"md_CancerCells.RData" ))

	xe_AT12Cancer@meta.data = md
	ct_colorz = c("Alveolar"="dodgerblue","Proliferative"="red","Hypoxic"="purple")
	# pdf(paste0(OutDir,"CancerCells_imageplot_CellStates_restrictive.pdf" ),12,12)
	# print(ImageDimPlot(xe_AT12Cancer,group.by="CellState",axes=T,cols=ct_colorz,border.size=NA))
	# dev.off()
	# pdf(paste0(OutDir,"CancerCells_imageplot_CellStates_restrictive.pdf"),4,2.5,pointsize=6)
    # coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	# plot=ImageDimPlot(xe_AT12Cancer,group.by="CellState",axes=T,size=0.07,border.size=NA,cols=ct_colorz)+labs(y="x (µm)",x="y (µm)",fill="Cell state")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	# plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
    # print(plot)
    # dev.off()
    xe_AT12Cancer = RunPCA(xe_AT12Cancer, npcs = 30, features = rownames(xe_AT12Cancer))
	xe_AT12Cancer = RunUMAP(xe_AT12Cancer, dims = 1:30)

	coords = GetTissueCoordinates(xe_AT12Cancer, which = "centroids")
	cells_zoomin = coords[((coords$x<9000) & (coords$x>500)) & ((coords$y<8000)),'cell' ]
    subxe = subset(xe_AT12Cancer,cells=cells_zoomin)

    ct_colorz = c("Alveolar"="#2B3990","Proliferative"="#FBB040","Hypoxic"="#BE1E2D")
    pdf(paste0(OutDir,"CancerCells_imageplot_CellStates_restrictive_zoomin.pdf"),3.2,2,pointsize=6)
    coords = GetTissueCoordinates(subxe, which = "centroids")
	plot=ImageDimPlot(subxe,group.by="CellState",axes=T,size=0.08,border.size=NA,cols=ct_colorz)+labs(y="x (µm)",x="y (µm)",fill="Cell state")+coord_flip() +scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
    print(plot)
    dev.off()

    for (cs in names(ct_colorz)){
    	subxe = nn_smoothing(subxe, feature = cs, radius = 50)
		# pdf(paste0(OutDir,"CancerCells_imagefeatureplot_CellStates_restrictive_",cs,".pdf" ),2.1,2.5)
		# coords = GetTissueCoordinates(subxe, which = "centroids")
		# plot=ImageFeaturePlot(subxe,features=cs,axes=T,size=0.05,border.size=NA)+NoLegend()+labs(y="x (µm)",x="y (µm)")+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )
		# plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(4, 'pt'))
		# print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		# dev.off()
    	# subxe
		pdf(paste0(OutDir,"CancerCells_imagefeatureplot_CellStates_restrictive_",cs,"_smoothed50_unrastered.pdf" ),1.8,1.8)
		coords = GetTissueCoordinates(subxe, which = "centroids")
		plot=ImageFeaturePlot(subxe,features=paste0("smoothed_",cs),axes=T,size=0.05,border.size=NA,max.cutoff='q99',min.cutoff='q1')+labs(y="x (µm)",x="y (µm)",fill=cs)+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
		# plot = plot + theme(legend.position = "bottom", legend.direction = "horizontal")
		print(plot+NoLegend() & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		dev.off()

		pdf(paste0(OutDir,"CancerCells_imagefeatureplot_CellStates_restrictive_",cs,"_smoothed50_unrastered_WithLegend.pdf" ),2.5,2.1)
		coords = GetTissueCoordinates(subxe, which = "centroids")
		plot=ImageFeaturePlot(subxe,features=paste0("smoothed_",cs),axes=T,size=0.05,border.size=NA,max.cutoff='q99',min.cutoff='q1')+labs(y="x (µm)",x="y (µm)",fill=cs)+coord_flip()+scale_x_reverse() + theme(aspect.ratio = (max(coords$y)-min(coords$y))/(max(coords$x)-min(coords$x)) )+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		plot = plot + theme(axis.line.x = element_line(size = 0.2),axis.line.y = element_line(size = 0.1), text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
		# plot = plot + theme(legend.position = "bottom", legend.direction = "horizontal")
		print(plot & scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(11, 'RdYlBu'))(100))))
		dev.off()
    }
    

	save(md,file = paste0( OutDir,"md_CancerCells.RData" ))
	load(file = paste0( OutDir,"md_CancerCells.RData" ))

	levelz = c("Alveolar","Proliferative","Hypoxic")
	colorz = ct_colorz
	pdf( paste0(OutDir,"distance_normalAlveolar_byCellState.pdf"),2,2 )
	for (tp in c( "nearest_dist" ) ){
	   x = md$CellState
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md[,tp]
	   # ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   # plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0("Distance from nearest \nnon-malignant alveolar cell"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0("Distance from nearest alveolar cell (µm)"), plotTitle = "",signifTest = "kruskal", labelycoo = ylimRight, ylimRight = ylimRight, ylimLeft=0, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	pdf( paste0(OutDir,"distance_normalAlveolar_byCellState_ViolinPlots.pdf"),6,4 )
	for (tp in c( "nearest_dist" ) ){
	   x = md$CellState
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md[,tp]
	   plot=dan.violinplots.multipages( x, y, xlab = "", ylab = "Distance from nearest\nalveolar cell (µm)", plotTitle = "",signifTest = "kruskal", labelycoo=max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
	dan.save(md,paste0(OutDir,"distance_normalAlveolar_byCellState.pdf"))
	# density plots
	for (tp in c( "nearest_dist" ) ){
	   x = md$CellState
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md[,tp]
	   dan.densityPlot( paste0(OutDir,"distance_normalAlveolar_byCellState_DensityPlot.pdf"), y, x, groupinglab = "", xlab = "Distance from nearest alveolar cell (µm)", ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = these_colorz, fileWidth = 3, fileHeight = 2 )
	}

	### spatial segregation
	rdf = dan.df(0,c( "CellState","Side","nearest_dist" ))
	md_cancer = md[,c( "x","y" )]
	dd_CancerAt2 = cdist(md_cancer,md_cancer)
	rownames(dd_CancerAt2) = rownames(md_cancer)
	colnames(dd_CancerAt2) = rownames(md_cancer)
	diag(dd_CancerAt2) = NA

	md_shuffled = md
	set.seed(123)
	for (cs in names(cs_signatures)){
		dcat(cs)
		these_cells = rownames(md[md$CellState==cs,])
		min_dist = apply(dd_CancerAt2[these_cells,these_cells],1,min,na.rm=T)
		rdf = rbind(rdf,data.frame(CellState=cs,Side="Real",nearest_dist=as.numeric(min_dist),stringsAsFactors=F))
		for (nshuf in c(1:10)){
			dcat(nshuf,1)
			md_shuffled$CellState = sample(md_shuffled$CellState,nrow(md_shuffled))
			these_cells = rownames(md_shuffled[md_shuffled$CellState==cs,])
			min_dist = apply(dd_CancerAt2[these_cells,these_cells],1,min,na.rm=T)
			rdf = rbind(rdf,data.frame(CellState=cs,Side="Shuffled",nearest_dist=as.numeric(min_dist),stringsAsFactors=F))
		}
	}
	dan.save(rdf,paste0(OutDir,"SpatialSegregation_densities.pdf"))
	# removing long range distances
	rdf = rdf[rdf$nearest_dist<=as.numeric(quantile(rdf$nearest_dist,seq(0,1,0.01))["95%"]),]
	rdf$CellState = factor(rdf$CellState,levels=c("Alveolar","Proliferative","Hypoxic"))
	plot = ggplot(rdf)+aes(x=nearest_dist,fill=CellState,alpha=Side,linetype=Side) + geom_density() + scale_alpha_manual(values = c("Real" = 0.5, "Shuffled" = 0.1)) + scale_linetype_manual(values = c("Shuffled" = "dashed", "Real" = "solid")) +  scale_fill_manual(values=ct_colorz) + facet_wrap(~CellState,ncol=1) + labs(x="Distance to nearest cell from the same state",y="Density") + theme(legend.title = element_blank()) + theme_classic()
	pdf( paste0(OutDir,"SpatialSegregation_densities.pdf"),7,4 )
	print(plot)
	dev.off()

	# rdf = dan.df(0,c( "CellState","Side","nearest_dist" ))
	# for (cs in names(cs_signatures)){
	# 	dcat(cs)
	# 	these_cells = rownames(md[md$CellState==cs,])
	# 	min_dist = apply(dd_CancerAt2[these_cells,these_cells],1,min,na.rm=T)
	# 	rdf = rbind(rdf,data.frame(CellState=cs,Side="Intra-state",nearest_dist=as.numeric(min_dist),stringsAsFactors=F))
	# 	min_dist = apply(dd_CancerAt2[these_cells,!(colnames(dd_CancerAt2) %in% these_cells)],1,min,na.rm=T)
	# 	rdf = rbind(rdf,data.frame(CellState=cs,Side="Inter-state",nearest_dist=as.numeric(min_dist),stringsAsFactors=F))
	# }
	# quantile(rdf$nearest_dist,seq(0,1,0.01))
	# # removing long range distances
	# rdf = rdf[rdf$nearest_dist<=as.numeric(quantile(rdf$nearest_dist,seq(0,1,0.01))["95%"]),]
	# rdf$CellState = factor(rdf$CellState,levels=c("Alveolar","Proliferative","Hypoxic"))
	# plot = ggplot(rdf)+aes(x=nearest_dist,fill=CellState,alpha=Side,linetype=Side) + geom_density() + scale_alpha_manual(values = c("Intra-state" = 0.5, "Inter-state" = 0.1)) + scale_linetype_manual(values = c("Inter-state" = "dashed", "Intra-state" = "solid")) +  scale_fill_manual(values=ct_colorz) + facet_wrap(~CellState,ncol=1) + labs(x="Distance to nearest cancer cell from the same or a different state",y="Density") + theme(legend.title = element_blank()) + theme_classic()
	# pdf( paste0(OutDir,"SpatialSegregation_IntraInter_densities.pdf"),7,4 )
	# print(plot)
	# dev.off()

	# ### tps 
	# load( paste0(OutDir,"../../tps_discovery/tps.RData") )
	# genes_xenium = rownames(xe_AT12Cancer)
	# tps = lapply(tps, function(x) intersect( x,genes_xenium ))
	# tps = tps[sapply(tps,length)>1]
	# # xe_AT12Cancer = AddModuleScore(xe_AT12Cancer,tps,pool = NULL,nbin = 10,ctrl = 20,k = FALSE,assay = 'SCT',name = "ams")
	# # colnames(xe_AT12Cancer@meta.data)[substr(colnames(xe_AT12Cancer@meta.data),1,3)=="ams"] = names(tps)
	# dcat( "AmsCentered" )
	# db_rand = dan.barkley_MakeRand_data(xe_AT12Cancer,tps, 3)
	# scores = dan.Barkley_GeneToEnrichment_AmsCentered( xe_AT12Cancer, tps, db_rand)
	# xe_AT12Cancer@meta.data = scores
	# md = xe_AT12Cancer@meta.data
	# md$highest_tps = names(tps)[apply(md[,names(tps)],1,which.max)]
	# xe_AT12Cancer@meta.data = md
	# pdf(paste0(OutDir,"CancerCells_imageplot_highest_tps.pdf" ),12,12)
	# print(ImageDimPlot(xe_AT12Cancer,group.by="highest_tps",axes=T,border.size=NA))
	# dev.off()

	### neighbourhood analysis
	load(file=paste0(OutDir,"../../xenium_preprocessing_",whichDataset,"/xe.RData" ))
	md$annd_level_2 = md$CellState
	dtable(xe@meta.data[rownames(md),"annd_level_1"])
	dtable(xe@meta.data[rownames(md),"annd_level_2"])
	xe@meta.data[rownames(md),"annd_level_2"] = md$CellState
	md = xe@meta.data
	coords = GetTissueCoordinates(xe, which = "centroids")
	rownames(coords) = coords$cell
	coords$annd_level_1 = md[rownames(coords),"annd_level_1"]
	coords$annd_level_2 = md[rownames(coords),"annd_level_2"]
	in_file_squidpy = paste0(OutDir,"in_file_squidpy.txt")
	out_file_squidpy = paste0(OutDir,"out_file_squidpy.txt")
	dan.write(coords,in_file_squidpy)
	# run python nei_squidpy.py in_file_squidpy out_file_squidpy
	# run python nei_squidpy.py --infile "processed/cvNMF_pipeline/run_05/xenium_analyses_xenium_05February2024/cancer_cells/in_file_squidpy.txt" --outfile "processed/cvNMF_pipeline/run_05/xenium_analyses_xenium_05February2024/cancer_cells/out_file_squidpy.txt"
	if (file.exists(out_file_squidpy) ){
		cellstates_spatialneighbourhood_singlecells( md,out_file_squidpy,OutDir )
	}
}

xenium_gridding_cancercells = function( OutDir, whichDataset, grid_spacing = 200, density_threshold = 0.002 ){
	library(RColorBrewer)
	celltype_map = data.frame(row.names=c( "AT1","AT2","Cancer" ),colorz = c("chocolate","dodgerblue2","firebrick") ,stringsAsFactors=F )
	load( file=paste0(OutDir,"../../CellStates/signatures/wilcox_cs_signatures_restrictive.RData" ) )
	this_OutDir = paste0(OutDir,"griddingCC_",grid_spacing,"/")
	dir.create(this_OutDir)
	load(file = paste0( OutDir,"md_CancerCells.RData" ))
	xmax = ceiling(max(md$x)/grid_spacing)
  	ymax = ceiling(max(md$y)/grid_spacing)
  	# Cancer
  	md_cancer = md[md$annd_level_2=="Cancer",]
  	md_cancer$xg = factor(ceiling(md_cancer$x/grid_spacing), levels = c(1:xmax))
  	md_cancer$yg = factor(ceiling(md_cancer$y/grid_spacing), levels = c(1:ymax))
  	meaMat_cancer = t(table(md_cancer$xg, md_cancer$yg))
  	md_cancer$grid_coo = paste0( "xg",md_cancer$xg,"_yg",md_cancer$yg )
  	md_cancer_grid = aggregate( .~grid_coo,data=md_cancer[,c( "grid_coo",names(cs_signatures) )],FUN='mean' )
  	md_cancer_grid$xg = as.numeric(gsub("xg","",gsub("_.*","",md_cancer_grid$grid_coo)))
  	md_cancer_grid$yg = as.numeric(gsub("yg","",gsub(".*_","",md_cancer_grid$grid_coo)))
  	rownames(md_cancer_grid) = md_cancer_grid$grid_coo
  	md_cancer_grid$n_Cancer = table(md_cancer$grid_coo)[rownames(md_cancer_grid)]
	mmelt = melt(meaMat_cancer)
	colnames(mmelt) = c( "yg","xg","value" )
	mmelt$grid_coo = paste0( "xg",mmelt$xg,"_yg",mmelt$yg )
	mmelt = mmelt[!(mmelt$grid_coo %in% md_cancer_grid$grid_coo),]
	md_cancer_grid = rbind(md_cancer_grid,data.frame(row.names=mmelt$grid_coo,grid_coo=mmelt$grid_coo,Alveolar=NA,Proliferative=NA,Hypoxic=NA,xg=mmelt$xg,yg=mmelt$yg,n_Cancer=0,stringsAsFactors=F))

	md_cancer_grid$spatial_density_Cancer = md_cancer_grid$n_Cancer/(grid_spacing*grid_spacing)

	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_cancer_spatialdensity.pdf")
	pdf(fileName, 9, 6, useDingbats = F )
	dfmelt = melt(meaMat_cancer)
	dfmelt$value = 0
	rownames(dfmelt) = paste0( "xg",dfmelt$Var2,"_yg",dfmelt$Var1 )
	dfmelt[rownames(md_cancer_grid),"value"] = md_cancer_grid$spatial_density_Cancer
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="Density",colours=colorRampPalette(c("white",as.character(celltype_map["Cancer","colorz"]) ))(n = 100))+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("Density (N/µm²) of cancer cells, grid spacing = ",grid_spacing,"µm"), x="",y="")+
	  coord_fixed()
	print(p)
	dev.off()

	md_cancer_grid$cancer_poor = 0
	md_cancer_grid[(md_cancer_grid$spatial_density_Cancer<=density_threshold),"cancer_poor"] = 1
	fileName = paste0(this_OutDir,"Heatmap_tileSizeMicrons",grid_spacing,"_is_cancer_poor_onlyDensity.pdf")
	pdf(fileName, 9, 6, useDingbats = F )
	dfmelt = melt(meaMat_cancer)
	dfmelt$value = 0
	rownames(dfmelt) = paste0( "xg",dfmelt$Var2,"_yg",dfmelt$Var1 )
	dfmelt[rownames(md_cancer_grid),"value"] = md_cancer_grid$cancer_poor
	levelz = levels(factor(dfmelt$Var2))
	p=ggplot(dfmelt) +
	  geom_tile(aes(Var1,ordered(Var2, levels=(levelz)),fill=value))+
	  scale_fill_gradientn(name="Poor",colours=colorRampPalette(c("white",as.character(celltype_map["Cancer","colorz"]) ))(n = 100))+
	  # scale_fill_gradientn(name="hlca_AT1",colours=colorRampPalette(c("white",as.character(celltype_map["AT1","colorz"]) ))(n = 100))+
	  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	  labs(title=paste0("Is cancer-poor, grid spacing = ",grid_spacing,"µm"), x="",y="")+
	  coord_fixed()
	print(p)
	dev.off()

	md_cancer_poor = md_cancer_grid[(md_cancer_grid$spatial_density_Cancer<=density_threshold),c( "xg","yg" )]
	if (nrow(md_cancer_poor)<2){ next }

	library(rdist)
	dd_CancerAt2 = cdist(md_cancer_poor,md_cancer_grid[,c( "xg","yg" )])
	rownames(dd_CancerAt2) = rownames(md_cancer_poor)
	colnames(dd_CancerAt2) = rownames(md_cancer_grid)
	min_dist = apply(dd_CancerAt2,2,min)
	md_cancer_grid[colnames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)*grid_spacing
	
	pdf( paste0(this_OutDir,"Cancer_cs_signatures_vsDistance_tileSizeMicrons",grid_spacing,".pdf"),6,5 )
	for (tp in intersect(names(cs_signatures),colnames(md_cancer_grid))){
	   x = md_cancer_grid$nearest_dist
	   y = md_cancer_grid[,tp]
	   plotTitle = paste0( "Spearman R = ",signif(cor(x,y,method="spearman",use='complete.obs'),2)," p-value = ",signif(cor.test(x,y,method="spearman")$p.value,2)  )
	   plot=dan.scatterplots.multipages( x, y, xlab = "Distance from nearest cancer-poor tile (µm)", ylab = paste0(tp," score"), plotTitle = plotTitle, dotSize = 0.5, plotFitLine = "yes", FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
	   print(plot)
	}
	dev.off()

	md_cancer_grid$dist_discrete = ">1mm"
	md_cancer_grid[md_cancer_grid$nearest_dist<=250,"dist_discrete"] = "<0.25mm"
	md_cancer_grid[(md_cancer_grid$nearest_dist>250) & (md_cancer_grid$nearest_dist<=500),"dist_discrete"] = "0.25-0.5mm"
	md_cancer_grid[(md_cancer_grid$nearest_dist>500) & (md_cancer_grid$nearest_dist<=750),"dist_discrete"] = "0.5-0.75mm"
	md_cancer_grid[(md_cancer_grid$nearest_dist>750) & (md_cancer_grid$nearest_dist<=1000),"dist_discrete"] = "0.75-1mm"
	this_md_cancer_grid = md_cancer_grid[!is.na(md_cancer_grid$Alveolar),]

	levelz = c( "<0.25mm","0.25-0.5mm","0.5-0.75mm","0.75-1mm",">1mm" )
	colorz = c( "firebrick4","tomato3","orange","gray66","gray33" )
	pdf( paste0(this_OutDir,"Cancer_cs_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots.pdf"),6,5 )
	for (tp in intersect(names(cs_signatures),colnames(this_md_cancer_grid)) ){
	   x = this_md_cancer_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = this_md_cancer_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-poor tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	library(viridis)

	md_cancer_grid$dist_discrete = ">=1mm"
	md_cancer_grid[md_cancer_grid$nearest_dist<grid_spacing,"dist_discrete"] = paste0("<",grid_spacing/1000,"mm")
	levelz = c(paste0("<",grid_spacing/1000,"mm"))
	for (it in 1:(ceiling(1000/grid_spacing)-1)  ){
		lowerbound = grid_spacing*it
		upperbound = grid_spacing*(it+1)
		md_cancer_grid[(md_cancer_grid$nearest_dist>=lowerbound) & (md_cancer_grid$nearest_dist<upperbound),"dist_discrete"] = paste0(lowerbound/1000,"-",upperbound/1000,"mm")
		levelz = c(levelz,paste0(lowerbound/1000,"-",upperbound/1000,"mm"))
	}
	this_md_cancer_grid = md_cancer_grid[!is.na(md_cancer_grid$Alveolar),]
	levelz = c(levelz,">=1mm")
	colorz = (viridis::plasma(length(levelz)))
	pdf( paste0(this_OutDir,"Cancer_cs_signatures_vsDistance_tileSizeMicrons",grid_spacing,"_Boxplots_SameGridSpacing.pdf"),6,5 )
	for (tp in intersect(names(cs_signatures),colnames(this_md_cancer_grid)) ){
	   x = this_md_cancer_grid$dist_discrete
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = this_md_cancer_grid[,tp]
	   ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-poor tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()

	md$nearest_dist = NA
	md_cancer = md[md$annd_level_2=="Cancer",c( "x","y" )]
	md_cancer_poor$x = md_cancer_poor$xg*grid_spacing
	md_cancer_poor$y = md_cancer_poor$yg*grid_spacing
	md_at2 = md_cancer_poor[,c( "x","y" )]
	dd_CancerAt2 = cdist(md_cancer,md_at2)
	rownames(dd_CancerAt2) = rownames(md_cancer)
	colnames(dd_CancerAt2) = rownames(md_at2)
	min_dist = apply(dd_CancerAt2,1,min)
	md[rownames(dd_CancerAt2),"nearest_dist"] = as.numeric(min_dist)

	levelz = c("Alveolar","Proliferative","Hypoxic")
	colorz = c( "dodgerblue","red","purple" )
	pdf( paste0(OutDir,"distance_CancerCellPoorTiles_byCellState_othercolors.pdf"),8,6 )
	for (tp in c( "nearest_dist" ) ){
	   x = md$CellState
	   these_colorz = colorz[levelz %in% unique(x)]
	   x = factor(x, levels = levelz[levelz %in% unique(x)])
	   y = md[,tp]
	   # ylimLeft = min(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["25%"]-(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   ylimRight = max(sapply( unique(x), function(lev) {quantile(y[x==lev],seq(0,1,0.25))["75%"]+(quantile(y[x==lev],seq(0,1,0.25))["75%"]-quantile(y[x==lev],seq(0,1,0.25))["25%"])*1.5} ))
	   # plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0("Distance from nearest \nnon-malignant alveolar cell"), plotTitle = "",signifTest = NULL, labelycoo = max(y), ylimLeft=ylimLeft, ylimRight=ylimRight, xColors = these_colorz, includeJitters = F)
	   plot=dan.boxplots.multipages( x, y, xlab = "", ylab = paste0("Distance from nearest cancer cell-poor tile"), plotTitle = "",signifTest = "kruskal", labelycoo = ylimRight, ylimRight = ylimRight, ylimLeft=0, xColors = these_colorz, includeJitters = F)
	   # plot=dan.boxplots.multipages( x, y, xlab = "Distance from nearest cancer-rich tile", ylab = paste0(tp," score"), plotTitle = "",signifTest = "kruskal", labelycoo = max(y), xColors = these_colorz, includeJitters = F)
	   print(plot)
	}
	dev.off()
}

TME_formatted_umaps = function( OutDir ){
	CellTypeDir = paste0("data/cell_typing/non_malignant/")

	ctdf = dan.df(0,c( "annd_level_1","annd_level_2","annd_level_3","stiched" ))
	# for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu")){
	for (dataset in c("wang")){ # wang as an example
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_annot_annd.RData"))
		annot$stiched = paste0(annot$annd_level_1,"_",annot$annd_level_2,"_",annot$annd_level_3)
		ctdf = rbind(ctdf,annot[,c("annd_level_1","annd_level_2","annd_level_3","stiched" )])
		ctdf = ctdf[!duplicated(ctdf$stiched),]	
	}
	ctdf = ctdf[!duplicated(ctdf$stiched),]
	ctdf = ctdf[!(ctdf$annd_level_1=="unclear"), ]
	ctdf = ctdf[order(ctdf$annd_level_1,ctdf$annd_level_2,ctdf$annd_level_3),]
	all_ct_1 = unique(ctdf$annd_level_1)[!(unique(ctdf$annd_level_1)=="unclear")]
	all_ct_2 = unique(ctdf$annd_level_2)[!(unique(ctdf$annd_level_2)=="unclear")]
	all_ct_3 = unique(ctdf$annd_level_3)[!(unique(ctdf$annd_level_3)=="unclear")]
	# ct_colorz_1 = dan.colors(length(all_ct_1))
	ct_colorz_1 = c( "goldenrod1","firebrick2","dodgerblue3","gray" )
	names(ct_colorz_1) = all_ct_1
	# ct_colorz_2 = dan.colors(length(all_ct_2))
	ct_colorz_2 = c( "gold","goldenrod1",
							"firebrick3","darkorange1","orangered","sienna4","tan3","sandybrown","wheat2",
							"forestgreen","deeppink2","purple","purple4","mediumpurple1","midnightblue","dodgerblue2",
							"gray22","gray55","gray88")
	names(ct_colorz_2) = all_ct_2
	# ct_colorz_3 = dan.colors(length(all_ct_3))
	ct_colorz_3 = c( "gold","goldenrod1","firebrick3","darkorange1","orangered","sienna4","tan3","sandybrown","wheat2",
							"green4","green3","yellowgreen","lightgreen",
							"deeppink3","deeppink1","pink2","pink",
							"purple2","mediumpurple3","purple4","mediumpurple2","mediumorchid",
							"midnightblue","dodgerblue4","dodgerblue3","dodgerblue2","dodgerblue1","lightblue3","lightblue2","lightblue1",
							"turquoise4","turquoise3","turquoise2","turquoise1","aquamarine3","aquamarine2","aquamarine1",
							"gray22","gray55","gray88")
	names(ct_colorz_3) = all_ct_3

	# for (dataset in c("qian","kim","wu","laughney","xing","bischoff","yang","zhu","salcher","he","wang","hu")){
	for (dataset in c("wang")){
		dcat(dataset)
		load(file = paste0(CellTypeDir,"final_seus_annots/",dataset,"_seu.RData"))

		# level1
		subseu = subset(seu,cells=rownames(seu@meta.data[seu@meta.data$annd_level_1!="unclear",]))
    	plot=(DimPlot(subseu,group.by = "annd_level_1", label = F,raster=F,cols=ct_colorz_1) + coord_fixed()+NoLegend())
    	ggsave(file=paste0(OutDir,dataset,"_celltyping_level1_assigned.png"),plot,width=6,height=6,dpi=1200 )
    	pdf(paste0(OutDir,dataset,"_celltyping_level1_assigned_rastered.pdf"),2.7,2.7,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_1", label = T,raster=T,cols=ct_colorz_1) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level1_assigned_unrastered.pdf"),2.7,2.7,pointsize=6)
		plot=(dan.DimPlot(subseu,pt.size=0.25,label.size=6/.pt, group.by = "annd_level_1", label = F,raster=F,cols=ct_colorz_1) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level1_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_1", label = T,raster=T,cols=ct_colorz_1) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot)
    	dev.off()

    	# level2
		subseu = subset(seu,cells=rownames(seu@meta.data[seu@meta.data$annd_level_2!="unclear",]))
    	plot=(DimPlot(subseu,group.by = "annd_level_2", label = F,raster=F,cols=ct_colorz_2) + coord_fixed()+NoLegend())
    	ggsave(file=paste0(OutDir,dataset,"_celltyping_level2_assigned.png"),plot,width=6,height=6,dpi=1200 )
    	pdf(paste0(OutDir,dataset,"_celltyping_level2_assigned_rastered.pdf"),2.7,2.7,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_2", label = T,raster=T,cols=ct_colorz_2) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level2_assigned_unrastered.pdf"),2.7,2.7,pointsize=6)
		plot=(dan.DimPlot(subseu,pt.size=0.25,label.size=6/.pt, group.by = "annd_level_2", label = F,raster=F,cols=ct_colorz_2) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level2_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_2", label = T,raster=T,cols=ct_colorz_2) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot)
    	dev.off()

    	# level3
		subseu = subset(seu,cells=rownames(seu@meta.data[seu@meta.data$annd_level_3!="unclear",]))
    	plot=(DimPlot(subseu,group.by = "annd_level_3", label = F,raster=F,cols=ct_colorz_3) + coord_fixed()+NoLegend())
    	ggsave(file=paste0(OutDir,dataset,"_celltyping_level3_assigned.png"),plot,width=6,height=6,dpi=1200 )
    	pdf(paste0(OutDir,dataset,"_celltyping_level3_assigned_rastered.pdf"),2.7,2.7,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_3", label = T,raster=T,cols=ct_colorz_3) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level3_assigned_unrastered.pdf"),2.7,2.7,pointsize=6)
		plot=(dan.DimPlot(subseu,pt.size=0.25,label.size=6/.pt, group.by = "annd_level_3", label = F,raster=F,cols=ct_colorz_3) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level3_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_3", label = T,raster=T,cols=ct_colorz_3) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot)
    	dev.off()

    	for (thislev in all_ct_2){
    		# level3 split
			subseu = subset(seu,cells=rownames(seu@meta.data[(seu@meta.data$annd_level_3!="unclear") & (seu@meta.data$annd_level_2==thislev),]))
			these_ct_colorz_3 = ct_colorz_3[unique(as.character(subseu@meta.data$annd_level_3))]
	    	plot=(DimPlot(subseu,group.by = "annd_level_3", label = F,raster=F,cols=these_ct_colorz_3) + coord_fixed()+NoLegend())
	    	ggsave(file=paste0(OutDir,dataset,"_celltyping_level3_sub_",thislev,"_assigned.png"),plot,width=6,height=6,dpi=1200 )
	    	pdf(paste0(OutDir,dataset,"_celltyping_level3_sub_",thislev,"_assigned_rastered.pdf"),2.7,2.7,pointsize=6)
			plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_3", label = T,raster=T,cols=these_ct_colorz_3) + coord_fixed())
			plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	    	print(plot+NoLegend())
	    	dev.off()
	    	pdf(paste0(OutDir,dataset,"_celltyping_level3_sub_",thislev,"_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
			plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_3", label = T,raster=T,cols=these_ct_colorz_3) + coord_fixed())
			plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	    	print(plot)
	    	dev.off()
    	}

    	# epithelial cells
		subseu = subset(seu,cells=rownames(seu@meta.data[(seu@meta.data$annd_level_2!="unclear") & (seu@meta.data$annd_level_1=="epithelial"),]))
		these_ct_colorz_2 = ct_colorz_2[unique(as.character(subseu@meta.data$annd_level_2))]
    	plot=(DimPlot(subseu,group.by = "annd_level_2", label = F,raster=F,cols=ct_colorz_2) + coord_fixed()+NoLegend())
    	ggsave(file=paste0(OutDir,dataset,"_celltyping_level2_sub_epithelial_assigned.png"),plot,width=6,height=6,dpi=1200 )
    	pdf(paste0(OutDir,dataset,"_celltyping_level3_sub_epithelial_assigned_rastered.pdf"),2.7,2.7,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_3", label = T,raster=T,cols=these_ct_colorz_2) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot+NoLegend())
    	dev.off()
    	pdf(paste0(OutDir,dataset,"_celltyping_level2_sub_epithelial_assigned_rastered_withLegend.pdf"),6,6,pointsize=6)
		plot=(DimPlot(subseu,pt.size=1/.pt,label.size=6/.pt, group.by = "annd_level_3", label = T,raster=T,cols=these_ct_colorz_2) + coord_fixed())
		plot = plot + theme_classic(base_size=6) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
    	print(plot)
    	dev.off()
	}
}

facs_data_analyses = function( OutDir ){
	df = dan.read(file = "data/facs_data_siKeap1_siStk11/facs_mhcii_percentages.txt")
	meltdf = melt(df)
	meltdf = meltdf[!is.na(meltdf$value),]
	fileName = paste0(OutDir,"barplot.pdf")
	x = factor(meltdf$variable,levels=c( "parental","siCtrl","siKeap1","siStk11" ))
	y = meltdf$value
	y = apply(df,2,function(x) median(x[!is.na(x)]))
	sd = apply(df,2,function(x) sd(x[!is.na(x)]))
	se = apply(df,2,function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)])))
	low25 = apply(df,2,function(x) quantile(x,seq(0,1,0.25),na.rm=T)['25%'] )
	top75 = apply(df,2,function(x) quantile(x,seq(0,1,0.25),na.rm=T)['75%'] )
	# dan.barplot(fileName,levels(x),y,sd=se[levels(x)],ylab="% of MHC-II+ cells",xlab="",fileWidth=2, fileHeight=2)
	# computing pvals
	library(ggpubr)
	library(rstatix)
	pval_df = data.frame(row.names=paste0("test",1:3),group1=c( "parental","siCtrl","siCtrl" ),group2=c( "siCtrl","siKeap1","siStk11" ),
		y.position=c(12.7,11.4,10.1),p=NA,p.adj=NA )
	pval_df["test1","p"] = wilcox.test( df$parental,df$siCtrl )$p.value
	pval_df["test2","p"] = wilcox.test( df$siCtrl,df$siKeap1 )$p.value
	pval_df["test3","p"] = wilcox.test( df$siCtrl,df$siStk11 )$p.value
	pval_df[,"p.adj"] = paste0("p = ", signif(p.adjust(pval_df[,"p"],method="BH"),1))
	meltdf$color = "gray77"
	meltdf[meltdf$variable=="siCtrl","color"] = "gray44"
	meltdf[meltdf$variable=="siKeap1","color"] = "forestgreen"
	meltdf[meltdf$variable=="siStk11","color"] = "dodgerblue3"
	se = se[levels(x)]
	pdf(fileName,1.9,1.9)
	plot = ggplot(mapping = aes(y = y[levels(x)], x = levels(x))) + geom_bar(size=0.25,stat="identity", position=position_dodge2(width = 0.9, preserve = "single"),fill="white",color="black")
  	plot = plot + geom_errorbar(size=0.25,aes(ymin=y-se, ymax=y+se), width=.4) + geom_jitter(aes(x = x, y = meltdf$value),width=.2,height=0, size = 1, alpha = 1, fill = meltdf$color,pch=21,colour="black") #geom_point(aes(x = x, y = meltdf$value),size = 1,position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0))
  	plot = plot + xlab( "" ) + ylab( "% of MHC-II+ cells" )
  	plot = plot + theme_classic(base_size=6) + theme(axis.title.x=element_blank(),text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6,angle = 45, hjust = 1), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
  	plot = plot + stat_pvalue_manual( pval_df, label = "p.adj",label.size = 6/.pt)
  	print(plot)
	dev.off()
}
