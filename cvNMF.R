
#################### Main functions

cvNMF_preprocessing = function( annot, counts, OutDir ){
	dcat( "Preprocessing" )
	dcat( "Normalizing",1 )
	seu = CreateSeuratObject(counts = counts, project = "qkwlxb", meta.data = annot)
	seu[["percent.mt"]] = PercentageFeatureSet(seu, pattern = "^MT-")
	seu = SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"), verbose = FALSE)
	ribo = read.csv(file = paste0(DataDir,"hgnc_ribosomial_proteins.csv"))
	mito.genes = rownames(seu)[grepl(pattern = "^MT-", x = rownames(seu))]
	ribo.genes = ribo$Approved.symbol
	seu = subset(x = seu, cells = rownames(annot), features = rownames(seu)[!(rownames(seu) %in% c(mito.genes,ribo.genes))])
	seuNorm = seu[['SCT']]@data
	annot$CellID = rownames(annot)
	snp = data.frame(seuNorm)
	snp0 = Matrix(as.matrix(snp[(rowSums(snp))>0,]),sparse=T)
	snp0 = snp0[!(rownames(snp0) %in% c(mito.genes,ribo.genes)),]
	colnames(snp0) = gsub("\\.","-",colnames(snp0))
	dcat( "Centering",1 )
	snp0 = snp0-rowMeans(snp0) # centering
	snp0[snp0<0] = 0 # setting to zero all negative values
	gene_universe = rownames(snp0)
	dcat( "Saving",1 )
	save(annot,file=paste0(OutDir,"annot.RData"))
	save(snp0,file=paste0(OutDir,"preprocessed_expression_mat.RData"))
	save(gene_universe,file=paste0(OutDir,"gene_universe.RData"))
}

cvNMF_partitioned_runs = function( annot, snp0, chosen_rank, ntrials, unbalanced_size_margin, cor_thresh, this_OutDir ){
	dcat("Computing partitioned patients models and matching LFs",2)
	model_p1 = list()
	model_p2 = list()
	all_patients = unique(annot$Sample)
	nt = 1
	while (nt <= ntrials)
	{
		p1 = sample( all_patients, size = length(all_patients)/2 )
		p2 = all_patients[!(all_patients %in% p1)]
		len_p1 = sum(annot$Sample %in% p1)
		len_p2 = sum(annot$Sample %in% p2)
		if ( abs(len_p1-len_p2) < nrow(annot)*unbalanced_size_margin )
		{
			# dcat(nt, 3)
			model = RcppML::nmf(snp0[,annot[annot$Sample %in% p1,"CellID"]], k = chosen_rank, maxit = 10000, tol = 1e-10, verbose = F)
			model_p1[[ nt ]] = model
			model = RcppML::nmf(snp0[,annot[annot$Sample %in% p2,"CellID"]], k = chosen_rank, maxit = 10000, tol = 1e-10, verbose = F)
			model_p2[[ nt ]] = model
			nt = nt+1
		}	
	}

	### Comparing LFs
	m_list = list()
	wpool_list = list()
	assigned_list = list()
	for (nt in 1:ntrials)
	{
		mod = model_p1[[nt]]
		w1 = mod$w
		mod = model_p2[[nt]]
		w2 = mod$w
		m = matrix(NA, nrow = chosen_rank, ncol = chosen_rank)
		for (lf1 in 1:chosen_rank)
		{
			for (lf2 in 1:chosen_rank)
			{
				m[lf1,lf2] = cor(w1[,lf1],w2[,lf2])
			}
		}
		mplot = m
		mplot[mplot<cor_thresh] = 0
		# pdf(paste0(this_OutDir,"Intra_corrplot_w_",nt,".pdf"),15,15)
		# corrplot(mplot,addCoef.col = 'white')
		# dev.off()
		mplot = data.frame(mplot)
		m_list[[nt]] = mplot
		assigned_df = AssignLfs(m,0)
		assigned_list[[nt]] = assigned_df
		assigned_df = assigned_df[assigned_df$cor > cor_thresh,]
		wpool = matrix(NA, nrow = nrow(w1), ncol = nrow(assigned_df))
		if (ncol(wpool)>0)
		{
			for (cn in 1:ncol(wpool))
			{
				wpool[,cn] = rowMeans(cbind(w1[,assigned_df[cn,"p1"]],w2[,assigned_df[cn,"p2"]]))

			}
		}

		wpool_list[[nt]] = wpool
	}

	save(m_list, file = paste0(this_OutDir,"m_list.RData"))
	save(assigned_list, file = paste0(this_OutDir,"assigned_list.RData"))
	save(wpool_list, file = paste0(this_OutDir,"wpool_list.RData"))
}

cvNMF_allpatients_run = function(snp0, chosen_rank, ntrials, cor_thresh, recurrence_ratio, this_OutDir){

	### Computing master model and comparing
	dcat("Computing all-patients model and comparing",2)

	load(file = paste0(this_OutDir,"wpool_list.RData"))
	master_mod = RcppML::nmf(snp0, k = chosen_rank, maxit = 10000, tol = 1e-10, verbose = F, seed = 1)
	master_robust = c()
	for (nt in 1:ntrials)
	{
		# dcat(nt,3)
		wpool = wpool_list[[nt]]
		if (ncol(wpool)==0) { next }
		wmaster = master_mod$w
		m = matrix(NA, nrow = chosen_rank, ncol = ncol(wpool))
		for (lf1 in 1:chosen_rank)
		{
			for (lf2 in 1:ncol(wpool))
			{
				m[lf1,lf2] = cor(wmaster[,lf1],wpool[,lf2])
			}
		}
		mplot = m
		mplot[mplot<cor_thresh] = 0
		pdf(paste0(this_OutDir,"AllPatients_Intra_corrplot_w_",nt,".pdf"),15,15)
		corrplot(mplot,addCoef.col = 'white')
		dev.off()
		mplot = data.frame(mplot)
		# dcat(which(rowSums(mplot)>0))
		master_robust = c(master_robust,which(rowSums(mplot)>0))
	}
	master_robust = as.numeric(names(which(table(master_robust)/ntrials >= recurrence_ratio)))
	if (length(master_robust)==0) { return() }
	w = master_mod$w
	wrobust = data.frame(w[,master_robust])
	colnames(wrobust) = paste0("lf",master_robust)
	h = master_mod$h
	hrobust = data.frame(h[master_robust,])
	if (length(master_robust)==1) { hrobust = data.frame(t(hrobust)) }
	rownames(hrobust) = paste0("lf",master_robust)
	save(wrobust, file = paste0(this_OutDir,"wrobust.RData"))
	save(hrobust, file = paste0(this_OutDir,"hrobust.RData"))
}

cvNMF_adjmat_construction = function(OutDir, kvec, exclude_genes = NULL){
	load(file=paste0(OutDir,"gene_universe.RData"))
	flat_list = list()
	flat_index = 1
	for (k in kvec)
	{
	   this_file = paste0(OutDir,"cvNMF_runs/","rank",k,"/","wrobust.RData")
	   if (!file.exists(this_file)) { next }
	   dcat(k)
	   load(this_file)
	   w = wrobust
	   a = ExtractGenes(w, gene_universe, n = 100, method = "max_NMF_package", exclude_genes = exclude_genes)
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

cvNMF_adjmat_clustering = function(OutDir, res_range, minimum_n_lfs = NULL){
	load(file=paste0(OutDir,"gene_universe.RData"))
	this_OutDir = paste0(OutDir,"louvain_clustering/")
	dir.create(this_OutDir)
	load(file = paste0(OutDir,"full_adjmat.RData"))
	dcat( paste0( "Matrix size: ", nrow(adjmat)) )
	if (!is.null(minimum_n_lfs)){
		adjmat2 = adjmat
		adjmat2[adjmat2<minimum_n_lfs] = 0
		namez = sort(names(which(rowSums(adjmat2)>0)))
		adjmat = adjmat2[namez,namez]
		dcat( paste0( "Filtered matrix size: ", nrow(adjmat)) )
	}
	set.seed(123)
	gw = graph_from_adjacency_matrix(adjmat, mode = "undirected", weighted = T)

	commdf = dan.df(rownames(adjmat), paste0("ig", res_range,"_sol1"), data = NA, as_df = T)
	for (res in res_range){
	   dcat(res,2)
	   nt = cluster_louvain(gw,resolution=res)
	   for (rnn in 1:nrow(nt$memberships)){ 
	      commdf[nt$names,paste0("ig", res,"_sol",rnn)] = as.numeric(nt$memberships[rnn,])
	   }
	}	
	### clustree analysis
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

	save(commdf, file = paste0(this_OutDir,"commdf.RData"))
	load( file = paste0(this_OutDir,"commdf.RData"))
	commdf$gene = rownames(commdf)
	for (which_clustering in colnames(commdf)[colnames(commdf)!="gene"]){
	   dir.create(paste0(this_OutDir,which_clustering))
	   tps = list()
	   for (ncomm in unique(commdf[,which_clustering]) )
	   {
	      these_wrGenes = commdf[commdf[,which_clustering]==ncomm,"gene"]
	      tps[[paste0("c",ncomm )]] = sort(these_wrGenes)
	   }
	   save(tps,file = paste0(this_OutDir,which_clustering, "/tps.RData"))
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

gsea_enrichment = function(OutDir, tps, exclude_genes){
	dir.create(paste0(OutDir,"gsea_enrichment/"))
	load(file=paste0(OutDir,"gene_universe.RData"))
	gene_universe = gene_universe[!(gene_universe %in% exclude_genes)]
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	for (f in names(tps)){
		genez = tps[[f]]
		fgRes = fgsea::fora(pathways = msig_list,
                       genes = genez,
                       universe = gene_universe)
		fgRes = fgRes[fgRes$padj <= 0.01,]
		a = as.data.frame(fgRes)
		a = apply(a,2,as.character)
		write.table(a,file=paste0( OutDir,"gsea_enrichment/",f,"_MSigDB_complete.txt" ),col.names=T,row.names=F,quote=F,sep="\t")
	}
}

############# Utils

AssignLfs = function( m1, cor_thresh ){
	el = melt(m1)
	colnames(el) = c("p1","p2","value" )
	el = el[order(-el$value),]
	el = el[el$value>cor_thresh,]
	assigned_df = data.frame(matrix(nrow=0,ncol=3,dimnames=list(NULL,c("p1","p2","cor"))))
	while (nrow(el) > 0)
	{
		this_p1 = el[1,"p1"]
		this_p2 = el[1,"p2"]
		this_cor = el[1,"value"]
		assigned_df = rbind(assigned_df,c("p1"=this_p1,"p2"=this_p2, "cor"=this_cor))
		el = el[(el$p1!=this_p1) & (el$p2!=this_p2),]
	}
	colnames(assigned_df) = c("p1","p2","cor")
	return( assigned_df )
}

dcat = function( string, tabb = 0 ){
	cat("\n", paste0(rep(".",tabb*5),collapse=""), string, "\n" )
}

ExtractGenes = function(w, snp0_genes, n = 100, method = "max_NMF_package", exclude_genes = NULL){
   genes = list()
   rownames(w) = snp0_genes
   single_lf = (ncol(w)==1)
   if (!is.null(exclude_genes)){
   	w = w[!(rownames(w) %in% exclude_genes),]
   	snp0_genes = snp0_genes[!(snp0_genes %in% exclude_genes)]
   }
   if (single_lf){
   	w = data.frame(row.names=snp0_genes,lf1=w)
   	genes[["lf1"]] = rownames(w)[order(w[,1],decreasing=T)][1:n] # in case there is only one lf, max_NMF_package does not make sense - so we keep the first 'n' genes by NMF loadings
   } else if (method == "max_NMF_package") # aka 'barkley'
   {
      res <- lapply(1:ncol(w), 
            function(i){
               mat <- w
               vect <- mat[,i]
               #order by decreasing contribution to factor i
               index.sort <- order(vect, decreasing=TRUE)      
               
               for( k in seq_along(index.sort) )
               {
                  index <- index.sort[k]
                  #if the feature contributes more to any other factor then return the features above it
                  if( any(mat[index,-i] >= vect[index]) )
                  {
                     if( k == 1 ) return(as.integer(NA))
                     else return( index.sort[1:(k-1)] )
                  }
               }
               
               # all features meet the criteria
               seq_along(vect)
            }
      )
      names(res) = colnames(w)
      for (lf in colnames(w))
      {
         genes[[lf]] = rownames(w)[res[[lf]]]
      }
   }
   return( genes )
}

dan.df = function(rownames, colnames, data = NA, as_df = T){
	if ((length(rownames)==1)) { mat = matrix(nrow = 0, ncol = length(colnames), dimnames = list(NULL,colnames)) }
	if ((length(rownames)!=1)) { mat = matrix(nrow = length(rownames), ncol = length(colnames), dimnames = list(rownames,colnames), data = data) }
	mat = data.frame(mat, stringsAsFactors = F)
	colnames(mat) = colnames
	return( mat )
}

###############################################
