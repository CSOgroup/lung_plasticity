
library(ggplot2)
library(ggpubr)
library(gplots)
library(plyr)
library(ggrepel)
library(ggbreak)

dan.save = function( object, file )
{
	formatted_filename = paste0(gsub("\\..*","",file),".RData")
	save(object, file = formatted_filename)
	formatted_filename = paste0(gsub("\\..*","",file),".txt")
	dan.write(object,formatted_filename)
}

dan.memory = function(){
	print(sort( sapply(ls(),function(x){object.size(get(x))})))
}

dan.read = function( file, row.names = NULL, header = T )
{
	rt = read.table(file = file, header = header, row.names = row.names, stringsAsFactors = F, sep = "\t", quote = '')
	return( rt )
}

dan.write = function( table, file, row.names = F )
{
	write.table(table, file = file, row.names = row.names, col.names = T, sep = "\t", quote = F)
	return()
}

peek = function( table, nrow = 5, ncol = 5 )
{
	nrow = min(c(nrow,nrow(table)))
	ncol = min(c(ncol,ncol(table)))
	print(table[1:nrow,1:ncol])
	print(paste0(nrow(table)," x ",ncol(table) ))
}

dan.colors = function( n, ggplot2_stile = T ){
	if (length(n)>1) { n = length(unique(n)) }
	if (ggplot2_stile){
		hues = seq(15, 375, length = n + 1)
		hcl(h = hues, l = 65, c = 100)[1:n]	
	}
}

dcat = function( string, tabb = 0 )
{
	cat("\n", paste0(rep(".",tabb*5),collapse=""), string, "\n" )
}

ddup = function( df, column )
{
	dupl = unique(df[,column][duplicated(df[,column])] )
	df = df[df[,column] %in% dupl,]
	df = df[order(df[,column]),]
	return(df)
}

dtable = function(...){
   return(table(...,useNA = "ifany"))
}

dan.expand_colors = function( vec, levelz, colorz ){
	colorz_vec = rep(NA,length(vec))
	for (l in levelz){
		colorz_vec[vec==l] = colorz[levelz==l]
	}
	return(colorz_vec)
}

dan.rowMedians = function( table, na.rm = T )
{
	return( apply(table,1,median, na.rm = na.rm) )
}

dan.colMedians = function( table, na.rm = T )
{
	return( apply(table,2,median, na.rm = na.rm) )
}

dan.boxplots = function( fileName, x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", signifLabel = "p.format", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, jitterDotSize = 4.5, fileWidth = 4, fileHeight = 3, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL, jitterShape = NULL )
{
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(fill))
	{
		if (!("default" %in% fillColors) )
		{
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,fatten=3) + scale_fill_manual(values = fillColors) + geom_point(pch = 16, size = jitterDotSize, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=6),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		} else {
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,fatten=3) + geom_point(pch = 16, size = jitterDotSize, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=6),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		}		
	} else {
		# cat(xlab,ylab) label = "p.format",
		p = ggplot(mapping = aes(y = y, x = x)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,fatten=3) +
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
			theme_classic(base_size=6) + theme(text = element_text(size=6),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
		
		if (!is.null(labelJitteredPoints))
		{
			p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 6/.pt, fontface = "bold", position = position_jitter(width = 0.2, height = 0), color = jitterColors)
		} else {
			if (!is.null(jitterShape)){
				js=jitterShape
				p = p + geom_jitter(aes(shape=js),width=.2,height=0,stroke=1,size = jitterDotSize, alpha = 1, color = jitterColors) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
			} else {
				p = p + geom_jitter(width=.2,height=0, size = jitterDotSize, alpha = 0.5, color = jitterColors) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
			}
		}
	}
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo, size = 6/.pt, label = signifLabel)}
	if (!is.null(hlines_coo))
	{
		p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" )# + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1, hjust = -1), color = "gray44")
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(legend_position)){
		p = p + theme(legend.position=legend_position)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.boxplots.multipages = function( x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	if (!is.null(fill))
	{
		if (!("default" %in% fillColors) )
		{
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,linewidth = 0.2) + scale_fill_manual(values = fillColors) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		} else {
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,linewidth = 0.2) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		}
		if (includeJitters){
				p = p + geom_point(position=position_jitterdodge(jitter.width = 0.1,jitter.height = 0),size=0.2)
			}
	} else {
		# cat(xlab,ylab) label = "p.format",
		p = ggplot(mapping = aes(y = y, x = x)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,linewidth = 0.2) +
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
			theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
		
		if (!is.null(labelJitteredPoints))
		{
			p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 6/.pt, position = position_jitter(width = 0.2, height = 0), color = jitterColors)
		} else {
			if (includeJitters){
				p = p + geom_jitter(width=.1, color = jitterColors,size=0.2)	
			}
		}
	}
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo, label = "p.format", size = 6/.pt)}
	if (!is.null(hlines_coo))
	{
		p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" ) + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1), color = "gray44",size = 6/.pt)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(legend_position)){
		p = p + theme(legend.position=legend_position)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.violinplots.multipages = function( x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	if (!is.null(fill))
	{
		if (!("default" %in% fillColors) )
		{
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_violin(linewidth=0.2, outlier.shape = NA) + scale_fill_manual(values = fillColors) + scale_color_manual(values=xColors) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		} else {
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_violin(linewidth=0.2, outlier.shape = NA) + scale_color_manual(values=xColors) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		}
		if (includeJitters){
				p = p + geom_point(position=position_jitterdodge(jitter.width = 0.2,jitter.height = 0))
			}
	} else {
		# cat(xlab,ylab) label = "p.format",
		p = ggplot(mapping = aes(y = y, x = x, color=x, show.legend = FALSE)) + geom_violin(linewidth=0.2, outlier.shape = NA) + scale_color_manual(values=xColors) +
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + 
			theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
		
		if (!is.null(labelJitteredPoints))
		{
			p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 6/.pt, position = position_jitter(width = 0.2, height = 0), color = jitterColors)
		} else {
			if (includeJitters){
				p = p + geom_jitter(width=.1, color = jitterColors)	
			}
		}
	}
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo, label = "p.format")}
	if (!is.null(hlines_coo))
	{
		p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" ) + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1), color = "gray44",size = 6/.pt)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	p = p + stat_summary(fun=median, geom="point", size=.1, color=xColors,position = position_dodge(0.9) )
	if (!is.null(legend_position)){
		p = p + theme(legend.position=legend_position)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}


dan.barplots.multipages = function( x, y, fill = NULL, sd = NULL, xlab = "default", ylab = "default", ylimLeft = NULL, ylimRight = NULL, filllab = "default", fillColors = NULL, labelBarsBottom = NULL, plotTitle = "")
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1))
		# p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
	} else {
		p = ggplot(mapping = aes(y = y, x = factor(x))) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45)) + theme_classic(base_size=6)
	}
	if (!is.null(fillColors))
	{
		p = p + scale_fill_manual(values = fillColors)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(sd))
	{
		p = p + geom_errorbar(aes(ymin=y-sd, ymax=y+sd), size = 0.2, alpha = 0.5, stat="identity", position=position_dodge2(width = 0.9, preserve = "single"))+ theme(text = element_text(size=20))
	}
	if (!is.null(labelBarsBottom))
	{
		p = p + geom_text(label = labelBarsBottom, y = y + 1*sign(y),size = 6/.pt)#, position = position_dodge2(width=0.9), size=2.7) # (0-(max(y)-min(y))/10)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.scatterplots.multipages = function( x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", dotSize = 1, fillColors = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
{
	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	} else {
		p = ggplot(mapping = aes(y = y, x = x)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + theme_classic(base_size=6)
	}
	if (!is.null(fillColors))
	{
		p = p + scale_color_manual( values=fillColors,drop = FALSE )
	}
	p = p + geom_point(size = dotSize, alpha = 1)
	if (!is.null(plotFitLine)){
		p = p + geom_smooth(method=FitLineMethod, color = FitLineColor, se = plotFitLine_se)	
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}


dan.densityPlot = function( fileName, x, grouping, groupinglab = "default", xlab = "default", ylab = "default", show_means = F, show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = "firebrick", fileWidth = 4, fileHeight = 3 )
{

	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = "Density" }
	if (groupinglab=="default") { groupinglab = deparse(substitute(grouping)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	p = ggplot(mapping = aes(x = x, color = grouping)) + geom_density() + scale_color_manual(values = groupingColors, name = groupinglab ) + 
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
			theme_classic(base_size=6) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	if (show_means)
	{
		mu = data.frame(grp = levels(grouping), grp.mean = NA)
		for (g in mu$grp)
		{
			mu[mu$grp==g,"grp.mean"] = mean(x[as.character(grouping)==g])
		}
		p = p + geom_vline(data=mu, aes(xintercept=grp.mean, color=grp), linetype="dashed")
	}
	if (show_medians)
	{
		mu = data.frame(grp = levels(grouping), grp.mean = NA)
		for (g in mu$grp)
		{
			mu[mu$grp==g,"grp.median"] = median(x[as.character(grouping)==g])
		}
		p = p + geom_vline(data=mu, aes(xintercept=grp.median, color=grp), linetype="dashed")
	}
	if (!is.null(xlimLeft))
	{
		p = p + xlim(xlimLeft, xlimRight)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.GOBarPlot = function( fileName, df, plotTitle = "", topN = NULL, topFDR = NULL, withGOid = F, barColor = "steelblue", fileWidth = 4, fileHeight = 5 )
{
	df = df[,c("GO.biological.process.complete","upload_1..FDR.")]
	colnames(df) = c("go","fdr")
	df$fdr = as.numeric(df$fdr)
	df = df[order(df$fdr),]
	if (!(is.null(topN)) ) { 
		if (length(topN)==1) {df = df[1:min(c(topN,nrow(df))),] 
			} else {
				df = df[topN,]
			}
	}
	if (!(is.null(topFDR)) ) { df = df[df$fdr < topFDR,] }
	df = df[order(df$fdr, decreasing = T),]
	if (!(withGOid)) { df$go = substr(df$go, 1, nchar(df$go)-13 ) }
	df$fdr_score = -log10(df$fdr)
	df$go = factor(df$go, levels = as.character(df$go))
	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	p = ggplot(mapping = aes(y = fdr_score, x = go), data = df) + geom_bar(stat="identity", fill = barColor) + scale_y_continuous(position = "right") + ggtitle( plotTitle ) + xlab( "" ) + ylab( "-log10(FDR)" ) + theme_classic(base_size=6) + theme(text = element_text(size=10), axis.text.y = element_text(size = 10)) + coord_flip() 
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.GOBarPlotReadFile = function( fileName, panther_file, plotTitle = "", topN = NULL, topFDR = NULL, withGOid = F, rightToLeft = F, verticalLine_color = "firebrick3", verticalLine_thresh = 0.01, barColor = "gray", fileWidth = 4, fileHeight = 5 )
{
	df = read.table(file = panther_file, sep = "\t", skip = 11, header = T,stringsAsFactors = F)
	df = df[,c("GO.biological.process.complete","upload_1..FDR.")]
	colnames(df) = c("go","fdr")
	df$fdr = as.numeric(df$fdr)
	rownames(df) = c(1:nrow(df))
	if (!(is.null(topN)) ) { 
		if (length(topN)==1) {
			df = df[1:min(c(topN,nrow(df))),] 
			df = df[order(df$fdr),]
			} else {
				topN = topN-10
				df = df[topN,]
			}
	}
	if (!(is.null(topFDR)) ) { df = df[df$fdr < topFDR,] }
	df = df[order(df$fdr, decreasing = T),]
	if (!(withGOid)) { df$go = substr(df$go, 1, nchar(df$go)-13 ) }
	df$fdr_score = -log10(df$fdr)
	df$go = factor(df$go, levels = as.character(df$go))
	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	p = ggplot(mapping = aes(y = fdr_score, x = go), data = df) + geom_bar(stat="identity", fill = barColor, colour="black") + geom_hline(yintercept = -log10(verticalLine_thresh), linetype = "dashed", color = verticalLine_color) + scale_y_continuous(position = "right") + ggtitle( plotTitle ) + xlab( "" ) + ylab( "-log10(FDR)" ) + ylim(0,max(df$fdr_score)) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.line.y = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_text(size = 12), ) + coord_flip() 
	if (rightToLeft) { p = p + scale_y_reverse()}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	ggsave(filename = paste0(substr(fileName,1,nchar(fileName)-4),"_ggsave.pdf"), plot = p, device = "pdf")
	return()
}

dan.msigdb_BarPlotReadFile = function( fileName, msigdb_file, plotTitle = "", topN = NULL, topFDR = NULL, rightToLeft = F, verticalLine_color = "firebrick3", verticalLine_thresh = 0.01, barColor = "gray", fileWidth = 4, fileHeight = 5 )
{
	lines = readLines(con = msigdb_file)
	lines = lines[10:(which(lines=="Gene/Gene Set Overlap Matrix")-3)]
	df = data.frame(t(sapply(strsplit(lines, split = "\t"),c)), stringsAsFactors=F)
	colnames(df) = c("go","x2","x3","x4","x5","x6","fdr")
	df = df[2:nrow(df),c("go","fdr")]
	df$fdr = as.numeric(df$fdr)
	dfgo_processed = tolower(gsub("_"," ",gsub("[^_]*_(.*)", "\\1", df$go)))
	substr(dfgo_processed, 1, 1) = toupper(substr(dfgo_processed, 1, 1))
	df$go = dfgo_processed
	rownames(df) = c(1:nrow(df))
	if (!(is.null(topN)) ) { 
		if (length(topN)==1) {
			df = df[order(df$fdr),]
			df = df[1:min(c(topN,nrow(df))),] 
			} else {
				topN = topN-10
				df = df[topN,]
			}
	}
	if (!(is.null(topFDR)) ) { df = df[df$fdr < topFDR,] }
	df = df[order(df$fdr, decreasing = T),]
	df$fdr_score = -log10(df$fdr)
	maxx = max(df$fdr_score)
	df$go = factor(df$go, levels = as.character(df$go))
	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	p = ggplot(mapping = aes(y = fdr_score, x = go), data = df) + geom_bar(stat="identity", fill = barColor, colour="black") + geom_hline(yintercept = -log10(verticalLine_thresh), linetype = "dashed", color = verticalLine_color) + scale_y_continuous(position = "right", breaks = seq(0,maxx,by = max(2,round(maxx/6)) )) + ggtitle( plotTitle ) + xlab( "" ) + ylab( "-log10(FDR)" ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.line.y = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_text(size = 12), ) + coord_flip() 
	if (rightToLeft) { p = p + scale_y_reverse()}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	# ggsave(filename = paste0(substr(fileName,1,nchar(fileName)-4),"_ggsave.pdf"), plot = p, device = "pdf")
	return()
}

dan.pairedDotPlot = function( fileName, x, y, xlab = "default", ylab = "default", plotTitle = "", labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = 1, jitterColors = "black", labelPoints = NULL, fileWidth = 4, fileHeight = 3, linesPairings = NULL )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
		# cat(xlab,ylab)
	p = ggplot(mapping = aes(y = y, x = x)) + geom_point(pch = 20, size = 7, alpha = 0.5, colour = jitterColors) +
		ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
		theme_classic(base_size=6) + theme(text = element_text(size=14),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo)}
	if (!is.null(labelPoints))
	{
		p = p + geom_text(aes(label=labelPoints),hjust=-0.5, vjust=-0.9, size = 6/.pt, color = jitterColors)
	} else {
		p = p + geom_jitter(width=.1, color = jitterColors)
	}
	if (!(is.null(linesPairings)))
	{
		p = p + stat_summary(fun.y=median, geom="line", aes(group=linesPairings), color = "gray44", size = 0.3)# + stat_summary(fun.y=median, geom="point", pch = "-", size = 16, aes(group=linesPairings), color = "gray44")  linetype="dashed",
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.pairedDotPlot.multipages = function( x, y, xlab = "default", ylab = "default", plotTitle = "", labelLines = NULL, signifTest = NULL, comparisons = NULL, ylimLeft = NULL, ylimRight = NULL, labelycoo = 1, jitterColors = "black", labelPoints = NULL, linesPairings = NULL )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
		# cat(xlab,ylab)
	p = ggplot(mapping = aes(y = y, x = x)) + geom_point(pch = 20, size = 7, alpha = 0.5, colour = jitterColors) +
		ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
		theme_classic(base_size=6) + theme(text = element_text(size=14),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo)}
	if (!is.null(labelPoints))
	{
		p = p + geom_text(aes(label=labelPoints),hjust=-0.5, vjust=-0.9, size = 6/.pt, color = jitterColors)
	} else {
		p = p + geom_jitter(width=.1, color = jitterColors)
	}
	if (!(is.null(linesPairings)))
	{
		p = p + stat_summary(fun.y=median, geom="line", aes(group=linesPairings), color = "gray44", size = 0.3)# + stat_summary(fun.y=median, geom="point", pch = "-", size = 16, aes(group=linesPairings), color = "gray44")  linetype="dashed",
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.heatmap = function( fileName, table, topVar = 1, Rowv = NULL, Colv = NULL, scale = 'none', dendrogram = 'none', colors = colorRampPalette(c("blue","white","red"))(n = 100), ColSideColors = NULL, cellnote = NULL, xlab = "", ylab = "", plotTitle = "", key.xlab = "", labRow = "default", labCol = "default", fileWidth = 4, fileHeight = 3, getRearrangedHeatmap = F  )
{
	
	if (topVar<1)
	{
		varS = apply(table,1,function(x) var(x))
		threshValue = quantile(varS, (1-topVar), na.rm = T)
		table = table[varS>=threshValue,]
	}

	if (labRow=='default') { labRow = rownames(table) }
	if (labCol=='default') { labCol = colnames(table) }
	if (is.null(cellnote))
	{
		cellnote = matrix(NA,nrow = nrow(table),ncol = ncol(table))
	}

	if (!is.null(ColSideColors))
	{
		pdf( fileName ,width = fileWidth, height = fileHeight , useDingbats = F )
		a = heatmap.2(as.matrix(table), Rowv = Rowv, Colv = Colv, dendrogram = dendrogram, scale = scale,na.color = "gray", col = colors, tracecol = NA, cellnote = cellnote,notecol="black", ColSideColors = ColSideColors, xlab = xlab, ylab = ylab, main = plotTitle, labRow = labRow, labCol = labCol, cexRow = 0.6,cexCol = 0.6, margins = c(6,10), 
				key.title = NA, key.xlab = key.xlab, key.ylab = NA, key.ytickfun = NA, symm=F,symkey=F,symbreaks=F,density.info = "none", keysize = 1)#, lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 ), )

	# lmat = cbind(c(0, 0, 0, 3), c(5, 4, 1, 2)), lhei = c(1,1,0.5,7), lwid = c(1,7)		 
		dev.off()
	} else {
		pdf( fileName ,width = fileWidth, height = fileHeight , useDingbats = F )
		a = heatmap.2(as.matrix(table), Rowv = Rowv, Colv = Colv, dendrogram = dendrogram, scale = scale,na.color = "gray", col = colors, tracecol = NA, cellnote = cellnote,notecol="black", xlab = xlab, ylab = ylab, main = plotTitle, labRow = labRow, labCol = labCol, cexRow = 0.6,cexCol = 0.6,symm=F,symkey=F,symbreaks=F, margins = c(8,8))				# key.title = NA, key.xlab = key.xlab, key.ylab = NA, density.info = "none", lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 ), keysize = 2 )
		dev.off()
	}
	if (getRearrangedHeatmap) { 
		tabRea = table[rev(a$rowInd), a$colInd]
		save( tabRea, file = paste0(substr(fileName,1,nchar(fileName)-4),"_ReArranged.RData") ) 
	}
	
	return()
}

dan.ScatterDensity_2Dplot = function( fileName, x, y, fill, fillColors=NULL, dotSize=0.1, xlimLeft=NULL, xlimRight=NULL, ylimLeft=NULL, ylimRight=NULL, plotTitle="", xlab="default", ylab="default", filllab="default", fileWidth = 4, fileHeight = 3){
	library(ggExtra)
	if (xlab=="default") { xlab = gsub(".*\\$","",deparse(substitute(x))) }
	if (ylab=="default") { ylab = gsub(".*\\$","",deparse(substitute(y))) }
	if (filllab=="default") { filllab = gsub(".*\\$","",deparse(substitute(fill))) }
	if (is.null(xlimLeft)){ xlimLeft=min(x) }
	if (is.null(xlimRight)){ xlimRight=max(x) }
	if (is.null(ylimLeft)){ ylimLeft=min(y) }
	if (is.null(ylimRight)){ ylimRight=max(y) }
	plot = ggplot(mapping = aes(x=x, y=y, color=fill)) + geom_point(size=dotSize,alpha=0.8) + ggtitle(plotTitle) + xlab(xlab) + ylab(ylab) + labs(color=filllab) + xlim(xlimLeft,xlimRight) + ylim(ylimLeft,ylimRight) + theme_classic(base_size=6) + theme(legend.position="bottom") + coord_fixed()
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	if (!is.null(fillColors)){ plot = plot + scale_color_manual( values=fillColors,drop = FALSE ) }
	plot2 = ggMarginal(plot,type="density",groupFill = TRUE, groupColour = TRUE)
	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	print(plot2)
	dev.off()
}


dan.ScatterDensity_2Dplot_2fills = function( fileNameRoot, x, y, fill1, fillColors1=NULL, filllab1="default", fill2=NULL, fillColors2=NULL, filllab2="default",splitByFill1=FALSE,splitByFill2=FALSE, dotSize=0.1, xlimLeft=NULL, xlimRight=NULL, ylimLeft=NULL, ylimRight=NULL, plotTitle="", xlab="default", ylab="default", hline=NULL,vline=NULL,fileWidth = 4, fileHeight = 4){
	library(ggExtra)
	if (xlab=="default") { xlab = gsub(".*\\$","",deparse(substitute(x))) }
	if (ylab=="default") { ylab = gsub(".*\\$","",deparse(substitute(y))) }
	if (filllab1=="default") { filllab1 = gsub(".*\\$","",deparse(substitute(fill1))) }
	if (filllab2=="default") { filllab2 = gsub(".*\\$","",deparse(substitute(fill2))) }
	if (is.null(xlimLeft)){ xlimLeft=min(x) }
	if (is.null(xlimRight)){ xlimRight=max(x) }
	if (is.null(ylimLeft)){ ylimLeft=min(y) }
	if (is.null(ylimRight)){ ylimRight=max(y) }
	# plotting with fill1
	plot = ggplot(mapping = aes(x=x, y=y, color=fill1)) + geom_point(stroke=0,size=dotSize,alpha=0.8,) + ggtitle(plotTitle) + xlab(xlab) + ylab(ylab) + labs(color=filllab1) + xlim(xlimLeft,xlimRight) + ylim(ylimLeft,ylimRight) + theme_classic(base_size=6) + theme(legend.position="bottom") + coord_fixed()
	if (!is.null(fillColors1)){ plot = plot + scale_color_manual( values=fillColors1,drop = FALSE ) }
	if (!is.null(hline)){ 
		plot = plot + geom_hline(yintercept=hline,linetype="dashed") + geom_vline(xintercept=vline,linetype="dashed") 
		bottomleft = paste0(round(100*sum( (x<hline) & (y<vline) )/length(x),1),"%")
		topleft = paste0(round(100*sum( (x<hline) & (y>=vline) )/length(x),1),"%")
		topright = paste0(round(100*sum( (x>=hline) & (y>=vline) )/length(x),1),"%")
		bottomright = paste0(round(100*sum( (x>=hline) & (y<vline) )/length(x),1),"%")
		annotations = data.frame( xpos=c(-Inf,-Inf,Inf,Inf),ypos=c(-Inf,Inf,-Inf,Inf), annotateText=c(bottomleft,topleft,bottomright,topright), hjustvar=c(-0.5,-0.5,1.5,1.5), vjustvar=c(-1,2,-1,2))
		plot = plot + geom_text(fontface = "bold",data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),colour='black',size = 6/.pt)
	}
	plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	plot2 = ggMarginal(plot,type="density",groupFill = TRUE, groupColour = TRUE)
	pdf( paste0(fileNameRoot,"_colorBy1",gsub(" ","",filllab1),".pdf"), width = fileWidth, height = fileHeight , useDingbats = F )
	print(plot2)
	dev.off()
	# plotting with fill2
	if (!is.null(fill2)){
		plot = ggplot(mapping = aes(x=x, y=y, color=fill2)) + geom_point(stroke=0,size=dotSize,alpha=0.8) + ggtitle(plotTitle) + xlab(xlab) + ylab(ylab) + labs(color=filllab2) + xlim(xlimLeft,xlimRight) + ylim(ylimLeft,ylimRight) + theme_classic(base_size=6) + theme(legend.position="bottom") + coord_fixed()
		if (!is.null(fillColors2)){ plot = plot + scale_color_manual( values=fillColors2,drop = FALSE ) }
		if (!is.null(hline)){ 
			plot = plot + geom_hline(yintercept=hline,linetype="dashed") + geom_vline(xintercept=vline,linetype="dashed") 
			bottomleft = paste0(round(100*sum( (x<hline) & (y<vline) )/length(x),1),"%")
			topleft = paste0(round(100*sum( (x<hline) & (y>=vline) )/length(x),1),"%")
			topright = paste0(round(100*sum( (x>=hline) & (y>=vline) )/length(x),1),"%")
			bottomright = paste0(round(100*sum( (x>=hline) & (y<vline) )/length(x),1),"%")
			annotations = data.frame( xpos=c(-Inf,-Inf,Inf,Inf),ypos=c(-Inf,Inf,-Inf,Inf), annotateText=c(bottomleft,topleft,bottomright,topright), hjustvar=c(-0.5,-0.5,1.5,1.5), vjustvar=c(-1,2,-1,2))
			plot = plot + geom_text(fontface = "bold",data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),colour='black',size = 6/.pt)
		}
		plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		plot2 = ggMarginal(plot,type="density",groupFill = TRUE, groupColour = TRUE)
		pdf( paste0(fileNameRoot,"_colorBy2",filllab2,".pdf"), width = fileWidth, height = fileHeight , useDingbats = F )
		print(plot2)
		dev.off()
	}
	if (splitByFill1){
		for (stratus in unique(as.character(fill1))){
			thisx = x[as.character(fill1)==stratus]
			thisy = y[as.character(fill1)==stratus]
			thisfill2 = fill2[as.character(fill1)==stratus]
			plot = ggplot(mapping = aes(x=thisx, y=thisy, color=thisfill2)) + geom_point(stroke=0,size=dotSize,alpha=0.8) + ggtitle(plotTitle) + xlab(xlab) + ylab(ylab) + labs(color=filllab2) + xlim(xlimLeft,xlimRight) + ylim(ylimLeft,ylimRight) + theme_classic(base_size=6) + theme(legend.position="bottom") + coord_fixed()
			if (!is.null(fillColors2)){ plot = plot + scale_color_manual( values=fillColors2,drop = FALSE ) }
			if (!is.null(hline)){ 
				plot = plot + geom_hline(yintercept=hline,linetype="dashed") + geom_vline(xintercept=vline,linetype="dashed") 
				bottomleft = paste0(round(100*sum( (thisx<hline) & (thisy<vline) )/length(thisx),1),"%")
				topleft = paste0(round(100*sum( (thisx<hline) & (thisy>=vline) )/length(thisx),1),"%")
				topright = paste0(round(100*sum( (thisx>=hline) & (thisy>=vline) )/length(thisx),1),"%")
				bottomright = paste0(round(100*sum( (thisx>=hline) & (thisy<vline) )/length(thisx),1),"%")
				annotations = data.frame( xpos=c(-Inf,-Inf,Inf,Inf),ypos=c(-Inf,Inf,-Inf,Inf), annotateText=c(bottomleft,topleft,bottomright,topright), hjustvar=c(-0.5,-0.5,1.5,1.5), vjustvar=c(-1,2,-1,2))
				plot = plot + geom_text(fontface = "bold",data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),colour='black',size = 6/.pt)
			}
			plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
			plot2 = ggMarginal(plot,type="density",groupFill = TRUE, groupColour = TRUE)
			pdf( paste0(fileNameRoot,"_only",stratus,"_colorBy2",filllab2,".pdf"), width = fileWidth, height = fileHeight , useDingbats = F )
			print(plot2)
			dev.off()
		}
	}
	if (splitByFill2){
		for (stratus in unique(as.character(fill2))){
			thisx = x[as.character(fill2)==stratus]
			thisy = y[as.character(fill2)==stratus]
			thisfill1 = fill1[as.character(fill2)==stratus]
			plot = ggplot(mapping = aes(x=thisx, y=thisy, color=thisfill1)) + geom_point(stroke=0,size=dotSize,alpha=0.8) + ggtitle(plotTitle) + xlab(xlab) + ylab(ylab) + labs(color=filllab1) + xlim(xlimLeft,xlimRight) + ylim(ylimLeft,ylimRight) + theme_classic(base_size=6) + theme(legend.position="bottom") + coord_fixed()
			if (!is.null(fillColors1)){ plot = plot + scale_color_manual( values=fillColors1,drop = FALSE ) }
			if (!is.null(hline)){ 
				plot = plot + geom_hline(yintercept=hline,linetype="dashed") + geom_vline(xintercept=vline,linetype="dashed") 
				bottomleft = paste0(round(100*sum( (thisx<hline) & (thisy<vline) )/length(thisx),1),"%")
				topleft = paste0(round(100*sum( (thisx<hline) & (thisy>=vline) )/length(thisx),1),"%")
				topright = paste0(round(100*sum( (thisx>=hline) & (thisy>=vline) )/length(thisx),1),"%")
				bottomright = paste0(round(100*sum( (thisx>=hline) & (thisy<vline) )/length(thisx),1),"%")
				annotations = data.frame( xpos=c(-Inf,-Inf,Inf,Inf),ypos=c(-Inf,Inf,-Inf,Inf), annotateText=c(bottomleft,topleft,bottomright,topright), hjustvar=c(-0.5,-0.5,1.5,1.5), vjustvar=c(-1,2,-1,2))
				plot = plot + geom_text(fontface = "bold",data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),colour='black',size = 6/.pt)
			}
			plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
			plot2 = ggMarginal(plot,type="density",groupFill = TRUE, groupColour = TRUE)
			pdf( paste0(fileNameRoot,"_only",stratus,"_colorBy1",filllab1,".pdf"), width = fileWidth, height = fileHeight , useDingbats = F )
			print(plot2)
			dev.off()
		}
	}
}


dan.scatterplot = function( fileName, x, y, fill = NULL, shape = NULL, xlab = "default", ylab = "default", xlimLeft = NULL, xlimRight = NULL, ylimLeft = NULL, ylimRight = NULL, filllab = "default", fillColors = NULL, fillColors_continuous = NULL, plotTitle = "", dotLabels = NULL, dotSize = 1, plotVline = NULL, plotHline = NULL, plotBisector = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, legend_position = NULL, xbreaks = NULL,ybreaks=NULL,fileWidth = 4, fileHeight = 3 )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	} else {
		if (!is.null(shape)){
			p = ggplot(mapping = aes(y = y, x = x)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + theme_classic(base_size=6)	
		} else {
			p = ggplot(mapping = aes(y = y, x = x)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + theme_classic(base_size=6)	
		}
	}
	if (!is.null(dotLabels))
	{
		p = p + geom_point(shape=NA) + geom_text( label = dotLabels, size = 6/.pt, show.legend = F ) + guides(color = guide_legend(override.aes = list(shape = 19, size = 5)))
	} else {
		if (!is.null(shape)){
			p = p + geom_point(aes(shape=shape),size = dotSize)
		} else {
			p = p + geom_point(size = dotSize, alpha = 1) # 0.7
		}
	}
	if (!is.null(plotBisector))
	{
		if (plotBisector) { p = p + geom_abline( slope=1,intercept=0 ) }
	}
	if (!is.null(plotVline))
	{
		p = p + geom_vline( xintercept = plotVline,linetype="dashed",size=0.2 )
	}
	if (!is.null(plotHline))
	{
		p = p + geom_hline( yintercept = plotHline,linetype="dashed",size=0.2 )
	}
	if (!is.null(plotFitLine))
	{
		if (plotFitLine & (is.null(fill))) { p = p + geom_smooth(method=FitLineMethod, color = FitLineColor, se = plotFitLine_se, linewidth=fileWidth/7) }
		if (plotFitLine & (!is.null(fill))) { p = p + geom_smooth(method=FitLineMethod, se = plotFitLine_se, linewidth=fileWidth/7) }
	}
	if (!is.null(fillColors))
	{
		p = p + scale_color_manual( values=fillColors,drop = FALSE )
	}
	if (!is.null(fillColors_continuous))
	{
		p = p + scale_colour_gradientn( colours=fillColors_continuous )
	}
	if (!is.null(xlimLeft))
	{
		p = p + xlim(xlimLeft, xlimRight)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(repel_labels))
	{
		p = p + geom_text_repel(segment.size=0.2,aes(label = repel_labels), box.padding = 0, label.padding = 0.15, min.segment.length = 0,force=15,size = 6/.pt,max.overlaps=20000,show.legend=FALSE )
	}
	p = p+ theme(text = element_text(size=14))
	if (coord_fixed){
		p = p + coord_fixed()
	}
	if (coord_flipped){
		p = p + coord_flip()
	}
	if (!is.null(legend_position)){
		p = p + theme(legend.position = legend_position)
	}
	if (!is.null(xbreaks)){
		dcat( xbreaks )
		p = p + scale_x_break(breaks = xbreaks)
	}
	if (!is.null(ybreaks)){
		dcat( ybreaks )
		p = p + scale_y_break(breaks = ybreaks)
	}
	# p = p + guides(color=guide_legend(ncol=1,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.lineplot = function( fileName, x, y, group = NULL, xlab = "default", ylab = "default", xlimLeft = NULL, xlimRight = NULL, ylimLeft = NULL, ylimRight = NULL, grouplab = "default", groupColors = NULL, lineShowLegend = T, dotColors = NULL, plotTitle = "", fileWidth = 4, fileHeight = 3 )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (grouplab=="default") { grouplab = deparse(substitute(group)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(group))
	{
		if (!is.null(dotColors))
		{
			p = ggplot(mapping = aes(y = y, x = x, group = group)) + geom_line(aes(color=group), show.legend = lineShowLegend) + geom_point(color=dotColors, size = 5) + geom_text(aes(label=group),color = "gray77",hjust=-0.5, vjust=-0.5,size = 6/.pt) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = grouplab ) + theme_classic(base_size=6)
		} else {
			p = ggplot(mapping = aes(y = y, x = x, group = group)) + geom_line(aes(color=group), show.legend = lineShowLegend) + geom_point(aes(color=group)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = grouplab ) + theme_classic(base_size=6)
		}
		
	} else {
		p = ggplot(mapping = aes(y = y, x = x, group = 1)) + geom_line() + geom_point() + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = grouplab ) + theme_classic(base_size=6)
	}
	if (!is.null(groupColors))
	{
		p = p + scale_color_manual( values=groupColors )
	}
	if (!is.null(xlimLeft))
	{
		p = p + xlim(xlimLeft, xlimRight)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	p = p+ theme(text = element_text(size=14))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.barplot = function( fileName, x, y, fill = NULL, sd = NULL, xlab = "default", ylab = "default", ylimLeft = NULL, ylimRight = NULL, filllab = "default", fillColors = NULL, labelBarsBottom = NULL, textOnTop = NULL, noxtext=FALSE, y_horizontalLine = NULL, color_horizontalLine = "gray", plotTitle = "", fileWidth = 4, fileHeight = 3 )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1))
		if (noxtext){ p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position="identity") + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) }
	} else {
		p = ggplot(mapping = aes(y = y, x = factor(x))) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1))
	}
	if (!is.null(fillColors))
	{
		p = p + scale_fill_manual(values = fillColors)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(sd))
	{
		p = p + geom_errorbar(aes(ymin=y-sd, ymax=y+sd), size = 0.2, alpha = 0.5, stat="identity", position=position_dodge2(width = 0.9, preserve = "single"))+ theme(text = element_text(size=20))
	}
	if (!is.null(labelBarsBottom))
	{
		p = p + geom_text(label = labelBarsBottom, y = y + 1*sign(y),size = 6/.pt)#, position = position_dodge2(width=0.9), size=2.7) # (0-(max(y)-min(y))/10)
	}
	if (!is.null(y_horizontalLine)){
		p = p + geom_hline(yintercept=y_horizontalLine,color=color_horizontalLine,linetype='dashed',linewidth=0.3)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}


dan.arrowplot = function( fileName, xStart, yStart, xEnd, yEnd, fillStart = NULL, xlab = "default", ylab = "default", filllab = "default", fillColors_continuous = NULL, fillColors = NULL, xlimLeft = NULL, xlimRight = NULL, ylimLeft = NULL, ylimRight = NULL, dotSize = 3, plotTitle = "", fileWidth = 4, fileHeight = 3 )
{
   if (xlab=="default") { xlab = deparse(substitute(x)) }
   if (ylab=="default") { ylab = deparse(substitute(y)) }
   if (filllab=="default") { filllab = deparse(substitute(fillStart)) }

   d = data.frame(xs = xStart, ys = yStart, xe = xEnd, ye = yEnd)

   pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
   if (!is.null(fill))
   {
      # p = ggplot(mapping = aes(y = y, x = x, group = group, color = fill)) + geom_line(aes(color=group), show.legend = lineShowLegend) + geom_point(size = dotSize, alpha = 0.7) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
      p = ggplot(mapping=aes(color=fillStart) ) + geom_point(data=d, mapping=aes(x=xs, y=ys), size=dotSize, alpha = 0.7) + geom_segment(data=d, mapping=aes(x=xs, y=ys, xend=xe, yend=ye), lineend = "round", linejoin = "round", arrow=arrow(length = unit(0.05, "inches")), size=0.2, color="black") + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
   }
   if (!is.null(fillColors))
   {
      p = p + scale_color_manual( values=fillColors,drop = FALSE )
   }
   if (!is.null(fillColors_continuous))
   {
      p = p + scale_colour_gradientn( colours=fillColors_continuous )
   }
   if (!is.null(xlimLeft))
   {
      p = p + xlim(xlimLeft, xlimRight)
   }
   if (!is.null(ylimLeft))
   {
      p = p + ylim(ylimLeft, ylimRight)
   }
   p = p+ theme(text = element_text(size=14))
   p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
   print(p)
   dev.off()
   return()
}


dan.CellCycleAnalysis = function(seu, OutDir){
	dir.create(OutDir)
	seu = CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
	tseu = seu
	md = tseu@meta.data
	### Cells in the s/g2m coordinates
	fileName = paste0(OutDir, "/S_G2M_coordinates_ScatterPlot.pdf")
	x = md$S.Score
	y = md$G2M.Score
	fill = factor(md$Phase, levels = c("G1","G2M","S"))
	dan.scatterplot( fileName, x, y, fill, xlab = "S score", ylab = "G2M score", filllab = "Phase", 
		fillColors = c("lightskyblue3","mediumorchid","orangered3" ), dotSize = 5000/length(x), fileWidth = 7, fileHeight = 7 )
	### PCA on cell cycle features only
	tseu = RunPCA(tseu, assay = 'SCT', features = c(cc.genes$s.genes, cc.genes$g2m.genes))
	pdf(paste0(OutDir,"/PCAcc_ScatterPlot.pdf"))
	print(DimPlot(tseu),reduction = 'pca')
	dev.off()
	### Regress-out everything, and see how it goes 
	tseu = seu
	tseu = SCTransform(
	  tseu,
	  assay = 'RNA',
	  new.assay.name = 'SCT',
	  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score')
	)
	### PCA on cell cycle features only
	tseu = RunPCA(tseu, assay = 'SCT', features = c(cc.genes$s.genes, cc.genes$g2m.genes))
	pdf(paste0(OutDir,"/RegressAll_PCAcc_ScatterPlot.pdf"))
	print(DimPlot(tseu),reduction = 'pca')
	dev.off()
	tseu = seu
	tseu$CC.Difference = tseu$S.Score - tseu$G2M.Score
	tseu = SCTransform(
	  tseu,
	  assay = 'RNA',
	  new.assay.name = 'SCT',
	  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference')
	)
	### PCA on cell cycle features only
	tseu = RunPCA(tseu, assay = 'SCT', features = c(cc.genes$s.genes, cc.genes$g2m.genes))
	pdf(paste0(OutDir,"/RegressDiff_PCAcc_ScatterPlot.pdf"))
	print(DimPlot(tseu),reduction = 'pca')
	dev.off()
}


dan.CellTypeAnalysis = function(seu, OutDir = "CellType_Analysis/", refs = c("BlueprintEncodeData","HumanPrimaryCellAtlasData","ImmGenData" ), 
	SubMainCellTypes = NULL, labelType = "label.main"){
	dir.create(OutDir)
	library(SingleR)
	library(celldex)
	library(pheatmap)
	seuMat = (seu[['SCT']]@data)
	for (r in refs)
	{
		if (r=="BlueprintEncodeData"){ 
			ref = celldex::BlueprintEncodeData()
			if (!is.null(SubMainCellTypes)){ ref = ref[, ref$label.main %in% SubMainCellTypes] }
		}
		if (r=="HumanPrimaryCellAtlasData"){ ref = celldex::HumanPrimaryCellAtlasData() 
			if (!is.null(SubMainCellTypes)){ ref = ref[, ref$label.main %in% SubMainCellTypes] }
		}
		if (r=="ImmGenData"){ ref = celldex::ImmGenData() 
			if (!is.null(SubMainCellTypes)){ ref = ref[, ref$label.main %in% SubMainCellTypes] }
		}
		pred = SingleR(test=seuMat, ref=ref, labels=ref[[labelType]])
		pdf(file = paste0(OutDir,"/Pred_Heatmap_",r,".pdf"), 10,10)
		plotScoreHeatmap(pred)
		dev.off()
		mt = pred[,"labels"]
		names(mt) = rownames(pred)
		seu = AddMetaData(seu,metadata = mt, col.name="CellType")
		tab = table(Assigned=seu@meta.data$CellType, Cluster=seu@meta.data$seurat_clusters)
		pdf(file = paste0(OutDir,"/AssignedClusters_Heatmap_",r,".pdf"), 10,10)
		pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
		dev.off()
		pdf(paste0(OutDir,"/umap_",r,".pdf"),8,8)
		print(DimPlot(seu, group.by = "CellType"))
		dev.off()
	}
}

dan.CellTypeAnalysisSingleR = function(seu, OutPrefix, ref, SubMainCellTypes = NULL, labelType = "label.main"){
	library(SingleR)
	library(celldex)
	library(pheatmap)
	seuMat = (seu[['SCT']]@data)
	r = ref # e.g. "BlueprintEncodeData","HumanPrimaryCellAtlasData"
	if (r=="BlueprintEncodeData"){ 
		ref = celldex::BlueprintEncodeData()
		if (!is.null(SubMainCellTypes)){ ref = ref[, ref$label.main %in% SubMainCellTypes] }
	}
	if (r=="HumanPrimaryCellAtlasData"){ ref = celldex::HumanPrimaryCellAtlasData() 
		if (!is.null(SubMainCellTypes)){ ref = ref[, ref$label.main %in% SubMainCellTypes] }
	}
	if (r=="ImmGenData"){ ref = celldex::ImmGenData() 
		if (!is.null(SubMainCellTypes)){ ref = ref[, ref$label.main %in% SubMainCellTypes] }
	}
	pred = SingleR(test=seuMat, ref=ref, labels=ref[[labelType]])
	predsc = data.frame(pred)
	mt = pred[,"labels"]
	names(mt) = rownames(pred)
	seu = AddMetaData(seu,metadata = mt, col.name="CellType")
	tab = table(Assigned=seu@meta.data$CellType, Cluster=seu@meta.data$seurat_clusters)
	pdf(paste0(OutPrefix,"_SingleR_",r,"_umap.pdf"),16,8)
	print(DimPlot(seu, group.by = "CellType") + coord_fixed())
	dev.off()
	pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = NA, stringsAsFactors = F)
	for (cl in pred$cluster) { 
		tabb = table( predsc[rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,]),"labels"] )
		pred[pred$cluster==cl,"ct"] = names(which.max(tabb))
		pred[pred$cluster==cl,"confidence"] = as.numeric(tabb[names(which.max(tabb))])/sum(tabb) 
	}
	save(pred, file = paste0(OutPrefix,"_pred_SingleR.RData"))
	return(pred)
}
# OutPrefix = "positron_29Jan"
# dan.CellTypeAnalysis_IntegrateHlca(seu,OutPrefix)
dan.CellTypeAnalysis_IntegrateHlca = function(seu, OutPrefix,level = 'all', ct1 = NULL, ct2 = NULL){
	reference = readRDS(paste0("/mnt/ndata/daniele/lung_multiregion/sc/Data/HLCA2022/ref.Rds"))
	if (!is.null(ct1)){
		reference = subset(reference,cells = rownames(reference@meta.data[tolower(reference@meta.data$ann_level_1)==ct1,]))
	}
	if (!is.null(ct2)){
		reference = subset(reference,cells = rownames(reference@meta.data[tolower(reference@meta.data$ann_level_3)==ct2,])) # Macrophages Monocytes
		level = c("level_4","level_5","finest_level" )
	}
	common_features = lapply(list(reference,seu), row.names) %>% Reduce(intersect, .)
	lung.anchors = FindTransferAnchors(reference = reference, query = seu, dims = 1:30, reference.reduction = "refDR", normalization.method="SCT", features = common_features)
	if (level[1]=='all') { level = c( "level_1","level_2","level_3","level_4","level_5","finest_level" ) }
	first = T
	for (levelz in level)
	{
		predictions = TransferData(anchorset = lung.anchors, refdata = reference@meta.data[,paste0("ann_",levelz)], dims = 1:30)
		mt = predictions$predicted.id
		names(mt) = rownames(predictions)
		seu = AddMetaData(seu, metadata = mt, col.name = paste0("ct_",levelz))
		pdf(paste0(OutPrefix,"_Integration_Hlca_",levelz,"_umap.pdf"),16,8)
		print(DimPlot(seu, group.by = paste0("ct_",levelz)) + coord_fixed())
		dev.off()
		pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = NA, stringsAsFactors = F)
		for (cl in pred$cluster) { 
			tabb = table( predictions[rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,]),"predicted.id"] )
			pred[pred$cluster==cl,"ct"] = names(which.max(tabb))
			pred[pred$cluster==cl,"confidence"] = as.numeric(tabb[names(which.max(tabb))])/sum(tabb) 
		}
		pred$method = "IntegrateHlca"
		pred$ref = levelz
		if (first) { 
			pred_all = pred
			first = F
		} else {
			pred_all = rbind(pred_all, pred)
		}
	}
	save(pred_all, file = paste0(OutPrefix,"_pred_IntegrateHlca.RData"))
	return( pred_all )
}

dan.CellTypeAnalysis_ScType = function(seu, OutPrefix){
	lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
	# load gene set preparation function
	source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
	# load cell type annotation function
	source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
	# DB file
	db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
	# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 	
	gs_list = gene_sets_prepare(db_, "Lung")
	gs_list2 = gene_sets_prepare(db_, "Immune system")
	gs_list$gs_positive = c(gs_list$gs_positive,gs_list2$gs_positive)
	gs_list$gs_negative = c(gs_list$gs_negative,gs_list2$gs_negative)
	es.max = sctype_score(scRNAseqData = seu[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
	cL_resutls = do.call("rbind", lapply(unique(seu@meta.data$seurat_clusters), function(cl){
	    es.max.cl = sort(rowSums(es.max[ ,rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
		    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$seurat_clusters==cl)), 10)
		}))
	sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
	sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
	pred = as.data.frame(sctype_scores)
	seu@meta.data$ScType = NA
	for (cl in pred$cluster) { seu@meta.data[seu@meta.data$seurat_clusters==cl,"ScType"] = pred[pred$cluster==cl,"type"] }
	pdf(paste0(OutPrefix,"_ScType_ScType_umap.pdf"),16,8)
	print(DimPlot(seu, group.by = "ScType") + coord_fixed())
	dev.off()
	pred = data.frame(cluster = pred$cluster, ct = pred$type, confidence = 1, stringsAsFactors = F)
	save(pred, file = paste0(OutPrefix,"_pred_ScType.RData"))
	return( pred )
}

dan.CellTypeAnalysis_UCell = function(seu, OutPrefix, marker_list){
	library(UCell)
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/",marker_list,"_markers.txt"))
	for (cn in colnames(aa)){ GeneSets[[paste0(cn)]] = aa[,cn][nchar(aa[,cn])>0] }
	seuNorm = as.matrix(seu[['SCT']]@data)
	GeneSets = GeneSets[sapply(GeneSets,length)<1500]
	scores = ScoreSignatures_UCell(seuNorm, features=GeneSets)	
	scores = data.frame(scores)
	ct = colnames(scores)[apply(scores,1,which.max)]
	predsc = data.frame(row.names = rownames(scores), cell_id = rownames(scores), ct = substr(ct,1,nchar(ct)-6))
	mt = predsc$ct
	names(mt) = rownames(predsc)
	seu = AddMetaData(seu, metadata = mt, col.name = "UCell")
	pdf(paste0(OutPrefix,"_UCell_",marker_list,"_umap.pdf"),16,8)
	print(DimPlot(seu, group.by = "UCell")+ coord_fixed())
	dev.off()
	pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = NA, stringsAsFactors = F)
	for (cl in pred$cluster) { 
		tabb = table( predsc[rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,]),"ct"] )
		pred[pred$cluster==cl,"ct"] = names(which.max(tabb))
		pred[pred$cluster==cl,"confidence"] = as.numeric(tabb[names(which.max(tabb))])/sum(tabb) 
	}
	save(pred, file = paste0(OutPrefix,"_pred_UCell.RData"))
	return( pred )
}

dan.CellTypeAnalysis_AUCell = function(seu, OutPrefix, marker_list){
	library(AUCell)
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/",marker_list,"_markers.txt"))
	for (cn in colnames(aa)){ GeneSets[[paste0(cn)]] = aa[,cn][nchar(aa[,cn])>0] }
	seuNorm = as.matrix(seu[['SCT']]@data)
	cells_rankings = AUCell_buildRankings(seuNorm)
	cells_auc = AUCell_calcAUC(GeneSets, cells_rankings,aucMaxRank = ceiling(5/100 * nrow(cells_rankings)))
	scores = data.frame(t(getAUC(cells_auc)))
	ct = colnames(scores)[apply(scores,1,which.max)]
	predsc = data.frame(row.names = rownames(scores), cell_id = rownames(scores), ct = substr(ct,1,nchar(ct)))
	mt = predsc$ct
	names(mt) = rownames(predsc)
	seu = AddMetaData(seu, metadata = mt, col.name = "AUCell")
	pdf(paste0(OutPrefix,"_AUCell_",marker_list,"_umap.pdf"),16,8)
	print(DimPlot(seu, group.by = "AUCell") + coord_fixed())
	dev.off()
	pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = NA, stringsAsFactors = F)
	for (cl in pred$cluster) { 
		tabb = table( predsc[rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,]),"ct"] )
		pred[pred$cluster==cl,"ct"] = names(which.max(tabb))
		pred[pred$cluster==cl,"confidence"] = as.numeric(tabb[names(which.max(tabb))])/sum(tabb) 
	}
	save(pred, file = paste0(OutPrefix,"_pred_AUCell.RData"))
	return( pred )
}

dan.CellTypeAnalysis_AddModuleScore = function(seu, OutPrefix, marker_list){
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/",marker_list,"_markers.txt"))
	seuNorm = as.matrix(seu[['SCT']]@data)
	for (cn in colnames(aa)){ 
		genes_ok = intersect(aa[,cn][nchar(aa[,cn])>0],rownames(seuNorm))
		if (length(genes_ok)<2){ next }
		GeneSets[[paste0(cn)]] = genes_ok
	}
	seu = AddModuleScore(seu, GeneSets,pool = NULL,nbin = 15,ctrl = 100,k = FALSE,assay = 'SCT',name = "ams")
	scores = seu@meta.data[,substr(colnames(seu@meta.data),1,3)=="ams"]
	colnames(scores) = names(GeneSets)
	ct = colnames(scores)[apply(scores,1,which.max)]
	predsc = data.frame(row.names = rownames(scores), cell_id = rownames(scores), ct = substr(ct,1,nchar(ct)))
	mt = predsc$ct
	names(mt) = rownames(predsc)
	seu = AddMetaData(seu, metadata = mt, col.name = "ams")
	pdf(paste0(OutPrefix,"_ams_",marker_list,"_umap.pdf"),16,8)
	print(DimPlot(seu, group.by = "ams") + coord_fixed())
	dev.off()
	pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = NA, stringsAsFactors = F)
	for (cl in pred$cluster) { 
		tabb = table( predsc[rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,]),"ct"] )
		pred[pred$cluster==cl,"ct"] = names(which.max(tabb))
		pred[pred$cluster==cl,"confidence"] = as.numeric(tabb[names(which.max(tabb))])/sum(tabb) 
	}
	pred = pred[order(as.numeric(pred$cluster)),]
	save(pred, file = paste0(OutPrefix,"_pred_AddModuleScore.RData"))
	return( pred )
}

dan.CellTypeAnalysis_AddModuleScore_xenium = function(seu, OutPrefix, marker_list){
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/",marker_list,"_markers.txt"))
	seuNorm = as.matrix(seu[['SCT']]@data)
	for (cn in colnames(aa)){ 
		genes_ok = intersect(aa[,cn][nchar(aa[,cn])>0],rownames(seuNorm))
		if (length(genes_ok)<2){ next }
		GeneSets[[paste0(cn)]] = genes_ok
	}
	seu = AddModuleScore(seu, GeneSets,pool = NULL,nbin = 10,ctrl = 20,k = FALSE,assay = 'SCT',name = "ams")
	scores = seu@meta.data[,substr(colnames(seu@meta.data),1,3)=="ams"]
	colnames(scores) = names(GeneSets)
	ct = colnames(scores)[apply(scores,1,which.max)]
	predsc = data.frame(row.names = rownames(scores), cell_id = rownames(scores), ct = substr(ct,1,nchar(ct)))
	mt = predsc$ct
	names(mt) = rownames(predsc)
	seu = AddMetaData(seu, metadata = mt, col.name = "ams")
	pdf(paste0(OutPrefix,"_ams_",marker_list,"_umap.pdf"),16,8)
	print(DimPlot(seu, group.by = "ams") + coord_fixed())
	dev.off()
	pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = NA, stringsAsFactors = F)
	for (cl in pred$cluster) { 
		tabb = table( predsc[rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,]),"ct"] )
		pred[pred$cluster==cl,"ct"] = names(which.max(tabb))
		pred[pred$cluster==cl,"confidence"] = as.numeric(tabb[names(which.max(tabb))])/sum(tabb) 
	}
	pred = pred[order(as.numeric(pred$cluster)),]
	save(pred, file = paste0(OutPrefix,"_pred_AddModuleScore.RData"))
	return( pred )
}

dan.CellTypeAnalysis_GSVA = function(seu, OutPrefix, marker_list){
	library(GSVA)
	GeneSets = list()
	aa = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/cell_typing_resources/",marker_list,"_markers.txt"))
	for (cn in colnames(aa)){ GeneSets[[paste0(cn)]] = aa[,cn][nchar(aa[,cn])>0] }
	seuNorm = as.matrix(seu[['SCT']]@data)
	ge = matrix(0, nrow = nrow(seuNorm), ncol = length(unique(seu@meta.data$seurat_clusters)), dimnames = list(rownames(seuNorm),paste0("cl",unique(seu@meta.data$seurat_clusters))) )
	for (cl in unique(seu@meta.data$seurat_clusters)) { 
		cellz = rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl,])
		ge[,paste0("cl",cl )] = rowMeans(seuNorm[,cellz])
	}
	scores = gsva(ge, GeneSets)
	pred = data.frame(cluster = unique(as.character(seu@meta.data$seurat_clusters)), ct = NA, confidence = 1, stringsAsFactors = F)
	for (cl in pred$cluster) { pred[pred$cluster==cl,"ct"] = names(which.max(scores[,paste0("cl",cl)])) }
	save(pred, file = paste0(OutPrefix,"_pred_GSVA.RData"))
	return( pred )
}

dan.CellTypeAnalysis_Full = function(seu, OutPrefix, level = 'all'){
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = NULL, labelType = "label.main")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.main")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))
	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level)
	pred_full = rbind(pred_full,pred_all)
	pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	pred$method = "ScType"
	pred$ref = "ScType"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	dcat("here?")
	for (marker_list in c( "Bischoff","HLCA","Travaglini" )) # Guo
	{
		pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
		pred$method = "ams"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		save(pred_full, file = paste0(OutPrefix,"_pred_full_until_Bischoff.RData"))
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_Full_NoAms = function(seu, OutPrefix, level = 'all'){
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = NULL, labelType = "label.main")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.main")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))
	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level)
	pred_full = rbind(pred_full,pred_all)
	pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	pred$method = "ScType"
	pred$ref = "ScType"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	dcat("here?")
	for (marker_list in c( "Bischoff","HLCA","Travaglini" )) # Guo
	{
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		save(pred_full, file = paste0(OutPrefix,"_pred_full_until_Bischoff.RData"))
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_Full_NoScType = function(seu, OutPrefix, level = 'all'){
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = NULL, labelType = "label.main")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.main")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))
	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level)
	pred_full = rbind(pred_full,pred_all)
	# pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	# pred$method = "ScType"
	# pred$ref = "ScType"
	# pred_full = rbind(pred_full,pred)
	# save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	dcat("here?")
	for (marker_list in c( "Bischoff","HLCA","Travaglini" )) # Guo
	{
		pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
		pred$method = "ams"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		save(pred_full, file = paste0(OutPrefix,"_pred_full_until_Bischoff.RData"))
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_Full_SubCellTypes = function(seu, OutPrefix, level = 'all', ct1 = NULL){
	if (ct1=="immune") { SubMainCellTypes = c(  "B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Erythrocytes","HSC","Macrophages","Monocytes","Neutrophils","NK cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Mesangial cells","Skeletal muscle","Smooth muscle") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial cells") }
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	if (ct1=="immune") { SubMainCellTypes = c(  "B_cell","BM","BM & Prog.","CMP","DC","Embryonic_stem_cells","Erythroblast","GMP","HSC_-G-CSF","HSC_CD34+","Macrophage","MEP","Monocyte","Neutrophils","NK_cell","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Platelets","T_cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Smooth_muscle_cells") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial_cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial_cells") }
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.fine")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))

	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,ct1)
	pred_full = rbind(pred_full,pred_all)
	pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	pred$method = "ScType"
	pred$ref = "ScType"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	if (ct1=="epithelial"){
		this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini","Han","HanNoKac" ),"_",ct1))
	} else {
		this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini" ),"_",ct1),"Salcher_neutrophils")
	}
	for (marker_list in this_marker_list ) # Guo
	{
		pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
		pred$method = "ams"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_Full_SubCellTypes_NoAms_NoGSVA = function(seu, OutPrefix, level = 'all', ct1 = NULL){
	if (ct1=="immune") { SubMainCellTypes = c(  "B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Erythrocytes","HSC","Macrophages","Monocytes","Neutrophils","NK cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Mesangial cells","Skeletal muscle","Smooth muscle") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial cells") }
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	if (ct1=="immune") { SubMainCellTypes = c(  "B_cell","BM","BM & Prog.","CMP","DC","Embryonic_stem_cells","Erythroblast","GMP","HSC_-G-CSF","HSC_CD34+","Macrophage","MEP","Monocyte","Neutrophils","NK_cell","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Platelets","T_cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Smooth_muscle_cells") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial_cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial_cells") }
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.fine")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))

	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,ct1)
	pred_full = rbind(pred_full,pred_all)
	pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	pred$method = "ScType"
	pred$ref = "ScType"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	# for (marker_list in paste0(c( "Bischoff","HLCA","Travaglini" ),"_",ct1) ) # Guo
	# {
	# 	# pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
	# 	# pred$method = "ams"
	# 	# pred$ref = marker_list
	# 	# pred_full = rbind(pred_full,pred)
	# 	pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
	# 	pred$method = "GSVA"
	# 	pred$ref = marker_list
	# 	pred_full = rbind(pred_full,pred)
	# }
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}



dan.CellTypeAnalysis_Full_SubCellTypes_NoAms = function(seu, OutPrefix, level = 'all', ct1 = NULL){
	if (ct1=="immune") { SubMainCellTypes = c(  "B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Erythrocytes","HSC","Macrophages","Monocytes","Neutrophils","NK cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Mesangial cells","Skeletal muscle","Smooth muscle") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial cells") }
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	if (ct1=="immune") { SubMainCellTypes = c(  "B_cell","BM","BM & Prog.","CMP","DC","Embryonic_stem_cells","Erythroblast","GMP","HSC_-G-CSF","HSC_CD34+","Macrophage","MEP","Monocyte","Neutrophils","NK_cell","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Platelets","T_cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Smooth_muscle_cells") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial_cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial_cells") }
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.fine")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))

	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,ct1)
	pred_full = rbind(pred_full,pred_all)
	pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	pred$method = "ScType"
	pred$ref = "ScType"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	if (ct1=="epithelial"){
		this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini","Han" ),"_",ct1))
	} else {
		this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini" ),"_",ct1),"Salcher_neutrophils")
	}
	for (marker_list in this_marker_list ) # Guo
	{
		# pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
		# pred$method = "ams"
		# pred$ref = marker_list
		# pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}


dan.CellTypeAnalysis_Full_SubCellTypes_NoScType = function(seu, OutPrefix, level = 'all', ct1 = NULL){
	if (ct1=="immune") { SubMainCellTypes = c(  "B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Erythrocytes","HSC","Macrophages","Monocytes","Neutrophils","NK cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Mesangial cells","Skeletal muscle","Smooth muscle") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial cells") }
	pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
	pred_full$method = "SingleR"
	pred_full$ref = "BlueprintEncodeData"
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_BlueprintEncodeData.RData"))
	if (ct1=="immune") { SubMainCellTypes = c(  "B_cell","BM","BM & Prog.","CMP","DC","Embryonic_stem_cells","Erythroblast","GMP","HSC_-G-CSF","HSC_CD34+","Macrophage","MEP","Monocyte","Neutrophils","NK_cell","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Platelets","T_cells" ) }
	if (ct1=="stroma") { SubMainCellTypes = c(  "Fibroblasts","Smooth_muscle_cells") }
	if (ct1=="epithelial") { SubMainCellTypes = c(  "Epithelial_cells") }
	if (ct1=="endothelial") { SubMainCellTypes = c(  "Endothelial_cells") }
	pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = NULL, labelType = "label.fine")
	pred$method = "SingleR"
	pred$ref = "HumanPrimaryCellAtlasData"
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_both.RData"))

	pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,ct1)
	pred_full = rbind(pred_full,pred_all)
	# pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
	# pred$method = "ScType"
	# pred$ref = "ScType"
	# pred_full = rbind(pred_full,pred)
	# save(pred_full, file = paste0(OutPrefix,"_pred_full_SingleR_ScType.RData"))
	if (ct1=="epithelial"){
		this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini","Han" ),"_",ct1))
	} else {
		this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini" ),"_",ct1),"Salcher_neutrophils")
	}
	for (marker_list in this_marker_list ) # Guo
	{
		pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
		pred$method = "ams"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.ClusterQualityAnalysis = function(seu, OutDir = "ClusterQuality_Analysis/"){
	
	dir.create(OutDir)
	pdf(paste0(OutDir,"Clustering_vs_CellCyclePhase.pdf"),12,5)
	print(DimPlot(seu,label = TRUE, split.by = "Phase")  + NoLegend())
	dev.off()
	metrics =  c("nCount_RNA", "nFeature_RNA", "percent.mt")
	pdf(paste0(OutDir,"Clustering_vs_Covariates.pdf"),12,10)
	FeaturePlot(seu, features = metrics, pt.size = 0.4, order = TRUE,min.cutoff = 'q10',label = TRUE)
	dev.off()

}


dan.ExtractGenes = function(w, n = 100, method = "max_NMF_package"){
   genes = list()
   if (method == "max_NMF_package") # aka 'barkley'
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
         genes[[lf]] = snp0_genes[res[[lf]]]
      }
   }
   return( genes )
}

dan.ExtractGenes2 = function(w, snp0_genes, n = 100, method = "max_NMF_package", exclude_genes = NULL){
   genes = list()
   rownames(w) = snp0_genes
   if (!is.null(exclude_genes)){
   	w = w[!(rownames(w) %in% exclude_genes),]	
   }
   if (method == "max_NMF_package") # aka 'barkley'
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

dan.compare_published_mp = function(wrComm_genes, OutFolder, snp0_genes, printAllGenes = T, printTopIntersections = F){
   dir.create(OutFolder)  
   kinker_mp = list()
   aa = dan.read(file = "data/kinker_mp.txt")
   for (cn in colnames(aa)){
      kinker_mp[[paste0("kinker_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   gavish_mp = list()
   aa = dan.read(file = "data/gavish_mp.txt")
   for (cn in colnames(aa)){
      gavish_mp[[paste0("gavish_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   lungatlas_mp = list()
   barkley_mp = list()
   aa = dan.read(file = "data/barkley_mp.txt")
   for (cn in colnames(aa)){
      barkley_mp[[paste0("barkley_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   cookEMT_mp = list()
   aa = dan.read(file = "data/cookEMT_mp.txt")
   for (cn in colnames(aa)){
   	if ( (grepl("du145" ,cn) | grepl("mcf7",cn)) | (grepl("GOBP" ,cn) | grepl("ovca420",cn)) ){ next }
      cookEMT_mp[[paste0("cookEMT_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   travaglini_mp = list()
   aa = dan.read(file = "data/cell_typing_resources/Travaglini_epithelial_markers.txt")
   for (cn in colnames(aa)){
      travaglini_mp[[paste0("trav_",tolower(gsub("_CELL","",cn)) )]] = aa[,cn][nchar(aa[,cn])>0]
   }
   hlca_mp = list()
   aa = dan.read(file = "data/cell_typing_resources/HLCA_epithelial_markers.txt")
   for (cn in colnames(aa)){
   	hlca_mp[[paste0("hlca_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   han_mp = list()
   aa = dan.read(file = "data/cell_typing_resources/Han_epithelial_markers.txt")
   for (cn in colnames(aa)){
   	# if (cn %in% c( "KRAS_signature","KAC_signature" )){ next }
   	han_mp[[paste0("han_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   he_mp = list() # from He suppl table, C1 tab
   aa = dan.read(file = "data/cell_typing_resources/HeFetal_epithelial_markers.txt")
   for (cn in colnames(aa)){
   	he_mp[[paste0("heFetal_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   cao_mp = list() # from Cao 2020 https://www.science.org/doi/full/10.1126/science.aba7721 suppl table 4. Filtered for FDR<0.05 as in their methods.
   aa = dan.read(file = "data/cell_typing_resources/CaoFetal_markers_degs.txt")
   aa = aa[aa$qval<0.05,]
   aa$max.cluster = gsub(" ",".",aa$max.cluster )
   aa$gene_short_name = gsub("'","",aa$gene_short_name )
   unique_cts = unique(aa$max.cluster)
   for (cn in unique_cts){
   	genez = aa[aa$max.cluster==cn,"gene_short_name"]
   	if (length(genez)<10){ next }
   	cao_mp[[paste0("caoFetal_",cn)]] = genez
   }
   cm_mp = list() # from CellMarkers2.0
   aa = dan.read(file = "data/cell_typing_resources/CellMarkers_markers.txt")
   aa = aa[aa$cell_type=="Normal cell",]
   aa$cell_name = gsub(" ",".",aa$cell_name )
   unique_cts = unique(aa$cell_name)
   for (cn in unique_cts){
   	genez = aa[aa$cell_name==cn,"Symbol"]
   	if (length(genez)<10){ next }
   	cm_mp[[paste0("cm_",cn)]] = genez
   }
   
   load( "data/marjnmf_mp.RData" )
   load( "data/marjcluComplete_mp.RData" )
   load( "data/marjcluTop100_mp.RData" )
   aext = readRDS("data/cell_typing_resources/fetalAtlas_eleni/epithelial_subset_sc_de_cell_types_positive_markers_filtered_present_LUAD_extended.rds")
   names(aext) = paste0( "ele_",names(aext) )

   all_mp = c(kinker_mp,gavish_mp,lungatlas_mp,barkley_mp,cookEMT_mp,travaglini_mp,hlca_mp,han_mp,marjnmf_mp,marjcluComplete_mp,marjcluTop100_mp,he_mp,cao_mp,cm_mp,aext)
   load(paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoS_SST/ClassicalSignature_29Oct_3PatientsSplitted/GeneSets/GeneSets_qval0.1_logfc1.RData"))
   all_mp[["lepidic_classic"]] = GeneSets[["lepidic_up"]]
   all_mp[["solid_classic"]] = GeneSets[["solid_up"]]
   load(paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoS_SST/ClassicalSignature_29Oct_3PatientsSplitted/GeneSets/gsaug_10.RData"))
   all_mp[["lepidic_augm"]] = gsaug[["lepidic_up"]]
   all_mp[["solid_augm"]] = gsaug[["solid_up"]]

   # Removing genes from mp not in snp0
   for (n in names(all_mp) )
   {
      all_mp[[n]] = intersect(all_mp[[n]],snp0_genes)
   }
   pvalz_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   jaccard_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   intersection_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   union_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   universe = length(snp0_genes)
   colnames(pvalz_df) = names(wrComm_genes)
   rownames(pvalz_df) = names(all_mp)
   for (cn in colnames(pvalz_df))
   {
   	overlap_list = list()
      for (rn in rownames(pvalz_df))
      {
         mymp_genes = length(wrComm_genes[[cn]])
         overlap = length(intersect(all_mp[[rn]],wrComm_genes[[cn]]))
         overlap_list[[ rn ]] = intersect(all_mp[[rn]],wrComm_genes[[cn]])
         othermp_genes = length(all_mp[[rn]])
         pvalz_df[rn,cn] = phyper(overlap-1, mymp_genes, universe-mymp_genes, othermp_genes, lower.tail = FALSE, log.p = FALSE)
         # pvalz_fisher_df[rn,cn] = fisher.test(rbind(c(overlap,othermp_genes-overlap),c(mymp_genes-overlap,universe-mymp_genes- othermp_genes+overlap)), alternative = "greater")$p.value
         intersection_df[rn,cn] = overlap
         union_df[rn,cn] = length(union(all_mp[[rn]],wrComm_genes[[cn]]))
         jaccard_df[rn,cn] = length(intersect(all_mp[[rn]],wrComm_genes[[cn]]))/length(union(all_mp[[rn]],wrComm_genes[[cn]]))
      }
      save(overlap_list,file=paste0(OutFolder,cn,"_vs_PublishedMps_overlaps.RData") )
   }
   for (cn in colnames(pvalz_df))
   {
      thisdf = data.frame(mp = rownames(pvalz_df), minusLog10Pvalz = -log10(pvalz_df[,cn]), qvals = p.adjust(pvalz_df[,cn]), jaccard = jaccard_df[,cn], label = paste0(intersection_df[,cn],"/",union_df[,cn]), intersection_n = intersection_df[,cn] )
      save(thisdf,file=paste0(OutFolder,cn,"_vs_PublishedMps_All.RData"))
      thisdf = thisdf[thisdf$qvals<0.01,]
      thisdf = thisdf[thisdf$jaccard>0,]
      thisdf = thisdf[thisdf$intersection_n>2,]
      thisdf = thisdf[order(-thisdf$minusLog10Pvalz),]
      thisdf$mp = factor(thisdf$mp,levels = as.character(thisdf$mp))
      plotTitle = paste0(cn,", N genes = ",length(wrComm_genes[[cn]]))
      splittedString = wrComm_genes[[cn]]
      for (index in 1:length(splittedString)){ if (index%%15==0) splittedString[index] = paste0(splittedString[index],"\n") }
      if (printAllGenes) { plotTitle = paste0(plotTitle, "\n",paste(splittedString,collapse=", "))  }
      pdf(paste0(OutFolder,cn,"_vs_PublishedMps.pdf"),max(c(4,ceiling(nrow(thisdf)*0.6))),4)
      p1 = ggplot(data=thisdf, aes(x=mp, y=minusLog10Pvalz)) + geom_bar(stat="identity")+ xlab("") + ylab("\n\n\n\n-log10(p-value), hypergeometric test") + ggtitle( plotTitle ) +
         geom_text(aes(label=label), vjust=1.6, color="white", size = 6/.pt)+geom_hline(yintercept=2, colour="red")+geom_label(aes(y=2, label="p-value = 0.01", x=1,hjust=0,vjust=0,label.size = 6/.pt), colour="red", fill = "white")+theme_minimal() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 10))
      p1 = p1 + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
      print(p1)
      dev.off()
      if (printTopIntersections){
      	dcat(cn)
      	dcat( paste0("Intersection with ",thisdf[1,"mp"] ),1 )
      	dcat(sort(intersect(all_mp[[as.character(thisdf$mp)[1]]],wrComm_genes[[cn]])),1)
      	dcat( paste0("Intersection with ",thisdf[2,"mp"] ),1 )
      	dcat(sort(intersect(all_mp[[as.character(thisdf$mp)[2]]],wrComm_genes[[cn]])),1)
      	dcat( paste0("Intersection with ",thisdf[3,"mp"] ),1 )
      	dcat(sort(intersect(all_mp[[as.character(thisdf$mp)[3]]],wrComm_genes[[cn]])),1)
      }
      save(thisdf,file=paste0(OutFolder,cn,"_vs_PublishedMps.RData"))
   }
   return()
}


dan.compare_published_celltypes = function(wrComm_genes, OutFolder, snp0_genes, printAllGenes = T, printTopIntersections = F){
   dir.create(OutFolder)
   travaglini_mp = list()
   aa = dan.read(file = "data/cell_typing_resources/Travaglini_markers.txt")
   for (cn in colnames(aa)){
      travaglini_mp[[paste0("trav_",tolower(gsub("_CELL","",cn)) )]] = aa[,cn][nchar(aa[,cn])>0]
   }
   hlca_mp = list()
   aa = dan.read(file = "data/cell_typing_resources/HLCA_markers.txt")
   for (cn in colnames(aa)){
   	hlca_mp[[paste0("hlca_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   pang_mp = list()
   aa = dan.read(file = "data/cell_typing_resources/PangalaoDb_markers.txt")
   for (cn in colnames(aa)){
   	pang_mp[[paste0("pang_",cn)]] = aa[,cn][nchar(aa[,cn])>0]
   }
   
   all_mp = c(travaglini_mp,hlca_mp,pang_mp)

   # Removing genes from mp not in snp0
   for (n in names(all_mp) )
   {
      all_mp[[n]] = intersect(all_mp[[n]],snp0_genes)
   }
   pvalz_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   jaccard_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   intersection_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   union_df = data.frame(matrix(nrow = length(names(all_mp)), ncol = length(names(wrComm_genes)), dimnames = list(names(all_mp), names(wrComm_genes))))
   universe = length(snp0_genes)
   colnames(pvalz_df) = names(wrComm_genes)
   rownames(pvalz_df) = names(all_mp)
   for (cn in colnames(pvalz_df))
   {
   	overlap_list = list()
      for (rn in rownames(pvalz_df))
      {
         mymp_genes = length(wrComm_genes[[cn]])
         overlap = length(intersect(all_mp[[rn]],wrComm_genes[[cn]]))
         overlap_list[[ rn ]] = intersect(all_mp[[rn]],wrComm_genes[[cn]])
         othermp_genes = length(all_mp[[rn]])
         pvalz_df[rn,cn] = phyper(overlap-1, mymp_genes, universe-mymp_genes, othermp_genes, lower.tail = FALSE, log.p = FALSE)
         # pvalz_fisher_df[rn,cn] = fisher.test(rbind(c(overlap,othermp_genes-overlap),c(mymp_genes-overlap,universe-mymp_genes- othermp_genes+overlap)), alternative = "greater")$p.value
         intersection_df[rn,cn] = overlap
         union_df[rn,cn] = length(union(all_mp[[rn]],wrComm_genes[[cn]]))
         jaccard_df[rn,cn] = length(intersect(all_mp[[rn]],wrComm_genes[[cn]]))/length(union(all_mp[[rn]],wrComm_genes[[cn]]))
      }
      save(overlap_list,file=paste0(OutFolder,cn,"_vs_PublishedMps_overlaps.RData") )
   }
   for (cn in colnames(pvalz_df))
   {
      thisdf = data.frame(mp = rownames(pvalz_df), minusLog10Pvalz = -log10(pvalz_df[,cn]), qvals = p.adjust(pvalz_df[,cn]), jaccard = jaccard_df[,cn], label = paste0(intersection_df[,cn],"/",union_df[,cn]), intersection_n = intersection_df[,cn] )
      thisdf = thisdf[thisdf$qvals<0.01,]
      thisdf = thisdf[thisdf$jaccard>0,]
      thisdf = thisdf[thisdf$intersection_n>2,]
      thisdf = thisdf[order(-thisdf$minusLog10Pvalz),]
      thisdf$mp = factor(thisdf$mp,levels = as.character(thisdf$mp))
      plotTitle = paste0(cn,", N genes = ",length(wrComm_genes[[cn]]))
      splittedString = wrComm_genes[[cn]]
      for (index in 1:length(splittedString)){ if (index%%15==0) splittedString[index] = paste0(splittedString[index],"\n") }
      if (printAllGenes) { plotTitle = paste0(plotTitle, "\n",paste(splittedString,collapse=", "))  }
      pdf(paste0(OutFolder,cn,"_vs_PublishedMps.pdf"),max(c(10,ceiling(nrow(thisdf)*0.6))),10)
      p1 = ggplot(data=thisdf, aes(x=mp, y=minusLog10Pvalz)) + geom_bar(stat="identity")+ xlab("") + ylab("\n\n\n\n-log10(p-value), hypergeometric test") + ggtitle( plotTitle ) +
         geom_text(aes(label=label), vjust=1.6, color="white", size = 6/.pt)+geom_hline(yintercept=2, colour="red")+geom_label(aes(y=2, label="p-value = 0.01", x=1,hjust=0,vjust=0), colour="red", fill = "white")+theme_minimal() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 10))
      p1 = p1 + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
      print(p1)
      dev.off()
      if (printTopIntersections){
      	dcat(cn)
      	dcat( paste0("Intersection with ",thisdf[1,"mp"] ),1 )
      	dcat(sort(intersect(all_mp[[as.character(thisdf$mp)[1]]],wrComm_genes[[cn]])),1)
      	dcat( paste0("Intersection with ",thisdf[2,"mp"] ),1 )
      	dcat(sort(intersect(all_mp[[as.character(thisdf$mp)[2]]],wrComm_genes[[cn]])),1)
      	dcat( paste0("Intersection with ",thisdf[3,"mp"] ),1 )
      	dcat(sort(intersect(all_mp[[as.character(thisdf$mp)[3]]],wrComm_genes[[cn]])),1)
      }
      save(thisdf,file=paste0(OutFolder,cn,"_vs_PublishedMps.RData"))
   }
   return()
}



dan.randomGeneSets = function( p, tgs, nbin=24, ctrl=1000 ){
   expr.levels = rowMeans(p)[order(rowMeans(p))]
   # data.cut = ggplot2::cut_number(x = expr.levels + rnorm(n = length(expr.levels))/1e+30,n = nbin,labels = FALSE,right = FALSE)
   data.cut = rep(1:nbin, each = ceiling(length(expr.levels)/nbin))
   n_excess = length(data.cut)-length(expr.levels)
   ids_toremove = length(data.cut)-c(1:n_excess)*ceiling(length(expr.levels)/nbin)+1
   data.cut = data.cut[-ids_toremove]
   names(data.cut) = names(expr.levels)
   randf = data.frame(matrix(NA,nrow=length(tgs),ncol=ctrl,dimnames = list(tgs,1:ctrl)))
   for (rn in rownames(randf)){
      randf[rn,] = names(sample(which(x = data.cut == data.cut[rn]),size=ctrl,replace=T))      
   }
   random_gs = list()
   for (cn in colnames(randf)){
      this_gs = sort(randf[,cn])
      if (!(paste(this_gs,collapse="_") %in% sapply(random_gs,function(x) paste(x,collapse = "_")) ))
      {
         random_gs[[cn]] = this_gs   
      }
   }
   return( random_gs )
}


dan.overdispersion = function( p, tgs, random_gs ){
   tabb = t(p[tgs,])
   pca = prcomp(tabb)
   var_pc1_real = summary(pca)$importance["Proportion of Variance","PC1"]
   coord_real = var_pc1_real/summary(pca)$importance["Proportion of Variance","PC2"]
   var_pc1_random = c()
   coord_random = c()
   for (rr in names(random_gs))
   {
      tabb = t(p[random_gs[[rr]],])
      pca = prcomp(tabb)
      var_pc1_random = c(var_pc1_random,summary(pca)$importance["Proportion of Variance","PC1"])
      coord_random = c(coord_random,summary(pca)$importance["Proportion of Variance","PC1"]/summary(pca)$importance["Proportion of Variance","PC2"])
   }
   pval_od = sum(var_pc1_random>var_pc1_real)/length(random_gs)
   pval_coord = sum(coord_random>coord_real)/length(random_gs)
   pvals = c(pval_od,pval_coord)
   names(pvals) = c("pval_od","pval_coord")
   return( pvals )
}


dan.df = function(rownames, colnames, data = NA, as_df = T){
	if ((length(rownames)==1)) { mat = matrix(nrow = 0, ncol = length(colnames), dimnames = list(NULL,colnames)) }
	if ((length(rownames)!=1)) { mat = matrix(nrow = length(rownames), ncol = length(colnames), dimnames = list(rownames,colnames), data = data) }
	mat = data.frame(mat, stringsAsFactors = F)
	colnames(mat) = colnames
	return( mat )
}

dan.barkley_MakeRand = function(
  seu,
  db,
  nrand = 3,
  nbin = 25,
  seed = 123
){
	set.seed(seed)
  data = as.matrix(seu@assays$SCT@data)
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                        n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  db_rand = lapply(names(db), function(m){
    lapply(1:10^nrand, function(i){
      used = vector()
      unused = binned
      for (g in db[[m]]){
        pool = data.cut[g]
        new = sample(unused[[pool]], 1)
        used = c(used, new)
        unused[[pool]] = setdiff(unused[[pool]], new)
      }
      return(used)
    })
  })
  names(db_rand) = names(db)
  return(db_rand)
}

dan.barkley_MakeRand_data = function(
  data,
  db,
  nrand = 3,
  nbin = 25,
  seed = 123
){
	set.seed(seed)
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                        n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  db_rand = lapply(names(db), function(m){
    lapply(1:10^nrand, function(i){
      used = vector()
      unused = binned
      for (g in db[[m]]){
        pool = data.cut[g]
        new = sample(unused[[pool]], 1)
        used = c(used, new)
        unused[[pool]] = setdiff(unused[[pool]], new)
      }
      return(used)
    })
  })
  names(db_rand) = names(db)
  return(db_rand)
}

dan.Barkley_GeneToEnrichment_AmsCentered = function(
  seu,
  db = NULL,
  db_rand = NULL
){
	 data = as.matrix(seu@assays$SCT@data)
	 data = t(scale(t(data)))
	 db = lapply(db, intersect, rownames(seu))
	   nrand = log10(length(db_rand[[1]]))
	 enrichment.profile = sapply(names(db), function(m){
	   ra = sapply(db_rand[[m]], function(i){
	     colMeans(data[i, ], na.rm = TRUE)
	   })
	   re = colMeans(data[db[[m]], ], na.rm = TRUE)
	   p = re-rowMeans(ra)
	   return(p)
	 })
	 scores = seu@meta.data
	 commonz = intersect(rownames(scores),rownames(enrichment.profile))
	 if (length(commonz)!=nrow(scores)){ dcat("Careful: not the same number of cells" ) }
	 scores = cbind(scores[commonz,],enrichment.profile[commonz,])
    return(scores)
}

dan.Barkley_GeneToEnrichment_AmsCentered_data = function(
  data, annot,
  db = NULL,
  db_rand = NULL
){
	 data = t(scale(t(as.matrix(data) )))
	 db = lapply(db, intersect, rownames(data))
	   nrand = log10(length(db_rand[[1]]))
	 enrichment.profile = sapply(names(db), function(m){
	 	dcat(m)
	   ra = sapply(db_rand[[m]], function(i){
	     colMeans(data[intersect(rownames(data),i), ], na.rm = TRUE)
	   })
	   re = colMeans(data[intersect(rownames(data),db[[m]]), ], na.rm = TRUE)
	   p = re-rowMeans(ra)
	   return(p)
	 })
	 scores = annot
	 commonz = intersect(rownames(scores),rownames(enrichment.profile))
	 if (length(commonz)!=nrow(scores)){ dcat("Careful: not the same number of cells" ) }
	 scores = cbind(scores[commonz,],enrichment.profile[commonz,])
    return(scores)
}

dan.Barkley_GeneToEnrichment_AmsCentered_data_corrected = function(
  data, annot,
  db = NULL,
  db_rand = NULL,
  by = NULL,
  correction_factors = NULL){
	 data = t(scale(t(as.matrix(data) )))
	 db = lapply(db, intersect, rownames(data))
	   nrand = log10(length(db_rand[[1]]))
	 enrichment.profile = sapply(names(db), function(m){
	 	dcat(m)
	   ra = sapply(db_rand[[m]], function(i){
	     colMeans(data[intersect(rownames(data),i), ], na.rm = TRUE)
	   })
	   re = colMeans(data[intersect(rownames(data),db[[m]]), ], na.rm = TRUE)
	   p = re-rowMeans(ra)
	   return(p)
	 })
	 scores = annot
	 commonz = intersect(rownames(scores),rownames(enrichment.profile))
	 if (length(commonz)!=nrow(scores)){ dcat("Careful: not the same number of cells" ) }
	 scores = cbind(scores[commonz,],enrichment.profile[commonz,])
    return(scores)
}

dan.CellTypeAnalysis_level3 = function(seu, OutPrefix, level = 'all', ct2 = NULL){
	if (ct2=="monocytes") { 
		SubMainCellTypes = c(  "Monocyte" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "HLCA","Travaglini" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="macrophages") { 
		SubMainCellTypes = c(  "Macrophage" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "Macrophages" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "Bischoff","HLCA","pooled","CasanovaAcebes" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="neutrophils") { 
		SubMainCellTypes = c(  "Neutrophils" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "Neutrophils" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,"immune",NULL)
		pred_full = rbind(pred_full,pred_all)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "Salcher" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="DC") { 
		SubMainCellTypes = c(  "DC" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "DC" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,"immune",NULL)
		pred_full = rbind(pred_full,pred_all)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "Zilionis" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="Bcells") { 
		SubMainCellTypes = c(  "B_cell" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "B-cells" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in  paste0(c( "Patil","Hao","pooled" ),"_",ct2)) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_level3_noScType = function(seu, OutPrefix, level = 'all', ct2 = NULL){
	if (ct2=="monocytes") { 
		SubMainCellTypes = c(  "Monocyte" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		for (marker_list in paste0(c( "HLCA","Travaglini","CasanovaAcebes" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="macrophages") { 
		SubMainCellTypes = c(  "Macrophage" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "Macrophages" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		for (marker_list in paste0(c( "Bischoff","HLCA","pooled","CasanovaAcebes" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="neutrophils") { 
		SubMainCellTypes = c(  "Neutrophils" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "Neutrophils" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,"immune",NULL)
		pred_full = rbind(pred_full,pred_all)
		for (marker_list in paste0(c( "Salcher" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="Bcells") { 
		SubMainCellTypes = c(  "B_cell" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "B-cells" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		for (marker_list in  paste0(c( "Patil","Hao","pooled" ),"_",ct2)) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="DC") { 
		SubMainCellTypes = c(  "DC" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "DC" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,"immune",NULL)
		pred_full = rbind(pred_full,pred_all)
		for (marker_list in paste0(c( "Zilionis" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
			pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
			pred$method = "GSVA"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_level3_NoGSVA = function(seu, OutPrefix, level = 'all', ct2 = NULL){
	if (ct2=="monocytes") { 
		SubMainCellTypes = c(  "Monocyte" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "HLCA","Travaglini" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="macrophages") { 
		SubMainCellTypes = c(  "Macrophage" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "Macrophages" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "Bischoff","HLCA","pooled" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="Bcells") { 
		SubMainCellTypes = c(  "B_cell" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "B-cells" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		pred$method = "ScType"
		pred$ref = "ScType"
		pred_full = rbind(pred_full,pred)
		for (marker_list in  paste0(c( "Patil","Hao","pooled" ),"_",ct2)) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_level3_NoGSVA_NoScType = function(seu, OutPrefix, level = 'all', ct2 = NULL){
	if (ct2=="monocytes") { 
		SubMainCellTypes = c(  "Monocyte" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		# pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		# pred$method = "ScType"
		# pred$ref = "ScType"
		# pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "HLCA","Travaglini" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="macrophages") { 
		SubMainCellTypes = c(  "Macrophage" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "Macrophages" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		pred_all = dan.CellTypeAnalysis_IntegrateHlca(seu, OutPrefix,level,NULL,ct2)
		pred_full = rbind(pred_full,pred_all)
		# pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		# pred$method = "ScType"
		# pred$ref = "ScType"
		# pred_full = rbind(pred_full,pred)
		for (marker_list in paste0(c( "Bischoff","HLCA","pooled" ),"_",ct2) ) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	if (ct2=="Bcells") { 
		SubMainCellTypes = c(  "B_cell" ) 
		pred_full = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "HumanPrimaryCellAtlasData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred_full$method = "SingleR"
		pred_full$ref = "HumanPrimaryCellAtlasData"
		SubMainCellTypes = c(  "B-cells" ) 
		pred = dan.CellTypeAnalysisSingleR(seu, OutPrefix, ref = "BlueprintEncodeData", SubMainCellTypes = SubMainCellTypes, labelType = "label.fine")
		pred$method = "SingleR"
		pred$ref = "BlueprintEncodeData"
		pred_full = rbind(pred_full,pred)
		# pred = dan.CellTypeAnalysis_ScType(seu, OutPrefix)
		# pred$method = "ScType"
		# pred$ref = "ScType"
		# pred_full = rbind(pred_full,pred)
		for (marker_list in  paste0(c( "Patil","Hao","pooled" ),"_",ct2)) # Guo
		{
			pred = dan.CellTypeAnalysis_AddModuleScore(seu, OutPrefix, marker_list)
			pred$method = "ams"
			pred$ref = marker_list
			pred_full = rbind(pred_full,pred)
		}
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),15,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_Full_xenium = function(seu, OutPrefix, ct1 = NULL, ct2 = NULL){
	pred_full = dan.df(0,c( "cluster","ct","confidence","method","ref" ))
	this_marker_list = "Travaglini"
	if (!is.null(ct1)){ 
		this_marker_list = paste0( this_marker_list,"_",ct1 ) 
		if (ct1=="epithelial"){
			this_marker_list = c(this_marker_list, "HLCAext_epithelial" )
		}
	}
	if (!is.null(ct2)){ this_marker_list = paste0( this_marker_list,"_",ct2 ) }
	for (marker_list in this_marker_list ) # Guo
	{
		pred = dan.CellTypeAnalysis_AddModuleScore_xenium(seu, OutPrefix, marker_list)
		pred$method = "ams"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
	}
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),6,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.CellTypeAnalysis_Full_cosmx = function(seu, OutPrefix, ct1 = NULL, ct2 = NULL){
	pred_full = dan.df(0,c( "cluster","ct","confidence","method","ref" ))
	this_marker_list = c( "Bischoff","HLCA","Travaglini" )
	if (!is.null(ct1)){ 
		if (ct1=="epithelial"){
			this_marker_list = c(paste0(c( "Bischoff","HLCA","HLCAext","Travaglini","Han","HanNoKac" ),"_",ct1))
		} else {
			this_marker_list = c(paste0(c( "Bischoff","HLCA","Travaglini" ),"_",ct1),"Salcher_neutrophils")
		}
	}
	if (!is.null(ct2)){ this_marker_list = paste0( this_marker_list,"_",ct2 ) }
	for (marker_list in this_marker_list ) # Guo
	{
		pred = dan.CellTypeAnalysis_AddModuleScore_xenium(seu, OutPrefix, marker_list)
		pred$method = "ams"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
		pred = dan.CellTypeAnalysis_GSVA(seu, OutPrefix, marker_list)
		pred$method = "GSVA"
		pred$ref = marker_list
		pred_full = rbind(pred_full,pred)
	}
	md = seu@meta.data
	tabb = dtable(md$seurat_clusters,md$CellType_theirpaper)/rowSums(dtable(md$seurat_clusters,md$CellType_theirpaper))
	pred = data.frame(row.names=paste0( "c",as.character(rownames(tabb))),cluster=as.character(rownames(tabb)), ct = colnames(tabb)[apply(tabb,1,which.max)], confidence = apply(tabb,1,max), method="paper",ref="paper" )
	pred_full = rbind(pred_full,pred)
	save(pred_full, file = paste0(OutPrefix,"_pred_full.RData"))
	pdf(file = paste0(OutPrefix,"_summary_barplot.pdf"),6,5,onefile = T)
	for (cl in sort(as.numeric(unique(pred_full$cluster))))
	{
		this = pred_full[pred_full$cluster==cl,]
		this = this[order(this$method),]
		this$ct = paste0(substr(this$ref,1,3),".",this$ct)
		if (any(duplicated(this$ct))) { this$ct = paste0(this$ct,".",1:length(this$ct) ) }
		this$ct = factor(this$ct, levels = as.character(this$ct))
		plot = ggplot(this, aes(x=ct,y=confidence,fill=method)) + geom_bar(stat="identity") + theme_minimal() + xlab("") + ylab(paste0("\n\nConfidence, cl ",cl)) + theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1))
		print(plot)
	}
	dev.off()
	return(pred_full)
}

dan.DimPlot <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  stroke = 0,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = FALSE,
  label.size = 4,
  label.color = 'black',
  label.box = FALSE,
  repel = FALSE,
  alpha = 1,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  ncol = NULL,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
) {
  if (!rlang::is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims > 0L)) {
    abort(message = "'dims' must be a two-length integer vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  # cells <- cells %||% colnames(x = object)
  ##### Cells for all cells in the assay.
  #### Cells function should not only get default layer
  cells <- cells %||% Cells(
    x = object,
    assay = DefaultAssay(object = object[[reduction]])
  )
  # data <- Embeddings(object = object[[reduction]])[cells, dims]
  # data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data <- FetchData(
    object = object,
    vars = c(dims, group.by),
    cells = cells,
    clean = 'project'
  )
  # cells <- rownames(x = object)
  # object[['ident']] <- Idents(object = object)
  # orig.groups <- group.by
  # group.by <- group.by %||% 'ident'
  # data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    split <- FetchData(object = object, vars = split.by, clean=TRUE)[split.by]
    data <- data[rownames(split),]
    data[, split.by] <- split
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- dan.SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        stroke = stroke,
        shape.by = shape.by,
        order = order,
        alpha = alpha,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      plot <- if (is.null(x = orig.groups)) {
        plot + labs(title = NULL)
      } else {
        plot + CenterTitle()
      }
    }
  )
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

dan.SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  stroke = 0,
  shape.by = NULL,
  alpha = 1,
  alpha.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  raster = NULL,
  raster.dpi = NULL
) {
  if ((nrow(x = data) > 1e5) & is.null(x = raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)

  if (!is.null(x = cells.highlight) && pt.size != AutoPointSize(data = data, raster = raster) && sizes.highlight != pt.size && isTRUE(x = raster)) {
    warning("When `raster = TRUE` highlighted and non-highlighted cells must be the same size. Plot will use the value provided to 'sizes.highlight'.")
  }

  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    if (inherits(x = cells.highlight, what = "data.frame")) {
      stop("cells.highlight cannot be a dataframe. ",
           "Please supply a vector or list")
    }
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size,
      raster = raster
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(!is.na(x = data[, col.by]), data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning(
      "Cannot find alpha variable ",
      alpha.by,
      " in data, setting to NULL",
      call. = FALSE,
      immediate. = TRUE
    )
    alpha.by <- NULL
  }

  plot <- ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + geom_scattermore(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      pointsize = pt.size,
      alpha = alpha,
      pixels = raster.dpi
    )
  } else {
    plot + geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      size = pt.size,
      stroke = stroke,
      alpha = alpha
    )
  }
  plot <- plot +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(color = NULL, title = col.by) +
    CenterTitle()
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + cowplot::theme_cowplot()
  return(plot)
}


