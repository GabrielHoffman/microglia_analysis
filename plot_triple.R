
plot_triple = function(SNP, ensGene, peakID,window, fileOut){

	chrom = as.numeric(unique(df_eqtl[Variant == SNP,'Chr']))
	center = as.numeric(unique(df_eqtl[Variant == SNP,'Position']))

	wh = GRanges(paste0("chr", chrom), IRanges(center-window, center+window))

	geneInfo = select(EnsDb.Hsapiens.v75, keys=ensGene, keytype="GENEID", column=c('GENENAME', 'GENEBIOTYPE'))
	Symbol = geneInfo$GENENAME

	edges = c(start(wh), end(wh))

	# AD GWAS
	df = AD_GWAS_summary[CHR == chrom,]
	df = df[(BP >= center-window) & (BP <= center+window),]

	df2 = merge(df, df_ad_finemap[,c('Variant', 'FINEMAP')], by.x="SNP", by.y="Variant", all.x=TRUE)
	df2[,isCandSet :=!is.na(FINEMAP)]

	ymax = max( 5, -log10(min(df$P)))

	fig_AD_GWAS = ggplot(df2, aes(BP, -log10(P), color=!is.na(FINEMAP), size=isCandSet+1.0 )) + geom_point() + scale_x_continuous(limits=edges, expand=c(0,0)) + scale_y_continuous(expand=c(0, 0), limits=c(0,ymax*1.03)) + scale_color_manual(values=c("black", "red")) + ylab(bquote(-log[10]~P)) + scale_size(range=c(.8, 2), breaks=c(1.1))


	# eQTL
	######
	df = df_eqtl[Chr == chrom,]
	df = df[(Position >= center-window) & (Position <= center+window),]
	df = df[Gene==ensGene,]
	df[,isCandSet:=!is.na(PIP)]

	sizeRng = c(.8)
	if( any(df$isCandSet) ){
		sizeRng = append(sizeRng, 2)
	}

	# Get CLPP from eQTL
	clpp_e = df_fm_eQTL[(Variant==SNP)&(Gene==ensGene),data.table(Author=`Consortium/author`, CLPP=PIP.prod)]
	clpp_e = paste(with(clpp_e, paste0('CLPP=', format(CLPP, digits=2), ' ')), collapse="\n")

	if(clpp_e == "CLPP= "){
		clpp_e = ''
	}

	ymax = max( 5, max(df$log10.p.value))
	fig_eQTL = ggplot(df, aes(Position, log10.p.value, color=!is.na(PIP), size=isCandSet+1.0)) + geom_point() + scale_x_continuous(limits=edges, expand=c(0,0)) + scale_y_continuous(expand=c(0, 0), limits=c(0,ymax*1.03)) + scale_color_manual(values=c("black", "red")) + ylab(bquote(-log[10]~P)) + scale_size(range=sizeRng, breaks=c(1.1)) + annotate("text",  x=edges[2], y = ymax, label = clpp_e, vjust=1, hjust=1)

	# caQTL
	#######
	df = df_caQTL[Chr == chrom,]
	df = df[(Position >= center-window) & (Position <= center+window),]
	df = df[Gene==peakID,]

	if( nrow(df) > 0){
		df[,isCandSet:=!is.na(PIP)]

		sizeRng = c(.8)
		if( any(df$isCandSet) ){
			sizeRng = append(sizeRng, 2)
		}

		# Get CLPP from caQTL
		clpp_ca = df_fm_caQTL[(Variant==SNP)&(Gene==peakID),data.table(Author=`Consortium/author`, CLPP=PIP.prod)]
		clpp_ca = paste(with(clpp_ca, paste0('CLPP=', format(CLPP, digits=2), ' ')), collapse="\n")

		if(clpp_ca == "CLPP= "){
			clpp_ca = ''
		}

		ymax = max( 5, max(df$log10.p.value))
		fig_caQTL = ggplot(df, aes(Position, log10.p.value, color=!is.na(PIP), size=isCandSet+1)) + geom_point() + scale_x_continuous(limits=edges, expand=c(0,0)) + scale_y_continuous(expand=c(0, 0), limits=c(0,ymax*1.03)) + scale_color_manual(values=c("black", "red")) + ylab(bquote(-log[10]~P)) + scale_size(range=sizeRng, breaks=c(1.1)) + annotate("text",  x=edges[2], y = ymax, label = clpp_ca, vjust=1, hjust=1)
		}else{
			fig_caQTL = ggplot()
		}
	# pkgload::unload("EnsDb.Hsapiens.v75")
	# library(EnsDb.Hsapiens.v75)

	fig_genebody = plotEnsGenes_gg( EnsDb.Hsapiens.v75, start(wh), end(wh), seqnames(wh), splice_variants=FALSE, non_coding=TRUE)

	# ATAC
	######

	gr_loc = gr_peak_location[findOverlaps(wh, gr_peak_location)@to]

	fig_ATAC = ggplot() + scale_x_continuous(expand=c(0,0), limits=edges) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + scale_color_manual(values = c("black", "red")) + theme_bw() + theme(legend.position = "none", axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.grid=element_blank(), axis.text.y=element_blank(),panel.background = element_blank(), strip.background = element_blank(), rect = element_rect(fill="white", linetype=0))
	if( length(gr_loc) > 0){
		df_loc = as.data.frame(gr_loc)

		col = c("black", "red")
		df_loc$status = 0
		df_loc$status[rownames(df_loc) == peakID] = 1
		df_loc$status = factor(as.character(df_loc$status), c(0, 1))

		col = c()
		if( any(df_loc$status==0)){
			col = append(col, "black")
		}
		if( any(df_loc$status==1)){
			col = append(col, "red")
		}

		df_loc$start = pmax(df_loc$start, edges[1])
		df_loc$end = pmin(df_loc$end, edges[2])


		fig_ATAC = ggplot(df_loc) + geom_segment(aes(x=start, y=1, xend = end, yend=1, color=status), size=4) + scale_x_continuous(expand=c(0,0), limits=edges) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + scale_color_manual(values = col) + theme_bw() + theme(legend.position = "none", axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), panel.grid=element_blank(), axis.text.y=element_blank(),panel.background = element_blank(), strip.background = element_blank(), rect = element_rect(fill="white", linetype=0))
	}

	 # ,
	fig_track = ggbio::tracks( "AD GWAS" = fig_AD_GWAS,
							"eQTL" 		= fig_eQTL,
							"caQTL" 	= fig_caQTL ,
							"Genes"   = fig_genebody ,
							"OCR" 		= fig_ATAC + ggbio::scale_x_sequnit(),
							xlim = wh,
							padding = unit(-0.65, "lines"),
							label.bg.fill="navy", label.text.color="white",
							heights=c(1, 1, 1, .3, .1),
							label.text.cex = c(1,1,1, .8, .1),
							theme = theme_bw(8) + theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank()),
							title=Symbol )
	fig_track@mutable['Genes'] = FALSE
	fig_track@mutable['OCR'] = FALSE
	ggbio::ggsave(fig_track, file=fileOut, width=7, height=7)
}