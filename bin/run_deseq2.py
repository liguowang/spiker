#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:46:12 2021

@author: m102324
"""
import sys,os
import subprocess
from optparse import OptionParser

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():

	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile",action="store",type="string",dest="in_file",help="Input data matrix file. First column is gene ID, the first row is sample ID.")
	parser.add_option("-o","--outfile",action="store",type="string",dest="out_file",help="Prefix of the output files.")
	parser.add_option("-g","--group",action="store",type='string', dest="group_id", default=2, help="Comma separated group IDs (such as \"Mut,Mut,WT,WT\" represents four samples of two group). 'The number of group IDs' must match to 'the number of samples' defined in the data matrix file." )
	parser.add_option("-s","--scaling",action="store",type='string',dest="sf",help="Comma separated scaling factors. Such as \"1.1,0.95,1.5,1.3\". The number of scaling factors also need to match to the number of samples defined in the data matrix file.")

	(options,args)=parser.parse_args()

	if not (options.in_file):
		parser.print_help()
		sys.exit(0)

	if not (options.out_file):
		parser.print_help()
		sys.exit(0)

	all_sample_ids = []    #ALL sample IDs, this will be the first row of data matrix file
	all_group_ids = []     #ALL group ID
	all_sfs = [] #ALL scale factors

	if options.group_id:
		all_group_ids = options.group_id.replace(' ','').split(',')
		#all_group_ids = options.group_id.replace(' ','').split(';')

	if options.sf:
		all_sfs = options.sf.replace(' ','').split(',')
		#all_sfs = options.sf.replace(' ','').split(';')

	# get all samples names 
	SAMP = open(options.in_file,'r')
	for l in SAMP:
		l = l.strip()
		f = l.split()
		all_sample_ids = f[1:]
		break
	SAMP.close()

	if options.sf:
		assert len(all_sample_ids) == len(all_group_ids) == len(all_sfs), "Numbe of sample IDs, group IDs, and scaling factors are not equal!"
	else:
		assert len(all_sample_ids) == len(all_group_ids), "Numbe of sample IDs, group IDs are not equal!"

	# get unwanted sample names
	unwanted = [] #unwanted sample name. Must be subset of samples[]
	used_group = []
	used_sf = []

	assert len(used_group) == len(used_sf), "Numbe of group IDs and number of scaling factors are not equal!"

	if options.sf:
		for s,g,f in zip(all_sample_ids, all_group_ids, all_sfs):
			if g == "NA":
				print ("%s\t%s\t(removed)" % (s,g), file=sys.stderr)
				tmp = s.replace('-','.')
				if tmp[0].isdigit():
					tmp = 'X' + tmp
				unwanted.append(tmp)
			else:
				print ("%s\t%s\t%s" % (s,g,f), file=sys.stderr)
				used_group.append(g)
				used_sf.append(f)
				
	else:
		for s,g in zip(all_sample_ids, all_group_ids):
			if g == "NA":
				print ("%s\t%s\t(removed)" % (s,g), file=sys.stderr)
				tmp = s.replace('-','.')
				if tmp[0].isdigit():
					tmp = 'X' + tmp
				unwanted.append(tmp)
			else:
				print ("%s\t%s" % (s,g), file=sys.stderr)
				used_group.append(g)
				

	name = options.out_file 
	ROUT = open(name + '.r','w')
	cwd = os.getcwd()
	
	print ('setwd("%s")' % cwd, file=ROUT)
	print ('library(RColorBrewer)', file=ROUT)
	print ('library(gplots)', file=ROUT)
	print ('library(scales)', file=ROUT)
	print ('require(genefilter)', file=ROUT)
	print ('require(calibrate)', file=ROUT)
	print ('require(RColorBrewer)', file=ROUT)

	r_functions='''
	### Volcano plot function
	
	
	volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", labelsig=FALSE, textcx=1, ...) {
	  res$padj = replace(res$padj, res$pvalue==0, 300);
	  y_min = 0;
	  y_max = max(-log10(res$padj),na.rm=TRUE);
	  x_min = min(res$log2FoldChange,na.rm=TRUE);
	  x_max = max(res$log2FoldChange,na.rm=TRUE);
	  with(subset(res, padj>sigthresh), plot(log2FoldChange, -log10(pvalue), pch=20, col=alpha("grey",0.3), main=main, ylim=c(y_min,y_max),xlim=c(x_min, x_max), xlab="Log2(Fold change)", ylab="-Log10(adjusted p-value)"))
	  with(subset(res, padj<=sigthresh & log2FoldChange <= -lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col=alpha("blue",0.3)))
	  with(subset(res, padj<=sigthresh & log2FoldChange>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col=alpha("red",0.3)))
	  with(subset(res, padj<=sigthresh & abs(log2FoldChange) < lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col=alpha("orange",0.3)))
	  if (labelsig) {
	    require(calibrate)
	    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
	  }
	  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p.adj>",sigthresh,sep=""), paste("p.adj<=", sigthresh, ", ", "LogFC<= -",lfcthresh, sep=""), paste("p.adj<=", sigthresh, ", ", "LogFC>",lfcthresh, sep=""),  paste("p.adj<=", sigthresh, ", ", "|LogFC|<",lfcthresh, sep="")), pch=20, cex=1.2, col=c("grey","blue","red", "orange"))
	}
	
	
	
	### MA plot
	
	maplot <- function (res, lfcthresh=1, sigthresh=0.05, labelsig=TRUE, legendpos="topleft", textcx=1, ...) {
	  tmp = log10(resdata$baseMean+1)
	  x_min = min(tmp, na.rm=TRUE)
	  x_max = max(tmp, na.rm=TRUE)
	  y_min = min(resdata$log2FoldChange, na.rm=TRUE)
	  y_max = max(resdata$log2FoldChange, na.rm=TRUE)
	  
	  X_1 = tmp[resdata$padj <= sigthresh & resdata$log2FoldChange > lfcthresh]
	  Y_1 = resdata$log2FoldChange[resdata$padj <= sigthresh & resdata$log2FoldChange > lfcthresh]
	
	  X_2 = tmp[resdata$padj <= sigthresh & resdata$log2FoldChange <= -lfcthresh]
	  Y_2 = resdata$log2FoldChange[resdata$padj <= sigthresh & resdata$log2FoldChange <= -lfcthresh]
	
	  X_3 = tmp[resdata$padj > sigthresh ]
	  Y_3 = resdata$log2FoldChange[resdata$padj > sigthresh]
	
	  X_4 = tmp[resdata$padj <= sigthresh & abs(resdata$log2FoldChange) <= lfcthresh]
	  Y_4 = resdata$log2FoldChange[resdata$padj <= sigthresh & abs(resdata$log2FoldChange) <= lfcthresh]
	
	  plot(X_1, Y_1, xlim=c(x_min, x_max), ylim=c(y_min, y_max), col=alpha("red",0.4), pch=20, cex=1, xlab="Log10(baseMean)", ylab="Log2(FoldChange)")
	  points(X_2, Y_2,  col=alpha("blue",0.4), pch=20, cex=1)
	  points(X_3, Y_3,  col=alpha("orange",0.3), pch=20, cex=0.5)
	  points(X_4, Y_4,  col=alpha("orange",0.3), pch=20, cex=0.5)
	  
	  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p.adj<=", sigthresh, ", ", "LogFC<= -",lfcthresh, sep=""),  paste("p.adj<=", sigthresh, ", ", "LogFC>",lfcthresh, sep=""), "Others"), pch=20, cex=1.2, col=c("blue","red","orange"))
	}
	
	
	rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
	  rv = rowVars(assay(rld))
	  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	  pca = prcomp(t(assay(rld)[select, ]))
	  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
	  if (is.null(colors)) {
	    if (nlevels(fac) >= 3) {
	      colors = brewer.pal(nlevels(fac), "Paired")
	    }   else {
	      colors = c("black", "red")
	    }
	  }
	  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
	  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
	  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
	  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
	  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
	  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
	  #legend(legendpos, legend=levels(fac), col=colors, pch=20)
	  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
	  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
	  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
	}
	'''

	print (r_functions, file=ROUT)
	print ("\n### read table", file=ROUT)
	print ('countdata <- read.table("%s", header=TRUE, row.names=1)' % options.in_file, file=ROUT)
	print ('ids <- rownames(countdata)', file=ROUT)

	if len(unwanted)>0:
		print ('### remove unwanted columns', file=ROUT)
		print ('countdata <- subset(countdata, select=-c(%s))' % (','.join(unwanted)), file=ROUT)

	#rescale count data
	if len(used_sf) > 0:
		print ('\n### rescale count data frame using spike-in factors', file=ROUT)
		print ('MultVec <- c(%s)' % ','.join(used_sf), file=ROUT)
		print ('countdata <- mapply(FUN = `*`, countdata, MultVec)', file=ROUT)
		print ('rownames(countdata) <- ids', file=ROUT)
		print ('write.table(countdata, file="%s", sep="\\t", quote=FALSE, row.names=TRUE)' % (name + '.scaled_counts.tsv'), file=ROUT)
	else:
		print ('write.table(countdata, file="%s", sep="\\t", quote=FALSE, row.names=TRUE)' % (name + '.counts.tsv'), file=ROUT)

	print ('\n### convert to matrix', file=ROUT)
	print ('countdata <- as.matrix(countdata)', file=ROUT)
	if len(used_sf) > 0:
		print ('countdata[] <- mapply(as.integer, countdata)', file=ROUT)
		print ('rownames(countdata) <- ids', file=ROUT)
	
	print ('\n### show head', file=ROUT)
	print ('head(countdata)', file=ROUT)


	print ('\n### Assign groups', file=ROUT)
	print ('condition <- factor(c(%s))' % ','.join(['"' + i + '"' for i in used_group]), file=ROUT)
	
	
	print ('\n### Run DESeq2 pipleline', file=ROUT)
	print ('library(DESeq2)', file=ROUT)
	print ('(coldata <- data.frame(row.names=colnames(countdata), condition))', file=ROUT)
	print ('dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)', file=ROUT)
	if len(used_sf) > 0:
		print ('### Turn off DESeq scalling', file=ROUT)
		print ('sizeFactors(dds) <-  c(rep(1.0,length(condition)))', file=ROUT)
	print ('dds <- DESeq(dds)', file=ROUT)

	
	print ('\n### Generate dispersion plot', file=ROUT)
	print ('pdf("%s")' % (name + '.dispersion.pdf'), file=ROUT)
	print ('plotDispEsts(dds, main="Dispersion plot")', file=ROUT)
	print ('dev.off()', file=ROUT)
	
	print ('\n### Log transformation', file=ROUT)
	print ('rld <- rlogTransformation(dds)', file=ROUT)
	print ('head(assay(rld))', file=ROUT)
	
	
	print ("\n### Generate samples' distance plot", file=ROUT)
	print ('(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])', file=ROUT)
	print ('sampleDists <- as.matrix(dist(t(assay(rld))))', file=ROUT)
	print ('pdf("%s")' % (name + '.sample_distance.pdf'), file=ROUT)
	print ('heatmap.2(as.matrix(sampleDists), key=F, trace="none",col=colorpanel(100, "black", "white"), ColSideColors=mycols[condition], RowSideColors=mycols[condition],margin=c(10, 10), main="Sample Distance Matrix")', file=ROUT)
	print ('dev.off()', file=ROUT)
	
	
	print ('\n### Generate PCA plot', file=ROUT)
	print ('pdf("%s")' % (name + '.sample_PCA.pdf'), file = ROUT)
	print ('rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))', file=ROUT)
	print ('dev.off()', file=ROUT)
	
	
	print ('\n### Differential gene expression analysis', file=ROUT) 
	print ('res <- results(dds)', file=ROUT)
	print ('table(res$padj<0.05)', file=ROUT)
	print ('res <- res[order(res$padj), ]', file=ROUT)
	print ('resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)', file=ROUT)
	print ('names(resdata)[1] <- "ID"', file=ROUT)
	print ('head(resdata)', file=ROUT)
	print ('write.table(resdata, file="%s", sep="\\t", quote=FALSE, row.names=TRUE)' % (name + '.DE.tsv'), file=ROUT)
	print ('pdf("%s")' % (name + '.pvalue_distribution.pdf'), file=ROUT)
	print ('hist(res$pvalue, breaks=50, col="grey",xlab="P values")', file=ROUT)
	print ('dev.off()', file=ROUT)
	
	
	print ('\n### Generate MA plot', file=ROUT)
	print ('pdf("%s")' % (name + '.MA_plot.pdf'), file=ROUT)
	print ('maplot(resdata, main="MA Plot")', file=ROUT)
	print ('dev.off()', file=ROUT)
	
	
	print ('\n### Generate Volcano plot', file=ROUT)
	print ('pdf("%s")' % (name + '.volcano_plot.pdf'), file=ROUT)
	print ('volcanoplot(resdata, main="Volcano Plot")', file=ROUT)
	print ('dev.off()', file=ROUT)
	
	ROUT.close()
	subprocess.call('Rscript %s' % (name + '.r') ,shell=True)

if __name__=='__main__':
	main()
