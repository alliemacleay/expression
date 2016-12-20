library(matrixStats)
library(reshape2)
library(MASS)
library('ggplot2')
library(ggrepel)  # relies on ggplot2 version 2.2.0
library(plyr)

## ~~ Functions ~~
getGene <- function(str){
  return(strsplit(toString(str), "_")[[1]][1])
}
getGene2 <- function(str){
  return(strsplit(toString(str), "_")[[1]][2])
}

getSamples <- function(str){
  sname <- strsplit(as.character(str), ".", fixed=TRUE)[[1]][4]
  if(substr(sname, 1, 2) == "NA"){
    # R converts "-" to "." in sample names so it truncates NA ids
    sname <- paste(sname, "-", strsplit(as.character(str), ".", fixed=TRUE)[[1]][5], sep="")
  }
  return(sname)
}

getLabel <- function(sample, fusion_info){
  label <- ifelse(dim(fusion_info[fusion_info["ID."]==sample,])[1]==0, '', fusion_info[fusion_info["ID."]==sample,]$Fusion..LP.RP[1])
  # if no label check NA ID too
  label <- ifelse(label == '' & dim(fusion_info[fusion_info["NA."]==sample,])[1]!=0, fusion_info[fusion_info["NA."]==sample,]$Fusion..LP.RP[1], label)
  return(label)
}

fixLabels <- function(df){
  labels <- c()
  for(i in c(1:dim(df)[1])){
    row <- df[i,]
    labels[i] <- ifelse(grepl(paste("^", row['Gene'], "-|-", row['Gene'], "$", sep=''), row['label']), row['label'], '')  # get indeces of labels with this gene in the fusion
  }
  return(labels)
}

getFusionGenes <- function(genes){
  gene_array <- c()
  i <- 1
  for(gs in genes){
    partners <- strsplit(toString(gs), "-")
    gene_array[i] <- partners[[1]][1]
    gene_array[i+1] <- partners[[1]][2]
    i <- i + 2
  }
  gene_array <- gene_array[!is.na(gene_array) & gene_array != '']
  return(unique(gene_array))
}

## ~~ Main ~~
# Read command line arguments
options(echo=TRUE) # to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
setwd("/Users/Admin/Documents/mgh/projects/pipeline_files/cnv/rna_to_dna")
#
# ~~ Data prep - this was run on the cluster to create the coverage file ~~ ##
# 1. create a bed file with positions 3 bases into the exon and 3 bases into the intron after each primer
# python cidtools/tools/multicov_bed_file_generator.py --exome-file gencodeV24lift37.bed --primer-bed targets.gtf --output-file expression_161214.bed
# 2. calculate coverage at these loci
# bedtools multicov -bams `cat rna_bams_list_161214.txt` -bed expression_161214.bed > rna_coverage_121416.txt
#
rna_covfile <- 'rna_coverage_161216.txt'
fusion_info_file <- 'sample_fusion_info.txt' # NA# ID# FISH+ AMP+ Fusion LPInfo Chromosome#(LP) RPInfo Chromosome#(RP) inFrame Notes Type
rna_coverage <- read.delim(rna_covfile, header=TRUE, stringsAsFactors = FALSE)
rna_coverage <- rna_coverage[, colSums(is.na(rna_coverage)) < nrow(rna_coverage)]  # Remove blank columns

fusion_info <- read.delim(fusion_info_file, header=TRUE, stringsAsFactors = FALSE)
fusion_genes <- getFusionGenes(fusion_info$Fusion..LP.RP.)

## Zongli's script used (cdna - gcorr)/gdna
#   where:
#     cnda is 5 minus the next exon end (in the exon)
#     gcorr is 5 plus the next exon end (in the intron)
#     gnda is the 5 before the next adjacent exon start (in the intron)
#  He updated it to (cdna - gcorr) / gcorr
#     * I am now labeling gcorr as gdna
# Arrange data like Zongli
zli <- melt(rna_coverage, id=names(rna_coverage)[1:6])
genes <- unlist(lapply(zli$primer, getGene))
samples <- unlist(lapply(zli$variable, getSamples))
source <- paste(samples, genes, sep="_")

zli$source <- source
cdna <-zli[zli$region == 'exon',] # exon
gdna <- zli[zli$region == 'intron', ]  # intron
if(! (all(cdna$source == gdna$source) & all(cdna$source == gdna$source))){
  paste("Warning! Source columns don't match!")
  quit()
}
cdna.adj <- cdna$value - gdna$value
paste('Negative CDNA values: ', length(cdna.adj[cdna.adj<0]), ' out of ', length(cdna.adj))
cdna.adj[cdna.adj<0] <- 0  # Set negative values to 0
# raw_ratios <- cdna.adj / gdna$value  # UPDATE: This makes no sense - this was in Zongli's original implementation

my_ratios <- cdna.adj/gdna$value  # This makes way more sense

all.3 <- data.frame(source=source,
                    cDNA.consolidated.adj=unlist(cdna$value - gdna$value),
                    gDNA.consolidated=unlist(gdna$value),
                    r2d0 <-unlist(my_ratios),
                    samples <- samples,
                    AP7.gene <- source,
                    Gene <- unlist(lapply(source, getGene2))
)


## ~~ Zongli's Code modified by Allie ~~ ##
all.4 <- all.3  # all.4 will hold all raw data
all.3$r2d <- all.3$r2d0
all.3 <- all.3[!is.na(all.3$r2d), ]  # Remove columns with NA - denominator was 0 (no dna coverage)

##=== summary by gene, median
gex = data.frame('gex' = tapply(all.3$r2d, all.3$AP7.gene, median))  # median by gene
head(gex)
gex$gex[gex$gex <0.001] = 0.001
gex$gex[gex$gex >10000] = 10000
gex$AP7 = sub('_.*', '', row.names(gex))
gex$Gene = sub('.*_', '', row.names(gex))

gex$label <- unlist(lapply(gex$AP7, getLabel, fusion_info))  # labels for fusions in each sample
gex$label <- unlist(fixLabels(gex))  # only label genes the fusion in present in

# Full boxplot
plot0 = ggplot(data=gex, aes_string(x='Gene', y='gex',  color='Gene')) +
  geom_boxplot() +
  theme_bw() +
  geom_jitter(alpha=0.2) +
  geom_text_repel(data=gex[gex$label!="",], aes(label=label), col=1, hjust=runif(1, 0, 1), vjust=0.5, size=3, angle=-30) +
  ylab('Ratio of cDNA-to-DNA coverage') +
  scale_y_log10(limits=c(0.01, 1000), breaks=c(0.01,0.1,1,10,100,1000),labels=c('<= 0.01',0.1, 1, 10, 100, '>=1000')) +
  theme(axis.text.x=element_text(angle=-90)) +
  ggtitle(paste('RNA/DNA expression ratio (n=', length(unique(gex$AP7)), ')', sep=''))

plot0
ggsave(file=paste('RNA_to_DNA_.raw.full.pdf', sep=''), plot=plot0, width = 22.3, height = 8.7)

genelist <- fusion_genes[which(fusion_genes %in% gex$Gene)]
for(igene in genelist){
  print(igene)
  pd <- gex[gex$Gene == igene,]
  plot0 = ggplot(data=pd, aes_string(x='Gene', y='gex',  color='Gene', label='label')) +
    geom_boxplot() +
    theme_bw() +
    geom_jitter() +
    geom_text(col=2, hjust=runif(1, 0, 1), vjust=0.5, size=3) +
    ylab('Ratio of cDNA-to-DNA coverage') +
    scale_y_log10(limits=c(0.01, 100), breaks=c(0.01,0.1,1,10,100),labels=c('<= 0.01',0.1, 1, 10, '>=100')) +
    ggtitle(paste('RNA/DNA expression ratio: ', igene, ' (n=', length(unique(gex$AP7)), ')', sep=''))
  ggsave(file=paste('RNA_to_DNA_', igene, '.raw.pdf', sep=''), plot=plot0, width = 22.3, height = 8.7)
  plot0
  
  nd <- data.frame(gex=log10(pd[order(pd$gex),]$gex), label=pd[order(pd$gex),]$label)
  idx <- grep(paste("^", igene, "-|-", igene, "$", sep=''), nd$label)  # get indeces of labels with this gene in the fusion
  if(length(idx) > 0){
    print(paste(igene, ' has data'))
  }
  nd$label[-idx] <- ''
  nd <- nd[!is.na(nd$gex),]
  fusion = factor(ifelse(nd$label=="", " WT", "fusion"))
  plot <- ggplot(nd, aes(x=c(1:dim(nd)[1]), y=nd$gex, label='label', data=nd$gex)) + 
    geom_point(position="jitter", aes(colour=nd$gex, size=fusion)) +
    geom_text(aes(label=nd$label, vjust=-1.5, hjust=1.3, angle=-25)) + 
    theme_bw() +
    ylab('RNA to DNA ratio, log10') +
    xlab('Sample') +
    scale_color_gradient(limits=c(0,4), low="blue", high="red", guide_legend(title="Expression")) +
    ggtitle(paste('Gene expression of identified fusions:', igene, sep='\n'))
  plot
  ggsave(file=paste('RNA_to_DNA_', igene, '_exp.pdf', sep=''), plot=plot, width=22.3, height=8.7)
}

## ~~ Calculate variance in genes by sample ~~

##=== r2d normalized by housekeeping
head(gex)
house = subset(gex, Gene %in% c('CTBP1', 'B2M', 'GAPDH'))  # Zongli's
house = subset(gex, Gene %in% c('MSMB', 'ERG', 'TFEB'))  # Picked random genes for now
head(house)
house.median = ddply(house, .(AP7), summarize, 'house'=median(gex))  # Calculate median of gex for all housekeeping genes by sample
head(house.median)

gex.house = merge(gex, house.median, by='AP7', all=T)  # add median of housekeeping gene expression for all samples
head(gex.house)
gex.house$gex.h = gex.house$gex / gex.house$house

# Full boxplot normalized by housekeeping genes
pd = gex.house
pd$label = ''
pd$label <- unlist(lapply(pd$AP7, getLabel, fusion_info))  # labels for fusions in each sample
pd$label <- unlist(fixLabels(pd))  # only label genes the fusion in present in
plot1 = ggplot(data=pd, aes_string(x='Gene', y='gex.h',  color='Gene')) +
  geom_boxplot() +
  theme_bw() +
  geom_jitter(alpha=0.2) +
  geom_text_repel(data=pd[gex$label!="",], aes(label=label), col=1, size=3, angle=-30) +
  ylab('Ratio of cDNA-to-DNA coverages, normalized by median of 3 housekeeping genes') +
  scale_y_log10(limits=c(0.01, 100), breaks=c(0.01,0.1,1,10,100),labels=c('<= 0.01',0.1, 1, 10, '>=100')) +
  theme(axis.text.x=element_text(angle=-90)) +
  ggtitle(paste('RNA/DNA expression ratio (n=', length(unique(gex.house$AP7)), '), normalized by median of 3 housekeeping genes', sep=''))

plot1
ggsave(file=paste('RNA_to_DNA_normalized_by_housekeeping.full.pdf', sep=''), plot=plot1, width = 22.3, height = 8.7)

# Plot just for a subset of genes - this panel only include ALK, BRD4, and ROS1 from below
pd = subset(gex.house, Gene %in% c('GAPDH','CTBP1','B2M','BRD4','JAK1','ERBB2','PRKACA','ALK','ROS1'))
plot1 = ggplot(data=pd, aes_string(x='Gene', y='gex.h',  color='Gene', label='label')) +
  geom_boxplot() +
  theme_bw() +
  geom_jitter() +
  geom_text(col=2, hjust=runif(1, 0, 1), vjust=0.5, size=3) +
  ylab('Ratio of cDNA-to-DNA coverages, normalized by median of 3 housekeeping genes') +
  scale_y_log10(limits=c(0.01, 100), breaks=c(0.01,0.1,1,10,100),labels=c('<= 0.01',0.1, 1, 10, '>=100')) +
  ggtitle(paste('RNA/DNA expression ratio, (n=', length(unique(gex.house$AP7)), '), normalized by median of 3 housekeeping genes', sep=''))

ggsave(file=paste('RNA_to_DNA_normalized_by_housekeeping.pdf', sep=''), plot=plot1, width = 11.3, height = 8.7)

##=== RNA only, normalized by housekeeping
rna = data.frame('rna' = tapply(all.4$cDNA.consolidated.adj, all.4$AP7.gene, median))
head(rna)
rna$rna[rna$rna <0.01] = 0.01
rna$AP7 = sub('_.*', '', row.names(rna))
rna$Gene = sub('.*_', '', row.names(rna))


rna$label = unlist(lapply(rna$AP7, getLabel, fusion_info))
rna$label = unlist(fixLabels(rna))
head(rna)
#house = subset(rna, Gene %in% c('CTBP1', 'B2M', 'GAPDH'))  # Zongli's
house = subset(rna, Gene %in% c('MSMB', 'ERG', 'TFEB'))  # Picked random genes for now
head(house)
house.median = ddply(house, .(AP7), summarize, 'house'=median(rna))
head(house.median)

rna.house = merge(rna, house.median, by='AP7', all=T)
head(rna.house)
rna.house$rna.h = rna.house$rna / rna.house$house

pd = rna.house
plotrna = ggplot(data=pd, aes_string(x='Gene', y='rna.h',  color='Gene', label='label')) +
  geom_boxplot() +
  theme_bw() +
  geom_jitter(alpha=0.2) +
  geom_text(col=2, hjust=runif(1, 0, 1), vjust=0.5, size=3) +
  ylab('RNA counts, normalized by median of 3 housekeeping genes') +
  scale_y_log10(limits=c(0.01, 100), breaks=c(0.01,0.1,1,10,100),labels=c('<= 0.01',0.1, 1, 10, '>=100')) +
  theme(axis.text.x=element_text(angle=-90)) +
  ggtitle(paste('RNA expression (n=', length(unique(gex.house$AP7)), '), normalized by median of 3 housekeeping genes', sep=''))
ggsave(file=paste('RNA_normalized_by_housekeeping.full.pdf', sep=''), plot=plotrna, width = 22.3, height = 8.7)

pd = subset(rna.house, Gene %in% c('GAPDH','CTBP1','B2M','BRD4','JAK1','ERBB2','PRKACA','ALK','ROS1'))
plotrna = ggplot(data=pd, aes_string(x='Gene', y='rna.h',  color='Gene', label='label')) +
  geom_boxplot() +
  theme_bw() +
  geom_jitter() +
  geom_text(col=2, hjust=runif(1, 0, 1), vjust=0.5, size=3) +
  ylab('RNA counts, normalized by median of 3 housekeeping genes') +
  scale_y_log10(limits=c(0.01, 100), breaks=c(0.01,0.1,1,10,100),labels=c('<= 0.01',0.1, 1, 10, '>=100')) +
  ggtitle(paste('RNA expression (n=', length(unique(gex.house$AP7)), '), normalized by median of 3 housekeeping genes', sep=''))
ggsave(file=paste('RNA_normalized_by_housekeeping_subset.pdf', sep=''), plot=plotrna, width = 11.3, height = 8.7)
