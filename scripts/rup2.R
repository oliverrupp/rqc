library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(rjson)
library(tximeta)
library(BiocFileCache)
library(DESeq2)
library(gtools)
library(getopt)

min_assigned_read_count <- 10e6
min_sample_correlation <- 0.9

spec = matrix(c(
  'datafolder' , 'd', 1, "character",
  'reads', 'r', 1, 'integer',
  'correlation', 'c', 1, "double",
#  'make', 'm', 0, "",
#  'threads', 't', 1, "integer",
#  'jobs', 'j', 1, "integer",
  'help', 'h', 0, ""
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if(!is.null(opt$datafolder)) setwd(opt$datafolder)
if(!is.null(opt$reads)) min_assigned_read_count <- opt$reads
if(!is.null(opt$correlation)) min_sample_correlation <- opt$correlation

### setwd("~/Projects/rnaseq_ranomics/CR")

outfile="report.pdf"

if(exists("snakemake")) {
  setwd(dirname(snakemake@output[["report"]]))
  min_assigned_read_count = snakemake@params[["reads"]]
  outfile = basename(snakemake@output[["report"]])
}

sample_info <- read.table("reference/samples.tsv", header=F, sep="\t", row.names = 2)
colnames(sample_info) <- "condition"

sname = mixedsort(unique(sample_info$condition))
sample_colors = rainbow(length(sname))
names(sample_colors) = sname






sample_info$files = paste0("results/salmon/", rownames(sample_info), "/quant.sf")
sample_info$names = rownames(sample_info)

bfcDir <- file.path(getwd(),"results/index/BFC")
indexDir <- file.path(getwd(),"results/index/salmon")
fastaPath <- file.path(getwd(),"reference/genome.fa")
gtfPath <- file.path(getwd(),"reference/annotation.gtf")

#bfcloc <- getTximetaBFC()
#message(bfcloc)
#unlink(paste0(bfcloc,"/*"))

setTximetaBFC(bfcDir)
makeLinkedTxome(indexDir=indexDir,
                source="denovo",
                organism=getwd(),
                release="1",
                genome=getwd(),
                fasta=fastaPath,
                gtf=gtfPath)

se <- tximeta(sample_info)
gse <- summarizeToGene(se)

TPM = assay(gse, "abundance")
gene_count_matrix = assay(gse, "counts")
dds <- DESeqDataSet(gse, design = ~condition)
vsd <- vst(dds, blind = FALSE)
cor_data = assay(vsd)
#cor_data = log2(TPM+1)





samples <- row.names(sample_info)

deg_res <- data.frame(sample=character(), quantile=character(), pct=integer())
results <- data.frame(sample=character(), status=character(), reference=character(), reads=integer())

for(sample in samples) {
  message(sample)
  jfile <- paste0("results/trimmed/",sample,".json")
  jdat <- fromJSON(file = jfile)

  total_reads <- jdat$summary$before_filtering$total_reads
  stats <- data.frame(status=c("LowQuality", "Nreads", "TooShort", "TooLong"),
                      reads=c(jdat$filtering_result$low_quality_reads/2,
                              jdat$filtering_result$too_many_N_reads/2,
                              jdat$filtering_result$too_short_reads/2,
                              jdat$filtering_result$too_long_reads/2))
  stats$sample <- sample
  stats$reference <- "Gene"
  results <- rbind(results, stats)

  jfile <- paste0("results/salmon/",sample,"/aux_info/meta_info.json")
  jdat <- fromJSON(file = jfile)

  stats <- data.frame(status=c("Unmapped", "NoFeature", "Assigned"),
                      reads=c(jdat$num_processed-(jdat$num_decoy_fragments+jdat$num_mapped),
                              jdat$num_decoy_fragments, jdat$num_mapped))

  stats$sample <- sample
  stats$reference <- "Gene"
  results <- rbind(results, stats)

  jfile <- paste0("results/salmon_rrna/",sample,"/aux_info/meta_info.json")
  if(file.exists(jfile)) {
    jdat <- fromJSON(file = jfile)

    stats <- data.frame(status=c("Assigned"),
                        reads=c(jdat$num_mapped))

    stats$sample <- sample
    stats$reference <- "rRNA"
    results <- rbind(results, stats)
  }

  jfile <- paste0("results/salmon_quantiles/", sample, "/quant.sf")

  if(file.exists(jfile)) {
    dc <- read.table(jfile, header=T, sep="\t", row.names=1)
    dc$quantile = gsub(".*_(q[0-9]+)$", "\\1", row.names(dc))
    dc$ids = gsub("(.*)_(q[0-9]+)$", "\\1", row.names(dc))

    dc <- dc[,c("ids", "quantile", "NumReads")]

    dcd <- as.data.frame(pivot_wider(dc, names_from = quantile, values_from = NumReads))
    dcd <- column_to_rownames(dcd, var="ids")

    dcd <- dcd[rowSums(dcd) > 0,]

    dcd <- dcd / rowSums(dcd) * 100

    pct = colSums(dcd) / sum(colSums(dcd)) * 100

    deg_res <- rbind(deg_res, data.frame(sample=sample, quantile=names(pct), pct = pct))
  }
}

deg_tmp = pivot_wider(deg_res, names_from = quantile, values_from = pct)

deg_matrix = as.matrix(column_to_rownames(deg_tmp, var="sample"))
deg_matrix = deg_matrix[,paste0("q", 1:10)]



results$status = factor(results$status, levels=rev(c("Assigned", "NoFeature", "Unmapped",
                                                     "Nreads", "TooShort", "TooLong", "LowQuality")))
reads_plot <- ggplot(results, aes(x=sample, y=reads, fill=status)) + geom_bar(stat="identity") +
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) + facet_wrap(~reference, ncol=1)





df = data.frame("G0"=colSums(gene_count_matrix == 0),
                "G1"=colSums(gene_count_matrix >= 1 & gene_count_matrix < 10),
                "G10"=colSums(gene_count_matrix >= 10 & gene_count_matrix < 100),
                "G100"=colSums(gene_count_matrix >= 100 & gene_count_matrix < 1000),
                "G1000"=colSums(gene_count_matrix >= 1000))

df$Replicate = rownames(df)
df$Sample = sample_info[df$Replicate,]$condition

dfm = melt(df, id.vars = c("Replicate", "Sample"))
for (sample_name in unique(sample_info$condition)) {
  message(sample_name)
  if(sum(sample_info$condition == sample_name) > 1) {
  	max_counts <- apply(gene_count_matrix[,sample_info[sample_info$condition == sample_name,]$names], 1, max, na.rm=TRUE)

  	dfm <- rbind(dfm, list(Replicate = sample_name, Sample = sample_name, variable="G0", value=sum(max_counts == 0)))
  	dfm <- rbind(dfm, list(Replicate = sample_name, Sample = sample_name, variable="G1", value=sum(max_counts >= 1 & max_counts < 10)))
  	dfm <- rbind(dfm, list(Replicate = sample_name, Sample = sample_name, variable="G10", value=sum(max_counts >= 10 & max_counts < 100)))
  	dfm <- rbind(dfm, list(Replicate = sample_name, Sample = sample_name, variable="G100", value=sum(max_counts >= 100 & max_counts < 1000)))
  	dfm <- rbind(dfm, list(Replicate = sample_name, Sample = sample_name, variable="G1000", value=sum(max_counts >= 1000)))
  }
}

dfm$variable = as.character(dfm$variable)

dfm[dfm$variable == "G0",]$variable = "No reads"
dfm[dfm$variable == "G1",]$variable = "1 to 10 reads"
dfm[dfm$variable == "G10",]$variable = "10 to 100 reads"
dfm[dfm$variable == "G100",]$variable = "100 to 1000 reads"
dfm[dfm$variable == "G1000",]$variable = "more than 1000 reads"

dfm$variable = factor(dfm$variable, levels=c("No reads", "1 to 10 reads", "10 to 100 reads", "100 to 1000 reads", "more than 1000 reads"))

dfm$Replicate = as.character(dfm$Replicate)
for(s in unique(dfm$Sample)) {
  dfm[dfm$Sample == s,]$Replicate = gsub(paste0(s, "_"), "", dfm[dfm$Sample == s,]$Replicate)
  dfm[dfm$Sample == s,]$Replicate = gsub(paste0(s, "-"), "", dfm[dfm$Sample == s,]$Replicate)
  dfm[dfm$Sample == s,]$Replicate = gsub(s, "S", dfm[dfm$Sample == s,]$Replicate)
}

dfm$Replicate = factor(dfm$Replicate, levels=c(mixedsort(unique(dfm$Replicate)), mixedsort(unique(dfm$Sample))))

gene_coverage_plot <- ggplot(dfm, aes(x=Replicate, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  ylab("Number of Genes") + xlab("Samples") +
  ggtitle("Number of Reads per Gene") +
  guides(fill=guide_legend(title="Number of assigned reads")) +
  theme(text = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  facet_wrap(~Sample, scale="free_x")

sample_cor <- cor(cor_data, method='pearson', use='pairwise.complete.obs')
sample_dist <- as.matrix(dist(t(cor_data)))



remove_by_reads <- results[results$status == "Assigned" & results$reference == "Gene" & results$reads <= min_assigned_read_count,]$sample
remove_by_correlation <- c()

sample_cor_filtered <- sample_cor[! rownames(sample_cor) %in% remove_by_reads, ! colnames(sample_cor) %in% remove_by_reads]

for(sample_name in unique(sample_info$condition)) {
  message(sample_name)

  if(sum(sample_info$condition == sample_name) > 2) {
  sample_replicates <- sample_info[sample_info$condition == sample_name,]$names
  sample_replicates <- sample_replicates[! sample_replicates %in% remove_by_reads]

  if(length(sample_replicates) > 1) {
    rep_num <- sum(sample_info$condition == sample_name)
    low_cor <- rowSums(sample_cor_filtered[sample_replicates, sample_replicates] <= min_sample_correlation)
    if(sum(low_cor/rep_num > 0.67)) {
      low_ids = names(low_cor[low_cor/rep_num > 0.67])
      remove_by_correlation <- c(remove_by_correlation, low_ids)
    }
  }
}
}

low_q_samples <- c(remove_by_reads, remove_by_correlation)

plot_rem <- F
if(length(low_q_samples) > 0) {
  lowq <- data.frame(sample = low_q_samples, reason = "low read number")
  if(length(remove_by_correlation) > 0) {
    lowq[lowq$sample %in% remove_by_correlation,]$reason = "low in-sample correlation"
  }

  removed_plot <- ggplot(lowq, aes(y=sample, x=reason, col=I("red"))) + geom_point(size=4) + # + facet_wrap(~reason, scale="free_x") +
    theme(text = element_text(size = 18)) +
      ylab("Replicates") + xlab(NULL) +
      ggtitle("Low Quality Samples")
  plot_rem <- T
}


message(low_q_samples)

if(length(low_q_samples) < ncol(cor_data)) {
cor_data_filtered = cor_data[,!colnames(cor_data) %in% low_q_samples]

sample_cor_filtered <- cor(cor_data_filtered, method='pearson', use='pairwise.complete.obs')
sample_dist_filtered <- as.matrix(dist(t(cor_data_filtered)))
}

library(PCAtools)


pdf(outfile, w=18, h=12)

theme_set(theme_bw())

print(reads_plot)
print(gene_coverage_plot)

pheatmap(deg_matrix, cluster_rows = T, cluster_cols = F,
         treeheight_row=0, fontsize = 8, border = F,
         heatmap_legend_param=list(title="coverage (%)"))

detach("package:ComplexHeatmap", unload=TRUE)

a = data.frame(samples = sample_info[colnames(cor_data),1])
row.names(a) = colnames(cor_data)

pheatmap(sample_cor, fontsize = 10, annotation_col=a, main="Sample Correlation Heatmap")
pheatmap(sample_dist, fontsize = 10, annotation_col=a, main="Sample Distance Heatmap")

PCA = pca(cor_data, metadata = colData(se))
biplot(PCA, x = "PC2", y = "PC1", colby="condition")

if(plot_rem) { print(removed_plot) }

if(!is.null(sample_cor_filtered)) {
  pheatmap(sample_cor_filtered, fontsize = 10,
           annotation_col=a, main="Sample Correlation Heatmap")

  pheatmap(sample_dist_filtered, fontsize = 10, annotation_col=a, main="Sample Distance Heatmap")

  PCA = pca(cor_data_filtered, metadata = colData(se)[! rownames(colData(se)) %in% low_q_samples,])
  bip = biplot(PCA, x = "PC2", y = "PC1", colby="condition")
  pair = pairsplot(PCA, colby="condition")
  scree = screeplot(PCA, axisLabSize = 18, titleLabSize = 22)

  print(bip)
  print(pair)
  print(scree)
}

dev.off()




