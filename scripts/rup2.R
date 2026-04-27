options(java.parameters = "-Xmx2048m")

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(GGally))
suppressMessages(library(edgeR))
suppressMessages(library(reshape2))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(rjson))
suppressMessages(library(tximeta))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(BiocFileCache))
suppressMessages(library(DESeq2))
suppressMessages(library(gtools))
suppressMessages(library(getopt))
suppressMessages(library(moments))
suppressMessages(library(PCAtools))
suppressMessages(library(xlsx))
suppressMessages(library(filelock))
suppressMessages(library(scales))

#### setwd("/home/orupp/Projects/trqc/CS")

##################### INIT ####################################################

### DEBUG: rm(list=ls())
### DEBUG: setwd("/vol/ranomics/rnaseq/QC/TT")

min_assigned_read_count <- 100
min_sample_correlation <- 0.9

spec <- matrix(c(
  "datafolder", "d", 1, "character",
  "reads", "r", 1, "integer",
  "correlation", "c", 1, "double",
  "help", "h", 0, ""
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

if (!is.null(opt$datafolder)) setwd(opt$datafolder)
if (!is.null(opt$reads)) min_assigned_read_count <- opt$reads
if (!is.null(opt$correlation)) min_sample_correlation <- opt$correlation

outfile <- "report.html"

if (exists("snakemake")) {
  setwd(dirname(snakemake@output[["report"]]))
  min_assigned_read_count <- snakemake@params[["reads"]]
  outfile <- basename(snakemake@output[["report"]])
}

dir.create("results/rds", recursive = TRUE)

# pdf(outfile, w = 24, h = 16)

theme_set(theme_bw())

###############################################################################





##################### SAMPLES #################################################

sample_info <- read.table("reference/samples.tsv",
                          header = TRUE, sep = "\t", row.names = 2)

sname <- mixedsort(unique(sample_info$condition))
sample_colors <- rainbow(length(sname))
names(sample_colors) <- sname

sample_info$files <-
  paste0("results/salmon/", rownames(sample_info), "/quant.sf")
sample_info$names <- rownames(sample_info)

samples <- row.names(sample_info)

message(sample_info)

###############################################################################





##################### READ COUNTS #############################################

results <- data.frame(sample = character(),
                      status = character(),
                      reference = character(),
                      reads = integer())

for (sample in samples) {
  message(sample)
  jfile <- paste0("results/trimmed/", sample, ".json")
  jdat <- fromJSON(file = jfile)

  type <- jdat$summary$sequencing
  div <- 1

  if (regexpr("paired end", type) >= 0) {
    div <- 2
  }
  total_reads <- jdat$summary$before_filtering$total_reads
  stats <- data.frame(status = c("LowQuality", "Nreads", "TooShort", "TooLong"),
                      reads = c(jdat$filtering_result$low_quality_reads / div,
                                jdat$filtering_result$too_many_N_reads / div,
                                jdat$filtering_result$too_short_reads / div,
                                jdat$filtering_result$too_long_reads / div))
  stats$sample <- sample
  stats$reference <- "Gene"
  results <- rbind(results, stats)

  jfile <- paste0("results/salmon/", sample, "/aux_info/meta_info.json")
  jdat <- fromJSON(file = jfile)

  stats <- data.frame(status = c("Unmapped", "NoFeature", "Assigned"),
                      reads = c(jdat$num_processed -
                                  (jdat$num_decoy_fragments + jdat$num_mapped),
                                jdat$num_decoy_fragments, jdat$num_mapped))

  stats$sample <- sample
  stats$reference <- "Gene"
  results <- rbind(results, stats)

  jfile <- paste0("results/salmon_rrna/", sample, "/aux_info/meta_info.json")
  if (file.exists(jfile)) {
    jdat <- fromJSON(file = jfile)

    stats <- data.frame(status = "rRNA",
                        reads = -1*jdat$num_mapped)

    stats$sample <- sample
    stats$reference <- "Gene"
    results <- rbind(results, stats)
    #selection <- results$status == "Unmapped" & results$sample == sample
    #results[selection, ]$reads <- results[selection, ]$reads  - jdat$num_mapped
  }
}

results$status <- factor(results$status,
                         levels = rev(c("Assigned", "rRNA",
                                        "NoFeature", "Unmapped",
                                        "Nreads", "TooShort",
                                        "TooLong", "LowQuality")))

# reads_plot <- ggplot(results, aes(x = sample, y = reads, fill = status)) +
#   geom_hline(yintercept=c(min_assigned_read_count, 50e6)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(
#     	breaks = function(x) {
#       		default_breaks <- extended_breaks()(x)
#       		sort(c(default_breaks, min_assigned_read_count))
#     	},
#     	labels = function(breaks) {
#      		ifelse(breaks == min_assigned_read_count, "minimal read count", breaks)
#     	}
#   ) + 
#   theme(text = element_text(size = 18)) + ggtitle("Read count analysis") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0))

saveRDS(results, "results/rds/reads_plot.rds")

# print(reads_plot)

###############################################################################

#message(sample_info)
sink(file.path(getwd(), "filtered_samples.txt"))
print(paste0("Low Read Count [<, ", min_assigned_read_count, "]:"))
print(results[results$status == "Assigned" & results$reads < min_assigned_read_count,]$sample)  
sink()

samples <- results[results$status == "Assigned" & results$reads >= min_assigned_read_count,]$sample

sample_info <- sample_info[samples,]



##################### COLLECT COUNTS ##########################################

mtx <- lock("~/.tximeta.lock")

bcf_dir <- file.path(getwd(), "results/index/BFC")
index_dir <- file.path(getwd(), "results/index/salmon")
fasta_path <- file.path(getwd(), "reference/genome.fa")
gtf_path <- file.path(getwd(), "reference/annotation.gtf")

if (!dir.exists(bcf_dir)) {
  dir.create(bcf_dir)
}

setTximetaBFC(bcf_dir)
makeLinkedTxome(indexDir = index_dir,
                source = "denovo",
                organism = getwd(),
                release = "1",
                genome = getwd(),
                fasta = fasta_path,
                gtf = gtf_path)

se <- tximeta(sample_info)
gse <- summarizeToGene(se)

transcript_count_matrix <- assay(se, "counts")
gene_count_matrix <- assay(gse, "counts")

##### edgeR #####

# TMM
dge <- DGEList(counts=gene_count_matrix, group=factor(sample_info[colnames(gene_count_matrix),]$condition))
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
TMM_counts <- cpm(dge, normalized.lib.sizes = TRUE)

tr_dge <- DGEList(counts=transcript_count_matrix, group=factor(sample_info[colnames(transcript_count_matrix),]$condition))
tr_dge <- tr_dge[filterByExpr(tr_dge), , keep.lib.sizes = FALSE]
tr_dge <- calcNormFactors(tr_dge, method = "TMM")
tr_TMM_counts <- cpm(tr_dge, normalized.lib.sizes = TRUE)

# geTMM
tr_rpk <- transcript_count_matrix / assay(se, "length")
gene_rpk <- gene_count_matrix / assay(gse, "length")

tr_len <- assay(se, "length")
tr_len_m <- melt(tr_len)
colnames(tr_len_m) <- c("id", "sample", "length")

# real_len_plot <- ggplot(tr_len_m, aes(x = sample, y = length)) + 
#	geom_boxplot() + scale_y_log10() + 
#	theme(axis.text.x = element_text(angle = 90, hjust = 0))
# print(real_len_plot)

tr_norm_edger <- DGEList(counts=tr_rpk,group=colData(se)$condition)
gene_norm_edger <- DGEList(counts=gene_rpk,group=colData(gse)$condition)

#tr_norm_edger <- tr_norm_edger[filterByExpr(tr_norm_edger), , keep.lib.sizes = FALSE]
#gene_norm_edger <- gene_norm_edger[filterByExpr(gene_norm_edger), , keep.lib.sizes = FALSE]

tr_norm_edger <- calcNormFactors(tr_norm_edger, method = "TMM")
gene_norm_edger <- calcNormFactors(gene_norm_edger, method = "TMM")

tr_geTMM_counts <- cpm(tr_norm_edger)
gene_geTMM_counts <- cpm(gene_norm_edger)

####### edgeR ######


sample_info %>% group_by(condition) %>% summarize(rep_num = n()) -> rep_num
message(rep_num)

#### DESeq2 ####
dds <- DESeqDataSet(gse, design = ~condition)
tr_dds <- DESeqDataSet(se, design = ~condition)
vsd <- vst(dds, blind = sum(rep_num$rep_num > 1) == 0) # blind = TRUE if no replicates!
cor_data <- assay(vsd)

transcript_TPM_matrix <- assay(tr_dds, "abundance")
gene_TPM_matrix <- assay(dds, "abundance")

#### DESeq2 ####

# correlation

sample_cor <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")
sample_dist <- as.matrix(dist(t(cor_data)))

unlock(mtx)

###############################################################################



##################### READS/GENE COUNT COMPARISON #############################

cols <- colnames(gene_count_matrix)
n  <- length(cols)

cutoff <- 25
# count rows where BOTH col_i and col_j exceed cutoff
mat <- matrix(0L, nrow = n, ncol = n, dimnames = list(cols, cols))
for (i in 1:n) {
  for (j in 1:n) {
    mat[i, j] <- sum(gene_count_matrix[,i] >= cutoff & gene_count_matrix[,j] >= cutoff)
  }
}

detected <- sum(rowMax(gene_count_matrix) >= cutoff)
detected2 <- sum(rowMax(gene_count_matrix) >= 10)

mat <- mat/detected

saveRDS(mat, file="results/rds/read_count_comp.rds")

###############################################################################



##################### READS/GENE COUNTS #######################################

df <- data.frame("G0" = colSums(gene_count_matrix == 0),
                 "G1" = colSums(gene_count_matrix > 0 &
                                  gene_count_matrix <= 10),
                 "G10" = colSums(gene_count_matrix > 10 &
                                   gene_count_matrix < 25),
                 "G100" = colSums(gene_count_matrix >= 25))

df$Replicate <- rownames(df)
df$Sample <- sample_info[df$Replicate, ]$condition

dfm <- melt(df, id.vars = c("Replicate", "Sample"))
for (sample_name in unique(sample_info$condition)) {
  message(sample_name)
  if (sum(sample_info$condition == sample_name) > 1) {
    scols <- sample_info[sample_info$condition == sample_name, ]$names
    max_counts <- apply(gene_count_matrix[, scols], 1, max, na.rm = TRUE)

    dfm <- rbind(dfm, list(Replicate = sample_name,
                           Sample = sample_name,
                           variable = "G0",
                           value = sum(max_counts == 0)))

    dfm <- rbind(dfm, list(Replicate = sample_name,
                           Sample = sample_name,
                           variable = "G1",
                           value = sum(max_counts > 0 & max_counts <= 10)))

    dfm <- rbind(dfm, list(Replicate = sample_name,
                           Sample = sample_name,
                           variable = "G10",
                           value = sum(max_counts > 10 & max_counts < 25)))

    dfm <- rbind(dfm, list(Replicate = sample_name,
                           Sample = sample_name,
                           variable = "G100",
                           value = sum(max_counts >= 25)))
  }
}

dfm$variable <- as.character(dfm$variable)

dfm[dfm$variable == "G0", ]$variable <- "No reads"
dfm[dfm$variable == "G1", ]$variable <- "Low read count"
dfm[dfm$variable == "G10", ]$variable <- "Medium read count"
dfm[dfm$variable == "G100", ]$variable <- "Good read count"

dfm$variable <- factor(dfm$variable,
                       levels = c("No reads", "Low read count",
                                  "Medium read count", "Good read count"))

dfm$Replicate <- as.character(dfm$Replicate)
for (s in unique(dfm$Sample)) {
  dfm[dfm$Sample == s, ]$Replicate <-
    gsub(paste0(s, "_"), "", dfm[dfm$Sample == s, ]$Replicate)
  dfm[dfm$Sample == s, ]$Replicate <-
    gsub(paste0(s, "-"), "", dfm[dfm$Sample == s, ]$Replicate)
  dfm[dfm$Sample == s, ]$Replicate <-
    gsub(s, "S", dfm[dfm$Sample == s, ]$Replicate)
}


sorted_names <- c(mixedsort(unique(dfm$Replicate)),
                  mixedsort(unique(dfm$Sample)))
dfm$Replicate <- factor(dfm$Replicate, levels = sorted_names)

# gcp_colors <- c("#e0e0e0", "#faa39d", "#dbc160", "#00BA38")
# names(gcp_colors) <- levels(dfm$variable)
# 
# gene_coverage_plot <- ggplot(dfm, aes(x = Replicate,
#                                       y = value,
#                                       fill = variable)) +
#   geom_bar(stat = "identity") +
#   ylab("Number of Genes") + xlab("Samples") +
#   ggtitle("Number of Reads per Gene") +
#   guides(fill = guide_legend(title = "Number of assigned reads")) +
#   theme(text = element_text(size = 18)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
#   facet_wrap(~Sample, scale = "free_x") +
#   scale_fill_manual(values = gcp_colors)


saveRDS(list(dfm=dfm, detected=detected, detected2=detected2), file="results/rds/gene_coverage_plot.rds")

# print(gene_coverage_plot)

###############################################################################





##################### GENE BODY COVERAGE ######################################

deg_res <- data.frame(sample = character(),
                      quantile = character(),
                      pct = integer())

for (sample in samples) {
  message(sample)

  jfile <- paste0("results/salmon_quantiles/", sample, "/quant.sf")

  if (file.exists(jfile)) {
    dc <- read.table(jfile, header = TRUE, sep = "\t", row.names = 1)
    dc$quantile <- gsub(".*_(q[0-9]+)$", "\\1", row.names(dc))
    dc$ids <- gsub("(.*)_(q[0-9]+)$", "\\1", row.names(dc))

    dc <- dc[, c("ids", "quantile", "NumReads")]

    dcd <- as.data.frame(pivot_wider(dc,
                                     names_from = quantile,
                                     values_from = NumReads))
    dcd <- column_to_rownames(dcd, var = "ids")

    dcd <- dcd[rowSums(is.na(dcd)) == 0, ]
    dcd <- dcd[rowSums(dcd) > 0, ]

    dcd <- dcd / rowSums(dcd) * 100

    pct <- colSums(dcd) / sum(colSums(dcd)) * 100

    deg_res <- rbind(deg_res, data.frame(sample = sample,
                                         quantile = names(pct),
                                         pct = pct))
  }
}

deg_tmp <- pivot_wider(deg_res, names_from = quantile, values_from = pct)

deg_matrix <- as.matrix(column_to_rownames(deg_tmp, var = "sample"))
deg_matrix <- deg_matrix[, paste0("q", 1:10)]

skewness <- apply(deg_matrix, 1, function(data) {
  mx <- max(data[2:9])
  mi <- min(data[2:9])
  px <- which(data == mx)
  pi <- which(data == mi)

  (px - pi) * abs(mx-mi) / 10
})

skewclass <- rep("normal", length(skewness))
names(skewclass) <- names(skewness)

skewclass[names(skewness)[(skewness < -1)]] <- "light 3-prime"
skewclass[names(skewness)[(skewness < -3)]] <- "strong 3-prime"
skewclass[names(skewness)[(skewness >  1)]] <- "light 5-prime"
skewclass[names(skewness)[(skewness >  3)]] <- "strong 5-prime"

skewclass <- factor(skewclass, levels = c("strong 5-prime", "light 5-prime",
                                          "normal",
                                          "light 3-prime", "strong 3-prime"))

ccols <- list(degradation = c("light 3-prime" = "#ff8080",
                              "strong 3-prime" = "#800000",
                              "normal" = "#00ff00",
                              "light 5-prime" = "#8080ff",
                              "strong 5-prime" = "#000080"))

skewclass <- skewclass[rev(order(skewness))]

cclass <- data.frame(degradation = skewclass)

s <- read.table("reference/samples.tsv", header=T, sep="\t", row.names = 2)

s <- s[samples,]

if (length(colnames(s)) > 1) {
	faktors <- colnames(s)
	message(faktors)
	faktors <- faktors[!faktors %in% c("condition")]

	for(f in faktors) {
            flist = s[,f]
            names(flist) = rownames(s)

            flist <- factor(flist)
            cclass[[f]] <- flist[rownames(cclass)]

            fcols <- rainbow(length(unique(flist)))
            names(fcols) <- unique(flist)

            ccols[[f]] <- fcols
	}
}

deg_matrix <- deg_matrix[rev(order(skewness)),]

saveRDS(list(matrix = deg_matrix, colors = ccols, skewclass = skewclass, cclass = cclass), "results/rds/genebody_coverage.rds")

# ComplexHeatmap::pheatmap(as.matrix(deg_matrix), cluster_rows = FALSE, cluster_cols = FALSE,
#                          annotation_row = cclass, annotation_colors = ccols,
#                          treeheight_row = 0, fontsize = 12, border = FALSE,
#                          row_split = skewclass,
#                          heatmap_legend_param = list(title = "coverage (%)"),
#                          main = "RNA degradation (gene body coverage)")

###############################################################################





##################### CHECK FACTORS ###########################################
# https://www.bioconductor.org/packages/release/bioc/vignettes/DEGreport/inst/doc/DEGreport.html

geoMeanNZ <- function(x) {
  if (all(x == 0)) { 0 }
  else {
    exp(sum(log(x[x > 0])) / length(x[x > 0]))
  }
}
geoMeans <- apply(cor_data, 1, geoMeanNZ)
loggeomeans <- log(geoMeans)

df <- lapply(1:ncol(cor_data), function(smple) {
  cnts <- cor_data[,smple]
  r <- (log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]
  smple_name <- colnames(cor_data)[smple]
  data.frame(ratios = r, sample = smple_name, stringsAsFactors = FALSE)
}) %>% bind_rows()

df$replicate = df$sample
df$sample = sample_info[df$replicate,]$condition

# factor_plot <- ggplot(df, aes(ratios, col = sample, group = replicate)) + 
#     geom_vline(xintercept=0) + geom_density() + theme_bw()
# 
# print(factor_plot)
# 
# split_factor_plot <- ggplot(df, aes(ratios, col = sample, group = replicate)) + 
#     geom_vline(xintercept=0) + geom_density() + theme_bw() + facet_wrap(~sample)
# 
# print(split_factor_plot)

saveRDS(df, file="results/rds/split_factor_plot.rds")

###############################################################################






##################### SAMPLE CORRELATION ######################################

a <- data.frame(samples = sample_info[colnames(cor_data), 1])
row.names(a) <- colnames(cor_data)

s <- read.table("reference/samples.tsv", header=T, sep="\t", row.names = 2)
s <- s[samples, , drop = FALSE]

ccols <- NA

if (length(colnames(s)) > 1) {
    conds <- sort(unique(s$condition))
    condcols <- rainbow(length(conds))
    names(condcols) <- conds

    ccols = list(condition = condcols)

    faktors <- colnames(s)
    message(faktors)
    faktors <- faktors[!faktors %in% c("condition")]

    for(f in faktors) {
        flist = s[,f]
        names(flist) = rownames(s)
            
        flist <- factor(flist)
        a[[f]] <- flist[rownames(a)]
            
        fcols <- rainbow(length(unique(flist)))
        names(fcols) <- unique(flist)

        ccols[[f]] <- fcols
    }
}


# pheatmap::pheatmap(sample_cor,
#                    fontsize = 12,
#                    annotation_col = a,
#                    annotation_colors = ccols,
#                    main = "Sample Correlation Heatmap of VST counts")

saveRDS(list(sample_cor = sample_cor, a = a, ccols = ccols), "results/rds/correlation_hm.rds")

#######################################


### PCA ####

pca_res <- pca(cor_data, metadata = colData(se))
saveRDS(pca_res, "results/rds/pca_res.rds")

# biplot(pca_res, x = "PC2", y = "PC1",
#        colby = "condition", title = "PCA of VST counts")

##################



##### PCA OUTLIER #####



group_by_distance <- function(v, max_dist) {
    if(length(v) >= 2) {
        d <- dist(v)
        hc <- hclust(d, method = "single")
        cutree(hc, h = max_dist)
    } else {
        c(1)
    }
}


pc95 <- min(length(pca_res$variance), 6) #sum(cumsum(pca_res$variance) < 75)
outliers <- data.frame(sample = character(), pc = integer())
dns_values <- data.frame(x = double(), y = double(), pc = character())
peaks <- data.frame(peak = double(), pc = character)

for(pc in 1:pc95) {
    dst <- c()
    for(condition in unique(colData(se)$condition)) {
        if(sum(colData(se)$condition == condition) > 0) {
            points <- pca_res$rotated[colnames(se)[colData(se)$condition == condition],][pc]
            dst <- c(dst, as.vector(dist(points[,1])))
        }
    }

    peak <- 1000

    if(length(dst) > 1) {
        dns <- density(dst)
        peak <- dns$x[which.max(dns$y)]

        pcid = sprintf("PC%d [%.2f]", pc, peak)
        dns_values <- rbind(dns_values, data.frame(x = dns$x, y = dns$y, pc = pcid))
        peaks <- rbind(peaks, data.frame(peak = peak, pc = pcid))
    }
    
    message(peak)
    
    for(condition in unique(colData(se)$condition)) {
        if(sum(colData(se)$condition == condition) > 0) {
            points <- pca_res$rotated[colnames(se)[colData(se)$condition == condition],][pc]
        
            clstr <- group_by_distance(points[,1], 3*peak)

            names(clstr) <- rownames(points)
        
            outlier <- names(clstr)[clstr != names(which.max(table(clstr)))]

            if(length(outlier) > 0) {
                outliers <- rbind(outliers, data.frame(sample = outlier, pc = pc))
            }
        }
    }
}

if(nrow(outliers) > 0) {
    outliers$val <- "+"
    outlier_df <- as.data.frame(pivot_wider(outliers, id_cols = sample, names_from = pc, values_from = val, values_fill = "-"))

    outlier_df <- cbind(outlier_df, s[outlier_df$sample,])
    
    write.table(outlier_df, "results/rds/PCA_outlier_PC1-5.tsv", sep="\t", quote=F, row.names = F)

    # g <- ggplot(dns_values, aes(x=x, y=y)) +
    #     geom_vline(data = peaks, mapping = aes(xintercept = peak), col = "gray") +
    #     geom_line() +
    #     xlab("inner sample PCA distance") + 
    #     facet_wrap(~pc)
    # print(g)
    saveRDS(list(dns_values = dns_values, peaks = peaks), file="results/rds/pca_outlier.rds")
    
    outlier_ids <- unique(outliers$sample)
    all_samples <- rownames(s[s$condition %in% unique(s[outlier_ids,,drop=F]$condition),,drop=F])
    
    pca_data = as.data.frame(pca_res$rotated[all_samples,1:pc95])

    rownames_to_column(pca_data, "samples") %>%
        pivot_longer(cols=-samples, names_to = "PC", values_to = "rot") -> pca_data

    pca_data$condition <- s[pca_data$samples,,drop=FALSE]$condition

    outliers$pc = paste0("PC", outliers$pc)
    colnames(outliers) = c("samples", "PC")

    pca_data %>% semi_join(outliers, by = c("samples", "PC")) -> outlier_data

    variance <- sprintf("%.2f%%", pca_res$variance)[1:pc95]
    names(variance) = paste0("PC", 1:pc95)

    pca_data$PC = paste(pca_data$PC, " [", variance[pca_data$PC], "]", sep="")
    outlier_data$PC = paste(outlier_data$PC, " [", variance[outlier_data$PC], "]", sep="")
    
    # g <- ggplot(pca_data, aes(x = condition, y = rot)) +
    #     geom_point() + ggtitle("PCA outlier") + 
    #     geom_label(data=outlier_data, aes(label=samples), hjust = 0) +
    #     geom_point(data=outlier_data, col="red", size=2) +
    #     facet_wrap(~PC)
    # print(g)
    saveRDS(list(pca_data=pca_data, outlier_data=outlier_data), file="results/rds/pca_outlier_2.rds")
}


###############################################################################




if (length(colnames(s)) > 1) {
  pc <- pca_res$rotated[,1:5]
  pc <- rownames_to_column(pc, "sample")

  pcl <- pivot_longer(pc, cols=c(paste0("PC", 1:5)), names_to = "PC", values_to = "rot")

  pcl <- cbind(pcl, s[pcl$sample,])

  pcm <- pivot_longer(pcl, cols = colnames(s))

  compl <- list()
  for(batch in colnames(s)[!(colnames(s) %in% "condition")]) {
    batchs <- combn(unique(s[[batch]]), 2, simplify = F)
    compl <- append(compl, batchs)
  }

  # batch_pca_plot <- ggplot(pcm, aes(x = value, y = rot)) + geom_boxplot() + geom_jitter(width = 0.2) +
  #   stat_compare_means(comparison=compl, method="wilcox.test") +
  #   facet_grid(PC~name, scale = "free_x") + ggtitle("PC by batch") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # print(batch_pca_plot)
  saveRDS(list(pcm=pcm, compl=compl), file="results/rds/batch_pca_plot.rds")
}

###############################################################################





##################### FILTER SAMPLES ##########################################

## if (FALSE) {
## remove_by_reads <- results[results$status == "Assigned" &
##                              results$reference == "Gene" &
##                              results$reads <= min_assigned_read_count, ]$sample

## remove_by_correlation <- c()

## sample_cor_filtered <- sample_cor[! rownames(sample_cor) %in% remove_by_reads,
##                                   ! colnames(sample_cor) %in% remove_by_reads]

## for (sample_name in unique(sample_info$condition)) {
##   message(sample_name)

##   if (sum(sample_info$condition == sample_name) > 2) {
##     sr <- sample_info[sample_info$condition == sample_name, ]$names
##     sr <- sr[! sr %in% remove_by_reads]

##     if (length(sr) > 1) {
##       rep_num <- sum(sample_info$condition == sample_name)
##       low_cor <-
##         rowSums(sample_cor_filtered[sr, sr] <= min_sample_correlation)

##       if (sum(low_cor / rep_num > 0.67)) {
##         low_ids <- names(low_cor[low_cor / rep_num > 0.67])
##         remove_by_correlation <- c(remove_by_correlation, low_ids)
##       }
##     }
##   }
## }

## low_q_samples <- c(remove_by_reads, remove_by_correlation)

## message(low_q_samples)

## if (length(low_q_samples) > 0) {
##   lowq <- data.frame(sample = low_q_samples, reason = "low read number")
##   if (length(remove_by_correlation) > 0) {
##     remove_rows <- lowq$sample %in% remove_by_correlation
##     lowq[remove_rows, ]$reason <- "low in-sample correlation"
##   }

##   removed_plot <- ggplot(lowq, aes(y = sample, x = reason, col = I("red"))) +
##     geom_point(size = 4) +
##     theme(text = element_text(size = 18)) +
##     ylab("Replicates") + xlab(NULL) +
##     ggtitle("Low Quality Samples")

##   print(removed_plot)
## }
## }

## low_q_samples <- c()

###############################################################################





##################### FILTERED SAMPLES PLOT ###################################

## if (length(low_q_samples) > 0 && length(low_q_samples) < ncol(cor_data)) {
##   message("FILTERED:")
##   message(length(low_q_samples))
##   message(ncol(cor_data))
##   cor_data_filtered <- cor_data[, !colnames(cor_data) %in% low_q_samples]

##   sample_cor_filtered <- cor(cor_data_filtered,
##                              method = "pearson",
##                              use = "pairwise.complete.obs")

##   sample_dist_filtered <- as.matrix(dist(t(cor_data_filtered)))

##   pheatmap::pheatmap(sample_cor_filtered, fontsize = 10,
##                      annotation_col = a,
##                      main = "Filtered Sample Correlation Heatmap")

##   metadata <- colData(se)[!rownames(colData(se)) %in% low_q_samples, ]
##   pca_res <- pca(cor_data_filtered, metadata = metadata)

##   biplot(pca_res, x = "PC2", y = "PC1",
##          colby = "condition",
##          title = "PCA of VST counts (filtered samples)")
## }

###############################################################################





##################### (FILTERED) PCA PLOTS ####################################

# components <- getComponents(pca_res, seq_len(5))
# components <- components[!is.na(components)]
# 
# pair <- pairsplot(pca_res, colby = "condition", components = components)
# 
# scree <- screeplot(pca_res, axisLabSize = 18, titleLabSize = 22)
# 
# print(pair)
# print(scree)

### TODO add filtered PCA

###############################################################################





##################### WRITE PDF ###############################################

# dev.off()

###############################################################################


##################### HTML Report #############################################

folder <- getwd()

outfile_global = paste0(folder, "/", outfile)

message(outfile_global)

# file.create(outfile_global)

# q(save="no")

template <- paste0(snakemake@scriptdir, "/QCReport.Rmd")

rmarkdown::render(
  input = template,
  params = list(results = folder, min_assigned_read_count = min_assigned_read_count),
  output_file = outfile_global
)

###############################################################################


##################### WRITE DATA ##############################################

outfile <- gsub(".html", "", outfile)

deb <- 1
message(paste0("output ", deb)); deb <- deb + 1

xsheet1 <- pivot_wider(results[results$reference == "Gene", ],
                       id_cols = "sample",
                       names_from = "status",
                       values_from = "reads")
message(paste0("output ", deb)); deb <- deb + 1

xsheet1$trimmed <- (xsheet1$LowQuality +
                      xsheet1$Nreads +
                      xsheet1$TooShort +
                      xsheet1$TooLong)
message(paste0("output ", deb)); deb <- deb + 1


xsheet1 <- xsheet1 %>%
  dplyr::select(-one_of("LowQuality", "Nreads", "TooShort", "TooLong"))
message(paste0("output ", deb)); deb <- deb + 1

rrnar <- results[results$reference == "rRNA", ]$reads
message(paste0("output ", deb)); deb <- deb + 1

names(rrnar) <- results[results$reference == "rRNA", ]$sample
message(paste0("output ", deb)); deb <- deb + 1

xsheet1 <- cbind(xsheet1, rRNA = rrnar[xsheet1$sample])
message(paste0("output ", deb)); deb <- deb + 1

skewn <- skewness(t(deg_matrix))
message(paste0("output ", deb)); deb <- deb + 1

#xsheet1 <- cbind(xsheet1, skewness = skewn[xsheet1$sample])
message(paste0("output ", deb)); deb <- deb + 1

xsheet2 <- df

xsheet3 <- sample_cor

write.table(xsheet1, paste0(outfile, ".tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

#write.xlsx2(xsheet1, paste0(outfile, ".report.xlsx"),
#            sheetName = "sequencing results",
#            col.names = TRUE, row.names = FALSE, append = FALSE)
#gc()

#write.xlsx2(xsheet2, paste0(outfile, ".report.xlsx"),
#            sheetName = "reads per gene",
#            col.names = TRUE, row.names = FALSE, append = TRUE)
#gc()

#write.xlsx2(xsheet3, paste0(outfile, ".report.xlsx"),
#            sheetName = "sample correlation",
#            col.names = TRUE, row.names = TRUE, append = TRUE)
#gc()

write.table(gene_count_matrix, paste0(outfile, ".genes.raw.counts.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

write.table(gene_TPM_matrix, paste0(outfile, ".genes.TPM.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

write.table(TMM_counts, paste0(outfile, ".genes.TMM.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

write.table(gene_geTMM_counts, paste0(outfile, ".genes.geTMM.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1


write.table(transcript_count_matrix, paste0(outfile, ".transcripts.raw.counts.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

write.table(transcript_TPM_matrix, paste0(outfile, ".transcripts.TPM.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

write.table(tr_TMM_counts, paste0(outfile, ".transcripts.TMM.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

write.table(tr_geTMM_counts, paste0(outfile, ".transcripts.geTMM.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)
message(paste0("output ", deb)); deb <- deb + 1

###############################################################################

