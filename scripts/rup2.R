suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
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
suppressMessages(library(xlsx))
suppressMessages(library(moments))
suppressMessages(library(PCAtools))



#### setwd("/vol/ranomics/rnaseq/QC/AC/")

##################### INIT ####################################################

min_assigned_read_count <- 10e6
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

outfile <- "report.pdf"

if (exists("snakemake")) {
  setwd(dirname(snakemake@output[["report"]]))
  min_assigned_read_count <- snakemake@params[["reads"]]
  outfile <- basename(snakemake@output[["report"]])
}

pdf(outfile, w = 18, h = 12)

theme_set(theme_bw())

###############################################################################





##################### SAMPLES #################################################

sample_info <- read.table("reference/samples.tsv",
                          header = FALSE, sep = "\t", row.names = 2)
colnames(sample_info) <- "condition"

sname <- mixedsort(unique(sample_info$condition))
sample_colors <- rainbow(length(sname))
names(sample_colors) <- sname

sample_info$files <-
  paste0("results/salmon/", rownames(sample_info), "/quant.sf")
sample_info$names <- rownames(sample_info)

samples <- row.names(sample_info)

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

reads_plot <- ggplot(results, aes(x = sample, y = reads, fill = status)) +
  geom_hline(yintercept=c(30e6, 50e6)) +
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 18)) + ggtitle("Read count analysis") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0))

print(reads_plot)

###############################################################################





##################### COLLECT COUNTS ##########################################

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

gene_count_matrix <- assay(gse, "counts")
dds <- DESeqDataSet(gse, design = ~condition)
vsd <- vst(dds, blind = FALSE)
cor_data <- assay(vsd)

sample_cor <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")
sample_dist <- as.matrix(dist(t(cor_data)))

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

gcp_colors <- c("#e0e0e0", "#faa39d", "#dbc160", "#00BA38")
names(gcp_colors) <- levels(dfm$variable)

gene_coverage_plot <- ggplot(dfm, aes(x = Replicate,
                                      y = value,
                                      fill = variable)) +
  geom_bar(stat = "identity") +
  ylab("Number of Genes") + xlab("Samples") +
  ggtitle("Number of Reads per Gene") +
  guides(fill = guide_legend(title = "Number of assigned reads")) +
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  facet_wrap(~Sample, scale = "free_x") +
  scale_fill_manual(values = gcp_colors)

print(gene_coverage_plot)

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

deg_matrix <- deg_matrix[rev(order(skewness)),]

ComplexHeatmap::pheatmap(as.matrix(deg_matrix), cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = cclass, annotation_colors = ccols, 
                         treeheight_row = 0, fontsize = 12, border = FALSE,
                         row_split = skewclass, 
                         heatmap_legend_param = list(title = "coverage (%)"),
                         main = "RNA degradation (gene body coverage)")

###############################################################################





##################### SAMPLE CORRELATION ######################################

a <- data.frame(samples = sample_info[colnames(cor_data), 1])
row.names(a) <- colnames(cor_data)

pheatmap::pheatmap(sample_cor,
                   fontsize = 12,
                   annotation_col = a,
                   main = "Sample Correlation Heatmap of VST counts")

pca_res <- pca(cor_data, metadata = colData(se))

biplot(pca_res, x = "PC2", y = "PC1",
       colby = "condition", title = "PCA of VST counts")

###############################################################################





##################### FILTER SAMPLES ##########################################

if (FALSE) {
remove_by_reads <- results[results$status == "Assigned" &
                             results$reference == "Gene" &
                             results$reads <= min_assigned_read_count, ]$sample

remove_by_correlation <- c()

sample_cor_filtered <- sample_cor[! rownames(sample_cor) %in% remove_by_reads,
                                  ! colnames(sample_cor) %in% remove_by_reads]

for (sample_name in unique(sample_info$condition)) {
  message(sample_name)

  if (sum(sample_info$condition == sample_name) > 2) {
    sr <- sample_info[sample_info$condition == sample_name, ]$names
    sr <- sr[! sr %in% remove_by_reads]

    if (length(sr) > 1) {
      rep_num <- sum(sample_info$condition == sample_name)
      low_cor <-
        rowSums(sample_cor_filtered[sr, sr] <= min_sample_correlation)

      if (sum(low_cor / rep_num > 0.67)) {
        low_ids <- names(low_cor[low_cor / rep_num > 0.67])
        remove_by_correlation <- c(remove_by_correlation, low_ids)
      }
    }
  }
}

low_q_samples <- c(remove_by_reads, remove_by_correlation)

message(low_q_samples)

if (length(low_q_samples) > 0) {
  lowq <- data.frame(sample = low_q_samples, reason = "low read number")
  if (length(remove_by_correlation) > 0) {
    remove_rows <- lowq$sample %in% remove_by_correlation
    lowq[remove_rows, ]$reason <- "low in-sample correlation"
  }

  removed_plot <- ggplot(lowq, aes(y = sample, x = reason, col = I("red"))) +
    geom_point(size = 4) +
    theme(text = element_text(size = 18)) +
    ylab("Replicates") + xlab(NULL) +
    ggtitle("Low Quality Samples")

  print(removed_plot)
}
}

low_q_samples <- c()

###############################################################################





##################### FILTERED SAMPLES PLOT ###################################

if (length(low_q_samples) > 0 && length(low_q_samples) < ncol(cor_data)) {
  message("FILTERED:")
  message(length(low_q_samples))
  message(ncol(cor_data))
  cor_data_filtered <- cor_data[, !colnames(cor_data) %in% low_q_samples]

  sample_cor_filtered <- cor(cor_data_filtered,
                             method = "pearson",
                             use = "pairwise.complete.obs")

  sample_dist_filtered <- as.matrix(dist(t(cor_data_filtered)))

  pheatmap::pheatmap(sample_cor_filtered, fontsize = 10,
                     annotation_col = a,
                     main = "Filtered Sample Correlation Heatmap")

  metadata <- colData(se)[!rownames(colData(se)) %in% low_q_samples, ]
  pca_res <- pca(cor_data_filtered, metadata = metadata)

  biplot(pca_res, x = "PC2", y = "PC1",
         colby = "condition",
         title = "PCA of VST counts (filtered samples)")
}

###############################################################################





##################### (FILTERED) PCA PLOTS ####################################

pair <- pairsplot(pca_res, colby = "condition")
scree <- screeplot(pca_res, axisLabSize = 18, titleLabSize = 22)

print(pair)
print(scree)

###############################################################################





##################### WRITE PDF ###############################################

dev.off()

###############################################################################





##################### WRITE DATA ##############################################

xsheet1 <- pivot_wider(results[results$reference == "Gene", ],
                       id_cols = "sample",
                       names_from = "status",
                       values_from = "reads")

xsheet1$trimmed <- (xsheet1$LowQuality +
                      xsheet1$Nreads +
                      xsheet1$TooShort +
                      xsheet1$TooLong)


xsheet1 <- xsheet1 %>%
  dplyr::select(-one_of("LowQuality", "Nreads", "TooShort", "TooLong"))

rrnar <- results[results$reference == "rRNA", ]$reads

names(rrnar) <- results[results$reference == "rRNA", ]$sample

xsheet1 <- cbind(xsheet1, rRNA = rrnar[xsheet1$sample])

skewn <- skewness(t(deg_matrix))

xsheet1 <- cbind(xsheet1, skewness = skewn[xsheet1$sample])

xsheet2 <- df

xsheet3 <- sample_cor

xsheet4 <- assay(dds, "counts")

xsheet5 <- assay(dds, "abundance")

write.table(xsheet1, paste0(outfile, ".report.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(xsheet1, paste0(outfile, ".report.xlsx"),
            sheetName = "sequencing results",
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(xsheet2, paste0(outfile, ".report.xlsx"),
            sheetName = "reads per gene",
            col.names = TRUE, row.names = FALSE, append = TRUE)

write.xlsx2(xsheet3, paste0(outfile, ".report.xlsx"),
            sheetName = "sample correlation",
            col.names = TRUE, row.names = TRUE, append = TRUE)

write.table(xsheet4, paste0(outfile, ".raw.counts.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = TRUE)

write.table(xsheet5, paste0(outfile, ".TMP.normalized.tsv"),
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE, append = FALSE)

###############################################################################
