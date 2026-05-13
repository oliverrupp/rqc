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
suppressMessages(library(stringr))

#### setwd("/home/orupp/Projects/trqc/CS")

##################### INIT ####################################################
message(" ----- INIT")
### DEBUG: rm(list=ls())
### DEBUG: setwd("/vol/ranomics/rnaseq/QC/TT")

min_assigned_read_count <- 100
min_sample_correlation <- 0.9

spec <- matrix(c(
  "datafolder", "d", 1, "character",
  "reads", "r", 1, "integer",
  "correlation", "c", 1, "double",
  "help", "h", 0, "logical"
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

if( !dir.exists("results/rds") ) {
  dir.create("results/rds", recursive = TRUE)
}

# pdf(outfile, w = 24, h = 16)

theme_set(theme_bw())

###############################################################################





##################### SAMPLES #################################################
message(" ----- SAMPLES")

sample_info <- read.table("reference/samples.tsv",
                          header = TRUE, sep = "\t", row.names = 2)

sname <- mixedsort(unique(sample_info$condition))
sample_colors <- rainbow(length(sname))
names(sample_colors) <- sname

sample_info$files <-
  paste0("results/salmon/", rownames(sample_info), "/quant.sf")
sample_info$names <- rownames(sample_info)

samples <- row.names(sample_info)



# message(sample_info)

###############################################################################


##################### FALCO QC REPORT #########################################
message(" ----- FALCO")

# ── PARSERS ───────────────────────────────────────────────────────────────────
read_falco <- function(path) {
  con   <- gzcon(file(path, "rb"))
  lines <- readLines(con, warn = FALSE); close(con)
  sec_idx <- grep("^>>", lines)
  secs    <- list()
  for (i in seq_along(sec_idx)) {
    hdr <- lines[sec_idx[i]]
    if (hdr == ">>END_MODULE") next
    end <- if (i < length(sec_idx)) sec_idx[i+1]-1 else length(lines)
    body <- lines[(sec_idx[i]+1):end]
    body <- body[!grepl("^>>END_MODULE", body)]
    nm   <- str_remove(hdr, "^>>") |> str_remove("\\s+(pass|warn|fail)$") |> str_trim()
    st   <- str_extract(hdr, "(pass|warn|fail)$")
    secs[[nm]] <- list(status = st %||% "unknown", lines = body)
  }
  secs
}

parse_basic <- function(sec) {
  ls <- sec$lines[!grepl("^#", sec$lines)]
  m  <- str_split_fixed(ls, "\t", 2)
  setNames(m[,2], m[,1])
}

parse_table <- function(sec) {
  hi <- tail(grep("^#", sec$lines), 1); if (is.na(hi)) return(NULL)
  hd <- str_remove(sec$lines[hi], "^#") |> str_split("\t") |> unlist()
  dl <- sec$lines[(hi+1):length(sec$lines)]; dl <- dl[nzchar(dl)]
  if (!length(dl)) return(NULL)
  mat <- str_split_fixed(dl, "\t", length(hd))
  df  <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- make.names(hd); df
}


# ── LOAD DATA ─────────────────────────────────────────────────────────────────

falco_dir <- "results/falco"
falco_files <- list.files(falco_dir, "\\.gz$", recursive = TRUE, full.names = TRUE)
sample_names <- gsub("\\.f(ast|)q\\.gz$", "", basename(falco_files))
all_data <- setNames(lapply(falco_files, read_falco), sample_names)

saveRDS(all_data, file="results/rds/falco_qc.rds")

###############################################################################



##################### READ COUNTS #############################################
message(" ----- READ COUNTS")

results <- data.frame(sample = character(),
                      status = character(),
                      reference = character(),
                      reads = integer())

tr <- c()
duprate <- c()

for (sample in samples) {
  # message(sample)
  jfile <- paste0("results/trimmed/", sample, ".json")
  jdat <- fromJSON(file = jfile)

  duprate <- c(duprate, jdat$duplication$rate)
  
  type <- jdat$summary$sequencing
  div <- 1

  if (regexpr("paired end", type) >= 0) {
    div <- 2
  }
  total_reads <- jdat$summary$before_filtering$total_reads / div
  
  tr <- c(tr, total_reads)
  
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
saveRDS(results, "results/rds/reads_plot.rds")

###############################################################################


######################## ADD SUMMARY ##########################################
message(" ----- SUMMARY 1")

PASS_FAIL <- data.frame(sample = samples)

results |> filter(status == "Assigned") |>  
  mutate(ASSIGNED_READS = case_when(reads > 50e6 ~ "PASS", reads < 20e6 ~ "FAIL", .default  = "ATTENTION")) |>
  dplyr::select(c(sample, ASSIGNED_READS)) |> inner_join(PASS_FAIL, by="sample") -> PASS_FAIL

results |> filter(status == "Unmapped") |> arrange(match(sample, samples)) |> mutate(tr = tr) |> mutate(pct_assigned = (tr-reads)/tr*100) |>
  mutate(PCT_MAPPED = case_when(pct_assigned > 85 ~ "PASS", pct_assigned < 70 ~ "FAIL", .default  = "ATTENTION")) |>
  dplyr::select(c(sample, PCT_MAPPED)) |> inner_join(PASS_FAIL, by="sample") -> PASS_FAIL


results |> filter(status == "rRNA") |> arrange(match(sample, samples)) |> mutate(tr = tr) |> mutate(pct_rrna = -1*reads/tr*100) |>  
  mutate(PCT_RRNA = case_when(pct_rrna < 5 ~ "PASS", pct_rrna > 20 ~ "FAIL", .default  = "ATTENTION")) |>
  dplyr::select(c(sample, PCT_RRNA)) |> inner_join(PASS_FAIL, by="sample") -> PASS_FAIL


names(duprate) <- samples

as.data.frame(duprate) |> 
  mutate(DUPRATE = case_when(duprate < 0.4 ~ "PASS", duprate > 0.7 ~ "FAIL", .default  = "ATTENTION")) |>
  rownames_to_column("sample") |> dplyr::select(c(sample, DUPRATE)) |> 
  inner_join(PASS_FAIL, by="sample") -> PASS_FAIL
  



###############################################################################


# #message(sample_info)
# sink(file.path(getwd(), "filtered_samples.txt"))
# print(paste0("Low Read Count [<, ", min_assigned_read_count, "]:"))
# print(results[results$status == "Assigned" & results$reads < min_assigned_read_count,]$sample)
# sink()
# 
# samples <- results[results$status == "Assigned" & results$reads >= min_assigned_read_count,]$sample
# 
# sample_info <- sample_info[samples,]



##################### COLLECT COUNTS ##########################################
message(" ----- READ SALMON")
suppressMessages({
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

#### DESeq2 ####
colData(gse)$condition <- factor(colData(gse)$condition)
colData(se)$condition <- factor(colData(se)$condition)

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
})

###############################################################################


##################### Saturation ##############################################
message(" ----- SATURATION")

counts <- counts(dds)
fractions <- seq(0.05, 1, by=0.05)

results <- list()

for(sample in colnames(counts)) {
  
  x <- counts[,sample]
  
  for(frac in fractions) {
    
    # probabilistic downsampling
    subsampled <- rbinom(
      length(x),
      size = round(x),
      prob = frac
    )
    
    detected <- sum(subsampled >= 5)
    
    results[[length(results)+1]] <- data.frame(
      sample = sample,
      fraction = frac,
      detected = detected
    )
  }
}

df <- do.call(rbind, results)

saveRDS(df, "results/rds/saturation.rds")

###############################################################################





##################### Gene Complexity #########################################
message(" ----- COMPLEXITY")

counts <- as.data.frame(counts(dds))

gene_fraction <- counts %>% mutate(across(everything(), ~ .x / sum(.x)))

gene_fraction_cumsum <- as.data.frame(apply(gene_fraction, 2, function(x) cumsum(sort(x, decreasing = T))))

gene_fraction_cumsum$gene = 1:nrow(gene_fraction_cumsum)

gfdf <- gene_fraction_cumsum |> pivot_longer(!gene)

gene_fraction_cumsum[c(10, 50, 100),] |> dplyr::select(-gene) * 100 -> top_x_genes

as.data.frame(t(top_x_genes)) |> 
  mutate(TOP100 = case_when(`100` < 50 ~ "PASS", `100`> 70 ~ "FAIL", .default  = "ATTENTION")) |>
  rownames_to_column("sample") |> dplyr::select(c(sample, TOP100)) |> 
  inner_join(PASS_FAIL, by="sample") -> PASS_FAIL

top_x_genes <- t(top_x_genes)

top_x_genes <- as.data.frame(top_x_genes) |> mutate(across(where(is.numeric), ~ sprintf("%.2f %%", .x)))

shannon_df <- apply(counts, 2, function(x) {
  
  p <- x / sum(x)
  p <- p[p > 0]
  
  -sum(p * log(p))
})

shannon_df <- data.frame(shannon=shannon_df, samples = names(shannon_df))
saveRDS(list(gfdf=gfdf, topx=top_x_genes, shannon=shannon_df), "results/rds/complexity.rds")

###############################################################################





##################### DISPERSION ##############################################
message(" ----- DISPERSION")

suppressMessages({
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
})

# Extract values
df <- data.frame(
  mean = mcols(dds)$baseMean,
  geneEst = mcols(dds)$dispGeneEst,
  fit = mcols(dds)$dispFit,
  final = dispersions(dds)
)

# Remove invalid rows
df <- df[
  is.finite(df$mean) &
    is.finite(df$geneEst) &
    df$mean > 0 &
    df$geneEst > 0,
]

saveRDS(df, file="results/rds/dispersion.rds")

###############################################################################



##################### READS/GENE COUNT COMPARISON #############################
message(" ----- READ COUNT COMP")

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
message(" ----- GENE COUNTS")

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
  # message(sample_name)
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


saveRDS(list(dfm=dfm, detected=detected, detected2=detected2), file="results/rds/gene_coverage_plot.rds")

# print(gene_coverage_plot)

###############################################################################





##################### GENE BODY COVERAGE ######################################
message(" ----- GENE BODY COVERAGE")

deg_res <- data.frame(sample = character(),
                      quantile = character(),
                      pct = integer())

for (sample in samples) {
  # message(sample)

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

as.data.frame(skewness) |> 
  mutate(BIAS53 = case_when(abs(skewness) > 3 ~ "FAIL", abs(skewness) > 1 ~ "ATTENTION", .default  = "PASS")) |>
  rownames_to_column("sample") |> dplyr::select(c(sample, BIAS53)) |> 
  inner_join(PASS_FAIL, by="sample") -> PASS_FAIL



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
	# message(faktors)
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

###############################################################################





##################### CHECK FACTORS ###########################################
message(" ----- FACTORS")

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


saveRDS(df, file="results/rds/split_factor_plot.rds")

###############################################################################






##################### SAMPLE CORRELATION ######################################
message(" ----- SAMPLE CORRELATION")

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
    # message(faktors)
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


saveRDS(list(sample_cor = sample_cor, a = a, ccols = ccols), "results/rds/correlation_hm.rds")

#######################################


############### PCA ###################
message(" ----- PCA")

pca <- prcomp(t(cor_data), scale. = FALSE)
saveRDS(pca, "results/rds/pca_res.rds")

n_pcs   <- min(ncol(pca$x), 5)
scores  <- as.data.frame(pca$x[, 1:n_pcs])
scores$condition <- colData(vsd)$condition

outliers <- lapply(unique(scores$condition), function(g) {
  
  mat    <- as.matrix(scores[scores$condition == g, 1:n_pcs])
  n_pcs_use <- min(ncol(mat), nrow(mat) - 1, 5)  # conservative for small n
  mat    <- mat[, 1:n_pcs_use, drop = FALSE]
  
#  mcd    <- covMcd(mat)
#  d2     <- mahalanobis(mat, center = mcd$center, cov = mcd$cov)
  
  centroid <- colMeans(mat)
  d2        <- apply(mat, 1, \(x) sqrt(sum((x - centroid)^2))) # Euclidean, not squared

  cutoff   <- median(d2) + 3 * mad(d2)

  data.frame(
    sample    = rownames(mat),
    condition = g,
    d2        = d2,
    outlier   = d2 > cutoff
  )
}) |> do.call(rbind, args = _)

outliers |> 
  mutate(PCA_OUTLIER = case_when(outlier ~ "FAIL", .default  = "PASS")) |> 
  dplyr::select(c(sample, PCA_OUTLIER)) |> 
  inner_join(PASS_FAIL, by="sample") -> PASS_FAIL


##################



# ##### PCA OUTLIER #####
# 
# group_by_distance <- function(v, max_dist) {
#     if(length(v) >= 2) {
#         d <- dist(v)
#         hc <- hclust(d, method = "single")
#         cutree(hc, h = max_dist)
#     } else {
#         c(1)
#     }
# }
# 
# 
# pc95 <- min(length(pca_res$variance), 6) #sum(cumsum(pca_res$variance) < 75)
# outliers <- data.frame(sample = character(), pc = integer())
# dns_values <- data.frame(x = double(), y = double(), pc = character())
# peaks <- data.frame(peak = double(), pc = character)
# 
# for(pc in 1:pc95) {
#     dst <- c()
#     for(condition in unique(colData(se)$condition)) {
#         if(sum(colData(se)$condition == condition) > 0) {
#             points <- pca_res$rotated[colnames(se)[colData(se)$condition == condition],][pc]
#             dst <- c(dst, as.vector(dist(points[,1])))
#         }
#     }
# 
#     peak <- 1000
# 
#     if(length(dst) > 1) {
#         dns <- density(dst)
#         peak <- dns$x[which.max(dns$y)]
# 
#         pcid = sprintf("PC%d [%.2f]", pc, peak)
#         dns_values <- rbind(dns_values, data.frame(x = dns$x, y = dns$y, pc = pcid))
#         peaks <- rbind(peaks, data.frame(peak = peak, pc = pcid))
#     }
#     
#     message(peak)
#     
#     for(condition in unique(colData(se)$condition)) {
#         if(sum(colData(se)$condition == condition) > 0) {
#             points <- pca_res$rotated[colnames(se)[colData(se)$condition == condition],][pc]
#         
#             clstr <- group_by_distance(points[,1], 3*peak)
# 
#             names(clstr) <- rownames(points)
#         
#             outlier <- names(clstr)[clstr != names(which.max(table(clstr)))]
# 
#             if(length(outlier) > 0) {
#                 outliers <- rbind(outliers, data.frame(sample = outlier, pc = pc))
#             }
#         }
#     }
# }
# 
# if(nrow(outliers) > 0) {
#     outliers$val <- "+"
#     outlier_df <- as.data.frame(pivot_wider(outliers, id_cols = sample, names_from = pc, values_from = val, values_fill = "-"))
# 
#     outlier_df <- cbind(outlier_df, s[outlier_df$sample,])
#     
#     write.table(outlier_df, "results/rds/PCA_outlier_PC1-5.tsv", sep="\t", quote=F, row.names = F)
# 
#     # g <- ggplot(dns_values, aes(x=x, y=y)) +
#     #     geom_vline(data = peaks, mapping = aes(xintercept = peak), col = "gray") +
#     #     geom_line() +
#     #     xlab("inner sample PCA distance") + 
#     #     facet_wrap(~pc)
#     # print(g)
#     saveRDS(list(dns_values = dns_values, peaks = peaks), file="results/rds/pca_outlier.rds")
#     
#     outlier_ids <- unique(outliers$sample)
#     all_samples <- rownames(s[s$condition %in% unique(s[outlier_ids,,drop=F]$condition),,drop=F])
#     
#     pca_data = as.data.frame(pca_res$rotated[all_samples,1:pc95])
# 
#     rownames_to_column(pca_data, "samples") %>%
#         pivot_longer(cols=-samples, names_to = "PC", values_to = "rot") -> pca_data
# 
#     pca_data$condition <- s[pca_data$samples,,drop=FALSE]$condition
# 
#     outliers$pc = paste0("PC", outliers$pc)
#     colnames(outliers) = c("samples", "PC")
# 
#     pca_data %>% semi_join(outliers, by = c("samples", "PC")) -> outlier_data
# 
#     variance <- sprintf("%.2f%%", pca_res$variance)[1:pc95]
#     names(variance) = paste0("PC", 1:pc95)
# 
#     pca_data$PC = paste(pca_data$PC, " [", variance[pca_data$PC], "]", sep="")
#     outlier_data$PC = paste(outlier_data$PC, " [", variance[outlier_data$PC], "]", sep="")
#     
#     # g <- ggplot(pca_data, aes(x = condition, y = rot)) +
#     #     geom_point() + ggtitle("PCA outlier") + 
#     #     geom_label(data=outlier_data, aes(label=samples), hjust = 0) +
#     #     geom_point(data=outlier_data, col="red", size=2) +
#     #     facet_wrap(~PC)
#     # print(g)
#     saveRDS(list(pca_data=pca_data, outlier_data=outlier_data), file="results/rds/pca_outlier_2.rds")
# }


###############################################################################




# if (length(colnames(s)) > 1) {
#   pc <- pca_res$rotated[,1:5]
#   pc <- rownames_to_column(pc, "sample")
# 
#   pcl <- pivot_longer(pc, cols=c(paste0("PC", 1:5)), names_to = "PC", values_to = "rot")
# 
#   pcl <- cbind(pcl, s[pcl$sample,])
# 
#   pcm <- pivot_longer(pcl, cols = colnames(s))
# 
#   compl <- list()
#   for(batch in colnames(s)[!(colnames(s) %in% "condition")]) {
#     batchs <- combn(unique(s[[batch]]), 2, simplify = F)
#     compl <- append(compl, batchs)
#   }
# 
#   # batch_pca_plot <- ggplot(pcm, aes(x = value, y = rot)) + geom_boxplot() + geom_jitter(width = 0.2) +
#   #   stat_compare_means(comparison=compl, method="wilcox.test") +
#   #   facet_grid(PC~name, scale = "free_x") + ggtitle("PC by batch") +
#   #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   # print(batch_pca_plot)
#   saveRDS(list(pcm=pcm, compl=compl), file="results/rds/batch_pca_plot.rds")
# }

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
message(" ----- REPORT")

saveRDS(PASS_FAIL, "results/rds/summary.rds")

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
message(" ----- WRITE DATA")

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

