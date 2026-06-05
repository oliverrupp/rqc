message(" == LIBRARY")

suppressMessages(library(getopt))
suppressMessages(library(gtools))
suppressMessages(library(tidyverse))
suppressMessages(library(filelock))
suppressMessages(library(tximeta))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(pracma))
suppressMessages(library(jsonlite))
suppressMessages(library(rmarkdown))
suppressMessages(library(DT))
suppressMessages(library(plotly))
suppressMessages(library(heatmaply))
suppressMessages(library(RColorBrewer))
suppressMessages(library(GGally))





##################### INIT ####################################################
init_project <- function() {
    message(" == INIT")

    spec <- matrix(c(
        "datafolder", "d", 1, "character",
        "samples", "s", 1, "character",
        "help", "h", 0, "logical"
    ), byrow = TRUE, ncol = 4)
    opt <- getopt(spec)

    samples_file <<- NULL

    if (exists("snakemake")) {
        samples_file <<- normalizePath(snakemake@input[["samples"]], mustWork = TRUE)
        setwd(dirname(snakemake@output[["report"]]))
        setwd(dirname(getwd()))
    } else {
        if (!is.null(opt$samples)) samples_file <<- normalizePath(opt$samples, mustWork = TRUE)
        if (!is.null(opt$datafolder)) setwd(opt$datafolder)
    }

    if(is.null(samples_file)) samples_file <<- normalizePath("reference/samples.tsv", mustWork = TRUE)

    project_name <<- basename(getwd())

    outdir <<- gsub(".tsv", "", basename(samples_file))
}

prj_message <- function(msg, lvl) {
    menu <- strrep(ifelse(lvl == 4, "-", "="), lvl)
    time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(paste0(time, " [", outdir, " - ", project_name, "] ", menu, " ", msg))
}
##################### INIT ####################################################





##################### SAMPLES #################################################
read_sample_info <- function() {
    prj_message("SAMPLES", 3)

    sample_info <<- read.table(samples_file, header = TRUE, sep = "\t") |> column_to_rownames("sample")
    
    sname <- mixedsort(unique(sample_info$condition))
    
    sample_info$files <<- paste0("results/salmon/", rownames(sample_info), "/quant.sf")
    sample_info$names <<- rownames(sample_info)
    
    replicate_counts <<- sample_info |> group_by(condition) |> summarize(number = n())
    
    prj_message(paste0("sample file: ", samples_file), 4)
    prj_message(paste0(nrow(sample_info), " samples"), 4)
    prj_message(paste0(length(unique(sample_info$condition)), " conditions"), 4)
}
##################### SAMPLES #################################################





##################### COLLECT COUNTS ##########################################
get_counts_from_salmon <- function() {
    prj_message("READ SALMON", 3)

    mtx <- lock("~/.tximeta.lock")

#### tximeta ####
    prj_message("tximeta", 4)
    suppressMessages({
        bcf_dir    <- file.path(getwd(), "results/index/BFC")
        index_dir  <- file.path(getwd(), "results/index/salmon")
        fasta_path <- file.path(getwd(), "reference/genome.fa")
        gtf_path   <- file.path(getwd(), "reference/annotation.gtf")
        
        if (!dir.exists(bcf_dir)) {
            dir.create(bcf_dir)
        }

        setTximetaBFC(bcf_dir)
        makeLinkedTxome(indexDir = index_dir,
                        source   = "denovo",
                        organism = getwd(),
                        release  = "1",
                        genome   = getwd(),
                        fasta    = fasta_path,
                        gtf      = gtf_path)
        
        t_se <- tximeta(sample_info)
        g_se <- summarizeToGene(t_se)
        
        t_count_matrix <<- assay(t_se, "counts")
        g_count_matrix <<- assay(g_se, "counts")

        check <- unlock(mtx)
    })
#### tximeta ####
    
    
##### edgeR #####
    prj_message("edgeR", 4)
    suppressMessages({
        ## TMM
        g_dge <- DGEList(counts=g_count_matrix, group=factor(sample_info[colnames(g_count_matrix),]$condition))
        g_dge <- g_dge[filterByExpr(g_dge), , keep.lib.sizes = FALSE]
        g_dge <- calcNormFactors(g_dge, method = "TMM")
        g_TMM <<- cpm(g_dge, normalized.lib.sizes = TRUE)
        
        t_dge <- DGEList(counts=t_count_matrix, group=factor(sample_info[colnames(t_count_matrix),]$condition))
        t_dge <- t_dge[filterByExpr(t_dge), , keep.lib.sizes = FALSE]
        t_dge <- calcNormFactors(t_dge, method = "TMM")
        t_TMM <<- cpm(t_dge, normalized.lib.sizes = TRUE)
        
        ## geTMM
        t_rpk <- t_count_matrix / assay(t_se, "length")
        g_rpk <- g_count_matrix / assay(g_se, "length")
        
        t_len <- assay(t_se, "length")
        t_len_m <- as.data.frame(t_len) |> rownames_to_column("id") |> 
            pivot_longer(cols=-id, names_to="sample", values_to = "length")
        
        t_norm_edger <- DGEList(counts=t_rpk,group=colData(t_se)$condition)
        g_norm_edger <- DGEList(counts=g_rpk,group=colData(g_se)$condition)
        
        t_norm_edger <- calcNormFactors(t_norm_edger, method = "TMM")
        g_norm_edger <- calcNormFactors(g_norm_edger, method = "TMM")
        
        t_geTMM <<- cpm(t_norm_edger)
        g_geTMM <<- cpm(g_norm_edger)
    })
####### edgeR ######
    
#### DESeq2 ####
    prj_message("DESeq2", 4)
    suppressMessages({
        colData(g_se)$condition <- factor(colData(g_se)$condition)
        colData(t_se)$condition <- factor(colData(t_se)$condition)
        
        g_dds <- DESeqDataSet(g_se, design = ~condition)
        t_dds <- DESeqDataSet(t_se, design = ~condition)
        
        g_dds <- estimateSizeFactors(g_dds)
        t_dds <- estimateSizeFactors(t_dds)

        
        if(sum(replicate_counts$number > 1) > 0) {
            g_dds <- estimateDispersions(g_dds)
            t_dds <- estimateDispersions(t_dds)

            g_dispersion <<- data.frame(
                mean = mcols(g_dds)$baseMean,
                geneEst = mcols(g_dds)$dispGeneEst,
                fit = mcols(g_dds)$dispFit,
                final = dispersions(g_dds)
            )
            
            t_dispersion <<- data.frame(
                mean = mcols(t_dds)$baseMean,
                geneEst = mcols(t_dds)$dispGeneEst,
                fit = mcols(t_dds)$dispFit,
                final = dispersions(t_dds)
            )

            rownames(g_dispersion) <<- rownames(g_dds)
            rownames(t_dispersion) <<- rownames(t_dds)
        } else {
            g_dispersion <<- data.frame(mean = 0, geneEst = 0, fit = 0, final = 0)
            t_dispersion <<- data.frame(mean = 0, geneEst = 0, fit = 0, final = 0)
        }

                                        # blind = TRUE if no replicates!
        t_vsd <- vst(t_dds, blind = sum(replicate_counts$number > 1) == 0) 
        t_vst <<- assay(t_vsd)

                                        # blind = TRUE if no replicates!
        g_vsd <- vst(g_dds, blind = sum(replicate_counts$number > 1) == 0) 
        g_vst <<- assay(g_vsd)

        
        ## TPM
        t_TPM <<- assay(t_dds, "abundance")
        g_TPM <<- assay(g_dds, "abundance")
    })
#### DESeq2 ####

####n correlation ####
    g_sample_correlation <- cor(g_vst, method = "pearson", use = "pairwise.complete.obs")
    g_sample_dist <- as.matrix(dist(t(g_vst)))
#### correlation ####

    for(cnts in c("count_matrix", "TPM", "TMM", "geTMM", "vst", "dispersion")) {
        for(ref in c("t", "g")) {
            v <- paste0(ref, "_", cnts)
            write_tsv(get(v), v)
        }
    }
}
##################### COLLECT COUNTS ##########################################





##################### TOP GENES ###############################################
compute_gene_detection <- function() {
    prj_message("TOP GENE READS", 3)
    read_fraction <- as.data.frame(g_count_matrix) |> mutate(across(everything(), ~ .x / sum(.x)))
    read_fraction_cumsum <- as.data.frame(apply(read_fraction, 2, function(x) cumsum(sort(x, decreasing = T))))
    read_fraction_cumsum$transcripts = 1:nrow(read_fraction_cumsum)

    write_tsv(read_fraction_cumsum, "gene_detection")
    
    top_x_transcripts <<- t(read_fraction_cumsum[c(10, 50, 100),] |> 
                            dplyr::select(-transcripts) * 100) |> as.data.frame() |>
        rename_with(~ paste0("top_", .x, "_transcripts"))
}
##################### TOP GENES ###############################################




##################### READS TO GENES ##########################################
reads_to_genes <- function() {
    prj_message("READS TO GENES", 3)

    fractions <- seq(0, 1, by = 0.01)

    gene_detection_curve <- function(counts) {
        sapply(fractions, function(f) { 
            sum(1 - (1 - f)^counts)
        })
    }

    reads_to_genes <<- apply(g_count_matrix, 2, gene_detection_curve)

    write_tsv(reads_to_genes, "reads_to_genes")
}
##################### READS TO GENES ##########################################





##################### LIB COMPLEXITY ##########################################
library_complexity <- function() {
    prj_message("LIB COMPLEXITY", 3)
    
    shannon_entropy <<- apply(t_count_matrix, 2, function(x) {
        p <- x / sum(x)
        p <- p[p > 0]
        
        -sum(p * log(p))
    })
}
##################### LIB COMPLEXITY ##########################################




##################### JACCARD TOP 1000 #########################################
jaccard_similarity <- function() {
    top_genes <- apply(g_count_matrix, 2, function(x) {
        names(sort(x, decreasing=TRUE)[1:1000])
    })

    samples <- colnames(top_genes)

    jaccard_mat <- matrix(
        NA,
        nrow = length(samples),
        ncol = length(samples),
        dimnames = list(samples, samples)
    )

    for(i in seq_along(samples)) {
        for(j in seq_along(samples)) {
            
            a <- top_genes[,i]
            b <- top_genes[,j]
    
            intersection <- length(intersect(a, b))
            union <- length(union(a,b))
    
            jaccard_mat[i,j] <- intersection / union
        }
    }

    jaccard_df <<- as.data.frame(jaccard_mat)

    write_tsv(jaccard_df, "jaccard")
}


###############################################################################




##################### SATURATION ##############################################
get_saturation <- function() {
    prj_message("SATURATION", 3)

    fractions <- seq(0.05, 1, by=0.05)

    results <- list()

    for(sample in colnames(t_count_matrix)) {
        x <- t_count_matrix[,sample]
        
        max_detected <- sum(x >= 10)
        
        for(frac in fractions) {
            subsampled <- rbinom(
                length(x),
                size = round(x),
                prob = frac
            )
            
            detected <- sum(subsampled >= 10)
            
            results[[length(results)+1]] <- data.frame(
                sample = sample,
                fraction = frac,
                detected = detected / max_detected
            )
        }
    }

    saturation_auc_df <<- do.call(rbind, results) |>
        group_by(sample) |> 
        arrange(fraction, .by_group = TRUE) |>
        summarise(
            saturation_auc = trapz(fraction, detected)
        ) |> column_to_rownames("sample")
}
##################### SATURATION ##############################################





##################### READ MAPPING ############################################
read_mapping <- function() {
    prj_message("READ MAPPING", 3)
    
    
    results <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(results) <- c("Duprate", "TotalReads", 
                           "LowQuality", "Nreads", "TooShort", "TooLong", 
                           "Unmapped", "NoFeature", "Assigned",
                           "rRNA")
    
    lib_type <- c()
    e_type <- c()
    
    for (sample in sample_info$names) {
        ## FastP
        jfile <- paste0("results/trimmed/", sample, ".json")
        jdat <- fromJSON(jfile)
        
        type <- jdat$summary$sequencing
        div <- 1
        e_type <- c(e_type, type)
        
        if (regexpr("paired end", type) >= 0) {
            div <- 2
        }

        res_v <- c(jdat$duplication$rate * 100,
                   jdat$summary$before_filtering$total_reads / div,
                   jdat$filtering_result$low_quality_reads / div,
                   jdat$filtering_result$too_many_N_reads / div,
                   jdat$filtering_result$too_short_reads / div,
                   jdat$filtering_result$too_long_reads / div)
        
        
        ## Salmon
        jfile <- paste0("results/salmon/", sample, "/aux_info/meta_info.json")
        jdat <- fromJSON(jfile)
        
        res_v <- c(res_v, 
                   jdat$num_processed - (jdat$num_decoy_fragments + jdat$num_mapped),
                   jdat$num_decoy_fragments, 
                   jdat$num_mapped)

        lib_type <- c(lib_type, jdat$library_types[1])
        
        ## rRNA
        jfile <- paste0("results/salmon_rrna/", sample, "/aux_info/meta_info.json")
        if (file.exists(jfile)) {
            jdat <- fromJSON(jfile)
            
            res_v <- c(res_v, jdat$num_mapped)
        } else { res_v <- c(res_v, -1) }
        
        results[sample,] <- res_v
    }

    lib_type <- data.frame(lib_type=lib_type, e_type=e_type)
    rownames(lib_type) = sample_info$names
    write_tsv(lib_type, "lib_type")

    read_mapping_df <<- results
}
##################### READ MAPPING ############################################





##################### 5'/3' BIAS ##############################################
get_coverage_bias <- function() {
    prj_message("5'/3' BIAS", 3)
    
    deg_res <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(deg_res) <- paste0("q", 1:10)
    
    lapply(sample_info$names, function(sample) {
        jfile <- file.path("results/salmon_quantiles", sample, "quant.sf")
        
        if (!file.exists(jfile)) { return(NULL) }
        
        dc <- read_tsv(jfile, show_col_types = FALSE) |>
            dplyr::select(Name, NumReads) |>
            mutate(ids      = str_remove(Name, "_q[0-9]+$"),
                   quantile = str_extract(Name, "q[0-9]+$"))
        
        dcd <- dc |> pivot_wider(id_cols = ids, names_from = quantile, values_from = NumReads) |> 
            column_to_rownames("ids")
        
        dcd[is.na(dcd)] <- 0
        dcd <- dcd[rowSums(dcd) >= 100,]
        
        mat <- as.matrix(dcd)
        rs <- rowSums(mat)
        mat <- mat / rs * 100
        
        pct <- colSums(mat)
        pct <- pct / sum(pct) * 100
        
        deg_res[sample, ] <<- pct[colnames(deg_res)]
    })

    write_tsv(deg_res, "gene_body_coverage")
    
    
    cov_bias_skewness <- apply(deg_res, 1, function(data) {
        mx <- max(data[2:9])
        mi <- min(data[2:9])
        px <- which(data == mx)
        pi <- which(data == mi)

        ((px - pi) / 7) * (abs(mx-mi))
    }) |> as.data.frame() |> rename_with(~ "skewness")

    write_tsv(cov_bias_skewness, "coverage_skewness")
}
##################### 5'/3' BIAS ##############################################





##################### PCA OUTLIER  ############################################
pca_outlier_global <- function() {
    prj_message("PCA", 3)

    pca <- prcomp(t(g_vst), scale. = FALSE, center = TRUE)

    threshold <- 0.9
    var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    n_pcs <- min(which(var_explained >= threshold)[1], ncol(pca$x))
    
    scores <- as.data.frame(pca$x[, 1:n_pcs])
    
    loo_outliers <- function(scores, labels, n_pcs = 2, threshold = 1.5) {
        scores_mat <- as.matrix(scores[, 1:n_pcs])
        groups     <- unique(labels)
        results    <- list()
        
        for (grp in groups) {
            idx <- which(labels == grp)
            X   <- scores_mat[idx, , drop = FALSE]
            n   <- nrow(X)
            
            if (n < 2) {
                results[[as.character(grp)]] <- data.frame(
                    sample         = rownames(scores)[idx],
                    n_pcs          = n_pcs,
                    label          = grp,
                    influence      = 0
                )
            } else {
                centroid  <- colMeans(X)
                
                influence <- sapply(seq_len(n), function(i) {
                    loo_centroid <- colMeans(X[-i, , drop = FALSE])
                    sqrt(sum((centroid - loo_centroid)^2))
                })
                
                results[[as.character(grp)]] <- data.frame(
                    sample         = rownames(scores)[idx],
                    n_pcs          = n_pcs,
                    label          = grp,
                    influence      = influence
                )
            }
        }
        
        all_results <- bind_rows(results)
        
        med     <- median(all_results$influence)
        mad_val <- median(abs(all_results$influence - med))
        cutoff  <- med + threshold * (mad_val + 1e-10)
        
        all_results |>
            mutate(
                cutoff     = cutoff,
                is_outlier = influence > cutoff
            )
    }
    

    result <- loo_outliers(scores, 
                           labels = sample_info[rownames(scores),]$condition, 
                           n_pcs = n_pcs, 
                           threshold = 3)

    result$var_explained = var_explained[n_pcs]
    
    n_pcs <- min(5, nrow(pca$x))
    scores  <- as.data.frame(pca$x[, 1:n_pcs])
    scores["__VAR",] <- (100*(pca$sdev^2 / sum(pca$sdev^2)))[1:n_pcs]
    
    write_tsv(scores, "pca")
    write_tsv(result, "pca_outlier")
}
##################### PCA OUTLIER ##############################################





##################### CONDITION PCA ############################################
condition_pca <- function() {
    prj_message("CONDITION PCA", 3)
    
    pca <- prcomp(t(g_vst), center = TRUE, scale. = FALSE)
    n_pcs <- min(nrow(pca$x), 5)
    
    pc_df <- as.data.frame(pca$x[, 1:n_pcs])
    
    pc_scaled <- pc_df
    
    factor_cols <- setdiff(colnames(sample_info), c("files", "names"))
    
    condition_qc <- do.call(rbind, lapply(factor_cols, function(fac) {
        
        do.call(rbind, lapply(unique(sample_info[[fac]]), function(cond) {
            
            idx <- sample_info[[fac]] == cond
            
            pcs <- pc_scaled[idx, , drop = FALSE]
            
            if (nrow(pcs) < 2) { return(data.frame(factor = fac, name = cond, qc_score = 0)) }
            
            centroid <- apply(pcs, 2, median)
            
            d <- apply(pcs, 1, function(x) { sqrt(sum((x - centroid)^2)) })
            
            data.frame(factor = fac, name = cond, qc_score = median(d))
        }))
    }))
    
    condition_pca_df <<- condition_qc
    
    write_tsv(condition_pca_df, "condition_df")
}
##################### CONDITION PCA ############################################





##################### CHECK FACTORS ###########################################
size_factor_qc <- function() {
    prj_message("SIZE FACTORS", 3)

                                        # https://www.bioconductor.org/packages/release/bioc/vignettes/DEGreport/inst/doc/DEGreport.html

    geoMeanNZ <- function(x) {
        if (all(x == 0)) { 0 }
        else {
            exp(sum(log(x[x > 0])) / length(x[x > 0]))
        }
    }
    geoMeans <- apply(g_vst, 1, geoMeanNZ)
    loggeomeans <- log(geoMeans)

    df <- lapply(1:ncol(g_vst), function(smple) {
        cnts <- g_vst[,smple]
        r <- (log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]
        smple_name <- colnames(g_vst)[smple]
        data.frame(ratios = r, sample = smple_name, stringsAsFactors = FALSE)
    }) %>% bind_rows()

    df$replicate = df$sample
    df$sample = sample_info[df$replicate,]$condition

    write_tsv(df, "size_factor_qc")  
}
##################### CHECK FACTORS ###########################################





##################### WRITE TSV RESULTS #######################################
write_tsv <- function(df, file) {
    dir <- paste0("tsv/", outdir)
    
    if( !dir.exists(dir) ) {
        dir.create(dir, recursive = TRUE)
    }

    prj_message(paste0(file, " [", nrow(df), "]"), 4)
    
    write.table(df, paste0(dir, "/", file, ".tsv"), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
}
##################### WRITE TSV RESULTS #######################################





##################### RENDER HTML #############################################
render_html <- function() {
    prj_message("RENDER HTML", 3)

    script_dir <- if (exists("snakemake")) {
                      snakemake@scriptdir
                  } else if (rstudioapi::isAvailable()) {
                      dirname(rstudioapi::getActiveDocumentContext()$path)
                  } else {
                      dirname(normalizePath(commandArgs(trailingOnly = FALSE) |>
                                            grep("--file=", x = _, value = TRUE) |>
                                            sub("--file=", "", x = _)))
                  }
    
    template <- paste0(script_dir, "/rqc.Rmd")
    
    if(!dir.exists(paste0(getwd(), "/report"))) dir.create(paste0(getwd(), "/report"), recursive = TRUE)

    html_report_file <- paste0(getwd(), "/report/", outdir, ".report.html")

    prjdir <- getwd()
    
    rmarkdown::render(
                   input = template,
                   params = list(project = prjdir, samples_file = samples_file),
                   output_file = html_report_file)

}
##################### RENDER HTML #############################################




##################### RUN ANALYSES ############################################
if(0) {
    setwd("~/Projects/trqc/art")
    setwd("~/Projects/trqc/CS")
    setwd("~/Projects/trqc/PP")
}

### init paths ####
init_project()
### init paths ####

### read results ####
read_sample_info()
get_counts_from_salmon()
### read results ####

### complexity ####
compute_gene_detection()
library_complexity()
reads_to_genes()
get_saturation()
read_mapping()
jaccard_similarity()

complexity_df <<- cbind(read_mapping_df, top_x_transcripts, shannon_entropy)

write_tsv(complexity_df, "complexity")
### complexity ####

### degradation ####
get_coverage_bias()
### degradation ####

### pca outlier ####
pca_outlier_global()
condition_pca()
### pca outlier ####

### size factor ####
size_factor_qc()
### size factor ####

### render HTML ####
render_html()
### render HTML ####

prj_message("DONE", 2)

