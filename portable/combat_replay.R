#!/usr/bin/env Rscript
# Paper 8 / Stage 6 / frozen full-matrix ComBat replay.
# Inputs read-only. Metrics and environment files written only to supplied Downloads directory.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop("STAGE6_EXPECTED_FOUR_ARGUMENTS")
pre_file <- args[1]; post_file <- args[2]; meta_file <- args[3]; out_dir <- args[4]
if (!requireNamespace("sva", quietly=TRUE)) stop("HOLD_SVA_MISSING")
if (as.character(utils::packageVersion("sva")) != "3.54.0") stop("HOLD_SVA_VERSION_NOT_LOCKED_3_54_0")
message("[R_START] P8 Stage6 numeric replay: full 44341 x 216, no models")
message("[R_PACKAGE] sva ", as.character(utils::packageVersion("sva")))

# Deliberately do not use row names: HGNC symbols can occur more than once.
meta <- utils::read.delim(meta_file, sep="\t", header=TRUE,
                          check.names=FALSE, stringsAsFactors=FALSE,
                          quote="", comment.char="")
if (nrow(meta) != 216L || !identical(names(meta), c("sample_name","condition","phase","batch"))) {
    stop("HOLD_METADATA_SCHEMA")
}
if (anyNA(meta) || anyDuplicated(meta$sample_name) != 0L) stop("HOLD_METADATA_MISSING_OR_DUPLICATE")
pre <- utils::read.delim(pre_file, sep="\t", header=TRUE,
                         check.names=FALSE, stringsAsFactors=FALSE,
                         quote="", comment.char="")
if (!identical(dim(pre), c(44341L,217L))) stop("HOLD_PRECOMBAT_SHAPE")
if (names(pre)[1] != "hgnc_symbol") stop("HOLD_PRECOMBAT_SYMBOL_COLUMN")
if (!identical(as.character(names(pre)[-1]), as.character(meta$sample_name))) stop("HOLD_PRECOMBAT_GSM_ORDER")
if (anyNA(pre[[1]])) stop("HOLD_PRECOMBAT_GENE_IDS_MISSING")
# Avoid character coersion across the entire 216-column data frame.
x <- as.matrix(pre[,-1,drop=FALSE])
if (!is.numeric(x) || any(!is.finite(x))) stop("HOLD_PRECOMBAT_NONNUMERIC_OR_NONFINITE")
colnames(x) <- as.character(meta$sample_name)
gene_ids <- as.character(pre[[1]])
rm(pre); invisible(gc(verbose=FALSE))

# Match the previously locked original full-cohort design.
# Both training and test samples enter ComBat; TRUE test conditions are supplied.
# This is not leakage-free prospective preprocessing.
if (!all(as.character(meta$batch) %in% c("1","2","3","4"))) {
    message("[INFO] batch levels ", paste(unique(as.character(meta$batch)),collapse=","))
}
if (!("1" %in% unique(as.character(meta$batch)))) stop("HOLD_REFERENCE_BATCH_1_MISSING")
mod <- stats::model.matrix(~condition, data=meta)
if (nrow(mod) != 216L || any(!is.finite(mod))) stop("HOLD_DESIGN_INVALID")
if (qr(mod)$rank != ncol(mod)) stop("HOLD_DESIGN_RANK_DEFICIENT")
message("[R_DESIGN] ",paste(colnames(mod), collapse=","),
        " | batch_counts=",paste(as.character(table(meta$batch)),collapse=","))
message("[R_COMBAT_START] all genes; outcome-informed; ref.batch=1")
start <- proc.time()[[3]]
corrected <- sva::ComBat(dat=x, batch=meta$batch, mod=mod,
                         ref.batch=1, par.prior=TRUE, prior.plots=FALSE)
spent <- proc.time()[[3]] - start
if (!identical(dim(corrected), c(44341L,216L)) || any(!is.finite(corrected))) {
    stop("HOLD_REPLAY_NONFINITE_OR_SHAPE")
}
rm(x); invisible(gc(verbose=FALSE))
message("[R_COMBAT_DONE] elapsed_sec=",format(spent,digits=8))

# The archived numeric values are compared cell by cell. Text-file SHA256
# equality is NOT required because R numeric serialization can differ.
reference <- utils::read.delim(post_file, sep="\t", header=TRUE,
                               check.names=FALSE, stringsAsFactors=FALSE,
                               quote="", comment.char="")
if (!identical(dim(reference),c(44341L,217L)) || names(reference)[1]!="symbol") {
    stop("HOLD_ARCHIVED_POSTCOMBAT_SCHEMA")
}
if (!identical(names(reference)[-1], as.character(meta$sample_name))) {
    stop("HOLD_POSTCOMBAT_GSM_ORDER")
}
if (!identical(as.character(reference[[1]]),gene_ids)) stop("HOLD_GENE_ROW_ORDER")

max_abs <- 0.0; sum_abs <- 0.0; sum_sq <- 0.0
above_1e6 <- 0; above_1e8 <- 0; above_1e10 <- 0
n_numeric <- as.double(nrow(corrected)) * ncol(corrected)
worst_gene <- NA_character_; worst_gsm <- NA_character_; worst_gene_row <- NA_integer_
for (j in seq_len(ncol(corrected))) {
    archived <- reference[[j+1L]]
    if (!is.numeric(archived) || any(!is.finite(archived))) stop("HOLD_REFERENCE_NONNUMERIC_OR_NONFINITE")
    err <- abs(corrected[,j] - archived)
    local_max <- max(err)
    if (local_max > max_abs) {
        idx <- which.max(err)
        max_abs <- local_max
        worst_gene <- gene_ids[[idx]]
        worst_gsm <- colnames(corrected)[[j]]
        worst_gene_row <- idx
    }
    sum_abs <- sum_abs + sum(err)
    sum_sq <- sum_sq + sum(err*err)
    above_1e6 <- above_1e6 + sum(err > 1e-6)
    above_1e8 <- above_1e8 + sum(err > 1e-8)
    above_1e10 <- above_1e10 + sum(err > 1e-10)
    if ((j %% 54L)==0L) message("[R_COMPARE] ",j,"/216 samples")
}
rm(reference,corrected); invisible(gc(verbose=FALSE))
metrics <- data.frame(
    key=c("status","n_genes","n_samples","n_values","sva_version",
          "R_version","design_columns","reference_batch","computation_seconds",
          "max_abs_error","mean_abs_error","rmse", "n_abs_gt_1e-6", "n_abs_gt_1e-8",
          "n_abs_gt_1e-10","worst_gene","worst_gsm","worst_gene_row"),
    value=as.character(c(if (above_1e6 == 0L) "PASS_NUMERIC_COMBAT_REPLAY_1E6" else "HOLD_NUMERIC_DIFFERENCE",
          "44341","216",format(n_numeric,scientific=FALSE),
          as.character(utils::packageVersion("sva")),R.version.string,
          paste(colnames(mod),collapse=";"),"1",format(spent,digits=16),
          format(max_abs,digits=16,scientific=TRUE),format(sum_abs/n_numeric,digits=16,scientific=TRUE),
          format(sqrt(sum_sq/n_numeric),digits=16,scientific=TRUE),
          format(above_1e6,scientific=FALSE),format(above_1e8,scientific=FALSE),
          format(above_1e10,scientific=FALSE),worst_gene,worst_gsm,worst_gene_row)),
    stringsAsFactors=FALSE
)
utils::write.table(metrics, file.path(out_dir,"numeric_metrics.tsv"),
                   sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
message("[R_RESULT] ",metrics$value[[1]]," | max_abs=",format(max_abs,digits=16),
        " | mismatch_1e-6=",above_1e6)
