source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(RcisTarget)
  library(AUCell)
  library(data.table)
  library(dplyr)
  library(tibble)
})

must_getenv <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || identical(value, "")) {
    stop(sprintf("Missing environment variable: %s", name), call. = FALSE)
  }
  value
}

patch_import_rankings <- function() {
  ns <- asNamespace("RcisTarget")
  original <- get("importRankings", envir = ns)

  patched <- function(dbFile, indexCol = "features", ...) {
    args <- list(...)

    postprocess_result <- function(res) {
      if (is.null(res)) {
        return(res)
      }

      if (isS4(res) && methods::.hasSlot(res, "rankings")) {
        rankings <- res@rankings
        if (!"features" %in% names(rankings)) {
          if ("motifs" %in% names(rankings)) {
            names(rankings)[names(rankings) == "motifs"] <- "features"
          } else if (!is.null(rownames(rankings))) {
            rankings$features <- rownames(rankings)
          }
          res@rankings <- rankings
        }
        return(res)
      }

      if (!"features" %in% names(res)) {
        if ("motifs" %in% names(res)) {
          names(res)[names(res) == "motifs"] <- "features"
        } else if (!is.null(rownames(res))) {
          res$features <- rownames(res)
        }
      }
      res
    }

    tryCatch(
      {
        postprocess_result(do.call(original, c(list(dbFile = dbFile, indexCol = indexCol), args)))
      },
      error = function(e) {
        if (!grepl("index column|features|columns are missing", conditionMessage(e))) {
          stop(e)
        }

        retry_args <- args
        if ("columns" %in% names(retry_args) && !is.null(retry_args$columns)) {
          retry_args$columns <- gsub("features", "motifs", retry_args$columns, fixed = TRUE)
        }

        postprocess_result(do.call(original, c(list(dbFile = dbFile, indexCol = "motifs"), retry_args)))
      }
    )
  }

  environment(patched) <- environment(original)
  unlockBinding("importRankings", ns)
  assign("importRankings", patched, envir = ns)
  lockBinding("importRankings", ns)
  invisible(TRUE)
}

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)
patch_import_rankings()

input_rds <- file.path(cfg$checkpoint_dir, "03_after_annotation.rds")
adj_file <- file.path(cfg$scenic_output_dir, "adjacencies.tsv")
expr_csv <- file.path(cfg$scenic_input_dir, "expr_mat_human.csv")
ortholog_csv <- file.path(cfg$scenic_input_dir, "chicken_to_human_orthologs.csv")

scenic_db_500bp <- must_getenv("SCENIC_DB_500BP")
scenic_db_10kb <- must_getenv("SCENIC_DB_10KB")
scenic_motif_ann <- must_getenv("SCENIC_MOTIF_ANN")

if (!file.exists(input_rds)) {
  stop("Missing checkpoint: 03_after_annotation.rds", call. = FALSE)
}
if (!file.exists(adj_file)) {
  stop("Missing SCENIC network output: adjacencies.tsv", call. = FALSE)
}
if (!file.exists(expr_csv)) {
  stop("Missing SCENIC expression export: expr_mat_human.csv", call. = FALSE)
}
if (!file.exists(scenic_db_500bp) || !file.exists(scenic_db_10kb) || !file.exists(scenic_motif_ann)) {
  stop("Missing SCENIC database files", call. = FALSE)
}

module_top_n <- as.integer(Sys.getenv("SCENIC_MODULE_TOP_N", "50"))
module_min_genes <- as.integer(Sys.getenv("SCENIC_MODULE_MIN_GENES", "10"))
regulon_min_targets <- as.integer(Sys.getenv("SCENIC_REGULON_MIN_TARGETS", "5"))
nes_threshold <- as.numeric(Sys.getenv("SCENIC_NES_THRESHOLD", "3"))
auc_rank_fraction <- as.numeric(Sys.getenv("SCENIC_AUC_MAX_RANK_FRACTION", "0.05"))
n_cores <- max(1L, as.integer(Sys.getenv("SCENIC_THREADS", "4")))

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)

expr_cells_x_genes <- fread(expr_csv, data.table = FALSE)
expr_cells <- expr_cells_x_genes[[1]]
expr_cells_x_genes[[1]] <- NULL
expr_cells_x_genes <- as.matrix(expr_cells_x_genes)
rownames(expr_cells_x_genes) <- expr_cells
expr_genes_x_cells <- t(expr_cells_x_genes)
storage.mode(expr_genes_x_cells) <- "numeric"

common_cells <- intersect(colnames(obj), colnames(expr_genes_x_cells))
if (length(common_cells) == 0) {
  stop("No overlapping cells between Seurat object and SCENIC expression matrix", call. = FALSE)
}
expr_genes_x_cells <- expr_genes_x_cells[, common_cells, drop = FALSE]
obj <- subset(obj, cells = common_cells)

adj <- fread(adj_file)
if (!all(c("TF", "target", "importance") %in% colnames(adj))) {
  adj <- fread(adj_file, sep = "\t")
}
if (!all(c("TF", "target", "importance") %in% colnames(adj))) {
  stop("adjacencies.tsv must contain TF, target and importance columns", call. = FALSE)
}

adj <- adj[target %in% rownames(expr_genes_x_cells)]
setorderv(adj, cols = c("TF", "importance"), order = c(1L, -1L))

adj_split <- split(adj, by = "TF", keep.by = FALSE)
modules <- lapply(adj_split, function(tf_dt) {
  unique_targets <- unique(tf_dt$target)
  unique_targets[seq_len(min(length(unique_targets), module_top_n))]
})
modules <- modules[lengths(modules) >= module_min_genes]
names(modules) <- paste0(names(modules), "_top", module_top_n)

if (length(modules) == 0) {
  stop("No valid SCENIC modules were generated from adjacencies.tsv", call. = FALSE)
}

motif_annotations <- importAnnotations(scenic_motif_ann)
db_files <- c("500bp" = scenic_db_500bp, "10kb" = scenic_db_10kb)

motif_results <- lapply(names(db_files), function(db_name) {
  message("Running RcisTarget on database: ", db_name)
  rankings <- importRankings(db_files[[db_name]])
  enrichment <- cisTarget(modules, rankings, motif_annotations, nCores = n_cores)
  enrichment$db <- db_name
  as.data.frame(enrichment)
})
names(motif_results) <- names(db_files)

motif_df <- rbindlist(motif_results, fill = TRUE)
fwrite(motif_df, file.path(cfg$scenic_output_dir, "motif_enrichment.tsv"), sep = "\t")
saveRDS(motif_results, file.path(cfg$scenic_output_dir, "motif_enrichment.rds"))

sig_motifs <- motif_df %>% filter(NES > nes_threshold)
if (nrow(sig_motifs) == 0) {
  stop("No significant motifs passed the configured NES threshold", call. = FALSE)
}

split_targets <- function(x) {
  values <- unlist(strsplit(x, ";|,", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  values[nzchar(values)]
}

regulons_list <- list()
for (idx in seq_len(nrow(sig_motifs))) {
  row <- sig_motifs[idx, ]
  tf_name <- sub(paste0("_top", module_top_n, "$"), "", row$geneSet)
  enriched_genes <- split_targets(row$enrichedGenes)
  if (length(enriched_genes) == 0) {
    next
  }
  current_targets <- regulons_list[[tf_name]]
  if (is.null(current_targets)) {
    current_targets <- character(0)
  }
  regulons_list[[tf_name]] <- unique(c(current_targets, enriched_genes))
}

regulons_list <- regulons_list[lengths(regulons_list) >= regulon_min_targets]
regulons_list <- lapply(regulons_list, intersect, rownames(expr_genes_x_cells))
regulons_list <- regulons_list[lengths(regulons_list) >= regulon_min_targets]

if (length(regulons_list) == 0) {
  stop("No regulons remained after motif filtering and expression matching", call. = FALSE)
}

regulon_summary <- tibble(
  TF = names(regulons_list),
  TargetCount = lengths(regulons_list)
) %>%
  arrange(desc(TargetCount), TF)
write.csv(regulon_summary, file.path(cfg$scenic_output_dir, "regulons.csv"), row.names = FALSE)

regulon_targets_long <- bind_rows(lapply(names(regulons_list), function(tf_name) {
  tibble(TF = tf_name, Target = regulons_list[[tf_name]])
}))
write.csv(regulon_targets_long, file.path(cfg$scenic_output_dir, "regulon_targets_long.csv"), row.names = FALSE)

gmt_file <- file.path(cfg$scenic_output_dir, "regulons.gmt")
gmt_conn <- file(gmt_file, open = "w")
on.exit(close(gmt_conn), add = TRUE)
for (tf_name in names(regulons_list)) {
  writeLines(paste(c(tf_name, "NA", regulons_list[[tf_name]]), collapse = "\t"), gmt_conn)
}

cells_rankings <- AUCell_buildRankings(expr_genes_x_cells, plotStats = FALSE, nCores = n_cores)
auc_max_rank <- max(1L, as.integer(ceiling(nrow(expr_genes_x_cells) * auc_rank_fraction)))
cells_auc <- AUCell_calcAUC(regulons_list, cells_rankings, aucMaxRank = auc_max_rank)
auc_features_x_cells <- getAUC(cells_auc)
auc_cells_x_regulons <- t(auc_features_x_cells)
write.csv(as.data.frame(auc_cells_x_regulons), file.path(cfg$scenic_output_dir, "auc_matrix.csv"), row.names = TRUE)

saveRDS(
  list(
    regulons = regulons_list,
    motif_enrichment = motif_results,
    auc_matrix = auc_features_x_cells,
    ortholog_map = if (file.exists(ortholog_csv)) read.csv(ortholog_csv, stringsAsFactors = FALSE) else NULL
  ),
  file.path(cfg$scenic_output_dir, "scenic_core_outputs.rds")
)

message(sprintf("SCENIC regulons generated: %d", length(regulons_list)))
