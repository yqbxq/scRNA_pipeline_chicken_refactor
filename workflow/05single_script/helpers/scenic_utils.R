patch_import_rankings_08 <- function() {
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

read_scenic_expr_csv_08 <- function(expr_csv) {
  expr_cells_x_genes <- data.table::fread(expr_csv, data.table = FALSE)
  expr_cells <- expr_cells_x_genes[[1]]
  expr_cells_x_genes[[1]] <- NULL
  expr_cells_x_genes <- as.matrix(expr_cells_x_genes)
  rownames(expr_cells_x_genes) <- expr_cells
  expr_genes_x_cells <- t(expr_cells_x_genes)
  storage.mode(expr_genes_x_cells) <- "numeric"
  expr_genes_x_cells
}

build_scenic_modules_from_adj_08 <- function(adj, expr_genes, module_top_n, module_min_genes) {
  if (!all(c("TF", "target", "importance") %in% colnames(adj))) {
    stop("adjacencies.tsv must contain TF, target and importance columns", call. = FALSE)
  }
  adj <- adj[target %in% expr_genes]
  data.table::setorderv(adj, cols = c("TF", "importance"), order = c(1L, -1L))

  adj_split <- split(adj, by = "TF", keep.by = FALSE)
  modules <- lapply(adj_split, function(tf_dt) {
    unique_targets <- unique(tf_dt$target)
    unique_targets[seq_len(min(length(unique_targets), module_top_n))]
  })
  modules <- modules[lengths(modules) >= module_min_genes]
  names(modules) <- paste0(names(modules), "_top", module_top_n)
  modules
}

split_scenic_targets_08 <- function(x) {
  values <- unlist(strsplit(x, ";|,", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  values[nzchar(values)]
}

regulons_from_motif_df_08 <- function(motif_df, expr_genes, module_top_n, regulon_min_targets, nes_threshold) {
  sig_motifs <- motif_df %>% dplyr::filter(NES > nes_threshold)
  if (nrow(sig_motifs) == 0) {
    stop("No significant motifs passed the configured NES threshold", call. = FALSE)
  }

  regulons_list <- list()
  for (idx in seq_len(nrow(sig_motifs))) {
    row <- sig_motifs[idx, ]
    tf_name <- sub(paste0("_top", module_top_n, "$"), "", row$geneSet)
    enriched_genes <- split_scenic_targets_08(row$enrichedGenes)
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
  regulons_list <- lapply(regulons_list, intersect, expr_genes)
  regulons_list <- regulons_list[lengths(regulons_list) >= regulon_min_targets]
  regulons_list
}

write_regulons_gmt_08 <- function(regulons_list, gmt_file) {
  ensure_dir(dirname(gmt_file))
  gmt_conn <- file(gmt_file, open = "w")
  on.exit(close(gmt_conn), add = TRUE)
  for (tf_name in names(regulons_list)) {
    writeLines(paste(c(tf_name, "NA", regulons_list[[tf_name]]), collapse = "\t"), gmt_conn)
  }
}

named_discrete_palette_08 <- function(levels, base_colors = NULL) {
  levels <- unique(as.character(levels))
  levels <- levels[nzchar(levels) & !is.na(levels)]
  if (length(levels) == 0) {
    return(setNames(character(0), character(0)))
  }

  if (is.null(base_colors) || length(base_colors) == 0) {
    base_colors <- c(
      "#79BB7D", "#C4A9D6", "#FDBE6F", "#2F6CB3", "#E30073",
      "#A54F37", "#69706A", "#0E9B84", "#D95F02", "#6B63B7"
    )
  }
  if (length(levels) > length(base_colors)) {
    base_colors <- c(base_colors, scales::hue_pal()(length(levels) - length(base_colors)))
  }
  setNames(base_colors[seq_along(levels)], levels)
}

calculate_csi_long_08 <- function(auc_features_x_cells, vendor_path) {
  if (requireNamespace("scFunctions", quietly = TRUE)) {
    auc_rankings <- AUCell::AUCell_buildRankings(auc_features_x_cells, plotStats = FALSE, nCores = 1)
    auc_rankings@assays@data@listData$AUC <- auc_features_x_cells
    return(list(
      csi_long = scFunctions::calculate_csi(auc_rankings, calc_extended = FALSE, verbose = FALSE),
      source = "scFunctions"
    ))
  }

  source(vendor_path, encoding = "UTF-8")
  auc_rankings <- AUCell::AUCell_buildRankings(auc_features_x_cells, plotStats = FALSE, nCores = 1)
  auc_rankings@assays@data@listData$AUC <- auc_features_x_cells
  list(
    csi_long = calculate_csi(auc_rankings, calc_extended = FALSE, verbose = FALSE),
    source = "vendor"
  )
}
