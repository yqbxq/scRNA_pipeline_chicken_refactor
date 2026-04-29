# Original calculate_csi implementation adapted from scFunctions source usage.
calculate_csi <- function(regulonAUC, calc_extended = FALSE, verbose = FALSE) {
  compare_pcc <- function(vector_of_pcc, pcc) {
    pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
    if (pcc_larger == length(vector_of_pcc)) {
      return(0)
    }
    length(vector_of_pcc)
  }

  calc_csi_pair <- function(reg, reg2, pearson_cor) {
    test_cor <- pearson_cor[reg, reg2]
    total_n <- ncol(pearson_cor)
    pearson_cor_sub <- subset(pearson_cor, rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)
    sums <- apply(pearson_cor_sub, MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
    length(sums[sums == nrow(pearson_cor_sub)]) / total_n
  }

  regulonAUC_sub <- regulonAUC@assays@data@listData$AUC

  if (calc_extended) {
    regulonAUC_sub <- subset(regulonAUC_sub, grepl("extended", rownames(regulonAUC_sub)))
  } else {
    regulonAUC_sub <- subset(regulonAUC_sub, !grepl("extended", rownames(regulonAUC_sub)))
  }

  regulonAUC_sub <- t(regulonAUC_sub)
  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    tidyr::gather(regulon_2, pcc, -regulon_1) %>%
    dplyr::mutate(regulon_pair = paste(regulon_1, regulon_2, sep = "_"))

  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names) * length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow = num_of_calculations, ncol = 3), stringsAsFactors = FALSE)
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")

  idx <- 0
  for (reg in regulon_names) {
    if (verbose) {
      message(reg)
    }
    for (reg2 in regulon_names) {
      idx <- idx + 1
      csi_regulons[idx, ] <- c(reg, reg2, calc_csi_pair(reg, reg2, pearson_cor))
    }
  }

  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  csi_regulons
}
