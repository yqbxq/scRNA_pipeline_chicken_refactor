normalize_key <- function(x) {
  toupper(trimws(as.character(x)))
}

# Coarse orthology grouping used in summaries and quality grading.
orthology_bucket <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("one2one$", x) ~ "one2one",
    grepl("one2many$", x) ~ "one2many",
    grepl("many2many$", x) ~ "many2many",
    TRUE ~ "other"
  )
}

# Stable priority used only for best-hit sorting within one chicken gene.
orthology_rank <- function(x) {
  dplyr::case_when(
    x == "ortholog_one2one" ~ 1L,
    x == "apparent_ortholog_one2one" ~ 2L,
    x == "ortholog_one2many" ~ 3L,
    x == "apparent_ortholog_one2many" ~ 4L,
    x == "ortholog_many2many" ~ 5L,
    TRUE ~ 99L
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

classify_pair_quality <- function(df) {
  bucket <- orthology_bucket(df$orthology_type)
  df$pair_quality <- dplyr::case_when(
    bucket == "one2one" & !is.na(df$orthology_confidence) & df$orthology_confidence == 1 ~ "gold",
    bucket == "one2one" ~ "silver",
    bucket %in% c("one2many", "many2many") ~ "ambiguous",
    TRUE ~ "ambiguous"
  )
  df
}
