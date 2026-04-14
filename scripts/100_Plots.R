#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(stringr)
  library(RColorBrewer)
})

# -------------------- argument parsing --------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  if (name %in% args) {
    args[which(args == name) + 1]
  } else {
    default
  }
}

metrics_path      <- get_arg("--metrics")
samples_path      <- get_arg("--samples")
bamdam_dir        <- get_arg("--bamdam_dir")
mmseqs_dir        <- get_arg("--mmseqs_dir")
metadata_path     <- get_arg("--metadata")
outdir            <- get_arg("--outdir", ".")
bamdam_plot_mode  <- tolower(get_arg("--bamdam_plot", "heatmap"))  # heatmap | bubble | both
min_reads_arg     <- get_arg("--min_reads", "1")

min_reads <- suppressWarnings(as.numeric(min_reads_arg))
if (is.na(min_reads) || min_reads < 0) min_reads <- 1

# damage threshold in percent (default 5)
damage_thr_arg <- get_arg("--damage_threshold", Sys.getenv("PLOTS_DAMAGE_THRESHOLD", "5"))
damage_threshold <- suppressWarnings(as.numeric(damage_thr_arg))
if (is.na(damage_threshold) || damage_threshold < 0) damage_threshold <- 5

# taxa per bamdam plot panel (0 = no splitting)
taxa_per_plot_arg <- get_arg("--taxa_per_plot", Sys.getenv("BAMDAM_TAXA_PER_PLOT", "0"))
taxa_per_plot <- suppressWarnings(as.integer(taxa_per_plot_arg))
if (is.na(taxa_per_plot) || taxa_per_plot < 0) taxa_per_plot <- 0

message("[100_Plots.R] damage_threshold(%): ", damage_threshold)
message("[100_Plots.R] taxa_per_plot: ", taxa_per_plot)
# --- NEW: taxa exclusion list for HEATMAP only ---
# Accepts comma / semicolon / whitespace separated names (exact matches).
exclude_taxa_arg <- get_arg("--exclude_taxa", Sys.getenv("PLOTS_EXCLUDE_TAXA", ""))
exclude_taxa <- stringr::str_split(exclude_taxa_arg, "[,;[:space:]]+")[[1]]
exclude_taxa <- exclude_taxa[nzchar(exclude_taxa)]
exclude_taxa <- stringr::str_trim(exclude_taxa)

message("[100_Plots.R] exclude_taxa (heatmap only) count: ", length(exclude_taxa))


# --- NEW: taxa evolution plot (optional) ---
# A text file with one taxon (genus name) per line.
taxa_trend_file <- get_arg("--taxa_trend_file", Sys.getenv("PLOTS_TAXA_TREND_FILE", ""))
taxa_trend_file <- stringr::str_trim(taxa_trend_file)

if (nzchar(taxa_trend_file)) {
  message("[100_Plots.R] taxa_trend_file: ", taxa_trend_file)
} else {
  message("[100_Plots.R] taxa_trend_file: NONE")
}

if (!bamdam_plot_mode %in% c("heatmap", "bubble", "both")) {
  message("[100_Plots.R] WARNING: invalid --bamdam_plot='", bamdam_plot_mode,
          "'. Using 'heatmap'. Allowed: heatmap, bubble, both.")
  bamdam_plot_mode <- "heatmap"
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("[100_Plots.R] outdir: ", outdir)
message("[100_Plots.R] bamdam_plot_mode: ", bamdam_plot_mode)
message("[100_Plots.R] min_reads: ", min_reads)

# -------------------- helper: read sample list --------------------
read_sample_list <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(character(0))
  }
  samp <- readr::read_lines(path)
  samp <- samp[samp != ""]

  # If the sample list contains FASTQ paths, convert to sample IDs.
  normalize_one <- function(x) {
    if (is.na(x) || !nzchar(x)) return(x)
    b <- basename(x)

    # strip compression + fastq suffixes
    b <- sub("\\.(fastq|fq)\\.gz$", "", b, ignore.case = TRUE)
    b <- sub("\\.(fastq|fq)$", "", b, ignore.case = TRUE)

    # strip common read mate / chunk suffixes
    b <- sub("(_R[12]_[0-9]+)$", "", b, ignore.case = TRUE)  # SAMPLE_R1_001
    b <- sub("(_R[12])$", "", b, ignore.case = TRUE)         # SAMPLE_R1

    b
  }

  unique(vapply(samp, normalize_one, character(1)))
}

# ================================================================
# 1) METRICS DOT PLOT (always)
# ================================================================
make_metrics_plot <- function(metrics_path, samples_path, outdir) {
  message("[metrics] starting metrics plot")

  if (is.null(metrics_path) || !file.exists(metrics_path)) {
    message("[metrics] metrics file not found: ", metrics_path, " - skipping metrics plot.")
    return(invisible(NULL))
  }

  primary_samples <- read_sample_list(samples_path)
  message("[metrics] samples in primary list: ",
          ifelse(length(primary_samples) == 0, "NONE",
                 paste(primary_samples, collapse = ", ")))

  m <- suppressMessages(read_tsv(metrics_path, show_col_types = FALSE))
  message("[metrics] metrics rows: ", nrow(m), " cols: ", ncol(m))

  if (!"sample" %in% names(m)) {
    message("[metrics] WARNING: metrics TSV has no 'sample' column - skipping metrics plot.")
    return(invisible(NULL))
  }


  # Optional: restrict metrics to the selected mapping/database tag
  map_tag <- Sys.getenv("MAP_LAST_DB_TAG", "")
  if (nzchar(map_tag) && "Database_name" %in% names(m)) {
    m <- m %>% filter(.data$Database_name == map_tag)
  }

  if (length(primary_samples) > 0) {
    m <- m %>% filter(sample %in% primary_samples)
  }
  if (nrow(m) == 0L) {
    message("[metrics] no samples left after filtering - skipping.")
    return(invisible(NULL))
  }

  # Case 1: already long format with a step column
  step_col  <- intersect(c("step", "Steps", "pipeline_step", "Step", "STEPS"), names(m))[1]
  reads_col <- intersect(c("reads", "read_count", "n_reads", "Reads", "READS"), names(m))[1]

  if (!is.na(step_col) && !is.na(reads_col)) {
    df <- m %>%
      mutate(
        step  = as.factor(.data[[step_col]]),
        reads = as.numeric(.data[[reads_col]])
      ) %>%
      select(sample, step, reads)
  } else {
    # Case 2: wide format (one row per sample; step columns are numeric)
    num_cols <- names(m)[vapply(m, is.numeric, logical(1))]
    num_cols <- setdiff(num_cols, "sample")

    if (length(num_cols) == 0) {
      message("[metrics] WARNING: couldn't find step+reads columns and no numeric step columns to pivot - skipping metrics plot.")
      message("[metrics] columns: ", paste(names(m), collapse = ", "))
      return(invisible(NULL))
    }

    # Keep the column order as the step order
    df <- m %>%
      select(sample, all_of(num_cols)) %>%
      pivot_longer(cols = all_of(num_cols), names_to = "step", values_to = "reads") %>%
      mutate(
        step  = factor(step, levels = num_cols),
        reads = as.numeric(reads)
      )
  }

  # Drop missing/zero reads safely
  df <- df %>% filter(!is.na(reads))
  if (nrow(df) == 0L) {
    message("[metrics] no numeric reads to plot - skipping.")
    return(invisible(NULL))
  }

  # Order samples along x-axis as in the primary sample list (when provided)
  if (length(primary_samples) > 0) {
    df <- df %>% mutate(sample = factor(sample, levels = primary_samples))
  } else {
    df <- df %>% mutate(sample = factor(sample))
  }



  p <- ggplot(df, aes(x = sample, y = reads, colour = step)) +
    geom_point(size = 2, alpha = 0.8, position = position_dodge(width = 0.7)) +
    scale_y_log10() +
    coord_flip() +
    labs(x = "Sample", y = "Number of reads", colour = "Step") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x       = element_text(size = 9),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      legend.key.height = unit(0.6, "cm")
    )

  pdf_file <- file.path(outdir, "reads_per_step_dots.pdf")
  png_file <- file.path(outdir, "reads_per_step_dots.png")

  message("[metrics] saving: ", pdf_file)
  ggsave(pdf_file, p, width = 7, height = 5)

  message("[metrics] saving: ", png_file)
  ggsave(png_file, p, width = 7, height = 5, dpi = 300)

  message("[metrics] metrics plot written.")
}

# ================================================================
# 2) BAMDAM ABUNDANCE PLOTS (heatmap OR bubble)
# ================================================================
make_bamdam_abundance_plots <- function(bamdam_dir, samples_path, metadata_path, outdir, min_reads = 1, plot_mode = "heatmap", taxa_per_plot = 0) {
  message("[bamdam] starting bamdam plot(s)")

  map_tag <- Sys.getenv("MAP_LAST_DB_TAG", "")

  if (is.null(bamdam_dir) || !dir.exists(bamdam_dir)) {
    message("[bamdam] directory not found at: ", bamdam_dir, " - skipping.")
    return(invisible(NULL))
  }
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    message("[bamdam] metadata file not found at: ", metadata_path, " - skipping.")
    return(invisible(NULL))
  }

  primary_samples <- read_sample_list(samples_path)
  message("[bamdam] samples in primary list: ",
          ifelse(length(primary_samples) == 0, "NONE",
                 paste(primary_samples, collapse = ", ")))

  # bamdam files: *.tsv with TaxName / TotalReads / taxpath
  files <- list.files(bamdam_dir, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)

  # If MAP_LAST_DB_TAG is set, pre-filter to files inside <SAMPLE>_<TAG>/ to avoid mixing runs
  if (nzchar(map_tag)) {
    tag_needle <- paste0("_", map_tag, .Platform$file.sep)
    files <- files[grepl(tag_needle, files, fixed = TRUE)]
  }

  message("[bamdam] found ", length(files), " bamdam *.tsv files")
  if (length(files) > 0) {
    message("[bamdam] example tsv: ", files[[1]])
  }
  if (length(files) == 0L) {
    message("[bamdam] no bamdam files found, skipping.")
    return(invisible(NULL))
  }

  # ---- metadata ----
  meta <- suppressMessages(read_tsv(metadata_path, show_col_types = FALSE))
  if (!"sample" %in% names(meta)) stop("[bamdam] ERROR: metadata TSV must contain a 'sample' column.")

  age_col <- intersect(c("age_ka", "age", "ka"), names(meta))[1]
  if (is.na(age_col)) stop("[bamdam] ERROR: metadata TSV must have an age column named one of: age_ka, age, ka")

  # helper to get rank from the first element of taxpath
  get_rank <- function(tp) {
    if (is.na(tp)) return(NA_character_)
    s <- gsub('^"|"$', "", tp)
    first <- strsplit(s, ";", fixed = TRUE)[[1]][1]
    parts <- strsplit(first, ":", fixed = TRUE)[[1]]
    if (length(parts) >= 3) parts[3] else NA_character_
  }

  # helper: assign clade for sorting (Viridiplantae first, then Opisthokonta, then other)
  get_clade <- function(tp) {
    if (is.na(tp) || !nzchar(tp)) return("Other")
    if (grepl("Viridiplantae", tp, ignore.case = TRUE)) return("Viridiplantae")
    if (grepl("Opisthokonta", tp, ignore.case = TRUE))  return("Opisthokonta")
    "Other"
  }

  # ---- read bamdam tsvs ----
  bam_list <- purrr::map(files, function(f) {
    x <- suppressMessages(read_tsv(f, show_col_types = FALSE))

    needed <- c("TaxName", "TotalReads", "taxpath", "Damage+1")
    missing <- setdiff(needed, names(x))
    if (length(missing) > 0) {
      stop("Bamdam file ", f, " is missing columns: ", paste(missing, collapse = ", "))
    }

    # Derive sample + database tag from the folder structure:
    # bamdam_dir/<sample>/<sample>_<tag>/.../*.tsv
    f_norm  <- normalizePath(f, winslash = "/", mustWork = FALSE)
    bd_norm <- normalizePath(bamdam_dir, winslash = "/", mustWork = FALSE)
    rel     <- sub(paste0("^", bd_norm, "/?"), "", f_norm)
    parts   <- strsplit(rel, "/", fixed = TRUE)[[1]]

    sample_id <- if (length(parts) >= 1) parts[1] else NA_character_
    tag_dir   <- if (length(parts) >= 2) parts[2] else NA_character_

    db_tag <- NA_character_
    if (!is.na(sample_id) && !is.na(tag_dir) && startsWith(tag_dir, paste0(sample_id, "_"))) {
      db_tag <- sub(paste0("^", sample_id, "_"), "", tag_dir)
    }

    # If MAP_LAST_DB_TAG is set, keep only matching-tag files
    if (nzchar(map_tag) && (is.na(db_tag) || db_tag != map_tag)) {
      return(tibble())
    }

    ranks <- vapply(x$taxpath, get_rank, character(1))

    tibble(
      sample  = sample_id,
      taxon   = as.character(x$TaxName),
      rank    = ranks,
      reads   = as.numeric(x$TotalReads),
      damage_p1 = as.numeric(x$`Damage+1`),
      damage_m1 = as.numeric(x$`Damage-1`),
      taxpath = as.character(x$taxpath)
    ) %>%
      filter(rank %in% c("genus", "family"))
  })

  dat_all <- bind_rows(bam_list) %>% filter(reads > 0)
  message("[bamdam] rows (genus/family, all bamdam samples): ", nrow(dat_all))
  if (nrow(dat_all) == 0L) {
    message("[bamdam] all counts are zero; skipping.")
    return(invisible(NULL))
  }

  # ---- sample matching ----
  bam_samples  <- unique(dat_all$sample)
  meta_samples <- unique(meta$sample)

  overlap3 <- intersect(bam_samples, intersect(meta_samples, primary_samples))
  overlap2 <- intersect(bam_samples, meta_samples)

  if (length(overlap3) > 0) {
    used_samples <- overlap3
  } else if (length(overlap2) > 0) {
    used_samples <- overlap2
  } else {
    used_samples <- bam_samples
  }

  meta_filt <- meta %>% filter(sample %in% used_samples)
  if (nrow(meta_filt) == 0L) {
    message("[bamdam] WARNING: no metadata rows match bamdam samples; proceeding without ages (sample-only ordering).")
    meta_filt <- tibble(sample = used_samples)
    meta_filt[[age_col]] <- NA_real_
  }

  meta_filt <- meta_filt %>% arrange(.data[[age_col]])
  sample_order <- meta_filt$sample


  # ---- auto-tune top sample labels based on number of samples ----
  n_samples <- length(sample_order)

  sample_name_size <- dplyr::case_when(
    n_samples <= 20 ~ 2.8,
    n_samples <= 30 ~ 2.4,
    n_samples <= 40 ~ 2.0,
    n_samples <= 60 ~ 1.4,
    TRUE            ~ 0.8
  )
  sample_age_size <- dplyr::case_when(
    n_samples <= 20 ~ 2.6,
    n_samples <= 30 ~ 2.2,
    n_samples <= 40 ~ 1.9,
    n_samples <= 60 ~ 1.6,
    TRUE            ~ 1.3
  )

  sample_angle <- dplyr::case_when(
    n_samples <= 20 ~ 0,
    n_samples <= 50 ~ 45,
    TRUE            ~ 60
  )
  sample_hjust <- ifelse(sample_angle == 0, 0.5, 0)

  # Keep sample names just above the plotting area, and ages below it.
  # Using small fixed offsets avoids the left-shift seen with large negative vjust on rotated labels.
  name_vjust <- dplyr::case_when(
    sample_angle == 0  ~ -0.35,
    sample_angle == 45 ~ -0.20,
    TRUE               ~ -0.12
  )
  age_vjust <- 1.55
  age_angle <- 0
  age_hjust <- 0.5

  top_margin_lines <- dplyr::case_when(
    sample_angle == 0  ~ 3.2,
    sample_angle == 45 ~ 4.8,
    TRUE               ~ 5.6
  )
  bottom_margin_lines <- 3.4
  plot_width <- min(18, max(9, 6 + 0.18 * n_samples))

  dat_all <- dat_all %>% filter(sample %in% sample_order)
  if (nrow(dat_all) == 0L) {
    message("[bamdam] no data left after sample filtering; skipping.")
    return(invisible(NULL))
  }

  # ---- sample_type + palette ----
  sample_type_col <- intersect(c("sample_type", "sampleType", "type", "sample_category", "sample_group"),
                               names(meta_filt))[1]

  if (!is.na(sample_type_col)) {
    sample_types <- factor(meta_filt[[sample_type_col]][match(sample_order, meta_filt$sample)])
  } else {
    sample_types <- factor(rep("unknown", length(sample_order)))
  }

  ntypes <- length(levels(sample_types))
  base_n <- max(3, min(8, ntypes))
  base_cols <- RColorBrewer::brewer.pal(base_n, "Set2")
  sample_type_cols <- if (ntypes <= length(base_cols)) base_cols[seq_len(ntypes)] else grDevices::colorRampPalette(base_cols)(ntypes)

  #heat_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(256)
  heat_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu")[1:8])(256)
  #heat_colors[256] <- "#69aedf"

  ages <- meta_filt %>%
    dplyr::distinct(sample, .keep_all = TRUE) %>%
    dplyr::transmute(sample, age = .data[[age_col]])

  label_df_name <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = sample_order,
    sample_type = sample_types
  )
  age_labels <- ages$age[match(sample_order, ages$sample)]
  age_labels <- dplyr::if_else(is.na(age_labels), "", as.character(age_labels))

  label_df_age <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = age_labels
  )

  header_x <- -0.1
  header_hjust <- 1

  legend_counts <- c(1, 10, 100, 1000, 10000, 100000)
  legend_breaks <- log10(legend_counts)

  # ---- shared: build long df + taxon ordering for a rank ----
  build_rank_df <- function(rank_sel) {
    dat_rank <- dat_all %>% filter(rank == rank_sel)
    if (nrow(dat_rank) == 0L) return(NULL)

    mat <- dat_rank %>%
      group_by(taxon, sample) %>%
      summarise(reads = sum(reads), .groups = "drop") %>%
      mutate(sample = factor(sample, levels = sample_order)) %>%
      pivot_wider(names_from = sample, values_from = reads, values_fill = 0) %>%
      column_to_rownames("taxon") %>%
      as.matrix()

    tax_totals <- rowSums(mat)
    keep_taxa  <- tax_totals >= min_reads
    mat <- mat[keep_taxa, , drop = FALSE]
    if (nrow(mat) == 0L) return(NULL)

    tax_info <- dat_rank %>%
      group_by(taxon) %>%
      summarise(taxpath = dplyr::first(na.omit(taxpath)), .groups = "drop") %>%
      mutate(clade = vapply(taxpath, get_clade, character(1))) %>%
      filter(taxon %in% rownames(mat))

    tax_info$clade <- factor(tax_info$clade, levels = c("Viridiplantae", "Opisthokonta", "Other"))

    # alphabetic within clade block
    taxon_order <- tax_info %>%
      arrange(clade, taxon) %>%
      pull(taxon)

    df_long <- as.data.frame(mat) %>%
      rownames_to_column("taxon") %>%
      pivot_longer(-taxon, names_to = "sample", values_to = "reads") %>%
      left_join(ages, by = "sample") %>%
      mutate(
        sample    = factor(sample, levels = sample_order),
        taxon     = factor(taxon, levels = rev(taxon_order)),  # top-to-bottom in plot
        reads     = tidyr::replace_na(reads, 0),
        log_reads = log10(reads + 1)
      )


# ---- damage class per cell (sample-dependent) ----
# Use Damage+1 directly from TSV (fraction). Convert to percent for thresholding.
damage_cell <- dat_rank %>%
  filter(sample %in% sample_order) %>%
  group_by(taxon, sample) %>%
  summarise(
    reads = sum(reads, na.rm = TRUE),
    dmg   = mean(c(dplyr::first(damage_p1), dplyr::first(damage_m1)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(dmg_pct = 100 * dmg)


# ---- NEW (heatmap-only): per-taxon max damage across samples (percent) ----
taxon_max_damage <- damage_cell %>%
  group_by(taxon) %>%
  summarise(max_dmg_pct = max(dmg_pct, na.rm = TRUE), .groups = "drop")

# For each sample: mean damage (percent) of the 3 most abundant taxa with >5% damage
ref_by_sample <- damage_cell %>%
  filter(dmg_pct > damage_threshold) %>%
  group_by(sample) %>%
  arrange(desc(reads), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  summarise(ref_mean = mean(dmg_pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(ref_low = pmax(damage_threshold, ref_mean - 5), ref_high = ref_mean + 5)

damage_cell <- damage_cell %>%
  left_join(ref_by_sample, by = "sample") %>%
  mutate(
    cond_damage_gt5 = dmg_pct > damage_threshold,
    cond_near_ref   = ifelse(is.na(ref_mean), FALSE, (dmg_pct >= ref_low & dmg_pct <= ref_high)),
    damage_class = dplyr::case_when(
      cond_damage_gt5 & cond_near_ref ~ "green",
      xor(cond_damage_gt5, cond_near_ref) ~ "orange",
      TRUE ~ "red"
    )
  ) %>%
  select(taxon, sample, damage_class)

df_long <- df_long %>%
  mutate(taxon_chr = as.character(taxon), sample_chr = as.character(sample)) %>%
  left_join(
    damage_cell %>% transmute(taxon_chr = taxon, sample_chr = sample, damage_class),
    by = c("taxon_chr", "sample_chr")
  ) %>%
  mutate(
    damage_class = factor(tidyr::replace_na(damage_class, "red"),
                          levels = c("red", "orange", "green"))
  ) %>%
  select(-taxon_chr, -sample_chr)

    list(df = df_long, max_log = max(df_long$log_reads, na.rm = TRUE), taxon_max_damage = taxon_max_damage)
  }

  # ---- plot: heatmap ----
  plot_heatmap <- function(df_long, max_log, out_prefix) {
    max_log <- max(max_log, max(legend_breaks))
    min_log <- 0
    label_df_name_plot <- label_df_name %>% dplyr::mutate(y = Inf)
    label_df_age_plot  <- label_df_age  %>% dplyr::mutate(y = -Inf)

    df_long <- df_long %>%
      mutate(label = formatC(as.integer(round(reads)), format = "f", digits = 0, big.mark = ","))

    p <- ggplot(df_long, aes(x = sample, y = taxon, fill = log_reads)) +
      geom_tile(width = 1, height = 1) +
      geom_text(aes(label = label, colour = reads >= 1e6), size = 1.4, show.legend = FALSE) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
      scale_x_discrete(
        position = "top",
        limits   = sample_order,
        labels   = rep("", length(sample_order)),
        expand   = expansion(mult = 0, add = 0)
      ) +
      scale_y_discrete(name = NULL, expand = expansion(mult = 0, add = 0)) +
      scale_fill_gradientn(
        colours = heat_colors,
        limits  = c(min_log, max_log),
        breaks  = legend_breaks,
        labels  = formatC(legend_counts, format = "fg", big.mark = ","),
        name    = "Unique Reads"
      ) +
      labs(x = NULL, title = NULL) +
      coord_cartesian(clip = "off") +
      theme_bw(base_size = 8) +
      theme(
        panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.title.x.bottom = element_blank(),
        axis.ticks.x        = element_blank(),
        panel.grid          = element_blank(),
        legend.position     = "right",
        legend.title        = element_text(size = 7),
        legend.text         = element_text(size = 6),
        plot.margin         = grid::unit(c(top_margin_lines, 0.5, bottom_margin_lines, 2.6), "lines")
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               label = "Samples",
               hjust = 1.1, vjust = -0.05, size = 2.5) +
      annotate("text",
               x = -Inf, y = -Inf,
               label = "Age",
               hjust = 1.4, vjust = 1.25, size = 2.5)

    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_label(
          data = label_df_name_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label, fill = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size,
          linewidth = 0.15,
          label.r = grid::unit(0.08, "lines"),
          key_glyph = ggplot2::draw_key_rect
        ) +
        geom_text(
          data = label_df_age_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust,
          size = sample_age_size, colour = "black"
        ) +
        scale_fill_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 3, title.position = "top")
        )
    } else {
      message("[100_Plots.R] NOTE: package 'ggnewscale' not found - drawing colored NAME text only.")
      p <- p +
        geom_text(
          data = label_df_name_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label, colour = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size, show.legend = TRUE
        ) +
        geom_text(
          data = label_df_age_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust,
          size = sample_age_size, colour = "black"
        ) +
        scale_colour_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 2, title.position = "top")
        )
    }

    pdf_file <- file.path(outdir, paste0(out_prefix, ".pdf"))
    png_file <- file.path(outdir, paste0(out_prefix, ".png"))

    message("[bamdam] saving: ", pdf_file)
    ggsave(pdf_file, p, width = plot_width, height = 8)

    message("[bamdam] saving: ", png_file)
    ggsave(png_file, p, width = plot_width, height = 8, dpi = 300)
  }

  # ---- plot: bubble ----
  plot_bubble_damage <- function(df_long, max_log, out_prefix) {
    max_log <- max(max_log, max(legend_breaks))
    min_log <- 0
    label_df_name_plot <- label_df_name %>% dplyr::mutate(y = Inf)
    label_df_age_plot  <- label_df_age  %>% dplyr::mutate(y = -Inf)

    # bubble size based on log10(reads+1)
    df_plot <- df_long %>%
      dplyr::filter(reads > 0) %>%
      dplyr::mutate(label = formatC(as.integer(round(reads)), format = "f", digits = 0, big.mark = ","))

    p <- ggplot(df_plot, aes(x = sample, y = taxon)) +
      geom_point(aes(size = log_reads, fill = damage_class), shape = 21, colour = "black", stroke = 0.15, alpha = 0.9) +
      geom_text(aes(label = label), colour = "black", size = 1.2) +
      scale_x_discrete(
        position = "top",
        limits   = sample_order,
        labels   = rep("", length(sample_order)),
        expand   = expansion(add = c(0.6, 0.6))
      ) +
      scale_y_discrete(name = NULL) +
      scale_fill_manual(
        name   = "Damage",
        values = c(red = "#D55E00", orange = "#F0E442", green = "#009E73"),
        breaks = c("green","orange","red"),
        labels = c("Confident", "Require investigation", "Insufficiently damaged"),
        guide  = guide_legend(order = 1, title.position = "top")
      ) +
      scale_size_continuous(
        range  = c(0.2, 7),
        limits = c(min_log, max_log),
        breaks = legend_breaks,
        labels = formatC(legend_counts, format = "fg", big.mark = ","),
        name   = "Unique Reads",
        guide  = guide_legend(order = 2, title.position = "top")
      ) +
      labs(x = NULL, title = NULL) +
      coord_cartesian(clip = "off") +
      theme_bw(base_size = 8) +
      theme(
        panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.title.x.bottom = element_blank(),
        axis.ticks.x        = element_blank(),
        panel.grid          = element_blank(),
        legend.position     = "right",
        legend.title        = element_text(size = 8),
        legend.text         = element_text(size = 7),
        plot.margin         = grid::unit(c(top_margin_lines, 0.5, bottom_margin_lines, 2.6), "lines")
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               label = "Samples",
               hjust = 1.1, vjust = -0.05, size = 2.5) +
      annotate("text",
               x = -Inf, y = -Inf,
               label = "Age",
               hjust = 1.4, vjust = 1.25, size = 2.5)

    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_label(
          data = label_df_name_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label, fill = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size,
          linewidth = 0.15,
          label.r = grid::unit(0.08, "lines"),
          key_glyph = ggplot2::draw_key_rect
        ) +
        geom_text(
          data = label_df_age_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust,
          size = sample_age_size, colour = "black"
        ) +
        scale_fill_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 3, title.position = "top")
        )
    } else {
      message("[100_Plots.R] NOTE: package 'ggnewscale' not found - drawing colored NAME text only.")
      p <- p +
        geom_text(
          data = label_df_name_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label, colour = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size, show.legend = TRUE
        ) +
        geom_text(
          data = label_df_age_plot, inherit.aes = FALSE,
          aes(x = sample, y = y, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust,
          size = sample_age_size, colour = "black"
        ) +
        scale_colour_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 3, title.position = "top")
        )
    }

    pdf_file <- file.path(outdir, paste0(out_prefix, ".pdf"))
    png_file <- file.path(outdir, paste0(out_prefix, ".png"))

    message("[bamdam] saving: ", pdf_file)
    ggsave(pdf_file, p, width = plot_width, height = 8)

    message("[bamdam] saving: ", png_file)
    ggsave(png_file, p, width = plot_width, height = 8, dpi = 300)
  }

plot_bubble_reads <- function(df_long, max_log, out_prefix) {
  legend_counts <- c(1, 10, 100, 1000, 10000, 100000)
  legend_breaks <- log10(legend_counts)

  max_log <- max(max_log, max(legend_breaks))
  min_log <- 0
  label_df_name_plot <- label_df_name %>% dplyr::mutate(y = Inf)
  label_df_age_plot  <- label_df_age  %>% dplyr::mutate(y = -Inf)

  df_plot <- df_long %>%
    dplyr::filter(reads > 0) %>%
    dplyr::mutate(
      label = formatC(as.integer(round(reads)), format = "f", digits = 0, big.mark = ","),
      high_read_label = reads >= 1e6
    )

  p <- ggplot(df_plot, aes(x = sample, y = taxon)) +
    geom_point(
      aes(size = log_reads, fill = log_reads),
      shape = 21, colour = "black", stroke = 0.15, alpha = 0.9
    ) +
    geom_text(
      aes(label = label, colour = high_read_label),
      size = 1.2,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "black")) +
    scale_x_discrete(
      position = "top",
      limits   = sample_order,
      labels   = rep("", length(sample_order)),
      expand   = expansion(add = c(0.6, 0.6))
    ) +
    scale_y_discrete(name = NULL) +
    scale_fill_gradientn(
      colours = heat_colors,
      limits  = c(min_log, max_log),
      breaks  = legend_breaks,
      labels  = formatC(legend_counts, format = "fg", big.mark = ","),
      name    = "Unique Reads",
      guide   = guide_colorbar(order = 1, title.position = "top")
    ) +
    scale_size_continuous(
      range  = c(0.2, 7),
      limits = c(min_log, max_log),
      breaks = legend_breaks,
      labels = formatC(legend_counts, format = "fg", big.mark = ","),
      name   = "Unique Reads",
      guide  = guide_legend(order = 2, title.position = "top")
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 9) +
    theme(
      panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.title.x.bottom = element_blank(),
      axis.ticks.x        = element_blank(),
      panel.grid          = element_blank(),
      legend.position     = "right",
      legend.title        = element_text(size = 8),
      legend.text         = element_text(size = 7),
      plot.margin         = grid::unit(c(top_margin_lines, 0.5, bottom_margin_lines, 2.6), "lines")
    ) +
    annotate("text",
      x = -Inf, y = Inf,
      label = "Samples",
      hjust = 1.1, vjust = -0.05, size = 2.5
    ) +
    annotate("text",
      x = -Inf, y = -Inf,
      label = "Age",
      hjust = 1.4, vjust = 1.25, size = 2.5
    )

  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_label(
        data = label_df_name_plot, inherit.aes = FALSE,
        aes(x = sample, y = y, label = label, fill = sample_type),
        vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size,
        linewidth = 0.15,
        label.r = grid::unit(0.08, "lines"),
        key_glyph = ggplot2::draw_key_rect
      ) +
      geom_text(
        data = label_df_age_plot, inherit.aes = FALSE,
        aes(x = sample, y = y, label = label),
        vjust = age_vjust, angle = age_angle, hjust = age_hjust,
        size = sample_age_size, colour = "black"
      ) +
      scale_fill_manual(
        name = "Sample type",
        values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
        guide = guide_legend(order = 3, title.position = "top")
      )
  } else {
    p <- p +
      geom_text(
        data = label_df_name_plot, inherit.aes = FALSE,
        aes(x = sample, y = y, label = label),
        vjust = name_vjust, angle = sample_angle, hjust = sample_hjust,
        size = sample_name_size, colour = "black", show.legend = FALSE
      ) +
      geom_text(
        data = label_df_age_plot, inherit.aes = FALSE,
        aes(x = sample, y = y, label = label),
        vjust = age_vjust, angle = age_angle, hjust = age_hjust,
        size = sample_age_size, colour = "black"
      )
  }

  pdf_file <- file.path(outdir, paste0(out_prefix, ".pdf"))
  png_file <- file.path(outdir, paste0(out_prefix, ".png"))

  message("[bamdam] saving: ", pdf_file)
  ggsave(pdf_file, p, width = plot_width, height = 8)

  message("[bamdam] saving: ", png_file)
  ggsave(png_file, p, width = plot_width, height = 8, dpi = 300)
}


  # ---- run genus + family ----
  for (rank_sel in c("genus", "family")) {
    built <- build_rank_df(rank_sel)
    if (is.null(built)) {
      message("[bamdam] no data for rank=", rank_sel, " after filtering; skipping.")
      next
    }
    df_long <- built$df
    max_log <- built$max_log
    taxon_max_damage <- built$taxon_max_damage

    # ---- split taxa into multiple panels if requested ----
    tax_levels <- levels(df_long$taxon)
    if (is.null(tax_levels) || length(tax_levels) == 0L) {
      tax_levels <- unique(as.character(df_long$taxon))
    }

    if (taxa_per_plot > 0 && length(tax_levels) > taxa_per_plot) {
      idx <- seq_along(tax_levels)
      chunk_ids <- ceiling(idx / taxa_per_plot)
      chunks <- split(tax_levels, chunk_ids)
    } else {
      chunks <- list(tax_levels)
    }

    for (ci in seq_along(chunks)) {
      taxa_chunk <- chunks[[ci]]

      df_chunk <- df_long %>%
        filter(.data$taxon %in% taxa_chunk) %>%
        mutate(taxon = factor(as.character(taxon), levels = taxa_chunk))

      part_lab <- NULL
      if (taxa_per_plot > 0) {
        start_i <- (ci - 1L) * taxa_per_plot + 1L
        end_i   <- start_i + length(taxa_chunk) - 1L
        part_lab <- sprintf("part%02d_%04d-%04d", ci, start_i, end_i)
      }

      if (plot_mode %in% c("heatmap", "both")) {
        out_prefix <- paste0("bamdam_", rank_sel, "_heatmap", if (!is.null(part_lab)) paste0(".", part_lab) else "")

        # ---- NEW: HEATMAP ONLY filter ----
        # Keep taxa whose max damage across samples is >= damage_threshold.
        keep_taxa <- taxon_max_damage %>%
          filter(max_dmg_pct >= damage_threshold) %>%
          pull(taxon) %>%
          as.character()

        df_heat <- df_chunk %>%
          filter(as.character(taxon) %in% keep_taxa)

        # Optional: user-specified exclusion list (exact matches on TaxName)
        if (length(exclude_taxa) > 0) {
          df_heat <- df_heat %>% filter(!(as.character(taxon) %in% exclude_taxa))
        }

        if (nrow(df_heat) == 0L) {
          message("[bamdam] heatmap: no taxa left after damage-threshold filter (", damage_threshold, "%); skipping: ", out_prefix)
        } else {
          # Drop unused factor levels to keep the y-axis clean.
          df_heat <- df_heat %>%
            mutate(taxon = factor(as.character(taxon), levels = levels(df_chunk$taxon))) %>%
            droplevels()

          plot_heatmap(df_heat, max_log, out_prefix)
        }
      }
      if (plot_mode %in% c("bubble", "both")) {
        out_prefix <- paste0("bamdam_", rank_sel, "_bubbleplot", if (!is.null(part_lab)) paste0(".", part_lab) else "")
        plot_bubble_damage(df_chunk, max_log, paste0(out_prefix, "_damage"))
        plot_bubble_reads(df_chunk, max_log, paste0(out_prefix, "_reads"))
      }
    }
  }

  message("[bamdam] done.")
}


# ================================================================
# 3) MMSEQS2 EVALUATION BUBBLE PLOT (optional)
# ================================================================
make_mmseqs_evaluation_bubbleplot <- function(mmseqs_dir, samples_path, metadata_path, outdir, taxa_per_plot = 0) {
  message("[mmseqs_eval] starting mmseqs evaluation plot")

  map_tag <- Sys.getenv("MAP_LAST_DB_TAG", "")

  if (is.null(mmseqs_dir) || !dir.exists(mmseqs_dir)) {
    message("[mmseqs_eval] mmseqs_dir not found at: ", mmseqs_dir, " - skipping.")
    return(invisible(NULL))
  }
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    message("[mmseqs_eval] metadata file not found at: ", metadata_path, " - skipping.")
    return(invisible(NULL))
  }

  primary_samples <- read_sample_list(samples_path)
  message("[mmseqs_eval] samples in primary list: ",
          ifelse(length(primary_samples) == 0, "NONE",
                 paste(primary_samples, collapse = ", ")))

  files <- list.files(mmseqs_dir, pattern = "bamdam_mmseqs\\.evaluation\\.summary\\.tsv$", full.names = TRUE, recursive = TRUE)

  # If MAP_LAST_DB_TAG is set, pre-filter to files inside <SAMPLE>_<TAG>/ to avoid mixing runs
  if (nzchar(map_tag)) {
    tag_needle <- paste0("_", map_tag, .Platform$file.sep)
    files <- files[grepl(tag_needle, files, fixed = TRUE)]
  }

  message("[mmseqs_eval] found ", length(files), " evaluation *.tsv files")
  if (length(files) > 0) message("[mmseqs_eval] example: ", files[[1]])
  if (length(files) == 0L) {
    message("[mmseqs_eval] no evaluation files found; skipping.")
    return(invisible(NULL))
  }

  meta <- suppressMessages(read_tsv(metadata_path, show_col_types = FALSE))
  if (!"sample" %in% names(meta)) stop("[mmseqs_eval] ERROR: metadata TSV must contain a 'sample' column.")
  meta <- meta %>% mutate(sample = trimws(as.character(sample)))

  age_col <- intersect(c("age_ka", "age", "ka"), names(meta))[1]
  if (is.na(age_col)) stop("[mmseqs_eval] ERROR: metadata TSV must have an age column named one of: age_ka, age, ka")

  # ---- helper: clade from taxpath (matches bamdam ordering blocks) ----
  get_clade <- function(tp) {
    if (is.na(tp) || !nzchar(tp)) return("Other")
    if (grepl("Viridiplantae", tp, ignore.case = TRUE)) return("Viridiplantae")
    if (grepl("Opisthokonta", tp, ignore.case = TRUE))  return("Opisthokonta")
    "Other"
  }

  # ---- read per-sample evaluation files ----
  eval_list <- purrr::map(files, function(f) {
    x <- suppressMessages(read_tsv(f, show_col_types = FALSE))

    needed <- c("exp_kingdom_name", "exp_genus_name", "n_queries", "same_genus", "same_family", "no_hit")
    missing <- setdiff(needed, names(x))
    if (length(missing) > 0) {
      stop("[mmseqs_eval] file ", f, " is missing columns: ", paste(missing, collapse = ", "))
    }

    # Derive sample + database tag from folder structure:
    # mmseqs_dir/<sample>/<sample>_<tag>/.../*.tsv
    f_norm  <- normalizePath(f, winslash = "/", mustWork = FALSE)
    ms_norm <- normalizePath(mmseqs_dir, winslash = "/", mustWork = FALSE)
    rel     <- sub(paste0("^", ms_norm, "/?"), "", f_norm)
    parts   <- strsplit(rel, "/", fixed = TRUE)[[1]]

    sample_id <- if (length(parts) >= 1) parts[1] else NA_character_
    sample_id <- trimws(sample_id)
    tag_dir   <- if (length(parts) >= 2) parts[2] else NA_character_

    db_tag <- NA_character_
    if (!is.na(sample_id) && !is.na(tag_dir) && startsWith(tag_dir, paste0(sample_id, "_"))) {
      db_tag <- sub(paste0("^", sample_id, "_"), "", tag_dir)
    }

    # If MAP_LAST_DB_TAG is set, keep only matching-tag files
    if (nzchar(map_tag) && (is.na(db_tag) || db_tag != map_tag)) {
      return(tibble())
    }

    # Ensure missing numeric columns exist (older files)
    if (!"same_order" %in% names(x))   x$same_order   <- 0
    if (!"same_kingdom" %in% names(x)) x$same_kingdom <- 0
    if (!"other" %in% names(x))        x$other        <- 0

    x %>%
      mutate(sample = sample_id) %>%
      mutate(exp_genus_taxid = if ("exp_genus_taxid" %in% names(.)) as.character(exp_genus_taxid) else NA_character_) %>%
      select(sample, exp_kingdom_name, exp_genus_taxid, exp_genus_name, n_queries, same_genus, same_family,
             same_order, same_kingdom, other, no_hit)
  })

  eval_raw <- bind_rows(eval_list) %>%
    mutate(sample = trimws(as.character(sample)))

  if (nrow(eval_raw) == 0L) {
    message("[mmseqs_eval] no rows read; skipping.")
    return(invisible(NULL))
  }

  # ---- sample matching (STRICT: do not fall back to plotting all eval samples) ----
  eval_samples <- unique(eval_raw$sample)
  meta_samples <- unique(meta$sample)

  overlap3 <- intersect(eval_samples, intersect(meta_samples, primary_samples))
  overlap2 <- intersect(eval_samples, meta_samples)

  if (length(overlap3) > 0) {
    used_samples <- overlap3
  } else if (length(overlap2) > 0) {
    used_samples <- overlap2
  } else {
    used_samples <- if (length(primary_samples) > 0) intersect(meta_samples, primary_samples) else meta_samples
  }

  meta_filt <- meta %>% filter(sample %in% used_samples)
  if (nrow(meta_filt) == 0L) {
    message("[mmseqs_eval] WARNING: no metadata rows match selected samples; skipping.")
    return(invisible(NULL))
  }
  meta_filt <- meta_filt %>% arrange(.data[[age_col]], .data$sample)
  sample_order <- meta_filt$sample

  # ---- auto-tune top sample labels based on number of samples ----
  n_samples <- length(sample_order)

  sample_name_size <- dplyr::case_when(
    n_samples <= 20 ~ 2.8,
    n_samples <= 30 ~ 2.4,
    n_samples <= 40 ~ 2.0,
    n_samples <= 60 ~ 1.4,
    TRUE            ~ 0.8
  )
  sample_age_size <- dplyr::case_when(
    n_samples <= 20 ~ 2.6,
    n_samples <= 30 ~ 2.2,
    n_samples <= 40 ~ 1.9,
    n_samples <= 60 ~ 1.6,
    TRUE            ~ 1.3
  )

  sample_angle <- dplyr::case_when(
    n_samples <= 20 ~ 0,
    n_samples <= 50 ~ 45,
    TRUE            ~ 60
  )
  sample_hjust <- ifelse(sample_angle == 0, 0.5, 0)

  name_vjust <- dplyr::case_when(
    sample_angle == 0  ~ -0.35,
    sample_angle == 45 ~ -0.20,
    TRUE               ~ -0.12
  )
  age_vjust <- 1.55
  age_angle <- 0
  age_hjust <- 0.5

  top_margin_lines <- dplyr::case_when(
    sample_angle == 0  ~ 3.2,
    sample_angle == 45 ~ 4.8,
    TRUE               ~ 5.6
  )
  bottom_margin_lines <- 3.4
  plot_width <- min(18, max(9, 6 + 0.18 * n_samples))

  ages <- meta_filt %>% dplyr::select(sample, age = dplyr::all_of(age_col))
  age_vec <- ages$age; names(age_vec) <- ages$sample

  # ---- sample_type + palette ----
  sample_type_col <- intersect(c("sample_type", "sampleType", "type", "sample_category", "sample_group"),
                               names(meta_filt))[1]

  if (!is.na(sample_type_col)) {
    sample_types <- factor(meta_filt[[sample_type_col]][match(sample_order, meta_filt$sample)])
  } else {
    sample_types <- factor(rep("unknown", length(sample_order)))
  }

  ntypes <- length(levels(sample_types))
  base_n <- max(3, min(8, ntypes))
  base_cols <- RColorBrewer::brewer.pal(base_n, "Set2")
  sample_type_cols <- if (ntypes <= length(base_cols)) base_cols[seq_len(ntypes)] else grDevices::colorRampPalette(base_cols)(ntypes)

  label_df_name <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = sample_order,
    sample_type = sample_types
  )
  label_df_age <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = ifelse(is.na(age_vec[sample_order]), "", format(age_vec[sample_order], trim = TRUE, scientific = FALSE))
  )

  # ---- attach bamdam genus read counts (defines plotted taxa) ----
  extract_rank <- function(tp) {
    if (is.na(tp) || !nzchar(tp)) return(NA_character_)
    first <- strsplit(tp, ";", fixed = TRUE)[[1]][1]
    parts <- strsplit(first, ":", fixed = TRUE)[[1]]
    if (length(parts) >= 3) parts[3] else NA_character_
  }

  get_bamdam_genus_reads <- function(sample_id) {
    if (!exists("bamdam_dir", inherits = TRUE) || is.null(bamdam_dir) || !dir.exists(bamdam_dir)) {
      return(tibble(sample = character(), exp_genus_taxid = character(), exp_genus_name = character(), bamdam_reads = numeric(), bamdam_taxpath = character()))
    }
    base_dir <- file.path(bamdam_dir, sample_id)
    if (!dir.exists(base_dir)) {
      return(tibble(sample = character(), exp_genus_taxid = character(), exp_genus_name = character(), bamdam_reads = numeric(), bamdam_taxpath = character()))
    }

    run_dir <- if (nzchar(map_tag)) {
      file.path(base_dir, paste0(sample_id, "_", map_tag))
    } else {
      cand <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
      cand <- cand[grepl(paste0("^", sample_id, "_"), basename(cand))]
      if (length(cand) > 0) cand[[1]] else base_dir
    }

    tsv <- file.path(run_dir, paste0(sample_id, ".tsv"))
    if (!file.exists(tsv)) {
      alt <- list.files(run_dir, pattern = "\\.tsv$", full.names = TRUE)
      if (length(alt) == 0) {
        return(tibble(sample = character(), exp_genus_taxid = character(), exp_genus_name = character(), bamdam_reads = numeric(), bamdam_taxpath = character()))
      }
      tsv <- alt[[1]]
    }

    bd <- tryCatch(
      readr::read_tsv(tsv, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
    if (is.null(bd)) {
      return(tibble(sample = character(), exp_genus_taxid = character(), exp_genus_name = character(), bamdam_reads = numeric(), bamdam_taxpath = character()))
    }

    need <- c("TaxNodeID", "TaxName", "TotalReads", "taxpath")
    if (!all(need %in% colnames(bd))) {
      return(tibble(sample = character(), exp_genus_taxid = character(), exp_genus_name = character(), bamdam_reads = numeric(), bamdam_taxpath = character()))
    }

    bd %>%
      mutate(
        TaxNodeID   = as.character(TaxNodeID),
        TaxName     = as.character(TaxName),
        TotalReads  = suppressWarnings(as.numeric(TotalReads)),
        rank        = vapply(taxpath, extract_rank, character(1))
      ) %>%
      filter(rank == "genus") %>%
      transmute(
        sample = sample_id,
        exp_genus_taxid = TaxNodeID,
        exp_genus_name  = TaxName,
        bamdam_reads    = TotalReads,
        bamdam_taxpath  = as.character(taxpath)
      )
  }

  bam_reads <- purrr::map_dfr(sample_order, get_bamdam_genus_reads) %>%
    filter(sample %in% sample_order) %>%
    mutate(
      exp_genus_taxid = as.character(exp_genus_taxid),
      exp_genus_name  = as.character(exp_genus_name),
      bamdam_reads    = suppressWarnings(as.numeric(bamdam_reads))
    ) %>%
    filter(!is.na(bamdam_reads) & bamdam_reads > 0) %>%
    group_by(sample, exp_genus_taxid, exp_genus_name) %>%
    summarise(
      bamdam_reads = sum(bamdam_reads, na.rm = TRUE),
      bamdam_taxpath = dplyr::first(na.omit(bamdam_taxpath)),
      .groups = "drop"
    )

  if (nrow(bam_reads) == 0L) {
    message("[mmseqs_eval] no bamdam genus reads found (or bamdam_dir missing); skipping.")
    return(invisible(NULL))
  }

  # Optional: apply min_reads filter as bamdam (across all samples)
  if (exists("min_reads", inherits = TRUE) && is.numeric(min_reads) && min_reads > 1) {
    keep_taxa <- bam_reads %>%
      group_by(exp_genus_name) %>%
      summarise(total = sum(bamdam_reads, na.rm = TRUE), .groups = "drop") %>%
      filter(total >= min_reads) %>%
      pull(exp_genus_name) %>%
      as.character()
    bam_reads <- bam_reads %>% filter(exp_genus_name %in% keep_taxa)
  }

  bam_reads <- bam_reads %>% mutate(exp_genus_name_key = tolower(trimws(exp_genus_name)))

  # ---- classify evaluation rows (MMSeqs2) ----
  eval_classed <- eval_raw %>%
    filter(sample %in% sample_order) %>%
    mutate(
      exp_genus_taxid     = as.character(exp_genus_taxid),
      exp_genus_name      = trimws(as.character(exp_genus_name)),
      exp_genus_name_key  = tolower(exp_genus_name),
      n_queries           = suppressWarnings(as.numeric(n_queries)),
      same_genus          = suppressWarnings(as.numeric(same_genus)),
      same_family         = suppressWarnings(as.numeric(same_family)),
      other               = suppressWarnings(as.numeric(other)),
      no_hit              = suppressWarnings(as.numeric(no_hit))
    ) %>%
    mutate(
      n_queries   = tidyr::replace_na(n_queries, 0),
      same_genus  = tidyr::replace_na(same_genus, 0),
      same_family = tidyr::replace_na(same_family, 0),
      other       = tidyr::replace_na(other, 0),
      no_hit      = tidyr::replace_na(no_hit, 0),
      hits                 = pmax(n_queries - no_hit, 0),
      hits_frac_total       = if_else(n_queries > 0, hits / n_queries, 0),
      same_genus_frac_hits  = if_else(hits > 0, same_genus / hits, 0),
      same_family_frac_hits = if_else(hits > 0, same_family / hits, 0),
      same_genus_frac_total  = if_else(n_queries > 0, same_genus / n_queries, 0),
      same_family_frac_total = if_else(n_queries > 0, same_family / n_queries, 0),
      eval_class = case_when(
        n_queries == 0 ~ "Require specific WGS genera MMSeqs2 investigation",
        hits_frac_total < 0.50 ~ "Require specific WGS genera MMSeqs2 investigation",
        hits > 0 & same_genus_frac_hits >= 0.40 & same_genus_frac_total >= 0.15 ~ "Confident for genus",
        hits > 0 & same_family_frac_hits >= (1/3) & same_family_frac_total >= 0.15 ~ "Confident for family\n but require investigation for genus",
        TRUE ~ "False assignation\n potential db bias or contamination"
      )
    ) %>%
    select(sample, exp_kingdom_name, exp_genus_taxid, exp_genus_name, exp_genus_name_key,
           n_queries, same_genus, same_family, other, no_hit, eval_class)

  # ---- join eval class onto bamdam genera; missing => Not tested ----
  df <- bam_reads %>%
    left_join(eval_classed %>% select(sample, exp_genus_taxid, eval_class),
              by = c("sample", "exp_genus_taxid"))

  df <- df %>%
    left_join(eval_classed %>% select(sample, exp_genus_name_key, eval_class) %>% rename(eval_class_by_name = eval_class),
              by = c("sample", "exp_genus_name_key")) %>%
    mutate(eval_class = dplyr::coalesce(eval_class, eval_class_by_name)) %>%
    select(-eval_class_by_name)

  df <- df %>%
    mutate(
      eval_class = if_else(is.na(eval_class) | !nzchar(eval_class), "Not tested for mmseqs", eval_class),
      eval_class = factor(eval_class, levels = c(
        "Confident for genus",
        "Confident for family\n but require investigation for genus",
        "False assignation\n potential db bias or contamination",
        "Require specific WGS genera MMSeqs2 investigation",
        "Not tested for mmseqs"
      )),
      clade = factor(vapply(bamdam_taxpath, get_clade, character(1)),
                     levels = c("Viridiplantae", "Opisthokonta", "Other"))
    )

  # ---- write merged file in plot folder (recreated each run) ----
  merged_tsv <- file.path(outdir, "mmseqs_evaluation_merged.tsv")
  suppressWarnings(write_tsv(df, merged_tsv))
  message("[mmseqs_eval] wrote merged file: ", merged_tsv)

  # ---- ordering of taxa ----
  taxon_order <- df %>%
    distinct(exp_genus_name, clade) %>%
    arrange(clade, exp_genus_name) %>%
    pull(exp_genus_name)

  df_plot <- df %>%
    mutate(
      sample = factor(sample, levels = sample_order),
      taxon  = factor(exp_genus_name, levels = rev(taxon_order)),
      log_reads = log10(bamdam_reads + 1),
      label = formatC(as.integer(round(bamdam_reads)), format = "f", digits = 0, big.mark = ",")
    )

  legend_counts <- c(1, 10, 100, 1000, 10000, 100000)
  legend_breaks <- log10(legend_counts)

  # ---- split taxa into multiple panels if requested ----
  tax_levels <- levels(df_plot$taxon)
  if (taxa_per_plot > 0 && length(tax_levels) > taxa_per_plot) {
    idx <- seq_along(tax_levels)
    chunk_ids <- ceiling(idx / taxa_per_plot)
    chunks <- split(tax_levels, chunk_ids)
  } else {
    chunks <- list(tax_levels)
  }

  for (ci in seq_along(chunks)) {
    taxa_chunk <- chunks[[ci]]

    df_chunk <- df_plot %>%
      dplyr::filter(.data$taxon %in% taxa_chunk) %>%
      dplyr::mutate(taxon = factor(as.character(taxon), levels = taxa_chunk)) %>%
      droplevels()

    part_lab <- NULL
    if (taxa_per_plot > 0) {
      start_i <- (ci - 1L) * taxa_per_plot + 1L
      end_i   <- start_i + length(taxa_chunk) - 1L
      part_lab <- sprintf("part%02d_%04d-%04d", ci, start_i, end_i)
    }

    p <- ggplot(df_chunk, aes(x = sample, y = taxon)) +
      geom_point(aes(size = log_reads, fill = eval_class), shape = 21, colour = "black", stroke = 0.15, alpha = 0.9) +
      geom_text(aes(label = label), colour = "black", size = 1.2) +
      scale_x_discrete(
        position = "top",
        limits   = sample_order,
        labels   = rep("", length(sample_order)),
        expand   = expansion(add = c(0.6, 0.6))
      ) +
      scale_y_discrete(name = NULL) +
      scale_fill_manual(
        name   = "MMSeqs2",
        values = c(
          "Confident for genus" = "#009E73",
          "Confident for family\n but require investigation for genus" = "#F0E442",
          "False assignation\n potential db bias or contamination" = "#D55E00",
          "Require specific WGS genera MMSeqs2 investigation" = "#0072B2",
          "Not tested for mmseqs" = "#999999"
        ),
        breaks = c(
          "Confident for genus",
          "Confident for family\n but require investigation for genus",
          "False assignation\n potential db bias or contamination",
          "Require specific WGS genera MMSeqs2 investigation",
          "Not tested for mmseqs"
        ),
        guide = guide_legend(order = 1, title.position = "top", override.aes = list(size = 3))
      ) +
      scale_size_continuous(
        range  = c(0.2, 7),
        limits = c(0, max(max(df_chunk$log_reads, na.rm = TRUE), max(legend_breaks))),
        breaks = legend_breaks,
        labels = formatC(legend_counts, format = "fg", big.mark = ","),
        name   = "Unique Reads",
        guide  = guide_legend(order = 2, title.position = "top")
      ) +
      labs(x = NULL, title = NULL) +
      coord_cartesian(clip = "off") +
      theme_bw(base_size = 9) +
      theme(
        panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.title.x.bottom = element_blank(),
        axis.ticks.x        = element_blank(),
        panel.grid          = element_blank(),
        legend.position     = "right",
        legend.title        = element_text(size = 8),
        legend.text         = element_text(size = 5),
        axis.text.y         = element_text(size = 6),
        axis.text.x.top     = element_blank(),
        plot.margin         = grid::unit(c(top_margin_lines, 0.5, bottom_margin_lines, 1.0), "lines")
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               label = "Samples",
               hjust = 1.1, vjust = -0.05, size = 2.5) +
      annotate("text",
               x = -Inf, y = -Inf,
               label = "Age",
               hjust = 1.4, vjust = 1.25, size = 2.5)

    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_label(
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, fill = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size,
          linewidth = 0.15,
          label.r = grid::unit(0.08, "lines"),
          key_glyph = ggplot2::draw_key_rect
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = -Inf, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust, size = sample_age_size, colour = "black",
          show.legend = FALSE
        ) +
        scale_fill_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 3, title.position = "top")
        )
    } else {
      p <- p +
        geom_text(
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, colour = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size, show.legend = TRUE
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = -Inf, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust, size = sample_age_size, colour = "black",
          show.legend = FALSE
        ) +
        scale_colour_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 3, title.position = "top")
        )
    }

    base_name <- paste0("mmseqs_evaluation_bubbleplot", if (!is.null(part_lab)) paste0(".", part_lab) else "")
    pdf_file <- file.path(outdir, paste0(base_name, ".pdf"))
    png_file <- file.path(outdir, paste0(base_name, ".png"))

    message("[mmseqs_eval] saving: ", pdf_file)
    ggsave(pdf_file, p, width = plot_width, height = 8.5)

    message("[mmseqs_eval] saving: ", png_file)
    ggsave(png_file, p, width = plot_width, height = 8.5, dpi = 300)
  }

  message("[mmseqs_eval] done.")
}


# ================================================================
# 4) TAXA EVOLUTION (optional; line plot per requested taxon)
# ================================================================
make_taxa_evolution_plots <- function(bamdam_dir, samples_path, metadata_path, metrics_path, outdir, taxa_file, min_reads = 1) {
  message("[taxa_evolution] starting taxa evolution plot(s)")

  if (is.null(taxa_file) || !nzchar(taxa_file) || !file.exists(taxa_file)) {
    message("[taxa_evolution] taxa file not provided / not found - skipping.")
    return(invisible(NULL))
  }
  if (is.null(bamdam_dir) || !dir.exists(bamdam_dir)) {
    message("[taxa_evolution] bamdam_dir not found - skipping.")
    return(invisible(NULL))
  }
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    message("[taxa_evolution] metadata file not found - skipping.")
    return(invisible(NULL))
  }

  # Read taxa list (one per line)
  taxa <- readr::read_lines(taxa_file)
  taxa <- taxa[!grepl("^\\s*#", taxa)]  # allow comments
  taxa <- stringr::str_trim(taxa)
  taxa <- taxa[nzchar(taxa)]
  taxa <- unique(taxa)

  if (length(taxa) == 0L) {
    message("[taxa_evolution] no taxa found in file - skipping.")
    return(invisible(NULL))
  }
  message("[taxa_evolution] taxa requested: ", paste(taxa, collapse = ", "))

  primary_samples <- read_sample_list(samples_path)

  meta <- suppressMessages(read_tsv(metadata_path, show_col_types = FALSE))
  if (!"sample" %in% names(meta)) stop("[taxa_evolution] ERROR: metadata TSV must contain a 'sample' column.")
  age_col <- intersect(c("age_ka", "age", "ka"), names(meta))[1]
  if (is.na(age_col)) stop("[taxa_evolution] ERROR: metadata TSV must have an age column named one of: age_ka, age, ka")


  # Optional: filter to one site (SITE_TAG), like other site-specific plots
  site_tag <- Sys.getenv("SITE_TAG", "")
  if (nzchar(site_tag)) {
    site_col <- intersect(
      c("site_tag", "SITE_TAG", "site", "Site", "site_name", "siteName", "Site_name", "site_id", "siteid"),
      names(meta)
    )[1]

    if (!is.na(site_col)) {
      meta <- meta %>% dplyr::filter(trimws(.data[[site_col]]) == site_tag)
      message("[taxa_evolution] SITE_TAG filter: ", site_tag, " (", nrow(meta), " metadata rows kept)")
    } else {
      message("[taxa_evolution] SITE_TAG is set (", site_tag, ") but metadata has no recognized site column; not filtering.")
    }
  }

  # Determine sample order (same strategy as other plots)
  meta_samples <- unique(meta$sample)
  used_samples <- meta_samples
  if (length(primary_samples) > 0) {
    overlap <- intersect(primary_samples, meta_samples)
    if (length(overlap) > 0) used_samples <- overlap
  }

  meta_filt <- meta %>% filter(sample %in% used_samples)
  meta_filt <- meta_filt %>% arrange(.data[[age_col]], .data$sample)
  sample_order <- meta_filt$sample

  if (length(sample_order) == 0L) {
    message("[taxa_evolution] no samples after filtering - skipping.")
    return(invisible(NULL))
  }

  # Auto-tune label sizes/rotation and output width (same scheme as your other plots)
  n_samples <- length(sample_order)

  sample_name_size <- dplyr::case_when(
    n_samples <= 20 ~ 2.8,
    n_samples <= 30 ~ 2.4,
    n_samples <= 40 ~ 2.0,
    n_samples <= 60 ~ 1.4,
    TRUE            ~ 0.8
  )
  sample_age_size <- dplyr::case_when(
    n_samples <= 20 ~ 2.6,
    n_samples <= 30 ~ 2.2,
    n_samples <= 40 ~ 1.9,
    n_samples <= 60 ~ 1.6,
    TRUE            ~ 1.3
  )

  sample_angle <- dplyr::case_when(
    n_samples <= 20 ~ 0,
    n_samples <= 50 ~ 45,
    TRUE            ~ 60
  )
  sample_hjust <- ifelse(sample_angle == 0, 0.5, 0)

  name_vjust <- dplyr::case_when(
    sample_angle == 0  ~ -0.35,
    sample_angle == 45 ~ -0.20,
    TRUE               ~ -0.12
  )
  age_vjust <- 1.55
  age_angle <- 0
  age_hjust <- 0.5

  top_margin_lines <- dplyr::case_when(
    sample_angle == 0  ~ 3.2,
    sample_angle == 45 ~ 4.8,
    TRUE               ~ 5.6
  )
  bottom_margin_lines <- 3.4
  plot_width <- min(18, max(9, 6 + 0.18 * n_samples))

  # Ages (as character labels)
  ages <- meta_filt %>% dplyr::select(sample, age = dplyr::all_of(age_col))
  age_vec <- ages$age
  names(age_vec) <- ages$sample

  label_df_name <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = sample_order
  )
  label_df_age <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = ifelse(is.na(age_vec[sample_order]), "", format(age_vec[sample_order], trim = TRUE, scientific = FALSE))
  )

  # ---- sample_type + palette (same strategy as bamdam/mmseqs top labels) ----
  sample_type_col <- intersect(
    c("sample_type", "sampleType", "type", "sample_category", "sample_group"),
    names(meta_filt)
  )[1]

  if (!is.na(sample_type_col)) {
    sample_types <- factor(meta_filt[[sample_type_col]][match(sample_order, meta_filt$sample)])
  } else {
    sample_types <- factor(rep("unknown", length(sample_order)))
  }

  base_n <- min(8, max(3, length(levels(sample_types))))
  base_cols <- RColorBrewer::brewer.pal(base_n, "Set2")
  ntypes <- length(levels(sample_types))
  sample_type_cols <- if (ntypes <= length(base_cols)) base_cols[seq_len(ntypes)] else grDevices::colorRampPalette(base_cols)(ntypes)

  label_df_name <- label_df_name %>% mutate(sample_type = sample_types)

  # ---- read bamdam genus counts for these samples (re-uses same parsing as MMseqs eval plot) ----
  map_tag <- Sys.getenv("MAP_LAST_DB_TAG", "")

  extract_rank <- function(tp) {
    if (is.na(tp) || !nzchar(tp)) return(NA_character_)
    first <- strsplit(tp, ";", fixed = TRUE)[[1]][1]
    parts <- strsplit(first, ":", fixed = TRUE)[[1]]
    if (length(parts) >= 3) parts[3] else NA_character_
  }

  get_bamdam_genus_reads <- function(sample_id) {
    base_dir <- file.path(bamdam_dir, sample_id)
    if (!dir.exists(base_dir)) {
      return(tibble(sample = character(), exp_genus_name = character(), bamdam_reads = numeric()))
    }

    run_dir <- if (nzchar(map_tag)) {
      file.path(base_dir, paste0(sample_id, "_", map_tag))
    } else {
      cand <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
      cand <- cand[grepl(paste0("^", sample_id, "_"), basename(cand))]
      if (length(cand) > 0) cand[[1]] else base_dir
    }

    tsv <- file.path(run_dir, paste0(sample_id, ".tsv"))
    if (!file.exists(tsv)) {
      alt <- list.files(run_dir, pattern = "\\.tsv$", full.names = TRUE)
      if (length(alt) == 0) {
        return(tibble(sample = character(), exp_genus_name = character(), bamdam_reads = numeric()))
      }
      tsv <- alt[[1]]
    }

    bd <- tryCatch(
      readr::read_tsv(tsv, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
    if (is.null(bd)) {
      return(tibble(sample = character(), exp_genus_name = character(), bamdam_reads = numeric()))
    }

    need <- c("TaxName", "TotalReads", "taxpath")
    if (!all(need %in% colnames(bd))) {
      return(tibble(sample = character(), exp_genus_name = character(), bamdam_reads = numeric()))
    }

    bd %>%
      mutate(
        TaxName    = as.character(TaxName),
        TotalReads = suppressWarnings(as.numeric(TotalReads)),
        rank       = vapply(taxpath, extract_rank, character(1))
      ) %>%
      filter(rank == "genus") %>%
      transmute(
        sample = sample_id,
        exp_genus_name = TaxName,
        bamdam_reads = TotalReads
      )
  }

  bam_reads <- purrr::map_dfr(sample_order, get_bamdam_genus_reads) %>%
    mutate(
      exp_genus_name = as.character(exp_genus_name),
      bamdam_reads   = suppressWarnings(as.numeric(bamdam_reads))
    ) %>%
    filter(!is.na(bamdam_reads)) %>%
    group_by(sample, exp_genus_name) %>%
    summarise(bamdam_reads = sum(bamdam_reads, na.rm = TRUE), .groups = "drop")

  if (nrow(bam_reads) == 0L) {
    message("[taxa_evolution] no bamdam genus reads found - skipping.")
    return(invisible(NULL))
  }

  # ---- raw reads from metrics (optional; needed only for normalized plot) ----
  raw_reads_tbl <- NULL

  # metrics may be needed for normalized plot; try to recover if path not provided
  metrics_path_final <- metrics_path
  if (is.null(metrics_path_final) || !file.exists(metrics_path_final)) {
    cand <- c(Sys.getenv("METRICS_TSV", ""), Sys.getenv("PIPELINE_METRICS", ""))
    cand <- cand[nzchar(cand)]

    guess_dirs <- unique(c(outdir, dirname(outdir), dirname(dirname(outdir))))
    cand2 <- unlist(lapply(guess_dirs, function(d) {
      if (is.null(d) || !nzchar(d)) return(character(0))
      c(file.path(d, "metrics.tsv"),
        file.path(d, "pipeline_metrics.tsv"),
        file.path(d, "reads_per_step.tsv"))
    }))

    cand <- unique(c(cand, cand2))
    cand_exist <- cand[file.exists(cand)]

    if (length(cand_exist) > 0) {
      metrics_path_final <- cand_exist[1]
      message("[taxa_evolution] using metrics file: ", metrics_path_final)
    } else {
      metrics_path_final <- NULL
      message("[taxa_evolution] metrics file not found; normalized plots will be skipped.")
    }
  }

  if (!is.null(metrics_path_final) && file.exists(metrics_path_final)) {
    m <- suppressMessages(read_tsv(metrics_path_final, show_col_types = FALSE))
    if ("sample" %in% names(m)) {

      # ---- wide metrics table (one row per sample; columns are steps) ----
      # Your pipeline metrics.tsv typically looks like: sample | Database_name | raw_reads | after_fastp | ...
      raw_col <- names(m)[grepl("^raw[-_ ]?reads$", names(m), ignore.case = TRUE)][1]
      if (!is.na(raw_col)) {
        m_wide <- m %>% dplyr::filter(sample %in% sample_order)

        raw_reads_tbl <- m_wide %>%
          dplyr::transmute(
            sample    = sample,
            raw_reads = suppressWarnings(as.numeric(.data[[raw_col]]))
          )

      } else {
        # ---- long metrics table (sample/step/reads) ----
        if (nzchar(map_tag) && "Database_name" %in% names(m)) {
          m_tag <- m %>% dplyr::filter(.data$Database_name == map_tag)
          # If the tag filter removes everything (common when Database_name is NA for some samples), don't filter.
          if (nrow(m_tag) > 0) m <- m_tag
        }

        m <- m %>% dplyr::filter(sample %in% sample_order)

        step_col  <- intersect(c("step", "Steps", "pipeline_step", "Step", "STEPS"), names(m))[1]
        reads_col <- intersect(c("reads", "read_count", "n_reads", "Reads", "READS"), names(m))[1]

        if (!is.na(step_col) && !is.na(reads_col)) {
          raw_reads_tbl <- m %>%
            dplyr::mutate(
              step  = as.character(.data[[step_col]]),
              reads = suppressWarnings(as.numeric(.data[[reads_col]]))
            ) %>%
            dplyr::filter(tolower(step) %in% c("raw-reads", "raw_reads", "raw reads", "rawreads", "raw")) %>%
            dplyr::group_by(sample) %>%
            dplyr::summarise(raw_reads = max(reads, na.rm = TRUE), .groups = "drop")
        }
      }
    }
  }

  if (!is.null(raw_reads_tbl) && nrow(raw_reads_tbl) > 0) {
    raw_reads_vec <- raw_reads_tbl$raw_reads
    names(raw_reads_vec) <- raw_reads_tbl$sample
  } else {
    raw_reads_vec <- NULL
    message("[taxa_evolution] raw-reads not found in metrics; normalized plots will be skipped.")
  }

  # helper for filenames
  safe_name <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    x
  }

  # Plot for each requested taxon (one line per plot)
  for (tx in taxa) {
    tx_key <- tolower(stringr::str_trim(tx))
    # Filter to requested taxon (case-insensitive match on genus name)
    df_tx <- bam_reads %>%
      mutate(tx_key2 = tolower(stringr::str_trim(exp_genus_name))) %>%
      filter(tx_key2 == tx_key) %>%
      select(sample, bamdam_reads)

    # Ensure all samples exist (fill 0)
    df_tx <- tibble(sample = sample_order) %>%
      left_join(df_tx, by = "sample") %>%
      mutate(
        bamdam_reads = ifelse(is.na(bamdam_reads), 0, bamdam_reads),
        sample = factor(sample, levels = sample_order)
      )

    # Apply min_reads filter as in other plots (but keep zeros for the line)
    # If *all* reads are below min_reads, still plot (otherwise user will think it's missing)
    max_reads <- max(df_tx$bamdam_reads, na.rm = TRUE)
    if (is.finite(max_reads) && max_reads < min_reads) {
      message("[taxa_evolution] NOTE: '", tx, "' max reads (", max_reads, ") < min_reads (", min_reads, "); plotting anyway.")
    }

    # Absolute plot
    p_abs <- ggplot(df_tx, aes(x = sample, y = bamdam_reads, group = 1)) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.5) +
      theme_bw(base_size = 10) +
      labs(
        x = NULL,
        y = "Bamdam reads (genus)",
        title = NULL
      ) +
      coord_cartesian(clip = "off") +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin  = grid::unit(c(top_margin_lines, 0.5, bottom_margin_lines, 1.0), "lines")
      ) +
      annotate("text",
               x = -Inf, y = Inf,
               label = "Samples",
               hjust = 1.1, vjust = -0.05, size = 2.5) +
      annotate("text",
               x = -Inf, y = -Inf,
               label = "Age",
               hjust = 1.4, vjust = 1.25, size = 2.5)

    # Top labels: same as bamdam/mmseqs plots
    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p_abs <- p_abs +
        ggnewscale::new_scale_fill() +
        geom_label(
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, fill = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size,
          label.size = 0.15,
          label.r = grid::unit(0.08, "lines"),
          key_glyph = ggplot2::draw_key_rect
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = -Inf, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust, size = sample_age_size, colour = "black",
          show.legend = FALSE
        ) +
        scale_fill_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 3, title.position = "top")
        )
    } else {
      p_abs <- p_abs +
        geom_text(
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, colour = sample_type),
          vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size, show.legend = TRUE
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = -Inf, label = label),
          vjust = age_vjust, angle = age_angle, hjust = age_hjust, size = sample_age_size, colour = "black",
          show.legend = FALSE
        ) +
        scale_colour_manual(
          name = "Sample type",
          values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
          guide = guide_legend(order = 2, title.position = "top")
        )
    }

    out_base <- file.path(outdir, paste0("taxa_evolution_", safe_name(tx)))
    pdf_abs  <- paste0(out_base, ".pdf")
    png_abs  <- paste0(out_base, ".png")

    message("[taxa_evolution] saving: ", pdf_abs)
    ggsave(pdf_abs, p_abs, width = plot_width, height = 4.8)

    message("[taxa_evolution] saving: ", png_abs)
    ggsave(png_abs, p_abs, width = plot_width, height = 4.8, dpi = 300)

    # Normalized plot (bamdam / raw reads)
    if (!is.null(raw_reads_vec)) {
      df_norm <- df_tx %>%
        mutate(
          raw_reads = as.numeric(raw_reads_vec[as.character(sample)]),
          frac_pct  = ifelse(is.na(raw_reads) | raw_reads <= 0, NA_real_, (bamdam_reads / raw_reads) * 1000000)
        )

      p_norm <- ggplot(df_norm, aes(x = sample, y = frac_pct, group = 1)) +
        geom_line(linewidth = 0.6) +
        geom_point(size = 1.5) +
        theme_bw(base_size = 10) +
        labs(
          x = NULL,
          y = "Bamdam reads / raw reads x 1 000 000",
          title = NULL
        ) +
        coord_cartesian(clip = "off") +
        theme(
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin  = grid::unit(c(top_margin_lines, 0.5, bottom_margin_lines, 1.0), "lines")
        ) +
        annotate("text",
               x = -Inf, y = Inf,
               label = "Samples",
               hjust = 1.1, vjust = -0.05, size = 2.5) +
        annotate("text",
               x = -Inf, y = -Inf,
               label = "Age",
               hjust = 1.4, vjust = 1.25, size = 2.5)

      # Top labels: same as bamdam/mmseqs plots
      if (requireNamespace("ggnewscale", quietly = TRUE)) {
        p_norm <- p_norm +
          ggnewscale::new_scale_fill() +
          geom_label(
            data = label_df_name, inherit.aes = FALSE,
            aes(x = sample, y = Inf, label = label, fill = sample_type),
            vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size,
            label.size = 0.15,
            label.r = grid::unit(0.08, "lines"),
            key_glyph = ggplot2::draw_key_rect
          ) +
          geom_text(
            data = label_df_age, inherit.aes = FALSE,
            aes(x = sample, y = -Inf, label = label),
            vjust = age_vjust, angle = age_angle, hjust = age_hjust, size = sample_age_size, colour = "black",
            show.legend = FALSE
          ) +
          scale_fill_manual(
            name = "Sample type",
            values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
            guide = guide_legend(order = 3, title.position = "top")
          )
      } else {
        p_norm <- p_norm +
          geom_text(
            data = label_df_name, inherit.aes = FALSE,
            aes(x = sample, y = Inf, label = label, colour = sample_type),
            vjust = name_vjust, angle = sample_angle, hjust = sample_hjust, size = sample_name_size, show.legend = TRUE
          ) +
          geom_text(
            data = label_df_age, inherit.aes = FALSE,
            aes(x = sample, y = -Inf, label = label),
            vjust = age_vjust, angle = age_angle, hjust = age_hjust, size = sample_age_size, colour = "black",
            show.legend = FALSE
          ) +
          scale_colour_manual(
            name = "Sample type",
            values = setNames(sample_type_cols, levels(label_df_name$sample_type)),
            guide = guide_legend(order = 2, title.position = "top")
          )
      }

      pdf_norm <- paste0(out_base, "_normalized.pdf")
      png_norm <- paste0(out_base, "_normalized.png")

      message("[taxa_evolution] saving: ", pdf_norm)
      ggsave(pdf_norm, p_norm, width = plot_width, height = 4.8)

      message("[taxa_evolution] saving: ", png_norm)
      ggsave(png_norm, p_norm, width = plot_width, height = 4.8, dpi = 300)
    }
  }

  message("[taxa_evolution] done.")
}

# ================================================================
# RUN
# ================================================================
make_metrics_plot(metrics_path, samples_path, outdir)
make_bamdam_abundance_plots(bamdam_dir, samples_path, metadata_path, outdir, min_reads, bamdam_plot_mode, taxa_per_plot)
make_mmseqs_evaluation_bubbleplot(mmseqs_dir, samples_path, metadata_path, outdir, taxa_per_plot)
make_taxa_evolution_plots(bamdam_dir, samples_path, metadata_path, metrics_path, outdir, taxa_trend_file, min_reads)

message("[100_Plots.R] all done.")
