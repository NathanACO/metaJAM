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

message("[100_Plots.R] damage_threshold(%): ", damage_threshold)

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

  p <- ggplot(df, aes(x = step, y = reads, colour = sample)) +
    geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0.15, height = 0)) +
    scale_y_log10() +
    labs(x = NULL, y = "Reads (log10)", colour = "Sample") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1),
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
make_bamdam_abundance_plots <- function(bamdam_dir, samples_path, metadata_path, outdir, min_reads = 1, plot_mode = "heatmap") {
  message("[bamdam] starting bamdam plot(s)")

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
  message("[bamdam] found ", length(files), " bamdam *.tsv files")
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

    folder_name <- basename(dirname(f))
    file_base   <- basename(f)
    file_base   <- sub("\\.tsv$", "", file_base)

    sample_id <- if (!is.na(folder_name) && nzchar(folder_name) && folder_name != "bamdam") {
      folder_name
    } else {
      file_base
    }

    ranks <- vapply(x$taxpath, get_rank, character(1))

    tibble(
      sample  = sample_id,
      taxon   = as.character(x$TaxName),
      rank    = ranks,
      reads   = as.numeric(x$TotalReads),
      damage_p1 = as.numeric(x$`Damage+1`),
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
  if (nrow(meta_filt) == 0L) meta_filt <- meta

  meta_filt <- meta_filt %>% arrange(.data[[age_col]])
  sample_order <- meta_filt$sample

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

  heat_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(256)

  ages <- meta_filt %>% select(sample, !!age_col)
  colnames(ages)[2] <- "age"
  age_vec <- ages$age; names(age_vec) <- ages$sample

  label_df_name <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = sample_order,
    sample_type = sample_types
  )
  label_df_age <- tibble(
    sample = factor(sample_order, levels = sample_order),
    label  = as.character(age_vec[sample_order])
  )

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
    dmg   = dplyr::first(damage_p1),
    .groups = "drop"
  ) %>%
  mutate(dmg_pct = 100 * dmg)

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

    list(df = df_long, max_log = max(df_long$log_reads, na.rm = TRUE))
  }

  # ---- plot: heatmap ----
  plot_heatmap <- function(df_long, max_log, out_prefix) {
    max_log <- max(max_log, max(legend_breaks))
    min_log <- 0

    df_long <- df_long %>%
      mutate(label = formatC(as.integer(round(reads)), format = "f", digits = 0, big.mark = ","))

    p <- ggplot(df_long, aes(x = sample, y = taxon, fill = log_reads)) +
      geom_tile() +
      geom_text(aes(label = label), size = 2) +
      scale_x_discrete(
        position = "top",
        limits   = sample_order,
        labels   = rep("", length(sample_order)),
        expand   = expansion(add = c(0.5, 0))
      ) +
      scale_y_discrete(name = NULL) +
      scale_fill_gradientn(
        colours = heat_colors,
        limits  = c(min_log, max_log),
        breaks  = legend_breaks,
        labels  = formatC(legend_counts, format = "fg", big.mark = ","),
        name    = "Unique Reads"
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
        legend.text         = element_text(size = 7),
        plot.margin         = grid::unit(c(1.2, 0.5, 0.5, 1.0), "lines")
      ) +
      annotate("text",
               x = 0.5, y = Inf,
               label = "Samples\nThousand years ago",
               hjust = 1, vjust = -0.3, size = 2.5)

    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_label(
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, fill = sample_type),
          vjust = -1.3, size = 2.7,
          label.size = 0.15,
          label.r = grid::unit(0.08, "lines"),
          key_glyph = ggplot2::draw_key_rect
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label),
          vjust = -0.7, size = 2.6, colour = "black",
          show.legend = FALSE
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
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, colour = sample_type),
          vjust = 1.35, size = 2.7, show.legend = TRUE
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label),
          vjust = 2.60, size = 2.6, colour = "black",
          show.legend = FALSE
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
    ggsave(pdf_file, p, width = 9, height = 8)

    message("[bamdam] saving: ", png_file)
    ggsave(png_file, p, width = 9, height = 8, dpi = 300)
  }

  # ---- plot: bubble ----
  plot_bubble_damage <- function(df_long, max_log, out_prefix) {
    max_log <- max(max_log, max(legend_breaks))
    min_log <- 0

    # bubble size based on log10(reads+1)
    df_plot <- df_long %>%
      dplyr::filter(reads > 0) %>%
      dplyr::mutate(label = formatC(as.integer(round(reads)), format = "f", digits = 0, big.mark = ","))

    p <- ggplot(df_plot, aes(x = sample, y = taxon)) +
      geom_point(aes(size = log_reads, fill = damage_class), shape = 21, colour = "black", stroke = 0.15, alpha = 0.9) +
      geom_text(aes(label = label), colour = "black", size = 2) +
      scale_x_discrete(
        position = "top",
        limits   = sample_order,
        labels   = rep("", length(sample_order)),
        expand   = expansion(add = c(0.6, 0.6))
      ) +
      scale_y_discrete(name = NULL) +
      scale_fill_manual(
        name   = "Damage",
        values = c(red = "red", orange = "orange", green = "green"),
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
      theme_bw(base_size = 9) +
      theme(
        panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.title.x.bottom = element_blank(),
        axis.ticks.x        = element_blank(),
        panel.grid          = element_blank(),
        legend.position     = "right",
        legend.title        = element_text(size = 8),
        legend.text         = element_text(size = 7),
        plot.margin         = grid::unit(c(1.2, 0.5, 0.5, 1.0), "lines")
      ) +
      annotate("text",
               x = 0.5, y = Inf,
               label = "Samples\nThousand years ago",
               hjust = 1, vjust = -0.3, size = 2.5)

    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        geom_label(
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, fill = sample_type),
          vjust = -1.3, size = 2.7,
          label.size = 0.15,
          label.r = grid::unit(0.08, "lines"),
          key_glyph = ggplot2::draw_key_rect
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label),
          vjust = -0.7, size = 2.6, colour = "black",
          show.legend = FALSE
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
          data = label_df_name, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label, colour = sample_type),
          vjust = 1.35, size = 2.7, show.legend = TRUE
        ) +
        geom_text(
          data = label_df_age, inherit.aes = FALSE,
          aes(x = sample, y = Inf, label = label),
          vjust = 2.60, size = 2.6, colour = "black",
          show.legend = FALSE
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
    ggsave(pdf_file, p, width = 9, height = 8)

    message("[bamdam] saving: ", png_file)
    ggsave(png_file, p, width = 9, height = 8, dpi = 300)
  }

plot_bubble_reads <- function(df_long, max_log, out_prefix) {
  legend_counts <- c(1, 10, 100, 1000, 10000, 100000)
  legend_breaks <- log10(legend_counts)

  max_log <- max(max_log, max(legend_breaks))
  min_log <- 0

  df_plot <- df_long %>%
    dplyr::filter(reads > 0) %>%
    dplyr::mutate(label = formatC(as.integer(round(reads)), format = "f", digits = 0, big.mark = ","))

  p <- ggplot(df_plot, aes(x = sample, y = taxon)) +
    geom_point(aes(size = log_reads, fill = log_reads),
               shape = 21, colour = "black", stroke = 0.15, alpha = 0.9) +
    geom_text(aes(label = label), colour = "black", size = 2) +
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
      name  = "Unique Reads",
      range = c(0.5, 6),
      guide = guide_legend(order = 2, title.position = "top")
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 9) +
    theme(
      axis.title.x.bottom = element_blank(),
      axis.ticks.x        = element_blank(),
      panel.grid          = element_blank(),
      legend.position     = "right",
      legend.title        = element_text(size = 8),
      legend.text         = element_text(size = 7),
      plot.margin         = grid::unit(c(1.2, 0.5, 0.5, 1.0), "lines")
    ) +
    annotate("text",
             x = 0.5, y = Inf,
             label = "Samples\nThousand years ago",
             hjust = 1, vjust = -0.3, size = 2.5)

  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_label(
        data = label_df_name, inherit.aes = FALSE,
        aes(x = sample, y = Inf, label = label, fill = sample_type),
        vjust = -1.3, size = 2.7,
        linewidth = 0.15,
        label.r = grid::unit(0.08, "lines"),
        key_glyph = ggplot2::draw_key_rect
      ) +
      geom_text(
        data = label_df_age, inherit.aes = FALSE,
        aes(x = sample, y = Inf, label = label),
        vjust = -0.7, size = 2.6, colour = "black",
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
        aes(x = sample, y = Inf, label = label),
        vjust = 1.35, size = 2.7, colour = "black", show.legend = FALSE
      ) +
      geom_text(
        data = label_df_age, inherit.aes = FALSE,
        aes(x = sample, y = Inf, label = label),
        vjust = 2.60, size = 2.6, colour = "black",
        show.legend = FALSE
      )
  }

  pdf_file <- file.path(outdir, paste0(out_prefix, ".pdf"))
  png_file <- file.path(outdir, paste0(out_prefix, ".png"))

  message("[bamdam] saving: ", pdf_file)
  ggsave(pdf_file, p, width = 9, height = 8)

  message("[bamdam] saving: ", png_file)
  ggsave(png_file, p, width = 9, height = 8, dpi = 300)
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

    if (plot_mode %in% c("heatmap", "both")) {
      out_prefix <- paste0("bamdam_", rank_sel, "_heatmap")
      plot_heatmap(df_long, max_log, out_prefix)
    }
    if (plot_mode %in% c("bubble", "both")) {
      out_prefix <- paste0("bamdam_", rank_sel, "_bubbleplot")
      plot_bubble_damage(df_long, max_log, paste0(out_prefix, "_damage"))
      plot_bubble_reads(df_long, max_log, paste0(out_prefix, "_reads"))
    }
  }

  message("[bamdam] done.")
}

# ================================================================
# RUN
# ================================================================
make_metrics_plot(metrics_path, samples_path, outdir)
make_bamdam_abundance_plots(bamdam_dir, samples_path, metadata_path, outdir, min_reads, bamdam_plot_mode)

message("[100_Plots.R] all done.")
