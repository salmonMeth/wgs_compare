library(tidyverse)

# Paths
input_base <- "/scratch/project_2019524/pinksalmon_alignments"
fai_file <- "/scratch/project_2019524/pinksalmon_reference/GCF_021184085.1_OgorEven_v1.0_genomic.fna.fai"

output_dir <- "/scratch/project_2019524/pinksalmon_depth_evenness_binned_approx"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c("lowinput", "lowinput_P100", "normalinput")

# Old -> desired sample IDs
sample_map <- c(
  "8167" = "8036",
  "8168" = "8038"
)

# The 27 actual chromosomes
chromosomes <- tibble(
  chromosome = paste0("NC_", sprintf("%06d", 60173:60199), ".1"),
  chr_number = 1:27
)

# Read chromosome lengths from FAI
fai <- read_tsv(
  fai_file,
  col_names = c("chromosome", "chr_length", "offset", "line_bases", "line_width"),
  show_col_types = FALSE
) %>%
  select(chromosome, chr_length)

# Keep only the 27 chromosomes
chromosomes <- chromosomes %>%
  left_join(fai, by = "chromosome")

# Check that all 27 chromosomes were found
if (any(is.na(chromosomes$chr_length))) {
  missing_chr <- chromosomes$chromosome[is.na(chromosomes$chr_length)]
  stop(
    "Missing chromosome(s) in FAI: ",
    paste(missing_chr, collapse = ", ")
  )
}

# Calculate weighted statistics for one chromosome
weighted_stats <- function(df, chr_length) {
  
  # Correct bin weight:
  # actual number of bases represented by each bin
  df <- df %>%
    mutate(
      weight = pmin(window_end, chr_length) - window_start + 1
    ) %>%
    filter(weight > 0, !is.na(mean_depth))
  
  total_weight <- sum(df$weight)
  
  # Weighted mean
  mean_depth <- sum(df$mean_depth * df$weight) / total_weight
  
  # Weighted population SD
  sd_depth <- sqrt(
    sum(df$weight * (df$mean_depth - mean_depth)^2) / total_weight
  )
  
  # CV
  cv_depth <- ifelse(
    mean_depth == 0,
    NA_real_,
    sd_depth / mean_depth
  )
  
  # Weighted empirical quantiles
  df_sorted <- df %>%
    arrange(mean_depth) %>%
    mutate(cum_weight = cumsum(weight))
  
  weighted_quantile <- function(p) {
    threshold <- p * total_weight
    df_sorted$mean_depth[
      which(df_sorted$cum_weight >= threshold)[1]
    ]
  }
  
  tibble(
    Mean_depth = mean_depth,
    SD_depth = sd_depth,
    CV_depth = cv_depth,
    Q25_depth = weighted_quantile(0.25),
    Median_depth = weighted_quantile(0.50),
    Q75_depth = weighted_quantile(0.75),
    Min_depth = min(df$mean_depth),
    Max_depth = max(df$mean_depth)
  )
}

# Process each sample and method
for (method in methods) {
  
  depth_dir <- file.path(
    input_base,
    method,
    "coverage_depth",
    "depth_windows"
  )
  
  files <- list.files(
    depth_dir,
    pattern = "\\.depth_100000bp\\.txt$",
    full.names = TRUE
  )
  
  for (file in files) {
    
    filename <- basename(file)
    
    # Extract sample ID from filename
    sample_id <- sub("\\.depth_100000bp\\.txt$", "", filename)
    
    # Rename old IDs
    if (sample_id %in% names(sample_map)) {
      sample_id <- sample_map[[sample_id]]
    }
    
    # Only process desired samples
    if (!sample_id %in% c("8034", "8035", "8036", "8038")) {
      next
    }
    
    message("Processing ", sample_id, " ", method)
      
    depth <- read_tsv( file, col_names = TRUE, show_col_types = FALSE ) %>% mutate( window_start = as.numeric(window_start), window_end = as.numeric(window_end), mean_depth = as.numeric(mean_depth) )
    
    # Keep only the 27 chromosomes and attach chromosome lengths
    depth_chr <- depth %>%
      inner_join(chromosomes, by = "chromosome")
    
    # Calculate statistics chromosome by chromosome
    results <- depth_chr %>%
      group_by(chr_number, chromosome, chr_length) %>%
      group_modify(~ weighted_stats(.x, .y$chr_length)) %>%
      ungroup() %>%
      mutate(
        Sample = sample_id,
        Method = method,
        .before = 1
      ) %>%
      select(
        Sample,
        Method,
        chr_number,
        chromosome,
        Mean_depth,
        SD_depth,
        CV_depth,
        Q25_depth,
        Median_depth,
        Q75_depth,
        Min_depth,
        Max_depth
      ) %>%
      arrange(chr_number)
    
    # Output contains only the desired sample ID
    output_file <- file.path(
      output_dir,
      paste0(sample_id, "_", method, "_depth_evenness_weighted_chr.txt")
    )
    
    write_tsv(results, output_file)
    
    message("Written: ", output_file)
  }
}
