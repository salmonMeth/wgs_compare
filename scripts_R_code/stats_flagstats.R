### the code to extract some stats from the stats.txt and flagstats.txt files
### and save them as csv files

library(tidyverse)

## change the paths and sample IDs accordingly

BAMDIR <- "/scratch/project_2019524/pinksalmon_alignments/lowinput"

SAMPLES <- c("8034", "8035", "8036", "8038")

##

OUTDIR <- file.path(BAMDIR, "flagstat_QC")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# get infor from the flagstat file

parse_flagstat <- function(sample) {
  
  file <- file.path(
    BAMDIR,
    paste0(sample, "_align"),
    paste0(sample, ".flagstat.txt")
  )
  
  if (!file.exists(file)) {
    stop("Flagstat file not found: ", file)
  }
  
  x <- readLines(file)
  
  # helper
  
  extract_number <- function(pattern) {
    
    line <- x[grepl(pattern, x)]
    
    if (length(line) == 0) {
      return(NA_real_)
    }
    
    as.numeric(
      sub("^([0-9]+).*", "\\1", line[1])
    )
  }
  
  ####################
  # extracting the stats we are interested in
  #total number of reads, mapped reads, properly paired reads
  # mate mapped to different chromosome
  #
  total_reads <- extract_number(
    "in total"
  )
  mapped_line <- x[
    grepl(" mapped \\(", x) &
      !grepl("primary mapped", x)
  ]
  
  if (length(mapped_line) != 1) {
    stop(
      "Could not uniquely identify mapped line for sample ",
      sample,
      "\nFound:\n",
      paste(mapped_line, collapse = "\n")
    )
  }
  
  mapped_reads <- as.numeric(
    sub("^([0-9]+).*", "\\1", mapped_line)
  )
  
  mapped_percent <- as.numeric(
    sub(
      ".*\\(([0-9.]+)%.*",
      "\\1",
      mapped_line
    )
  )


  #
  paired_line <- x[grepl("properly paired", x)]
  
  properly_paired_reads <- as.numeric(
    sub("^([0-9]+).*", "\\1", paired_line[1])
  )
  
  properly_paired_percent <- as.numeric(
    sub(".*\\(([0-9.]+)%.*","\\1",paired_line
    )
  )

  #
  diff_chr_line <- x[
    grepl(
      "with mate mapped to a different chr$",
      x
    )
  ]
  
  mate_diff_chr <- as.numeric(
    sub("^([0-9]+).*", "\\1", diff_chr_line[1])
  )
  
  
  data.frame(
    sample = sample,
    total_reads = total_reads,
    mapped_reads = mapped_reads,
    mapped_percent = mapped_percent,
    properly_paired_reads = properly_paired_reads,
    properly_paired_percent = properly_paired_percent,
    mate_diff_chr = mate_diff_chr
  )
}

# go through all the samples

flagstat <- map_dfr(
  SAMPLES,
  parse_flagstat
)


# add some percentages

flagstat <- flagstat %>%
  mutate(
    mate_diff_chr_percent =
      mate_diff_chr / total_reads * 100
  )


write_csv(
  flagstat,
  file.path(OUTDIR, "flagstat_summary.csv")
)

##########get some stats from the stats.txt file
OUTDIR <- file.path(BAMDIR, "stats_QC")

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)

#various heper functions
read_stats <- function(sample) {
  file <- file.path(
    BAMDIR,
    paste0(sample, "_align"),
    paste0(sample, ".stats.txt")
  )
  if (!file.exists(file)) {
    stop("Stats file not found: ", file)
  }
  readLines(file)
}
get_SN <- function(lines, name) {
  line <- lines[
    grepl(
      paste0("^SN\\t", name, ":"),
      lines
    )
  ]
  
  if (length(line) == 0) {
    return(NA_character_)
  }
    value <- sub(
    paste0("^SN\\t", name, ":\\t"),
    "",
    line[1]
  )
  value
}
get_SN_numeric <- function(lines, name) {
  value <- get_SN(lines, name)
  if (is.na(value)) {
    return(NA_real_)
  }
    value <- sub(
    "\\s+#.*$",
    "",
    value
  )
    value <- trimws(value)
  as.numeric(value)
}
extract_sample_stats <- function(sample) {
  
  lines <- read_stats(sample)
  
  data.frame(
    
    sample = sample,
    
    # Error rate NOTE this is #mismatches/#mapped bases
    error_rate =
      get_SN_numeric(
        lines,
        "error rate"
      ),
    
    # Average read length
    average_length =
      get_SN_numeric(
        lines,
        "average length"
      ),
    
    # Average quality
    average_quality =
      get_SN_numeric(
        lines,
        "average quality"
      ),
    
    # Insert size average
    insert_size_average =
      get_SN_numeric(
        lines,
        "insert size average"
      ),
    
    # Insert size standard deviation
    insert_size_sd =
      get_SN_numeric(
        lines,
        "insert size standard deviation"
      )
  )
}

stats_summary <- map_dfr(
  SAMPLES,
  extract_sample_stats
)
write_csv(
  stats_summary,
  file.path(
    OUTDIR,
    "samtools_stats_summary.csv"
  )
)

