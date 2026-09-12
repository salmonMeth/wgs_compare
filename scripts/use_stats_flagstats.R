library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

#change the directory accordingly 
work_dir <- "/scratch/project_2019524/pinksalmon_alignments"

flagstat_files <- list.files(
  work_dir,
  pattern = "\\.flagstat\\.txt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

stats_files <- list.files(
  work_dir,
  pattern = "\\.stats\\.txt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("Found", length(flagstat_files), "flagstat files\n")
cat("Found", length(stats_files), "stats files\n")

# get sample and method name from the file path

get_sample <- function(x) {
  sub(
    "_align$",
    "",
    basename(dirname(x))
  )
}

get_method <- function(x) {
  basename(dirname(dirname(x)))
}

# parsing helpers

parse_flagstat <- function(file) {
  
  x <- readLines(file, warn = FALSE)
  
  get_count <- function(pattern) {
    line <- x[grepl(pattern, x)]
    
    if (length(line) == 0)
      return(NA_real_)
    
    as.numeric(
      str_extract(line[1], "^[0-9]+")
    )
  }
  
  total <- get_count("\\+ .*in total")
  primary <- get_count("\\+ .*primary$")
  secondary <- get_count("\\+ .*secondary$")
  supplementary <- get_count("\\+ .*supplementary$")
  mapped <- get_count("\\+ .*mapped \\(")
  
  different_chr_line <- x[
    grepl(
      "with mate mapped to a different chr",
      x,
      ignore.case = TRUE
    )
  ]
  
  different_chr <- if (length(different_chr_line) > 0) {
    as.numeric(
      str_extract(
        different_chr_line[1],
        "^[0-9]+"
      )
    )
  } else {
    NA_real_
  }
  
  data.frame(
    Sample = get_sample(file),
    Method = get_method(file),
    Total_reads = total,
    Primary_reads = primary,
    Secondary_alignments = secondary,
    Supplementary_alignments = supplementary,
    Mapped_reads = mapped,
    Different_chr_mate = different_chr
  )
}

#
parse_stats <- function(file) {
  
  x <- readLines(file, warn = FALSE)
  
  sn <- x[grepl("^SN\\t", x)]
  
  get_sn <- function(pattern) {
    
    line <- sn[
      grepl(
        pattern,
        sn,
        ignore.case = TRUE
      )
    ]
    
    if (length(line) == 0)
      return(NA_real_)
    
    parts <- strsplit(
      line[1],
      "\t",
      fixed = TRUE
    )[[1]]
    
    if (length(parts) < 3)
      return(NA_real_)
    
    as.numeric(parts[3])
  }
  
  mq0 <- get_sn("^SN\\treads MQ0:")
  error_rate <- get_sn("^SN\\terror rate:")
  
  properly_paired_pct <- get_sn(
    "^SN\\tpercentage of properly paired reads"
  )
  
  properly_paired_reads <- get_sn(
    "^SN\\treads properly paired:"
  )
  
  # mapq
  
  mapq_lines <- x[
    grepl("^MAPQ\\t", x)
  ]
  
  mapq <- lapply(
    mapq_lines,
    function(line) {
      
      parts <- strsplit(
        line,
        "\t",
        fixed = TRUE
      )[[1]]
      
      tibble(
        MAPQ = as.numeric(parts[2]),
        MAPQ_count = as.numeric(parts[3])
      )
    }
  ) %>%
    bind_rows()
  
  if (nrow(mapq) > 0) {
    
    mapq_total <- sum(
      mapq$MAPQ_count
    )
    
    mapq <- mapq %>%
      mutate(
        MAPQ_percent =
          MAPQ_count /
          mapq_total *
          100
      )
    
    mapq_sorted <- mapq %>%
      arrange(MAPQ) %>%
      mutate(
        cumulative =
          cumsum(MAPQ_count) /
          mapq_total
      )
    
    median_mapq <- mapq_sorted$MAPQ[
      which(
        mapq_sorted$cumulative >= 0.5
      )[1]
    ]
    
    mean_mapq <- sum(
      mapq$MAPQ *
        mapq$MAPQ_count
    ) / mapq_total
    
  } else {
    
    mapq_total <- NA_real_
    median_mapq <- NA_real_
    mean_mapq <- NA_real_
  }
  average_length <- get_sn("^SN\\taverage length:")
  
  average_quality <- get_sn("^SN\\taverage quality:")
  
  insert_size_average <- get_sn("^SN\\tinsert size average:")
  
  insert_size_sd <- get_sn("^SN\\tinsert size standard deviation:")
  
  data.frame(
    Sample = get_sample(file),
    Method = get_method(file),
    Average_length = average_length,
    Average_quality = average_quality,
    Insert_size_average = insert_size_average,
    Insert_size_SD = insert_size_sd,
    Properly_paired_reads = properly_paired_reads,
    Properly_paired_percent = properly_paired_pct,
    MQ0_reads = mq0,
    Error_rate = error_rate,
    MAPQ_total = mapq_total,
    Mean_MAPQ = mean_mapq,
    Median_MAPQ = median_mapq
  )
  
  
  
}

# parse the stats and flagstat files
flagstat_data <- bind_rows(
  lapply(
    flagstat_files,
    parse_flagstat
  )
)

stats_data <- bind_rows(
  lapply(
    stats_files,
    parse_stats
  )
)

# put the flagstat and stats info together

alignment_summary <- full_join(
  flagstat_data,
  stats_data,
  by = c("Sample", "Method")
)

# get percents

alignment_summary <- alignment_summary %>%
  mutate(
    Mapped_percent =
      Mapped_reads /
      Total_reads *
      100,
    
    Unmapped_reads =
      Total_reads -
      Mapped_reads,
    
    Unmapped_percent =
      Unmapped_reads /
      Total_reads *
      100,
    
    Primary_percent =
      Primary_reads /
      Total_reads *
      100,
    
    Secondary_percent =
      Secondary_alignments /
      Total_reads *
      100,
    
    Supplementary_percent =
      Supplementary_alignments /
      Total_reads *
      100,
    
    MQ0_percent =
      MQ0_reads /
      Mapped_reads *
      100,
    
    Different_chr_percent =
      Different_chr_mate /
      Total_reads *
      100
  )

# save

output_dir <- file.path(
  work_dir,
  "alignment_comparison"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

write_csv(
  alignment_summary,
  file.path(
    output_dir,
    "alignment_summary.csv"
  )
)


### create the shortened file to make a latex table

input_file <- file.path(
  output_dir,
  "alignment_summary.csv"
)
# get the mapping info
report_alignment <- read_csv(
  input_file,
  show_col_types = FALSE
) %>%
  select(
    Sample,
    Method,
    Total_reads,
    Primary_percent,
    Properly_paired_percent,
    Unmapped_percent
  ) %>%
  mutate(
    Total_reads = round(Total_reads),
    Primary_percent = round(Primary_percent, 3),
    Properly_paired_percent = round(Properly_paired_percent, 3),
    Unmapped_percent = round(Unmapped_percent, 3)
  )

output_file <- file.path(
  output_dir,
  "alignment_summary_for_latex.csv"
)

write_csv(
  report_alignment,
  output_file
)

#get the read length + insert length info

# get the mapping info
report_insert <- read_csv(
  input_file,
  show_col_types = FALSE
) %>%
  select(
    Sample,
    Method,
    Average_length,
    Insert_size_average,
    Insert_size_SD,
  ) 

output_file <- file.path(
  output_dir,
  "insert_read_lengths_for_latex.csv"
)

write_csv(
  report_insert,
  output_file
)

# ------------------------------------------------------------
# Create MAPQ distribution data
# ------------------------------------------------------------

mapq_data <- bind_rows(
  lapply(
    stats_files,
    function(file) {
      
      x <- readLines(
        file,
        warn = FALSE
      )
      
      lines <- x[
        grepl("^MAPQ\\t", x)
      ]
      
      if (length(lines) == 0)
        return(NULL)
      
      out <- lapply(
        lines,
        function(line) {
          
          parts <- strsplit(
            line,
            "\t",
            fixed = TRUE
          )[[1]]
          
          tibble(
            Sample =
              get_sample(file),
            Method =
              get_method(file),
            MAPQ =
              as.numeric(parts[2]),
            Count =
              as.numeric(parts[3])
          )
        }
      ) %>%
        bind_rows()
      
      total <- sum(out$Count)
      
      out %>%
        mutate(
          Percent =
            Count / total * 100
        )
    }
  )
)

write_csv(
  mapq_data,
  file.path(
    output_dir,
    "MAPQ_distribution.csv"
  )
)

cat(
  "\nSaved files to:\n",
  output_dir,
  "\n"
)
