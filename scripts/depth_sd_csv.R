library(tidyverse)

#change to the path of the files created using depth_sd.slurm
input_dir <- "/scratch/project_2019524/depth_evenness_stats"

output_file <- file.path(
  input_dir,
  "genome_wide_depth_evenness_stats.csv"
)

#change the sample IDs accordingly
desired_samples <- c(
  "8034",
  "8035",
  "8036",
  "8038"
)

files <- list.files(
  input_dir,
  pattern = "_depth_evenness\\.txt$",
  full.names = TRUE
)

genome_wide <- map_dfr(
  files,
  function(file) {
    
    read_tsv(
      file,
      col_types = cols(
        Sample = col_character(),
        Method = col_character(),
        Level = col_character(),
        chr_number = col_character(),
        chromosome = col_character(),
        Mean_depth = col_character(),
        SD_depth = col_character(),
        CV_depth = col_character(),
        Min_depth = col_character(),
        Max_depth = col_character()
      ),
      show_col_types = FALSE
    ) %>%
      filter(
        Level == "Genome-wide",
        chromosome == "Genome-wide"
      )
  }
)

# this fixes the sample ID issue with 8167 <-> 8036 and 8168 <->  8038
genome_wide$Sample[
  genome_wide$Sample == "8167"
] <- "8036"

genome_wide$Sample[
  genome_wide$Sample == "8168"
] <- "8038"

genome_wide <- genome_wide %>%
  filter(
    Sample %in% desired_samples
  ) %>%
  select(
    Sample,
    Method,
    Level,
    Mean_depth,
    SD_depth,
    CV_depth,
    Min_depth,
    Max_depth
  ) %>%
  arrange(
    Sample,
    Method
  )

write_csv(
  genome_wide,
  output_file
)

message("Written: ", output_file)
