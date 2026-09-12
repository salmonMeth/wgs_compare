library(ggplot2)
library(dplyr)
library(zoo)
library(ggrepel)

#plots the deduplicated sequence counts and the quality score boxplots

################################################
################################################
# plot the deduplicated sequence counts
library(jsonlite)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

methods <- c(
  "Low input" = "/scratch/project_2019524/PinkSalmon_lowinput",
  "Low input P100"  = "/scratch/project_2019524/PinkSalmon_lowinput_P100",
  "Normal input"  = "/scratch/project_2019524/PinkSalmon_normalinput"
)

extract_unique_reads <- function(method_dir, method_name) {
  
  json_file <- file.path(
    method_dir,
    "qc",
    "multiqc",
    "multiqc_data",
    "multiqc_data.json"
  )
  
  if (!file.exists(json_file)) {
    stop("File not found: ", json_file)
  }
  
  message("Reading: ", json_file)
  
  mq <- fromJSON(
    json_file,
    simplifyVector = FALSE
  )
  
  # MultiQC FastQC Sequence Counts plot
  seq_counts <- mq$report_plot_data$fastqc_sequence_counts_plot
  
  dataset <- seq_counts$datasets[[1]]
  
  # Sample names
  samples <- unlist(dataset$samples)
  
  # Find "Unique Reads"
  unique_category <- dataset$cats[
    which(
      sapply(
        dataset$cats,
        function(x) x$name == "Unique Reads"
      )
    )
  ][[1]]
  
  # Unique / deduplicated read counts
  unique_reads <- unlist(unique_category$data)
  
  # Return data frame
  tibble(
    sample = samples,
    method = method_name,
    deduplicated_sequences = as.numeric(unique_reads)
  )
}


df <- bind_rows(
  lapply(
    seq_along(methods),
    function(i) {
      extract_unique_reads(
        method_dir = methods[[i]],
        method_name = names(methods)[i]
      )
    }
  )
)

df <- df %>%
  mutate(
    plot_label = str_extract(
      sample,
      "(8034|8035|8036|8167|8038|8168)_[12]$"
    )
  )

#order the samples by method
method_order <- unique(df$method)

df <- df %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )
#order the samples within the methods

sample_order <- unique(df$plot_label)

df <- df %>%
  mutate(
    plot_label = factor(
      plot_label,
      levels = sample_order
    )
  ) %>%
  arrange(method, plot_label)


df <- df %>%
  mutate(
    x_position = paste(method, plot_label, sep = "__")
  )

df$x_position <- factor(
  df$x_position,
  levels = unique(df$x_position)
)
#plot
p <- ggplot(
  df,
  aes(
    x = x_position,
    y = deduplicated_sequences,
    fill = method
  )
) +
  geom_col(
    width = 0.8
  ) +
  scale_x_discrete(
    labels = function(x) {
      sub(
        ".*__",
        "",
        x
      )
    }
  ) +
  scale_y_continuous(
    labels = label_number(
      scale = 1e-6,
      suffix = "M"
    ),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    title = "Deduplicated Sequence Counts",
    x = "Sample ID",
    y = "Deduplicated sequences (millions)",
    fill = "Sequencing method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )
#save
ggsave(
  filename = "/scratch/project_2019524/plots/fastq_multiq/deduplicated_sequence_counts_by_method.png",
  plot = p,
  width = 14,
  height = 7,
  dpi = 300
)

##################################################################
##################################################################
#plot the quality scores
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

files <- c(
  "/scratch/project_2019524/PinkSalmon_lowinput/qc/multiqc/multiqc_data/multiqc_data.json",
  "/scratch/project_2019524/PinkSalmon_lowinput_P100/qc/multiqc/multiqc_data/multiqc_data.json",
  "/scratch/project_2019524/PinkSalmon_normalinput/qc/multiqc/multiqc_data/multiqc_data.json"
)

methods <- c(
  "Low input",
  "Low input P100",
  "Normal input"
)

extract_quality <- function(file, method_name) {
  
  mqc <- fromJSON(
    file,
    simplifyVector = FALSE
  )
  
  quality <- mqc$report_plot_data$fastqc_per_sequence_quality_scores_plot
  
  map_dfr(
    quality$datasets[[1]]$lines,
    function(line) {
      
      sample_name <- line$name
      
      map_dfr(
        line$pairs,
        function(pair) {
          
          data.frame(
            sample = sample_name,
            phred = pair[[1]],
            count = pair[[2]]
          )
        }
      )
    }
  ) %>%
    mutate(
      method = method_name,
      
      sample_label = str_extract(
        sample,
        "(8034|8035|8036|8167|8038|8168)_[12]$"
      )
    )
}

quality_df <- map2_dfr(files,methods,extract_quality)


weighted_quantile <- function(x, w, probs) {
  
  ord <- order(x)
  
  x <- x[ord]
  w <- w[ord]
  
  cumulative <- cumsum(w)
  total <- sum(w)
  
  sapply(
    probs,
    function(p) {
      x[which(cumulative >= p * total)[1]]
    }
  )
}

quality_box <- quality_df %>%
  group_by(
    method,
    sample_label
  ) %>%
  summarise(
    ymin = weighted_quantile(phred, count, 0),
    lower = weighted_quantile(phred, count, 0.25),
    middle = weighted_quantile(phred, count, 0.50),
    upper = weighted_quantile(phred, count, 0.75),
    ymax = weighted_quantile(phred, count, 1),
    .groups = "drop"
  )

quality_box <- quality_box %>%
  mutate(
    # Extract just 8034, 8035, etc. from 8034_1
    sample_id = str_extract(
      sample_label,
      "8034|8035|8036|8038|8167|8168"
    ),
    # Make sure methods stay in the desired order
    method = factor(
      method,
      levels = methods
    )
  ) %>%
  arrange(method, sample_id, sample_label) %>%
  mutate(
    # Unique position for every method + sample
    x = factor(
      paste(method, sample_label, sep = "__"),
      levels = paste(method, sample_label, sep = "__")
    )
  )

#define the colors
blue_shades <- c(
  "8034" = "#C6DBEF",
  "8035" = "#9ECAE1",
  "8036" = "#6BAED6",
  "8038" = "#2171B5"
)

orange_shades <- c(
  "8034" = "#FDD49E",
  "8035" = "#FDBB84",
  "8036" = "#FC8D59",
  "8038" = "#D7301F"
)

green_shades <- c(
  "8167" = "#C7E9C0",
  "8168" = "#238B45"
)

#assign the colors
quality_box <- quality_box %>%
  mutate(
    method_number = as.integer(method),
    
    fill_colour = case_when(
      
      method_number == 1 ~
        unname(blue_shades[as.character(sample_id)]),
      
      method_number == 2 ~
        unname(orange_shades[as.character(sample_id)]),
      
      method_number == 3 ~
        unname(green_shades[as.character(sample_id)]),
      
      TRUE ~ "grey70"
    )
  )
#legend
method_colours <- c(
  "Low input" = "#377EB8",
  "Low input P100" = "#E69F00",
  "Normal input" = "#009E73"
)
#plot
p_quality <- ggplot(
  quality_box,
  aes(
    x = x,
    ymin = ymin,
    lower = lower,
    middle = middle,
    upper = upper,
    ymax = ymax
  )
) +

  geom_boxplot(
    stat = "identity",
    aes(fill = fill_colour),
    width = 0.7,
    colour = "grey25"
  ) +

  geom_point(
    aes(
      x = x,
      y = middle,
      colour = method
    ),
    alpha = 0
  ) +
  
  scale_fill_identity() +
  
  scale_colour_manual(
    name = "Sequencing method",
    values = c(
      "Low input" = "#377EB8",
      "Low input P100" = "#E69F00",
      "Normal input" = "#009E73"
    )
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        alpha = 1,
        size = 5,
        shape = 15
      )
    )
  ) +
  
  scale_x_discrete(
    labels = quality_box$sample_label
  ) +
  
  labs(
    title = "Per-Sequence Quality Scores",
    x = "Sample",
    y = "Mean sequence quality (Phred)"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave(filename = "/scratch/project_2019524/plots/fastq_multiq/quality_score_boxplots_by_method.png",
       plot = p_quality,
       width = 14,
       height = 7,
       dpi = 300
)