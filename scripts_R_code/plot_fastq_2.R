library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)

root="/scratch/project_2019524"

methods=c("PinkSalmon_lowinput",
          "PinkSalmon_lowinput_P100",
          "PinkSalmon_normalinput")

files=file.path(root,methods,
                "qc/multiqc/multiqc_data/multiqc_fastqc.txt")

#get all the multiqc files and put them in a dataframe
dat <- bind_rows(lapply(seq_along(files), function(i) {
  
  x <- read.delim(
    files[i],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  x$Method <- methods[i]
  
  x
}))

dat$Read <- ifelse(
  grepl("_1", dat$Filename),
  "R1",
  ifelse(grepl("_2", dat$Filename), "R2", NA)
)

dat$Sample <- sub("_([12]).*$", "", dat$Filename)

#we get the genome size using the fasta file  but we use a bash terminal command due to file size
#execute the following from the terminal
#awk '!/^>/ {total += length($0)} END {print total}' /scratch/project_2019524/pinksalmon_reference/GCF_021184085.1_OgorEven_v1.0_genomic.fna
#result is 2690283663
genome_size=2690283663

dat$TotalBases_numeric <- as.numeric(
  sub(" .*", "", dat$`Total Bases`)
) * 1e9

summary_table <- dat %>%
  group_by(Method, Sample) %>%
  summarise(
    TotalReads = first(`Total Sequences`),
    TotalBases = sum(TotalBases_numeric),
    UniquePercent = mean(
      `total_deduplicated_percentage`,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    DuplicatePercent = 100 - UniquePercent,
    TheoreticalCoverage = TotalBases / genome_size,
    UniqueBases = TotalBases * UniquePercent / 100,
    DeduplicatedCoverage = TheoreticalCoverage * UniquePercent / 100
  )

summary_table <- summary_table %>%
  mutate(
    TotalBases_Gb = round(TotalBases / 1e9, 2),
    UniquePercent = round(UniquePercent, 1),
    DuplicatePercent = round(DuplicatePercent, 1),
    UniqueBases_Gb = round(UniqueBases / 1e9, 2),
    TheoreticalCoverage = round(TheoreticalCoverage, 2),
    DeduplicatedCoverage = round(DeduplicatedCoverage, 2)
  ) %>%
  select(
    Method,
    Sample,
    TotalReads,
    TotalBases_Gb,
    TheoreticalCoverage,
    UniquePercent,
    DuplicatePercent,
    UniqueBases_Gb,
    DeduplicatedCoverage
  )

write.csv(
  summary_table,
  "/scratch/project_2019524/plots/fastq_multiq/pink_salmon_coverage_duplication_summary.csv",
  row.names = FALSE
)

##########################
##########################
#plot the theoretical vs deduplicated theoretical coverage
base_colors <- c(
  "PinkSalmon_lowinput"      = "#0072B2",  # blue
  "PinkSalmon_lowinput_P100" = "#D55E00",  # orange/red
  "PinkSalmon_normalinput"   = "#009E73"   # green/teal
)


# create ids for the samples

sample_info <- summary_table %>%
  distinct(Method, Sample) %>%
  group_by(Method) %>%
  arrange(Sample) %>%
  mutate(
    SampleNumber = row_number(),
    SampleID = paste(Method, Sample, sep = "_")
  ) %>%
  ungroup()

# assign shades to the samples

sample_colors <- c()

for (m in names(base_colors)) {
  
  samples <- sample_info %>%
    filter(Method == m)
  
  n <- nrow(samples)
  
  # Different vibrant shades within each method
  if (m == "PinkSalmon_lowinput") {
    
    shades <- grDevices::colorRampPalette(
      c("#6BAED6", base_colors[m])
    )(n)
    
  } else if (m == "PinkSalmon_lowinput_P100") {
    
    shades <- grDevices::colorRampPalette(
      c("#FC8D59", base_colors[m])
    )(n)
    
  } else if (m == "PinkSalmon_normalinput") {
    
    shades <- grDevices::colorRampPalette(
      c("#66C2A5", base_colors[m])
    )(n)
  }
  
  names(shades) <- samples$SampleID
  
  sample_colors <- c(
    sample_colors,
    shades
  )
}


# get the coverage data ready

coverage_long <- summary_table %>%
  left_join(
    sample_info,
    by = c("Method", "Sample")
  ) %>%
  pivot_longer(
    cols = c(
      TheoreticalCoverage,
      DeduplicatedCoverage
    ),
    names_to = "CoverageType",
    values_to = "Coverage"
  ) %>%
  mutate(
    
    CoverageType = factor(
      CoverageType,
      levels = c(
        "TheoreticalCoverage",
        "DeduplicatedCoverage"
      ),
      labels = c(
        "Theoretical",
        "Deduplicated"
      )
    ),
    
    MethodGroup = case_when(
      Method == "PinkSalmon_lowinput" ~ 1,
      Method == "PinkSalmon_lowinput_P100" ~ 2,
      Method == "PinkSalmon_normalinput" ~ 3
    ),
    
    x = MethodGroup +
      ifelse(
        CoverageType == "Theoretical",
        -0.15,
        0.15
      )
  )
# plot

p <- ggplot(
  coverage_long,
  aes(
    x = x,
    y = Coverage,
    group = SampleID,
    color = SampleID
  )
) +
  
  geom_vline(
    xintercept = c(1.5, 2.5),
    linewidth = 0.7,
    color = "grey75"
  ) +
  
  geom_line(
    linewidth = 0.9,
    alpha = 0.7
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  # Sample ID labels
  geom_text_repel(
    data = subset(
      coverage_long,
      CoverageType == "Deduplicated"
    ),
    aes(
      label = substr(
        Sample,
        nchar(Sample) - 3,
        nchar(Sample)
      )
    ),
    direction = "y",
    hjust = 0,
    nudge_x = 0.05,
    size = 3,
    alpha = 0.8,
    segment.color = "grey70",
    min.segment.length = 0
  )  +
  
  scale_color_manual(
    values = sample_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = 1:3,
    labels = c(
      "Low input",
      "Low input P100",
      "Normal input"
    ),
    limits = c(0.5, 3.5)
  ) +
  
  labs(
    x = "Sequencing method",
    y = "Coverage (×)",
    title = "Theoretical vs deduplicated theoretical coverage"
  ) +
  
  theme_classic() +
  
  theme(
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    panel.grid.minor.y = element_line(
      color = "grey93",
      linewidth = 0.3
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  "/scratch/project_2019524/plots/fastq_multiq/tot_coverage_dedup_coverage.png",
  plot = p,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)
