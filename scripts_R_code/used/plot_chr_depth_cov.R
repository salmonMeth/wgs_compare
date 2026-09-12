library(tidyverse)

# dirs where the coverage stats are of the different libraries

method_dirs <- c(
  lowinput_P100 = "/scratch/project_2019524/pinksalmon_alignments/lowinput_P100/coverage_depth",
  lowinput       = "/scratch/project_2019524/pinksalmon_alignments/lowinput/coverage_depth",
  normalinput       = "/scratch/project_2019524/pinksalmon_alignments/normalinput/coverage_depth"
)
read_method <- function(method, directory) {
  
  files <- list.files(
    path = directory,
    pattern = "\\.coverage\\.txt$",
    full.names = TRUE
  )
  
  message(
    method, ": found ",
    length(files),
    " coverage files"
  )
  
  map_dfr(files, function(file) {
    
    sample <- str_extract(
      basename(file),
      "^\\d+"
    )

    
    dat <- read.delim(
      file,
      header = TRUE,
      sep = "\t",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    # Change "#rname" to "rname"
    names(dat)[1] <- "rname"

    
    dat %>%
      mutate(
        Sample = sample,
        Method = method
      )
  })
}


coverage <- map2_dfr(
  names(method_dirs),
  method_dirs,
  read_method
)


coverage <- coverage %>%
  mutate(
    Sample = case_when(
      Sample == "8167" ~ "8036",
      Sample == "8168" ~ "8038",
      TRUE ~ Sample
    )
  )


coverage <- coverage %>%
  mutate(
    chr = as.numeric(
      str_extract(rname, "(?<=NC_060)\\d+")
    ) - 172
  )


coverage <- coverage %>% filter(!is.na(chr))



plot_data <- coverage %>% group_by(Sample, Method, chr) %>% summarise( coverage = mean(coverage, na.rm = TRUE), .groups = "drop" )


output_dir <- "/scratch/project_2019524/pinksalmon_alignments/coverage_comparison/chr_coverages"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_data %>%
  group_split(Sample) %>%
  walk(function(dat) {
    
    plot_id <- unique(dat$Sample)
    
    # Dynamic y-axis
    min_coverage <- min(dat$coverage, na.rm = TRUE)
    max_coverage <- max(dat$coverage, na.rm = TRUE)
    
    y_min <- max(0, floor(min_coverage - 2))
    y_max <- ceiling(max_coverage + 1)
    
    p <- ggplot(
      dat,
      aes(
        x = chr,
        y = coverage,
        color = Method,
        group = Method
      )
    ) +
      geom_line(linewidth = 1) +
      geom_point(
        aes(shape = Method),
        size = 2.5
      ) +
      scale_x_continuous(
        breaks = sort(unique(dat$chr))
      ) +
      scale_y_continuous(
        limits = c(y_min, y_max),
        breaks = scales::pretty_breaks(n = 8)
      ) +
      labs(
        title = paste("Chromosome coverage —", plot_id),
        x = "Chromosome",
        y = "Coverage (%)",
        color = "Library method",
        shape = "Library method"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    print(p)
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0("Chromosome_coverage_", plot_id, ".png")
      ),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  })



###########################
###########################
#create plots for the mean depth


plot_data_depth <- coverage %>%
  group_by(Sample, Method, chr) %>%
  summarise(
    meandepth = mean(meandepth, na.rm = TRUE),
    .groups = "drop"
  )

output_dir <- "/scratch/project_2019524/pinksalmon_alignments/coverage_comparison/chr_depths"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_data_depth %>% 
  group_split(Sample) %>% 
  walk(function(dat) { 
    
    plot_id <- unique(dat$Sample) 
    
    # Dynamic y-axis 
    min_depth <- min(dat$meandepth, na.rm = TRUE) 
    max_depth <- max(dat$meandepth, na.rm = TRUE) 
    
    y_min <- max(0, floor(min_depth - 2)) 
    y_max <- ceiling(max_depth + 1) 
    
    p <- ggplot( 
      dat, 
      aes( 
        x = chr, 
        y = meandepth, 
        color = Method,
        linetype = Method,
        group = Method 
      ) 
    ) + 
      geom_line(linewidth = 1) + 
      geom_point( 
        aes(shape = Method), 
        size = 2.5 
      ) + 
      scale_x_continuous( 
        breaks = sort(unique(dat$chr)) 
      ) + 
      scale_y_continuous( 
        limits = c(y_min, y_max),
        breaks = scales::pretty_breaks(n = 8)
      ) + 
      labs( 
        title = paste("Chromosome mean depth —", plot_id), 
        x = "Chromosome", 
        y = "Mean depth", 
        color = "Library method",
        linetype = "Library method",
        shape = "Library method" 
      ) + 
      theme_classic() + 
      theme( 
        plot.title = element_text(face = "bold"), 
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )

    print(p) 
    
    ggsave( 
      filename = file.path( 
        output_dir, 
        paste0("Chromosome_mean_depth_", plot_id, ".png") 
      ), 
      plot = p, 
      width = 10, 
      height = 6, 
      dpi = 300 
    ) 
  })
