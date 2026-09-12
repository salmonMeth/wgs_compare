library(tidyverse)

df <- read_csv("/scratch/project_2019524/pinksalmon_alignments/alignment_comparison/MAPQ_distribution.csv") %>%
  mutate(
    MAPQ = as.numeric(MAPQ),
    Count = as.numeric(Count),
    Sample = as.character(Sample),
    
    # Samples that should be plotted together
    PlotSample = case_when(
      Sample == "8167" ~ "8036",
      Sample == "8168" ~ "8038",
      TRUE ~ Sample
    )
  )

df %>%
  group_split(PlotSample) %>%
  walk(function(dat) {
    
    plot_id <- unique(dat$PlotSample)
    
    # Combine counts from samples belonging to the same plot
    # for each library method and MAPQ
    density_data <- dat %>%
      group_by(Method, MAPQ) %>%
      summarise(
        Count = sum(Count),
        .groups = "drop"
      ) %>%
      
      # Calculate weighted KDE for each method
      group_by(Method) %>%
      group_modify(~ {
        
        x <- .x$MAPQ
        w <- .x$Count
        
        # KDE bandwidth
        bw <- bw.nrd0(x)
        
        grid <- seq(0, 60, length.out = 512)
        
        dens <- sapply(grid, function(g) {
          sum(w * dnorm(g, mean = x, sd = bw)) / sum(w)
        })
        
        tibble(
          MAPQ = grid,
          Density = dens
        )
      }) %>%
      ungroup()
    
    p <- ggplot(
      density_data,
      aes(
        x = MAPQ,
        y = Density,
        color = Method,
        group = Method
      )
    ) +
      geom_line(linewidth = 1.2) +
      scale_x_continuous(
        limits = c(0, 60),
        breaks = seq(0, 60, 5)
      ) +
      labs(
        title = paste("MAPQ probability density —", plot_id),
        x = "MAPQ",
        y = "Probability density",
        color = "Library method"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    
    print(p)
    
    # Save plot
    ggsave(
      filename = paste0("/scratch/project_2019524/pinksalmon_alignments/alignment_comparison/MAPQ_densities_", plot_id, ".png"),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  })

