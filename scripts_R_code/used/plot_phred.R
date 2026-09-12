#gets the phred score plots, puts R1 and R2 plots together and saves them
#as pdf and png


library(magick)

#change to the path of the file that stores all the fastqc data files
#work_dir <- "C:/Users/ezele/Downloads/lowinput_fastqc"

files <- list.files(
    work_dir,
    pattern = "_fastqc\\.zip$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
)

cat("Found", length(files), "FastQC ZIPs\n")

get_fastqc_plot <- function(zip_file) {
    contents <- unzip(zip_file, list = TRUE)$Name
    plot_name <- contents[
        grepl("Images/per_base_quality\\.png$", contents)
    ]
    
    if (length(plot_name) != 1)
        stop("Could not find per_base_quality.png in ", basename(zip_file))
    
    tmp <- tempfile("plot_")
    dir.create(tmp)
    
    unzip(zip_file, files = plot_name, exdir = tmp)
    
    plot <- list.files(
        tmp,
        pattern = "per_base_quality\\.png$",
        recursive = TRUE,
        full.names = TRUE
    )
    
    image_read(plot[1])
}

info <- data.frame(
    FileName = basename(files),
    stringsAsFactors = FALSE
)

info$Read <- case_when(
    grepl("_1_fastqc\\.zip$", info$FileName, ignore.case = TRUE) ~ "R1",
    grepl("_2_fastqc\\.zip$", info$FileName, ignore.case = TRUE) ~ "R2",
    TRUE ~ NA_character_
)

info$Sample <- sub(
    "_[12]_fastqc\\.zip$",
    "",
    info$FileName
)

# Extract the sample ID to add to the plot
info$SampleID <- sub(
    ".*-",
    "",
    info$Sample
)

print(info[, c("SampleID", "Read")])

# Create phred_scores folder INSIDE work_dir
output_dir <- file.path(
    work_dir,
    "phred_scores"
)

dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
)

samples <- unique(info$Sample)

for (sample in samples) {
    
    sample_info <- info[info$Sample == sample, ]
    
    r1_file <- sample_info$FileName[sample_info$Read == "R1"]
    r2_file <- sample_info$FileName[sample_info$Read == "R2"]
    
    if (length(r1_file) != 1 || length(r2_file) != 1) {
        warning("Missing R1 or R2 for ", sample)
        next
    }
    
    r1_path <- files[basename(files) == r1_file]
    r2_path <- files[basename(files) == r2_file]
    
    r1 <- get_fastqc_plot(r1_path)
    r2 <- get_fastqc_plot(r2_path)
    
    # Add labels to mark R1 and R2
    add_read_label <- function(img, label) {
        
        info_img <- image_info(img)
        
        img <- image_extent(
            img,
            geometry = paste0(
                info_img$width, "x",
                info_img$height + 70
            ),
            gravity = "south",
            color = "white"
        )
        
        image_annotate(
            img,
            label,
            size = 32,
            weight = 700,
            color = "black",
            gravity = "north",
            location = "+0+15"
        )
    }
    
    r1 <- add_read_label(r1, "R1")
    r2 <- add_read_label(r2, "R2")
    
    combined <- image_append(
        image_join(r1, r2),
        stack = FALSE
    )
    
    sample_id <- sample_info$SampleID[1]
    
    info_combined <- image_info(combined)
    
    combined <- image_extent(
        combined,
        geometry = paste0(
            info_combined$width, "x",
            info_combined$height + 70
        ),
        gravity = "north",
        color = "white"
    )
    
    combined <- image_annotate(
        combined,
        sample_id,
        size = 28,
        weight = 700,
        color = "black",
        gravity = "south",
        location = "+0+15"
    )
    
    image_write(
        combined,
        file.path(
            output_dir,
            paste0(sample_id, "_FastQC_per_base_quality.png")
        ),
        format = "png"
    )
    
    image_write(
        combined,
        file.path(
            output_dir,
            paste0(sample_id, "_FastQC_per_base_quality.pdf")
        ),
        format = "pdf"
    )
    
    cat("Created:", sample_id, "\n")
}

cat("\nFinished. Files saved to:\n", output_dir, "\n")
