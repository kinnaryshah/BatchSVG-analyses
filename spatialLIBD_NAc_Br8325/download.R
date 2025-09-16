# download an spe and clean it 
# run this directly on JHPCE since data is stored there

# start with Br8667 since it has 4 samples across 2 slides
setwd("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc")

library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)
library(DropletUtils)
library(scuttle)
library(spatialLIBD)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

sample_info_path <- here("raw-data", "sample_key_spatial_NAc.csv")
sample_info_path2 <- here(
    "processed-data", "02_image_stitching", "sample_info_clean.csv"
)
transformed_dir <- here("processed-data", "04_VisiumStitcher")
raw_out_path <- here("processed-data", "05_harmony_BayesSpace", "01-build_spe", "spe_raw.rds")

################################################################################
#   Read in the two sources of sample info and merge
################################################################################

message("Gathering sample info")
sample_info <- read.csv(sample_info_path) |>
    as_tibble() |>
    rename(sample_id = Slide) |>
    select(c(sample_id, Age, Sex, Diagnosis, In.analysis, Refined.transforms))

sample_info <- read.csv(sample_info_path2) |>
    as_tibble() |>
    rename(sample_id = X) |>
    select(-In.analysis) |>
    right_join(sample_info, by = "sample_id") |>
    filter(In.analysis == "Yes") |>
    select(-In.analysis) |>
    mutate(spaceranger_dir = dirname(normalizePath(spaceranger_dir))) |>
    rename(slide_num = Slide.., array_num = Array.., donor = Brain)

#   We only need certain columns in the colData
sample_info <- sample_info |>
    select(
        c(
            sample_id, donor, slide_num, array_num, spaceranger_dir,
            raw_image_path, Age, Sex, Diagnosis
        )
    )

all_donors <- unique(sample_info$donor)

################################################################################
#   Build the SPE object using [slide]_[capture area] as samples, to start
################################################################################

#   'read10xVisiumWrapper' silently breaks when multiple validly named tissue
#   positions files exist for the same sample. A script for VisiumStitcher
#   supposedly depends on always having a 'tissue_positions_list.csv', which
#   creates conflicting requirements. As a temporary workaround, just rename
#   duplicate files while read10xVisiumWrapper is running
new_paths <- file.path(
    sample_info$spaceranger_dir, "spatial", "tissue_positions.csv"
)
old_paths <- file.path(
    sample_info$spaceranger_dir, "spatial", "tissue_positions_list.csv"
)
bad_old_paths <- old_paths[file.exists(old_paths) & file.exists(new_paths)]
temp_old_paths <- sub("/tissue_", "/.temp_tissue_", bad_old_paths)

stopifnot(all(file.rename(bad_old_paths, temp_old_paths)))

message("Building SpatialExperiment using [slide]_[capture area] as sample ID")
spe <- read10xVisiumWrapper(
    samples = sample_info$spaceranger_dir,
    sample_id = sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    verbose = TRUE
)

stopifnot(all(file.rename(temp_old_paths, bad_old_paths)))

save(spe, file="~/BatchSVG-analyses/spatialLIBD_NAc_Br8325/results/spatialLIBD_NAc_Br8325_spe_all.Rdata")