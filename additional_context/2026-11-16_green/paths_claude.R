# Green split-half study: the paths and the design every step shares
# Drafted by Claude for Kevin Z. Lin, 2026-11-16

# Sourced by every numbered script of this folder. The study repeats the
# Gabitto split-half study (../2026-11-15_gabitto/) on the Green (ROSMAP)
# microglia cohort. The two external locations are Kevin's laptop paths
# (`GREEN_COUNTS` and `GREEN_SPLIT_OUT`); see `CLAUDE_kevin.md`, External
# Locations.

# Paths ------------------------------------------------------------------------

data_file <- file.path("~", "Library", "CloudStorage", "Dropbox",
                       "Collaboration-and-People", "archive", "tati", "data",
                       "primary_datasets", "green_2024",
                       "Writeup43_green_cleaned_only-counts.RData")
out_dir <- file.path("~", "Library", "CloudStorage", "Dropbox",
                     "Collaboration-and-People", "Kathryn Roeder - private",
                     "eSVD2", "out", "2026-11-16_green")
# A trial run of the later steps on small fits points this elsewhere.
if(Sys.getenv("GREEN_OUT_DIR") != "") out_dir <- Sys.getenv("GREEN_OUT_DIR")

# The design of Tati's Writeup47 fit -------------------------------------------

categorical_vars <- c("msex", "APOEe4_status")
numerical_vars <- c("pmi", "age_death")
case_control_var <- "AD_path_cerad"
case_control_levels <- c("FALSE", "TRUE") # control first
id_var <- "individualID"

fit_label_vec <- c("split1_half1", "split1_half2", "split2_half1",
                   "split2_half2")
