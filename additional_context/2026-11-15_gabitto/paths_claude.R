# Gabitto split-half study: the paths and the design every step shares
# Drafted by Claude for Kevin Z. Lin, 2026-11-15

# Sourced by every numbered script of this folder. The two external
# locations are Kevin's laptop paths (`TATI_ESVD_OUT`'s sibling `data/`
# folder, and the study's output folder beside the repository clone); see
# `CLAUDE_kevin.md`, External Locations.

# Paths ------------------------------------------------------------------------

data_file <- file.path("~", "Library", "CloudStorage", "Dropbox",
                       "Collaboration-and-People", "archive", "tati", "data",
                       "primary_datasets", "gabitto_2024",
                       "Writeup43_gabitto_cleaned_only-counts.RData")
out_dir <- file.path("~", "Library", "CloudStorage", "Dropbox",
                     "Collaboration-and-People", "Kathryn Roeder - private",
                     "eSVD2", "out", "2026-11-15_gabitto")
# A trial run of the later steps on small fits points this elsewhere.
if(Sys.getenv("GABITTO_OUT_DIR") != "") out_dir <- Sys.getenv("GABITTO_OUT_DIR")

# The design of Tati's Writeup45 fit -------------------------------------------

categorical_vars <- c("sex", "APOE4status")
numerical_vars <- c("PMI", "Ageatdeath")
case_control_var <- "CERAD_binary"
case_control_levels <- c("Control", "Case") # control first
id_var <- "donor_id"

fit_label_vec <- c("split1_half1", "split1_half2", "split2_half1",
                   "split2_half2")
