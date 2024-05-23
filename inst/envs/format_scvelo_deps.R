library(dplyr)
library(tidyr)

pkg_info <- read.table("M1_20240522.txt")
pkg_specs <- pkg_info %>%
    as_tibble() %>%
    unite("spec", V1, V2, sep = "==") %>%
    pull(spec)

cat(paste0(
    "c(",
    "\n",
    paste0("  ", "'", pkg_specs, "'", collapse = ",\n"),
    "\n",
    ")",
    "\n"
))

# paste the result in R/basilisk.R