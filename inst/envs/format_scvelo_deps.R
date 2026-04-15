library(dplyr)
library(tidyr)

pkg_info <- read.table("/Users/kevin/Downloads/test-scvelo/pip.list.txt", skip = 2)
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