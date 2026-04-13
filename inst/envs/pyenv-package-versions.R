packages_info <- basilisk::listPackages(velociraptor:::velo.env)$full
cat(paste0(
  "c\n",
  paste0(sprintf("  '%s'", packages_info), collapse = ",\n"),
  "\n)"
))
