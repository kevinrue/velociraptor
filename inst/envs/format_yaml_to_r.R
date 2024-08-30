library(yaml)
yaml_data <- yaml::read_yaml("Windows.yaml")

which_pip_dependencies <- which(sapply(yaml_data$dependencies, function(x) any(names(x) %in% "pip")))

conda_dependencies <- yaml_data$dependencies[-c(which_pip_dependencies)]

conda_dependencies <- gsub("([[:alnum:]_]+)=([[:alnum:].]+)=.+", "\\1==\\2" , conda_dependencies)


message("Conda dependencies:")
cat(paste0(
  "c(",
  "\n",
  paste0("  ", "'", conda_dependencies, "'", collapse = ",\n"),
  "\n",
  ")",
  "\n"
))

message("pip dependencies:")
if (length(which_pip_dependencies)) {
  pip_dependencies <- yaml_data$dependencies[[which_pip_dependencies]]$pip
  cat(paste0(
    "c(",
    "\n",
    paste0("  ", "'", pip_dependencies, "'", collapse = ",\n"),
    "\n",
    ")",
    "\n"
  ))
}
