library(yaml)
yaml_data <- yaml::read_yaml("inst/envs/M1.yaml")
r_data <- gsub("([[:alnum:]_]+)=([[:alnum:].]+)=.+", "\\1==\\2" , yaml_data$dependencies)
cat(paste0(
  "c(",
  "\n",
  paste0("  ", "'", r_data, "'", collapse = ",\n"),
  "\n",
  ")",
  "\n"
))
