library(yaml)
yaml_data <- yaml::read_yaml("inst/envs/Linux_x86_64_anaconda.yaml")
r_data <- gsub("([[:alnum:]_]+)=([[:alnum:].]+)=.+", "\\1==\\2" , yaml_data$dependencies)
cat(paste0(
  "c(",
  "\n",
  paste0("  ", "'", r_data, "'", collapse = ",\n"),
  "\n",
  ")",
  "\n"
))
