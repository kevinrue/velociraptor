m1_info <- read.table("M1_20240522.txt")[, 1:2]
intel_info <- read.table("Intel_20240601.txt")[, 1:2]
merged_info <- merge(m1_info, intel_info, by = "V1", suffixes = c(".m1", ".intel"))
merged_info[merged_info$V2.m1 != merged_info$V2.intel, ]
