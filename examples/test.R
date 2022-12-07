library(devtools)
install_github("W-Holtz/R4scSHARP")
library(R4scSHARP)

# input vars
data_path <- "splat_0.8/query_counts.csv.gz"
out_path <- "output"
marker_path <- "splat_0.8/markers.txt"
ref_path <- "splat_0.8/ref_counts.csv.gz"
ref_label_path <- "splat_0.8/ref_labels.csv"

output <- run_r4scsharp(data_path, out_path, marker_path, ref_path, ref_label_path)
