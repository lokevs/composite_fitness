#### Fetch data ----

# Downloads the data and fitted-model files archived on Zenodo into data/
# and models/, verifying each file against the MD5 checksum reported by
# Zenodo. Files already present with the correct checksum are skipped.
#
# Usage, from the repository root:
#   Rscript fetch_data.R            # published record on zenodo.org
#   Rscript fetch_data.R sandbox    # test record on sandbox.zenodo.org

record_ids <- c(zenodo = "22297476",
                sandbox = "597737")

# Files the scripts require (the parameter estimates are derived from the fitted
# models but are needed to source R/functions.R before step 1 regenerates them);
# anything else in the record is downloaded too
required_files <- c("thermal_performance.txt", "ovi_rate_conversion.txt", "era5_thermal_regimes.txt",
                    "parameter_estimates_LRS.txt", "parameter_estimates_dev_rate.txt", "parameter_estimates_growth_rate.txt",
                    "LRS_model.rds", "devrate_model.rds", "growthrate_model.rds")

args <- commandArgs(trailingOnly = TRUE)
site <- if (length(args) > 0 && args[1] == "sandbox") "sandbox" else "zenodo"
host <- if (site == "sandbox") "https://sandbox.zenodo.org" else "https://zenodo.org"
record_url <- paste0(host, "/api/records/", record_ids[site])

options(timeout = max(3600, getOption("timeout")))
for (d in c("data", "models", "output")) dir.create(d, showWarnings = FALSE)

##########

#### Read the file list from the Zenodo record ----

json_file <- tempfile(fileext = ".json")
download.file(record_url, json_file, quiet = TRUE)
record <- paste(readLines(json_file, warn = FALSE), collapse = "")

# Each file entry lists its name, size and checksum before its download links
file_keys <- regmatches(record, gregexpr('"key": ?"[^"]+"', record))[[1]]
file_names <- sub('"key": ?"([^"]+)"', "\\1", file_keys)
key_positions <- gregexpr('"key": ?"[^"]+"', record)[[1]]
file_info <- do.call(rbind, lapply(seq_along(file_names), function(i) {
  block_end <- if (i < length(file_names)) key_positions[i + 1] else nchar(record)
  block <- substr(record, key_positions[i], block_end)
  data.frame(name = file_names[i],
             size = as.numeric(sub('.*"size": ?([0-9]+).*', "\\1", block)),
             md5 = sub('.*"checksum": ?"md5:([0-9a-f]+)".*', "\\1", block),
             url = paste0(host, "/api/records/", record_ids[site], "/files/", file_names[i], "/content"),
             stringsAsFactors = FALSE)
}))

missing <- setdiff(required_files, file_info$name)
if (length(missing) > 0) stop("Record ", record_ids[site], " lacks required files: ", paste(missing, collapse = ", "))
cat("Record", record_ids[site], "on", host, "lists", nrow(file_info), "files\n")

##########

#### Download and verify ----

for (i in seq_len(nrow(file_info))) {

  f <- file_info[i, ]
  destination <- file.path(ifelse(grepl("[.]rds$", f$name), "models", "data"), f$name)

  if (file.exists(destination) && unname(tools::md5sum(destination)) == f$md5) {
    cat("  - Present, checksum OK:", destination, "\n")
    next
  }

  cat("  - Downloading", f$name, sprintf("(%.1f MB)", f$size / 1e6), "\n")
  download.file(f$url, destination, mode = "wb", quiet = TRUE)

  if (unname(tools::md5sum(destination)) != f$md5) stop("Checksum mismatch after downloading ", f$name)
  cat("    Checksum OK\n")
}

cat("All files present and verified.\n")

##########
