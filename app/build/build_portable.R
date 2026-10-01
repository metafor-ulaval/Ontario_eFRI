# 🟡🟡 Build the app folder of the deliverable 🟡🟡 ----
# Run with the R version to ship (R 4.4.x), on a computer where the app runs:
#   Rscript build_portable.R <deliverable app folder> <eFRItools source folder>
#
# 1. Copy the app files (app.R, R/, Launch_eFRI.bat)
# 2. Copy this R installation with only its base and recommended packages
# 3. Copy the packages needed by the app, and their dependencies, from the libraries of this R
# 4. Install eFRItools from its source folder
# 5. Write the list of shipped packages





# 🟡🟡 Parameters 🟡🟡 ----
args <- commandArgs(trailingOnly = TRUE)
dest <- normalizePath(args[1], winslash = "/", mustWork = FALSE)
efritools_src <- normalizePath(args[2], winslash = "/", mustWork = TRUE)

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
app_src <- dirname(dirname(normalizePath(script, winslash = "/")))

app_packages <- c("shiny", "shinyjs", "shinyBS", "leaflet", "dplyr", "purrr", "stringr", "tibble",
                  "magrittr", "sf", "terra", "mapedit", "callr", "ps", "processx")

r_portable <- paste0(dest, "/R-portable")
lib_portable <- paste0(r_portable, "/library")





# 🟡🟡 Functions 🟡🟡 ----
copy_dir <- function(from, to) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(from, full.names = TRUE, all.files = TRUE, no.. = TRUE), to, recursive = TRUE, copy.date = TRUE)
  if (!all(ok)) stop("Copy failed: ", from)
}





# 🟡🟡 1. App files 🟡🟡 ----
cat("Copy app files\n")
dir.create(paste0(dest, "/R"), recursive = TRUE, showWarnings = FALSE)
file.copy(paste0(app_src, c("/app.R", "/Launch_eFRI.bat")), dest, overwrite = TRUE)
file.copy(list.files(paste0(app_src, "/R"), full.names = TRUE), paste0(dest, "/R"), overwrite = TRUE)





# 🟡🟡 2. R installation 🟡🟡 ----
cat("Copy R installation\n")
if (dir.exists(r_portable)) stop("R-portable already exists, delete it first: ", r_portable)

for (d in c("bin", "etc", "modules", "share", "Tcl")) copy_dir(paste0(R.home(), "/", d), paste0(r_portable, "/", d))

ip <- installed.packages()
base_packages <- unique(rownames(ip)[ip[, "Priority"] %in% "base"])
recommended_packages <- unique(rownames(ip)[ip[, "Priority"] %in% "recommended"])
for (pkg in base_packages) copy_dir(paste0(.Library, "/", pkg), paste0(lib_portable, "/", pkg))





# 🟡🟡 3. App packages and dependencies 🟡🟡 ----
cat("Copy app packages\n")
efritools_deps <- read.dcf(paste0(efritools_src, "/DESCRIPTION"), fields = "Imports")[1, 1]
efritools_deps <- trimws(sub("\\(.*", "", strsplit(efritools_deps, ",")[[1]]))

ip <- ip[!duplicated(ip[, "Package"]), ] # First library of .libPaths() wins, like library()

needed <- unique(c(app_packages, efritools_deps, recommended_packages))
needed <- unique(c(needed, unlist(tools::package_dependencies(needed, db = ip, which = c("Depends", "Imports"), recursive = TRUE))))
needed <- setdiff(needed, c(base_packages, "eFRItools"))

missing <- setdiff(needed, rownames(ip))
if (length(missing) != 0) stop("Packages not installed: ", paste(missing, collapse = ", "))

for (pkg in needed) copy_dir(paste0(ip[pkg, "LibPath"], "/", pkg), paste0(lib_portable, "/", pkg))





# 🟡🟡 4. eFRItools 🟡🟡 ----
cat("Install eFRItools\n")
status <- system2(paste0(r_portable, "/bin/R.exe"),
                  c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", paste0("--library=", shQuote(lib_portable)), shQuote(efritools_src)))
if (status != 0) stop("eFRItools installation failed")





# 🟡🟡 5. List of shipped packages 🟡🟡 ----
shipped <- installed.packages(lib.loc = lib_portable)[, c("Package", "Version", "Priority")]
write.csv(shipped, paste0(r_portable, "/packages.csv"), row.names = FALSE)
cat("Done:", nrow(shipped), "packages in", lib_portable, "\n")
