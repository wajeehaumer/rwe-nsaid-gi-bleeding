# Downloads the Eunomia GiBleed dataset (synthetic OMOP CDM v5.3, ~2,700
# patients) and loads every table into a local DuckDB file.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

raw_dir <- here::here("data", "raw")
zip_file <- file.path(raw_dir, "GiBleed_5.3.zip")
csv_dir <- file.path(raw_dir, "GiBleed_5.3")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(zip_file)) download.file(dataset_url, zip_file, mode = "wb")
if (!dir.exists(csv_dir)) unzip(zip_file, exdir = raw_dir)

if (file.exists(paths$db)) file.remove(paths$db)
con <- connect_omop()

for (csv in list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)) {
  table <- tolower(tools::file_path_sans_ext(basename(csv)))
  dbExecute(con, sprintf(
    "CREATE TABLE %s AS SELECT * FROM read_csv('%s', header = true, normalize_names = true)",
    table, csv
  ))
}

message("Loaded ", length(dbListTables(con)), " OMOP tables; ",
        dbGetQuery(con, "SELECT COUNT(*) AS n FROM person")$n, " persons.")
dbDisconnect(con, shutdown = TRUE)
