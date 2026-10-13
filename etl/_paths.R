# Bundle-relative paths shared by the ETL scripts (works from any location / working directory).
# Sourced after R/utils.R, which defines ROOT.

RAW       <- path_root("raw")
SINAN_DIR <- file.path(RAW, "sinan")
INMET_DIR <- file.path(RAW, "inmet")
IBGE_DIR  <- file.path(RAW, "ibge")
REF       <- path_root("reference")
ORIG_DB   <- file.path(REF, "original_db")
DB        <- path_root("db")
CLIN      <- path_root("clinical_tables")

for (.d in c(DB, CLIN)) dir.create(.d, showWarnings = FALSE, recursive = TRUE)

DIM_LOCATION <- file.path(REF, "dim_location.csv")
