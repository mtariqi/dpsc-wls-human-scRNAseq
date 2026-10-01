library(GEOquery)

datasets <- c(
  "GSE164157",
  "GSE202476",
  "GSE185222"
)

for(ds in datasets){

  cat("\nDownloading:", ds, "\n")

  tryCatch({

    getGEOSuppFiles(
      GEO = ds,
      makeDirectory = TRUE,
      baseDir = "data/raw"
    )

  }, error=function(e){

    cat("Failed
