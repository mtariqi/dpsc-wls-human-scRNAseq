library(Seurat)

samples <- c(
  "GSM4998457_Pulp1",
  "GSM4998458_Pulp2",
  "GSM4998459_Pulp3",
  "GSM4998460_Pulp4",
  "GSM4998461_Pulp5"
)

objs <- list()

for(i in seq_along(samples)){

  prefix <- samples[i]

  counts <- ReadMtx(
    mtx = paste0("data/GSE164157/", prefix, "_matrix.mtx.gz"),
    features = paste0("data/GSE164157/", prefix, "_genes.tsv.gz"),
    cells = paste0("data/GSE164157/", prefix, "_barcodes.tsv.gz")
  )

  objs[[i]] <- CreateSeuratObject(
    counts = counts,
    project = prefix
  )

}

combined <- merge(
  objs[[1]],
  y = objs[-1]
)

saveRDS(
  combined,
  "data/processed/GSE164157_raw.rds"
)

print(combined)
