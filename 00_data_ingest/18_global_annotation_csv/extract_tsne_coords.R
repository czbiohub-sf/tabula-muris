library(tidyverse)
library(stringr)
library(Seurat)
library(here)

load(file=here("00_data_ingest", "11_global_robj", "FACS_all.Robj"))

df = FetchData(tiss_FACS, c('cell', 'tSNE_1', 'tSNE_2', 'cluster'))
write_csv(df, here("00_data_ingest", "18_global_annotation_csv", "tsne_facs.csv"))

load(file=here("00_data_ingest", "11_global_robj", "droplet_all.Robj"))

df_droplet = FetchData(tiss_droplet, c('cell', 'tSNE_1', 'tSNE_2', 'cluster'))
write_csv(df_droplet, here("00_data_ingest", "18_global_annotation_csv", "tsne_droplet.csv"))
