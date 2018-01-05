# tabula-muris


## Folder Organization


* FACS = SmartSeq2 on FACS-sorted plates
* Microfluidic = 10x droplet-based unique molecular identifier (UMI)-barcoded transcripts and cells

```
tabula_muris/
    00_data_ingest/               # How the data was processed from gene-cell tables
        README.md
        download_robj.Rmd
        Aorta_FACS.Rmd
        Brain-Non-microglia_FACS.Rmd
        Brain-Microglia_FACS.Rmd
        Colon_FACS.Rmd
        Heart_FACS.Rmd
        Heart_Microfluidic.Rmd
        ... more files ...
    01_figure1/                   # Overview + #cell barplots + #gene/#reads horizonplots
        README.md
        figure1{b-g}.ipynb
    02_figure2/                   # FACS TSNE plots + annotation barplots
        README.md
        figure2a.Rmd
        figure2b.Rmd
        figure2c.ipynb
    03_figure3/                   # All-cell clustering heatmap with dendrogram.
        figure3.Rmd
    04_figure4/                   # Analysis of all T cells sorted by FACS.
        figure4{a-d}.Rmd
    05_figure5/                   # Transcription factor expression analysis.
        figure5.Rmd
```
