# tabula-muris

## Folder Organization

* FACS = SmartSeq2 on FACS-sorted plates
* Microfluidic = 10x droplet-based unique molecular identifier (UMI)-barcoded transcripts and cells

```
tabula_muris/
    00_data_ingest/               # How the data was processed from gene-cell tables
        README.md
        download_robj.Rmd         # Download R objects for figures using this script
        tissues/                  # *Generate* R objects for figures yourself
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
    03_figure3/                   # All-cell clustering heatmap with dendrogram
        figure3.Rmd
    04_figure4/                   # Analysis of all T cells sorted by FACS
        figure4{a-d}.Rmd
    05_figure5/                   # Transcription factor expression analysis
        figure5.Rmd
    11_supplementary_figure1/     # Histograms of number of genes detected across tissues
    12_supplementary_figure2/     # FACS vs Microfluidics - # cells expressing a gene
    13_supplementary_figure3/     # FACS vs Microfluidics - # genes detected per cell
    14_supplementary_figure4/     # FACS vs Microfluidics - dynamic range
    15_supplementary_figure5/     # Microfluidics TSNE plots + annotation barplots
    16_supplementary_figure6/     # Analysis of dissociation-induced genes
    17_supplementary_figure7/     # Transcription factor enrichment in cell types
```
