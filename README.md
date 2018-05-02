# tabula-muris

This is the code for annotating the tissues and cells behind the [Transcriptomic characterization of 20 organs and tissues from mouse at single cell resolution creates a Tabula Muris](https://www.biorxiv.org/content/early/2017/12/20/237446) paper.

## Installation - Python

To install the Python dependencies, create a `tabula-muris-env` environment by using the `environment.yml` file provided:

```
conda env create -f environment.yml
```

Activate the environment and install it to your Jupyter notebook with:

```
source activate tabula-muris-env
python -m ipykernel install --user --name tabula-muris-env --display-name "Python 3.6 (tabula-muris-env)"
```

## Installation - R

Packages:

```
install.packages(c("here", "Seurat", "useful", "ontologyIndex", "tidyverse"))
```

## Getting started

### From "raw" gene-cell counts tables

If you want to start from the raw gene-cell counts tables, then first download the data from figshare. You can download manually from the links ([FACS](https://figshare.com/articles/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells/5715040) and [Droplet](https://figshare.com/articles/Single-cell_RNA-seq_data_from_microfluidic_emulsion/5715025)) or run a script we've prepared:

```
bash 00_data_ingest/download_data.sh
```

This will download two zip files, `droplet_raw_data.zip` and `facs_raw_data.zip` and unzip them into the folder structure described below. Then you'll have two folders in `00_data_ingest` (the location is important - everything here depends on the folder structure). 


#### FACS

The FACS folder should look like this:

```
00_facs_raw_data
├── FACS
│   ├── Aorta-counts.csv
│   ├── Bladder-counts.csv
│   ├── Brain_Myeloid-counts.csv
│   ├── Brain_Non-Myeloid-counts.csv
│   ├── Diaphragm-counts.csv
│   ├── Fat-counts.csv
│   ├── Heart-counts.csv
│   ├── Kidney-counts.csv
│   ├── Large_Intestine-counts.csv
│   ├── Limb_Muscle-counts.csv
│   ├── Liver-counts.csv
│   ├── Lung-counts.csv
│   ├── Mammary_Gland-counts.csv
│   ├── Marrow-counts.csv
│   ├── Pancreas-counts.csv
│   ├── Skin-counts.csv
│   ├── Spleen-counts.csv
│   ├── Thymus-counts.csv
│   ├── Tongue-counts.csv
│   └── Trachea-counts.csv
├── FACS.zip
├── annotations_FACS.csv
└── metadata_FACS.csv
```

#### Droplet

Now your droplet folders should look like this:

```
01_droplet_raw_data
├── annotations_droplet.csv
├── droplet
│   ├── Bladder-10X_P4_3
│   ├── Bladder-10X_P4_4
│   ├── Bladder-10X_P7_7
│   ├── Heart_and_Aorta-10X_P7_4
│   ├── Kidney-10X_P4_5
│   ├── Kidney-10X_P4_6
│   ├── Kidney-10X_P7_5
│   ├── Limb_Muscle-10X_P7_14
│   ├── Limb_Muscle-10X_P7_15
│   ├── Liver-10X_P4_2
│   ├── Liver-10X_P7_0
│   ├── Liver-10X_P7_1
│   ├── Lung-10X_P7_8
│   ├── Lung-10X_P7_9
│   ├── Lung-10X_P8_12
│   ├── Lung-10X_P8_13
│   ├── Mammary_Gland-10X_P7_12
│   ├── Mammary_Gland-10X_P7_13
│   ├── Marrow-10X_P7_2
│   ├── Marrow-10X_P7_3
│   ├── Spleen-10X_P4_7
│   ├── Spleen-10X_P7_6
│   ├── Thymus-10X_P7_11
│   ├── Tongue-10X_P4_0
│   ├── Tongue-10X_P4_1
│   ├── Tongue-10X_P7_10
│   ├── Trachea-10X_P8_14
│   └── Trachea-10X_P8_15
├── droplet.zip
└── metadata_droplet.csv
```

All of the `*-10X_*` folders contain a `barcodes.tsv`, `genes.tsv` and `matrix.mtx` file as output by cellranger from 10X genomics.

```
01_droplet_raw_data/droplet/Bladder-10X_P4_3
├── barcodes.tsv
├── genes.tsv
└── matrix.mtx

```



## Folder Organization

* FACS = SmartSeq2 on FACS-sorted plates
* Microfluidic = 10x droplet-based unique molecular identifier (UMI)-barcoded transcripts and cells

```
tabula_muris/
    00_data_ingest/               # How the data was processed from gene-cell tables
        README.md
        download_robj.Rmd         # Download R objects for figures using this script
        02_tissue_analysis_rmd/                  # *Generate* R objects for figures yourself
            Aorta_facs.Rmd
            Brain-Non-microglia_facs.Rmd
            Brain-Microglia_facs.Rmd
            Bladder_facs.Rmd
            Bladder_droplet.Rmd
            Colon_facs.Rmd
            Heart_facs.Rmd
            Heart_droplet.Rmd
            ... more files ...
        03_tissue_annotation_csv/
            Aorta_facs_annotation.csv
            Brain-Non-microglia_facs_annotation.csv
            Brain-Microglia_facs_annotation.csv
            Bladder_facs_annotation.csv
            Bladder_droplet_annotation.csv
            Colon_facs_annotation.csv
            Heart_facs_annotation.csv
            Heart_droplet_annotation.csv
            ... more files ...
        04_tissue_robj_generated/
        10_tissue_robj_downloaded/
        11_global_robj/
        12_extract_number_of_genes_cells/
        13_ngenes_ncells_facs/
        14_ngenes_ncells_droplet/
        15_color_palette/
        16_genes_for_tissue_tsne/
        20_dissociation_genes/
        All_Droplet_Notebook.Rmd
        All_FACS_Notebook.Rmd
        Droplet_Notebook.Rmd
        FACS_Notebook.Rmd
        README.md
        cell_order_FACS.txt
        cell_order_droplets.txt
        download_data.sh
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

## Tissues and platforms studied

|                     | FACS + SmartSeq2 | Microfluidic droplets (10x) |
|:--------------------|:-----------------|:----------------------------|
| Aorta               | Yes              | No                          |
| Bladder             | Yes              | Yes                         |
| Brain_Microglia     | Yes              | No                          |
| Brain_Non-microglia | Yes              | No                          |
| Colon               | Yes              | No                          |
| Diaphragm           | Yes              | No                          |
| Fat                 | Yes              | No                          |
| Heart               | Yes              | Yes                         |
| Kidney              | Yes              | Yes                         |
| Liver               | Yes              | Yes                         |
| Lung                | Yes              | Yes                         |
| Mammary             | Yes              | Yes                         |
| Marrow              | Yes              | Yes                         |
| Muscle              | Yes              | Yes                         |
| Pancreas            | Yes              | No                          |
| Skin                | Yes              | No                          |
| Spleen              | Yes              | Yes .                       |
| Thymus              | Yes              | Yes                         |
| Tongue              | Yes              | Yes                         |
| Trachea             | Yes              | Yes                         |


## Other mouse cell atlas efforts

### "Mouse Cell Atlas"

From Guoji Guo's lab at Zhejiang University School of Medicine in China

- [MCA Data Analysis](https://github.com/ggjlab/mca_data_analysis)
- [scMCA](https://github.com/ggjlab/scMCA) - An simple R package for large scale data (large DGE) from Mouse Cell Atlas, to alleviate burdens of our main Server

