#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Downloading FACS data zip file to $DIR"
curl https://ndownloader.figshare.com/articles/5829687/versions/1 > $DIR/00_facs_raw_data.zip
echo "Downloading Droplet data zip file to $DIR"
curl https://ndownloader.figshare.com/articles/5715025/versions/1 > $DIR/01_droplet_raw_data.zip

echo "Unzipping directories"

unzip 00_facs_raw_data.zip -d 00_facs_raw_data
unzip 01_droplet_raw_data.zip -d 01_droplet_raw_data
unzip 00_facs_raw_data/FACS.zip -d 00_facs_raw_data/
unzip 01_droplet_raw_data/droplet.zip -d 01_droplet_raw_data/

rm 00_facs_raw_data.zip
rm 01_droplet_raw_data.zip
