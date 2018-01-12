#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Downloading FACS data zip file to $DIR"
curl https://ndownloader.figshare.com/articles/5715040/versions/1 > $DIR/facs_raw_data.zip
echo "Downloading Droplet data zip file to $DIR"
curl https://ndownloader.figshare.com/articles/5715025/versions/1 > $DIR/droplet_raw_data.zip

echo "Unzipping directories"

unzip droplet_raw_data.zip -d droplet_raw_data
unzip facs_raw_data.zip -d facs_raw_data
unzip facs_raw_data/FACS.zip -d facs_raw_data/
unzip droplet_raw_data/droplet.zip -d droplet_raw_data/
