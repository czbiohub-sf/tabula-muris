environment:
	conda env export | head -n -2 > environment.yml

environment_no_versions: environment
	cut -f 1 -d '=' environment.yml > environment_no_versions.yml

connect_ndnd_figures:
	sshfs olga@ndnd.czbiohub.org:/home/olga/tabula-muris/30_tissue_supplement_figures $HOME/tabula-muris-supplemental

download_data:
	bash 00_data_ingest/download_data.sh