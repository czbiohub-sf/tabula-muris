
# Generate all Rmds from template
YAMLS := $(wildcard ../00_data_ingest/02_tissue_analysis_rmd/00_parameters/*yaml)


rmd:
	for YML in $(YAMLS); do \
		python ../00_data_ingest/02_tissue_analysis_rmd/generate_from_template.py $$YML ; \
	done
