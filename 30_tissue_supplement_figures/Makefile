
# Generate all Rmds from template
YAMLS := $(wildcard ../00_data_ingest/02_tissue_analysis_rmd/00_parameters/*yaml)

rmd:
	for YML in $(YAMLS); do \
		LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 python ../utilities/generate_from_template.py $$YML ; \
	done


all: rmd
	# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!
	../utilities/run_rmds.py
