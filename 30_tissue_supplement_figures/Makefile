
# Generate all Rmds from template
YAMLS := $(wildcard ../00_data_ingest/02_tissue_analysis_rmd/00_parameters/*yaml)

rmd:
	for YML in $(YAMLS); do \
		python ../utilities/generate_from_template.py $$YML ; \
	done


all: rmd
	# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!
	RMDS='$(wildcard *.Rmd)'
	@echo "got Rmds"
	for RMD in $(RMDS); do \
		@echo "got to loop"; \
		@echo $$RMD; \
		echo "rmarkdown::render('$$RMD', clean=TRUE)" | R --slave > $$RMD.out 2>$$RMD.err & \
	done

