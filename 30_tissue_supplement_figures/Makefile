
# Generate all Rmds from template
YAMLS := $(wildcard ../28_tissue_yamls_for_supplement/*yaml)

clean:
	rm -rf *_auto_generated.Rmd

rmd: clean
	for YML in $(YAMLS); do \
		LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 python \
			../utilities/generate_from_template.py \
			--suffix _auto_generated.Rmd $$YML ; \
	done


all: rmd
	# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!
	../utilities/run_rmds.py
