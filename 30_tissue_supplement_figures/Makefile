
# Generate all Rmds from template
YAMLS := $(wildcard ../28_tissue_yamls_for_supplement/*yaml)

clean_rmd:
	rm -rf *_auto_generated.Rmd

rmd: clean
	for YML in $(YAMLS); do \
		LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 python \
			../utilities/generate_from_template.py \
			--suffix _auto_generated.Rmd $$YML ; \
	done

clean:
	rm -rf *.out *.err
	rm -rf *html
	rm -rf */

run_rmd: clean
	cp -r ../27_tissue_additional_figures/* .
	LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 ../utilities/run_rmds.py


all: rmd clean
	# copy over microbiome figures
	cp -r ../27_tissue_additional_figures .
	# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!
	LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 ../utilities/run_rmds.py
