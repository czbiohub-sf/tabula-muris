# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!

RMDS := $(wildcard *.Rmd)
CORES=32

clean:
	rm -rf *.out *.err *_files

all: clean
	export LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 ../../utilities/run_rmds.py

all_oldschool: clean
	for RMD in $(RMDS); do \
		echo $$RMD; \
		echo "rmarkdown::render('$$RMD', clean=TRUE)" | R --slave > $$RMD.out 2>$$RMD.err & \
	done

create_R_files_from_Rmds: clean
	for RMD in $(RMDS); do \
		echo $$RMD; \
		echo "library(knitr); purl('$$RMD')" | R --vanilla; \
	done

generate_Robj_files: create_R_files_from_Rmds
	ls -1 *.R | grep -v boiler | grep -v Template | xargs -P ${CORES} -I{} bash -c "R -f {} >{}.out 2>{}.err"

skin_heart: create_R_files_from_Rmds
	ls -1 *.R | grep facs | grep -E 'Skin|Heart' | xargs -P ${CORES} -I{} bash -c "R -f {} >{}.out 2>{}.err"

trachea: create_R_files_from_Rmds
	ls -1 *.R | grep droplet | grep -E 'Trachea' | xargs -P ${CORES} -I{} bash -c "R -f {} >{}.out 2>{}.err"

mammary: create_R_files_from_Rmds
	ls -1 *.R | grep facs | grep -E 'Mammary' | xargs -P ${CORES} -I{} bash -c "R -f {} >{}.out 2>{}.err"
