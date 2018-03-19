# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!

RMDS := $(wildcard *.Rmd)

clean:
	rm -rf *.out *.err *html *_files

all: clean
	export LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 ../../utilities/run_rmds.py

all_oldschool: clean
	for RMD in $(RMDS); do \
		echo $$RMD; \
		echo "rmarkdown::render('$$RMD', clean=TRUE)" | R --slave > $$RMD.out 2>$$RMD.err & \
	done
