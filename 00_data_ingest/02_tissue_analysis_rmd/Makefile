# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!

RMDS := $(wildcard *.Rmd)

clean:
	rm -rf *.out *.err *html *_files

all: clean
	../../utilities/run_rmds.py

all_oldschool: clean
	for RMD in $(RMDS); do \
		echo $$RMD; \
		echo "rmarkdown::render('$$RMD', clean=TRUE)" | R --slave > $$RMD.out 2>$$RMD.err & \
	done
