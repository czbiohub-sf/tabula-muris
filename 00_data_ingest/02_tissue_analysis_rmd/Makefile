# Launches all notebooks as a background process. Only recommended to be used on high-powered computers!

RMDS := $(wildcard *.Rmd)

all:
	../../utilities/run_rmds.py

all_oldschool:
	for RMD in $(RMDS); do \
		echo $$RMD; \
		echo "rmarkdown::render('$$RMD', clean=TRUE)" | R --slave > $$RMD.out 2>$$RMD.err & \
	done
