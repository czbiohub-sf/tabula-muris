# Launches all notebooks as a background process. Only reccommed to be used on high-powered computers!

RMDS := $(wildcard *.Rmd)

all:
	for RMD in $(RMDS); do \
		echo $$RMD; \
		echo "rmarkdown::render('$$RMD', clean=TRUE)" | R --slave > $$RMD.out 2>$$RMD.err & \
	done
