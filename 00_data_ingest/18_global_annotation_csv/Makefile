clean:
	rm *csv

droplet:
	./concat_all_annotations_by_method.py --method droplet

facs:
	./concat_all_annotations_by_method.py --method facs

all: clean facs droplet
