environment:
	conda env export | head -n -2 > environment.yml

environment_no_versions: environment
	cut -f 1 -d '=' environment.yml > environment_no_versions.yml

