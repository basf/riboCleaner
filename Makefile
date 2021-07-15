
build_core: rm_core
	docker build -t ribocleaner-core ./workflow/envs

rm_core:
	-docker rm ribocleaner-core

build: rm build_core
	docker build -t ribocleaner-wf -f ./Dockerfile ./workflow

rm:
	-docker rm ribocleaner-wf

dryrun:
	docker run -v ${PWD}/inputs:/analysis/inputs -v ${PWD}/results:/analysis/results -v ${PWD}/.snakemake:/analysis/.snakemake ribocleaner-wf snakemake --cores 8 -npr all

run:
	docker run -v ${PWD}/inputs:/analysis/inputs -v ${PWD}/results:/analysis/results -v ${PWD}/.snakemake:/analysis/.snakemake ribocleaner-wf snakemake --cores 8 -pr all

report:
	docker run -v ${PWD}/inputs:/analysis/inputs -v ${PWD}/results:/analysis/results -v ${PWD}/.snakemake:/analysis/.snakemake ribocleaner-wf snakemake --cores 8 --report results/report.zip -pr all

all: build run report

# pushing the image to dockerhub (after login)
push:
        docker tag ribocleaner-wf:latest basfcontainers/ribocleaner:latest
        docker push basfcontainers/ribocleaner:latest
