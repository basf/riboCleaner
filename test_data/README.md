
# Running the Test Data

## Overview

This test data is from our published paper.

It is RNASeq data from several tissue types in field and greenhouse locations.

Before following these steps, make sure you have already followed the [instructions for Running riboCleaner using Docker](../README.md) to become familiar with how riboCleaner works and have built the *ribocleaner-wf* Docker image. 

For this, we used Docker, but the instructions from the Running riboCleaner using Singularity section of the main README can be used from within this directory to run this analysis using Singularity.

---

## Get the Soybean Reference files

If you have cloned the riboCleaner repo you should have all the files necessary to run the test workflow *except* the soybean reference files. If you haven't already cloned the riboCleaner workflow, we recommend doing that now.

```bash
git clone https://github.com/basf/riboCleaner.git
```

Once you have the riboCleaner repo which contains the test data, you must download the soybean reference files. We recommend using the Soybean A4V1 fasta and gene model annotations (gff) from [Phytozome](https://phytozome-next.jgi.doe.gov/info/Gmax_Wm82_a4_v1) (this will require a user account).

The files needed are:

- Gmax_508_v4.0.fa.gz
- Gmax_508_Wm82.a4.v1.gene.gff3.gz

(they will need to be unzipped to be used)

Place the downloaded, decompressed references in `inputs/references`

---

## Checking the input files


After downloading the references, the files in your `inputs` directory should look like this:

```
├── inputs
│   ├── config.yaml
│   ├── raw_data
│   │   ├── Flower_Field_1.R1.fastq.gz
│   │   ├── Flower_Field_1.R2.fastq.gz
│   │   ├── Flower_Field_2.R1.fastq.gz
│   │   ├── Flower_Field_2.R2.fastq.gz
│   │   ├── Flower_GH_1.R1.fastq.gz
│   │   ├── Flower_GH_1.R2.fastq.gz
│   │   ├── Flower_GH_2.R1.fastq.gz
│   │   ├── Flower_GH_2.R2.fastq.gz
│   │   ├── Leaf_Field_1.R1.fastq.gz
│   │   ├── Leaf_Field_1.R2.fastq.gz
│   │   ├── Leaf_Field_2.R1.fastq.gz
│   │   ├── Leaf_Field_2.R2.fastq.gz
│   │   ├── Leaf_GH_1.R1.fastq.gz
│   │   ├── Leaf_GH_1.R2.fastq.gz
│   │   ├── Leaf_GH_2.R1.fastq.gz
│   │   ├── Leaf_GH_2.R2.fastq.gz
│   │   ├── Seed_Field_1.R1.fastq.gz
│   │   ├── Seed_Field_1.R2.fastq.gz
│   │   ├── Seed_Field_2.R1.fastq.gz
│   │   ├── Seed_Field_2.R2.fastq.gz
│   │   ├── Seed_GH_1.R1.fastq.gz
│   │   ├── Seed_GH_1.R2.fastq.gz
│   │   ├── Seed_GH_2.R1.fastq.gz
│   │   └── Seed_GH_2.R2.fastq.gz
│   ├── references
│   │   ├── Gmax_508_v4.0.fa
│   │   └── Gmax_508_Wm82.a4.v1.gene.gff3
│   └── sample_data.tsv

```

As a final check, make sure the [config file](./inputs/config.yaml) reflects the paths of references and our sample data. If you haven't changed them since cloning the repository and you have named your references as we have above, they should be correct

---

## Running riboCleaner

Run the code below in this directory to identify putative rDNA genes and quantify them in the test data. 

NOTE: this code will run riboCleaner with up to 8 cores. To use more or less, please adjust the `--cores` argument.


```bash
# running the full workflow (on 8 cores)
docker run \
    -v ${PWD}/inputs:/analysis/inputs \
    -v ${PWD}/results:/analysis/results \
    -v ${PWD}/.snakemake:/analysis/.snakemake \
    ribocleaner-wf \
    snakemake --cores 8 -pr all

# generating a report
# to generate this report, Snakemake requires web access
# you may need to pass proxy variables if you are on a restricted network
#-e http_proxy=${http_proxy} -e https_proxy=${https_proxy} \
docker run \
    -v ${PWD}/inputs:/analysis/inputs \
    -v ${PWD}/results:/analysis/results \
    -v ${PWD}/.snakemake:/analysis/.snakemake \
    ribocleaner-wf \
    snakemake --cores 1 --report results/riboCleaner_full_report.zip -pr all

```

This will take some time to complete and when it is done, your results can be found in `results/riboCleaner_full_report.zip`. 
