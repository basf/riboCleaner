
# riboCleaner

riboCleaner is an open-source workflow designed to identify and quantify rRNA abundance in RNA-Seq data. riboCleaner has two main workflows: identify and quantify.

The identify workflow *identifies* putative rDNA genes in a genome of interest. First, it uses [barrnap](https://github.com/tseemann/barrnap) to identify coordinates that contain rDNA elements in a genome (fasta) and then bundles these elements into non-overlapping operons. Then it compares these regions to the predicted genes (gff) to identify genes that were predicted in overlapping coordinate space as a predicted rDNA operon. These genes are flagged as putative rDNA contaminants. Then, the rest of the genes are compared to the rDNA contaminant gene list using *blastn* and any that pass a sequence identity and coverage threshold are also flagged as rDNA.

The quantify workflow *quantifies* the abundance of rRNA transcripts in RNA-Seq samples. For this we use [salmon](https://combine-lab.github.io/salmon/) to quantify the abundance of transcripts in each sample (including the putative rDNA gene models). This allows any reads derived from rRNA to map to the rRNA transcripts and avoids biasing the counts of mRNA transcripts with homologous rRNA reads. Then, riboCleaner quantifies the abundance of the rRNA transcripts to give an estimate of rRNA prevalence in the RNA-seq data and visualizes the results. As a final step, riboCleaner removes the rRNA transcripts from the salmon counts table to yield a table of mRNA abundances that can be normalized and used for downstream analysis.

This [workflow diagram](./riboCleaner_figure_S1.png) shows how the different parts of the riboCleaner standard workflow fits together. You can also read a [more detailed description of the methods](./detailed_methods.md) or see the [Snakemake DAG](./riboCleaner.png). 

![workflow_diagram](./riboCleaner_figure_S1.png)


---

## Quick Start

If you want to see the output riboCleaner produces, we have provided the [report from running riboCleaner on a test dataset](test_data.example_report.zip)

To jump right in and run riboCleaner as fast as possible:


```bash
# clone this repo
git clone https://github.com/basf/riboCleaner.git
cd riboCleaner/test_data
```

Then follow the [instructions in the test data README](./test_data/README.md) to make sure everything is set up properly.

Once the workflow runs on the test data, replace the test data inputs with your own, delete the contents of the `results` directory, and rerun. 

---

Otherwise see our [detailed documentation for installing and running riboCleaner](./docs/manual.md).
