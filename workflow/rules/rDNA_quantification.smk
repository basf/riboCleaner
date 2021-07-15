
# rules associated with running salmon 
rule _salmon_index:
    """ build the index for mapping """
    input:
        features=rules._get_feature_fasta.output.features,
        genome=config["reference"]["genome"]
    threads: max_threads
    output: directory(os.path.join(RDNA_QUANT, "salmon.index"))
    threads: max_threads
    shell:
        """
        salmon index -p {threads} -t <(cat {input.features} {input.genome}) -d <(grep "^>" {input.genome} | cut -d " " -f 1 | sed 's/>//g') -i {output}
        """

ruleorder: _run_salmon > _run_salmon_single

rule _run_salmon:
    """ run salmon on each of the samples """
    input:
        r1=rules._subsample_reads.output.r1,
        r2=rules._subsample_reads.output.r2,
        index = rules._salmon_index.output
    threads: max_threads
    output:
        sf=os.path.join(RDNA_QUANT, "salmon", "{sample}.salmon.sf"),
        tmp=directory(os.path.join(RDNA_QUANT, "salmon", "{sample}.salmon")),
        meta_info=os.path.join(RDNA_QUANT, "salmon", "{sample}.meta_info.json")
    log: os.path.join(RDNA_QUANT, "salmon", "{sample}.salmon.log")
    shell:
        """
        salmon quant --threads {threads} -i {input.index} --libType A -1 {input.r1} -2 {input.r2} -o {output.tmp} 2>&1 | tee {log}
        cp -v {output.tmp}/quant.sf {output.sf}
        cp -v {output.tmp}/aux_info/meta_info.json {output.meta_info}
        """

rule _run_salmon_single:
    """ run salmon on each of the samples (single end mode) """
    input:
        r1=rules._subsample_reads_single.output.r1,
        index=rules._salmon_index.output
    threads: max_threads
    output:
        sf=os.path.join(RDNA_QUANT, "salmon", "{sample}.salmon.sf"),
        tmp=directory(os.path.join(RDNA_QUANT, "salmon", "{sample}.salmon")),
        meta_info=os.path.join(RDNA_QUANT, "salmon", "{sample}.meta_info.json")
    log: os.path.join(RDNA_QUANT, "salmon", "{sample}.salmon.log")
    shell:
        """
        salmon quant --threads {threads} -i {input.index} --libType A -r {input.r1} -o {output.tmp} 2>&1 | tee {log}
        cp -v {output.tmp}/quant.sf {output.sf}
        cp -v {output.tmp}/aux_info/meta_info.json {output.meta_info}
        """


rule quantify_rdna:
    """ Reports the proportion of counts that correspond to a feature identified as rDNA """
    input:
        sf=expand(rules._run_salmon.output.sf, sample=samples.index),
        meta_info=expand(rules._run_salmon.output.meta_info, sample=samples.index),
        rdna_names=rules.rDNA_prediction.input.names
    output:
        rdna_abundance=report(os.path.join(RDNA_QUANT, "rDNA_gene_TPM.tsv"), caption="../report/rDNA_gene_TPM.rst", category="Quantification"),
        summary=report(os.path.join(RDNA_QUANT, "rDNA_summary.tsv"), caption="../report/rDNA_summary.rst", category="Quantification")
    params:
        samples=samples.index
    threads: 1
    run:
        import pandas as pd
        import json

        # get the rDNA names
        with open(input.rdna_names, 'r') as IN:
            names = [line.strip() for line in IN]
        
        sample_summary = {}

        # gather the counts for each file
        all_counts = []
        for sf, meta_info, sample in zip(input.sf, input.meta_info, params.samples):
            with open(meta_info) as IN:
                info = json.load(IN)

            d = pd.read_csv(sf, sep="\t", header=0, index_col=0)
            counts = d["TPM"]

            counts.name = sample
            all_counts.append(counts)

            # get the sample summary info
            sample_summary[sample] = {
                'percent_unmapped': 100 - float(info['percent_mapped']),
                # scale percent rDNA to the mapped portion
                'percent_rDNA': float(info['percent_mapped']) * counts.filter(items=names).sum() / counts.sum()
                }

            # the usable percent is whatever is left over
            sample_summary[sample]['percent_mRNA'] = 100 -  sample_summary[sample]['percent_unmapped'] - sample_summary[sample]['percent_rDNA']

        # and combine into one df
        all_counts = pd.DataFrame(all_counts).fillna(0).transpose()

        # make the index be just the gene names (not the position)
        all_counts.index = [g.split(":")[0] for g in all_counts.index]

        # get just the rdna
        all_counts_rdna = all_counts.filter(items=names, axis=0)
        all_counts_rdna.to_csv(output.rdna_abundance, sep="\t", index_label="sample")

        # done as is above to write out more details in the future, for now, just the TPM
        summary = pd.DataFrame.from_dict(sample_summary, orient="index").round(2)

        # write sorted by sample name
        summary.sort_index().to_csv(output.summary, sep="\t", index_label="sample", float_format="%.2f")


rule clean_tpm:
    """ Generates a TPM abundance table with rDNA removed (and TPM rescaled to the remaining mRNA) """
    input:
        sf=expand(rules._run_salmon.output.sf, sample=samples.index),
        meta_info=expand(rules._run_salmon.output.meta_info, sample=samples.index),
        rdna_names=rules.rDNA_prediction.input.names
    output:
        clean_tpm=report("results/mRNA_tpm.tsv", caption="../report/mRNA_tpm.rst", category="Quantification")
    params:
        samples=samples.index
    threads: 1
    run:
        import pandas as pd
        import json

        # get the rDNA names
        with open(input.rdna_names, 'r') as IN:
            rdna = [line.strip() for line in IN]
        
        sample_summary = {}

        # gather the counts for each file
        all_counts = []
        for sf, sample in zip(input.sf, params.samples):
            d = pd.read_csv(sf, sep="\t", header=0, index_col=0)
            counts = d["TPM"]

            # get the total (including rDNA)
            total = counts.sum()

            # remove rDNA
            counts = counts.filter(items=[n for n in counts.index if n not in rdna])

            # rescale counts by the rDNA ratio (to get the TPM as if the rDNA was never there)
            counts = counts * (total / counts.sum())

            counts.name = sample
            all_counts.append(counts)

        # and combine into one df
        all_counts = pd.DataFrame(all_counts).fillna(0).transpose()

        # make the index be just the gene names (not the position)
        all_counts.index = [g.split(":")[0] for g in all_counts.index]

        all_counts.to_csv(output.clean_tpm, sep="\t", index_label="gene")

# vim: set syntax=python expandtab tabstop=4 shiftwidth=4 autoindent foldmethod=indent foldnestmax=1:
