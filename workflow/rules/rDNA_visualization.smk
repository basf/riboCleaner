
rule plot_abundance_summary:
    """ Plots the rDNA abundance at the sample level and at the experimental level """
    input:
        summary=rules.quantify_rdna.output.summary
    output:
        # output png for convenience and svg for publications
        experiment_summary=report(os.path.join(RDNA_QUANT, "experiment_summary.png"), caption="../report/experiment_summary.rst", category="Quantification"),
        experiment_summary_svg=report(os.path.join(RDNA_QUANT, "experiment_summary.svg"), caption="../report/experiment_summary.rst", category="Quantification"),
        sample_summary=report(os.path.join(RDNA_QUANT, "sample_summary.png"), caption="../report/sample_summary.rst", category="Quantification"),
        sample_summary_svg=report(os.path.join(RDNA_QUANT, "sample_summary.svg"), caption="../report/sample_summary.rst", category="Quantification")
    threads: 1
    run:
        import pandas as pd
        from matplotlib import pyplot as plt

        summary = pd.read_csv(input.summary, sep="\t", index_col=0)[["percent_mRNA", "percent_rDNA", "percent_unmapped"]]
        summary.columns = ["mRNA", "rDNA", "unmapped"]

        # make the experiment level plot
        # adapted from: https://www.dataforeverybody.com/matplotlib-seaborn-pie-charts/
        overall_summary = summary.mean(axis=0)
        pie, ax = plt.subplots(figsize=[10,6])
        labels = overall_summary.keys()
        plt.pie(x=overall_summary, autopct="%.2f%%", explode=[0.05]*3, labels=labels, pctdistance=0.5)
        plt.title("Experiment Overview", fontsize=14)
        pie.savefig(output.experiment_summary)
        pie.savefig(output.experiment_summary_svg)

        # make the sample level plot
        # limited to 5 samples per inch of plot (and 3 inches for overhead)
        width = 3 + int(len(summary) / 5)
        summary.plot(kind="bar", stacked=True, figsize=(width, 8))
        plt.legend(loc='lower left', bbox_to_anchor=(1.0, 0.5))
        plt.xticks(rotation=45, ha="right")
        plt.title("Sample Composition")
        plt.ylabel("% of Reads")
        plt.tight_layout()
        plt.gcf().savefig(output.sample_summary)
        plt.gcf().savefig(output.sample_summary_svg)


# vim: set syntax=python expandtab tabstop=4 shiftwidth=4 autoindent foldmethod=indent foldnestmax=1:
