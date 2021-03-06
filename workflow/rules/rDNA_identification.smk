
# rules associated with rDNA prediction

rule _barrnap:
    """ Runs barrnap to search for rDNA models """
    input:
        reference=config["reference"]["genome"]
    output:
        gff=os.path.join(RDNA_PRED, "barrnap_rdna.gff3") 
    params:
        kingdom=config["kingdom"]
    threads: max_threads
    shell:
        """
        barrnap --kingdom {params.kingdom} --threads {threads} {input.reference} > {output.gff}
        """

rule _infer_non_overlapping_repeats:
    """ 
    merge to infer non-overlapping repeats
    
    this will merge neighboring rDNA elements together into one feature and extract the features with a certain number of elements.
    """
    input:
        gff=rules._barrnap.output.gff,
    output:
        nor=os.path.join(RDNA_PRED, "barrnap_rdna.nor.bed")
    params:
        nor_distance=config.get("nor_distance", 5000),
        nor_counts=config.get("nor_counts", 5)
    threads: 1
    shell:
        """
        bedtools merge -c 1 -o count -d {params.nor_distance} -i {input.gff} | grep -v ^unmapped | awk -v c={params.nor_counts} '{{if ($4>c) print}}' > {output.nor}
        """

rule _get_rdna_models:
    """ 
    Make the final call of the rDNA genes by overlapping the gene models in the gff file with the predicted rDNA models and rDNA regions
    
    This will depend somewhat on the annotations (e.g. mRNA) and may not work for all gff files
    """
    input:
        gene_models=config["reference"]["gff"],
        barrnap=rules._barrnap.output.gff,
        nor=rules._infer_non_overlapping_repeats.output.nor
    output:
        rdna_models=os.path.join(RDNA_PRED, "putative_rdna_models.gff3"),
        tmp=temp(os.path.join(RDNA_PRED, "_putative_rdna_models.tmp"))
    params:
        gff_feature=config.get("gff_feature", "mRNA")
    threads: 1
    shell:
        """
        bedtools intersect -a {input.gene_models} -b {input.barrnap} {input.nor} > {output.tmp}
        
        # grep will exit with a 1 if it finds no matches
        set +e
        grep {params.gff_feature} {output.tmp} > {output.rdna_models}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        # no matches were found
        then
            exit 0
        else
            exit $exitcode
        fi
        """
        
rule _gff2bed:
    """ Converts the GFF to a BED so we can use the names """
    input:
        gff=config["reference"]["gff"]
    output:
        bed=os.path.join(RDNA_PRED, "reference.features.bed")
    params:
        # the attribute in the gff that corresponds to the name
        name="Name",
        # the feature to use to make the bed file
        feature=config.get("gff_feature", "mRNA")
    threads: 1
    run:
        # parse through the gff and extract the fields required for the BED
        # [chromosome, start, end, name]
        with open(input.gff, 'r') as IN, open(output.bed, 'w') as OUT:
            for line in IN:
                if line.startswith("#"):
                    continue

                columns = line.strip().split("\t")

                chrom = columns[0]
                feature_type = columns[2]
                start = columns[3]
                end = columns[4]

                # only grab the features we want
                if feature_type == params.feature:
                    # get a string like "ID=Glyma.13G020400.1.Wm82.a2.v1;Name=Glyma.13G020400.1;pacid=30501300;longest=1;Parent=Glyma.13G020400.Wm82.a2.v1"
                    attrs = columns[8]
                
                    # tokenize attrs
                    tokens = {attr.split("=")[0]: attr.split("=")[1] for attr in attrs.split(";")}
                    try:
                        name = tokens[params.name]
                    except KeyError:
                        raise ValueError("'{name}' missing in '{line}'".format(name=params.name, line=line))
        
                    # write the line in the BED file
                    OUT.write("\t".join((chrom, start, end, name)) + "\n")

rule _get_rdna_fasta:
    """ Extract fasta for the putative genes """
    input:
        reference=rules._barrnap.input.reference,
        models=rules._get_rdna_models.output.rdna_models
    output:
        fasta=os.path.join(RDNA_PRED, "putative_rdna_models.fasta")
    threads: 1
    shell:
        "bedtools getfasta -name -fi {input.reference} -bed {input.models} -s > {output.fasta}"

rule _dedupe_rdna_fasta:
    """ 
    Deduplicates the rDNA for more efficient searching 
    
    The output of this might also be interesting to report (e.g. number of duplicates)
    """
    input:
        fasta=rules._get_rdna_fasta.output.fasta
    output:
        deduped=os.path.join(RDNA_PRED, "putative_rdna_models.deduped.fasta")
    params:
        identity=97
    threads: max_threads
    shell:
        "dedupe.sh minidentity={params.identity} in={input.fasta} out={output.deduped}"

rule _get_barrnap_rdna_names:
    """ Extract the names of the sequences from the gff file """
    input:
        models=rules._get_rdna_models.output.rdna_models
    output:
        names=report(os.path.join(RDNA_PRED, "putative_rdna_models.txt"), caption="../report/predicted_rdna_models.rst", category="Identification")
    params:
        # the attribute in the gff that corresponds to the name
        name="ID"
    threads: 1
    run:
        import re

        # gff/gtf have different ways to specify attributes
        gff_regex = re.compile('(?P<key>[^=]+) *= *(?P<value>[^;]+);? *')
        gtf_regex = re.compile('(?P<key>[^ ]+) *\"(?P<value>[^;\"]+)\"?;? *')

        filetype = None
        names = set()
        with open(input.models, 'r') as IN:
            # read the first line to determine if we are using gff/gtf
            first_line = IN.readline()
            if first_line:
                attrs = first_line.split("\t")[8]
                if gff_regex.match(attrs):
                    attr_regex = gff_regex
                    name_key = "ID"
                elif gtf_regex.match(attrs):
                    attr_regex = gtf_regex
                    name_key = "gene_id"
                else:
                    raise ValueError("Couldn't filetype as either GFF/GTF")

            # reset to the beginning of the file
            IN.seek(0)

            for line in IN:
                # get a string like "ID=Glyma.13G020400.1.Wm82.a2.v1;Name=Glyma.13G020400.1;pacid=30501300;longest=1;Parent=Glyma.13G020400.Wm82.a2.v1"
                # or 'gene_id "FBgn0267517"; gene_symbol "5.8SrRNA-Psi:CR45857"';
                # see: http://uswest.ensembl.org/info/website/upload/gff.html
                attrs = line.split("\t")[8]
                
                # tokenize attrs
                #tokens = {attr.split("=")[0]: attr.split("=")[1] for attr in attrs.split(";")}

                for match in attr_regex.finditer(attrs):
                    if match.group("key") == name_key:
                        names.add(match.group("value"))
                        break
                else:
                    raise ValueError("Couldn't detect gene name in '{line}'".format(line=line))
        
        # write the names 1/line to output
        with open(output.names, 'w') as OUT:
            for n in names:
                OUT.write(n + "\n")

rule _get_feature_fasta:
    """
    Extract the features from the genome and the GFF file

    For prokaryotes, this would often just be the CDS, for eukaryotes, this is typically the mRNA.

    Some genomes have annotated pseudogenes with exons that will be written by gffread. We currently make no assumptions about pseudogenes.
    To remove pseudogenes in gffread, the `--no-pseudo` option can be given. 
    """
    input:
        fasta=rules._barrnap.input.reference,
        gff=rules._get_rdna_models.input.gene_models
    output:
        features=os.path.join(RDNA_PRED, "reference.features.fasta")
    shell:
        """ 
        # -w gets the full trascript (+ UTR, etc)
        gffread -w {output.features} -g {input.fasta} {input.gff}
        """

rule _blast_for_more_rdna:
    """ 
    Some of the models may not be predicted by barrnap (because they are partial or duplicates, etc)

    This BLAST search will also add genes with high homology and coverage to rDNA genes
    """
    input:
        subject=rules._get_rdna_fasta.output.fasta,
        query=rules._get_feature_fasta.output.features
    output:
        hits=os.path.join(RDNA_PRED, "homology_rdna_models.blast6.tsv")
    shell:
        "blastn -query {input.query} -subject {input.subject} -outfmt '6 qseqid sseqid qstart qend pident qcovs' -out {output.hits} || touch {output.hits}"

rule _get_names_from_blast:
    """ 
    Parse the blast results to get the names of genes with high homology to rDNA 
    """
    input:
        blast=rules._blast_for_more_rdna.output.hits,
        barrnap_names=rules._get_barrnap_rdna_names.output.names
    output:
        names=report(os.path.join(RDNA_PRED, "homology_rdna_models.txt"), caption="../report/homology_rdna_names.rst", category="Identification")
    params:
        pident=90,
        # the hits to a single subject must cover 90% of this query
        qcovs=90
    run:
        rdna = set()
        with open(input.blast, 'r') as IN:
            current_query = ""
            current_subject = ""
            current_weight_id = 0
            current_weight_n = 0
            current_cov = 0
            for line in IN:
                query, subject, start, end, pident, qcovs = line.strip().split("\t")

                # trim the bedtools names at the first ":"
                #query = query.split(":")[0]

                # gather the HSPs to each subject for each query and check them against our thresholds
                if query != current_query or subject != current_subject:
                    # check the coverage and weighted perc id of the hits
                    if current_query != "" and float(current_cov) >= params.qcovs and current_weight_id / current_weight_n >= params.pident:
                        rdna.add(current_query)

                    # reset
                    current_query = query
                    current_subject = subject
                    current_weight_id = 0
                    current_weight_n = 0
                    current_cov = qcovs

                # add this hit; keeping a running % identity for the HSPs in this query/subject pairing
                current_weight_id += float(pident) * abs(int(end) - int(start))
                current_weight_n += abs(int(end) - int(start))
                        
        # subtract the known names
        barrnap_names = set()
        with open(input.barrnap_names, 'r') as IN:
            for line in IN:
                barrnap_names.add(line.strip())
                
        with open(output.names, 'w') as OUT:
            for n in rdna.difference(barrnap_names):
                OUT.write(n + "\n")
            
rule _get_rdna_names:
    """ Combine the blast and barrnap names """
    input:
        barrnap=rules._get_barrnap_rdna_names.output.names,
        homology=rules._get_names_from_blast.output.names
    output:
        names=report(os.path.join(RDNA_PRED, "all_rdna_models.txt"), caption="../report/all_rdna_names.rst", category="Identification")
    shell:
        "cat {input.barrnap} {input.homology} > {output.names}"

rule _write_rdna_fasta:
    """ Write a FASTA of the identified rDNA names """
    input:
        names=rules._get_rdna_names.output.names,
        features=rules._get_feature_fasta.output.features
    output:
        fasta=report(os.path.join(RDNA_PRED, "all_rdna_models.fasta"), caption="../report/all_rdna_fasta.rst", category="Identification")
    shell:
        "/opt/bbmap/filterbyname.sh include=t substring=name in={input.features} names={input.names} out={output.fasta}"        


rule _summarize_identification:
    """ Sometimes riboCleaner will fail to detect any rDNA genes. Warn the user if this is the case """
    input:
        barrnap_predictions=rules._barrnap.output.gff,
        rdna_regions=rules._infer_non_overlapping_repeats.output.nor,
        homology_models=rules._get_names_from_blast.output.names,
        gene_models=rules._get_rdna_names.output.names
    params:
        gff_feature=rules._get_rdna_models.params.gff_feature
    output:
        summary=report(os.path.join(RDNA_PRED, "identification_SUMMARY.txt"), caption="../report/identification_SUMMARY.rst", category="Identification")
    run:
        from collections import defaultdict
        # count non-partial rDNA predicted by barrnap
        barrnap_results = defaultdict(lambda: defaultdict(int))
        with open(input.barrnap_predictions, 'r') as IN:
            for line in IN:
                if not line.startswith("#"):
                    # detect if this is partial
                    partial = "partial" in line

                    # quantify each rDNA type
                    if "5S_rRNA" in line:
                        barrnap_results["5S"][partial] += 1
                    elif "18S_rRNA" in line:
                        barrnap_results["18S"][partial] += 1
                    elif "28S_rRNA" in line:
                        barrnap_results["28S"][partial] += 1
                    elif "5_8S_rRNA" in line:
                        barrnap_results["5.8S"][partial] += 1
                    # needs to be expanded for prokaryotes 
                    else:
                        barrnap_results["other"][partial] += 1
        
        # count the rDNA regions
        rDNA_regions = 0
        with open(input.rdna_regions, 'r') as IN:
            for line in IN:
                rDNA_regions += 1

        # count the models added via BLAST homology
        homology_models = 0
        with open(input.homology_models, 'r') as IN:
            for line in IN:
                homology_models += 1

        # count the rDNA genes
        rDNA_gene_models = 0
        with open(input.gene_models, 'r') as IN:
            for line in IN:
                rDNA_gene_models += 1

        # output the summary
        with open(output.summary, 'w') as OUT:
            OUT.write("riboCleaner ran barrnap to detect the following counts of rDNA features:\n")

            for elem in ["18S", "5.8S", "28S", "5S", "other"]:
                OUT.write("    {elem}: {complete} complete, {partial} partial\n".format(elem=elem, complete=barrnap_results[elem][False], partial=barrnap_results[elem][True]))

            OUT.write("\n")
            OUT.write("These rDNA features were merged into {} rDNA-rich regions.\n".format(rDNA_regions))
            OUT.write("\n")
            OUT.write("The rDNA regions were found to overlap with {overlap} features with the `{tag}` tag in the GFF file.\n".format(overlap=rDNA_gene_models - homology_models, tag=params.gff_feature))
            OUT.write("Additionally, {homology} models were identified as having significant homology to these features for a total of {total} models flagged as putative rDNA.\n".format(homology=homology_models, total=rDNA_gene_models))

            OUT.write("\n")

            # begin the warning section
            if not rDNA_regions:
                OUT.write("WARNING: There were no rDNA regions detected. This may be because there are no rDNA sequences in your reference genome (detectable by barrnap) or because they do not form regular rDNA repeats. riboCleaner works best when it is able to detect rDNA regions.\n")

            if not rDNA_gene_models:
                OUT.write("WARNING: There were no rDNA gene models detected. This is often good news but means riboCleaner will detect no rDNA for this genome. However, since there are no false gene models originating from rDNA, your counts table resulting from standard transcriptome analysis won't be biased by false rDNA gene models!\n")

# vim: set syntax=python expandtab tabstop=4 shiftwidth=4 autoindent foldmethod=indent foldnestmax=1:
