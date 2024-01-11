locals().update(config)
rule run_homer:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        background = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz",
        genome = GENOME
    output:
        report = "output/homer/finemapped_results/{experiment_label}/homerResults.html"
    params:
        error_file = "stderr/{experiment_label}.run_homer.err",
        out_file = "stdout/{experiment_label}.run_homer.out",
        run_time = "00:40:00",
        memory = "8000",
        job_name = "run_homer"
    benchmark: "benchmarks/run_homer/{experiment_label}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:Homer_4.11"
    resources:
        mem_mb=8000
    shell:
        "findMotifsGenome.pl <(zcat {input.finemapped_windows} | awk -v OFS=\"\t\" '{{print $4 \":\"$9,$1,$2+1,$3,$6}}') "
            "{input.genome} output/homer/finemapped_results/{wildcards.experiment_label} -preparsedDir output/homer/preparsed -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 "
            "-bg <(zcat {input.background} | awk -v OFS=\"\t\" '{{print $4,$1,$2+1,$3,$6}}') "

rule consult_encode_reference:
    input:
        enriched_windows = lambda wildcards: expand("output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz", experiment_label = experiment_labels),
        enriched_re = lambda wildcards: expand("output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", experiment_label = experiment_labels),
        encode_references = lambda wildcards: expand(TOOL_DIR + "/{reference}.reference.tsv", reference = ["encode3_feature_summary", "encode3_eclip_enrichment", "encode3_class_assignment"])
    output:
        tsne_coordinates = "output/tsne/skipper.tsne_query.tsv",
        tsne_plot = "output/figures/tsne/skipper.tsne_query.pdf"
    params:
        error_file = "stderr/skipper.consult_encode_reference.err",
        out_file = "stdout/skipper.consult_encode_reference.out",
        run_time = "00:10:00",
        memory = "1000",
        job_name = "consult_encode_reference"
    benchmark: "benchmarks/consult_encode_reference/skipper.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=1000
    shell:
        "Rscript --vanilla {TOOL_DIR}/consult_encode_reference.R output/reproducible_enriched_windows output/reproducible_enriched_re {TOOL_DIR} skipper "

rule consult_term_reference:
    input:
        enriched_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        gene_sets = GENE_SETS,
        gene_set_reference = GENE_SET_REFERENCE,
        gene_set_distance = GENE_SET_DISTANCE
    output:
        enrichment_results = "output/gene_sets/{experiment_label}.enriched_terms.tsv.gz",
        enrichment_plot = "output/figures/gene_sets/{experiment_label}.clustered_top_terms.pdf"
    params:
        error_file = "stderr/{experiment_label}.consult_term_reference.err",
        out_file = "stdout/{experiment_label}.consult_term_reference.out",
        run_time = "00:15:00",
        memory = "1000",
        job_name = "consult_term_reference"
    benchmark: "benchmarks/consult_term_reference/{experiment_label}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"'
    resources:
        mem_mb=1000
    shell:
        "Rscript --vanilla {TOOL_DIR}/consult_term_reference.R {input.enriched_windows} {input.gene_sets} {input.gene_set_reference} {input.gene_set_distance} {wildcards.experiment_label} "
