locals().update(config)

rule sample_background_windows_by_region:
    input:
        enriched_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        all_windows = ancient(FEATURE_ANNOTATIONS)
    output:
        variable_windows = "output/homer/region_matched_background/variable/{experiment_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz"
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/sample_background_windows_by_region/{experiment_label}.sample_background_windows_by_region.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.sample_background_windows_by_region.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.sample_background_windows_by_region.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting sample_background_windows_by_region" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/sample_matched_background_by_region.R \
            {input.enriched_windows} \
            {input.all_windows} \
            75 \
            output/homer/region_matched_background \
            {wildcards.experiment_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished sample_background_windows_by_region" | tee -a {log.stdout}
        """

rule run_homer:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        background = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz",
        genome = GENOME
    output:
        report = "output/homer/finemapped_results/{experiment_label}/homerResults.html",
        pwm = "output/homer/finemapped_results/{experiment_label}/homerMotifs.all.motifs",
    resources:
        mem_mb=lambda wildcards, attempt: 12000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/run_homer/{experiment_label}.all_replicates.reproducible.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.run_homer.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.run_homer.err",
    conda:
        "envs/homer.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting run_homer" | tee  -a {log.stdout}

        # Prepare process-substitution inputs
        fg_input=<(zcat {input.finemapped_windows} \
            | awk -v OFS="\t" '{{print $4 ":" $9, $1, $2+1, $3, $6}}')

        bg_input=<(zcat {input.background} \
            | awk -v OFS="\t" '{{print $4, $1, $2+1, $3, $6}}')

        # Run HOMER motif analysis
        findMotifsGenome.pl "$fg_input" {input.genome} \
            output/homer/finemapped_results/{wildcards.experiment_label} \
            -preparsedDir output/homer/preparsed \
            -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 \
            -bg "$bg_input" \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished run_homer" | tee -a {log.stdout}
        """

rule consult_encode_reference:
    input:
        enriched_windows = lambda wildcards: expand(
            "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
            experiment_label = experiment_labels
        ),
        encode_references = lambda wildcards: expand(
            TOOL_DIR + "/{reference}.reference.tsv",
            reference=["encode3_feature_summary", "encode3_eclip_enrichment", "encode3_class_assignment"]
        )
    output:
        tsne_coordinates = "output/tsne/skipper.tsne_query.tsv",
        tsne_plot = "output/figures/tsne/skipper.tsne_query.pdf"
    resources:
        mem_mb = 1000,
        runtime = "30m"
    benchmark: "benchmarks/consult_encode_reference/skipper.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.consult_encode_reference.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.consult_encode_reference.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting consult_encode_reference" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/consult_encode_reference_windows.R \
            output/reproducible_enriched_windows \
            {TOOL_DIR} \
            skipper \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished consult_encode_reference" | tee -a {log.stdout}
        """

rule consult_encode_reference_re:
    input:
        enriched_re = lambda wildcards: expand(
            "output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz",
            experiment_label = experiment_labels
        ),
        encode_references = lambda wildcards: expand(
            TOOL_DIR + "/{reference}.reference.tsv",
            reference=["encode3_feature_summary", "encode3_eclip_enrichment", "encode3_class_assignment"]
        )
    output:
        tsne_coordinates = "output/tsne_re/skipper.tsne_re_query.tsv",
        tsne_plot = "output/figures/tsne_re/skipper.tsne_re_query.pdf"
    resources:
        mem_mb = 1000,
        runtime = "30m"
    benchmark: "benchmarks/consult_encode_reference_re/skipper.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.consult_encode_reference_re.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.consult_encode_reference_re.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting consult_encode_reference" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/consult_encode_reference_re.R \
            output/reproducible_enriched_re \
            {TOOL_DIR} \
            skipper \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished consult_encode_reference_re" | tee -a {log.stdout}
        """

rule consult_term_reference:
    input:
        enriched_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
    output:
        enrichment_results = "output/gene_sets/{experiment_label}.enriched_terms.tsv.gz",
        enrichment_plot = "output/figures/gene_sets/{experiment_label}.clustered_top_terms.pdf"
    params:
        gene_sets = config["GENE_SETS"],
        gene_set_reference = config["GENE_SET_REFERENCE"],
        gene_set_distance = config["GENE_SET_DISTANCE"]
    resources:
        mem_mb = 1000,
        runtime = "30m"
    benchmark: "benchmarks/consult_term_reference/{experiment_label}.all_replicates.reproducible.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.consult_term_reference.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.consult_term_reference.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting consult_term_reference" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/consult_term_reference.R \
            {input.enriched_windows} \
            {params.gene_sets} \
            {params.gene_set_reference} \
            {params.gene_set_distance} \
            {wildcards.experiment_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished consult_term_reference" | tee -a {log.stdout}
        """
