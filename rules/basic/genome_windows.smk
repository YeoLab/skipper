locals().update(config)

rule filter_gff:
    input:
        gff = ancient(GFF),
        rankings = ancient(ACCESSION_RANKINGS),
    output:
        gff_filt = "output/gff/filtered.gff3",
        rankings_filt = "output/gff/filtered_ranks.txt",
    params:
        source = config["GFF_SOURCE"]
    threads: 1
    resources:
        mem_mb = 32000,
        runtime = "1h"
    benchmark: "benchmarks/filter_gff.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/filter_gff.out",
        stderr = config["WORKDIR"] + "/stderr/filter_gff.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting filter_gff" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/filter_gff.R \
            {params.source} \
            {input.gff} \
            {input.rankings} \
            {output.gff_filt} \
            {output.rankings_filt} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished filter_gff" | tee -a {log.stdout}
        """

rule parse_gff:
    input:
        gff_filt = "output/gff/filtered.gff3",
        rankings = ancient(ACCESSION_RANKINGS),
    output:
        partition = PARTITION,
        feature_annotations = FEATURE_ANNOTATIONS,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 64000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 180 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/parse_gff.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/parse_gff.out",
        stderr = config["WORKDIR"] + "/stderr/parse_gff.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting parse_gff" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/parse_gff.R \
            {input.gff_filt} \
            {input.rankings} \
            {output.partition} \
            {output.feature_annotations} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished parse_gff" | tee -a {log.stdout}
        """

rule partition_bam_reads:
    input:
        chrom_sizes = config["CHROM_SIZES"],
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
        region_partition = PARTITION,
    output:
        counts = "output/secondary_results/counts/genome/vectors/{replicate_label}.counts",
    params:
        uninformative = config["UNINFORMATIVE_READ"]
    resources:
        mem_mb=lambda wildcards, attempt: 48000 * (1.5 ** (attempt - 1)),
        tmpdir = TMPDIR,
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/counts/unassigned_experiment.{replicate_label}.partition_bam_reads.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.partition_bam_reads.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.partition_bam_reads.err",
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting partition_bam_reads" | tee -a {log.stdout}

        bedtools bamtobed -i {input.bam} \
            | awk '($1 != "chrEBV") && ($4 !~ "/{params.uninformative}$")' \
            | bedtools flank -s -l 1 -r 0 -g {input.chrom_sizes} -i - \
            | bedtools shift -p 1 -m -1 -g {input.chrom_sizes} -i - \
            | LC_ALL=C sort -S 50% -T "{resources.tmpdir}" -k1,1 -k2,2n \
            | perl -lane 'print if @F==6' \
            | bedtools coverage -counts -s -a {input.region_partition} -b - \
            | cut -f 7 \
            | awk 'BEGIN {{print "{wildcards.replicate_label}"}} {{print}}' \
        >> {output.counts} 2> {log.stderr}

        echo "[`date`] Finished partition_bam_reads" | tee -a {log.stdout}
        """

rule calc_partition_nuc:
    input:
        partition = PARTITION,
        genome = GENOME
    output:
        nuc = PARTITION.replace(".bed", ".nuc")
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/partition_nuc.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/calc_partition_nuc.out",
        stderr = config["WORKDIR"] + "/stderr/calc_partition_nuc.err",
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting calc_partition_nuc" | tee -a {log.stdout}

        bedtools nuc -s \
            -fi {input.genome} \
            -bed {input.partition} \
            | gzip -c \
            > {output.nuc} \
        2> {log.stderr}

        echo "[`date`] Finished calc_partition_nuc" | tee -a {log.stdout}
        """

rule make_genome_count_table:
    input:
        partition = PARTITION.replace(".bed", ".nuc"),
        replicate_counts = lambda wildcards: expand(
            "output/secondary_results/counts/genome/vectors/{replicate_label}.counts",
            replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]
        )
    output:
        count_table = "output/secondary_results/counts/genome/tables/{experiment_label}.tsv.gz"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 30 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_genome_count_table.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.make_genome_count_table.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.make_genome_count_table.err",
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting make_genome_count_table" | tee -a {log.stdout}

        paste <(
            zcat {input.partition} \
                | awk -v OFS='\t' 'BEGIN {{print "chr\tstart\tend\tname\tscore\tstrand\tgc"}} NR > 1 {{print $1,$2,$3,$4,$5,$6,$8}}'
        ) {input.replicate_counts} \
            | gzip -c \
            > {output.count_table} \
        2> {log.stderr}

        echo "[`date`] Finished make_genome_count_table" | tee -a {log.stdout}
        """

rule fit_input_betabinomial_model:
    input:
        table = rules.make_genome_count_table.output.count_table
    output:
        coef = "output/secondary_results/input_model_coef/{experiment_label}.{input_replicate_label}.tsv"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/betabinomial/{experiment_label}.{input_replicate_label}.fit_input.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.{input_replicate_label}.fit_input_betabinomial_model.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.{input_replicate_label}.fit_input_betabinomial_model.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting fit_input_betabinomial_model" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/fit_input_betabinom.R \
            {input.table} \
            {wildcards.experiment_label} \
            {wildcards.input_replicate_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished fit_input_betabinomial_model" | tee -a {log.stdout}
        """

rule fit_clip_betabinomial_model:
    input:
        table = rules.make_genome_count_table.output.count_table
    output:
        coef = "output/clip_model_coef/{experiment_label}.{clip_replicate_label}.tsv"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/fit_clip_betabinomial_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_model.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_model.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting fit_clip_betabinomial_model" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/fit_clip_betabinom.R \
            {input.table} \
            {wildcards.experiment_label} \
            {wildcards.clip_replicate_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished fit_clip_betabinomial_model" | tee -a {log.stdout}
        """

rule call_enriched_windows:
    input:
        feature_annotations = ancient(FEATURE_ANNOTATIONS),
        accession_rankings = ancient(ACCESSION_RANKINGS),
        replicate = lambda wildcards: (
            "output/secondary_results/counts/genome/vectors/" +
            re.sub(r"IP_\d$", "IP_2", wildcards.clip_replicate_label) +
            ".counts"
        ),
        table = rules.make_genome_count_table.output.count_table,
        parameters = lambda wildcards: (
            "output/secondary_results/" + OVERDISPERSION_MODE + "_model_coef/{experiment_label}." +
            overdispersion_replicate_lookup[wildcards.clip_replicate_label] +
            ".tsv"
        )
    output:
        "output/secondary_results/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_data.tsv",
        "output/secondary_results/tested_windows/{experiment_label}.{clip_replicate_label}.tested_windows.tsv.gz",
        "output/secondary_results/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz",
        "output/secondary_results/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_feature_summary.tsv",
        "output/secondary_results/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_transcript_summary.tsv",
        "output/secondary_results/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_gene_summary.tsv",
        "output/secondary_results/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions_feature_data.tsv",
        "output/secondary_results/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_data.tsv",
        "output/secondary_results/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_transcript_data.tsv",
        "output/secondary_results/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_gc_data.tsv",
        "output/figures/secondary_figures/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_scan.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_coverage.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_rates.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.linear.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.log10.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.feature.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.all_transcript_types.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.select_transcript_types.pdf",
        "output/figures/secondary_figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.per_gene_feature.pdf",
        "output/figures/secondary_figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions.feature.pdf",
        "output/figures/secondary_figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature.pdf",
        "output/figures/secondary_figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.all_transcript_types.pdf",
        "output/figures/secondary_figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature_gc.pdf"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 24000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/call_enriched_windows/{experiment_label}.{clip_replicate_label}.call_enriched_windows.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.{clip_replicate_label}.call_enriched_windows.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.{clip_replicate_label}.call_enriched_windows.err",
    conda:
        "envs/skipper_R.yaml"
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label]
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting call_enriched_windows" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/call_enriched_windows.R \
            {input.table} \
            {input.accession_rankings} \
            {input.feature_annotations} \
            {input.parameters} \
            {params.input_replicate_label} \
            {wildcards.clip_replicate_label} \
            {wildcards.experiment_label}.{wildcards.clip_replicate_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished call_enriched_windows" | tee -a {log.stdout}
        """

rule check_window_concordance:
    input:
        windows = lambda wildcards: expand(
            "output/secondary_results/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz",
            clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label]
        )
    output:
        "output/figures/secondary_figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf",
        "output/secondary_results/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.tsv",
        "output/secondary_results/enrichment_reproducibility/{experiment_label}.odds_data.tsv"
    resources:
        mem_mb = 8000,
        runtime = "30m"
    benchmark: "benchmarks/check_window_concordance/{experiment_label}.all_replicates.concordance.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.check_window_concordance.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.check_window_concordance.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail
        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting check_window_concordance" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/check_window_concordance.R \
            output/secondary_results/tested_windows \
            {wildcards.experiment_label}
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished check_window_concordance" | tee -a {log.stdout}
        """

rule find_reproducible_enriched_windows:
    input:
        windows = lambda wildcards: expand(
            "output/secondary_results/enriched_windows/{{experiment_label}}.{clip_replicate_label}.enriched_windows.tsv.gz",
            clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label]
        )
    output:
        reproducible_windows = "output/secondary_results/unfiltered_reproducible_enriched_windows/{experiment_label}.unfiltered_reproducible_enriched_windows.tsv.gz",
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/find_reproducible_enriched_windows/{experiment_label}.all_replicates.reproducible.txt"
    params:
        blacklist = (BLACKLIST if BLACKLIST is not None else "")
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.find_reproducible_enriched_windows.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.find_reproducible_enriched_windows.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting find_reproducible_enriched_windows" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/identify_reproducible_windows.R \
            output/secondary_results/enriched_windows/ \
            {wildcards.experiment_label} \
            {params.blacklist} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished find_reproducible_enriched_windows" | tee -a {log.stdout}
        """

rule filter_reproducible_windows:
    input:
        unfiltered_enriched_windows = "output/secondary_results/unfiltered_reproducible_enriched_windows/{experiment_label}.unfiltered_reproducible_enriched_windows.tsv.gz",
        nt_coverage = "output/finemapping/nt_coverage/{experiment_label}.nt_coverage.bed"
    output:
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        linear_bar = "output/figures/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_window_counts.linear.pdf",
        log_bar = "output/figures/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_window_counts.log10.pdf",
        filtered_out = "output/secondary_results/filtered_out_windows/{experiment_label}.filtered_out_windows.tsv.gz"
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    params:
        filter = config["GINI_CUTOFF"]
    benchmark: "benchmarks/filter_reproducible_windows/{experiment_label}.all_replicates.reproducible.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.filter_reproducible_windows.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.filter_reproducible_windows.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting filter_reproducible_windows" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/reproducible_windows_filtration.R \
            {input.unfiltered_enriched_windows} \
            {input.nt_coverage} \
            {params.filter} \
            {wildcards.experiment_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished filter_reproducible_windows" | tee -a {log.stdout}
        """