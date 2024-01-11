locals().update(config)
rule parse_gff:
    input:
        gff = ancient(GFF),
        rankings = ancient(ACCESSION_RANKINGS),
    output:
        partition = PARTITION,
        feature_annotations = FEATURE_ANNOTATIONS,
    threads: 1
    params:
        error_file = "stderr/parse_gff.err",
        out_file = "stdout/parse_gff.out",
        run_time = "3:00:00",
        job_name = "parse_gff",
        memory = "48000"
    benchmark: "benchmarks/parse_gff.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=48000
    shell:
        "Rscript --vanilla {TOOL_DIR}/parse_gff.R {input.gff} {input.rankings} {output.partition} {output.feature_annotations}"

rule partition_bam_reads:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
        region_partition = PARTITION,
    output:
        counts = "output/counts/genome/vectors/{replicate_label}.counts",
    params:
        error_file = "stderr/{replicate_label}.partition_bam_reads.err",
        out_file = "stdout/{replicate_label}.partition_bam_reads.out",
        run_time = "12:00:00",
        memory = "60000",
        job_name = "partition_bam_reads"
    benchmark: "benchmarks/counts/unassigned_experiment.{replicate_label}.partition_bam_reads.txt"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    resources:
        mem_mb=60000
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
        "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
        "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
        "bedtools sort -i - | "
        "bedtools coverage -counts -s -a {input.region_partition} -b - | cut -f 7 | "
        "awk 'BEGIN {{print \"{wildcards.replicate_label}\"}} {{print}}' > {output.counts};"
        
rule calc_partition_nuc:
    input:
        partition = PARTITION,
        genome = GENOME
    output:
        nuc = PARTITION.replace(".bed", ".nuc")
    params:
        error_file = "stderr/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "2:00:00",
        memory = "16000",
        job_name = "calc_partition_nuc"
    benchmark: "benchmarks/partition_nuc.txt"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    resources:
        mem_mb=16000
    shell:
        "bedtools nuc -s -fi {input.genome} -bed {input.partition} | gzip -c > {output.nuc}"

rule make_genome_count_table:
    input:
        partition = rules.calc_partition_nuc.output.nuc,
        replicate_counts = lambda wildcards: expand("output/counts/genome/vectors/{replicate_label}.counts", replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        count_table = "output/counts/genome/tables/{experiment_label}.tsv.gz",
    threads: 4
    params:
        error_file = "stderr/{experiment_label}.make_count_table.err",
        out_file = "stdout/{experiment_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = "1000",
        job_name = "make_genome_count_table"
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_genome_count_table.txt"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    resources:
        mem_mb=1000
    shell:
        "paste <(zcat {input.partition} | awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\\tgc\"}} NR > 1 {{print $1,$2,$3,$4,$5,$6,$8}}' ) {input.replicate_counts} | gzip -c > {output.count_table};"

rule fit_input_betabinomial_model:
    input:
        table = rules.make_genome_count_table.output.count_table
    output:
        coef = "output/input_model_coef/{experiment_label}.{input_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/input_distributions/{{experiment_label}}.{{input_replicate_label}}.{other_label}.input_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    threads: 4
    params:
        error_file = "stderr/{experiment_label}.{input_replicate_label}.fit_input_betabinom.err",
        out_file = "stdout/{experiment_label}.{input_replicate_label}.fit_input_betabinom.out",
        run_time = "6:00:00",
        memory = "32000",
        job_name = "fit_input_betabinomial_model"
    benchmark: "benchmarks/betabinomial/{experiment_label}.{input_replicate_label}.fit_input.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=32000
    shell:
        "Rscript --vanilla {TOOL_DIR}/fit_input_betabinom.R {input.table} {wildcards.experiment_label} {wildcards.input_replicate_label}"

rule fit_clip_betabinomial_model:
    input:
        table = rules.make_genome_count_table.output.count_table
    output:
        coef = "output/clip_model_coef/{experiment_label}.{clip_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment_label}}.{{clip_replicate_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    threads: 2
    params:
        error_file = "stderr/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_model.err",
        out_file = "stdout/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_model.out",
        run_time = "6:00:00",
        memory = "32000",
        job_name = "fit_clip_betabinomial_model"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=32000
    shell:
        "Rscript --vanilla {TOOL_DIR}/fit_clip_betabinom.R {input.table} {wildcards.experiment_label} {wildcards.clip_replicate_label}"

rule call_enriched_windows:
    input:
        feature_annotations = ancient(FEATURE_ANNOTATIONS),
        accession_rankings = ancient(ACCESSION_RANKINGS),
        replicate = lambda wildcards: "output/counts/genome/vectors/" + re.sub("IP_\d$","IP_2",wildcards.clip_replicate_label) + ".counts",
        table = rules.make_genome_count_table.output.count_table,
        parameters = lambda wildcards: "output/" + OVERDISPERSION_MODE + "_model_coef/{experiment_label}." + overdispersion_replicate_lookup[wildcards.clip_replicate_label] + ".tsv",
        # parameters = lambda wildcards: "output/clip_model_coef/{experiment_label}.{wildcards.clip_replicate_label}.tsv",
    output:
        "output/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_data.tsv",
        "output/tested_windows/{experiment_label}.{clip_replicate_label}.tested_windows.tsv.gz",
        "output/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_feature_summary.tsv",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_transcript_summary.tsv",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_gene_summary.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions_feature_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_transcript_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_gc_data.tsv",
        "output/figures/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_scan.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_coverage.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_rates.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.linear.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.log10.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.feature.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.all_transcript_types.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.select_transcript_types.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.per_gene_feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions.feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.all_transcript_types.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature_gc.pdf"
    threads: 2
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label],
        error_file = "stderr/{experiment_label}.{clip_replicate_label}.call_enriched_windows.err",
        out_file = "stdout/{experiment_label}.{clip_replicate_label}.call_enriched_windows.out",
        run_time = "06:00:00",
        memory = "24000",
        job_name = "call_enriched_windows"
    benchmark: "benchmarks/call_enriched_windows/{experiment_label}.{clip_replicate_label}.call_enriched_windows.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=24000
    shell:
        "Rscript --vanilla {TOOL_DIR}/call_enriched_windows.R {input.table} {input.accession_rankings} {input.feature_annotations} {input.parameters} {params.input_replicate_label} {wildcards.clip_replicate_label} {wildcards.experiment_label}.{wildcards.clip_replicate_label}"

rule check_window_concordance:
    input:
        windows = lambda wildcards: expand("output/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        "output/figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf",
        "output/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.tsv"
    params:
        error_file = "stderr/{experiment_label}.check_window_concordance.err",
        out_file = "stdout/{experiment_label}.check_window_concordance.out",
        run_time = "0:15:00",
        memory = "8000",
        job_name = "check_window_concordance"
    benchmark: "benchmarks/check_window_concordance/{experiment_label}.all_replicates.concordance.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=8000
    shell:
        "Rscript --vanilla {TOOL_DIR}/check_window_concordance.R output/tested_windows {wildcards.experiment_label} " + (BLACKLIST if BLACKLIST is not None else "") 

rule find_reproducible_enriched_windows:
    input:
        windows = lambda wildcards: expand("output/enriched_windows/{{experiment_label}}.{clip_replicate_label}.enriched_windows.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        linear_bar = "output/figures/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_window_counts.linear.pdf",
        log_bar = "output/figures/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_window_counts.log10.pdf"
    params:
        error_file = "stderr/{experiment_label}.find_reproducible_enriched_windows.err",
        out_file = "stdout/{experiment_label}.find_reproducible_enriched_windows.out",
        run_time = "00:30:00",
        memory = "2000",
        job_name = "find_reproducible_enriched_windows"
    benchmark: "benchmarks/find_reproducible_enriched_windows/{experiment_label}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=2000
    shell:
        "Rscript --vanilla {TOOL_DIR}/identify_reproducible_windows.R output/enriched_windows/ {wildcards.experiment_label} " + (BLACKLIST if BLACKLIST is not None else "") 

rule sample_background_windows_by_region:
    input:
        enriched_windows = rules.find_reproducible_enriched_windows.output.reproducible_windows,
        all_windows = ancient(FEATURE_ANNOTATIONS),
    output:
        variable_windows = "output/homer/region_matched_background/variable/{experiment_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz"
    params:
        error_file = "stderr/{experiment_label}.sample_background_windows_by_region.err",
        out_file = "stdout/{experiment_label}.sample_background_windows_by_region.out",
        run_time = "00:30:00",
        memory = "16000",
        job_name = "sample_background_windows"
    benchmark: "benchmarks/sample_background_windows_by_region/{experiment_label}.sample_background_windows_by_region.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=16000
    shell:
        "Rscript --vanilla {TOOL_DIR}/sample_matched_background_by_region.R {input.enriched_windows} {input.all_windows} 75 output/homer/region_matched_background {wildcards.experiment_label};"
