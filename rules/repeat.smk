locals().update(config)
rule uniq_repeats:
    input:
        repeatmasker = ancient(REPEAT_TABLE),
        genome = ancient(GENOME)
    output:
        sorted_bed = temp("repeats.sort.temp.bed.gz"),
        unique_repeats = REPEAT_BED
    benchmark: "benchmarks/uniq_repeats.txt"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    resources:
        mem_mb=32000,
        runtime="4h"
    shell:
        "zcat {REPEAT_TABLE} | awk -v OFS=\"\\t\" '{{print $6,$7,$8,$11 \":\" name_count[$11]++, $2, $10,$11,$12,$13}} "
            "$13 == \"L1\" || $13 == \"Alu\" {{$11 = $11 \"_AS\"; $12 = $12 \"_AS\"; $13 = $13 \"_AS\"; "
            "if($10 == \"+\") {{$10 = \"-\"}} else {{$10 = \"+\"}}; print $6,$7,$8,$11 \":\" name_count[$11]++, $2, $10,$11,$12,$13}}' | "
            "tail -n +2 | bedtools sort -i - | gzip > {output.sorted_bed}; "
        "bedtools coverage -s -d -a {output.sorted_bed} -b {output.sorted_bed}  | awk -v OFS=\"\\t\" "
            "'$NF >1 {{print $1,$2+$(NF-1)-1,$2+$(NF-1),$4,$5,$6}}' | "
            "bedtools sort -i - | "
            "bedtools merge -c 4,5,6 -o distinct -s -i - | "
            "bedtools subtract -s -a {output.sorted_bed} -b - | "
            "bedtools nuc -s -fi {input.genome} -bed -  | awk -v OFS=\"\\t\" 'NR > 1 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}}' | "
            "gzip -c > {output.unique_repeats}"


rule quantify_repeats:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: replicate_label_to_bams[wildcards.replicate_label],
        repeats = REPEAT_BED
    output:
        counts = "output/counts/repeats/vectors/{replicate_label}.counts"
    benchmark: "benchmarks/repeats/unassigned_experiment.{replicate_label}.quantify_repeats.txt"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 40 * (2 ** (attempt - 1)),
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -s -counts -a {input.repeats} -b - | "
            "awk 'BEGIN {{print \"{wildcards.replicate_label}\"}} {{print $NF}}' > {output.counts}"

rule make_repeat_count_tables:
    input:
        unique_repeats = REPEAT_BED,
        replicate_counts = lambda wildcards: expand("output/counts/repeats/vectors/{replicate_label}.counts", replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        name_table = "output/counts/repeats/tables/name/{experiment_label}.tsv.gz",
        class_table = "output/counts/repeats/tables/class/{experiment_label}.tsv.gz",
        family_table = "output/counts/repeats/tables/family/{experiment_label}.tsv.gz",
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_repeat_count_table.txt"
    resources:
        mem_mb=2000,
        runtime=120
    shell:
        "echo \"repeat_name\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.name_table};"
        "echo \"repeat_class\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.class_table};"
        "echo \"repeat_family\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.family_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_name\";}} {{print $7}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.name_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_class\";}} {{print $8}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.class_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_family\";}} {{print $9}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.family_table};"

rule fit_clip_betabinomial_re_model:
    input:
        table = rules.make_repeat_count_tables.output.name_table,
    output:
        coef = "output/clip_model_coef_re/{experiment_label}.{clip_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment_label}}.{{clip_replicate_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    benchmark: "benchmarks/fit_clip_betabinomial_re_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=32000,
        runtime="6h"
    shell:
        "Rscript --vanilla {TOOL_DIR}/fit_clip_betabinom_re.R {input.table} {wildcards.experiment_label} {wildcards.clip_replicate_label}"

rule fit_input_betabinomial_re_model:
    input:
        table = rules.make_repeat_count_tables.output.name_table,
    output:
        coef = "output/input_model_coef_re/{experiment_label}.{input_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/input_distributions/{{experiment_label}}.{{input_replicate_label}}.{other_label}.input_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    benchmark: "benchmarks/fit_input_betabinomial_re_model/{experiment_label}.{input_replicate_label}.fit_input.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=32000,
        runtime=60
    shell:
        "Rscript --vanilla {TOOL_DIR}/fit_input_betabinom_re.R {input.table} {wildcards.experiment_label} {wildcards.input_replicate_label}"

rule call_enriched_re:
    input:
        table = rules.make_repeat_count_tables.output.name_table,
        replicate = lambda wildcards: "output/counts/repeats/vectors/" + re.sub("IP_\d$","IP_2",wildcards.clip_replicate_label) + ".counts",
        repeats = REPEAT_BED,
        parameters = lambda wildcards: "output/" + OVERDISPERSION_MODE + "_model_coef_re/{experiment_label}." + overdispersion_replicate_lookup[wildcards.clip_replicate_label] + ".tsv",
    threads: 2
    output:
        "output/figures/clip_scatter_re/{experiment_label}.{clip_replicate_label}.clip_test_distribution.pdf",
        "output/enriched_re/{experiment_label}.{clip_replicate_label}.enriched_re.tsv.gz"
    benchmark: "benchmarks/call_enriched_re/{experiment_label}.{clip_replicate_label}.call_enriched_re.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label]
    resources:
        mem_mb=48000,
        runtime="3h"
    shell:
        "Rscript --vanilla {TOOL_DIR}/call_enriched_re.R {input.table} {input.repeats} {input.parameters} {params.input_replicate_label} {wildcards.clip_replicate_label} {wildcards.experiment_label}.{wildcards.clip_replicate_label}"

rule find_reproducible_enriched_re:
    input:
        windows = lambda wildcards: expand("output/enriched_re/{{experiment_label}}.{clip_replicate_label}.enriched_re.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        reproducible_windows = "output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz",
    benchmark: "benchmarks/find_reproducible_enriched_re/{experiment_label}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=16000,
        runtime=120
    shell:
        "Rscript --vanilla {TOOL_DIR}/identify_reproducible_re.R output/enriched_re/ {wildcards.experiment_label}"
        
