locals().update(config)
rule uniq_repeats:
    input:
        repeatmasker = ancient(REPEAT_TABLE),
        genome = ancient(GENOME)
    output:
        sorted_bed = temp("repeats.sort.temp.bed.gz"),
        unique_repeats = REPEAT_BED
    benchmark: "benchmarks/uniq_repeats.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/uniq_repeats.out",
        stderr = config["WORKDIR"] + "/stderr/uniq_repeats.err",
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting uniq_repeats" | tee -a {log.stdout}

        (zcat {REPEAT_TABLE} \
            | awk -v OFS="\t" '{{print $6,$7,$8,$11 ":" name_count[$11]++, $2, $10,$11,$12,$13}} 
                $13 == "L1" || $13 == "Alu" {{$11 = $11 "_AS"; $12 = $12 "_AS"; $13 = $13 "_AS"; 
                if($10 == "+") {{$10 = "-"}} else {{$10 = "+"}}; 
                print $6,$7,$8,$11 ":" name_count[$11]++, $2, $10,$11,$12,$13}}' \
            | tail -n +2 \
            | bedtools sort -i - \
            | gzip \
            > {output.sorted_bed}) \
        >> {log.stdout} 2> {log.stderr}

        (bedtools coverage -s -d -a {output.sorted_bed} -b {output.sorted_bed} \
            | awk -v OFS="\t" '$NF > 1 {{print $1,$2+$(NF-1)-1,$2+$(NF-1),$4,$5,$6}}' \
            | bedtools sort -i - \
            | bedtools merge -c 4,5,6 -o distinct -s -i - \
            | bedtools subtract -s -a {output.sorted_bed} -b - \
            | bedtools nuc -s -fi {input.genome} -bed - \
            | awk -v OFS="\t" 'NR > 1 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}}' \
            | gzip -c \
            > {output.unique_repeats}) \
        >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished uniq_repeats" | tee -a {log.stdout}
        """

rule quantify_repeats:
    input:
        chrom_sizes = config["CHROM_SIZES"],
        bam = lambda wildcards: replicate_label_to_bams[wildcards.replicate_label],
        repeats = REPEAT_BED,
    output:
        counts = "output/secondary_results/counts/repeats/vectors/{replicate_label}.counts"
    params:
        uninformative = config["UNINFORMATIVE_READ"]
    benchmark: "benchmarks/repeats/unassigned_experiment.{replicate_label}.quantify_repeats.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.quantify_repeats.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.quantify_repeats.err",
    conda:
        "envs/bedbam_tools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting quantify_repeats for {wildcards.replicate_label}" | tee -a {log.stdout}

        (bedtools bamtobed -i {input.bam} \
            | awk '($1 != "chrEBV") && ($4 !~ "/{params.uninformative}$")' \
            | bedtools flank -s -l 1 -r 0 -g {input.chrom_sizes} -i - \
            | bedtools shift -p 1 -m -1 -g {input.chrom_sizes} -i - \
            | bedtools sort -i - \
            | bedtools coverage -s -counts -a {input.repeats} -b - \
            | awk 'BEGIN {{print "{wildcards.replicate_label}"}} {{print $NF}}' \
            > {output.counts}) \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished quantify_repeats for {wildcards.replicate_label}" | tee -a {log.stdout}
        """

rule make_repeat_count_tables:
    input:
        unique_repeats = REPEAT_BED,
        replicate_counts = lambda wildcards: expand(
            "output/secondary_results/counts/repeats/vectors/{replicate_label}.counts", 
            replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]
        ),
    output:
        name_table = "output/secondary_results/counts/repeats/tables/name/{experiment_label}.tsv.gz",
        class_table = "output/secondary_results/counts/repeats/tables/class/{experiment_label}.tsv.gz",
        family_table = "output/secondary_results/counts/repeats/tables/family/{experiment_label}.tsv.gz",
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_repeat_count_table.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.make_repeat_count_tables.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.make_repeat_count_tables.err",
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting make_repeat_count_tables for {wildcards.experiment_label}" | tee -a {log.stdout}

        (echo "repeat_name" \
            | paste - {input.replicate_counts} \
            | sed -n '1p' \
            | gzip \
            > {output.name_table}) \
        >> {log.stdout} 2> {log.stderr}

        (echo "repeat_class" \
            | paste - {input.replicate_counts} \
            | sed -n '1p' \
            | gzip \
            > {output.class_table}) \
        >> {log.stdout} 2> {log.stderr}

        (echo "repeat_family" \
            | paste - {input.replicate_counts} \
            | sed -n '1p' \
            | gzip \
            > {output.family_table}) \
        >> {log.stdout} 2> {log.stderr}

        (paste <(
            zcat {input.unique_repeats} \
                | awk -v OFS="\t" 'BEGIN {{print "repeat_name";}} {{print $7}}'
        ) {input.replicate_counts} \
            | awk -v OFS="\t" 'NR > 1 {{for(i=2;i<=NF;i++) {{tabulation[$1][i]+=$i}}}} END {{for(name in tabulation) {{printf name; for(i=2;i<=NF;i++) {{printf "\t" tabulation[name][i]}} print ""}}}}' \
            | sort -k1,1 \
            | gzip \
            >> {output.name_table}) \
        >> {log.stdout} 2> {log.stderr}

        (paste <(
            zcat {input.unique_repeats} \
                | awk -v OFS="\t" 'BEGIN {{print "repeat_class";}} {{print $8}}'
        ) {input.replicate_counts} \
            | awk -v OFS="\t" 'NR > 1 {{for(i=2;i<=NF;i++) {{tabulation[$1][i]+=$i}}}} END {{for(name in tabulation) {{printf name; for(i=2;i<=NF;i++) {{printf "\t" tabulation[name][i]}} print ""}}}}' \
            | sort -k1,1 \
            | gzip \
            >> {output.class_table}) \
        >> {log.stdout} 2> {log.stderr}

        (paste <(
            zcat {input.unique_repeats} \
                | awk -v OFS="\t" 'BEGIN {{print "repeat_family";}} {{print $9}}'
        ) {input.replicate_counts} \
            | awk -v OFS="\t" 'NR > 1 {{for(i=2;i<=NF;i++) {{tabulation[$1][i]+=$i}}}} END {{for(name in tabulation) {{printf name; for(i=2;i<=NF;i++) {{printf "\t" tabulation[name][i]}} print ""}}}}' \
            | sort -k1,1 \
            | gzip \
            >> {output.family_table}) \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished make_repeat_count_tables for {wildcards.experiment_label}" | tee -a {log.stdout}
        """

rule fit_clip_betabinomial_re_model:
    input:
        table = rules.make_repeat_count_tables.output.name_table,
    output:
        coef = "output/secondary_results/clip_model_coef_re/{experiment_label}.{clip_replicate_label}.tsv",
    benchmark: "benchmarks/fit_clip_betabinomial_re_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_re_model.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.{clip_replicate_label}.fit_clip_betabinomial_re_model.err",
    conda:
        "envs/skipper_R.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 180 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting fit_clip_betabinomial_re_model for {wildcards.experiment_label}.{wildcards.clip_replicate_label}" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/fit_clip_betabinom_re.R \
            {input.table} \
            {wildcards.experiment_label} \
            {wildcards.clip_replicate_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished fit_clip_betabinomial_re_model for {wildcards.experiment_label}.{wildcards.clip_replicate_label}" | tee -a {log.stdout}
        """

rule fit_input_betabinomial_re_model:
    input:
        table = rules.make_repeat_count_tables.output.name_table,
    output:
        coef = "output/secondary_results/input_model_coef_re/{experiment_label}.{input_replicate_label}.tsv",
    benchmark: "benchmarks/fit_input_betabinomial_re_model/{experiment_label}.{input_replicate_label}.fit_input.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.{input_replicate_label}.fit_input_betabinomial_re_model.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.{input_replicate_label}.fit_input_betabinomial_re_model.err",
    conda:
        "envs/skipper_R.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting fit_input_betabinomial_re_model for {wildcards.experiment_label}.{wildcards.input_replicate_label}" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/fit_input_betabinom_re.R \
            {input.table} \
            {wildcards.experiment_label} \
            {wildcards.input_replicate_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished fit_input_betabinomial_re_model for {wildcards.experiment_label}.{wildcards.input_replicate_label}" | tee -a {log.stdout}
        """

rule call_enriched_re:
    input:
        table = rules.make_repeat_count_tables.output.name_table,
        replicate = lambda wildcards: "output/secondary_results/counts/repeats/vectors/" + re.sub(r"IP_\d$","IP_2",wildcards.clip_replicate_label) + ".counts",
        repeats = REPEAT_BED,
        parameters = lambda wildcards: "output/secondary_results/" + OVERDISPERSION_MODE + "_model_coef_re/{experiment_label}." + overdispersion_replicate_lookup[wildcards.clip_replicate_label] + ".tsv",
    threads: 2
    output:
        "output/figures/secondary_figures/clip_scatter_re/{experiment_label}.{clip_replicate_label}.clip_test_distribution.pdf",
        "output/secondary_results/enriched_re/{experiment_label}.{clip_replicate_label}.enriched_re.tsv.gz"
    benchmark: "benchmarks/call_enriched_re/{experiment_label}.{clip_replicate_label}.call_enriched_re.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.{clip_replicate_label}.call_enriched_re.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.{clip_replicate_label}.call_enriched_re.err",
    conda:
        "envs/skipper_R.yaml"
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label]
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 180 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting call_enriched_re for {wildcards.experiment_label}.{wildcards.clip_replicate_label}" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/call_enriched_re.R \
            {input.table} \
            {input.repeats} \
            {input.parameters} \
            {params.input_replicate_label} \
            {wildcards.clip_replicate_label} \
            {wildcards.experiment_label}.{wildcards.clip_replicate_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished call_enriched_re for {wildcards.experiment_label}.{wildcards.clip_replicate_label}" | tee -a {log.stdout}
        """

rule find_reproducible_enriched_re:
    input:
        windows = lambda wildcards: expand(
            "output/secondary_results/enriched_re/{{experiment_label}}.{clip_replicate_label}.enriched_re.tsv.gz",
            clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label]
        )
    output:
        reproducible_windows = "output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz",
    benchmark: "benchmarks/find_reproducible_enriched_re/{experiment_label}.all_replicates.reproducible.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.find_reproducible_enriched_re.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.find_reproducible_enriched_re.err",
    conda:
        "envs/skipper_R.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting find_reproducible_enriched_re for {wildcards.experiment_label}" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/identify_reproducible_re.R \
            output/secondary_results/enriched_re/ \
            {wildcards.experiment_label} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished find_reproducible_enriched_re for {wildcards.experiment_label}" | tee -a {log.stdout}
        """
        
