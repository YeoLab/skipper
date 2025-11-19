locals().update(config)
rule make_genome_mega_table:
    input:
        feature = FEATURE_ANNOTATIONS,
        replicate_counts = lambda wildcards: expand(
            "output/counts/genome/vectors/{replicate_label}.counts", 
            replicate_label = replicate_labels),
    output:
        "output/counts/genome/megatables/megatable.tsv.gz",
    threads: 4
    log:
        stdout = config["WORKDIR"] + "/stdout/make_genome_mega_table.out",
        stderr = config["WORKDIR"] + "/stderr/make_genome_mega_table.err",
    resources:
        mem_mb=8000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log.stdout}
        echo "[`date`] Starting make_genome_mega_table" | tee {log.stdout}

        (
        paste <(zcat {input.feature}) {input.replicate_counts} \
            | gzip > {output} 
        ) >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished make_genome_mega_table" | tee -a {log.stdout}
        """

rule make_repeat_mega_tables:
    input:
        unique_repeats = REPEAT_BED,
        replicate_counts = lambda wildcards: expand(
            "output/counts/repeats/vectors/{replicate_label}.counts", 
            replicate_label = replicate_labels),
    output:
        name_table = "output/counts/repeats/megatables/name.tsv.gz",
        class_table = "output/counts/repeats/megatables/class.tsv.gz",
        family_table = "output/counts/repeats/megatables/family.tsv.gz",
    log:
        stdout = config["WORKDIR"] + "/stdout/make_repeat_mega_tables.out",
        stderr = config["WORKDIR"] + "/stderr/make_repeat_mega_tables.err",
    resources:
        mem_mb=8000,
        runtime="2h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting make_repeat_mega_tables" | tee -a {log.stdout}

        (
          echo "repeat_name" \
          | paste - {input.replicate_counts} \
          | sed -n '1p' \
          | gzip > {output.name_table}
        ) >> {log.stdout} 2> {log.stderr}

        (
        echo "repeat_class" \
            | paste - {input.replicate_counts} \
            | sed -n '1p' \
            | gzip > {output.class_table} 
        ) >> {log.stdout} 2>> {log.stderr}

        (
        echo "repeat_family" \
            | paste - {input.replicate_counts} \
            | sed -n '1p' \
            | gzip > {output.family_table} 
        ) >> {log.stdout} 2>> {log.stderr}

        {
        paste <(zcat {input.unique_repeats} \
            | awk -v OFS="\t" 'BEGIN {{print "repeat_name";}} {{print $7}}') \
            {input.replicate_counts} \
            | awk -v OFS="\t" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf "\t" tabulation[name][i]}} print ""}} }}' \
            | sort -k 1,1 \
            | gzip >> {output.name_table} 
        ) >> {log.stdout} 2>> {log.stderr}

        (
        paste <(zcat {input.unique_repeats} \
            | awk -v OFS="\t" 'BEGIN {{print "repeat_class";}} {{print $8}}') \
            {input.replicate_counts} \
            | awk -v OFS="\t" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf "\t" tabulation[name][i]}} print ""}} }}' \
            | sort -k 1,1 \
            | gzip >> {output.class_table}
        ) >> {log.stdout} 2>> {log.stderr}

        (
        paste <(zcat {input.unique_repeats} \
            | awk -v OFS="\t" 'BEGIN {{print "repeat_family";}} {{print $9}}') \
            {input.replicate_counts} \
            | awk -v OFS="\t" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf "\t" tabulation[name][i]}} print ""}} }}' \
            | sort -k 1,1 \
            | gzip >> {output.family_table} 
        ) >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished make_repeat_mega_tables" | tee -a {log.stdout}
        """

rule join_aligned_reads:
    input:
        expand("output/QC/{replicate_label}.aligned_reads", replicate_label = replicate_labels)
    output:
        "output/QC/aligned_reads.csv"
    log:
        stdout = config["WORKDIR"] + "/stdout/join_aligned_reads.out",
        stderr = config["WORKDIR"] + "/stderr/join_aligned_reads.err",
    resources:
        mem_mb=2000,
        runtime="1h"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting join_aligned_reads" | tee -a {log.stdout}

        (awk '{{print FILENAME "," $0}}' {input} > {output})  >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished join_aligned_reads" | tee -a {log.stdout}
        """

# summarize per transcript type and family type
rule summarize_genome_megatable:
    input:
        "output/counts/genome/megatables/megatable.tsv.gz",
    output:
        f="output/counts/genome/megatables/feature_type_top.tsv.gz",
        t="output/counts/genome/megatables/transcript_type_top.tsv.gz",
    log:
        stdout = config["WORKDIR"] + "/stdout/summarize_genome_megatable.out",
        stderr = config["WORKDIR"] + "/stderr/summarize_genome_megatable.err",
    conda:
        "envs/metadensity.yaml"
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting summarize_genome_megatable" | tee -a {log.stdout}

        python {TOOL_DIR}/group_genome_megatable.py \
            {input} \
            {output.f} \
            {output.t} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished summarize_genome_megatable" | tee -a {log.stdout}
        """

rule join_reproducible_enriched_re:
    input:
        expand(
            "output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz", 
            experiment_label=experiment_labels
        )
    output:
        binary="output/joined_reproducible_re/binary.csv",
        l2or="output/joined_reproducible_re/l2or.csv",
    resources:
        mem_mb=2000,
        runtime="1h"
    log:
        stdout = config["WORKDIR"] + "/stdout/join_reproducible_enriched_re.out",
        stderr = config["WORKDIR"] + "/stderr/join_reproducible_enriched_re.err",
    conda:
        "envs/metadensity.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting join_reproducible_enriched_re" | tee -a {log.stdout}

        python {TOOL_DIR}/join_reproducible_enriched_re.py \
            . \
            {output.binary} \
            {output.l2or} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished join_reproducible_enriched_re" | tee -a {log.stdout}
        """

rule join_reproducible_enriched_windows_binary:
    input:
        expand(
            "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
            experiment_label=experiment_labels
        )
    output:
        binary="output/joined_reproducible_windows/binary.{feature_type}.csv",
        l2or="output/joined_reproducible_windows/l2or.{feature_type}.csv",
    resources:
        mem_mb=2000,
        runtime="1h"
    log:
        stdout = config["WORKDIR"] + "/stdout/{feature_type}.join_reproducible_enriched_windows_binary.out",
        stderr = config["WORKDIR"] + "/stderr/{feature_type}.join_reproducible_enriched_windows_binary.err",
    conda:
        "envs/metadensity.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting join_reproducible_enriched_windows" | tee -a {log.stdout}

        python {TOOL_DIR}/join_reproducible_enriched_windows.py \
            . \
            {wildcards.feature_type} \
            {output.binary} \
            {output.l2or} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished join_reproducible_enriched_windows" | tee -a {log.stdout}
        """

def find_all_tested_windows(experiment_labels, experiment_to_clip_replicate_labels):
    tested_windows = []
    for e in experiment_labels:
        for rep in experiment_to_clip_replicate_labels[e]:
            tested_windows.append(f"output/tested_windows/{e}.{rep}.tested_windows.tsv.gz")
    return tested_windows

rule join_reproducible_enriched_windows:
    input:
        expand(
            "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
            experiment_label=experiment_labels
        ),
        tested_windows = lambda wildcards: find_all_tested_windows(
            experiment_labels, experiment_to_clip_replicate_labels
        )
    output:
        binary="output/joined_reproducible_windows/all.csv",
    resources:
        mem_mb=2000,
        runtime="1h"
    log:
        stdout = config["WORKDIR"] + "/stdout/join_reproducible_enriched_windows.out",
        stderr = config["WORKDIR"] + "/stderr/join_reproducible_enriched_windows.err",
    conda:
        "envs/metadensity.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting join_reproducible_enriched_windows" | tee -a {log.stdout}

        python {TOOL_DIR}/join_all_reproducible_enriched_windows.py \
            . \
            {output.binary} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished join_reproducible_enriched_windows" | tee -a {log.stdout}
        """
