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
    params:
        error_file = "stderr/make_mega_table.err",
        out_file = "stdout/make_mega_table.out",
        run_time = "00:25:00",
        cores = "1",
        memory = "200",
        job_name = "make_genome_mega_table"
    shell:
        """
        paste <(zcat {input.feature}) {input.replicate_counts} | gzip > {output}
        """

rule make_repeat_mega_tables:
    input:
        unique_repeats = REPEAT_BED,
        replicate_counts = lambda wildcards: expand("output/counts/repeats/vectors/{replicate_label}.counts", 
        replicate_label = replicate_labels),
    output:
        name_table = "output/counts/repeats/megatables/name.tsv.gz",
        class_table = "output/counts/repeats/megatables/class.tsv.gz",
        family_table = "output/counts/repeats/megatables/family.tsv.gz",
    params:
        error_file = "stderr/make_repeat_mega_tables.err",
        out_file = "stdout/make_repeat_mega_tables.out",
        run_time = "06:00:00",
        cores = "1",
        memory = "2000",
        job_name = "make_repeat_mega_tables"
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

# Unique fragment per library
rule join_unique_fragments:
    input:
        expand("output/QC/{replicate_label}.uniq_fragments", replicate_label = replicate_labels)
    output:
        "output/QC/unique_fragments.csv"
    params:
        error_file = "stderr/join_fragment_counts.err",
        out_file = "stdout/join_fragment_counts.out",
        run_time = "05:00",
        cores = "1",
        memory = "2000",
        job_name = "join_unique_fragments"
    shell:
        """
        awk '{{print FILENAME "," $0}}' {input} > {output}
        """

# summarize per transcript type and family type
rule summarize_genome_megatable:
    input:
        "output/counts/genome/megatables/megatable.tsv.gz",
    output:
        f="output/counts/genome/megatables/feature_type_top.tsv.gz",
        t="output/counts/genome/megatables/transcript_type_top.tsv.gz",
    params:
        error_file = "stderr/summarize_genome_megatable.err",
        out_file = "stdout/summarize_genome_megatable.out",
        run_time = "20:00",
        cores = "1",
        memory = "2000",
        job_name = "summarize_genome_megatable"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/group_genome_megatable.py {input} {output.f} {output.t}
        """
# tested windows # with nan ready for analysis?
# table containing QC statistics

# prepare deep learning training
# fetch variants