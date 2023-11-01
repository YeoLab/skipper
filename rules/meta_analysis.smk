locals().update(config)
rule make_genome_mega_table:
    input:
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
        "paste -d '\t' {input} > {output}"

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
