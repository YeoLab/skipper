locals().update(config)
rule get_nt_coverage:
    input:
        windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz",
        clip_bams = lambda wildcards: [config['replicate_label_to_bams'][clip_replicate_label] for clip_replicate_label in experiment_to_clip_replicate_labels[wildcards.experiment_label]],
        input_bams = lambda wildcards: [config['replicate_label_to_bams'][input_replicate_label] for input_replicate_label in experiment_to_input_replicate_labels[wildcards.experiment_label]],         
    output:
        nt_census = temp("output/finemapping/nt_coverage/{experiment_label}.nt_census.bed"),
        nt_input_counts = temp("output/finemapping/nt_coverage/{experiment_label}.nt_coverage.input.counts"),
        nt_clip_counts = temp("output/finemapping/nt_coverage/{experiment_label}.nt_coverage.clip.counts"),
        nt_coverage = "output/finemapping/nt_coverage/{experiment_label}.nt_coverage.bed"
    threads: 6
    params:
        error_file = "stderr/{experiment_label}.get_nt_coverage.err",
        out_file = "stdout/{experiment_label}.get_nt_coverage.out",
        run_time = "6:00:00",
        memory = "45000",
        job_name = "get_nt_coverage"
    benchmark: "benchmarks/get_nt_coverage/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "set +eu;"
        "module load samtools/1.16 bedtools;"
        "zcat {input.windows} | tail -n +2 | sort -k1,1 -k2,2n | awk -v OFS=\"\t\" '{{print $1, $2 -37, $3+37,$4,$5,$6}}' | "
            "bedtools merge -i - -s -c 6 -o distinct | awk -v OFS=\"\t\" '{{for(i=$2;i< $3;i++) {{print $1,i,i+1,\"MW:\" NR \":\" i - $2,0,$4, NR}} }}' > {output.nt_census}; "
        "samtools cat {input.input_bams} | bedtools intersect -s -wa -a - -b {output.nt_census} | "
            "bedtools bamtobed -i - | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -counts -s -a {output.nt_census} -b - | awk '{{print $NF}}' > {output.nt_input_counts};"
        "samtools cat {input.clip_bams} | bedtools intersect -s -wa -a - -b {output.nt_census} | "
            "bedtools bamtobed -i - | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -counts -s -a {output.nt_census} -b - | awk '{{print $NF}}' > {output.nt_clip_counts};"
        "paste {output.nt_census} {output.nt_input_counts} {output.nt_clip_counts} > {output.nt_coverage}"

rule finemap_windows:
    input:
        nt_coverage = rules.get_nt_coverage.output.nt_coverage,        
    output:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    threads: 6,
    params:
        error_file = "stderr/{experiment_label}.finemap_windows.err",
        out_file = "stdout/{experiment_label}.finemap_windows.out",
        run_time = "6:00:00",
        memory = "60000",
        job_name = "finemap_windows"
    benchmark: "benchmarks/finemap_windows/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/finemap_enriched_windows.R {input.nt_coverage} output/finemapping/mapped_sites/ {wildcards.experiment_label}"
