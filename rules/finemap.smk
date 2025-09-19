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
    resources:
        mem_mb=45000,
        runtime="2h"
    benchmark: "benchmarks/get_nt_coverage/{experiment_label}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    resources:
        mem_mb=45000
    shell:
        "zcat {input.windows} | tail -n +2 | sort -k1,1 -k2,2n | awk -v OFS=\"\t\" '{{start = $2-37; if(start < 0) {{start = 0}}; print $1, start, $3+37,$4,$5,$6}}' | "
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
    resources:
        mem_mb=45000,
        runtime="2h"
    benchmark: "benchmarks/finemap_windows/{experiment_label}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    resources:
        mem_mb=45000
    shell:
        "Rscript --vanilla {TOOL_DIR}/finemap_enriched_windows.R {input.nt_coverage} output/finemapping/mapped_sites/ {wildcards.experiment_label}"

rule annotate_finemap:
    input:
        finemapped = rules.finemap_windows.output.finemapped_windows,
        feature_annotations = FEATURE_ANNOTATIONS,
        ranking = ACCESSION_RANKINGS,
    output:
        "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv"
    threads:
        1
    resources:
        mem_mb=45000,
        runtime="1h"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        if [ -s {input.finemapped} ]; then
            python {TOOL_DIR}/annotate_finemapped_regions.py \
                {input.finemapped} \
                {input.ranking} \
                {input.feature_annotations} \
                {output}
        else
            touch {output}
        fi
        """

rule find_both_tested_windows:
    input:
        lambda wildcards: expand("output/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz",
         clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        tested_windows_in_2_rep = "output/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.bed",
        tested_windows_merged = "output/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.merged.bed"
    threads: 1
    resources:
        mem_mb=45000,
        runtime="1h"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/find_both_tested_windows.py \
            "{input}" \
            {output.tested_windows_in_2_rep} \
            {output.tested_windows_merged}
        """
