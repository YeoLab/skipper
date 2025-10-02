locals().update(config)

rule get_nt_coverage:
    input:
        windows = "output/unfiltered_reproducible_enriched_windows/{experiment_label}.unfiltered_reproducible_enriched_windows.tsv.gz",
        clip_bams = lambda wildcards: [
            config['replicate_label_to_bams'][clip_replicate_label]
            for clip_replicate_label in experiment_to_clip_replicate_labels[wildcards.experiment_label]
        ],
        input_bams = lambda wildcards: [
            config['replicate_label_to_bams'][input_replicate_label]
            for input_replicate_label in experiment_to_input_replicate_labels[wildcards.experiment_label]
        ],         
    output:
        nt_census = temp("output/finemapping/nt_coverage/{experiment_label}.nt_census.bed"),
        nt_input_counts = temp("output/finemapping/nt_coverage/{experiment_label}.nt_coverage.input.counts"),
        nt_clip_counts = temp("output/finemapping/nt_coverage/{experiment_label}.nt_coverage.clip.counts"),
        nt_coverage = "output/finemapping/nt_coverage/{experiment_label}.nt_coverage.bed"
    params:
        chrom_sizes = config["CHROM_SIZES"],
        uninformative = config["UNINFORMATIVE_READ"]
    threads: 6
    resources:
        mem_mb = 45000,
        runtime = "2h"
    benchmark: "benchmarks/get_nt_coverage/{experiment_label}.all_replicates.reproducible.txt"
    log: "logs/{experiment_label}.get_nt_coverage.log"
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting get_nt_coverage" 2>&1 | tee {log}

        zcat {input.windows} \
            | tail -n +2 \
            | sort -k1,1 -k2,2n \
            | awk -v OFS="\t" '{{start = $2-37; if(start < 0) {{start = 0}}; print $1, start, $3+37,$4,$5,$6}}' \
            | bedtools merge -i - -s -c 6 -o distinct \
            | awk -v OFS="\t" '{{for(i=$2;i<$3;i++) {{print $1,i,i+1,"MW:" NR ":" i - $2,0,$4, NR}} }}' \
            > {output.nt_census} 2>&1 | tee -a {log}

        samtools cat {input.input_bams} \
            | bedtools intersect -s -wa -a - -b {output.nt_census} \
            | bedtools bamtobed -i - \
            | awk '($1 != "chrEBV") && ($4 !~ "/{params.uninformative}$")' \
            | bedtools flank -s -l 1 -r 0 -g {params.chrom_sizes} -i - \
            | bedtools shift -p 1 -m -1 -g {params.chrom_sizes} -i - \
            | bedtools sort -i - \
            | bedtools coverage -counts -s -a {output.nt_census} -b - \
            | awk '{{print $NF}}' \
            > {output.nt_input_counts} 2>&1 | tee -a {log}

        samtools cat {input.clip_bams} \
            | bedtools intersect -s -wa -a - -b {output.nt_census} \
            | bedtools bamtobed -i - \
            | awk '($1 != "chrEBV") && ($4 !~ "/{params.uninformative}$")' \
            | bedtools flank -s -l 1 -r 0 -g {params.chrom_sizes} -i - \
            | bedtools shift -p 1 -m -1 -g {params.chrom_sizes} -i - \
            | bedtools sort -i - \
            | bedtools coverage -counts -s -a {output.nt_census} -b - \
            | awk '{{print $NF}}' \
            > {output.nt_clip_counts} 2>&1 | tee -a {log}

        paste {output.nt_census} {output.nt_input_counts} {output.nt_clip_counts} \
            > {output.nt_coverage} 2>&1 | tee -a {log}

        echo "[`date`] Finished get_nt_coverage" 2>&1 | tee -a {log}
        """

rule finemap_windows:
    input:
        nt_coverage = rules.get_nt_coverage.output.nt_coverage        
    output:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    threads: 6
    resources:
        mem_mb = 45000,
        runtime = "2h"
    benchmark: "benchmarks/finemap_windows/{experiment_label}.all_replicates.reproducible.txt"
    log: "logs/{experiment_label}.finemap_windows.log"
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting finemap_windows" 2>&1 | tee {log}

        Rscript --vanilla {TOOL_DIR}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            output/finemapping/mapped_sites/ \
            {wildcards.experiment_label} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished finemap_windows" 2>&1 | tee -a {log}
        """

rule annotate_finemap:
    input:
        finemapped = rules.finemap_windows.output.finemapped_windows,
        feature_annotations = FEATURE_ANNOTATIONS,
        ranking = ACCESSION_RANKINGS
    output:
        "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv"
    threads: 1
    log: "logs/{experiment_label}.annotate_finemap.log"
    resources:
        mem_mb = 45000,
        runtime = "1h"
    conda:
        "envs/metadensity.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting annotate_finemap" 2>&1 | tee {log}

        if [ -s {input.finemapped} ]; then
            python {TOOL_DIR}/annotate_finemapped_regions.py \
                {input.finemapped} \
                {input.ranking} \
                {input.feature_annotations} \
                {output} \
                2>&1 | tee -a {log}
        else
            touch {output} 2>&1 | tee -a {log}
        fi

        echo "[`date`] Finished annotate_finemap" 2>&1 | tee -a {log}
        """

rule find_both_tested_windows:
    input:
        lambda wildcards: expand(
            "output/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz",
            clip_replicate_label=experiment_to_clip_replicate_labels[wildcards.experiment_label]
        )
    output:
        tested_windows_in_2_rep = "output/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.bed",
        tested_windows_merged = "output/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.merged.bed"
    threads: 1
    log: "logs/{experiment_label}.find_both_tested_windows.log"
    resources:
        mem_mb = 45000,
        runtime = "1h"
    conda:
        "envs/metadensity.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting find_both_tested_windows" 2>&1 | tee {log}

        python {TOOL_DIR}/find_both_tested_windows.py \
            "{input}" \
            {output.tested_windows_in_2_rep} \
            {output.tested_windows_merged} \
            2>&1 | tee -a {log}

        echo "[`date`] Finished find_both_tested_windows" 2>&1 | tee -a {log}
        """
