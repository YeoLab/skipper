locals().update(config)

rule get_nt_coverage:
    input:
        windows = "output/secondary_results/unfiltered_reproducible_enriched_windows/{experiment_label}.unfiltered_reproducible_enriched_windows.tsv.gz",
        clip_bams = lambda wildcards: [
            config['replicate_label_to_bams'][clip_replicate_label]
            for clip_replicate_label in experiment_to_clip_replicate_labels[wildcards.experiment_label]
        ],
        input_bams = lambda wildcards: [
            config['replicate_label_to_bams'][input_replicate_label]
            for input_replicate_label in experiment_to_input_replicate_labels[wildcards.experiment_label]
        ],         
    output:
        nt_census = temp("output/secondary_results/finemapping/nt_coverage/{experiment_label}.nt_census.bed"),
        nt_input_counts = temp("output/secondary_results/finemapping/nt_coverage/{experiment_label}.nt_coverage.input.counts"),
        nt_clip_counts = temp("output/secondary_results/finemapping/nt_coverage/{experiment_label}.nt_coverage.clip.counts"),
        nt_coverage = "output/secondary_results/finemapping/nt_coverage/{experiment_label}.nt_coverage.bed"
    params:
        chrom_sizes = config["CHROM_SIZES"],
        uninformative = config["UNINFORMATIVE_READ"]
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 45000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 120 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/get_nt_coverage/{experiment_label}.all_replicates.reproducible.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.get_nt_coverage.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.get_nt_coverage.err",
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting get_nt_coverage" | tee -a {log.stdout}

        zcat {input.windows} \
            | tail -n +2 \
            | sort -k1,1 -k2,2n \
            | awk -v OFS="\t" '{{start = $2-37; if(start < 0) {{start = 0}}; print $1, start, $3+37,$4,$5,$6}}' \
            | bedtools merge -i - -s -c 6 -o distinct \
            | awk -v OFS="\t" '{{for(i=$2;i<$3;i++) {{print $1,i,i+1,"MW:" NR ":" i - $2,0,$4, NR}} }}' \
            > {output.nt_census} 2> {log.stderr}

        samtools cat {input.input_bams} \
            | bedtools intersect -s -wa -a - -b {output.nt_census} \
            | bedtools bamtobed -i - \
            | awk '($1 != "chrEBV") && ($4 !~ "/{params.uninformative}$")' \
            | bedtools flank -s -l 1 -r 0 -g {params.chrom_sizes} -i - \
            | bedtools shift -p 1 -m -1 -g {params.chrom_sizes} -i - \
            | bedtools sort -i - \
            | bedtools coverage -counts -s -a {output.nt_census} -b - \
            | awk '{{print $NF}}' \
            > {output.nt_input_counts} 2>> {log.stderr}

        samtools cat {input.clip_bams} \
            | bedtools intersect -s -wa -a - -b {output.nt_census} \
            | bedtools bamtobed -i - \
            | awk '($1 != "chrEBV") && ($4 !~ "/{params.uninformative}$")' \
            | bedtools flank -s -l 1 -r 0 -g {params.chrom_sizes} -i - \
            | bedtools shift -p 1 -m -1 -g {params.chrom_sizes} -i - \
            | bedtools sort -i - \
            | bedtools coverage -counts -s -a {output.nt_census} -b - \
            | awk '{{print $NF}}' \
            > {output.nt_clip_counts} 2>> {log.stderr}

        paste {output.nt_census} {output.nt_input_counts} {output.nt_clip_counts} \
            > {output.nt_coverage} 2>> {log.stderr}

        echo "[`date`] Finished get_nt_coverage" >> {log.stdout}
        """

rule finemap_windows:
    input:
        nt_coverage = rules.get_nt_coverage.output.nt_coverage        
    output:
        finemapped_windows = "output/secondary_results/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 45000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 180 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/finemap_windows/{experiment_label}.all_replicates.reproducible.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.finemap_windows.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.finemap_windows.err",
    conda:
        "envs/skipper_R.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting finemap_windows" | tee -a {log.stdout}

        Rscript --vanilla {TOOL_DIR}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            output/secondary_results/finemapping/mapped_sites/ \
            {wildcards.experiment_label} \
        >> {log.stdout} 2> {log.stderr}
        gzip -t {output.finemapped_windows}
        echo "[`date`] Finished finemap_windows" >> {log.stdout}
        """

rule annotate_finemap:
    input:
        finemapped = rules.finemap_windows.output.finemapped_windows,
        feature_annotations = ancient(FEATURE_ANNOTATIONS),
        ranking = ACCESSION_RANKINGS
    output:
        "output/secondary_results/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv"
    threads: 1
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.annotate_finemap.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.annotate_finemap.err",
    resources:
        mem_mb=lambda wildcards, attempt: 45000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    conda:
        "envs/metadensity.yaml"
    envmodules:
        "finemap/1.100.0"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting annotate_finemap" | tee -a {log.stdout}

        PY="$CONDA_PREFIX/bin/python"

        echo "Using python: $($PY -c 'import sys; print(sys.executable)')" | tee -a {log.stdout}
        echo "PATH=$PATH" | tee -a {log.stdout}

        if [ -s {input.finemapped} ]; then
            "$PY" {TOOL_DIR}/annotate_finemapped_regions.py \
                {input.finemapped} \
                {input.ranking} \
                {input.feature_annotations} \
                {output} \
            >> {log.stdout} 2> {log.stderr}
        else
            : > {output}  # Create empty output file.
        fi

        echo "[`date`] Finished annotate_finemap" | tee -a {log.stdout}
        """

rule find_both_tested_windows:
    input:
        lambda wildcards: expand(
            "output/secondary_results/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz",
            clip_replicate_label=experiment_to_clip_replicate_labels[wildcards.experiment_label]
        )
    output:
        tested_windows_in_2_rep = "output/secondary_results/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.bed",
        tested_windows_merged = "output/secondary_results/finemapping/both_tested_sites/{experiment_label}.both_tested_windows.merged.bed"
    threads: 1
    log:
        stdout = config["WORKDIR"] + "/stdout/{experiment_label}.find_both_tested_windows.out",
        stderr = config["WORKDIR"] + "/stderr/{experiment_label}.find_both_tested_windows.err",
    resources:
        mem_mb=lambda wildcards, attempt: 45000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    conda:
        "envs/metadensity.yaml"
    envmodules:
        "finemap/1.100.0"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting find_both_tested_windows" | tee -a {log.stdout}

        PY="$CONDA_PREFIX/bin/python"

        echo "Using python: $($PY -c 'import sys; print(sys.executable)')" | tee -a {log.stdout}
        echo "PATH=$PATH" | tee -a {log.stdout}

        "$PY" {TOOL_DIR}/find_both_tested_windows.py \
            "{input}" \
            {output.tested_windows_in_2_rep} \
            {output.tested_windows_merged} \
        >> {log.stdout} 2> {log.stderr}

        echo "[`date`] Finished find_both_tested_windows" | tee -a {log.stdout}
        """
