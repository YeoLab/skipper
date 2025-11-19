locals().update(config)

rule make_unscaled_bigwig:
    input:
        chrom_sizes = config["CHROM_SIZES"],
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/unscaled/plus/{replicate_label}.unscaled.plus.bg"),
        bg_minus = temp("output/bedgraphs/unscaled/minus/{replicate_label}.unscaled.minus.bg"),
        bw_plus = "output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw",
        bw_minus = "output/bigwigs/unscaled/minus/{replicate_label}.unscaled.minus.bw",
    resources:
        tmpdir = TMPDIR,
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.make_unscaled_bigwig.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.make_unscaled_bigwig.err",
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail
        
        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[$(date)] Starting make_unscaled_bigwig" | tee -a {log.stdout}
        
        samtools index {input.bam} >> {log.stdout} 2> {log.stderr}
        
        bedtools genomecov -split -strand + -bg -ibam {input.bam} \
          | sort -k1,1 -k2,2n \
          | grep -v EBV \
          > {output.bg_plus} 2>> {log.stderr}
        
        bedtools genomecov -split -strand - -bg -ibam {input.bam} \
          | sort -k1,1 -k2,2n \
          | grep -v EBV \
          > {output.bg_minus} 2>> {log.stderr}
        
        bedGraphToBigWig {output.bg_plus}  {input.chrom_sizes} {output.bw_plus}  >> {log.stdout} 2>> {log.stderr}
        bedGraphToBigWig {output.bg_minus} {input.chrom_sizes} {output.bw_minus} >> {log.stdout} 2>> {log.stderr}
        
        echo "[$(date)] Finished make_unscaled_bigwig" | tee -a {log.stdout}
        """

rule make_scaled_bigwig:
    input:
        chrom_sizes = config["CHROM_SIZES"],
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/scaled/plus/{replicate_label}.scaled.plus.bg"),
        bg_minus = temp("output/bedgraphs/scaled/minus/{replicate_label}.scaled.minus.bg"),
        bw_plus = "output/bigwigs/scaled/plus/{replicate_label}.scaled.plus.bw",
        bw_minus = "output/bigwigs/scaled/minus/{replicate_label}.scaled.minus.bw",
    resources:
        tmpdir = TMPDIR,
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.make_scaled_bigwig.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.make_scaled_bigwig.err",
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting make_scaled_bigwig" | tee -a {log.stdout}

        samtools index {input.bam} 2>&1 | tee -a {log.stdout}

        FACTOR=$(samtools idxstats {input.bam} \
            | cut -f 3 \
            | paste -sd+ \
            | bc \
            | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' \
            | bc)

        bedtools genomecov -scale $FACTOR -split -strand + -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_plus} 2> {log.stderr}

        bedtools genomecov -scale $FACTOR -split -strand - -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_minus} 2>> {log.stderr}

        bedGraphToBigWig {output.bg_plus} {input.chrom_sizes} {output.bw_plus} >> {log.stdout} 2>> {log.stderr}
        bedGraphToBigWig {output.bg_minus} {input.chrom_sizes} {output.bw_minus} >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished make_scaled_bigwig" | tee -a {log.stdout}
        """

rule make_scaled_bigwig_coverage:
    input:
        chrom_sizes = config["CHROM_SIZES"],
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/scaled/plus/{replicate_label}.scaled.cov.plus.bg"),
        bg_minus = temp("output/bedgraphs/scaled/minus/{replicate_label}.scaled.cov.minus.bg"),
        bw_plus = "output/bigwigs/scaled/plus/{replicate_label}.scaled.cov.plus.bw",
        bw_minus = "output/bigwigs/scaled/minus/{replicate_label}.scaled.cov.minus.bw",
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    log:
        stdout = config["WORKDIR"] + "/stdout/{replicate_label}.make_scaled_bigwig_coverage.out",
        stderr = config["WORKDIR"] + "/stderr/{replicate_label}.make_scaled_bigwig_coverage.err",
    resources:
        tmpdir = TMPDIR,
        mem_mb=lambda wildcards, attempt: 16000 * (1.5 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 60 * (2 ** (attempt - 1)),
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee {log.stdout}
        echo "[`date`] Starting make_scaled_bigwig_coverage" | tee -a {log.stdout}

        factor=$(samtools idxstats {input.bam} \
            | cut -f 3 \
            | paste -sd+ \
            | bc \
            | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' \
            | bc)

        bedtools genomecov -scale $factor -strand + -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
        > {output.bg_plus} 2> {log.stderr}

        bedtools genomecov -scale $factor -strand - -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
        > {output.bg_minus} 2>> {log.stderr}

        bedGraphToBigWig {output.bg_plus} {input.chrom_sizes} {output.bw_plus} >> {log.stdout} 2>> {log.stderr}
        bedGraphToBigWig {output.bg_minus} {input.chrom_sizes} {output.bw_minus} >> {log.stdout} 2>> {log.stderr}

        echo "[`date`] Finished make_scaled_bigwig_coverage" | tee -a {log.stdout}
        """