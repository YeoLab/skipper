locals().update(config)

rule make_unscaled_bigwig:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/unscaled/plus/{replicate_label}.unscaled.plus.bg"),
        bg_minus = temp("output/bedgraphs/unscaled/minus/{replicate_label}.unscaled.minus.bg"),
        bw_plus = "output/bigwigs/unscaled/plus/{replicate_label}.unscaled.plus.bw",
        bw_minus = "output/bigwigs/unscaled/minus/{replicate_label}.unscaled.minus.bw",
    resources:
        mem_mb = 12000,
        tmpdir = TMPDIR,
        runtime = "1h"
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    log: "logs/{replicate_label}.make_unscaled_bigwig.log"
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting make_unscaled_bigwig" | tee "{log}"

        samtools index {input.bam} 2>&1 | tee -a "{log}"

        bedtools genomecov -split -strand + -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_plus} 2>&1 | tee -a "{log}"

        bedtools genomecov -split -strand - -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_minus} 2>&1 | tee -a "{log}"

        bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus} 2>&1 | tee -a "{log}"
        bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus} 2>&1 | tee -a "{log}"

        echo "[`date`] Finished make_unscaled_bigwig" | tee -a "{log}"
        """

rule make_scaled_bigwig:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/scaled/plus/{replicate_label}.scaled.plus.bg"),
        bg_minus = temp("output/bedgraphs/scaled/minus/{replicate_label}.scaled.minus.bg"),
        bw_plus = "output/bigwigs/scaled/plus/{replicate_label}.scaled.plus.bw",
        bw_minus = "output/bigwigs/scaled/minus/{replicate_label}.scaled.minus.bw",
    resources:
        mem_mb = 12000,
        tmpdir = TMPDIR,
        runtime = "1h"
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    log: "logs/{replicate_label}.make_scaled_bigwig.log"
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting make_scaled_bigwig" | tee "{log}"

        samtools index {input.bam} 2>&1 | tee -a "{log}"

        FACTOR=$(samtools idxstats {input.bam} \
            | cut -f 3 \
            | paste -sd+ \
            | bc \
            | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' \
            | bc)

        bedtools genomecov -scale $FACTOR -split -strand + -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_plus} \
            2>&1 | tee -a "{log}"

        bedtools genomecov -scale $FACTOR -split -strand - -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_minus} \
            2>&1 | tee -a "{log}"

        bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus} 2>&1 | tee -a "{log}"
        bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus} 2>&1 | tee -a "{log}"

        echo "[`date`] Finished make_scaled_bigwig" | tee -a "{log}"
        """

rule make_scaled_bigwig_coverage:
    input:
        CHROM_SIZES,
        bam = lambda wildcards: config['replicate_label_to_bams'][wildcards.replicate_label],
    output:
        bg_plus = temp("output/bedgraphs/scaled/plus/{replicate_label}.scaled.cov.plus.bg"),
        bg_minus = temp("output/bedgraphs/scaled/minus/{replicate_label}.scaled.cov.minus.bg"),
        bw_plus = "output/bigwigs/scaled/plus/{replicate_label}.scaled.cov.plus.bw",
        bw_minus = "output/bigwigs/scaled/minus/{replicate_label}.scaled.cov.minus.bw",
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    log: "logs/{replicate_label}.make_unscaled_bigwig_coverage.log"
    resources:
        mem_mb = 10000,
        tmpdir = TMPDIR,
        runtime = "1h"
    conda:
        "envs/bedbam_tools.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Running on node: $(hostname)" | tee -a {log}
        echo "[`date`] Starting make_scaled_bigwig_coverage" | tee "{log}"

        factor=$(samtools idxstats {input.bam} \
            | cut -f 3 \
            | paste -sd+ \
            | bc \
            | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' \
            | bc)

        bedtools genomecov -scale $factor -strand + -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_plus} 2>&1 | tee -a "{log}"

        bedtools genomecov -scale $factor -strand - -bg -ibam {input.bam} \
            | sort -k1,1 -k2,2n \
            | grep -v EBV \
            > {output.bg_minus} 2>&1 | tee -a "{log}"

        bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus} 2>&1 | tee -a "{log}"
        bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus} 2>&1 | tee -a "{log}"

        echo "[`date`] Finished make_scaled_bigwig_coverage" | tee -a "{log}"
        """