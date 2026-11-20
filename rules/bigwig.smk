locals().update(config)
# Read 5' ends (truncation / crosslink sites), as in vbb63a25-tscc and the Her 2026 RBPNet models; RBPNet prep_data.py
# trains on these. f8f7c84 had switched this to coverage (-split); see YeoLab/skipper#69. Scaled bigwigs stay coverage.
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
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        tmpdir=TMPDIR,
        runtime=lambda wildcards, attempt: 30 * attempt
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        tmpdir=TMPDIR,
        runtime=lambda wildcards, attempt: 60 * attempt
    shell:
        "samtools index {input.bam};"
        "bedtools genomecov -5 -strand + -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_plus};"
        "bedtools genomecov -5 -strand - -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_minus};"
        "bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus};" 
        "bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus};" 

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
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        tmpdir=TMPDIR,
        runtime=lambda wildcards, attempt: 60 * attempt
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
        "samtools index {input.bam};"
        "FACTOR=$(samtools idxstats {input.bam} | cut -f 3 | paste -sd+ | bc | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' | bc);"
        "bedtools genomecov -scale $FACTOR -split -strand + -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_plus};"
        "bedtools genomecov -scale $FACTOR -split -strand - -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_minus};"
        "bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus};"
        "bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus};"

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
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        tmpdir=TMPDIR,
        runtime=lambda wildcards, attempt: 60 * attempt
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
        "factor=$(samtools idxstats {input.bam} | cut -f 3 | paste -sd+ | bc | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' | bc);"
        "bedtools genomecov -scale $factor -strand + -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_plus};"
        "bedtools genomecov -scale $factor -strand - -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_minus};"
        "bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus};" 
        "bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus};" 