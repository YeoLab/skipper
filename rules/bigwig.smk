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
    params:
        error_file = "stderr/{replicate_label}.make_bigwig.err",
        out_file = "stdout/{replicate_label}.make_bigwig.out",
        run_time = "4:00:00",
        memory = "10000",
        job_name = "make_bigwig"
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
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
    params:
        error_file = "stderr/{replicate_label}.make_bigwig.err",
        out_file = "stdout/{replicate_label}.make_bigwig.out",
        run_time = "40:00",
        memory = "1000",
        job_name = "make_bigwig"
    benchmark: "benchmarks/bigwigs/unassigned_experiment.{replicate_label}.make_bigwig.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
        "factor=$(samtools idxstats {input.bam} | cut -f 3 | paste -sd+ | bc | xargs -I {{}} echo 'scale=6; 10^6 / {{}}' | bc);"
        "bedtools genomecov -scale $factor -5 -strand + -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_plus};"
        "bedtools genomecov -scale $factor -5 -strand - -bg -ibam {input.bam} | sort -k1,1 -k2,2n | grep -v EBV > {output.bg_minus};"
        "bedGraphToBigWig {output.bg_plus} {CHROM_SIZES} {output.bw_plus};" 
        "bedGraphToBigWig {output.bg_minus} {CHROM_SIZES} {output.bw_minus};" 