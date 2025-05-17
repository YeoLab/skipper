rule ctk:
    input:
        bam_ip_umap = "output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam"
    output:
        mdtag = temp("output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.md.sam.gz"),
        mutation_file = "output/ctk/{replicate_label}.mutation.txt",
        tagbed = "output/ctk/{replicate_label}.tag.bed",
        peak = "output/ctk/{replicate_label}.uniq.peak.sig.bed",
        peak_bd = "output/ctk/{replicate_label}.uniq.peak.sig.boundary.bed",
        peak_PH = "output/ctk/{replicate_label}.uniq.peak.sig.halfPH.bed",
    params:
        error_out_file = "error_files/ctk.{replicate_label}.err",
        out_file = "stdout/ctk.{replicate_label}.out",
        run_time = "2:10:00",
        memory = 10000,
        cores = 1,
    conda:
        "envs/ctk.yaml"
    shell:
        """
        samtools fillmd {input.bam_ip_umap} {GENOMEFA} | gzip -c > {output.mdtag}
        parseAlignment.pl \
            -v --map-qual 1 \
            --min-len 18 \
            --mutation-file {output.mutation_file} \
            {output.mdtag} - > {output.tagbed}

        tag2peak.pl -big -ss \
            -v --valley-seeking -p 0.05 --valley-depth 0.9 \
            --multi-test --dbkey hg38 \
            {output.tagbed} \
            {output.peak} \
            --out-boundary {output.peak_bd} \
            --out-half-PH {output.peak_PH}
        """

rule mcross_get_kmer_seed:
    input:
        foreground = rules.fetch_sequence.output.finemapped_fa,
        background = rules.fetch_sequence.output.background_fa,
    output:
        kmer_enrichment = "output/ctk/mcross/{experiment_label}.kmer.txt",
        config = "output/ctk/mcross/{experiment_label}.config.txt",
        topn_kmer_matrix = "output/ctk/mcross/{experiment_label}.w7.zcore.mat.txt",
        top_peak = "output/ctk/mcross/top7mer/top.{experiment_label}.txt",
    params:
        error_file = "error_files/mcross_kmer.{experiment_label}.err",
        out_file = "stdout/mcross_kmer.{experiment_label}.out",
        run_time = "2:00:00",
        memory = 10000,
        cores = 1,
    conda:
        "envs/ctk.yaml"
    shell:
        """
        word_enrich.pl -w 7 \
            -test binom -v \
            {input.foreground} \
            {input.background} \
            {output.kmer_enrichment}
        
        # generate config
        echo '{output.kmer_enrichment}"\t\"{wildcards.experiment_label}' > {output.config}

        gen_word_enrich_matrix.pl  \
            {output.config}  {output.topn_kmer_matrix}

        topword.R {output.topn_kmer_matrix} {wildcards.experiment_label}_top7mer
        """