"""
snakemake -kps /tscc/nfs/home/s5xu/projects/skipper/rules/ctk_mcross.py -j 30 -w 30 --use-singularity --singularity-prefix /tscc/lustre/ddn/scratch/s5xu/singularity --singularity-args "--bind /tscc"
"""

OUTPUT = "/tscc/nfs/home/s5xu/scratch/skipper_output"
replicate_label = "RBFOX2_K562_ENCSR756CKJ_IP_1"
experiment_label = "RBFOX2_K562_ENCSR756CKJ"
INFORMATIVE_READ = 2
UNINFORMATIVE_READ = 1

rule all:
    input:
        finemapped_fa = f"{OUTPUT}/ml/sequence/{experiment_label}.foreground.fa",
        background_fa = f"{OUTPUT}/ml/sequence/{experiment_label}.background.fa",
        # mutation_file = f"{OUTPUT}/ctk/{replicate_label}.mutation.txt",
        # tagbed = f"{OUTPUT}/ctk/{replicate_label}.tag.bed",
        # peak = f"{OUTPUT}/ctk/{replicate_label}.uniq.peak.sig.bed",
        # peak_bd = f"{OUTPUT}/ctk/{replicate_label}.uniq.peak.sig.boundary.bed",
        # peak_PH = f"{OUTPUT}/ctk/{replicate_label}.uniq.peak.sig.halfPH.bed",
        kmer_enrichment = f"{OUTPUT}/ctk/mcross/{experiment_label}.kmer.txt",
        config = f"{OUTPUT}/ctk/mcross/{experiment_label}.config.txt",
        topn_kmer_matrix = f"{OUTPUT}/ctk/mcross/{experiment_label}.w7.zcore.mat.txt",
        top_peak = f"{OUTPUT}/ctk/mcross/top7mer/top.{experiment_label}.txt",


rule select_informative_read:
    input:
        bam_combined=f"{OUTPUT}/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.bam"
    output:
        bam_informative=f"{OUTPUT}/bams/genome_R1/{replicate_label}.genome.Aligned.sort.dedup.R2.bam"
    params:
        error_file = f"{OUTPUT}/stderr/{replicate_label}.select_informative_read.err",
        out_file = f"{OUTPUT}/stdout/{replicate_label}.select_informative_read.out",
        run_time = "00:30:00",
        memory = "10000",
        job_name = "select_informative_read",
    benchmark: f"{OUTPUT}/benchmarks/select/unassigned_experiment.{replicate_label}.select_informative_read.txt"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    resources:
        mem_mb=10000
    shell:
        "samtools view -bF 128 {input.bam_combined} > {output.bam_informative}"
        

rule uniquely_mapped_reads:
    input:
        bam = rules.select_informative_read.output.bam_informative
    output:
        bam_umap = f"{OUTPUT}/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam",
        bai_umap = f"{OUTPUT}/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam.bai",
    threads: 1
    params:
        error_out_file = f"{OUTPUT}/stderr/uniquemap.{replicate_label}.err",
        out_file = f"{OUTPUT}/stdout/uniquemap.{replicate_label}.out",
        run_time = "0:30:00",
        memory = 40000,
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        module load bamtools;
        bamtools filter -in {input.bam} -out {output.bam_umap} -mapQuality ">3"
        samtools index {output.bam_umap}
        """
        
        
rule fetch_sequence:
    input:
        finemapped_windows = f"{OUTPUT}/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        background = f"{OUTPUT}/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz",
    output:
        finemapped_fa = f"{OUTPUT}/ml/sequence/{experiment_label}.foreground.fa",
        background_fa = f"{OUTPUT}/ml/sequence/{experiment_label}.background.fa"
    params:
        error_file = f"{OUTPUT}/stderr/{experiment_label}.fetch_sequence.err",
        out_file = f"{OUTPUT}/stdout/{experiment_label}.fetch_sequence.out",
        run_time = "40:00",
        memory = "2000",
        job_name = "run_homer",
        fa = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    resources:
        mem_mb=2000
    shell:
        '''
        bedtools getfasta -fo {output.finemapped_fa} -fi {params.fa} -bed {input.finemapped_windows} -s
        bedtools getfasta -fo {output.background_fa} -fi {params.fa} -bed {input.background} -s
        '''
        
        
rule ctk:
    input:
        bam_ip_umap = f"{OUTPUT}/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam"
    output:
        mdtag = temp(f"{OUTPUT}/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.md.sam.gz"),
        mutation_file = f"{OUTPUT}/ctk/{replicate_label}.mutation.txt",
        tagbed = f"{OUTPUT}/ctk/{replicate_label}.tag.bed",
        peak = f"{OUTPUT}/ctk/{replicate_label}.uniq.peak.sig.bed",
        peak_bd = f"{OUTPUT}/ctk/{replicate_label}.uniq.peak.sig.boundary.bed",
        peak_PH = f"{OUTPUT}/ctk/{replicate_label}.uniq.peak.sig.halfPH.bed",
    params:
        fa = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        error_out_file = f"{OUTPUT}/error_files/ctk.{replicate_label}.err",
        out_file = f"{OUTPUT}/stdout/ctk.{replicate_label}.out",
        run_time = "2:00:00",
        memory = 10000,
        cores = 1,
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        module load ctk;
        samtools fillmd {input.bam_ip_umap} {params.fa} | gzip -c > {output.mdtag}
        /tscc/projects/ps-yeolab4/software/ctk/1.0.8/bin/ctk-1.0.8/parseAlignment.pl \
            -v --map-qual 1 \
            --min-len 18 \
            --mutation-file {output.mutation_file} \
            {output.mdtag} - > {output.tagbed}

        /tscc/projects/ps-yeolab4/software/ctk/1.0.8/bin/ctk-1.0.8/tag2peak.pl -big -ss \
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
        kmer_enrichment = f"{OUTPUT}/ctk/mcross/{experiment_label}.kmer.txt",
        config = f"{OUTPUT}/ctk/mcross/{experiment_label}.config.txt",
        topn_kmer_matrix = f"{OUTPUT}/ctk/mcross/{experiment_label}.w7.zcore.mat.txt",
        top_peak = f"{OUTPUT}/ctk/mcross/top7mer/top.{experiment_label}.txt",
    params:
        error_file = f"{OUTPUT}/error_files/mcross_kmer.{experiment_label}.err",
        out_file = f"{OUTPUT}/stdout/mcross_kmer.{experiment_label}.out",
        run_time = "2:00:00",
        memory = 10000,
        cores = 1,
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        module load ctk;
        /tscc/nfs/home/s5xu/projects/mCross/word_enrich.pl -w 7 \
            -test binom -v \
            {input.foreground} \
            {input.background} \
            {output.kmer_enrichment}
        
        # generate config
        echo '{output.kmer_enrichment}"\t\"{experiment_label}' > {output.config}

        /tscc/nfs/home/s5xu/projects/mCrossgen_word_enrich_matrix.pl  \
            {output.config}  {output.topn_kmer_matrix}

        /tscc/nfs/home/s5xu/projects/mCross/topword.R {output.topn_kmer_matrix} {experiment_label}_top7mer
        """