"""
mamba install -c bioconda perl-math-cdf perl-bio-featureio perl-env
export PERL5LIB=/tscc/nfs/home/s5xu/projects/czplib-v.1.1.1
export PATH=$PATH:/tscc/nfs/home/s5xu/projects/patternmatch

snakemake -kps /tscc/nfs/home/s5xu/projects/skipper/rules/ctk_mcross.py -j 30 -w 30 --use-singularity --singularity-prefix /tscc/lustre/ddn/scratch/s5xu/singularity --singularity-args "--bind /tscc" --rerun-incomplete
"""

OUTPUT = "/tscc/nfs/home/s5xu/scratch/skipper_output"
replicate_labels = ["RBFOX2_K562_ENCSR756CKJ_IP_1"]
experiment_labels = ["RBFOX2_K562_ENCSR756CKJ"]
INFORMATIVE_READ = 1
UNINFORMATIVE_READ = 2

mutation_types = ['del', 'ins', 'sub']
data_types = ['CIMS', 'CITS']

rule all:
    input:
        bam_informative = expand(os.path.join(OUTPUT, "bams/genome_R1/{replicate_label}.genome.Aligned.sort.dedup.R1.bam"), replicate_label=replicate_labels),
        bam_umap = expand(os.path.join(OUTPUT, "bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam"), replicate_label=replicate_labels),
        bai_umap = expand(os.path.join(OUTPUT, "bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam.bai"), replicate_label=replicate_labels),
        mdtag = expand(os.path.join(OUTPUT, "bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.md.sam"), replicate_label=replicate_labels),
        mutation_file = expand(os.path.join(OUTPUT, "ctk/{replicate_label}.mutation.txt"), replicate_label=replicate_labels),
        tagbed = expand(os.path.join(OUTPUT, "ctk/{replicate_label}.tag.bed"), replicate_label=replicate_labels),
        peak = expand(os.path.join(OUTPUT, "ctk/{replicate_label}.uniq.peak.sig.bed"), replicate_label=replicate_labels),
        peak_bd = expand(os.path.join(OUTPUT, "ctk/{replicate_label}.uniq.peak.sig.boundary.bed"), replicate_label=replicate_labels),
        peak_PH = expand(os.path.join(OUTPUT, "ctk/{replicate_label}.uniq.peak.sig.halfPH.bed"), replicate_label=replicate_labels),
        mutation_beds = expand(os.path.join(OUTPUT, "ctk", "{replicate_label}.uniq.{mutation_type}.bed"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CIMSs = expand(os.path.join(OUTPUT, "ctk/CIMS/{replicate_label}.uniq.{mutation_type}.CIMS.txt"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CIMS_stats = expand(os.path.join(OUTPUT, "ctk/CIMS/{replicate_label}.uniq.{mutation_type}.pos.stat.CIMS.txt"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CIMS_s30_txts = expand(os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.txt"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CIMS_s30_beds = expand(os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.bed"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CIMS_s30_21nt_beds = expand(os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.21nt.bed"), replicate_label=replicate_labels, mutation_type=mutation_types),
        finemapped_CIMS_fa = expand(os.path.join(OUTPUT, "sequence/CIMS.{replicate_label}.{mutation_type}.foreground.fa"), replicate_label=replicate_labels, mutation_type=mutation_types),
        background_CIMS_fa = expand(os.path.join(OUTPUT, "sequence/CIMS.{replicate_label}.{mutation_type}.background.fa"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CITS = expand(os.path.join(OUTPUT, "ctk/CITS/{replicate_label}.uniq.clean.CITS.s30.{mutation_type}.bed"), replicate_label=replicate_labels, mutation_type=mutation_types),
        mutation_CITS_s30_singleton_bed = expand(os.path.join(OUTPUT, "ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.{mutation_type}.bed"), replicate_label=replicate_labels, mutation_type=mutation_types),    
        mutation_CITS_s30_singleton_n21_bed = expand(os.path.join(OUTPUT, "ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.21nt.{mutation_type}.bed"), replicate_label=replicate_labels, mutation_type=mutation_types),
        finemapped_CITS_fa = expand(os.path.join(OUTPUT, "sequence/CITS.{replicate_label}.{mutation_type}.foreground.fa"), replicate_label=replicate_labels, mutation_type=mutation_types),
        background_CITS_fa = expand(os.path.join(OUTPUT, "sequence/CITS.{replicate_label}.{mutation_type}.background.fa"), replicate_label=replicate_labels, mutation_type=mutation_types),
        clip_finemapped_fa = expand(os.path.join(OUTPUT, "sequence/{experiment_label}.foreground.fa"), experiment_label=experiment_labels),
        clip_background_fa = expand(os.path.join(OUTPUT, "sequence/{experiment_label}.background.fa"), experiment_label=experiment_labels),
        ctk_kmer_enrichment = expand(os.path.join(OUTPUT, "ctk/ctk_mcross/{data_type}.{replicate_label}.{mutation_type}.kmer.txt"), data_type=data_types, replicate_label=replicate_labels, mutation_type=mutation_types),
        ctk_config = expand(os.path.join(OUTPUT, "ctk/ctk_mcross/{data_type}.{replicate_label}.{mutation_type}.config.txt"), data_type=data_types, replicate_label=replicate_labels, mutation_type=mutation_types),
        ctk_topn_kmer_matrix = expand(os.path.join(OUTPUT, "ctk/ctk_mcross/{data_type}.{replicate_label}.{mutation_type}.w7.zcore.mat.txt"), data_type=data_types, replicate_label=replicate_labels, mutation_type=mutation_types),
        ctk_top_peak = expand(os.path.join(OUTPUT, "ctk/ctk_mcross/top7mer/{data_type}.{replicate_label}.{mutation_type}/top.X.{data_type}.{replicate_label}.{mutation_type}.txt"), data_type=data_types, replicate_label=replicate_labels, mutation_type=mutation_types),
        ctk_mat = expand(os.path.join(OUTPUT, "ctk/ctk_mcross/mcross/{data_type}.{replicate_label}.{mutation_type}/{data_type}.{replicate_label}.{mutation_type}.00.mat"), data_type=data_types, replicate_label=replicate_labels, mutation_type=mutation_types),
        skipper_kmer_enrichment = expand(os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_label}.kmer.txt"), experiment_label=experiment_labels),
        skipper_clip_config = expand(os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_label}.config.txt"), experiment_label=experiment_labels),
        skipper_clip_topn_kmer_matrix = expand(os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_label}.w7.zcore.mat.txt"), experiment_label=experiment_labels),
        skipper_top_peak = expand(os.path.join(OUTPUT, "ctk/skipper_mcross/top7mer/{experiment_label}/top.X.{experiment_label}.txt"), experiment_label=experiment_labels),
        skipper_mat = expand(os.path.join(OUTPUT, "ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.mat"), experiment_label=experiment_labels),

        

rule select_informative_read:
    input:
        bam_combined=os.path.join(OUTPUT, "bams/genome/{replicate_labels}.genome.Aligned.sort.dedup.bam"),
    output:
        bam_informative=os.path.join(OUTPUT, "bams/genome_R1/{replicate_labels}.genome.Aligned.sort.dedup.R1.bam"),
    params:
        error_out_file = os.path.join(OUTPUT, "stderr/{replicate_labels}.select_informative_read.err"),
        out_file = os.path.join(OUTPUT, "stdout/{replicate_labels}.select_informative_read.out"),
        run_time = "00:30:00",
        memory = "20000",
        job_name = "select_informative_read {replicate_labels}",
    benchmark: os.path.join(OUTPUT, "benchmarks/select/unassigned_experiment.{replicate_labels}.select_informative_read.txt"),
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    threads: 8
    shell:
        "samtools view -bF 128 {input.bam_combined} > {output.bam_informative}"
        

rule uniquely_mapped_reads:
    input:
        bam = rules.select_informative_read.output.bam_informative
    output:
        bam_umap = os.path.join(OUTPUT, "bams/genome/{replicate_labels}.genome.Aligned.sort.dedup.umap.bam"),
        bai_umap = os.path.join(OUTPUT, "bams/genome/{replicate_labels}.genome.Aligned.sort.dedup.umap.bam.bai"),
    threads: 1
    params:
        error_out_file = os.path.join(OUTPUT, "stderr/uniquemap.{replicate_labels}.err"),
        out_file = os.path.join(OUTPUT, "stdout/uniquemap.{replicate_labels}.out"),
        run_time = "0:30:00",
        memory = "20000",
        job_name = "filter_map_quality {replicate_labels}",
    threads: 8
    shell:
        """
        module load bamtools;
        bamtools filter -in {input.bam} -out {output.bam_umap} -mapQuality ">3"
        samtools index {output.bam_umap}
        """
        
rule md_tag:
    input:
        bam_ip_umap = os.path.join(OUTPUT, "bams/genome/{replicate_labels}.genome.Aligned.sort.dedup.umap.bam"),
    output:
        mdtag = os.path.join(OUTPUT, "bams/genome/{replicate_labels}.genome.Aligned.sort.dedup.umap.md.sam"),
    params:
        fa = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        error_out_file = os.path.join(OUTPUT, "error_files/md_tag.{replicate_labels}.err"),
        out_file = os.path.join(OUTPUT, "stdout/md_tag.{replicate_labels}.out"),
        run_time = "1:00:00",
        memory = "5000",
        cores = 1,
        job_name = "fillmd {replicate_labels}",
    threads: 8
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        samtools fillmd {input.bam_ip_umap} {params.fa} > {output.mdtag}
        """
    
    
rule ctk_parse:
    input:
        mdtag = os.path.join(OUTPUT, "bams/genome/{replicate_labels}.genome.Aligned.sort.dedup.umap.md.sam"),
    output:
        mutation_file = os.path.join(OUTPUT, "ctk/{replicate_labels}.mutation.txt"),
        tagbed = os.path.join(OUTPUT, "ctk/{replicate_labels}.tag.bed"),
        peak = os.path.join(OUTPUT, "ctk/{replicate_labels}.uniq.peak.sig.bed"),
        peak_bd = os.path.join(OUTPUT, "ctk/{replicate_labels}.uniq.peak.sig.boundary.bed"),
        peak_PH = os.path.join(OUTPUT, "ctk/{replicate_labels}.uniq.peak.sig.halfPH.bed"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/ctk_parse.{replicate_labels}.err"),
        out_file = os.path.join(OUTPUT, "stdout/ctk_parse.{replicate_labels}.out"),
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "parse {replicate_labels}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        perl /tscc/nfs/home/s5xu/projects/ctk-1.1.5/parseAlignment.pl \
            -v --map-qual 1 \
            --min-len 18 \
            --mutation-file {output.mutation_file} {input.mdtag} {output.tagbed}

        /tscc/nfs/home/s5xu/projects/ctk-1.1.5/tag2peak.pl -big -ss \
            -v --valley-seeking -p 0.05 --valley-depth 0.9 \
            --multi-test --dbkey hg38 \
            {output.tagbed} \
            {output.peak} \
            --out-boundary {output.peak_bd} \
            --out-half-PH {output.peak_PH}
        """


rule ctk_get_mutation:
    wildcard_constraints:
        mutation_types="(del|ins|sub)",
    input:
        mutation = rules.ctk_parse.output.mutation_file,
    output:
        mutation_bed = os.path.join(OUTPUT, "ctk/{replicate_labels}.uniq.{mutation_types}.bed"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files", "ctk_get_mutation.{replicate_labels}.{mutation_types}.err"),
        out_file = os.path.join(OUTPUT, "stdout", "ctk_get_mutation.{replicate_labels}.{mutation_types}.out"),
        mutation_type = "{mutation_types}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "get_mutation {replicate_labels} {mutation_types}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        perl /tscc/nfs/home/s5xu/projects/ctk-1.1.5/getMutationType.pl -t {params.mutation_type} {input.mutation} {output.mutation_bed}
        """
        
rule ctk_CIMS:
    input:
        tagbed = rules.ctk_parse.output.tagbed,
        mutation = os.path.join(OUTPUT, "ctk/{replicate_labels}.uniq.{mutation_types}.bed"),
    output:
        mutation_CIMS = os.path.join(OUTPUT, "ctk/CIMS/{replicate_labels}.uniq.{mutation_types}.CIMS.txt"),
        mutation_CIMS_stat = os.path.join(OUTPUT, "ctk/CIMS/{replicate_labels}.uniq.{mutation_types}.pos.stat.CIMS.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/ctk_CIMS.{replicate_labels}.{mutation_types}.err"),
        out_file = os.path.join(OUTPUT, "stdout/ctk_CIMS.{replicate_labels}.{mutation_types}.out"),
        run_time = "3:00:00",
        memory = "20000",
        cores = 1,
        job_name = "CIMS {replicate_labels} {mutation_types}",
        tmp_dir = os.path.join(OUTPUT, "tmp", "CIMS.{replicate_labels}.{mutation_types}"),
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        perl /tscc/nfs/home/s5xu/projects/ctk-1.1.5/CIMS.pl -big -n 10 -p -c {params.tmp_dir} -outp {output.mutation_CIMS_stat} -v {input.tagbed} {input.mutation} {output.mutation_CIMS}
        """
        
rule ctk_CIMS_process:
    wildcard_constraints:
        mutation_types="(del|ins|sub)",
    input:
        mutation_CIMS = os.path.join(OUTPUT, "ctk/CIMS/{replicate_labels}.uniq.{mutation_types}.CIMS.txt"),
    output:
        mutation_CIMS_s30_txt = os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_labels}.uniq.{mutation_types}.CIMS.s30.txt"),
        mutation_CIMS_s30_bed = os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_labels}.uniq.{mutation_types}.CIMS.s30.bed"),
        mutation_CIMS_s30_21nt_bed = os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_labels}.uniq.{mutation_types}.CIMS.s30.21nt.bed"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/ctk_CIMS_process.{replicate_labels}.{mutation_types}.err"),
        out_file = os.path.join(OUTPUT, "stdout/ctk_CIMS_process.{replicate_labels}.{mutation_types}.out"),
        run_time = "10:00",
        memory = "5000",
        cores = 1,
        job_name = "CIMS process {replicate_labels} {mutation_types}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        awk '{{if($9<=0.001) {{print $0}}}}' {input.mutation_CIMS} | sort -k 9,9n -k 8,8nr -k 7,7n > {output.mutation_CIMS_s30_txt}
        cut -f 1-6  {output.mutation_CIMS_s30_txt} > {output.mutation_CIMS_s30_bed}
        awk '{{print $1"\\t"$2-10"\\t"$3+10"\\t"$4"\\t"$5"\\t"$6}}' {output.mutation_CIMS_s30_bed} > {output.mutation_CIMS_s30_21nt_bed}
        """
        
rule fetch_CIMS_sequence:
    input:
        mutation_CIMS_s30_21nt_bed = os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_labels}.uniq.{mutation_types}.CIMS.s30.21nt.bed"),
    output:
        finemapped_CIMS_fa = os.path.join(OUTPUT, "sequence/CIMS.{replicate_labels}.{mutation_types}.foreground.fa"),
        background_CIMS_fa = os.path.join(OUTPUT, "sequence/CIMS.{replicate_labels}.{mutation_types}.background.fa"),
    params:
        error_out_file = os.path.join(OUTPUT, "stderr/CIMS.{replicate_labels}.{mutation_types}.fetch_sequence.err"),
        out_file = os.path.join(OUTPUT, "stdout/CIMS.{replicate_labels}.{mutation_types}.fetch_sequence.out"),
        run_time = "30:00",
        memory = "5000",
        job_name = "fetch_CIMS_{replicate_labels}.{mutation_types}_sequence",
        fa = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        tmp_dir = os.path.join(OUTPUT, "tmp/CIMS_fetch_sequence"),
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        mkdir -p {params.tmp_dir}

        awk -v OFS="\\t" '{{center=$2+10; print $1, center-50, center+50, $4, $5, $6}}' \
            {input.mutation_CIMS_s30_21nt_bed} > {params.tmp_dir}/{params.job_name}.foreground.bed
        bedtools getfasta -fo {output.finemapped_CIMS_fa} -fi {params.fa} -bed {params.tmp_dir}/{params.job_name}.foreground.bed -s
        awk -v OFS="\\t" '{{center=$2+10; print $1, center-550, center-450, $4"_up", $5, $6; print $1, center+450, center+550, $4"_down", $5, $6}}' \
            {input.mutation_CIMS_s30_21nt_bed} > {params.tmp_dir}/{params.job_name}.background.bed
        bedtools getfasta -fo {output.background_CIMS_fa} -fi {params.fa} -bed {params.tmp_dir}/{params.job_name}.background.bed -s
        
        rm -rf {params.tmp_dir}
        """

        
        
rule ctk_CITS:
    input:
        tagbed = rules.ctk_parse.output.tagbed,
        mutation = os.path.join(OUTPUT, "ctk/{replicate_labels}.uniq.{mutation_types}.bed"),
    output:
        mutation_CITS = os.path.join(OUTPUT, "ctk/CITS/{replicate_labels}.uniq.clean.CITS.s30.{mutation_types}.bed"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/ctk_CITS.{replicate_labels}.{mutation_types}.err"),
        out_file = os.path.join(OUTPUT, "stdout/ctk_CITS.{replicate_labels}.{mutation_types}.out"),
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "CITS {replicate_labels} {mutation_types}",
        tmp_dir = os.path.join(OUTPUT, "tmp", "CITS.{replicate_labels}.{mutation_types}"),
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        perl /tscc/nfs/home/s5xu/projects/ctk-1.1.5/CITS.pl -big -p 0.001 --gap 25 -c {params.tmp_dir} -v {input.tagbed} {input.mutation} {output.mutation_CITS}
        """   
                       
        
rule ctk_CITS_process:
    wildcard_constraints:
        mutation_types="(del|ins|sub)",
    input:
        mutation_CITS = os.path.join(OUTPUT, "ctk/CITS/{replicate_labels}.uniq.clean.CITS.s30.{mutation_types}.bed"),
    output:
        mutation_CITS_s30_singleton_bed = os.path.join(OUTPUT, "ctk/CITS/sig/{replicate_labels}.uniq.clean.CITS.s30.singleton.{mutation_types}.bed"),    
        mutation_CITS_s30_singleton_n21_bed = os.path.join(OUTPUT, "ctk/CITS/sig/{replicate_labels}.uniq.clean.CITS.s30.singleton.21nt.{mutation_types}.bed"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/ctk_CITS_process.{replicate_labels}.{mutation_types}.err"),
        out_file = os.path.join(OUTPUT, "stdout/ctk_CITS_process.{replicate_labels}.{mutation_types}.out"),
        run_time = "10:00",
        memory = "5000",
        cores = 1,
        job_name = "CITS process {replicate_labels} {mutation_types}",
    threads: 8
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"     
    shell:
        """
        awk '{{if($3-$2==1) {{print $0}}}}' {input.mutation_CITS} > {output.mutation_CITS_s30_singleton_bed}
        awk '{{print $1"\\t"$2-10"\\t"$3+10"\\t"$4"\\t"$5"\\t"$6}}' {output.mutation_CITS_s30_singleton_bed} > {output.mutation_CITS_s30_singleton_n21_bed}
        """

rule fetch_CITS_sequence:
    input:
        mutation_CITS_s30_singleton_n21_bed = os.path.join(OUTPUT, "ctk/CIMS/sig/{replicate_labels}.uniq.{mutation_types}.CIMS.s30.21nt.bed"),
    output:
        finemapped_CITS_fa = os.path.join(OUTPUT, "sequence/CITS.{replicate_labels}.{mutation_types}.foreground.fa"),
        background_CITS_fa = os.path.join(OUTPUT, "sequence/CITS.{replicate_labels}.{mutation_types}.background.fa"),
    params:
        error_out_file = os.path.join(OUTPUT, "stderr/CITS.{replicate_labels}.{mutation_types}.fetch_sequence.err"),
        out_file = os.path.join(OUTPUT, "stdout/CITS.{replicate_labels}.{mutation_types}.fetch_sequence.out"),
        run_time = "30:00",
        memory = "20000",
        job_name = "fetch_CITS_{replicate_labels}.{mutation_types}_sequence",
        fa = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        tmp_dir = os.path.join(OUTPUT, "tmp/CITS_fetch_sequence"),
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        mkdir -p {params.tmp_dir}

        awk -v OFS="\\t" '{{center=$2+10; print $1, center-50, center+50, $4, $5, $6}}' \
            {input.mutation_CITS_s30_singleton_n21_bed} > {params.tmp_dir}/{params.job_name}.foreground.bed
        bedtools getfasta -fo {output.finemapped_CITS_fa} -fi {params.fa} -bed {params.tmp_dir}/{params.job_name}.foreground.bed -s
        awk -v OFS="\\t" '{{center=$2+10; print $1, center-550, center-450, $4"_up", $5, $6; print $1, center+450, center+550, $4"_down", $5, $6}}' \
            {input.mutation_CITS_s30_singleton_n21_bed} > {params.tmp_dir}/{params.job_name}.background.bed
        bedtools getfasta -fo {output.background_CITS_fa} -fi {params.fa} -bed {params.tmp_dir}/{params.job_name}.background.bed -s

        rm -rf {params.tmp_dir}
        """
         
rule clip_fetch_sequence:
    input:
        finemapped_windows = os.path.join(OUTPUT, "finemapping/mapped_sites/{experiment_labels}.finemapped_windows.bed.gz"),
        background = os.path.join(OUTPUT, "homer/region_matched_background/fixed/{experiment_labels}.sampled_fixed_windows.bed.gz"),
    output:
        clip_finemapped_fa = os.path.join(OUTPUT, "sequence/{experiment_labels}.foreground.fa"),
        clip_background_fa = os.path.join(OUTPUT, "sequence/{experiment_labels}.background.fa"),
    params:
        error_out_file = os.path.join(OUTPUT, "stderr/{experiment_labels}.fetch_sequence.err"),
        out_file = os.path.join(OUTPUT, "stdout/{experiment_labels}.fetch_sequence.out"),
        run_time = "1:00:00",
        memory = "20000",
        job_name = "run_homer",
        fa = "/tscc/projects/ps-yeolab4/genomes/GRCh38/chromosomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        '''
        bedtools getfasta -fo {output.clip_finemapped_fa} -fi {params.fa} -bed {input.finemapped_windows} -s
        bedtools getfasta -fo {output.clip_background_fa} -fi {params.fa} -bed {input.background} -s
        '''
        
        
rule ctk_mcross_kmer_enrichment:
    input:
        ctk_foreground = os.path.join(OUTPUT, "sequence/{data_types}.{replicate_labels}.{mutation_types}.foreground.fa"),
        ctk_background = os.path.join(OUTPUT, "sequence/{data_types}.{replicate_labels}.{mutation_types}.background.fa"),
    output:
        ctk_kmer_enrichment = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.kmer.txt"),
        ctk_config = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.config.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{data_types}.{replicate_labels}.{mutation_types}.ctk_mcross_kmer_enrich.err"),
        out_file = os.path.join(OUTPUT, "stdout/{data_types}.{replicate_labels}.{mutation_types}.ctk_mcross_kmer_enrich.out"),
        experiment_label = "{data_types}.{replicate_labels}.{mutation_types}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "ctk_mcross {data_types}.{replicate_labels}.{mutation_types}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/word_enrich.pl -w 7 \
            -test binom -v \
            {input.ctk_finemapped_fa} \
            {input.ctk_background_fa} \
            {output.ctk_kmer_enrichment}
        
        # generate config
        echo '{output.ctk_kmer_enrichment}\t\{params.experiment_label}' > {output.ctk_config}
        """
        
        
rule ctk_mcross_enrich_matrix:
    input:
        ctk_kmer_enrichment = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.kmer.txt"),
        ctk_config = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.config.txt"),
    output:
        ctk_topn_kmer_matrix = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.w7.zcore.mat.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{data_types}.{replicate_labels}.{mutation_types}.enrich_matrix.err"),
        out_file = os.path.join(OUTPUT, "stdout/{data_types}.{replicate_labels}.{mutation_types}.enrich_matrix.out"),
        experiment_label = "{data_types}.{replicate_labels}.{mutation_types}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "enrich_matrix {data_types}.{replicate_labels}.{mutation_types}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/gen_word_enrich_matrix.pl  \
            {input.ctk_config}  {output.ctk_topn_kmer_matrix}
        """
        
        
rule ctk_mcross_topword:
    input:
        ctk_topn_kmer_matrix = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.w7.zcore.mat.txt"),
    output:
        ctk_top_peak = os.path.join(OUTPUT, "ctk/ctk_mcross/top7mer/{data_types}.{replicate_labels}.{mutation_types}/top.X.{data_types}.{replicate_labels}.{mutation_types}.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{data_types}.{replicate_labels}.{mutation_types}.mcross_topword.err"),
        out_file = os.path.join(OUTPUT, "stdout/{data_types}.{replicate_labels}.{mutation_types}.mcross_topword.out"),
        output_dir = os.path.join(OUTPUT, "ctk/ctk_mcross/top7mer"),
        experiment_label = "{data_types}.{replicate_labels}.{mutation_types}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "mcross_topword {data_types}.{replicate_labels}.{mutation_types}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/topword.R {input.ctk_topn_kmer_matrix} {params.output_dir}/{params.experiment_label}
        """
        
        
rule ctk_mcross:
    input:
        clip_finemapped_fa = os.path.join(OUTPUT, "sequence/{data_types}.{replicate_labels}.{mutation_types}.foreground.fa"),
        clip_kmer_enrichment = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.kmer.txt"),
        clip_topn_kmer_matrix = os.path.join(OUTPUT, "ctk/ctk_mcross/{data_types}.{replicate_labels}.{mutation_types}.w7.zcore.mat.txt"),
        clip_top_peak = os.path.join(OUTPUT, "ctk/ctk_mcross/top7mer/{data_types}.{replicate_labels}.{mutation_types}/top.X.{data_types}.{replicate_labels}.{mutation_types}.txt"),
    output:
        mat = os.path.join(OUTPUT, "ctk/ctk_mcross/mcross/{data_types}.{replicate_labels}.{mutation_types}/{data_types}.{replicate_labels}.{mutation_types}.00.mat")
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{data_types}.{replicate_labels}.{mutation_types}.ctk_mcross.err"),
        out_file = os.path.join(OUTPUT, "stdout/{data_types}.{replicate_labels}.{mutation_types}.ctk_mcross.out"),
        output_dir = os.path.join(OUTPUT, "ctk/ctk_mcross/mcross/{data_types}.{replicate_labels}.{mutation_types}"),
        experiment_label = "{data_types}.{replicate_labels}.{mutation_types}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "ctk_mcross {data_types}.{replicate_labels}.{mutation_types}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/mCross.pl -l 10 -p 2 -N 10 -m 1 --cluster-seeds --seed {input.clip_top_peak} --prefix {params.experiment_label} --score-method sqrt {input.clip_finemapped_fa} {params.output_dir}/{params.experiment_label}
        """  
    
        
        
rule skipper_mcross_kmer_enrichment:
    input:
        clip_finemapped_fa = os.path.join(OUTPUT, "sequence/{experiment_labels}.foreground.fa"),
        clip_background_fa = os.path.join(OUTPUT, "sequence/{experiment_labels}.background.fa"),
    output:
        clip_kmer_enrichment = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.kmer.txt"),
        clip_config = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.config.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{experiment_labels}.skipper_mcross_kmer_enrich.err"),
        out_file = os.path.join(OUTPUT, "stdout/{experiment_labels}.skipper_mcross_kmer_enrich.out"),
        experiment_label = "{experiment_labels}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "skipper_mcross {experiment_labels}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/word_enrich.pl -w 7 \
            -test binom -v \
            {input.clip_finemapped_fa} \
            {input.clip_background_fa} \
            {output.clip_kmer_enrichment}
        
        # generate config
        echo '{output.clip_kmer_enrichment}\t\{params.experiment_label}' > {output.clip_config}
        """
        
        
rule skipper_mcross_enrich_matrix:
    input:
        clip_kmer_enrichment = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.kmer.txt"),
        clip_config = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.config.txt"),
    output:
        clip_topn_kmer_matrix = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.w7.zcore.mat.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{experiment_labels}.enrich_matrix.err"),
        out_file = os.path.join(OUTPUT, "stdout/{experiment_labels}.enrich_matrix.out"),
        experiment_label = "{experiment_labels}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "enrich_matrix {experiment_labels}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/gen_word_enrich_matrix.pl  \
            {input.clip_config}  {output.clip_topn_kmer_matrix}
        """
        

rule skipper_mcross_topword:
    input:
        clip_topn_kmer_matrix = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.w7.zcore.mat.txt"),
    output:
        clip_top_peak = os.path.join(OUTPUT, "ctk/skipper_mcross/top7mer/{experiment_labels}/top.X.{experiment_labels}.txt"),
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{experiment_labels}.mcross_topword.err"),
        out_file = os.path.join(OUTPUT, "stdout/{experiment_labels}.mcross_topword.out"),
        output_dir = os.path.join(OUTPUT, "ctk/skipper_mcross/top7mer"),
        experiment_label = "{experiment_labels}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "mcross_topword {experiment_labels}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/topword.R {input.clip_topn_kmer_matrix} {params.output_dir}/{params.experiment_label}
        """
        
        
rule skipper_mcross:
    input:
        clip_finemapped_fa = os.path.join(OUTPUT, "sequence/{experiment_labels}.foreground.fa"),
        clip_kmer_enrichment = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.kmer.txt"),
        clip_topn_kmer_matrix = os.path.join(OUTPUT, "ctk/skipper_mcross/{experiment_labels}.w7.zcore.mat.txt"),
        clip_top_peak = os.path.join(OUTPUT, "ctk/skipper_mcross/top7mer/{experiment_labels}/top.X.{experiment_labels}.txt"),
    output:
        mat = os.path.join(OUTPUT, "ctk/skipper_mcross/mcross/{experiment_labels}/{experiment_labels}.00.mat")
    params:
        error_out_file = os.path.join(OUTPUT, "error_files/{experiment_labels}.mcross.err"),
        out_file = os.path.join(OUTPUT, "stdout/{experiment_labels}.mcross.out"),
        output_dir = os.path.join(OUTPUT, "ctk/skipper_mcross/mcross/{experiment_labels}"),
        experiment_label = "{experiment_labels}",
        run_time = "2:00:00",
        memory = "20000",
        cores = 1,
        job_name = "mcross {experiment_labels}",
    threads: 8
    # container:
    #     "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        /tscc/nfs/home/s5xu/projects/mCross/mCross.pl -l 10 -p 2 -N 10 -m 1 --cluster-seeds --seed {input.clip_top_peak} --prefix {params.experiment_label} --score-method sqrt {input.clip_finemapped_fa} {params.output_dir}/{params.experiment_label}
        """
        
        
        
