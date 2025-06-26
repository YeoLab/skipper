"""
mamba install -c bioconda perl-math-cdf perl-bio-featureio perl-env
export PERL5LIB=/tscc/nfs/home/s5xu/projects/czplib-v.1.1.1
export PATH=$PATH:/tscc/nfs/home/s5xu/projects/patternmatch
# /tscc/nfs/home/s5xu/projects/ctk-1.1.5/CITS.pltscc/nfs/home/s5xu/projects/ctk-1.1.5/CITS.pl
snakemake -kps /tscc/nfs/home/s5xu/projects/skipper/rules/ctk_mcross.py -j 30 -w 30 --use-singularity --singularity-prefix /tscc/lustre/ddn/scratch/s5xu/singularity --singularity-args "--bind /tscc" --rerun-incomplete
"""
locals().update(config)

# OUTPUT = "/tscc/nfs/home/s5xu/scratch/skipper_output"
# replicate_label = ["RBFOX2_K562_ENCSR756CKJ_IP_1"]
# experiment_label = ["RBFOX2_K562_ENCSR756CKJ"]
# INFORMATIVE_READ = 1
# UNINFORMATIVE_READ = 2

# mutation_type = ['del', 'ins', 'sub']
# data_types = ['CIMS', 'CITS']

# rule all_mcross:
#     input:
#         bam_informative = expand("output/bams/genome_R1/{replicate_label}.genome.Aligned.sort.dedup.R1.bam"), replicate_label=replicate_label),
#         bam_umap = expand("output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam"), replicate_label=replicate_label),
#         bai_umap = expand("output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam.bai"), replicate_label=replicate_label),
#         mdtag = expand("output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.md.sam"), replicate_label=replicate_label),
#         mutation_file = expand("output/ctk/{replicate_label}.mutation.txt"), replicate_label=replicate_label),
#         tagbed = expand("output/ctk/{replicate_label}.tag.bed"), replicate_label=replicate_label),
#         peak = expand("output/ctk/{replicate_label}.uniq.peak.sig.bed"), replicate_label=replicate_label),
#         peak_bd = expand("output/ctk/{replicate_label}.uniq.peak.sig.boundary.bed"), replicate_label=replicate_label),
#         peak_PH = expand("output/ctk/{replicate_label}.uniq.peak.sig.halfPH.bed"), replicate_label=replicate_label),
#         mutation_beds = expand("output/ctk", "{replicate_label}.uniq.{mutation_type}.bed"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CIMSs = expand("output/ctk/CIMS/{replicate_label}.uniq.{mutation_type}.CIMS.txt"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CIMS_stats = expand("output/ctk/CIMS/{replicate_label}.uniq.{mutation_type}.pos.stat.CIMS.txt"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CIMS_s30_txts = expand("output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.txt"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CIMS_s30_beds = expand("output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.bed"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CIMS_s30_21nt_beds = expand("output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.21nt.bed"), replicate_label=replicate_label, mutation_type=mutation_type),
#         finemapped_CIMS_fa = expand("output/sequence/CIMS.{replicate_label}.{mutation_type}.foreground.fa"), replicate_label=replicate_label, mutation_type=mutation_type),
#         background_CIMS_fa = expand("output/sequence/CIMS.{replicate_label}.{mutation_type}.background.fa"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CITS = expand("output/ctk/CITS/{replicate_label}.uniq.clean.CITS.s30.{mutation_type}.bed"), replicate_label=replicate_label, mutation_type=mutation_type),
#         mutation_CITS_s30_singleton_bed = expand("output/ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.{mutation_type}.bed"), replicate_label=replicate_label, mutation_type=mutation_type),    
#         mutation_CITS_s30_singleton_n21_bed = expand("output/ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.21nt.{mutation_type}.bed"), replicate_label=replicate_label, mutation_type=mutation_type),
#         finemapped_CITS_fa = expand("output/sequence/CITS.{replicate_label}.{mutation_type}.foreground.fa"), replicate_label=replicate_label, mutation_type=mutation_type),
#         background_CITS_fa = expand("output/sequence/CITS.{replicate_label}.{mutation_type}.background.fa"), replicate_label=replicate_label, mutation_type=mutation_type),
#         clip_finemapped_fa = expand("output/ml/sequence/{experiment_label}.foreground.fa"), experiment_label=experiment_label),
#         clip_background_fa = expand("output/sequence/{experiment_label}.background.fa"), experiment_label=experiment_label),
#         ctk_kmer_enrichment = expand("output/ctk/ctk_mcross/{data_type}.{replicate_label}.{mutation_type}.kmer.txt"), data_type=data_types, replicate_label=replicate_label, mutation_type=mutation_type),
#         ctk_config = expand("output/ctk/ctk_mcross/{data_type}.{replicate_label}.{mutation_type}.config.txt"), data_type=data_types, replicate_label=replicate_label, mutation_type=mutation_type),
#         ctk_topn_kmer_matrix = expand("output/ctk/ctk_mcross/{data_type}.{replicate_label}.{mutation_type}.w7.zcore.mat.txt"), data_type=data_types, replicate_label=replicate_label, mutation_type=mutation_type),
#         ctk_top_peak = expand("output/ctk/ctk_mcross/top7mer/{data_type}.{replicate_label}.{mutation_type}/top.X.{data_type}.{replicate_label}.{mutation_type}.txt"), data_type=data_types, replicate_label=replicate_label, mutation_type=mutation_type),
#         ctk_mat = expand("output/ctk/ctk_mcross/mcross/{data_type}.{replicate_label}.{mutation_type}/{data_type}.{replicate_label}.{mutation_type}.00.mat"), data_type=data_types, replicate_label=replicate_label, mutation_type=mutation_type),
#         ctk_mat_plot = expand("output/ctk/ctk_mcross/mcross/{data_type}.{replicate_label}.{mutation_type}/{data_type}.{replicate_label}.{mutation_type}.00.pdf"), data_type=data_types, replicate_label=replicate_label, mutation_type=mutation_type),
#         skipper_kmer_enrichment = expand("output/ctk/skipper_mcross/{experiment_label}.kmer.txt"), experiment_label=experiment_label),
#         skipper_clip_config = expand("output/ctk/skipper_mcross/{experiment_label}.config.txt"), experiment_label=experiment_label),
#         skipper_clip_topn_kmer_matrix = expand("output/ctk/skipper_mcross/{experiment_label}.w7.zcore.mat.txt"), experiment_label=experiment_label),
#         skipper_top_peak = expand("output/ctk/skipper_mcross/top7mer/{experiment_label}/top.X.{experiment_label}.txt"), experiment_label=experiment_label),
#         skipper_mat = expand("output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.mat"), experiment_label=experiment_label),
#         skipper_plot = expand("output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.pdf"), experiment_label=experiment_label),
        

        
        
# SE: bam="output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam",
# PE: "output/bams/dedup/genome_R"+str(INFORMATIVE_READ)+"/{replicate_label}.genome.Aligned.sort.dedup.R"+str(INFORMATIVE_READ)+".bam"
rule uniquely_mapped_reads:
    input:
        bam = lambda wildcards: "output/bams/raw/genome/{replicate_label}.genome.Aligned.sort.bam" if protocol == "ENCODE4" else 
            "output/bams/dedup/genome_R"+str(INFORMATIVE_READ)+"/{replicate_label}.genome.Aligned.sort.dedup.R"+str(INFORMATIVE_READ)+".bam",
    output:
        bam_umap = "output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam",
        bai_umap = "output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam.bai",
    resources:
        runtime = 30,
        mem_mb = "20000",
    threads: 8
    shell:
        """
        module load bamtools;
        module load samtools;
        bamtools filter -in {input.bam} -out {output.bam_umap} -mapQuality ">3"
        samtools index {output.bam_umap}
        """
        
        
rule md_tag:
    input:
        bam_ip_umap = "output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.bam",
    output:
        mdtag = "output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.md.sam",
    resources:
        runtime = "1h",
        mem_mb = "5000",
    threads: 8
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        samtools fillmd {input.bam_ip_umap} {GENOME} |
                    awk 'BEGIN {{
            valid = "ATCG"
            }}
            /^@/ {{ print; next }}  # Print headers
            {{
            md = ""; found = 0
            for (i = 12; i <= NF; i++) {{
                if ($i ~ /^MD:Z:/) {{
                md = substr($i, 6)
                if (md ~ /[^0-9ATCG^]/) found = 1
                break
                }}
            }}
            if (!found) print
            }}' > {output.mdtag}
        """
    
    
rule ctk_parse:
    input:
        mdtag = "output/bams/genome/{replicate_label}.genome.Aligned.sort.dedup.umap.md.sam"
    output:
        mutation_file = "output/ctk/{replicate_label}.mutation.txt",
        tagbed = "output/ctk/{replicate_label}.tag.bed",
        peak = "output/ctk/{replicate_label}.uniq.peak.sig.bed",
        peak_bd = "output/ctk/{replicate_label}.uniq.peak.sig.boundary.bed",
        peak_PH = "output/ctk/{replicate_label}.uniq.peak.sig.halfPH.bed",
    resources:
        runtime = "2h",
        mem_mb = "20000",
        tmpdir = "output/ctktmp",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        parseAlignment.pl \
            -v --map-qual 1 \
            --min-len 18 \
            --mutation-file {output.mutation_file} {input.mdtag} {output.tagbed}
        tag2peak.pl -big -ss \
             -v --valley-seeking -p 0.05 --valley-depth 0.9 \
             --multi-test --dbkey hg38 \
             {output.tagbed} \
             {output.peak} \
             --out-boundary {output.peak_bd} \
             --out-half-PH {output.peak_PH}
        """



rule ctk_get_mutation:
    input:
        mutation = rules.ctk_parse.output.mutation_file,
    output:
        mutation_bed = "output/ctk/{replicate_label}.uniq.{mutation_type}.bed",
    wildcard_constraints:
        mutation_type="(del|ins|sub)",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        getMutationType.pl -t {wildcards.mutation_type} {input.mutation} {output.mutation_bed}
        """
        
        
rule ctk_CIMS:
    input:
        tagbed = rules.ctk_parse.output.tagbed,
        mutation = "output/ctk/{replicate_label}.uniq.{mutation_type}.bed",
    output:
        mutation_CIMS = "output/ctk/CIMS/{replicate_label}.uniq.{mutation_type}.CIMS.txt",
        mutation_CIMS_stat = "output/ctk/CIMS/{replicate_label}.uniq.{mutation_type}.pos.stat.CIMS.txt",
    resources:
        runtime = "3h",
        mem_mb = "20000",
    conda: "envs/ctk.yaml"
    shell:
        """
        CIMS.pl -big -n 10 -p -outp {output.mutation_CIMS_stat} -v {input.tagbed} {input.mutation} {output.mutation_CIMS}
        """
        
        
rule ctk_CIMS_process:
    wildcard_constraints:
        mutation_type="(del|ins|sub)",
    input:
        mutation_CIMS = "output/ctk/CIMS/{replicate_label}.uniq.{mutation_type}.CIMS.txt",
    output:
        mutation_CIMS_s30_txt = "output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.txt",
        mutation_CIMS_s30_bed = "output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.bed",
        mutation_CIMS_s30_21nt_bed = "output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.21nt.bed",
    resources:
        runtime = 10,
        mem_mb = "5000",
    threads: 8
    shell:
        """
        awk '{{if($9<=0.001) {{print $0}}}}' {input.mutation_CIMS} | sort -k 9,9n -k 8,8nr -k 7,7n > {output.mutation_CIMS_s30_txt}
        cut -f 1-6  {output.mutation_CIMS_s30_txt} > {output.mutation_CIMS_s30_bed}
        awk '{{print $1"\\t"$2-10"\\t"$3+10"\\t"$4"\\t"$5"\\t"$6}}' {output.mutation_CIMS_s30_bed} > {output.mutation_CIMS_s30_21nt_bed}
        """
        
        
rule fetch_CIMS_sequence:
    input:
        mutation_CIMS_s30_21nt_bed = lambda wildcards: expand("output/ctk/CIMS/sig/{replicate_label}.uniq.{mutation_type}.CIMS.s30.21nt.bed",
        replicate_label = config['experiment_to_replicate_labels'][wildcards.experiment_label],
        mutation_type = ['del', 'ins', 'sub']),
    output:
        allbed = temp("output/ctk/CIMS/sig/{experiment_label}.mutation_CIMS_s30_21nt_bed.bed"),
        finemapped_CIMS_bed = temp("output/sequence/CIMS.{experiment_label}.foreground.bed"),
        background_CIMS_bed = temp("output/sequence/CIMS.{experiment_label}.background.bed"),
        finemapped_CIMS_fa = "output/sequence/CIMS.{experiment_label}.foreground.fa",
        background_CIMS_fa = "output/sequence/CIMS.{experiment_label}.background.fa",
    params:
        fa=GENOME,
    resources:
        runtime = 30,
        mem_mb = "5000",
    threads: 8
    shell:
        """
        cat {input} > {output.allbed}
        awk -v OFS="\\t" '{{center=$2+10; print $1, center-50, center+50, $4, $5, $6}}' \
            {output.allbed} > {output.finemapped_CIMS_bed}
        bedtools getfasta -fo {output.finemapped_CIMS_fa} -fi {params.fa} -bed {output.finemapped_CIMS_bed} -s
        awk -v OFS="\\t" '{{center=$2+10; print $1, center-550, center-450, $4"_up", $5, $6; print $1, center+450, center+550, $4"_down", $5, $6}}' \
            {output.allbed} > {output.background_CIMS_bed}
        bedtools getfasta -fo {output.background_CIMS_fa} -fi {params.fa} -bed {output.background_CIMS_bed} -s
        """

        
rule ctk_CITS:
    input:
        tagbed = rules.ctk_parse.output.tagbed,
        mutation = "output/ctk/{replicate_label}.uniq.{mutation_type}.bed",
        peak = "output/ctk/{replicate_label}.uniq.peak.sig.bed",
        peak_bd = "output/ctk/{replicate_label}.uniq.peak.sig.boundary.bed",
        peak_PH = "output/ctk/{replicate_label}.uniq.peak.sig.halfPH.bed",
    output:
        mutation_CITS = "output/ctk/CITS/{replicate_label}.uniq.clean.CITS.s30.{mutation_type}.bed",
    resources:
        runtime = "6h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        CITS.pl -big -p 0.001 --gap 25 -v {input.tagbed} {input.mutation} {output.mutation_CITS}
        """   
                       
        
rule ctk_CITS_process:
    input:
        mutation_CITS = "output/ctk/CITS/{replicate_label}.uniq.clean.CITS.s30.{mutation_type}.bed",
    output:
        mutation_CITS_s30_singleton_bed = "output/ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.{mutation_type}.bed",    
        mutation_CITS_s30_singleton_n21_bed = "output/ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.21nt.{mutation_type}.bed",
    wildcard_constraints:
        mutation_type="(del|ins|sub)",
    resources:
        runtime = 10,
        mem_mb = "5000",
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
        mutation_CITS_s30_singleton_n21_bed = lambda wildcards: expand("output/ctk/CITS/sig/{replicate_label}.uniq.clean.CITS.s30.singleton.21nt.{mutation_type}.bed",
        replicate_label=config['experiment_to_replicate_labels'][wildcards.experiment_label],
        mutation_type = ['del', 'ins', 'sub'])
    output:
        allbed = temp("output/ctk/CITS/sig/{experiment_label}.mutation_CITS_s30_singleton_n21_bed.bed"),
        finemapped_CITS_bed = temp("output/sequence/CITS.{experiment_label}.foreground.bed"),
        background_CITS_bed = temp("output/sequence/CITS.{experiment_label}.background.bed"),
        finemapped_CITS_fa = "output/sequence/CITS.{experiment_label}.foreground.fa",
        background_CITS_fa = "output/sequence/CITS.{experiment_label}.background.fa",
    params:
        fa=GENOME,
    resources:
        runtime = 30,
        mem_mb = "20000",
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        cat {input} > {output.allbed};
        awk '!seen[$1 FS $2 FS $3]++' {output.allbed} \
            | egrep "^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]" \
            | awk -v OFS="\t" '{{center=$2+10; print $1, center-50, center+50, $4, $5, $6}}' \
            | awk '$2 >= 0' > {output.finemapped_CITS_bed};
            
        bedtools getfasta -fo {output.finemapped_CITS_fa} -fi {params.fa} -bed {output.finemapped_CITS_bed} -s;
        
        awk '!seen[$1 FS $2 FS $3]++' {output.allbed} \
            | egrep "^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]" \
            | awk -v OFS="\\t" '{{center=$2+10; print $1, center-550, center-450, $4"_up", $5, $6; print $1, center+450, center+550, $4"_down", $5, $6}}' \
            | awk '$2 >= 0' > {output.background_CITS_bed}
        
        bedtools getfasta -fo {output.background_CITS_fa} -fi {params.fa} -bed {output.background_CITS_bed} -s
        """
        
rule ctk_mcross_kmer_enrichment:
    input:
        ctk_foreground = "output/sequence/{data_types}.{experiment_label}.foreground.fa",
        ctk_background = "output/sequence/{data_types}.{experiment_label}.background.fa",
    output:
        ctk_kmer_enrichment = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.kmer.txt",
        ctk_config = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.config.txt",
    params:
        experiment_label = "{data_types}.{experiment_label}",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        word_enrich.pl -w 7 \
            -test binom -v \
            {input.ctk_foreground} \
            {input.ctk_background} \
            {output.ctk_kmer_enrichment}
        
        # generate config
        echo '{output.ctk_kmer_enrichment}\t\{params.experiment_label}' > {output.ctk_config}
        """
        
        
rule ctk_mcross_enrich_matrix:
    input:
        ctk_kmer_enrichment = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.kmer.txt",
        ctk_config = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.config.txt",
    output:
        ctk_topn_kmer_matrix = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.w7.zcore.mat.txt",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        gen_word_enrich_matrix.pl  \
            {input.ctk_config}  {output.ctk_topn_kmer_matrix}
        """
        
        
rule ctk_mcross_topword:
    input:
        ctk_topn_kmer_matrix = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.w7.zcore.mat.txt",
    output:
        ctk_top_peak = "output/ctk/ctk_mcross/top7mer/{data_types}.{experiment_label}/top.X.{data_types}.{experiment_label}.txt",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        topword.R {input.ctk_topn_kmer_matrix} \
            output/ctk/ctk_mcross/top7mer/{wildcards.data_types}.{wildcards.experiment_label}/
        """
        
        
rule ctk_mcross:
    input:
        clip_finemapped_fa = "output/sequence/{data_types}.{experiment_label}.foreground.fa",
        clip_kmer_enrichment = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.kmer.txt",
        clip_topn_kmer_matrix = "output/ctk/ctk_mcross/{data_types}.{experiment_label}.w7.zcore.mat.txt",
        clip_top_peak = "output/ctk/ctk_mcross/top7mer/{data_types}.{experiment_label}/top.X.{data_types}.{experiment_label}.txt",
    output:
        mat = "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.00.mat"
    params:
        output_dir = "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}",
        experiment_label = "{data_types}.{experiment_label}",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        export PATH=$PATH:/tscc/nfs/home/s5xu/projects/patternmatch
        mCross.pl -l 10 -p 2 -N 10 -m 1 \
            --cluster-seeds --seed {input.clip_top_peak} --prefix {params.experiment_label} \
            --score-method sqrt {input.clip_finemapped_fa} {params.output_dir}/{params.experiment_label}
        """  
    
        
rule skipper_mcross_kmer_enrichment:
    input:
        clip_finemapped_fa = "output/ml/sequence/{experiment_label}.foreground.fa",
        clip_background_fa = "output/ml/sequence/{experiment_label}.background.fa",
    output:
        clip_kmer_enrichment = "output/ctk/skipper_mcross/{experiment_label}.kmer.txt",
        clip_config = "output/ctk/skipper_mcross/{experiment_label}.config.txt",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        word_enrich.pl -w 7 \
            -test binom -v \
            {input.clip_finemapped_fa} \
            {input.clip_background_fa} \
            {output.clip_kmer_enrichment}
        
        # generate config
        echo '{output.clip_kmer_enrichment}\t\{wildcards.experiment_label}' > {output.clip_config}
        """
        
        
rule skipper_mcross_enrich_matrix:
    input:
        clip_kmer_enrichment = "output/ctk/skipper_mcross/{experiment_label}.kmer.txt",
        clip_config = "output/ctk/skipper_mcross/{experiment_label}.config.txt",
    output:
        clip_topn_kmer_matrix = "output/ctk/skipper_mcross/{experiment_label}.w7.zcore.mat.txt",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        gen_word_enrich_matrix.pl  \
            {input.clip_config}  {output.clip_topn_kmer_matrix}
        """
        

rule skipper_mcross_topword:
    input:
        clip_topn_kmer_matrix = "output/ctk/skipper_mcross/{experiment_label}.w7.zcore.mat.txt",
    output:
        clip_top_peak = "output/ctk/skipper_mcross/top7mer/{experiment_label}/top.X.{experiment_label}.txt",
    params:
        output_dir = "output/ctk/skipper_mcross/top7mer",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        topword.R {input.clip_topn_kmer_matrix} \
            {params.output_dir}/{wildcards.experiment_label}
        """
        
        
rule skipper_mcross:
    input:
        clip_finemapped_fa = "output/ml/sequence/{experiment_label}.foreground.fa",
        clip_kmer_enrichment = "output/ctk/skipper_mcross/{experiment_label}.kmer.txt",
        clip_topn_kmer_matrix = "output/ctk/skipper_mcross/{experiment_label}.w7.zcore.mat.txt",
        clip_top_peak = "output/ctk/skipper_mcross/top7mer/{experiment_label}/top.X.{experiment_label}.txt",
    output:
        mat = "output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.mat"
    params:
        output_dir = "output/ctk/skipper_mcross/mcross/{experiment_label}",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    conda: "envs/ctk.yaml"
    shell:
        """
        export PATH=$PATH:/tscc/nfs/home/s5xu/projects/patternmatch
        mCross.pl -l 10 -p 2 -N 10 -m 1 --cluster-seeds \
            --seed {input.clip_top_peak} --prefix {wildcards.experiment_label} \
            --score-method sqrt {input.clip_finemapped_fa} {params.output_dir}/{wildcards.experiment_label}
        """
        
        
        
rule plot_skipper_mcross:
    input:
        skipper_mat = "output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.mat",
    output:
        skipper_plot = "output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.pdf",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    container:
        "docker://howardxu520/skipper:R_4.1.3_2"
    shell:
        """
        export PERL5LIB=/tscc/nfs/home/s5xu/projects/czplib-v.1.1.1;
        export PATH=$PATH:/tscc/nfs/home/s5xu/projects/patternmatch;
        
        Rscript /tscc/nfs/home/s5xu/projects/mCross/mCross2logo.R \
            -i {input.skipper_mat} -o {output.skipper_plot} -s rna -v
        """
        
        
rule plot_ctk_mcross:
    input:
        ctk_mat = "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.00.mat"
    output:
        ctk_mat_plot = "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.00.pdf",
    resources:
        runtime = "2h",
        mem_mb = "20000",
    threads: 8
    container:
        "docker://howardxu520/skipper:R_4.1.3_2"
    shell:
        """
        export PERL5LIB=/tscc/nfs/home/s5xu/projects/czplib-v.1.1.1;
        export PATH=$PATH:/tscc/nfs/home/s5xu/projects/patternmatch;
        
        Rscript /tscc/nfs/home/s5xu/projects/mCross/mCross2logo.R \
            -i {input.ctk_mat} -o {output.ctk_mat_plot} -s rna -v
        """

rule mcross2homer_skipper:
    input:
        "output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.00.mat"
    output:
        "output/ctk/skipper_mcross/mcross/{experiment_label}/{experiment_label}.homer"
    resources:
        runtime = 10,
        mem_mb=20000,
    threads: 1
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/mcross2homer.py \
            output/ctk/skipper_mcross/mcross/{wildcards.experiment_label} \
            {output}
        """
        
rule mcross2homer_ctk:
    input:
        "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.00.mat"
    output:
        "output/ctk/ctk_mcross/mcross/{data_types}.{experiment_label}/{data_types}.{experiment_label}.homer"
    resources:
        runtime = 10,
        mem_mb=20000,
    threads: 1
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/mcross2homer.py \
            output/ctk/ctk_mcross/mcross/{wildcards.data_types}.{wildcards.experiment_label} \
            {output}
        """
        
