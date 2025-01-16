VARIANT='/tscc/nfs/home/hsher/ps-yeolab5/ENCODE_paper_tables/CHD_variants.tsv'
#'CHROM','POS','ID','REF','ALT' first 5 columns needs to be this, no index, no header


TABLE='/tscc/projects/ps-yeolab4/software/skipper/bb63a25/bin/skipper/annotations/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.Heart_Left_Ventricle.tiled_partition.bed.gz'
VEP_CACHEDIR='/tscc/nfs/home/hsher/scratch/vep_cache/'
import pandas as pd
locals().update(config)

tabledir = Path('/tscc/nfs/home/hsher/ps-yeolab5/ENCODE_paper_tables/')
metrics=pd.read_csv(tabledir/'model_performance.csv')
selected_models = metrics.loc[metrics['selected'], 'Experiment'].sort_values().tolist()

workdir: "/tscc/nfs/home/hsher/scratch/CHD_analysis"
from pathlib import Path
model_dir = Path('/tscc/nfs/home/hsher/scratch/')
models = model_dir.glob('ENCO*/output/ml/rbpnet_model/*')
model_dict = {}
for f in models:
    model_dict[f.name] = f


print(model_dict)
"""
snakemake -s utils/score_any_variants.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2 \
    -n
"""

rule all:
    input:
        expand("output/variants/{model_name}.score.csv",
        model_name = selected_models
        ),
        #"output/variants/vep.tsv"


rule slop_finemap:
    input:
        finemapped_windows = TABLE
    output:
        "output/table.slop.bed.gz"
    threads: 2
    params:
        error_file = "stderr/slop_finemap",
        out_file = "stdout/slop_finemap",
        run_time = "06:20:00",
        cores = 1,
        memory = 40000,
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
        """
        bedtools slop -i {input.finemapped_windows} -g {CHROM_SIZES} -b 100 | gzip -c > {output}
        """

rule fetch_peak_sequence:
    input:
        finemapped_windows = rules.slop_finemap.output,
    output:
        finemapped_fa = "output/sequence/table.slop.fa",
    params:
        error_file = "stderr/fetch_sequence.err",
        out_file = "stdout/fetch_sequence.out",
        run_time = "40:00",
        memory = "2000",
        job_name = "run_homer",
        fa = config['GENOME']
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    shell:
        '''
        bedtools getfasta -fo {output.finemapped_fa} -fi {params.fa} -bed {input.finemapped_windows} -s
        '''

rule fetch_variant_sequence:
    input:
        subset_vcf=VARIANT,
        seq_fa = rules.fetch_peak_sequence.output.finemapped_fa,
        finemapped_windows = rules.slop_finemap.output
    output:
        ref_fa = temp("output/variants/annotation.ref.fa"),
        alt_fa = temp("output/variants/annotation.alt.fa"),
        all_fa = "output/variants/annotation.all.fa",
        csv = "output/variants/annotation.csv"
    threads: 2
    params:
        error_file = "stderr/fetch_variant",
        out_file = "stdout/fetch_variant",
        run_time = "01:20:00",
        cores = 1,
        out_prefix = lambda wildcards, output: output.csv.replace('.csv', ''),
        memory = 80000,
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/metadensity.yaml"
    shell:
        """
        if [ -s {input.subset_vcf} ]; then
            python {TOOL_DIR}/generate_variant_sequence.py \
                {input.subset_vcf} \
                {input.seq_fa} \
                {input.finemapped_windows} \
                {params.out_prefix}
        else
            touch {output.csv}
            touch {output.ref_fa}
            touch {output.alt_fa}
        fi
        """

rule score_fa:
    input:
        fa = "output/variants/annotation.all.fa"
    output:
        score=temp("output/variants/{model_name}.score.csv"),
    threads: 1
    params:
        error_file = "stderr/score_fa.{model_name}",
        out_file = "stdout/score_fa.{model_name}",
        run_time = "00:20:00",
        cores = 1,
        memory = 80000,
        model_path = lambda wildcards: model_dict[wildcards.model_name]
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/eugene.yaml"
    shell:
        """
        if [ -s {input.fa} ]; then
            python {RBPNET_PATH}/score_fa.py \
                {params.model_path} \
                {input.fa} \
                {output.score} \
                /tscc/lustre/ddn/scratch/${{USER}}
        else
            touch {output.score}
        fi
        """

rule to_vcf:
    input:
        TABLE
    output:
        "output/variants/variants.vcf"
    container:
        "docker://brianyee/bcftools:1.17"
    params:
        error_file = "stderr/tovcf",
        out_file = "stdout/tovcf",
        run_time = "00:20:00",
        cores = 1,
        memory = 80000,
    shell:
        """
        bcftools convert -c CHROM,POS,REF,ALT,-,-,filter,- -f {GENOME} --tsv2vcf {input} -o {output}
        """
rule vep:
    input:
        "output/variants/variants.vcf"
    output:
        "output/variants/vep.tsv"
    threads: 2
    params:
        error_file = "stderr/vep",
        out_file = "stdout/vep",
        run_time = "1:20:00",
        cores = 1,
        memory = 40000,
        cache= VEP_CACHEDIR
    container:
        "docker://ensemblorg/ensembl-vep:latest"
    shell:
        """
        vep \
        -i {input} \
        --force_overwrite \
        -o {output} -offline --cache {params.cache}
        """

