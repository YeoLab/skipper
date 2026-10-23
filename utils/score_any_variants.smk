VARIANT='/tscc/nfs/home/hsher/ps-yeolab5/ENCODE_paper_tables/CHD_variants.tsv'
#'CHROM','POS','ID','REF','ALT' first 5 columns needs to be this, no index, no header


TABLE='/tscc/projects/ps-yeolab4/software/skipper/bb63a25/bin/skipper/annotations/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.Heart_Left_Ventricle.tiled_partition.bed.gz'
import pandas as pd
from pathlib import Path
locals().update(config)
UTILS_WORKDIR = config.get('UTILS_WORKDIR', '/tscc/lustre/ddn/scratch/bay001')

tabledir = Path('/tscc/nfs/home/hsher/ps-yeolab5/ENCODE_paper_tables/')
metrics=pd.read_csv(tabledir/'model_performance.csv')
selected_models = metrics.loc[metrics['selected'], 'Experiment'].sort_values().tolist()

workdir: f"{UTILS_WORKDIR}/CHD_analysis"
VEP_CACHEDIR = config.get('UTILS_VEP_CACHEDIR', '/tscc/lustre/ddn/scratch/bay001/vep_cache')
UTILS_DIR = workflow.basedir
model_dir = Path('/tscc/nfs/home/hsher/ps-yeolab5/rbpnet_models')
model_dict = {f.name: f for f in model_dir.glob('*')}


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
    resources:
        mem_mb=40000,
        runtime=40
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
    resources:
        mem_mb=2000,
        runtime=40
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    shell:
        '''
        bedtools getfasta -fo {output.finemapped_fa} -fi {GENOME} -bed {input.finemapped_windows} -s
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
    resources:
        mem_mb=80000,
        runtime=80
    params:
        out_prefix = lambda wildcards, output: output.csv.replace('.csv', ''),
    conda:
        "../rules/envs/metadensity.yaml"
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
    resources:
        mem_mb=80000,
        runtime=20
    params:
        model_path = lambda wildcards: model_dict[wildcards.model_name]
    container:
        "/tscc/nfs/home/bay001/bin/eugene-tools_0.1.2.sif"
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
    resources:
        mem_mb=80000,
        runtime=20
    shell:
        """
        bcftools convert -c CHROM,POS,REF,ALT,-,-,filter,- -f {GENOME} --tsv2vcf {input} -o {output}
        """
rule download_vep_cache:
    output:
        Path(VEP_CACHEDIR) / f'homo_sapiens/{VEP_CACHE_VERSION}_GRCh38/1/all_vars.gz'
    threads: 1
    resources:
        mem_mb=2000,
        runtime="3h"
    container:
        f"docker://ensemblorg/ensembl-vep:release_{VEP_CACHE_VERSION}.4"
    shell:
        """
        mkdir -p {VEP_CACHEDIR}
        cd {VEP_CACHEDIR}
        wget https://ftp.ensembl.org/pub/release-{VEP_CACHE_VERSION}/variation/indexed_vep_cache/homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz
        tar xzf homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz
        rm homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz
        """

rule vep:
    input:
        vcf = "output/variants/variants.vcf",
        cache = rules.download_vep_cache.output,
    output:
        "output/variants/vep.tsv"
    threads: 2
    resources:
        mem_mb=40000,
        runtime=20
    container:
        f"docker://ensemblorg/ensembl-vep:release_{VEP_CACHE_VERSION}.4"
    shell:
        """
        vep \
        -i {input.vcf} \
        --force_overwrite \
        -o {output} --offline --cache --dir_cache {VEP_CACHEDIR} --cache_version {VEP_CACHE_VERSION}
        """

