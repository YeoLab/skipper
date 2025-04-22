from pathlib import Path
VCFs = ROULETTE_DIR.glob('*rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz')

# Defines regions to query
TABLE='/tscc/nfs/home/hsher/bin/poison_exon_variant_analysis/data/Felker2023SuppelementaryTable2_hg38_PE_cassettes.strand.bed'
import pandas as pd
locals().update(config)

tabledir = Path('/tscc/nfs/home/hsher/ps-yeolab5/ENCODE_paper_tables/')
metrics=pd.read_csv(tabledir/'model_performance.csv')
selected_models = metrics.loc[metrics['selected'], 'Experiment'].sort_values().tolist()

workdir: "/tscc/nfs/home/hsher/scratch/poison_exon_analysis"
from pathlib import Path
model_dir = Path('/tscc/nfs/home/hsher/scratch/')
models = model_dir.glob('ENCO*/output/ml/rbpnet_model/*')
eclip_dict = {}
model_dict = {}
for f in models:
    model_dict[f.name] = f
    eclip_dict[f.name] = str(f).split('output')[0] # point to eclip output


print(model_dict)
"""
snakemake -s utils/score_any_vcf.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2_gpu \
    --until score_fa
"""

rule all:
    input:
        expand("output/variants/{model_name}.score.csv.gz",
        model_name = selected_models
        ),
        expand("output/variants/{model_name}.LoF_oe.csv",
        model_name = selected_models
        ),
        expand("output/wt_score/{model_name}.score",
        model_name = selected_models
        ),
        #"output/variants/vep.tsv"

rule fetch_SNP_from_gnomAD_and_roulette:
    ''' fetch gnomAD variants from database '''
    input:
        vcf=Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz',
        finemapped_windows = TABLE
    output:
        "output/variants/{chr_number}.vcf"
    threads: 2
    resources:
        mem_mb=40000,
        runtime="13h"
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        if [ -s {input.finemapped_windows} ]; then
            bcftools query -R {input.finemapped_windows} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' \
                {input.vcf} > {output}
        else
            touch {output}
        fi
        """

rule slop_finemap:
    input:
        finemapped_windows = TABLE
    output:
        "output/table.slop.bed.gz"
    threads: 2
    resources:
        mem_mb=40000,
        runtime="6h"
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

rule score_wt:
    input:
        "output/sequence/table.slop.fa"
    output:
        "output/wt_score/{model_name}.score"
    threads: 1
    resources:
        mem_mb=80000,
        runtime="2h"
    params:
        model_path = lambda wildcards: model_dict[wildcards.model_name]
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}

        python /tscc/nfs/home/hsher/projects/RBPNet/score_fa.py \
            {params.model_path} \
            {input}  {output} \
            /tscc/lustre/ddn/scratch/${{USER}}
        """



rule fetch_variant_sequence:
    input:
        subset_vcf="output/variants/{chr_number}.vcf",
        seq_fa = rules.fetch_peak_sequence.output.finemapped_fa,
        finemapped_windows = rules.slop_finemap.output
    output:
        ref_fa = temp("output/variants/annotation.{chr_number}.ref.fa"),
        alt_fa = temp("output/variants/annotation.{chr_number}.alt.fa"),
        all_fa = "output/variants/annotation.{chr_number}.all.fa",
        csv = "output/variants/annotation.{chr_number}.csv"
    threads: 2
    resources:
        mem_mb=80000,
        runtime="80:00"
    params:
        out_prefix = lambda wildcards, output: output.csv.replace('.csv', ''),
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
        fa = "output/variants/annotation.{chr_number}.all.fa"
    output:
        score="output/variants/{model_name}.{chr_number}.score.csv",
    threads: 1
    params:
        model_path = lambda wildcards: model_dict[wildcards.model_name]
    resources:
        mem_mb=80000,
        runtime="2h"
    container: 
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
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

rule combine_score:
    input:
        expand("output/variants/{model_name}.{chr_number}.score.csv",
        model_name = "{model_name}",
        chr_number = range(1, 23)
        ),
    output:
        "output/variants/{model_name}.score.csv.gz"
    threads: 1
    resources:
        mem_mb=80000,
        runtime="20:00"
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/metadensity.yaml"
    shell:
        """
        python /tscc/nfs/home/hsher/projects/skipper/utils/join_scores.py \
            {wildcards.model_name} output/variants \
            {output}
        """
    
rule analysis:
    input:
        score="output/variants/{model_name}.score.csv.gz",
        gnomad=lambda wildcards: Path(eclip_dict[wildcards.model_name]) / f'output/variants/gnomAD_roulette/{wildcards.model_name}.total.csv'
    output:
        "output/variants/{model_name}.global_oe.csv",
        "output/variants/{model_name}.global_maps.csv",
        "output/variants/{model_name}.MIM_oe.csv",
        "output/variants/{model_name}.SFARI_oe.csv",
        "output/variants/{model_name}.LoF_oe.csv"
    threads: 1
    resources:
        mem_mb=80000,
        runtime="20:00"
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/metadensity.yaml"
    shell:
        """
        python /tscc/nfs/home/hsher/projects/skipper/utils/gnomad_analysis_poison_exon.py \
            {input.score} \
            {input.gnomad} \
            output/variants/{wildcards.model_name}
        """
rule vep:
    input:
        "output/variants/variants.vcf"
    output:
        "output/variants/vep.tsv"
    threads: 2
    resources:
        mem_mb=40000,
        runtime="80:00"
    container:
        "docker://ensemblorg/ensembl-vep:latest"
    shell:
        """
        vep \
        -i {input} \
        --force_overwrite \
        -o {output} -offline --cache {VEP_CACHEDIR}
        """

