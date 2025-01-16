from pathlib import Path
ROULETTE_DIR=Path('/tscc/nfs/home/hsher/ps-yeolab5/roulette/')
VCF='/tscc/projects/ps-yeolab5/hsher/clinvar/clinvar.rename.vcf.gz'


# Defines regions to query
TABLE='/tscc/projects/ps-yeolab4/software/skipper/bb63a25/bin/skipper/annotations/gencode.v41.annotation.tiled_partition.bed.gz'
VEP_CACHEDIR='/tscc/nfs/home/hsher/scratch/vep_cache/'
import pandas as pd
locals().update(config)

tabledir = Path('/tscc/nfs/home/hsher/ps-yeolab5/ENCODE_paper_tables/')
metrics=pd.read_csv(tabledir/'model_performance.csv')
selected_models = metrics.loc[metrics['selected'], 'Experiment'].sort_values().tolist()

workdir: "/tscc/nfs/home/hsher/scratch/clinvar_analysis"
from pathlib import Path
model_dir = Path('/tscc/nfs/home/hsher/scratch/')
models = model_dir.glob('ENCO*/output/ml/rbpnet_model/*')
eclip_dict = {}
model_dict = {}
for f in models:
    model_dict[f.name] = f
    eclip_dict[f.name] = str(f).split('output')[0] # point to eclip output


"""
snakemake -s utils/score_all_clinvar.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2 \
    -n
"""

rule all:
    input:
        expand("output/variants/{model_name}.score.csv.gz",
        model_name = selected_models
        ),
        "output/variants/vep.tsv",
        expand("output/analysis/{model_name}.clinvar_variants.csv",
        model_name = selected_models
        ),
        expand("output/mondo/{model_name}.clinvar_mondo_tagged.csv",
        model_name = selected_models
        ),
        # expand("output/variants/{model_name}.LoF_oe.csv",
        # model_name = selected_models
        # ),
        #"output/variants/vep.tsv"

rule fetch_SNP:
    ''' fetch gnomAD variants from database '''
    input:
        vcf=VCF,
        finemapped_windows = TABLE
    output:
        "output/variants/fetched.vcf"
    threads: 2
    params:
        error_file = "stderr/fetch_snp",
        out_file = "stdout/fetch_snp",
        run_time = "13:20:00",
        cores = 1,
        memory = 40000,
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        if [ -s {input.finemapped_windows} ]; then
            bcftools query -R {input.finemapped_windows} \
                -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNDN\t%INFO/CLNVC\t%INFO/CLNSIG\t%INFO/CLNDISDB\t%INFO/AF_ESP\t%INFO/AF_EXAC\t%INFO/AF_TGP\t%INFO/ALLELEID\n' \
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
        subset_vcf="output/variants/fetched.vcf",
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
        score="output/variants/{model_name}.score.csv",
    threads: 1
    params:
        error_file = "stderr/score_fa.{model_name}",
        out_file = "stdout/score_fa.{model_name}",
        run_time = "02:00:00",
        cores = 1,
        memory = 80000,
        model_path = lambda wildcards: model_dict[wildcards.model_name]
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
#%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNDN\t%INFO/CLNVC\t%INFO/CLNSIG\t%INFO/CLNDISDB\t%INFO/AF_ESP\t%INFO/AF_EXAC\t%INFO/AF_TGP\t%INFO/ALLELEID\n
rule combine_score:
    input:
        score = "output/variants/{model_name}.score.csv",
        vcf = "output/variants/fetched.vcf"
    output:
        "output/variants/{model_name}.score.csv.gz"
    threads: 1
    params:
        error_file = "stderr/combine_score.{model_name}",
        out_file = "stdout/combine_score.{model_name}",
        run_time = "00:20:00",
        cores = 1,
        memory = 80000,
    run:
        print(input.vcf)
        variants = pd.read_csv(input.vcf, 
                        sep = '\t', names =['CHROM','POS','ID','REF','ALT',
                        'INFO/CLNDN','INFO/CLNVC','INFO/CLNSIG','INFO/CLNDSIGDB','INFO/AF_ESP','INFO/AF_EXAC',
                        'INFO/AF_TGP','INFO/ALLELEID'
                        ],
            na_values = '.')
        print(variants.columns)
        score = pd.read_csv(input.score,
                    index_col = 0)
        score[['CHROM', 'POS', 'NU', 'name']]=score['ID'].str.split('-', expand = True)
        score['POS']=score['POS'].astype(int)
        ref_score = score.merge(variants, left_on = ['CHROM', 'POS', 'NU'],
                right_on = ['CHROM', 'POS', 'REF'],
                        how = 'right')[['CHROM', 'POS', 'REF', 'ALT', 'dlogodds_pred']]
        alt_score = score.merge(variants, left_on = ['CHROM', 'POS', 'NU'],
                    right_on = ['CHROM', 'POS', 'ALT'],
                    suffixes = ('_score', ''),
                            how = 'right')[['CHROM', 'POS', 'ID','REF', 'ALT', 'dlogodds_pred',
                            'INFO/CLNDN','INFO/CLNVC','INFO/CLNSIG','INFO/CLNDSDB','INFO/AF_ESP','INFO/AF_EXAC',
                        'INFO/AF_TGP','INFO/ALLELEID', 'name']]
        all_scores = ref_score.merge(alt_score, 
                                    left_on = ['CHROM', 'POS', 'REF', 'ALT'],
                                    right_on = ['CHROM', 'POS','REF', 'ALT'],
                                    suffixes = ('_REF', '_ALT')
                                )
        all_scores['delta_score'] = all_scores['dlogodds_pred_ALT']-all_scores['dlogodds_pred_REF']
        all_scores.to_csv(output[0], index = False, compression = 'gzip')
    

rule vep:
    input:
        "output/variants/fetched.vcf"
    output:
        "output/variants/vep.tsv"
    threads: 2
    params:
        error_file = "stderr/vep",
        out_file = "stdout/vep",
        run_time = "6:20:00",
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

rule analysis:
    input:
        score="output/variants/{model_name}.score.csv.gz",
        gnomad=lambda wildcards: Path(eclip_dict[wildcards.model_name]) / f'output/variants/gnomAD_roulette/{wildcards.model_name}.total.csv',
        vep = "output/variants/vep.tsv",
        gnomadech=expand("output/gnomAD_roulette/{chr_number}.vcf",
        chr_number = range(1, 23)
        )
    output:
        "output/analysis/{model_name}.clinvar_variants.csv",
        "output/analysis/{model_name}.clinvar_CLINSIC_counts.csv",
        "output/analysis/{model_name}.clinvar_impact_enrichment.csv",
        "output/analysis/{model_name}.clinvar_impact_counts.csv",
    threads: 1
    params:
        error_file = "stderr/analysis.{model_name}",
        out_file = "stdout/analysis.{model_name}",
        run_time = "01:20:00",
        cores = 1,
        memory = 80000,
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/metadensity.yaml"
    shell:
        """
        python /tscc/nfs/home/hsher/projects/skipper/utils/clivnar_analysis.py \
            {input.score} \
            output/gnomAD_roulette/ \
            {input.gnomad} \
            {input.vep} \
            output/analysis/
        """

rule fetch_SNP_from_gnomAD_and_roulette:
    ''' fetch gnomAD variants from database '''
    input:
        vcf=ROULETTE_DIR/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz',
        finemapped_windows = TABLE
    output:
        temp("output/gnomAD_roulette/{chr_number}.vcf")
    threads: 2
    params:
        error_file = "stderr/fetch_snp.{chr_number}",
        out_file = "stdout/fetch_snp.{chr_number}",
        run_time = "13:20:00",
        cores = 1,
        memory = 40000,
    # container:
    #     "docker://brianyee/bcftools:1.17"
    shell:
        """
        module load bcftools
        if [ -s {input.finemapped_windows} ]; then
            bcftools query -R {input.finemapped_windows} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' \
                {input.vcf} > {output}
        else
            touch {output}
        fi
        """

rule aggregate_mondo:
    input:
        variants="output/analysis/{model_name}.clinvar_variants.csv",
        href="/tscc/nfs/home/hsher/ontology/mondo/mondo_href.csv"
    output:
        "output/mondo/{model_name}.clinvar_mondo_tagged.csv",
        "output/mondo/{model_name}.clinvar_mondo_cnts.csv",
        "output/mondo/{model_name}.clinvar_mondo_stats.csv",
    threads: 1
    params:
        error_file = "stderr/mondo.{model_name}",
        out_file = "stdout/mondo.{model_name}",
        run_time = "01:20:00",
        cores = 1,
        memory = 80000,
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/mondo.yaml"
    shell:
        """
        python /tscc/nfs/home/hsher/projects/skipper/utils/propagate_mondo.py \
            {input.variants} \
            output/mondo/{wildcards.model_name}
        """
