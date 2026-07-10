from pathlib import Path
import pandas as pd
locals().update(config)
workdir: "/tscc/nfs/home/hsher/scratch/gnomad_vep_store"
# VCF='/tscc/projects/ps-yeolab5/hsher/clinvar/clinvar.rename.vcf.gz'
# Defines regions to query
TABLE='/tscc/projects/ps-yeolab4/software/skipper/bb63a25/bin/skipper/annotations/gencode.v41.annotation.tiled_partition.bed.gz'


# load models
tabledir = Path('/tscc/nfs/home/hsher/projects/ENCODE/tables/')

"""
snakemake -s utils/annotate_entire_gnomAD_vep.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2_snakemake9 \
    -n
"""

rule all:
    input:
        expand("output/variants/vep/chr{chr_number}.vcf", chr_number=list(range(1,23)))

rule fetch_SNP_from_gnomAD_and_roulette:
    ''' fetch gnomAD variants from database '''
    input:
        vcf=Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz',
        finemapped_windows = TABLE,
    output:
        "output/variants/gnomAD_roulette/chr{chr_number}.vcf"
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
rule download_vep_cache:
    output:
        Path(VEP_CACHEDIR) / f'homo_sapiens/{VEP_CACHE_VERSION}_GRCh38/1/all_vars.gz'
    threads: 1
    resources:
        mem_mb=2000,
        runtime="1h"
    container:
        "docker://ensemblorg/ensembl-vep:release_113.4"
    shell:
        """
        cd {VEP_CACHEDIR}
        wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz
        tar xzf homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz
        """
rule vep:
    input:
        "output/variants/gnomAD_roulette/chr{chr_number}.vcf"
    output:
        "output/variants/vep/chr{chr_number}.vcf"
    threads: 2
    resources:
        mem_mb=40000,
        runtime="6h",
        cache= VEP_CACHEDIR
    container:
        "docker://ensemblorg/ensembl-vep:latest"
    shell:
        """
        vep \
        -i {input} \
        --force_overwrite \
        -o {output} -offline --cache {VEP_CACHEDIR}
        """

