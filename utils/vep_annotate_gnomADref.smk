
ROULETTE = '/tscc/nfs/home/hsher/ps-yeolab5/roulette/22_rate_v5.2_TFBS_correction_all.vcf.bgz'
HEADER = '/tscc/nfs/home/hsher/bin/Roulette/header.hr'
N_SPLIT=1000
JVARKIT_DIR='/tscc/nfs/home/hsher/bin/jvarkit/JVARKIT'
GENCODE='/tscc/nfs/home/hsher/gencode_coords/gencode.v40.primary_assembly.annotation.gff3'
MODEL_DIR='/tscc/nfs/home/hsher/ps-yeolab5/roulette/model'

locals().update(config)

workdir: '/tscc/nfs/home/hsher/scratch/annotate_roulette/'
"""
snakemake -s utils/vep_annotate_gnomADref.smk \
    --profile /tscc/nfs/home/hsher/projects/skipper/profiles/tscc2 \
    -n
"""

rule all:
    input:
        "22_rate_v5.2_TFBS_correction_all.vep.tsv",
        '22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.filtered.tsv.gz',
        MODEL_DIR+'/singleton.pickle'

rule split_vcf:
    input:
        Path(ROULETTE_DIR)/'22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz'
    output:
        temp(expand("output/{prefix}.{chunk}.vcf.gz", 
        chunk = [f"{i:05d}" for i in range(1, N_SPLIT+1)],
        prefix="{prefix}")),
        temp(expand("{prefix}.{chunk}.vcf.gz.tbi", 
       chunk = [f"{i:05d}" for i in range(1, N_SPLIT+1)],
       prefix="{prefix}"))
    threads: 1
    resources:
        mem_mb=80000,
        runtime="6h"
    conda:
        "envs/java.yaml"
    shell:
        """
        java -jar {JVARKIT_DIR}/jvarkit.jar vcfsplitnvariants \
             --vcf-count {N_SPLIT} \
             --prefix output/{wildcards.prefix} \
             {input}
        """

rule vep:
    input:
        "output/{prefix}.{chunk}.vcf.gz"
    output:
        temp("output/{prefix}.{chunk}.tsv"),
        temp("output/{prefix}.{chunk}.tsv_summary.html")
    threads: 2
    resources:
        mem_mb=40000,
        runtime="1h"
    container:
        "docker://ensemblorg/ensembl-vep:latest"
    shell:
        """
        vep \
        -i {input} \
        --force_overwrite \
        -o {output} -offline --cache {VEP_CACHEDIR}
        """

rule combine_vep:
    input:
        expand("output/{prefix}.{chunk}.tsv", 
       chunk = [f"{i:05d}" for i in range(1, N_SPLIT+1)],
       prefix="{prefix}")
    output:
        "{prefix}_all.vep.tsv"
    threads: 2
    resources:
        mem_mb=40000,
        runtime="1h"
    shell:
        """
        cat {input} | grep -v '#' > {output}
        """

rule annotate_RBP_site:
    input:
        vcf=Path(ROULETTE_DIR)/'22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz'
    output:
        allrbp="allrbpsite.bed.gz",
        allrbp_rename="allrbpsite.bed.rename.gz",
        vcf=Path(ROULETTE_DIR)/'22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.filtered.vcf.gz'
    threads: 1
    resources:
        mem_mb=40000,
        runtime="1h"
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        cat /tscc/nfs/home/hsher/scratch/ENCO*_*/output/finemapping/mapped_sites/*finemapped_windows.bed.gz > {output.allrbp}
        zcat {output.allrbp} | \
            sort -k1,1 -k2,2n -k3,3n | \
            awk -v OFS=\"\t\" '{{print $1, $2, $3, "RBP_site", $5, $6}}' | \
            bgzip > {output.allrbp_rename}
        tabix -p bed {output.allrbp_rename}

        bcftools annotate \
            -a {output.allrbp_rename} \
            -c CHROM,FROM,TO,RBPSITE \
            -h <(echo '##INFO=<ID=RBPSITE,Number=1,Type=String,Description="RBP site">') \
            -Oz -o {output.vcf} \
            {input.vcf}
        
        bcftools index {output.vcf}
        """

rule make_tsv:
    input:
        vcf=Path(ROULETTE_DIR)/'22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.filtered.vcf.gz'
    output:
        tsv='22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.filtered.tsv.gz'
    threads: 1
    resources:
        mem_mb=40000,
        runtime="1h"
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/RBPSITE\t%INFO/MR\t%INFO/AC\t%INFO/AN\n' \
            {input.vcf} \
            | gzip > {output.tsv} 
        """

rule fit_expected_singleton_rate:
    input:
        vcf=Path(ROULETTE_DIR)/'22_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.filtered.vcf.gz',
        vep='/tscc/nfs/home/hsher/scratch/annotate_roulette/22_rate_v5.2_TFBS_correction_all.vep.tsv'
    output:
        expected=MODEL_DIR+'/singleton.pickle'
    threads: 1
    resources:
        mem_mb=40000,
        runtime="1h"
    conda:
        "/tscc/nfs/home/hsher/projects/skipper/rules/envs/metadensity.yaml"
    shell:
        """
        python /tscc/nfs/home/hsher/projects/skipper/utils/expected_singleton_rate_regression.py \
            {input.vcf} \
            {input.vep} \
            {MODEL_DIR}
        """