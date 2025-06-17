
import pandas as pd
from pathlib import Path
locals().update(config)
VEP_CACHEDIR = config['VEP_CACHEDIR']
VEP_CACHE_VERSION = config['VEP_CACHE_VERSION']

rule filter_roulette_for_high:
    input:
        vcf=Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.vcf.bgz',
        header=HEADER
    output:
        reheader = temp(Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.vcf'),
        filtered = temp(Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.vcf'),
    container:
        "docker://brianyee/bcftools:1.17"
    threads: 1
    resources:
        mem_mb=40000,
        runtime="1h"
    shell:
        """
        bcftools reheader -h {input.header} \
            {input.vcf} > {output.reheader}
        bcftools filter -O z -o {output.filtered} -i 'FILTER=="high"' {output.reheader}
        """

rule annotate_roulette_w_gnomAD:
    input:
        rename=RENAME,
        gnomad=Path(GNOMAD_DIR) / 'gnomad.genomes.v4.1.sites.chr{chr_number}.vcf.bgz',
        roulette=Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.vcf'
    output:
        rename=temp(Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.rename.vcf.gz'),
        annotated=Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz'
    threads: 1
    resources:
        mem_mb=80000,
        runtime="6h"
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        bcftools annotate --rename-chrs {input.rename} \
            {input.roulette} -Oz -o {output.rename}
        bcftools index {output.rename}
        bcftools annotate -a {input.gnomad} \
            -c INFO/AC,INFO/AN \
            -k {output.rename} \
            -o {output.annotated}
        bcftools index {output.annotated}
        """

rule fetch_SNP_from_gnomAD_and_roulette:
    ''' fetch gnomAD variants from database '''
    input:
        vcf=Path(ROULETTE_DIR)/'{chr_number}_rate_v5.2_TFBS_correction_all.header.filtered.rename.annotated.vcf.gz',
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    output:
        "output/variants/gnomAD_roulette/{experiment_label}.chr{chr_number}.vcf"
    threads: 2
    resources:
        runtime=lambda wildcards, attempt: 480 * (1.5 ** (attempt - 1)),
        mem_mb=40000,
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
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    output:
        "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.slop.bed.gz"
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
        finemapped_fa = "output/ml/sequence/{experiment_label}.foreground.slop.fa",
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
        subset_vcf="output/variants/{subset}/{experiment_label_thing}.vcf",
        seq_fa = lambda wildcards: f"output/ml/sequence/{wildcards.experiment_label_thing}.foreground.slop.fa" 
            if 'chr' not in wildcards.experiment_label_thing
            else "output/ml/sequence/"+wildcards.experiment_label_thing.split('.')[0]+".foreground.slop.fa",
        finemapped_windows = lambda wildcards: f"output/finemapping/mapped_sites/{wildcards.experiment_label_thing}.finemapped_windows.slop.bed.gz" 
            if 'chr' not in wildcards.experiment_label_thing
            else "output/finemapping/mapped_sites/"+wildcards.experiment_label_thing.split('.')[0]+".finemapped_windows.slop.bed.gz"
    output:
        ref_fa = temp("output/variants/{subset}/{experiment_label_thing}.ref.fa"),
        alt_fa = temp("output/variants/{subset}/{experiment_label_thing}.alt.fa"),
        csv = "output/variants/{subset}/{experiment_label_thing}.csv"
    threads: 2
    resources:
        mem_mb=80000,
        runtime="1h"
    params:
        out_prefix = lambda wildcards, output: output.csv.replace('.csv', ''),
        tmpdir = config["TMPDIR"],
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        module purge;
        export TMPDIR={params.tmpdir};
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
        model=lambda wildcards: "output/ml/rbpnet_model/{experiment_label_thing}/training_done" if 'chr' not in wildcards.experiment_label_thing 
            else "output/ml/rbpnet_model/"+wildcards.experiment_label_thing.split('.')[0]+"/training_done",
        fa = "output/variants/{subset}/{experiment_label_thing}.{type}.fa"
    output:
        score=temp("output/variants/{subset}/{experiment_label_thing}.{type}.score.csv"),
        fai=temp("output/variants/{subset}/{experiment_label_thing}.{type}.fa.fai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 25000 * (2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 180 * (2 ** (attempt - 1)),
    params:
        exp =lambda wildcards: wildcards.experiment_label_thing.split('.')[0],
    container:
        "/tscc/nfs/home/bay001/eugene-tools_0.1.2.sif"
    shell:
        """
        export NUMBA_CACHE_DIR=/tscc/lustre/ddn/scratch/${{USER}} # TODO: HARCODED IS BAD
        export MPLCONFIGDIR=/tscc/lustre/ddn/scratch/${{USER}}
        
        if [ -s {input.fa} ]; then
            samtools faidx {input.fa};
            python {RBPNET_PATH}/score_fa.py \
                output/ml/rbpnet_model/{params.exp}/ \
                {input.fa} \
                {output.score} \
                /tscc/lustre/ddn/scratch/${{USER}}
        else
            touch {output.score}
            touch {output.fai}
        fi
        """

rule join_gnomAD_info:
    input:
        scores = expand("output/variants/gnomAD_roulette/{experiment_label}.chr{chr_number}.{type}.score.csv",
            experiment_label = ["{experiment_label}"],
            chr_number = list(range(1,23)),
            type = ["ref", "alt"]),
        vcf = expand("output/variants/gnomAD_roulette/{experiment_label}.chr{chr_number}.vcf",
            experiment_label = ["{experiment_label}"],
            chr_number = list(range(1,23))),
    output:
        "output/variants/gnomAD_roulette/{experiment_label}.total.csv"
    threads: 1
    resources:
        mem_mb=80000,
        runtime=20
    run:
        indir = Path('output/variants/gnomAD_roulette/')
        scores = []
        exp = wildcards.experiment_label
        for chr in range(1,23):
            if os.stat(indir / f'{exp}.chr{chr}.vcf').st_size == 0:
                continue
            print(f'handling chrom{chr}')
            alt_score = pd.read_csv(indir / f'{exp}.chr{chr}.alt.score.csv', index_col = 0)
            ref_score = pd.read_csv(indir / f'{exp}.chr{chr}.ref.score.csv', index_col = 0)

            
            vcf = pd.read_csv(indir / f'{exp}.chr{chr}.vcf', sep = '\t',
                            names = ['CHROM', 'POS', '.', 'REF', 'ALT',
                                    'INFO/AC', 'INFO/AN', 'INFO/MR', 'INFO/AR',
                                    'INFO/MG', 'INFO/MC'])

            ref_score[['CHROM', 'POS', 'REF', 'name']]=ref_score['ID'].str.split('-', expand = True)
            alt_score[['CHROM', 'POS', 'ALT', 'name']]=alt_score['ID'].str.split('-', expand = True)
            score = alt_score.drop('ID', axis = 1).merge(ref_score.drop('ID', axis = 1), 
                            left_on = ['CHROM', 'POS', 'name'],
                            right_on = ['CHROM', 'POS', 'name'],
                            suffixes = ('_ALT', '_REF')
                            )
            score['delta_score'] = score['dlogodds_pred_ALT']-score['dlogodds_pred_REF']
            score['POS'] = score['POS'].astype(int)
            score_m = score.merge(vcf.drop('.', axis = 1), left_on = ['CHROM', 'POS', 'REF', 'ALT'],
                        right_on = ['CHROM', 'POS', 'REF', 'ALT']
                    )
            scores.append(score_m)
        try:
            scores = pd.concat(scores, axis = 0)
            scores.to_csv(output[0], index = False)
        except:
            print('no gnomAD variants found')
            open(output[0], 'w').close()

rule fetch_Clinvar_SNP:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        vcf = CLINVAR_VCF.replace('.vcf.gz', '.rename.vcf.gz')
    output:
        "output/variants/clinvar/{experiment_label}.vcf"
    threads: 2
    resources:
        mem_mb=50000,
        runtime="3h"
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
        mkdir -p {VEP_CACHEDIR}/ && cd {VEP_CACHEDIR}/;        
        curl -o {VEP_CACHEDIR}/homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz;
        tar xzf homo_sapiens_vep_{VEP_CACHE_VERSION}_GRCh38.tar.gz
        """
rule vep:
    input:
        vcf = "output/variants/clinvar/{experiment_label}.vcf",
        vars = Path(VEP_CACHEDIR) / f'homo_sapiens/{VEP_CACHE_VERSION}_GRCh38/1/all_vars.gz'
    output:
        "output/variants/clinvar/{experiment_label}.vep.tsv"
    threads: 2
    resources:
        mem_mb=40000,
        runtime="1h"
    container:
        "docker://ensemblorg/ensembl-vep:release_113.4"
    shell:
        """
        if [ -s {input.vcf} ]; then
            vep \
            -i {input.vcf} \
            --force_overwrite \
            -o {output} -offline --cache {VEP_CACHEDIR}
        else
            touch {output}
        fi
        """

rule variant_analysis:
    input:
        clinvar = "output/variants/clinvar/{experiment_label}.vep.tsv",
        gnomAD = "output/variants/gnomAD_roulette/{experiment_label}.total.csv",
        annotated = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv"
    output:
        "output/variant_analysis/{experiment_label}.clinvar_variants.csv",
        "output/variant_analysis/{experiment_label}.annotated.csv.gz",
        "output/variant_analysis/{experiment_label}.feature_type_top.oe_stat.lofbins.csv",
        "output/variant_analysis/{experiment_label}.transcript_type_top.MAPS_stat.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top_spectrum_enrichment.csv",
        "output/variant_analysis/{experiment_label}.transcript_type_top.MAPS.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.oe.csv",
        "output/variant_analysis/{experiment_label}.transcript_type_top_spectrum_enrichment.csv",
        "output/variant_analysis/{experiment_label}.transcript_type_top.oe_stat.csv",
        "output/variant_analysis/{experiment_label}.transcript_type_top.oe.csv",
        "output/variant_analysis/{experiment_label}.clinvar_CLINSIC_counts.csv",
        "output/variant_analysis/{experiment_label}.global_MAPS.csv",
        "output/variant_analysis/{experiment_label}.global_spectrum_enrichment.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.MAPS.lofbins.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.oe_stat.csv",
        "output/variant_analysis/{experiment_label}.clinvar_impact_counts.pdf",
        "output/variant_analysis/{experiment_label}.clinvar_CLINSIC_counts.pdf",
        "output/variant_analysis/{experiment_label}.global_MAPS_stat.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.MAPS_stat.lofbins.csv",
        "output/variant_analysis/{experiment_label}.global_oe.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.MAPS_stat.csv",
        "output/variant_analysis/{experiment_label}.clinvar_impact_counts.csv",
        "output/variant_analysis/{experiment_label}.clinvar_variants_exploded.csv",
        "output/variant_analysis/{experiment_label}.global_oe_stat.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.MAPS.csv",
        "output/variant_analysis/{experiment_label}.feature_type_top.oe.lofbins.csv",
    threads: 1
    resources:
        mem_mb=40000,
        runtime=20
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        module purge;
        if [ -s {input.gnomAD} ]; then
            python {TOOL_DIR}/mega_variant_analysis.py \
                . \
                {wildcards.experiment_label} \
                {SINGLETON_REFERENCE} \
                {OE_RATIO_REFERENCE} \
                {GNOMAD_CONSTRAINT} \

        else
            touch {output}
        fi
        """

