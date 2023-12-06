import pandas as pd
locals().update(config)
# fetch GnomAD, ClinVar and cancer somatic mutations in RBP binding sites and prepare sequence for ML inference
CLINVAR_VCF='/projects/ps-yeolab5/hsher/clinvar/clinvar.vcf.gz'
COSMIC_CODING_VCF='/projects/ps-yeolab5/hsher/cosmic_data/Cosmic_GenomeScreensMutant_Normal_v98_GRCh38.vcf.gz'
COSMIC_NONCODING_VCF='/projects/ps-yeolab5/hsher/cosmic_data/Cosmic_NonCodingVariants_Normal_v98_GRCh38.vcf.gz'

rule fetch_gnomAD_SNP:
    ''' fetch gnomAD variants from database '''
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    output:
        temp("output/variants/gnomAD/{experiment_label}.{chr}.vcf")
    params:
        error_file = "stderr/fetch_snp.{experiment_label}.{chr}",
        out_file = "stdout/fetch_snp.{experiment_label}.{chr}",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools
        bcftools query -R {input.finemapped_windows} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
            https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.{wildcards.chr}.vcf.bgz > {output}
        """

rule combine_gnomAD_vcf:
    input:
        expand(rules.fetch_gnomAD_SNP.output,
            chr = [f'chr{i}' for i in list(range(1,23))+['X','Y']],
            experiment_label = ['{experiment_label}'])
    output:
        "output/variants/gnomAD/{experiment_label}.vcf"
    params:
        error_file = "stderr/combine_snp.{experiment_label}",
        out_file = "stdout/combine_snp.{experiment_label}",
        run_time = "20:00",
        cores = 1,
    shell:
        """
        cat {input} > {output}
        """
    

rule reannotate_vcf:
# contig 1 -> chr1
    input:
        vcf="{anything}.vcf.gz",
        rename="/home/hsher/projects/oligoCLIP/utils/rename_chr.txt"
    output:
        "{anything}.rename.vcf.gz",
    params:
        error_file = "stderr/rename_chr",
        out_file = "stdout/rename_chr",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools
        bcftools annotate --rename-chrs {input.rename} \
            {input.vcf} -Oz -o {output}
        bcftools index {output}
        """



rule fetch_Clinvar_SNP:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        vcf = CLINVAR_VCF.replace('.vcf.gz', '.rename.vcf.gz')
    output:
        "output/variants/clinvar/{experiment_label}.vcf"
    params:
        error_file = "stderr/fetch_clinvar_snp.{experiment_label}",
        out_file = "stdout/fetch_clinvar_snp.{experiment_label}",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools 
        bcftools query -R {input.finemapped_windows} \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNDN\t%INFO/CLNVC\t%INFO/CLNSIG\t%INFO/CLNDISDB\t%INFO/AF_ESP\t%INFO/AF_EXAC\t%INFO/AF_TGP\t%INFO/ALLELEID\n' \
            {input.vcf} > {output}
        """


rule fetch_COSMIC_SNP:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        vcf = COSMIC_CODING_VCF.replace('.vcf.gz', '.rename.vcf.gz')
    output:
        "output/variants/cosmic_coding/{experiment_label}.vcf"
    params:
        error_file = "stderr/fetch_cosmic_snp.{experiment_label}",
        out_file = "stdout/fetch_cosmic_snp.{experiment_label}",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools 
        bcftools query -R {input.finemapped_windows} \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AA\t%INFO/SAMPLE_COUNT\t%INFO/TIER\n' \
            {input.vcf} > {output}
        """

rule fetch_COSMIC_NONCODING_SNP:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        vcf = COSMIC_NONCODING_VCF.replace('.vcf.gz', '.rename.vcf.gz')
    output:
        "output/variants/cosmic_noncoding/{experiment_label}.vcf"
    params:
        error_file = "stderr/fetch_cosmic_snp.{experiment_label}",
        out_file = "stdout/fetch_cosmic_snp.{experiment_label}",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools 
        bcftools query -R {input.finemapped_windows} \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SAMPLE_COUNT\t%INFO/TIER\n' \
            {input.vcf} > {output}
        """


rule fetch_sequence:
    input:
        subset_vcf="output/variants/{subset}/{experiment_label}.vcf",
        seq_fa = "output/ml/sequence/{experiment_label}.foreground.fa",
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    output:
        ref_fa = "variants/{subset}/{experiment_label}.ref.fa",
        alt_fa = "variants/{subset}/{experiment_label}.alt.fa",
        csv = "variants/{subset}/{experiment_label}.csv"
    params:
        error_file = "stderr/fetch_sequence.{subset}.{experiment_label}",
        out_file = "stdout/fetch_sequence.{subset}.{experiment_label}",
        run_time = "03:20:00",
        cores = 1,
        out_prefix = lambda wildcards, output: output.csv.replace('.csv', '')
    conda:
        "/home/hsher/projects/oligoCLIP/rules/envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/generate_variant_sequence.py \
            {input.subset_vcf} \
            {input.seq_fa} \
            {input.finemapped_windows} \
            {params.out_prefix}
        """
rule score_variants:
    input:
        ref_fa = rules.fetch_sequence.output.ref_fa,
        alt_fa = rules.fetch_sequence.output.alt_fa,
        model = "output/ml/gkmsvm/{experiment_label}.model.txt",
    output:
        ref_score = "variants/{subset}/{experiment_label}.ref.score.txt",
        alt_score = "variants/{subset}/{experiment_label}.alt.score.txt"
    params:
        error_file = "stderr/score_variants.{subset}.{experiment_label}",
        out_file = "stdout/score_variants.{subset}.{experiment_label}",
        run_time = "04:20:00",
        cores = 1,
    shell:
        """
        /home/hsher/bin/lsgkm-0.1.1/bin/gkmpredict {input.ref_fa} {input.model} {output.ref_score}
        /home/hsher/bin/lsgkm-0.1.1/bin/gkmpredict {input.alt_fa} {input.model} {output.alt_score}
        """

def find_well_trained_model(wildcards):
    auprc_threshold = 0.6
    
    auprc_df =  pd.read_csv(checkpoints.ml_gkmsvm_AUPRC.get().output[0], index_col = 0)
    print('RBP with good enough model to do variant interpretation:')
    print(auprc_df.loc[auprc_df['mean AUPRC']>auprc_threshold])
    with_good_model = auprc_df.loc[auprc_df['mean AUPRC']>auprc_threshold, 'Experiment'].tolist()
    
    return expand("variants/{variant_set}/{experiment_label}.{type}.score.txt",
        experiment_label = with_good_model, 
        variant_set = ['gnomAD', 'clinvar', 'cosmic_coding', 'cosmic_noncoding'], type = ['alt', 'ref'])
        

rule variants_done:
    input:
        find_well_trained_model
    output:
        "variants_done.txt"
    params:
        error_file = "stderr/variants_done",
        out_file = "stdout/variants_done",
        run_time = "04:20:00",
        cores = 1,
    shell:
        """
        touch {output}
        """