import pandas as pd
locals().update(config)
# fetch GnomAD, ClinVar and cancer somatic mutations in RBP binding sites and prepare sequence for ML inference
CLINVAR_VCF='/tscc/projects/ps-yeolab5/hsher/clinvar/clinvar.vcf.gz'
COSMIC_CODING_VCF='/tscc/projects/ps-yeolab5/hsher/cosmic_data/Cosmic_GenomeScreensMutant_Normal_v98_GRCh38.vcf.gz'
COSMIC_NONCODING_VCF='/tscc/projects/ps-yeolab5/hsher/cosmic_data/Cosmic_NonCodingVariants_Normal_v98_GRCh38.vcf.gz'
ROULETTE_PATH='/tscc/projects/ps-yeolab5/hsher/roulette'

rule fetch_gnomAD_SNP:
    ''' fetch gnomAD variants from database '''
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz"
    output:
        "output/variants/gnomAD/{experiment_label}.{chr}.vcf"
    threads: 2
    params:
        error_file = "stderr/fetch_snp.{experiment_label}.{chr}",
        out_file = "stdout/fetch_snp.{experiment_label}.{chr}",
        run_time = "13:20:00",
        cores = 1,
        memory = 40000,
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        bcftools query -R {input.finemapped_windows} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
            https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.{wildcards.chr}.vcf.bgz > {output}
        """

# sometimes will miss 1 or 2 chromosomes... wtf TSCC?
rule combine_gnomAD_vcf:
    input:
        expand(rules.fetch_gnomAD_SNP.output,
            chr = [f'chr{i}' for i in list(range(1,23))+['X','Y']],
            experiment_label = ['{experiment_label}'])
    output:
        "output/variants/gnomAD/{experiment_label}.vcf"
    threads: 2
    params:
        error_file = "stderr/combine_snp.{experiment_label}",
        out_file = "stdout/combine_snp.{experiment_label}",
        run_time = "20:00",
        cores = 1,
        memory = "16000",
    resources:
        mem_mb=16000
    shell:
        """
        cat {input} > {output}
        """

# rule annotate_roulette_rate:
#     # for observed gnomAD mutations, annotate mutation rates
#     input: 
#         vcf = "output/variants/gnomAD/{experiment_label}.vcf"
#     output:
#         input_tsv=temp("output/variants/gnomAD/{experiment_label}.tsv"),
#         tsv = "output/variants/gnomAD/{experiment_label}.vcf.roulette.tsv"
#     threads: 1
#     params:
#         error_file = "stderr/annotate_roulette_rate.{experiment_label}",
#         out_file = "stdout/annotate_roulette_rate.{experiment_label}",
#         run_time = "3:20:00",
#         cores = 1,
#         memory = 60000,
#         roulette_path = ROULETTE_PATH
#     container:
#         "docker://brianyee/bcftools:1.17"
#     shell:
#         """
#         cut -f 1,2,4,5 {input.vcf} | awk 'BEGIN {{print "CHROM\tPOS\tREF\tALT"}} {{ gsub(/chr/,"", $1); print }} ' OFS=\"\t\" > {output.input_tsv}
#         bash /tscc/nfs/home/hsher/bin/Roulette/adding_mutation_rate/roulette_annotate.sh {output.input_tsv} {params.roulette_path}
#         """
# annotate ROulette rates
#cut -f 1,2,4,5  ~/scratch/ENCODE3_HepG2/output/variants/gnomAD/YBX3_HepG2_ENCSR735HOK.vcf | awk 'BEGIN {print "CHROM\tPOS\tREF\tALT"} { gsub(/chr/,"", $1); print } ' OFS='\t' > ~/scratch/YBX3_test.tsv
#bash ~/bin/Roulette/adding_mutation_rate/roulette_annotate.sh ~/scratch/YBX3_test.tsv ~/ps-yeolab5/roulette
rule reannotate_vcf:
# contig 1 -> chr1
    input:
        vcf="{anything}.vcf.gz",
        rename="/home/hsher/projects/oligoCLIP/utils/rename_chr.txt"
    output:
        "{anything}.rename.vcf.gz",
    threads: 2
    params:
        error_file = "stderr/rename_chr",
        out_file = "stdout/rename_chr",
        run_time = "3:20:00",
        cores = 1,
        memory = "16000",
    container:
        "docker://miguelpmachado/bcftools:1.9-01"
    resources:
        mem_mb=16000
    shell:
        """
        bcftools annotate --rename-chrs {input.rename} \
            {input.vcf} -Oz -o {output}
        bcftools index {output}
        """

rule mutation_rate_in_binding_site:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        vcf = lambda wildcards: Path(ROULETTE_PATH)/(wildcards.chr.replace('chr', '')+'_rate_v5.2_TFBS_correction_all.vcf.bgz')
    output:
        renamed_bed = temp("output/finemapping/mapped_sites/{experiment_label}.{chr}.finemapped_windows.rename.bed.gz"),
        final_vcf = "output/variants/roulette/{experiment_label}.{chr}.vcf"
    threads: 2
    params:
        error_file = "stderr/fetch_roulette_snp.{experiment_label}.{chr}",
        out_file = "stdout/fetch_roulette_snp.{experiment_label}.{chr}",
        run_time = "3:20:00",
        cores = 1,
        memory = 60000,
    container:
        "docker://brianyee/bcftools:1.17"
    shell:
        """
        zcat {input.finemapped_windows} |
        awk '{{ gsub(/chr/,"", $1); print }}' OFS=\"\t\"> {output.renamed_bed}

        bcftools query -R {output.renamed_bed} \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' \
            {input.vcf} \
            | awk '$1="chr"$1' OFS=\"\t\" > {output.final_vcf}
        """



rule fetch_Clinvar_SNP:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        vcf = CLINVAR_VCF.replace('.vcf.gz', '.rename.vcf.gz')
    output:
        "output/variants/clinvar/{experiment_label}.vcf"
    threads: 2
    params:
        error_file = "stderr/fetch_clinvar_snp.{experiment_label}",
        out_file = "stdout/fetch_clinvar_snp.{experiment_label}",
        run_time = "3:20:00",
        cores = 1,
        memory = "16000",
    container:
        "docker://miguelpmachado/bcftools:1.9-01"
    resources:
        mem_mb=16000
    shell:
        """
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
    threads: 2
    params:
        error_file = "stderr/fetch_cosmic_snp.{experiment_label}",
        out_file = "stdout/fetch_cosmic_snp.{experiment_label}",
        run_time = "6:20:00",
        cores = 1,
        memory = "16000",
    container:
        "docker://miguelpmachado/bcftools:1.9-01"
    resources:
        mem_mb=16000
    shell:
        """
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
        memory = "16000",
    container:
        "docker://miguelpmachado/bcftools:1.9-01"
    resources:
        mem_mb=16000
    shell:
        """
        bcftools query -R {input.finemapped_windows} \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SAMPLE_COUNT\t%INFO/TIER\n' \
            {input.vcf} > {output}
        """


rule fetch_sequence:
    input:
        subset_vcf="output/variants/{subset}/{experiment_label_thing}.vcf",
        seq_fa = lambda wildcards: f"output/ml/sequence/{wildcards.experiment_label_thing}.foreground.fa" 
            if 'chr' not in wildcards.experiment_label_thing
            else "output/ml/sequence/"+wildcards.experiment_label_thing.split('.')[0]+".foreground.fa",
        finemapped_windows = lambda wildcards: f"output/finemapping/mapped_sites/{wildcards.experiment_label_thing}.finemapped_windows.bed.gz" 
            if 'chr' not in wildcards.experiment_label_thing
            else "output/finemapping/mapped_sites/"+wildcards.experiment_label_thing.split('.')[0]+".finemapped_windows.bed.gz"
    output:
        ref_fa = "output/variants/{subset}/{experiment_label_thing}.ref.fa",
        alt_fa = "output/variants/{subset}/{experiment_label_thing}.alt.fa",
        csv = "output/variants/{subset}/{experiment_label_thing}.csv"
    threads: 2
    params:
        error_file = "stderr/fetch_sequence.{subset}.{experiment_label_thing}",
        out_file = "stdout/fetch_sequence.{subset}.{experiment_label_thing}",
        run_time = "06:20:00",
        cores = 1,
        out_prefix = lambda wildcards, output: output.csv.replace('.csv', ''),
        memory = "16000",
    conda:
        "/home/hsher/projects/oligoCLIP/rules/envs/metadensity.yaml"
    resources:
        mem_mb=16000
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
        
rule score_variants_ref: 
    input:
        fa = "output/variants/{subset}/{experiment_label_thing}.{label}.fa",
        model = lambda wildcards: f"output/ml/gkmsvm/{wildcards.experiment_label_thing}.model.txt"
        if 'chr' not in wildcards.experiment_label_thing
        else "output/ml/gkmsvm/"+wildcards.experiment_label_thing.split('.')[0]+".model.txt"
        ,
    output:
        score = "output/variants/{subset}/{experiment_label_thing}.{label}.score.txt",
    threads: 2
    params:
        error_file = "stderr/score_ref_variants.{subset}.{experiment_label_thing}.{label}",
        out_file = "stdout/score_ref_variants.{subset}.{experiment_label_thing}.{label}",
        run_time = "4:00:00",
        cores = 1,
        memory = "16000",
    container:
        "docker://algaebrown/lsgkm"
    resources:
        mem_mb=16000
    shell:
        """
        /bin/gkmpredict -T 16 {input.fa} {input.model} {output.score}
        """

rule combine_roulette_results:
    input:
        roulette_scores=expand("output/variants/roulette/{experiment_label}.{chr}.{label}.score.txt",
        chr=['chr'+str(i) for i in list(range(1,23))],
        experiment_label = ['{experiment_label}'],
        label = ['alt', 'ref']),
        roulette_csvs=expand("output/variants/roulette/{experiment_label}.{chr}.csv",
        chr=['chr'+str(i) for i in list(range(1,23))],
        experiment_label = ['{experiment_label}']
        ),
        gnomAD_csv = "output/variants/gnomAD/{experiment_label}.csv",
        finemap_annotation = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv",
    output:
        "output/variants/roulette/{experiment_label}.full.csv.gz",
        "output/variants/roulette/{experiment_label}.filtered.csv.gz"
    threads: 2
    params:
        error_file = "stderr/combine_roulette_results.{experiment_label}",
        out_file = "stdout/combine_roulette_results.{experiment_label}",
        run_time = "00:30:00",
        memory = 40000,
    shell:
        """
        python {TOOL_DIR}/annotate_roulette.py . {wildcards.experiment_label}
        """
rule gnomAD_variant_analysis:
    input:
        ref_score = "output/variants/gnomAD/{experiment_label}.ref.score.txt",
        alt_score = "output/variants/gnomAD/{experiment_label}.alt.score.txt",
        vcf = "output/variants/gnomAD/{experiment_label}.vcf",
        csv = "output/variants/gnomAD/{experiment_label}.csv",
        finemap_annotation = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.annotated.tsv",
    output:
        "output/variants/gnomAD_analysis/{experiment_label}.delta_score_distribution.pdf",
        "output/variants/gnomAD_analysis/{experiment_label}.MAF_vs_impact.pdf",
        "output/variants/gnomAD_analysis/{experiment_label}.feature_stat.tsv",
        "output/variants/gnomAD_analysis/{experiment_label}.transcript_stat.tsv",
        'output/variants/gnomAD/{experiment_label}.variants_scores.tsv'
    params:
        error_file = "stderr/gnomAD_analysis.{experiment_label}",
        out_file = "stdout/genomAD_analysis.{experiment_label}",
        run_time = "00:30:00",
        cores = 1,
        memory = 40000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/annotate_gnomAD.py . {wildcards.experiment_label}
        """

rule gnomAD_permutation_test:
    input:
        'output/variants/gnomAD/{experiment_label}.variants_scores.tsv'
    output:
        'output/variants/gnomAD_analysis/{experiment_label}.test_stats.csv',
        'output/variants/gnomAD_analysis/{experiment_label}.null.csv',
        'output/variants/gnomAD_analysis/{experiment_label}.qqplot.pdf',
    params:
        error_file = "stderr/gnomAD_permutation.{experiment_label}",
        out_file = "stdout/genomAD_permutation.{experiment_label}",
        run_time = "04:30:00",
        cores = 8,
        memory = 160000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/permutation_gnomAD.py {input} \
            gene_name \
            output/variants/gnomAD_analysis/{wildcards.experiment_label}
        """

rule clinvar_variant_analysis:
    input:
        rules.gnomAD_permutation_test.output,
        expand("output/variants/clinvar/{experiment_label}.{subtype}.score.txt",
        experiment_label = ['{experiment_label}'], subtype = ['alt', 'ref']),
        'output/variants/gnomAD/{experiment_label}.variants_scores.tsv'
    output:
        "output/variants/clinvar_analysis/{experiment_label}.constrain_disease_count.csv",
        "output/variants/clinvar_analysis/{experiment_label}.constrain_disease_gene_count.csv",
        "output/variants/clinvar_analysis/{experiment_label}.constrain_disease_test",
        'output/variants/clinvar/{experiment_label}.variants_scores.csv'
    params:
        error_file = "stderr/clinvar_analysis.{experiment_label}",
        out_file = "stdout/clinvar_analysis.{experiment_label}",
        run_time = "00:30:00",
        cores = 1,
        memory = 40000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/annotate_clinvar.py . {wildcards.experiment_label}
        """

rule cosmic_variant_analysis:
    input:
        rules.gnomAD_permutation_test.output,
        expand("output/variants/cosmic_coding/{experiment_label}.{subtype}.score.txt",
        experiment_label = ['{experiment_label}'], subtype = ['alt', 'ref']),
        expand("output/variants/cosmic_noncoding/{experiment_label}.{subtype}.score.txt",
        experiment_label = ['{experiment_label}'], subtype = ['alt', 'ref']),
        'output/variants/gnomAD/{experiment_label}.variants_scores.tsv'
    output:
        "output/variants/cosmic_analysis/{experiment_label}.pos_selection_tsg.csv",
        "output/variants/cosmic_analysis/{experiment_label}.pos_selection_oncogene.csv",
        "output/variants/cosmic_analysis/{experiment_label}.impact_vs_role_stat.csv",
        'output/variants/cosmic_analysis/{experiment_label}.impact_vs_tier_stat.csv',
        'output/variants/cosmic_noncoding/{experiment_label}.variants_scores.csv'
    params:
        error_file = "stderr/cosmic_analysis.{experiment_label}",
        out_file = "stdout/cosmic_analysis.{experiment_label}",
        run_time = "00:30:00",
        cores = 1,
        memory = 40000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/annotate_cosmic.py . {wildcards.experiment_label}
        """
def find_well_trained_model(wildcards):
    auprc_threshold = 0.6
    
    auprc_df =  pd.read_csv(checkpoints.ml_gkmsvm_AUPRC.get().output[0], index_col = 0)
    print('RBP with good enough model to do variant interpretation:')
    print(auprc_df.loc[auprc_df['mean AUPRC']>auprc_threshold])
    with_good_model = auprc_df.loc[auprc_df['mean AUPRC']>auprc_threshold, 'Experiment'].tolist()
    

    return expand("output/variants/{variant_set}/{experiment_label}.{type}.score.txt",
        experiment_label = with_good_model, 
        variant_set = ['gnomAD', 'clinvar', 'cosmic_coding', 'cosmic_noncoding'], 
        type = ['alt', 'ref']
        )+expand("output/variants/gnomAD_analysis/{experiment_label}.test_stats.csv",
        experiment_label = with_good_model
        )+expand('output/variants/clinvar/{experiment_label}.variants_scores.csv',
        experiment_label = with_good_model
        )+expand('output/variants/cosmic_noncoding/{experiment_label}.variants_scores.csv',
        experiment_label = with_good_model
        )+expand("output/variants/roulette/{experiment_label}.filtered.csv.gz",
        experiment_label = with_good_model
        )
        

rule variants_done:
    input:
        find_well_trained_model
    output:
        "variants_done.txt"
    params:
        error_file = "stderr/variants_done",
        out_file = "stdout/variants_done",
        run_time = "00:20:00",
        cores = 1,
        memory = "1000",
    resources:
        mem_mb=1000
    shell:
        """
        touch {output}
        """
