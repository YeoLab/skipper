workdir: config['WORKDIR']
import pandas as pd
"""
snakemake -kps rules/visualize_ml.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2_single -n

# making splice site predictions
snakemake -kps rules/visualize_ml.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_HepG2_20230620.yaml \
    --profile profiles/tscc2_single output/gkmsvm_explain/windows/{AQR_HepG2_ENCSR018WPY,BUD13_HepG2_ENCSR830BSQ,U2AF2_HepG2_ENCSR202BFN,SUGP2_HepG2_ENCSR506UPY,HNRNPC_HepG2_ENCSR550DVK,PTBP1_HepG2_ENCSR384KAN,HLTF_HepG2_ENCSR647HOX,GTF2F1_HepG2_ENCSR265ZIS,BCLAF1_HepG2_ENCSR876EYA,SF3A3_HepG2_ENCSR331MIC,QKI_HepG2_ENCSR570WLM,MATR3_HepG2_ENCSR290VLT,SF3B4_HepG2_ENCSR279UJF,ILF3_HepG2_ENCSR786TSC,HNRNPK_HepG2_ENCSR828ZID,UCHL5_HepG2_ENCSR490IEE,RBFOX2_HepG2_ENCSR987FTF,TAF15_HepG2_ENCSR841EQA,CSTF2_HepG2_ENCSR384MWO,DDX59_HepG2_ENCSR214BZA,FAM120A_HepG2_ENCSR987NYS,PRPF8_HepG2_ENCSR121NVA,HNRNPM_HepG2_ENCSR267UCX,KHSRP_HepG2_ENCSR366DGX,NKRF_HepG2_ENCSR277DEO,PCBP2_HepG2_ENCSR339FUY,AGGF1_HepG2_ENCSR543TPH,CSTF2T_HepG2_ENCSR919HSE,GRWD1_HepG2_ENCSR893NWB,EFTUD2_HepG2_ENCSR527DXF,HNRNPL_HepG2_ENCSR724RDN,SFPQ_HepG2_ENCSR965DLL,PPIG_HepG2_ENCSR097NEE,XPO5_HepG2_ENCSR921SXC,PRPF4_HepG2_ENCSR977OXG}.709.705.gkmexplain.txt \
    -n
# compare with dlogodds
snakemake -kps rules/visualize_ml.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_HepG2_20230620.yaml \
    --profile profiles/tscc2_single output/ml/sequence/reproducible_enriched_windows/{RBFOX2_HepG2_ENCSR987FTF,AQR_HepG2_ENCSR018WPY,BUD13_HepG2_ENCSR830BSQ,U2AF2_HepG2_ENCSR202BFN,PRPF4_HepG2_ENCSR977OXG}.txt \
    -n
# compare with dlogodds
snakemake -kps rules/visualize_ml.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_HepG2_20230620.yaml \
    --profile profiles/tscc2 \
    -n

# compare with biophysical
snakemake -kps rules/visualize_ml.smk \
    --configfile /tscc/nfs/home/hsher/projects/skipper/encode_configs/encode_pe_rules_config_K562_20230929.yaml \
    --profile profiles/tscc2_single output/ml/PUM_biophysical_model_enriched/{PUM1_K562_ENCSR308YNT,PUM2_K562_ENCSR661ICQ}/individual_scores.0.tsv \
    -n
"""
locals().update(config)
auprc = pd.read_csv('output/ml/gkmsvm/AUPRC.txt', index_col = 0)
pass_exp = auprc.loc[auprc['mean AUPRC']>0.65, 'Experiment'].tolist()


rule all:
    input:
        lambda wildcards: [f"output/gkmsvm_explain/{experiment_label}.foreground.{experiment_label}.gkmexplain.txt"
        for experiment_label in pass_exp
        ],
        expand("output/ml/sequence/reproducible_enriched_windows/{experiment_label}.txt", 
        experiment_label = pass_exp)


rule fix_model:
    input: 
        model = "output/ml/gkmsvm/{experiment_label}.model.txt"
    output: 
        "output/ml/gkmsvm/{experiment_label}.model.fix.txt"
    params:
        error_file = "stderr/{experiment_label}.fix_model.err",
        out_file = "stdout/{experiment_label}.fix_model.out",
        run_time = "10:00",
        memory = "1000",
        job_name = "run_homer",
    shell: 
        """
        sed -e 's/norc/gamma/g' {input.model} > {output}
        """


rule gkmsvm_explain:
    input: 
        fa="output/ml/sequence/{something}.fa",
        model = "output/ml/gkmsvm/{experiment_label}.model.fix.txt"
    output:
        output="output/gkmsvm_explain/{something}.{experiment_label}.gkmexplain.txt"
    container: 
        "docker://kundajelab/lsgkm:latest"
    params:
        error_file = "stderr/{experiment_label}.{something}.gkmexplain.err",
        out_file = "stdout/{experiment_label}.{something}.gkmexplin.out",
        run_time = "10:00",
        memory = "40000",
        job_name = "gkmsvm_explain",
    shell: 
        """
        head -n 40 {input.fa} > {input.fa}.head
        /opt/lsgkm/src/gkmexplain \
                {input.fa}.head \
                {input.model} \
                {output.output}
        """

rule window_to_sequence_with_range:
    input:
        config['PARTITION']
    output:
        "output/gkmsvm_explain/windows/sequence/{max}.{min}.fa"
    params:
        error_file = "stderr/{max}.{min}.fa.err",
        out_file = "stdout/{max}.{min}.fa.out",
        run_time = "10:00",
        memory = "1000",
        job_name = "window_to_sequence",
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        zcat {input} | awk '{{ if ($4 >= {wildcards.min} && $4 <= {wildcards.max}) print}}' |\
        bedtools getfasta -fi {GENOME} -bed stdin -s -fo {output}
        """

rule gkmsvm_explain_windows:
    input:
        fa="output/gkmsvm_explain/windows/sequence/{max}.{min}.fa",
        model = "output/ml/gkmsvm/{experiment_label}.model.fix.txt"
    output:
        output="output/gkmsvm_explain/windows/{experiment_label}.{max}.{min}.gkmexplain.txt"
    params:
        error_file = "stderr/{experiment_label}.{max}.{min}.gkmexplain.err",
        out_file = "stdout/{experiment_label}.{max}.{min}.gkmexplin.out",
        run_time = "10:00",
        memory = "40000",
        job_name = "gkmsvm_explain",
    container: 
        "docker://kundajelab/lsgkm:latest"
    shell: 
        """
        /opt/lsgkm/src/gkmexplain \
                {input.fa} \
                {input.model} \
                {output.output}
        """

rule sequence_from_reproducible_enriched_windows:
    input:
        "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz"
    output:
        temp("output/ml/sequence/reproducible_enriched_windows/{experiment_label}.fa")
    params:
        error_file = "stderr/{experiment_label}.fa.err",
        out_file = "stdout/{experiment_label}.fa.out",
        run_time = "10:00",
        memory = "1000",
        job_name = "sequence_from_reproducible_enriched_windows"
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """
        zcat {input} | tail -n +2 | bedtools getfasta -fi {GENOME} -bed stdin -s -fo {output}
        """

rule predict_enriched_windows:
    input:
        fa = "output/ml/sequence/reproducible_enriched_windows/{experiment_label}.fa",
        model = "output/ml/gkmsvm/{experiment_label}.model.txt",
    output:
        score = "output/ml/sequence/reproducible_enriched_windows/{experiment_label}.txt",
    threads: 2
    params:
        error_file = "stderr/score_ref_variants.{experiment_label}",
        out_file = "stdout/score_ref_variants.{experiment_label}",
        run_time = "4:00:00",
        cores = 1,
        memory = 40000, #320000 is ceiling of gold slow
    container:
        "docker://shl198/lsgkm:0.1.1"
    shell:
        """
        /bin/gkmpredict -T 16 {input.fa} {input.model} {output.score}
        """

rule cleaup_fa:
    input:
        "{anything}.fa"
    output:
        replacet=temp("{anything}.U.fa"),
        removen=temp("{anything}.UN.fa")
    threads: 1
    params:
        error_file = "stderr/score_ref_variants.{anything}",
        out_file = "stdout/score_ref_variants.{anything}",
        run_time = "00:00:05",
        cores = 1,
        memory = 40000, #320000 is ceiling of gold slow
    conda: "envs/metadensity.yaml"
    shell:
        """
        sed '/^[^>]/s/T/U/g' {input} > {output.replacet}
        python /tscc/nfs/home/hsher/projects/PUM2_biophysical/removeN.py {output.replacet} {output.removen}
        """

BIOPHYSICAL_MODEL_PATH='/tscc/nfs/home/hsher/projects/PUM2_biophysical/final_pum2_w_flips_coupling_eab_predict.py'
rule pum2_biophysical_model:
    input:
        foreground="output/ml/sequence/{experiment_label}.foreground.UN.fa",
        background="output/ml/sequence/{experiment_label}.background.UN.fa",
    output:
        "output/ml/PUM_biophysical_model/{experiment_label}/individual_scores.0.tsv",
    threads: 1
    params:
        error_file = "stderr/pum_biophysical.{experiment_label}",
        out_file = "stdout/pum_biophysical.{experiment_label}",
        run_time = "00:20:00",
        cores = 1,
        memory = 40000, #320000 is ceiling of gold slow
    conda: "envs/pum2model.yaml"
    shell:
        """
        python {BIOPHYSICAL_MODEL_PATH} \
            --pos_fasta {input.foreground} \
            --neg_fasta {input.background} \
            --out_dir output/ml/PUM_biophysical_model/{wildcards.experiment_label}
        """

rule pum2_biophysical_model_enriched_windows:
    input:
        foreground="output/ml/sequence/reproducible_enriched_windows/{experiment_label}.UN.fa",
        background="output/ml/sequence/{experiment_label}.background.UN.fa",
    output:
        "output/ml/PUM_biophysical_model_enriched/{experiment_label}/individual_scores.0.tsv",
    threads: 1
    params:
        error_file = "stderr/pum_biophysical.{experiment_label}",
        out_file = "stdout/pum_biophysical.{experiment_label}",
        run_time = "00:20:00",
        cores = 1,
        memory = 40000, #320000 is ceiling of gold slow
    conda: "envs/pum2model.yaml"
    shell:
        """
        python {BIOPHYSICAL_MODEL_PATH} \
            --pos_fasta {input.foreground} \
            --neg_fasta {input.background} \
            --out_dir output/ml/PUM_biophysical_model_enriched/{wildcards.experiment_label}
        """