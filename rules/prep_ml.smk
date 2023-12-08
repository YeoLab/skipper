import pandas as pd
locals().update(config)
rule fetch_sequence:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        background = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz",
    output:
        finemapped_fa = "output/ml/sequence/{experiment_label}.foreground.fa",
        background_fa = "output/ml/sequence/{experiment_label}.background.fa"
    params:
        error_file = "stderr/{experiment_label}.fetch_sequence.err",
        out_file = "stdout/{experiment_label}.fetch_sequence.out",
        run_time = "40:00",
        memory = "2000",
        job_name = "run_homer",
        fa = config['GENOME']
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    shell:
        '''
        bedtools getfasta -fo {output.finemapped_fa} -fi {params.fa} -bed {input.finemapped_windows} -s
        bedtools getfasta -fo {output.background_fa} -fi {params.fa} -bed {input.background} -s
        '''

rule train_gkmsvm:
    # jobs exceeding time limit are those that cannot converge
    input:
        foreground = rules.fetch_sequence.output.finemapped_fa,
        background = rules.fetch_sequence.output.background_fa,
    output:
        model = "output/ml/gkmsvm/{experiment_label}.model.txt",
    params:
        error_file = "stderr/{experiment_label}.train_gkmsvm.err",
        out_file = "stdout/{experiment_label}.train_gkmsvm.out",
        run_time = "3:10:00",
        memory = "2000",
        job_name = "run_homer",
        prefix = lambda wildcards, output: output.model.replace('.model.txt', '')
    container:
        "docker://kundajelab/lsgkm"
    shell:
        """
        timeout 3h gkmtrain \
            -R \
            -T 4 \
            -m 2000 \
            {input.foreground} \
            {input.background} \
            {params.prefix}  || true
        
        if ! test -f {output.model} ; then
            echo "Timeout. gkmsvm does not converge" > {output.model}
        fi
        """
# TODO: need to containerize lsgkm latest version with -R option. The ones on conda does not have it
# lsgkm does not have maximal iterations, leading to infinite training/cv for some dataset such as RBFOX3
rule cv_gkmsvm:
    input:
        foreground = rules.fetch_sequence.output.finemapped_fa,
        background = rules.fetch_sequence.output.background_fa,
    output:
        cv = "output/ml/gkmsvm/{experiment_label}.cvpred.txt",
    params:
        error_file = "stderr/{experiment_label}.fetch_sequence.err",
        out_file = "stdout/{experiment_label}.fetch_sequence.out",
        run_time = "6:10:00",
        memory = "2000",
        job_name = "gkmsvm_cv",
        prefix = lambda wildcards, output: output.cv.replace('.cvpred.txt', '')
    container:
        "docker://kundajelab/lsgkm"
    shell:
        """
        timeout 6h gkmtrain \
            -x 5 \
            -T 4 \
            -R \
            -m 2000 \
            {input.foreground} \
            {input.background} \
            {params.prefix} || true

        if ! test -f {output.cv} ; then
            echo "Timeout. gkmsvm does not converge" > {output.cv}
        fi
        """
checkpoint gkmsvm_AUPRC:
    input:
        cv_output = expand(rules.cv_gkmsvm.output.cv, experiment_label = config['manifest'].Experiment)
    output:
        "output/ml/gkmsvm/AUPRC.txt"
    params:
        error_file = "stderr/gkmsvm_AUPRC.err",
        out_file = "stdout/gkmsvm_AUPRC.out",
        run_time = "40:00",
        memory = "2000",
        job_name = "gkmsvm_cv",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {TOOL_DIR}/gather_gkmsvm_auprc.py {output}
        """


# rule make_prismnet_tsv:


