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

# TODO: optimize thread use
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
        run_time = "24:10:00",
        memory = "320000",
        job_name = "run_homer",
        prefix = lambda wildcards, output: output.model.replace('.model.txt', '')
    container:
        "docker://algaebrown/lsgkm"
    benchmark: "benchmarks/gkmsvm/train.{experiment_label}.txt"
    shell:
        """
        /usr/src/lsgkm/bin/gkmtrain \
            -R \
            -m 280000 \
            {input.foreground} \
            {input.background} \
            {params.prefix}  
        """
# TODO: optimize memory usage and timeoue.
# on TSCC 1.0 -m 2000 used to work
# but after containerize -m 1000 goes segmentation fault
# -T causes segmentation fault anyways
# -T 1, -m 1000 TARDBP uses 4 hr to CV on TSCC 1.0
# -T can only be 1, 4 or 16
# 8 CPU, 16 thread, 32GB mem, -m 31GB died
# 8 CPU, 16 thread, 64GB mem, -m 1.8BG died
# Please note that it will run SVM training N times, which can take time if training sets are large. In this case, you can perform CV analysis on a specific set by using -i <I> option for parallel runnings. The output will be <outprefix>.cvpred.<I>.txt
# tried -c in slurm, thread per CPU, doesn't work either
rule cv_gkmsvm:
    input:
        foreground = rules.fetch_sequence.output.finemapped_fa,
        background = rules.fetch_sequence.output.background_fa,
    output:
        cv = "output/ml/gkmsvm/{experiment_label}.cvpred.txt",
    params:
        error_file = "stderr/{experiment_label}.cv.err",
        out_file = "stdout/{experiment_label}.cv.out",
        run_time = "24:00:00",
        memory = "640000",
        job_name = "gkmsvm_cv",
        prefix = lambda wildcards, output: output.cv.replace('.cvpred.txt', '')
    threads: 1
    container:
        "docker://algaebrown/lsgkm"
    benchmark: "benchmarks/gkmsvm/cv.{experiment_label}.txt"
    shell:
        """
        if ! test -s {input.foreground}
            then
            echo "No reproducible finemapped sequences" > {output.cv}
        else
            timeout -t 86000 /usr/src/lsgkm/bin/gkmtrain \
                -x 5 \
                -R \
                -T 1 \
                -m 620000 \
                {input.foreground} \
                {input.background} \
                {params.prefix} 
        fi
        
        if ! test -f {output.cv} ; then
           echo "Timeout. gkmsvm does not converge" > {output.cv}
        fi
        """
        # https://stackoverflow.com/questions/50382603/how-to-timeout-with-exit0-from-bash

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


