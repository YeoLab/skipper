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
    resources:
        mem_mb=2000
    shell:
        '''
        bedtools getfasta -fo {output.finemapped_fa} -fi {params.fa} -bed {input.finemapped_windows} -s
        bedtools getfasta -fo {output.background_fa} -fi {params.fa} -bed {input.background} -s
        '''