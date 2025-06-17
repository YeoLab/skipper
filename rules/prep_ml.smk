import pandas as pd
locals().update(config)
rule fetch_sequence:
    input:
        finemapped_windows = "output/finemapping/mapped_sites/{experiment_label}.finemapped_windows.bed.gz",
        background = "output/homer/region_matched_background/fixed/{experiment_label}.sampled_fixed_windows.bed.gz",
    output:
        finemapped_fa = "output/ml/sequence/{experiment_label}.foreground.fa",
        background_fa = "output/ml/sequence/{experiment_label}.background.fa"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    resources:
        mem_mb=2000,
        runtime=40
    shell:
        '''
        bedtools getfasta -fo {output.finemapped_fa} -fi {GENOME} -bed {input.finemapped_windows} -s
        bedtools getfasta -fo {output.background_fa} -fi {GENOME} -bed {input.background} -s
        '''
