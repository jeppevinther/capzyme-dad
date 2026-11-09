rule collapse_ncrna:
    input:
        fasta=SETTINGS["ncrna"], legacy=COLLAPSER, core=CORE
    output:
        fasta=OUT + "/reference/ncrna.fa",
        membership=OUT + "/reference/ncrna.tsv"
    params:
        step="collapse", prefix=OUT + "/reference/ncrna",
        k=SETTINGS["k"], drops=SETTINGS["drops"]
    threads: 1
    resources: mem_mb=1024, runtime=30
    log: OUT + "/logs/collapse_ncrna.json"
    conda: "../envs/python.yaml"
    script: "../scripts/run_step.py"


rule combine_reference:
    input:
        # Preserve original experiment ordering: explicitly supplied viral records first.
        fastas=SETTINGS["viral"] + [OUT + "/reference/ncrna.fa"], core=CORE
    output: REF
    params: step="combine"
    threads: 1
    resources: mem_mb=1024, runtime=10
    log: OUT + "/logs/combine_reference.json"
    conda: "../envs/python.yaml"
    script: "../scripts/run_step.py"


rule index_reference:
    input: REF
    output: REF + ".fai"
    threads: 1
    resources: mem_mb=512, runtime=10
    log: OUT + "/logs/faidx.log"
    conda: "../envs/samtools.yaml"
    shell:
        "samtools --version > {log:q}; samtools faidx {input:q} 2>> {log:q}"


rule bowtie2_index:
    input: REF
    output:
        multiext(INDEX, ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
    threads: SETTINGS["threads"]
    resources: mem_mb=4000, runtime=60
    log: OUT + "/logs/bowtie2_build.log"
    conda: "../envs/bowtie2.yaml"
    shell:
        # A single deterministic suffix family also supports references over 4 Gb.
        "bowtie2-build --version > {log:q}; "
        "bowtie2-build --large-index --threads {threads} {input:q} {INDEX:q} >> {log:q} 2>&1"
