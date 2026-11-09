rule trim_fastq:
    input:
        fastq=lambda wc: SETTINGS["samples"][wc.sample],
        reference_check=REFERENCE_CHECK
    output:
        fastq=OUT + "/trimmed/{sample}.fastq.gz",
        json=OUT + "/qc/{sample}.cutadapt.json"
    params:
        adapter=SETTINGS["pre"]["adapter"], quality=SETTINGS["pre"]["quality"],
        minimum=SETTINGS["pre"]["minimum_length"], nextseq=("" if SETTINGS["pre"]["nextseq_trim"] is None else "--nextseq-trim="+str(SETTINGS["pre"]["nextseq_trim"]))
    threads: 1
    resources: mem_mb=1000, runtime=60
    log: OUT + "/logs/{sample}.cutadapt.log"
    conda: "../envs/cutadapt.yaml"
    shell:
        "cutadapt --version > {log:q}; "
        "cutadapt -j {threads} -a {params.adapter:q} -q {params.quality} "
        "-m {params.minimum} {params.nextseq} "
        "--json {output.json:q} -o {output.fastq:q} {input.fastq:q} >> {log:q} 2>&1"
