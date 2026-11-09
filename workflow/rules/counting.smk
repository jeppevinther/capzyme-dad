rule count_5ends:
    input:
        bams=expand(RNA_BAM, sample=SAMPLES),
        indexes=expand(RNA_BAM + '.bai', sample=SAMPLES),
        fai=SIZES, core=CORE, reference_check=REFERENCE_CHECK
    output: RAW
    params: step="matrix", names=SAMPLES, strand=SETTINGS["strand"]
    threads: 2 if GENOME else 1
    resources: mem_mb=4000, runtime=60
    log: OUT + "/logs/count_5ends.json"
    conda: "../envs/counting.yaml"
    script: "../scripts/run_step.py"


rule threshold_counts:
    input: matrix=RAW, core=CORE
    output: FILTERED
    threads: 1
    resources: mem_mb=256, runtime=10
    params: step="threshold", minimum=SETTINGS["minimum"]
    log: OUT + "/logs/threshold_counts.json"
    conda: "../envs/python.yaml"
    script: "../scripts/run_step.py"


if not GENOME:
  rule record_provenance:
    input:
        matrix=RAW,
        threshold=([FILTERED] if SETTINGS["minimum"] > 0 else []),
        manifest=SETTINGS["manifest"], code=CODE,
        tool_logs=[OUT + "/logs/" + SAMPLES[0] + ".cutadapt.log",
                   OUT + "/logs/" + SAMPLES[0] + ".bowtie2.log",
                   OUT + "/logs/" + SAMPLES[0] + ".filter.log",
                   OUT + "/logs/count_5ends.json"]
    output: OUT + "/provenance/run.json"
    params:
        step="provenance", settings=SETTINGS, repo=str(REPO),
        sources=[SETTINGS["manifest"], *SETTINGS["reference_sources"], *SETTINGS["samples"].values()]
    log: OUT + "/logs/provenance.json"
    threads: 1
    resources: mem_mb=1024, runtime=10
    # Coordinator environment captures actual host-tool versions; no separate Conda env.
    script: "../scripts/run_step.py"
else:
  rule record_sequential_provenance:
    input:
        matrix=RAW, genome_matrix=OUT + '/counts/genome.5ends.raw.tsv',
        threshold=([FILTERED] if SETTINGS['minimum'] > 0 else []),
        manifest=SETTINGS['manifest'], code=CODE,
        raw=list(SETTINGS['samples'].values()),
        trim=expand(OUT+'/qc/{sample}.cutadapt.json', sample=SAMPLES),
        accounting=expand(OUT+'/qc/{sample}.sequential.json', sample=SAMPLES),
        stages=expand(OUT+'/rna/{sample}/provenance.json', sample=SAMPLES),
        publications=expand(OUT+'/rna/{sample}/publication.json', sample=SAMPLES),
        genome=OUT+'/provenance/genome.json'
    output: OUT+'/provenance/run.json'
    params: step='provenance', settings=SETTINGS, repo=str(REPO), code=CODE
    threads: 1
    resources: mem_mb=1024, runtime=30
    log: OUT+'/logs/provenance.json'
    script: '../scripts/genome_step.py'
