if not GENOME:
  rule map_filter_sort:
    input:
        fastq=OUT + "/trimmed/{sample}.fastq.gz",
        index=INDEX_FILES, reference_check=REFERENCE_CHECK
    output:
        bam=OUT + "/bam/{sample}.accepted.bam",
        diagnostic=([OUT + "/diagnostics/{sample}.unfiltered.sam"] if SETTINGS["diagnostic"] else [])
    params:
        prefix=INDEX, orientation=("--norc" if SETTINGS["norc"] else ""),
        preset=SETTINGS["preset"], seed=SETTINGS["seed"],
        mapq=SETTINGS["mapq"], flags=SETTINGS["flags"],
        bowtie_threads=lambda wc, threads: max(1, threads - (3 if SETTINGS["diagnostic"] else 2)),
        diagnostic_pipe=lambda wc: ("tee " + shlex.quote(OUT + "/diagnostics/" + wc.sample + ".unfiltered.sam") + " | " if SETTINGS["diagnostic"] else "")
    threads: SETTINGS["threads"]
    resources: mem_mb=4000, runtime=120
    log:
        align=OUT + "/logs/{sample}.bowtie2.log",
        filter=OUT + "/logs/{sample}.filter.log",
        sort=OUT + "/logs/{sample}.sort.log"
    conda: "../envs/mapping.yaml"
    shell:
        "bowtie2 --version > {log.align:q}; samtools --version > {log.filter:q}; "
        "bowtie2 -p {params.bowtie_threads} --{params.preset} {params.orientation} "
        "--seed {params.seed} -x {params.prefix:q} -U {input.fastq:q} 2>> {log.align:q} "
        "| {params.diagnostic_pipe}samtools view -u -F {params.flags} -q {params.mapq} - 2>> {log.filter:q} "
        "| samtools sort -@ 0 -m 256M -T {output.bam:q}.sorttmp "
        "-o {output.bam:q} - 2> {log.sort:q}"


else:
  rule map_rna_handoff:
    input:
        fastq=OUT+'/trimmed/{sample}.fastq.gz', reference=REF,
        index=INDEX_FILES, verified=REFERENCE_CHECK, code=CODE
    output:
        bam=OUT+'/rna/{sample}/accepted.bam',
        all_bam=OUT+'/rna/{sample}/alignments.bam',
        unmapped=OUT+'/rna/{sample}/unmapped.fastq',
        qc=OUT+'/rna/{sample}/qc.json',
        provenance=OUT+'/rna/{sample}/provenance.json',
        publication=OUT+'/rna/{sample}/publication.json',
        commands=OUT+'/rna/{sample}/commands.jsonl',
        sizes=OUT+'/rna/{sample}/reference.sizes.tsv'
    params: step='rna', settings=SETTINGS, repo=str(REPO)
    threads: SETTINGS['threads']
    resources: mem_mb=4000, runtime=120
    log: OUT+'/logs/{sample}.rna.jsonl'
    conda: '../envs/sequential.yaml'
    script: '../scripts/genome_step.py'


rule index_bam:
    input: RNA_BAM
    output: RNA_BAM + '.bai'
    threads: 1
    resources: mem_mb=512, runtime=10
    log: OUT + "/logs/{sample}.index.log"
    conda: "../envs/samtools.yaml"
    shell: "samtools index {input:q} {output:q} 2> {log:q}"
