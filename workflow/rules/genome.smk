rule map_genome:
    input:
        fastq=(lambda wc: SETTINGS['samples'][wc.sample]) if GENOME_ONLY else OUT+'/rna/{sample}/unmapped.fastq',
        reference=GENOME['reference'], index=GENOME['index_files'],
        decoys=([GENOME['decoys']] if GENOME['decoys'] else []),
        policies=GENOME['policy_files'], code=CODE
    output:
        bam=OUT+'/genome/{sample}/accepted.bam',
        all_bam=OUT+'/genome/{sample}/alignments.bam',
        unmapped=OUT+'/genome/{sample}/unmapped.fastq',
        qc=OUT+'/genome/{sample}/qc.json',
        provenance=OUT+'/genome/{sample}/provenance.json',
        publication=OUT+'/genome/{sample}/publication.json',
        commands=OUT+'/genome/{sample}/commands.jsonl',
        sizes=OUT+'/genome/{sample}/reference.sizes.tsv',
        mapping=OUT+'/genome/{sample}/mapping_by_reference.tsv'
    params: step='genome', settings=SETTINGS, repo=str(REPO)
    threads: GENOME['threads']
    resources: mem_mb=GENOME['mem_mb'], runtime=240
    log: OUT+'/logs/{sample}.genome.jsonl'
    conda: '../envs/sequential.yaml'
    script: '../scripts/genome_step.py'


rule count_genome:
    input:
        bams=expand(OUT+'/genome/{sample}/accepted.bam', sample=SAMPLES),
        sizes=OUT+'/genome/'+SAMPLES[0]+'/reference.sizes.tsv', code=CODE
    output: OUT+'/counts/genome.5ends.raw.tsv'
    params: step='genome_count', settings=SETTINGS, repo=str(REPO)
    threads: GENOME['count_threads']
    resources: mem_mb=GENOME['count_mem_mb'], runtime=120
    log: OUT+'/logs/genome.count.jsonl'
    conda: '../envs/counting.yaml'
    script: '../scripts/genome_step.py'


rule genome_provenance:
    input:
        matrix=OUT+'/counts/genome.5ends.raw.tsv', manifest=SETTINGS['manifest'],
        sources=GENOME['sources'], code=CODE,
        stages=expand(OUT+'/genome/{sample}/provenance.json', sample=SAMPLES),
        publications=expand(OUT+'/genome/{sample}/publication.json', sample=SAMPLES)
    output: OUT+'/provenance/genome.json'
    params: step='provenance', settings=SETTINGS, repo=str(REPO), code=CODE
    threads: 1
    resources: mem_mb=1024, runtime=30
    log: OUT+'/logs/genome.provenance.jsonl'
    script: '../scripts/genome_step.py'


if not GENOME_ONLY:
  rule sequential_accounting:
    input:
        trim=OUT+'/qc/{sample}.cutadapt.json',
        rna=OUT+'/rna/{sample}/alignments.bam',
        accepted_rna=OUT+'/rna/{sample}/accepted.bam',
        genome=OUT+'/genome/{sample}/alignments.bam',
        accepted_genome=OUT+'/genome/{sample}/accepted.bam',
        unmapped=OUT+'/rna/{sample}/unmapped.fastq',
        qc=[OUT+'/rna/{sample}/qc.json', OUT+'/genome/{sample}/qc.json'],
        provenance=[OUT+'/rna/{sample}/provenance.json', OUT+'/genome/{sample}/provenance.json'],
        publications=[OUT+'/rna/{sample}/publication.json', OUT+'/genome/{sample}/publication.json'],
        code=CODE
    output: OUT+'/qc/{sample}.sequential.json'
    params:
        step='accounting', settings=SETTINGS, repo=str(REPO),
        rna=lambda wc: OUT+'/rna/'+wc.sample,
        genome=lambda wc: OUT+'/genome/'+wc.sample
    threads: 2
    resources: mem_mb=2048, runtime=60
    log: OUT+'/logs/{sample}.accounting.jsonl'
    conda: '../envs/sequential.yaml'
    script: '../scripts/genome_step.py'
