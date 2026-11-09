# No reference/index generation rules are registered in fixed mode.
rule verify_fixed_reference:
    input:
        fasta=REF, index=INDEX_FILES, core=CORE,
        verifier=str(REPO / "workflow/scripts/fixed_reference.py")
    output:
        report=OUT + "/reference/fixed_verified.json",
        sizes=SIZES
    threads: 1
    resources: mem_mb=512, runtime=10
    params: step="verify_fixed", reference=SETTINGS["fixed"]
    log: OUT + "/logs/fixed_reference.json"
    conda: "../envs/fixed.yaml"
    script: "../scripts/run_step.py"
