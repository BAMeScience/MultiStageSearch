rule getTargets:
    input:
        PSM_Report = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
        taxIDs = RESULT_DIR / "taxidMapping/protacc2taxids_virus.txt"
    output:
        RESULT_DIR / "{sample}/taxids/mapped_taxids.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/taxids/getTargets/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/taxids/getTargets/{sample}/stdout.log"
    params:
        query= RESULT_DIR / "{sample}/ref_query_accessions.txt",
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/taxIDMapping.py"