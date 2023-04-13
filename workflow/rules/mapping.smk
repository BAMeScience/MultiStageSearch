rule mapTaxIDs:
    input:
        PSM_Report = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
        taxIDs = RESULT_DIR / "taxidMapping/protacc2taxids_virus.txt"
    output:
        all_mapped_taxIDs = RESULT_DIR / "{sample}/taxids/mapped_taxids.txt",
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/taxids/mapTaxIDs/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/taxids/mapTaxIDs/{sample}/stdout.log"
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/taxIDMapping.py"