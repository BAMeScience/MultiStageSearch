rule filterUnidentified:
    input: 
        mgf = Filtering.get_input_filter_undefined()["mgf"],
        peptide_shaker_report = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
    output: 
        filtered_mgf = RESULT_DIR / "{sample}/RefFilter/Filtered_ref.mgf"
    log:
        stdout_log = RESULT_DIR / "logs/RefFilter/{sample}/stdout.log"
    conda: 
        '../envs/pandas.yml'
    threads: 1
    script: 
        "../scripts/ref_filtering.py"