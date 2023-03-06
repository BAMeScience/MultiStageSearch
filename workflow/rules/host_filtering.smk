# from PepGM, edited

rule AddDecoysHost:
    input: 
        host_fasta = HOST_FASTA
    output: 
        touch(RESULT_DIR / "{sample}/Database/AddDecoysHost.done")
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/HostDB/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/HostDB/{sample}/stdout.log" 
    conda:
        "../envs/host_filtering.yml"
    threads: 1
    shell: 
        "searchgui eu.isas.searchgui.cmd.FastaCLI -in {input.host_fasta} -decoy > {log.stdout_log} 2> {log.stderr_log}"

rule MoveHostDecoyDB:
    input:
        RESULT_DIR / "{sample}/Database/AddDecoysHost.done",
    output:
        moved_host_decoy_fasta = RESULT_DIR / "{sample}/Database/{sample}_protein_concatenated_target_decoy.fasta"
    params:
        res_dir = RESULT_DIR,
        sample_name = "{sample}"
    threads: 1
    script:
        "../scripts/rename_db.py"



rule SearchHostSpectra:
    input: 
        mgf = MGF_FILE,
        host_decoy_db = RESULT_DIR / "{sample}/Database/{sample}_protein_concatenated_target_decoy.fasta",
        par = PAR_FILE
    output:  
        out_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip"
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/SearchHostSpectra/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/SearchHostSpectra/{sample}/stdout.log"
    params:
        result_dir = str(RESULT_DIR / "{sample}/SpectraFilter"),
        hostname = "host",
        search_engine = config["host_filtering"]["search_engine"],
        psm_fdr = config["host_filtering"]["psm_fdr"],
        peptide_fdr = config["host_filtering"]["peptide_fdr"],
        protein_fdr = config["host_filtering"]["protein_fdr"]
    conda:
        "../envs/host_filtering.yml"
    threads: workflow.cores / 2
    shell: 
        "searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.host_decoy_db} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.hostname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule RunPeptideShakerHost:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip",
        mgf = MGF_FILE,
        host_decoy_db = RESULT_DIR / "{sample}/Database/{sample}_protein_concatenated_target_decoy.fasta",
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/SpectraFilter/host.psdb"
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/RunPeptideShakerHost/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/RunPeptideShakerHost/{sample}/stdout.log",
        peptide_shaker_log = RESULT_DIR / "logs/hostfiltering/RunPeptideShakerHost/{sample}/PeptideShaker.log",
    params:
        hostname = "host",
    conda:
        "../envs/host_filtering.yml"
    threads: workflow.cores / 2
    shell:
        "peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.hostname} -fasta_file {input.host_decoy_db} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule SimplePeptideListHost:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/SpectraFilter/host.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/SpectraFilter/host_Default_PSM_Report.txt",
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/SimplePeptideListHost/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/SimplePeptideListHost/{sample}/stdout.log"
    params:
        hostname = "host",
        out_dir = str(RESULT_DIR / "{sample}/SpectraFilter")
    conda:
        "../envs/host_filtering.yml"
    threads: 1
    shell:
        "peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 
   

rule FilterSpectra:
    input: 
        mgf = MGF_FILE,
        peptide_shaker_report = RESULT_DIR / "{sample}/SpectraFilter/host_Default_PSM_Report.txt",
    output: 
        filtered_mgf = RESULT_DIR / "{sample}/SpectraFilter/Filtered_host.mgf"
    log:
        stdout_log = RESULT_DIR / "logs/hostfiltering/FilterSpectra/{sample}/stdout.log"
    conda: 
        '../envs/pandas.yml'
    threads: 1
    script: 
        "../scripts/host_filtering.py"