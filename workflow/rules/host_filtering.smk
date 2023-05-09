# from PepGM, edited

rule AddHostandCrap:
    input:
        host_fasta = HOST_FASTA,
        crap = config["contaminants"]
    output:
        concat_db = RESULT_DIR / "{sample}/Database/Host_crap.fasta"
    shell:"cat {input} > {output}"  


rule AddDecoysHost:
    input: 
        concat_db = RESULT_DIR / "{sample}/Database/Host_crap.fasta"
    output: 
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/HostDB/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/HostDB/{sample}/stdout.log" 
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.concat_db} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule SearchHostSpectra:
    input: 
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta",
        par = PAR_FILE
    output:  
        out_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip"
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/SearchHostSpectra/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/SearchHostSpectra/{sample}/stdout.log"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/SpectraFilter"),
        hostname = "host",
        search_engine = config["host_filtering"]["search_engine"],
        psm_fdr = config["host_filtering"]["psm_fdr"],
        peptide_fdr = config["host_filtering"]["peptide_fdr"],
        protein_fdr = config["host_filtering"]["protein_fdr"]
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.decoy_crap_host} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.hostname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule RunPeptideShakerHost:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip",
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/SpectraFilter/host.psdb"
    log:
        stderr_log = RESULT_DIR / "logs/hostfiltering/RunPeptideShakerHost/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/hostfiltering/RunPeptideShakerHost/{sample}/stdout.log",
        peptide_shaker_log = RESULT_DIR / "logs/hostfiltering/RunPeptideShakerHost/{sample}/PeptideShaker.log",
    params:
        hostname = "host",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.hostname} -fasta_file {input.decoy_crap_host} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


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
        out_dir = str(RESULT_DIR / "{sample}/SpectraFilter"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 
   

rule FilterSpectra:
    input: 
        mgf = MGF_FILE,
        peptide_shaker_report = RESULT_DIR / "{sample}/SpectraFilter/host_Default_PSM_Report.txt",
    output: 
        filtered_mgf = RESULT_DIR / "{sample}/SpectraFilter/Filtered_host.mgf"
    log:
        stdout_log = RESULT_DIR / "logs/hostfiltering/FilterSpectra/{sample}/stdout.log"
    conda: 
        "../envs/pandas.yml"
    threads: 1
    script: 
        "../scripts/host_filtering.py"