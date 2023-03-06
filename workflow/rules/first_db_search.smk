rule AddDecoysRef:
    input: 
        ref = config["db_search"]["ref"],
    output: 
        ref_decoy_fasta = "resources/Database/refSeqViral_concatenated_target_decoy.fasta"
    log:
        stderr_log=RESULT_DIR / "logs/FirstSearch/RefDB/stderr.log",
        stdout_log=RESULT_DIR / "logs/FirstSearch/RefDB/stdout.log" 
    conda:
        "../envs/db_search.yml"
    threads: 1
    shell: 
        "searchgui eu.isas.searchgui.cmd.FastaCLI -in {input.ref} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule SearchSpectraAgainstReference:
    input: 
        ref_decoy_fasta = "resources/Database/refSeqViral_concatenated_target_decoy.fasta",
        mgf = SearchDB.get_input_SearchDB()["mfg"],
        par = PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/FirstSearch/SearchSpectraAgainstReference/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FirstSearch/SearchSpectraAgainstReference/{sample}/stdout.log"
    params:
        result_dir = str(RESULT_DIR / "{sample}/FirstSearch"),
        refname = "ref",
        search_engine = config["db_search"]["search_engine"],
        psm_fdr = config["db_search"]["psm_fdr"],
        peptide_fdr = config["db_search"]["peptide_fdr"],
        protein_fdr = config["db_search"]["protein_fdr"]
    conda:
        "../envs/db_search.yml"
    threads: workflow.cores / 2
    shell: 
        "searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.ref_decoy_fasta} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.refname}_searchgui_out -psm_fdr 1 -peptide_fdr 1 -protein_fdr 1 {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule RunPeptideShakerRef:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip",
        mgf = SearchDB.get_input_SearchDB()["mfg"],
        ref = config["db_search"]["ref"],
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/FirstSearch/ref.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/FirstSearch/RunPeptideShakerRef/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FirstSearch/RunPeptideShakerRef/{sample}/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/FirstSearch/RunPeptideShakerRef/{sample}/PeptideShaker.log",
    params:
        refname = "ref",
    conda:
        "../envs/db_search.yml"
    threads: workflow.cores / 2
    shell:
        "peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.ref} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule SimplePeptideListRef:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/FirstSearch/ref.psdb"
    output: 
        #TODO hostname varialbe: ResultsDir+SampleName+'/FirstSearch/'+HostName+'.psdb'
        peptide_shaker_report = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/FirstSearch/SimplePeptideListRef/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FirstSearch/SimplePeptideListRef/{sample}/stdout.log"
    params:
        out_dir = str(RESULT_DIR / "{sample}/FirstSearch")
    conda:
        "../envs/db_search.yml"
    threads: 1
    shell:
        "peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 
   
rule extractSearchGuiResults:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip",
    output:
        out_file = RESULT_DIR / "{sample}/FirstSearch/Filtered_host.t.xml.gz"
    params:
        out_dir = str(RESULT_DIR / "{sample}/FirstSearch")
    threads: 1
    shell:
        "unzip -u {input.searchgui_zip} -d {params.out_dir}"


rule createMS2RescoreConfig:
    input: 
        xTandem_out = RESULT_DIR / "{sample}/FirstSearch/Filtered_host.t.xml.gz"
    output:
        out_dir = RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/config.json"
    params:
        pipeline = config["ms2rescore"]["RescorePipeline"],
        features = config["ms2rescore"]["RescoreFeatures"],
        percolator = config["ms2rescore"]["RunPercolator"],
        fragmodel = config["ms2rescore"]["FragModel"]
    threads: 1
    script:
        "../scripts/ms2rescore_config.py"